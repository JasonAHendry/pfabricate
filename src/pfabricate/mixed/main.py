import uuid
import random
import numpy as np
import pandas as pd
from allel import read_vcf
from scipy.stats import dirichlet
from pfabricate.util.generic import produce_dir
from pfabricate.util.ibd import PairwiseIBDStatistics
from pfabricate.util.process_vcfs import VCFBuilder
from pfabricate.util.plotting import WSAFPlotter
from .meiosis import MeiosisEngine, OocystMaker, partition_strains
from .sequencing import convert_ad_to_haplotypes, simulate_read_data
from .chromosomes import ChromosomeFactory


def mixed(
    input_vcf,
    output_dir,
    coi,
    n_simulate,
    max_M,
    depth_mean,
    depth_shape,
    wsaf_shape,
    e_0,
    e_1,
    include_plots,
):
    """
    Stochastically simulate mixed infections of COI=`K` from an `input_vcf`

    Model:

    - Read depth from a negative binomial
    - WSAF from a betabinomial
    - Print statements


    TODO:
    - What is the best way to handle `M`
        - What are the requirements?
    - Provide access to dirichlet shape parameter
    - Vastly improve storage
    - Better handle cases where no IBD is generated

    """

    # Prepare parameters
    K = coi
    output_dir = produce_dir(output_dir)

    # Load VCF
    vcf = read_vcf(
        input=input_vcf,
        fields=[
            "samples",
            "variants/CHROM",
            "variants/POS",
            "variants/REF",
            "variants/ALT",
            "calldata/AD",
        ],
    )

    # Extract relevant fields
    samples = vcf["samples"]
    chroms = vcf["variants/CHROM"]
    pos = vcf["variants/POS"]
    ref = vcf["variants/REF"]
    alt = vcf["variants/ALT"]
    ad = vcf["calldata/AD"]

    if include_plots:
        plot_dir = produce_dir(output_dir, "plots")
        wsaf_plotter = WSAFPlotter()
        wsaf_plotter.set_snp_positions(chroms=chroms, pos=pos)

    # Create haplotypes from allelic depth
    # TODO: create them instead from GT
    haplotypes = convert_ad_to_haplotypes(ad)

    n_snps = pos.shape[0]
    n_samples = samples.shape[0]

    chrom_factory = ChromosomeFactory(chroms=chroms, pos=pos)
    chromosomes = chrom_factory.create_chromosomes()
    genome_kbp = sum([c.l_kbp for c in chromosomes])

    oocyst_maker = OocystMaker()

    # SAMPLE
    simulated_proportions = dirichlet.rvs(alpha=np.ones(K), size=n_simulate)
    simulated_proportions.sort(axis=1)  # sort for comparison later
    simulated_bites = np.random.choice(a=range(1, K + 1), size=n_simulate, replace=True)
    simulated_meioses = np.random.choice(
        a=range(1, max_M + 1), size=n_simulate, replace=True
    )
    simulated_meioses[simulated_bites == K] = 0

    # STORAGE
    # VCF
    vcf_builder = VCFBuilder()
    vcf_builder.populate_variant_columns(
        CHROM=chroms, POS=pos, REF=ref, ALT=alt[:, 0]  # assuming biallelic
    )
    vcf_builder.set_sample_format(AD=True, DP=True, GT=True)

    # Summary table
    sample_dt = {f"samp{i:02d}": [] for i in range(K)}
    prop_dt = {f"prop{i:02d}": simulated_proportions[:, i] for i in range(K)}
    statistic_dt = {
        "K": np.repeat(K, n_simulate),
        "Keff": 1 / (simulated_proportions**2).sum(1),
        "B": simulated_bites,
        "M": simulated_meioses,
        "f_ibd": np.zeros(n_simulate),
        "l_ibd": np.zeros(n_simulate),
        "n_ibd": np.zeros(n_simulate),
        "n50_ibd": np.zeros(n_simulate),
    }
    summary_dt = {
        "sample_id": [
            f"SMI{str(uuid.uuid4())[:6].upper()}-K{K:02d}-B{B:02d}-M{M:02d}"
            for i, (B, M) in enumerate(zip(simulated_bites, simulated_meioses))
        ]
    }
    summary_dt.update(statistic_dt)
    summary_dt.update(sample_dt)
    summary_dt.update(prop_dt)

    # IBD segments
    expected_columns = ["sample_id", "strain1", "strain2", "chrom", "start", "end", "length"]
    ibd_dfs = [pd.DataFrame(columns=expected_columns)]   # handles cases with no IBD
    ibd_pairwise_dfs = [pd.DataFrame(columns=["sample_id", "strain1", "strain2", "f_ibd", "l_ibd", "n50_ibd", "n_ibd"])]
    for i in range(n_simulate):

        # Sample K strains at random
        ixs = random.sample(range(n_samples), k=K)

        # Extract infection information
        sample_id = summary_dt["sample_id"][i]
        props = simulated_proportions[i]
        B = simulated_bites[i]
        M = simulated_meioses[i]
        infection_samples = samples[ixs]
        infection_haplotypes = haplotypes[:, ixs].transpose()

        # Store samples
        for j in range(K):
            summary_dt[f"samp{j:02d}"].append(infection_samples[j])

        if B < K:  # co-transmission

            bites = partition_strains(K, B)

            transmitted = []
            ibd_segs = []
            ibd_pairs = []
            for bite in bites:

                if len(bite) == 1:  # if one strain, no meiosis
                    transmitted.append(infection_haplotypes[bite])
                    continue

                # Run meioses 'M' times
                meiosis_engine = MeiosisEngine(
                    haplotypes=infection_haplotypes[bite],
                    chroms=chroms,
                    pos=pos,
                    oocyst_maker=oocyst_maker,
                )
                meiosis_engine.run(
                    n_rounds=M, n_select=len(bite), force_outbred=True  # same in as out
                )

                # Get haplotypes
                bite_haplotypes = meiosis_engine.get_progeny_haplotypes()
                transmitted.append(bite_haplotypes)

                # Get IBD
                meiosis_engine.detect_generated_ibd()
                ibd_segs.append(meiosis_engine.get_ibd_segment_dataframe())
                ibd_pairs.append(meiosis_engine.get_ibd_pairwise_dataframe())

            # Combine across bites
            infection_haplotypes = np.vstack(transmitted)
            # IBD Segments
            ibd_df = pd.concat(ibd_segs)
            ibd_df.insert(0, "sample_id", sample_id)
            ibd_dfs.append(ibd_df)
            # Ibd Pairwise
            ibd_pair_df = pd.concat(ibd_pairs)
            ibd_pair_df.insert(0, "sample_id", sample_id)
            ibd_pairwise_dfs.append(ibd_pair_df)

            # Compute summary statistics
            n_pairs = K * (K - 1) / 2
            total_genome_length = 10**3 * genome_kbp * n_pairs
            ibd_stats = PairwiseIBDStatistics.from_dataframe(ibd_df, total_genome_length)
            for s in ["f_ibd", "l_ibd", "n50_ibd", "n_ibd"]:
                summary_dt[s][i] = getattr(ibd_stats, s)

        # Simulate read data
        read_data = simulate_read_data(
            haplotypes=infection_haplotypes,
            proportions=props,
            depth_mean=depth_mean,
            depth_shape=depth_shape,
            alt_shape=wsaf_shape,
            e_0=e_0,
            e_1=e_1,
        )

        # TODO:
        # - I could also just specifically set the props / depth &c here
        # - Right here I need to include simulated genotype calls
        #  - These would be based on the input infection_haplotypes

        # Store
        vcf_builder.add_sample(
            sample_name=sample_id,
            DP=read_data["depth"],
            AD=[f"{r},{a}" for r, a in zip(read_data["ref"], read_data["alt"])],
            GT=read_data["gts"]
        )

        # Optionally plot
        if include_plots:
            wsaf = read_data["alt"] / read_data["depth"] + 0.01
            wsaf_plotter.plot(
                wsaf=wsaf,
                title=f"{sample_id} | $p=${', '.join([f'{p:.02f}' for p in props])} | $f=${summary_dt['f_ibd'][i]:02f}",
                output_path=f"{plot_dir}/plot.{sample_id}.wsaf.png",
            )

    # Write results
    # VCF
    vcf_builder.write_vcf(
        output_path=f"{output_dir}/simulated_infections.vcf", source="pfabricate"
    )
    # Summary
    summary_df = pd.DataFrame(summary_dt)
    summary_df.to_csv(f"{output_dir}/simulated_infections.summary.csv", index=False)
    # IBD
    combined_ibd = pd.concat(ibd_dfs)  # Get "no objects to concatenate, for COI=1"
    combined_ibd.to_csv(
        f"{output_dir}/simulated_infections.ibd_segments.csv", index=False
    )
    pd.concat(ibd_pairwise_dfs).to_csv(f"{output_dir}/simulated_infections.ibd_pairwise.csv", index=False)

