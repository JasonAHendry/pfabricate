import uuid
import random
import numpy as np
import pandas as pd
from allel import read_vcf
from scipy.stats import dirichlet
from pfabricate.util.generic import produce_dir
from pfabricate.util.process_vcfs import VCFBuilder
from pfabricate.util.plotting import WSAFPlotter
from pfabricate.mixed.infection import MalariaInfectionMaker
from pfabricate.mixed.sequencing import convert_ad_to_haplotypes, simulate_read_data


def mixed(
    input_vcf: str,
    output_dir: str,
    max_K: int,
    n_simulate: int,
    max_M: int,
    depth_mean: int,
    depth_shape: float,
    wsaf_shape: float,
    e_0: float,
    e_1: float,
    include_plots: bool,
    id_prefix
):
    """
    Stochastically simulate mixed infections from an `input_vcf`

    """

    print(f"Input VCF: {input_vcf}")
    print(f"Max COI (K): {max_K}")
    print(f"No. simulations per COI: {n_simulate}")
    print(f"Max meioses (M): {max_M}")

    # Prepare parameters
    output_dir = produce_dir(output_dir)

    # Load VCF
    vcf = read_vcf(
        input=input_vcf,
        fields=["samples", "variants/CHROM", "variants/POS", "variants/REF", "variants/ALT", "calldata/AD"],
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
    n_samples = samples.shape[0]

    # STORAGE
    # VCF
    vcf_builder = VCFBuilder()
    vcf_builder.populate_variant_columns(
        CHROM=chroms, POS=pos, REF=ref, ALT=alt[:, 0]  # assuming biallelic
    )
    vcf_builder.set_sample_format(AD=True, DP=True, GT=True)

    # IBD and Summary
    outputs = {
        "summary": [],
        "ibd_pairwise": [],
        "ibd_segments": []
    }

    # Iterate over K and replicates
    rep_ix = 0
    for K in range(1, max_K + 1):

        print(f" Simulating for COI={K}...")
        
        # Sample infection parameters for all samples of this COI
        simulated_proportions = dirichlet.rvs(alpha=np.ones(K), size=n_simulate)
        simulated_proportions.sort(axis=1)  # sort for comparison later
        simulated_bites = np.random.choice(a=range(1, K + 1), size=n_simulate, replace=True)
        simulated_meioses = np.random.choice(
            a=range(1, max_M + 1), size=n_simulate, replace=True
        )
        simulated_meioses[simulated_bites == K] = 0
        
        # Simulate an individual sample
        for i in range(n_simulate):

            # Sample K strains, and get parameters
            ixs = random.sample(range(n_samples), k=K)
            props = simulated_proportions[i]
            B = simulated_bites[i]
            M = simulated_meioses[i]
            infection_samples = samples[ixs]
            infection_haplotypes = haplotypes[:, ixs].transpose()

            # Simulate the infection
            infection = MalariaInfectionMaker(
                input_haplotypes=infection_haplotypes,
                chroms=chroms,
                pos=pos,
                B=B,
                M=M,
                props=props,
                input_names=infection_samples
            )
            infection.assign_sample_id(id_prefix=id_prefix, index=rep_ix)
            infection.generate()

            # Simulate read data
            read_data = simulate_read_data(
                haplotypes=infection.get_infection_haplotypes(),
                proportions=props,
                depth_mean=depth_mean,
                depth_shape=depth_shape,
                alt_shape=wsaf_shape,
                e_0=e_0,
                e_1=e_1,
            )

            # Store outputs
            outputs["summary"].append(infection.get_sample_summary())
            outputs["ibd_pairwise"].append(infection.get_ibd_pairwise())
            outputs["ibd_segments"].append(infection.get_ibd_segments())

            vcf_builder.add_sample(
                sample_name=infection.sample_id,
                DP=read_data["depth"],
                AD=[f"{r},{a}" for r, a in zip(read_data["ref"], read_data["alt"])],
                GT=read_data["gts"]
            )

            # Optionally plot
            if include_plots:
                wsaf = read_data["alt"] / read_data["depth"] + 0.01
                wsaf_plotter.plot(
                    wsaf=wsaf,
                    title=f"{infection.sample_id} | $p=${', '.join([f'{p:.02f}' for p in props])} | $f=${outputs['summary'][-1]['f_ibd']:02f}",
                    output_path=f"{plot_dir}/plot.{infection.sample_id}.wsaf.png",
                )
            rep_ix += 1

    # Write results
    print("Writing outputs...")
    print(f"  to: {output_dir}")
    # VCF
    vcf_builder.write_vcf(
        output_path=f"{output_dir}/simulated_infections.vcf", source="pfabricate"
    )
    # Summary
    pd.DataFrame(outputs["summary"]).to_csv(f"{output_dir}/simulated_infections.summary.csv", index=False)

    # IBD
    for ibd_summary in ["ibd_pairwise", "ibd_segments"]:
        pd.concat(outputs[ibd_summary]).to_csv(f"{output_dir}/simulated_infections.{ibd_summary}.csv", index=False)
    print("Done.")

