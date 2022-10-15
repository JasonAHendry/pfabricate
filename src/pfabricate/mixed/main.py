import os
import numpy as np
from allel import read_vcf
from scipy.stats import dirichlet
from dataclasses import dataclass


def mixed(input_vcf, K):
    """
    Simulate mixed infections of COI=`K` from an `input_vcf`


    TODO:
    - What is the best way to handle `M`
        - What are the requirements?
    - How am I handling multiple chromosomes

    """

    # PARSE INPUTS

    # LOAD DATA
    vcf = read_vcf(input=input_vcf)
    samples = vcf["samples"]
    chroms = vcf["variants/CHROM"]
    pos = vcf["variants/POS"]
    gt = vcf["calldata/GT"].sum(2)

    # PARAMETERS
    # Data
    n_snps = pos.shape[0]
    n_samples = samples.shape[0]

    # Model
    # TODO: will probably come from arguments
    depth_mean = 500
    depth_shape = 20
    e_0 = 0.0001  # error rate, ref -> alt
    e_1 = 0.005  # error rate, alt -> ref

    # Simulation
    n_simulate = 20
    K = 2
    max_M = 3

    # SAMPLE
    # Proportions
    simulated_proportions = dirichlet.rvs(alpha=np.ones(K), size=n_simulate)

    # Bites
    simulated_bites = np.random.choice(a=np.arange(1, K), size=n_simulate, replace=True)

    # Meioses
    # For now, we will draw one per infection
    # - Should be zero if B == K
    # - Probably want to draw this later...
    simulated_meioses = np.random.choice(
        a=np.arange(1, max_M), size=n_simulate, replace=True
    )

    # SIMULATE
    @dataclass
    class SimulatedInfection:
        """
        Well... this needs to expand somehow
        depending on K

        """

        pass

    for i in range(n_simulate):

        # Sample K strains
        ixs = np.random.choice(a=n_samples, size=K, replace=False)

        # Pull parameters
        B = simulated_bites[i]
        M = simulated_meioses[i]
        inf_samples = samples[ixs]
        haplotypes = gt[:, ixs]  # I may want to transpose

        # may need to *haploidise* here

        # --------------------------------------------------------------------------------
        # Create haplotypes within the infection, by simulating
        # meiosis
        # --------------------------------------------------------------------------------

        if B == K:
            # No meiosis is required
            continue

        # Partition strains into bites,
        # simulate meiosis for each bite
        bites = partition_strains(K, B)
        transmitted = []
        for bite in bites:

            # If there is one strain, you don't perform meiosis
            if len(bite) == 1:
                continue

            bite_haplotypes = haplotypes[bite]

            # Here we run meiosis `M` times
            meiosis_engine = MeosisEngine(haplotypes=bite_haplotypes)
            meiosis_engine.run(rounds=M)

            # Here are the output progeny
            progeny = meiosis_engine.select_progeny(K=b)
            transmitted.append(progeny)

        # Finalise haplotypes
        infection_haplotypes = np.vstack(transmitted)

        # --------------------------------------------------------------------------------
        # Create haplotypes within the infection, by simulating
        # meiosis
        # --------------------------------------------------------------------------------

        read_simulator = ReadSimulator(
            haplotypes=infection_haplotypes,
            proportions=infection_proportions,
            mean_depth=mean_depth,
            e_0=e_0,
            e_1=e_1,
        )
        read_data = read_simulator.simulate()

        # So, now we have our ALT and REF read counts under the model
        # I can also (relatively) easily make the genotype calls
        #   - based on the infection haplotypes
        # So from here, I have a single-sample VCF

        # --------------------------------------------------------------------------------
        # Outputting 'truths'
        #
        # --------------------------------------------------------------------------------
        #
        # - So I know K and proportions, those are easy
        # - What I need to add to this is the IBD results
        # - So:
        # - (1) Pairwise IBD profiles between all strains (boolean matrix)
        # - (2) IBD tracks [truth tracks]
        # - (3) - Mean fraction, Length, and N50

        # What is the best structure of these?
        # How does hmmIBD output it's pairwise IBD data?
        # Probably want to encode as a table of segments, but look at other tools
        # /output_dir
        #   /parameters.json (or .log)
        #   /summary.csv
        #   /simulated_infections.vcf
        #   /simulated_ibd
        # /<inf_name>
        #     pairwise_profiles.npy
        #     pairwise_statistics.npy
