import click


# --------------------------------------------------------------------------------
# DEFAULT PARAMTERS
#
# --------------------------------------------------------------------------------


# MIN_PLAF = 0.1
# MAX_PLAF = 0.9
# MIN_AVG_DEPTH = 20
# MAX_AVG_DEPTH = 150
# MAX_LINKAGE = 0.2
# WINDOW_SIZE = 4000
# MAX_FRAC_MISSING = 0.05


# --------------------------------------------------------------------------------
# COMMAND
#
# --------------------------------------------------------------------------------


@click.command(short_help="Simulate a mixed infection.")
@click.option("-i", "--input_vcf", type=str, required=True, help="Input VCF to filter.")
@click.option(
    "-o", "--output_dir", type=str, required=True, help="Path of output directory."
)
@click.option(
    "-K", "--COI", type=int, required=True, help="Complexity of infection to simulate."
)
@click.option(
    "-n",
    "--n_simulate",
    type=int,
    required=True,
    default=20,
    help="Number of mixed infections to simulate.",
)
@click.option(
    "-M",
    "--max_M",
    type=int,
    required=False,
    default=3,
    help="Maximum number of meioses to simulate per mixed infection. The higher this value, the more IBD will be present in the mixed infections.",
)
@click.option(
    "-d",
    "--depth_mean",
    type=float,
    required=False,
    default=50,
    help="Mean depth of simulated read data.",
)
@click.option(
    "-s",
    "--depth_shape",
    type=float,
    required=False,
    default=100,
    help="Shape parameter controlling noise in read depth and WSAF values. Feeds into betabinomial distribution; higher values reduce variance.",
)
@click.option(
    "--e_0",
    type=float,
    required=False,
    default=0.01,
    help="Error rate for conversion of reference calls to alternative.",
)
@click.option(
    "--e_1",
    type=float,
    required=False,
    default=0.01,
    help="Error rate for conversion of alternative calls to reference.",
)
def mixed(
    input_vcf, output_dir, K, n_simulate, max_M, depth_mean, depth_shape, e_0, e_1
):
    """
    Stochastically simulate mixed infections of a particular COI from a VCF

    """
    from .main import mixed

    mixed(
        input_vcf, output_dir, K, n_simulate, max_M, depth_mean, depth_shape, e_0, e_1
    )
