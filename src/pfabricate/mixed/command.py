import click


# --------------------------------------------------------------------------------
# DEFAULT PARAMTERS
#
# --------------------------------------------------------------------------------


N_SIMULATE = 20
MAX_M = 3
DEPTH_MEAN = 50
DEPTH_SHAPE = 100
E_0 = 0.0001
E_1 = 0.005


# --------------------------------------------------------------------------------
# COMMAND
#
# --------------------------------------------------------------------------------


@click.command(short_help="Simulate a mixed infection.")
@click.option("-i", "--input_vcf", type=str, required=True, help="Input VCF from which mixed infections simulated. Ideally contains only clonal samples.")
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
    required=False,
    default=N_SIMULATE,
    help=f"Number of mixed infections to simulate. Default={N_SIMULATE}.",
)
@click.option(
    "-M",
    "--max_M",
    type=int,
    required=False,
    default=MAX_M,
    help=f"Maximum number of meioses to simulate per mixed infection. The higher this value, the more IBD will be present in the mixed infections. Default={MAX_M}.",
)
@click.option(
    "-d",
    "--depth_mean",
    type=float,
    required=False,
    default=DEPTH_MEAN,
    help=f"Mean depth of simulated read data. Default={DEPTH_MEAN}.",
)
@click.option(
    "-s",
    "--depth_shape",
    type=float,
    required=False,
    default=DEPTH_SHAPE,
    help=f"Shape parameter controlling noise in read depth and WSAF values. Feeds into betabinomial distribution; higher values reduce variance. Default={DEPTH_SHAPE}.",
)
@click.option(
    "--e_0",
    type=float,
    required=False,
    default=E_0,
    help=f"Error rate for conversion of reference calls to alternative. Default={E_0}.",
)
@click.option(
    "--e_1",
    type=float,
    required=False,
    default=E_1,
    help=f"Error rate for conversion of alternative calls to reference. Default={E_1}.",
)
def mixed(
    input_vcf, output_dir, coi, n_simulate, max_m, depth_mean, depth_shape, e_0, e_1
):
    """
    Stochastically simulate mixed infections of a particular COI from a VCF

    """
    from .main import mixed

    mixed(
        input_vcf, output_dir, coi, n_simulate, max_m, depth_mean, depth_shape, e_0, e_1
    )
