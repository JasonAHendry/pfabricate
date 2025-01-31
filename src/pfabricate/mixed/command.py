import click


# --------------------------------------------------------------------------------
# DEFAULT PARAMTERS
#
# --------------------------------------------------------------------------------


N_SIMULATE = 20
MAX_K = 4
MAX_M = 3
DEPTH_MEAN = 50
DEPTH_SHAPE = 15
WSAF_SHAPE = 500
E_0 = 0.0001
E_1 = 0.005


# --------------------------------------------------------------------------------
# COMMAND
#
# --------------------------------------------------------------------------------


@click.command(short_help="Simulate a mixed infection.")
@click.option(
    "-i",
    "--input_vcf",
    type=str,
    required=True,
    help="Input VCF from which mixed infections simulated. Ideally contains only clonal samples.",
)
@click.option(
    "-o", "--output_dir", type=str, required=True, help="Path of output directory."
)
@click.option(
    "-K", "--max_K", type=int, default=MAX_K, show_default=True, help="Maximum complexity of infection to simulate."
)
@click.option(
    "-n",
    "--n_simulate",
    type=int,
    required=False,
    default=N_SIMULATE,
    show_default=True,
    help=f"Number of mixed infections to simulate for each COI.",
)
@click.option(
    "-M",
    "--max_M",
    type=int,
    required=False,
    default=MAX_M,
    show_default=True,
    help=f"Maximum number of meioses to simulate per mixed infection. The higher this value, the more IBD will be present in the mixed infections.",
)
@click.option(
    "-d",
    "--depth_mean",
    type=float,
    required=False,
    default=DEPTH_MEAN,
    show_default=True,
    help=f"Mean depth of simulated read data.",
)
@click.option(
    "-s",
    "--depth_shape",
    type=float,
    required=False,
    default=DEPTH_SHAPE,
    show_default=True,
    help=f"Shape parameter controlling variance in read depth. Feeds into negative binomial distribution; higher values reduce variance.",
)
@click.option(
    "-w",
    "--wsaf_shape",
    type=float,
    required=False,
    default=WSAF_SHAPE,
    show_default=True,
    help=f"Shape parameter controlling variance in WSAF. Feeds into betabinomial distribution; higher values reduce variance.",
)
@click.option(
    "--e_0",
    type=float,
    required=False,
    default=E_0,
    show_default=True,
    help=f"Error rate for conversion of reference calls to alternative.",
)
@click.option(
    "--e_1",
    type=float,
    required=False,
    default=E_1,
    show_default=True,
    help=f"Error rate for conversion of alternative calls to reference.",
)
@click.option(
    "--include_plots",
    is_flag=True,
    show_default=True,
    default=False,
    help="Create WSAF vs genomic position plots for every simulated infection.",
)
@click.option(
    "--id_prefix",
    type=str,
    required=False,
    default="SMI",
    help="Provide a custom ID prefix for simulated samples. [Optional]"
)
def mixed(
    input_vcf,
    output_dir,
    max_k,
    n_simulate,
    max_m,
    depth_mean,
    depth_shape,
    wsaf_shape,
    e_0,
    e_1,
    include_plots,
    id_prefix,
):
    """
    Stochastically simulate mixed infections from a VCF

    """
    from .main import mixed

    mixed(
        input_vcf,
        output_dir,
        max_k,
        n_simulate,
        max_m,
        depth_mean,
        depth_shape,
        wsaf_shape,
        e_0,
        e_1,
        include_plots,
        id_prefix
    )
