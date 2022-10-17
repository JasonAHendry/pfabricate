import click


# --------------------------------------------------------------------------------
# DEFAULT PARAMTERS
#
# --------------------------------------------------------------------------------


MIN_PLAF = 0.1
MAX_PLAF = 0.9
MIN_AVG_DEPTH = 20
MAX_AVG_DEPTH = 150
MAX_LINKAGE = 0.2
WINDOW_SIZE = 4000
MAX_FRAC_MISSING = 0.05


# --------------------------------------------------------------------------------
# COMMAND
#
# --------------------------------------------------------------------------------


@click.command(short_help="Filter to informative SNPs.")
@click.option("-i", "--input_vcf", type=str, required=True, help="Input VCF to filter.")
@click.option("-o", "--output_vcf", type=str, required=True, help="Name of output VCF.")
@click.option(
    "-s",
    "--samples",
    type=str,
    default=None,
    help=f"Comma separated list of samples to include in population. [optional]",
)
@click.option(
    "-S",
    "--samples_file",
    type=str,
    default=None,
    help=f"Text file listing samples to include in population. Same as `bcftools view -S`. [optional]",
)
@click.option(
    "-p",
    "--min_plaf",
    type=float,
    default=MIN_PLAF,
    help=f"Filter SNPs with population-level allele frequency less than this value. Default={MIN_PLAF}.",
)
@click.option(
    "-P",
    "--max_plaf",
    type=float,
    default=MAX_PLAF,
    help=f"Filter SNPs with population-level allele frequency greater than this value. Default={MAX_PLAF}.",
)
@click.option(
    "-d",
    "--min_avg_depth",
    type=float,
    default=MIN_AVG_DEPTH,
    help=f"Filter SNPs with a mean read depth across samples of < `--min_avg_depth`. Default={MIN_AVG_DEPTH}.",
)
@click.option(
    "-D",
    "--max_avg_depth",
    type=float,
    default=MAX_AVG_DEPTH,
    help=f"Filter SNPs with a mean read depth across samples of > `--max_avg_depth`. Default={MAX_AVG_DEPTH}.",
)
@click.option(
    "-L",
    "--max_linkage",
    type=float,
    default=MAX_LINKAGE,
    help=f"Filter SNPs with a r2 > `-L` within a window size `-w`. Default={MAX_LINKAGE}.",
)
@click.option(
    "-w",
    "--window_size",
    type=float,
    default=WINDOW_SIZE,
    help=f"Apply linkage filter to SNPs within `-w` bp of each other. Default={WINDOW_SIZE}.",
)
@click.option(
    "-F",
    "--max_frac_missing",
    type=float,
    default=MAX_FRAC_MISSING,
    help=f"Remove SNPs missing in a fraction > `-F` of samples. Default={MAX_FRAC_MISSING}.",
)
def filter(
    input_vcf,
    output_vcf,
    samples,
    samples_file,
    min_plaf,
    max_plaf,
    min_avg_depth,
    max_avg_depth,
    max_linkage,
    window_size,
    max_frac_missing,
):
    """
    Filter to a reduced set of biallelic SNPs most
    informative for COI estimation

    """
    from .main import filter

    filter(
        input_vcf=input_vcf,
        output_vcf=output_vcf,
        samples=samples,
        samples_file=samples_file,
        min_plaf=min_plaf,
        max_plaf=max_plaf,
        min_avg_depth=min_avg_depth,
        max_avg_depth=max_avg_depth,
        max_linkage=max_linkage,
        window_size=window_size,
        max_frac_missing=max_frac_missing,
    )
