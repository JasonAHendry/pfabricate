import click


@click.command(short_help="Simulate a mixed infection.")
@click.option("-i", "--input_vcf", type=str, required=True, help="Input VCF to filter.")
@click.option("-o", "--output_vcf", type=str, required=True, help="Name of output VCF.")
def mixed(input_vcf, output_vcf):
    """
    Simulate a mixed infection

    """
    from .main import main

    main(
        input_vcf=input_vcf,
        output_vcf=output_vcf,
    )
