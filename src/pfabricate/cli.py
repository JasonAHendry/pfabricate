import click
from .filter.command import filter
from .mixed.command import mixed


@click.group()
def cli():
    """
    Simulate P. falciparum mixed infections from VCF files

    """
    pass


cli.add_command(filter)
cli.add_command(mixed)

if __name__ == "__main__":
    cli()
