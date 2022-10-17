import click
import pandas as pd


# --------------------------------------------------------------------------------
# GLOBAL CONSTANTS
#
# --------------------------------------------------------------------------------


DATA_DIR = "data/pf3k"
METADATA_DIR = "data/pf3k/metadata"
METADATA_PATH = f"{METADATA_DIR}/pf3k_release_5_metadata_20190910_cleaned.csv"


# --------------------------------------------------------------------------------
# ENUMERATE COUTRIES
#
# --------------------------------------------------------------------------------


METADATA_DF = pd.read_csv(METADATA_PATH)
COUNTRY_COLLECTION = {
    country: country_df["sample"].values.tolist()
    for country, country_df in METADATA_DF.groupby("country")
}


# --------------------------------------------------------------------------------
# SCRIPT SELECTION
#
# --------------------------------------------------------------------------------


@click.command()
@click.option(
    "-c",
    "--country",
    type=click.Choice(COUNTRY_COLLECTION),
    required=True,
    help="Output a text file listing all samples belonging to a specific country.",
)
def main(country):
    """
    Create a text file containing a set of sample names
    for a given country in Pf6

    """
    samples = COUNTRY_COLLECTION[country]
    output_path = f"{METADATA_DIR}/samples.{country}.txt"
    with open(output_path, "w") as file:
        file.write("\n".join(samples))
    print(f"Getting samples from Pf6...")
    print(f"  Country: {country}")
    print(f"  No. samples: {len(samples)}")
    print(f"  Written to: {output_path}")
    print(f"Done.\n")


if __name__ == "__main__":
    main()
