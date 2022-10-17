import os
import click
import pandas as pd
from dataclasses import dataclass
from pfabricate.util.generic import produce_dir


SOURCE_COMMAND = "qsub cluster/bmrc_pfabricate.sh"
CHROMOSOMES = range(1, 15)


# --------------------------------------------------------------------------------
# DEFINE DATASETS
#
# --------------------------------------------------------------------------------


@dataclass
class DatasetVCF:
    name: str
    metadata: str
    vcf_dir: str
    vcf_template: str

    def get_chrom_vcf(self, chrom_int):
        return f"{self.vcf_dir}/{self.vcf_template.format(chrom=chrom_int)}"


dataset_pf3k = DatasetVCF(
    name="pf3k",
    metadata="data/pf3k/metadata/pf3k_release_5_metadata_20190910_cleaned.csv",
    vcf_dir="data/pf3k/vcfs",
    vcf_template="SNP_INDEL_Pf3D7_{chrom:02d}_v3.high_quality_biallelic_snps.vcf.gz",
)

dataset_pf6 = DatasetVCF(
    name="pf6",
    metadata="data/pf6/metadata/Pf_60_candidate_public_metadata_20170919_cleaned.csv",
    vcf_dir="data/pf6/vcfs",
    vcf_template="Pf_60_public_Pf3D7_{chrom:02d}_v3.high_quality_biallelic_snps.vcf.gz",
)

DATASET_COLLECTION = {"pf3k": dataset_pf3k, "pf6": dataset_pf6}


# --------------------------------------------------------------------------------
# Extracting samples
#
# --------------------------------------------------------------------------------


def get_samples_in_region(metadata, region_dt):
    qry = " and ".join([f"{k} == " + (f"{v}" if isinstance(v, int) else f"'{v}'")
                    for k, v in region_dt.items()])
    df = metadata.query(qry)
    return df["sample"]


# --------------------------------------------------------------------------------
# MAIN
#
# --------------------------------------------------------------------------------


@click.command()
@click.option(
    "-d",
    "--dataset",
    type=click.Choice(DATASET_COLLECTION),
    required=True,
    help="MalariaGEN dataset.",
)
@click.option("-p", "--population", type=str, default=None, help="Target population.")
@click.option("-c", "--country", type=str, default=None, help="Target country.")
@click.option("-s", "--site", type=str, default=None, help="Target site.")
@click.option("-y", "--year", type=int, default=None, help="Target year.")
def main(dataset, population, country, site, year):
    """
    Create a BMRC cluster submission script for pfabricate, grouping
    by sample metadata features

    """
    # Parse inputs
    vcf_dataset = DATASET_COLLECTION[dataset]
    print("Inputs")
    print(f"  Dataset: {dataset}")
    print(f"  Population: {population}")
    print(f"  Country: {country}")
    print(f"  Site: {site}")
    print(f"  Year: {year}")

    # Prepare the output directories
    ps = {"population": population, "country": country, "site": site, "year": year}
    ps = {p: v for p, v in ps.items() if v}
    if not ps:
        raise ValueError("At least one of region flags must be set.")
    region = "".join([v for _, v in ps.items()])
    grouped_by = f"by_{'-'.join([p for p in ps])}"
    region_dir = produce_dir("results", dataset, grouped_by, region)
    vcf_dir = produce_dir(region_dir, "vcfs")
    print("Directories")
    print(f"  Region: {region_dir}")
    print(f"  VCF: {vcf_dir}")
    print("Samples")
    metadata_df = pd.read_csv(vcf_dataset.metadata)
    samples = get_samples_in_region(metadata_df, ps)
    sample_str = ",".join(samples)
    print(f"  No: {len(samples)}")
    print(f"  E.g.: {','.join(samples[:5])}...")

    # Define submission script
    submit_fn = f"submit_pfabricate-filter-{region}.sh"
    print(f"Writing submission script to:\n  {submit_fn}")
    with open(submit_fn, "w") as file:
        for c in CHROMOSOMES:
            in_vcf = vcf_dataset.get_chrom_vcf(c)
            out_vcf = (
                f"{vcf_dir}/{os.path.basename(in_vcf).replace('.vcf', '.filter.vcf')}"
            )
            file.write(f"{SOURCE_COMMAND} filter -i {in_vcf} -o {out_vcf} -s {sample_str}\n")
    print("Done.\n")
    os.chmod(submit_fn, 0o775)


if __name__ == "__main__":
    main()
