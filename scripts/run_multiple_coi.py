import os
import re
import click
import pandas as pd
from dataclasses import dataclass
from pfabricate.util.generic import produce_dir
from pfabricate.mixed.main import mixed
from pfabricate.util.process_vcfs import  bcftools


# --------------------------------------------------------------------------------
# DEFAULT PARAMTERS
#
# --------------------------------------------------------------------------------


MAX_COI = 4
N_SIMULATE = 10
MAX_M = 3
DEPTH_MEAN = 50
DEPTH_SHAPE = 15
WSAF_SHAPE = 500
E_0 = 0.00001
E_1 = 0.0005


# --------------------------------------------------------------------------------
# MERGE FUNCTIONS
#
# --------------------------------------------------------------------------------


@dataclass(frozen=False)
class SimulatedFiles:
    output_dir: str
    vcf: str=""
    summary_csv: str=""
    segment_csv: str=""
    pairwise_csv: str=""

    def __post_init__(self):
        self.vcf = f"{self.output_dir}/simulated_infections.vcf.gz"
        self.summary_csv = f"{self.output_dir}/simulated_infections.summary.csv"
        self.segment_csv = f"{self.output_dir}/simulated_infections.ibd_segments.csv"
        self.pairwise_csv = f"{self.output_dir}/simulated_infections.ibd_pairwise.csv"

        for f in [self.output_dir, self.vcf, self.summary_csv, self.segment_csv, self.pairwise_csv]:
            assert os.path.exists(f)


# --------------------------------------------------------------------------------
# ENTRY POINT
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
    "-K", "--max_COI", type=int, required=True, default=MAX_COI, show_default=True, help="Complexity of infection to simulate."
)
@click.option(
    "-n",
    "--n_simulate",
    type=int,
    required=False,
    default=N_SIMULATE,
    show_default=True,
    help=f"Number of mixed infections to simulate.",
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
    "--merge",
    is_flag=True,
    show_default=True,
    default=False,
    help="Merge output files across COIs into a single set of files.",
)
def main(
    input_vcf,
    output_dir,
    max_coi,
    n_simulate,
    max_m,
    depth_mean,
    depth_shape,
    wsaf_shape,
    e_0,
    e_1,
    include_plots,
    merge
):
    """
    Simulate across multiple COI values, and optionally merge the final result
    
    """

    simulated_files = []
    for coi in range(1, max_coi + 1):
        print(f"Running for COI={coi}")
        
        coi_output_dir = f"{output_dir}/K{coi}"
        
        mixed(
            input_vcf,
            coi_output_dir,
            coi,
            n_simulate,
            max_m,
            depth_mean,
            depth_shape,
            wsaf_shape,
            e_0,
            e_1,
            include_plots,
        )

        simulated_files.append(SimulatedFiles(coi_output_dir))

    if merge:
        print("Merging...")
        merged_dir = produce_dir(output_dir, "merged")
        
        # VCF files
        bcftools(subcommand="merge",
                args=[f"-o {merged_dir}/simulated_infections.vcf.gz"],
                input_vcf=" ".join([f.vcf for f in simulated_files]),
                run=True
                )

        # IBD segments    
        merged_segment_df = pd.concat([pd.read_csv(f.segment_csv) for f in simulated_files])
        merged_segment_df.to_csv(f"{merged_dir}/simulated_infections.ibd_segments.csv", index=False)
        # IBD pairwise
        merged_pairwise_df = pd.concat([pd.read_csv(f.pairwise_csv) for f in simulated_files])
        merged_pairwise_df.to_csv(f"{merged_dir}/simulated_infections.ibd_pairwise.csv", index=False)

        # Summary CSVs
        summary_dfs = []
        for f in simulated_files:
            summary_df = pd.read_csv(f.summary_csv)
            merge_column_prefixes = ["samp[0-9]{2}", "prop[0-9]{2}"]
            for m in merge_column_prefixes:
                m_cols = [c for c in summary_df.columns if re.match(m, c)]
                m_df = summary_df[m_cols]
                summary_df[f"{m[:4]}s"] = [
                    ";".join([str(v) for v in row.values])
                    for _, row in m_df.iterrows()
                ]
                summary_df.drop(m_cols, axis=1, inplace=True)
            summary_dfs.append(summary_df)
        pd.concat(summary_dfs).to_csv(f"{merged_dir}/simulated_infections.summary.csv", index=False)
    print("Done.")


if __name__ == "__main__":
    main()
