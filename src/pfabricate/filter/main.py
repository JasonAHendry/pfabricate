# Filter a VCF down to SNPs useful for COI inference
# 2022/10/12, JHendry
#
# TODO:
# - Depth based on percentiles
# - How does +prune actually work? What sites get dropped

import os
import re
import uuid
import shutil
import subprocess
import pandas as pd
from dataclasses import dataclass
from pfabricate.util.generic import (
    print_header,
    print_footer,
    produce_dir,
)
from pfabricate.util.process_vcfs import bcftools


# --------------------------------------------------------------------------------
# RECORD VCF STATISTICS
#
# --------------------------------------------------------------------------------


@dataclass
class ReportVCF:
    info: str
    vcf: str
    n_samples: int
    n_snps: int

    @staticmethod
    def _collect_bcftools_stats_summary(result, summary):
        """
        Extract a specific `summary` from a bcftools stats
        `result`

        """

        lines = result.stdout.decode("utf8").split("\n")
        n = [
            int(re.search("[0-9]+", s).group(0)) for s in lines if s.startswith(summary)
        ]

        return n[0] if n else None

    @classmethod
    def from_vcf(cls, vcf, info):
        """
        Extract VCF statistics from an `input_vcf`

        """

        # Run `bcftools stats`
        cmd_stats = bcftools("stats", args=[], input_vcf=vcf, run=False)
        cmd = f"{cmd_stats} | grep ^SN | cut -f 3-4"
        result = subprocess.run(cmd, check=True, shell=True, capture_output=True)

        # Collect specific statistic from result
        n_snps = ReportVCF._collect_bcftools_stats_summary(result, "number of SNPs")
        n_samples = ReportVCF._collect_bcftools_stats_summary(
            result, "number of samples"
        )

        return cls(info=info, vcf=vcf, n_samples=n_samples, n_snps=n_snps)


# --------------------------------------------------------------------------------
# FILTER VCFs
#
# --------------------------------------------------------------------------------


class FilterVCF:
    def __init__(self, input_vcf, output_vcf):
        """
        Filter an `input_vcf` to create an `output_vcf`
        through a series of steps each using bcftools

        """
        self.input_vcf = input_vcf
        self.input_fn = os.path.basename(self.input_vcf)
        self.output_vcf = output_vcf
        self.output_dir = produce_dir(os.path.dirname(self.output_vcf))
        self.step = 0

    def create_intermediate_dir(self):
        self.input_dir = os.path.dirname(self.input_vcf)
        self.inter_dir = produce_dir(
            self.input_dir, f"intermediates_{str(uuid.uuid4())[:8]}"
        )

    @property
    def current_vcf(self):
        self.current_fn = self.input_fn.replace(".vcf", f".filter{self.step:02d}.vcf")
        return f"{self.inter_dir}/{self.current_fn}"

    def bcftools(self, subcommand, args, final_step=False):
        if self.step == 0:
            iv = self.input_vcf
        else:
            iv = self.current_vcf
        if final_step:
            ov = self.output_vcf
        else:
            self.step += 1
            ov = self.current_vcf

        bcftools(input_vcf=iv, subcommand=subcommand, args=args + ["-o", ov], run=True)

    def clean_up(self):
        shutil.rmtree(self.inter_dir)


# --------------------------------------------------------------------------------
# MAIN SCRIPT
#
# --------------------------------------------------------------------------------


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
    Filter a VCF file stepwise


    """

    # PARSE INPUTS
    t0 = print_header("PFABRICATE: Filter to SNPs informative for COI inference")
    print("Parameters")
    print(f"  Input VCF: {input_vcf}")
    print(f"  Output VCF: {output_vcf}")
    if samples_file is not None:
        print(f"  Target sample file: {samples_file}")
    elif samples is not None:
        print(f"  Target samples include: {samples[:20]}...")
    else:
        print("  Including all samples.")
    print(f"  Minimum PLAF: {min_plaf}")
    print(f"  Maximum PLAF: {max_plaf}")
    print(f"  Minimum depth: {min_avg_depth}")
    print(f"  Maximum depth: {max_avg_depth}")
    print(f"  Max linkage (r2): {max_linkage}")
    print(f"  Window size for linkage (bp): {window_size}")
    print(f"  Maximum frac. missing: {max_frac_missing}")
    print("Done.\n")

    # PREPARE
    print("Computing statistics of input VCF...")
    vcf_reports = [ReportVCF.from_vcf(vcf=input_vcf, info="input")]
    vcf_filter = FilterVCF(input_vcf=input_vcf, output_vcf=output_vcf)
    vcf_filter.create_intermediate_dir()

    # FILTERS
    # (1) Samples
    if samples_file is not None:
        print("Filtering to indicated samples...")
        vcf_filter.bcftools(subcommand="view", args=["-S", samples_file])
        vcf_reports.append(ReportVCF.from_vcf(vcf_filter.current_vcf, "samples"))
    elif samples is not None:
        print("Filtering to indicated samples...")
        vcf_filter.bcftools(subcommand="view", args=["-s", samples])
        vcf_reports.append(ReportVCF.from_vcf(vcf_filter.current_vcf, "samples"))

    # (2) Biallelic SNPs
    print("Filtering to high-quality biallelic SNPs...")
    vcf_filter.bcftools(
        "view",
        args=[
            "--apply-filters",
            "PASS",
            "--min-alleles",
            "2",
            "--max-alleles",
            "2",
            "--types",
            "snps",
            "--output-type",
            "z",
        ],
    )
    vcf_reports.append(ReportVCF.from_vcf(vcf_filter.current_vcf, "biallelic"))

    # (3) PLAF
    print("Filtering to SNPs within PLAF range...")
    vcf_filter.bcftools(
        subcommand="view",
        args=["--min-af", min_plaf, "--max-af", max_plaf],
    )
    vcf_reports.append(ReportVCF.from_vcf(vcf_filter.current_vcf, "plaf"))

    # (4) Depth
    print("Filtering to SNPs within average read depth range...")
    vcf_filter.bcftools(
        subcommand="view",
        args=[
            "-i",
            f"'MEAN(FORMAT/DP) >= {min_avg_depth} & MEAN(FORMAT/DP) <= {max_avg_depth}'",
        ],
    )
    vcf_reports.append(ReportVCF.from_vcf(vcf_filter.current_vcf, "depth"))

    # (5) Linkage disequilibrium
    print("Filtering linked and missing SNPs...")
    vcf_filter.bcftools(
        subcommand="+prune",
        args=[
            f"-m {max_linkage}",
            f"-w {window_size}",
            f"-e 'F_MISSING>={max_frac_missing}'",
        ],
        final_step=True,
    )
    vcf_reports.append(ReportVCF.from_vcf(vcf_filter.output_vcf, "linkage"))
    print("Done.\n")
    vcf_filter.clean_up()

    # PRINT RESULTS
    report_df = pd.DataFrame(vcf_reports)
    report_df.rename({"info": "filter_step"}, inplace=True)
    report_df.insert(0, "step", range(report_df.shape[0]))
    output_report = output_vcf.replace(".gz", "").replace(".vcf", ".report.csv")
    report_df.to_csv(output_report, index=False)

    print_footer(t0)
