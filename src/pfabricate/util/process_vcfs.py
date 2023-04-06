import subprocess
import pandas as pd
import numpy as np
from datetime import datetime
from dataclasses import dataclass


# --------------------------------------------------------------------------------
# bcftools interface
#
# --------------------------------------------------------------------------------


def bcftools(subcommand, args, input_vcf, run=False):
    """
    Construct `bcftools` commands from python

    """
    assert isinstance(subcommand, str)
    assert isinstance(args, list)
    assert isinstance(input_vcf, str)
    cmd = "bcftools"
    cmd += f" {subcommand}"
    cmd += f" {' '.join([str(a) for a in args])}"
    cmd += f" {input_vcf}"

    if not run:
        return cmd

    subprocess.run(cmd, shell=True, check=True)


# --------------------------------------------------------------------------------
# VCF building and writing
#
# --------------------------------------------------------------------------------


def convert_ad_array_to_string(ad):
    a = np.apply_along_axis(lambda vs: ",".join([str(v) for v in vs]), 1, ad)
    return a.tolist()


@dataclass
class VCFFormatType:
    ID: str
    Number: int
    Type: str
    Description: str

    def get_as_string(self):
        return f"##FORMAT=<ID={self.ID},Number={self.Number},Type={self.Type},Description='{self.Description}'>\n"


class VCFBuilder:
    variant_columns = [
        "CHROM",
        "POS",
        "ID",
        "REF",
        "ALT",
        "QUAL",
        "FILTER",
        "INFO",
        "FORMAT",
    ]
    variant_defaults = {
        "ID": ".",
        "QUAL": 10,
        "FILTER": "PASS",
        "INFO": ".",
        "FORMAT": ".",
    }
    format_types = [
        VCFFormatType("GT", 1, "String", "Genotype"),
        VCFFormatType("DP", 1, "Integer", "Read Depth"),
        VCFFormatType("AD", 2, "Integer", "Allelic depth"),
    ]

    def __init__(self):
        """
        Build a VCF in python by providing information about
        variants and samples

        """
        self.header = None
        self.variant_dict = {}
        self.sample_dict = {}

    def populate_variant_columns(self, **kwargs):
        """
        Populate values for columns that provide variant information

        These include `CHROM`, `POS`, `ID`, `REF`, `ALT`, `QUAL`
        `FILTER`, `INFO`, and `FORMAT`

        Invoke as follows:

            .populate_variant_colums(CHROM=<values>, POS=<values>, ...)

        """
        # Sanity check inputs
        assert "CHROM" in kwargs, "Must provide `CHROM` column."
        self.n_variants = len(kwargs["CHROM"])
        for key, values in kwargs.items():
            vs = np.array(values)
            assert np.array(values).shape == (
                self.n_variants,
            ), f"Column {key} has incorrect size. Check length and dimensions."
            assert (
                key in self.variant_columns
            ), f"Column {key} is not one of {','.join(self.variant_columns)}."

        # Populate `self.variant_dict`
        for column in self.variant_columns:
            if column in kwargs:
                self.variant_dict[column] = kwargs[column]
            elif column in self.variant_defaults:
                vs = [self.variant_defaults[column]] * self.n_variants
                self.variant_dict[column] = vs
            else:
                raise ValueError(f"Must provide column {column}.")

    def set_sample_format(self, **kwargs):
        """
        Define the format columns for each sample

        Invoke as follows:

            .set_sample_format(GT=True, DP=False, ...)

        """
        self._fmts = []
        for fmt in self.format_types:
            if fmt.ID in kwargs:
                if kwargs[fmt.ID]:
                    self._fmts.append(fmt.ID)
        self.sample_format_string = ":".join(self._fmts)
        self.variant_dict["FORMAT"] = [self.sample_format_string] * self.n_variants

    def add_sample(self, sample_name, **kwargs):
        """
        Add a sample to the VCF; ensuring it conforms to the expected
        FORMAT as defined by `.sample_format_string`

        Invoke as follows:

            .add_sample(sample_name="control1",
                        GT=["0/1", ...],
                        DP=[100, 110, ..]
                        )

        """
        assert (
            self.sample_format_string is not None
        ), "No sample format string. Run `.set_sample_format()`"

        sample_column = []
        for values in zip(*[kwargs[fmt] for fmt in self._fmts]):
            sample_column.append(":".join([str(v) for v in values]))
        self.sample_dict[sample_name] = sample_column

    def _create_header(self, source):
        """
        Write a VCF header string based on variant data and samples added

        """
        # Essential starting lines
        self.header = "##fileformat=VCFv4.2\n"
        self.header += f"##fileDate={datetime.now().strftime('%Y%m%d')}\n"
        self.header += f"##source={source}\n"

        # Compute some psuedo-contig information
        contig_names = set(self.variant_dict["CHROM"])
        contigs = []
        for n in contig_names:
            ca = np.array(self.variant_dict["CHROM"])
            pa = np.array(self.variant_dict["POS"])
            pac = pa[ca == n]
            l = pac.max() - pac.min()
            contigs.append((n, l))

        # Add contig information
        for n, l in contigs:
            self.header += f"##contig=<ID={n},length={l},nb='Computed from POS.'>\n"

        # Add format information
        for ID in self._fmts:
            fmt = [f for f in self.format_types if f.ID == ID][0]
            self.header += fmt.get_as_string()

        return self.header

    def get_vcf_as_dataframe(self):
        """
        Return the entire body of the VCF (i.e. without header)
        as a pandas dataframe

        """
        return pd.DataFrame({**self.variant_dict, **self.sample_dict})

    def write_vcf(self, 
                  output_path, 
                  source="VCFBuilder", 
                  compress=True,
                  index=True
                  ):
        """
        Write variant and sample into a VCF file at `output_path`

        """
        vcf_body = self.get_vcf_as_dataframe()
        with open(output_path, "w") as vcf:
            vcf.write(self._create_header(source=source))

            sep = "\t"
            vcf.write(f"#{sep.join(vcf_body.columns)}\n")
            for _, row in vcf_body.iterrows():
                vcf.write(f"{sep.join([str(v) for v in row.values])}\n")
        
        if compress:
            bcftools(subcommand="view", 
                     args=["-Ob", f"-o {output_path}.gz"],
                     input_vcf=output_path
                     )
            output_path += ".gz"
        
        if index:
            bcftools("index", args=[], input_vcf=output_path)

