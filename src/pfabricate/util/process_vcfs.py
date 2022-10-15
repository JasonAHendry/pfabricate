import subprocess


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
