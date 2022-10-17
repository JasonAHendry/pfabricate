# pfabricate

`pfabricate` is a command-line tool for filtering and simulating from malaria VCF files. The overarching aim is to provide a framework to easily prepare data for COI estimation and benchmark the performance of COI estimation algorithms.


## Install

`pfabricate` is implemented in python and has dependencies that can be installed using [conda](https://conda.io/projects/conda/en/latest/index.html):

With conda installed, run:

```

conda update conda
conda env create
conda activate pfabricate

```

Now you have installed and activated a virtual environment containing all necessary dependencies. Next, locally install `pfabricate`:

```
pip install .
```
The command `pfabricate` should now be accessible.

## Usage
```
Usage: pfabricate [OPTIONS] COMMAND [ARGS]...

  Simulate P. falciparum mixed infections from VCF files

Options:
  --help  Show this message and exit.

Commands:
  filter  Filter to informative SNPs.
  mixed   Simulate a mixed infection.
```

## TODO
- DP filtering by percentile
- Assumptions associated with simulations
- Incorporation of background IBD
- Handling cases where no IBD is generated
- Most accurate approach for IBD fragment length distances given SNP positions
- Reconciling K > 4 with recombinant progeny
