# pfabricate

`pfabricate` is a command-line tool for filtering and simulating from malaria VCF files. The overarching aim is to provide a framework to easily prepare data for COI estimation and benchmark the performance of COI estimation algorithms


## Install

`pfabricate` is implement in python and it's dependencies can be installed using conda. If haven't already, first install [conda](https://conda.io/projects/conda/en/latest/index.html). Then run:

```
conda update conda
conda env create
conda activate pfabricate
```
Now you have installed and activated a virtual environment containing all of `pfabricate`s dependencies. Finally, locally install `pfabricate`:
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

## TODO:
- DP filtering by percentile
- Assumptions associated with simulations
- Incorporation of background IBD
