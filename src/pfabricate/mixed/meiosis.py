import random
import numpy as np
import pandas as pd
from typing import List
from itertools import combinations
from pfabricate.mixed.chromosomes import ChromosomeFactory
from pfabricate.util.ibd import (
    create_pairwise_ibd_matrix, 
    get_ibd_segments,
    IBDSummaryStatistics
)


# TODO
# - Needs better documentation
# - Bonus would be to do a bit of typing stuff
# - And they are nicely broken up so should be quite easy to write some tests


# --------------------------------------------------------------------------------
# Trasmission function
#
# --------------------------------------------------------------------------------


def partition_strains(K: int, n_bites: int) -> List:
    """
    Partition `K` strains randomly into
    `n_bites`

    TODO: 
    - What is the affect of this, given that we are
    *sorting* the strains by proportion at the top?
    - Probably should start from an empty list

    params
        K : int
            Number of strains in the infection.
        n_bites : int
            Number of bites by which the `K`
            strains are delivered to the infection.

    returns
        bites : list of lists, int, len(n_bites)
            A list of `n_bites` lists, each sub-list
            gives the indices of the strains that were
            delivered in that bite.

    """
    bites = [[i] for i in np.arange(n_bites)]
    for i in np.arange(n_bites, K):
        bites[np.random.choice(n_bites)].append(i)
    return bites


# --------------------------------------------------------------------------------
# Individual steps in Meiosis
#
# --------------------------------------------------------------------------------


class OocystMaker:
    def __init__(self):
        """
        Sample oocysts from a bespoke distribution modelled
        from empirical data on wild-caught A.g.

        """

        # Approximate distribution
        self.oocysts = np.arange(1, 11)
        self.probs = 0.4 * (0.5) ** (self.oocysts - 1)
        self.probs[2:6] = self.probs[2:6] + 0.025
        self.probs[9] = 1 - self.probs.sum()
        self.probs[9] += 1 - self.probs.sum()

    def sample_oocysts(self, n):

        return random.choices(self.oocysts, k=n, weights=self.probs)


def sample_zygotes(n_parents, n_oocysts, force_outbred=True):
    """
    Sample zygotes (pairs of parental strains) based on a `n_parents`
    and `n_oocysts`

    """

    # Sample pairs of parents randomly
    if not force_outbred:
        return [random.choices(range(n_parents), k=2) for _ in range(n_oocysts)]

    # Enforce one outbred zygote
    outbred_zygote = [random.sample(range(n_parents), k=2)]
    inbred_zygotes = [
        random.choices(range(n_parents), k=2) for _ in range(n_oocysts - 1)
    ]

    return outbred_zygote + inbred_zygotes


def sample_n_co(mean_n_co):
    """
    Sample an integer number of cross-over events for
    a chromosome, given a mean cross-over rate for
    the chromosome

    """

    return max(1, np.random.poisson(lam=mean_n_co))


def sample_co_brks(n_co, n_snps, co_probs):
    """
    Sample break points for an integer number of
    crossovers `n_co`, a number of SNPs `n_snps`,
    and a probability of recombination between the
    SNPs `co_probs`

    """

    co_brks = random.choices(range(n_snps), k=n_co, weights=co_probs)
    co_brks.sort()

    return co_brks


def sample_co_pairs(n_co):
    """
    Sample pairs of chromatids (indexed by 0,1) for
    each of `n_co` cross-over events

    """

    return [random.choices(range(2), k=2) for _ in range(n_co)]


def create_bivalent(n_snps, parents, zygote, pos_mask):
    """
    Create a bivalent, consistent of a pair of
    sister chromatids (i.e. four chromosomes total),
    from a set of parental strains `parents` and `zygote`
    indices

    """

    bivalent = np.zeros((2, n_snps, 2), dtype="int8")
    bivalent[:, :, 0] = np.copy(parents[zygote, :][:, pos_mask])
    bivalent[:, :, 1] = np.copy(parents[zygote, :][:, pos_mask])

    return bivalent


def resolve_cos(bivalent, co_brks, co_pairs):
    """
    Resolve cross-overs in a bivalent in place

    """

    for brk, pair in zip(co_brks, co_pairs):
        bivalent[[0, 1], :brk, pair] = bivalent[[1, 0], :brk, pair[::-1]]


def segregate(bivalent):
    """
    Segregate the chromosomes in a `bivalent`

    """

    shuffle = random.sample(range(4), k=4)

    return np.vstack([bivalent[:, :, 0], bivalent[:, :, 1]])[shuffle]


def get_recombinant_selections(zygotes, n_select):
    """
    Given a set of zygotes, create a list of all
    possible recombinant pairs of progeny, under the assumption
    that four progeny are produced per zygote

    """

    total_progeny = len(zygotes) * 4

    # Label each progeny 0 if inbred, 1 if outbred
    ii = np.concatenate([np.zeros(4) if a == b else np.ones(4) for a, b in zygotes])

    # All posssible progeny selections of size `n_select`
    all_selections = combinations(range(total_progeny), n_select)

    # At least one of the progeny must come from an outbred oocyst
    recombinant_selections = [
        selection
        for selection in all_selections
        if ii[list(selection)].sum()
        > 0  # at least one of the progeny must be from an outbred oocyst
    ]

    return recombinant_selections


def select_transmitted(zygotes, progeny, n_select, force_outbred=True):
    """
    Randomly select a set of progeny to be transmitted to the host,
    without replacement

    """

    # Enumerate progeny and ensure not trying to select too many
    n_progeny = progeny.shape[0]
    progeny_ixs = list(range(n_progeny))
    assert n_select <= n_progeny

    # If not outbred, just sample from progeny indexes
    if not force_outbred or n_select == 1:
        transmitted = random.sample(progeny_ixs, n_select)
        return progeny[transmitted]

    # If forcing outbred, sample from selections containing
    # at least one recombinant progeny
    recombinant_selections = get_recombinant_selections(zygotes, n_select)
    transmitted = list(random.choice(recombinant_selections))

    return progeny[transmitted]


# --------------------------------------------------------------------------------
# Meiosis function
#
# --------------------------------------------------------------------------------


def meiosis(parents, chromosomes, n_oocysts, n_select=None, force_outbred=True):
    """
    Here is our simple meiosis function

    """

    n_parents = parents.shape[0]
    zygotes = sample_zygotes(n_parents, n_oocysts, force_outbred=force_outbred)
    progeny = []
    for zygote in zygotes:
        zyg_progeny = []
        for chrom in chromosomes:
            n_co = sample_n_co(chrom.mean_n_co)
            co_brks = sample_co_brks(n_co, chrom.n_snps, chrom.co_probs)
            co_pairs = sample_co_pairs(n_co)
            bivalent = create_bivalent(chrom.n_snps, parents, zygote, chrom.pos_mask)
            resolve_cos(bivalent, co_brks, co_pairs)
            segregants = segregate(bivalent)
            zyg_progeny.append(segregants)
        progeny.append(np.hstack(zyg_progeny))
    progeny = np.vstack(progeny)

    if n_select is None:
        return progeny

    return select_transmitted(zygotes, progeny, n_select, force_outbred)


# --------------------------------------------------------------------------------
# Meiosis engine, for running meiosis multiple times while tracking IBD
#
# --------------------------------------------------------------------------------


class MeiosisEngine:
    """
    Engine for running meiosis multiple times, while tracking
    the IBD that is generated
    
    """
    def __init__(self, 
                 haplotypes: np.ndarray,
                 chroms: np.ndarray, 
                 pos: np.ndarray, 
                 haplotype_ixs: List[int] = None
    ) -> None:
        """
        Parse input data and prepare key attributes

        params
            haplotypes : ndarray [n_haplotypes, n_sites]
                Binary array [0, 1] representing alleles at each
                site for each haplotype.
            chroms : ndarray [n_sites]
                Array of chromosome numbers.
            pos : ndarray [n_sites]
                Array of variant site positions.
            haplotype_ixs: List[int] (optional)
                List of integers providing *indexes* to associate
                with each haplotype. These indexes are what define
                the haplotype within the infection, for purposes
                of pairwise IBD comparison. E.g. The IBD between
                strains (0, 1) is defined from the indices in this 
                list. We allow it to be optionally passed for
                cases where we are simulating meiosis between
                a subset of haplotypes in a larger infection.

        """
        # Input data
        self.n_haplotypes, self.n_sites = haplotypes.shape
        self.haplotypes = haplotypes
        self.haplotype_ixs = haplotype_ixs
        if self.haplotype_ixs is None:
            self.haplotype_ixs = list(range(self.n_haplotypes))

        # Necessary for simulating meiosis and/or detecting IBD
        self.chroms, self.pos = chroms, pos
        self.chromosomes = ChromosomeFactory(chroms, pos).create_chromosomes()
        self.genome_length = sum([c.l_kbp*10**3 for c in self.chromosomes])
        self.oocyst_maker = OocystMaker()
        
        # Results
        self.parents = self._create_parents()
        self.progeny = None
        self.n_progeny = 0
        self.ibd_segments = {}

    def _create_parents(self) -> np.ndarray:
        """
        Assigns integers to each parental haplotype in
        order to track the generation of IBD


        Overview
        * Each parental haplotype gets a unique integer (0, 1, 2...)
        * These parental haplotypes are subject to meiosis
        * After meiosis, they are used to detect IBD between pairs of strains
            * i.e. are two regions derived from the same parental haplotype?
        * They are also used to create the progeny haplotypes
            * The indices are used to recover the correct SNPs

        Notes
        * In theory, this is where we would want to allow for background IBD

        """

        parents = np.zeros((self.n_haplotypes, self.n_sites))
        for i in np.arange(self.n_haplotypes):
            parents[i] = i
        return parents

    def run(self, n_rounds: int = 1, n_select: int = None, force_outbred: bool = True):
        """
        Run meiosis for `n_rounds`, selecting `n_select` progeny
        after each round

        """

        oocysts = self.oocyst_maker.sample_oocysts(n_rounds)
        for n in oocysts:
            self.progeny = meiosis(
                parents=self.parents,
                chromosomes=self.chromosomes,
                n_oocysts=n,
                n_select=n_select,
                force_outbred=force_outbred,
            )
            self.parents = self.progeny
            self.n_progeny = self.progeny.shape[0]

    def get_progeny_haplotypes(self) -> np.ndarray:
        """
        Get haplotypes of the current meiotic progeny

        The progeny are created from the parental index arrays,
        they contain integers (0, 1, 2...) marking what parent they are
        derived from for each site. These are used to construct the haplotypes.

        """
        return self.haplotypes[self.progeny, range(self.progeny.shape[1])]
    
    def detect_generated_ibd(self) -> None:
        """
        Detect all of the IBD that was generated during round of meiosis and is now
        present in the progeny

        """
        # Create a boolean matrix that indicates presence of IBD between progeny
        ibd_matrix = create_pairwise_ibd_matrix(self.progeny)

        # Iterate over pairs of progeny, detect IBD segments
        for i, j in combinations(range(self.n_progeny), 2):
            s1 = self.haplotype_ixs[i]
            s2 = self.haplotype_ixs[j]
            self.ibd_segments[(s1, s2)] = get_ibd_segments(ibd_matrix[i, j], self.chroms, self.pos)

