import random
import numpy as np
from itertools import combinations


# TODO
# - These functions are definitely subtle enough to require good documentation
# - Bonus would be to do a bit of typing stuff
# - And they are nicely broken up so should be quite easy to write some tests


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
