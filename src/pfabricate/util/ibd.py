import numpy as np
import pandas as pd
from itertools import combinations
from dataclasses import dataclass


def create_pairwise_ibd_matrix(progeny, background_pairwise_ibd=None):
    """
    Create a pairwise IBD matrix from a set of meiotic progeny

    params
        progeny : ndarray, int, shape(n_strains, n_sites)
            Progeny from meiosis. Haplotypes should be indicated
            by integer values indicating from what parental
            strain each site was derived.

    """

    n_progeny, n_snps = progeny.shape
    pairwise_ibd = np.ones((n_progeny, n_progeny, n_snps), dtype="bool")
    for i, j in combinations(range(n_progeny), 2):
        ibd = progeny[i] == progeny[j]
        pairwise_ibd[i, j] = ibd
        pairwise_ibd[j, i] = ibd

    if background_pairwise_ibd is not None:
        assert background_pairwise_ibd.dtype == bool
        assert pairwise_ibd.shape == background_pairwise_ibd.shape
        pairwise_ibd = pairwise_ibd | background_pairwise_ibd

    return pairwise_ibd


@dataclass
class IBDSegment:
    chrom: str
    start: int
    end: int
    length: int = None

    def __post_init__(self):
        self.length = self.end - self.start


def get_ibd_segment_dataframe(ibd, chroms, pos):
    """
    Return a dataframe of IBD segments from an `ibd`,
    `chrom` and `pos` arrays

    TODO
    - Some questions about how to assign IBD to the SNP
    positions -- probably want to use midpoints between
    SNPs
    - This function (minus the dataclass) could be
    quite easily jit'd

    """

    ibd_segs = []

    cc = chroms[0]
    start = pos[0]
    end = pos[0]

    for state, c, p in zip(ibd, chroms, pos):
        if state and c == cc:
            end = p
            continue

        if end - start > 0:
            ibd_segs.append(IBDSegment(chrom=cc, start=start, end=end))

        cc = c
        start = p
        end = p

    if state:
        ibd_segs.append(IBDSegment(chrom=cc, start=start, end=end))

    return pd.DataFrame(ibd_segs)


def calc_n50(segments):
    """
    Calculate the N50 of a set of `segments`

    params
        segments : ndarray, shape(n_segments,)
            Numpy array of segment lengths. Could be
            IBD segments, contigs, &c.
    returns
        _ : float
            A length X of the smallest segment in the set where
            at least 50% of all bases are in segments smaller than X.
            Alternatively, the median segment size weighted by size.

    """

    if not isinstance(segments, np.ndarray):
        segments = np.array(segments)
    assert (segments > 0).all()

    if len(segments) == 0:
        return np.nan
    segments.sort()
    total = segments.sum()
    cuml_frac = (segments / total).sum()
    ix = np.argmax(cuml_frac > 0.5)

    return segments[ix]
