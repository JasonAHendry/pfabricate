import numpy as np
from typing import List
from scipy.stats import nbinom, betabinom


def convert_ad_to_haplotypes(ad):
    """
    Convert allleic depth array `ad` into clonal
    haplotypes

    """
    if ad.shape[2] > 2:
        ad = ad[:, :, :2]  # biallelic

    haplotypes = np.argmax(ad, axis=2)
    assert ((haplotypes == 0) | (haplotypes == 1)).all()

    return haplotypes


def call_genotypes_perfectly(infection_haplotypes: np.ndarray) -> List[str]:
    """
    Given a set of haplotypes in an infection, perfectly call the overall
    genotype
    
    """

    k = infection_haplotypes.shape[0]
    alt_counts = infection_haplotypes.sum(0)

    assert k <= 20
    assert alt_counts.max() <= k
    assert alt_counts.min() >= 0

    gts = []
    for n_alt in alt_counts:
        if n_alt == 0:
            gts.append("0/0")
        elif n_alt == k:
            gts.append("1/1")
        else:
            gts.append("0/1")

    return gts


def simulate_read_data(
    haplotypes, proportions, depth_mean, depth_shape, alt_shape, e_0, e_1
):
    """
    Simulate read data haplotypes and model parameters

    params
        haplotypes : ndarray, int, shape (K, n_snps)
            Allelic states for each of the `K`
            strains across all SNPs.
        proportions : ndarray, float, shape (K)
            Proportions for all strains.
        depth_mean : float
            Specifies mean depth in negative binommial
            model.
        depth_shape : float
            Specifies dispersion in the negative
            binomial distribution that is used to
            model read depth.
        alt_shape : float
            Specifies the dispersion in the
            beta-binomial distribution used to
            compute the probability of different
            WSAF, given the expectation under
            the proportions and error.
        e_0 : float
            ref->alt error rate.
        e_1 : float
            alt->ref error rate.

    returns
        output_dt : dict
            "depth" : ndarray, int, shape (n_snps)
                Simulated read depth across SNPs. Equal
                to ref + alt.
            "alt" : ndarray, int, shape (n_snps)
                Simulated alternative allele counts across SNPs.
            "ref" : ndarray, int, shape (n_snps)
                Simulated reference allele counts across SNPs.
            "gts" : List[str, shape(n_snps)
                Genotype calls for infection haplotypes

    """

    np.testing.assert_almost_equal(proportions.sum(), 1.0)

    n_snps = haplotypes.shape[1]

    # READ DEPTH
    read_depths = nbinom.rvs(
        n=depth_shape, p=depth_shape / (depth_shape + depth_mean), size=n_snps
    )

    # EXPECTED WSAF (including error)
    expected_wsaf = proportions.dot(haplotypes)
    expected_wsaf += e_0 * (1 - expected_wsaf) - e_1 * expected_wsaf

    # ALTERNATIVE COUNTS
    alphas = alt_shape * expected_wsaf
    betas = alt_shape * (1 - expected_wsaf)

    # NUMBER OF ALTERNATIVE READS
    sim_alt = betabinom.rvs(n=read_depths, a=alphas, b=betas)
    sim_ref = read_depths - sim_alt

    # Compile output
    output_dt = {"depth": read_depths, 
                 "alt": sim_alt, 
                 "ref": sim_ref,
                 "gts": call_genotypes_perfectly(haplotypes)
    }

    return output_dt
