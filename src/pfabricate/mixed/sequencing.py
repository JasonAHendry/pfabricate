import numpy as np
from scipy.stats import nbinom, betabinom


def simulate_read_data(
    haplotypes, proportions, mean_depth, depth_shape, alt_shape, e_0, e_1
):
    """
    Simulate read data haplotypes and model parameters

    params
        haplotypes : ndarray, int, shape (K, n_snps)
            Allelic states for each of the `K`
            strains across all SNPs.
        proportions : ndarray, float, shape (K)
            Proportions for all strains.
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

    """

    np.testing.assert_almost_equal(proportions.sum(), 1.0)

    n_snps = haplotypes.shape[1]

    # READ DEPTH
    read_depths = nbinom.rvs(
        n=depth_shape, p=depth_shape / (depth_shape + mean_depth), size=n_snps
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
    output_dt = {"depth": read_depths, "alt": sim_alt, "ref": sim_ref}

    return output_dt
