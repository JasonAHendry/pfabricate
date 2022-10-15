import numpy as np
from dataclasses import dataclass
from allel import read_vcf


@dataclass
class Chromosome:
    name: str
    n_snps: int
    l_kbp: int
    mean_n_co: float
    co_probs: np.ndarray
    pos_mask: np.ndarray


class ChromosomeFactory:

    kbp_per_cM = 13.5

    def __init__(self, chroms, pos):
        """
        Create chromosome objects from genome-wide position and chromosome information

        """
        self.chroms = chroms
        self.pos = pos
        self.chrom_set = set(chroms)

    @classmethod
    def from_vcf(cls, vcf_path):
        vcf = read_vcf(input=vcf_path)
        return cls(chroms=vcf["variants/CHROMS"], pos=vcf["variants/POS"])

    def create_chromosomes(self):
        chromosomes = []
        for c in self.chrom_set:

            # Extract positions belonging to chromosome with a mask
            chrom_mask = self.chroms == c
            n_snps = chrom_mask.sum()

            # Compute mean number of cross-overs for this chromosome
            l_kbp = (self.pos[chrom_mask].max() - self.pos[chrom_mask].min()) / 10**3
            Morgans = l_kbp / (self.kbp_per_cM * 100)
            mean_n_co = 2 * Morgans

            # Compute the probability that a cross-over event occurs between
            # pair of SNPs, assuming a uniform recombination rate
            frac_chrom_bw_snps = np.append(np.diff(self.pos[chrom_mask]), 0) / (
                l_kbp * 10**3
            )

            chromosomes.append(
                Chromosome(
                    name=c,
                    n_snps=n_snps,
                    l_kbp=l_kbp,
                    mean_n_co=mean_n_co,
                    co_probs=frac_chrom_bw_snps,
                    pos_mask=chrom_mask,
                )
            )

        return chromosomes
