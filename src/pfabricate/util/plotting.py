import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np


class WSAFPlotter:

    EDGE_BUF = 5 * 10**5
    CHROM_DIFF = 10**4
    CHROM_COLS = sns.color_palette("Paired", 10) + sns.color_palette("Paired", 4)

    def __init__(self):
        pass

    def _set_chromosome_colors(self, chroms):
        self.uniq_chroms = np.unique(chroms)
        self.col_dt = dict(zip(self.uniq_chroms, self.CHROM_COLS))
        self.chrom_cols = [self.col_dt[c] for c in chroms]

    def set_snp_positions(self, chroms, pos):
        """Set positions of SNPs in plot based on chromosome and position array"""

        # Store as attributes
        self.chroms = chroms
        self.pos = pos

        # Colors
        self._set_chromosome_colors(chroms)

        # Compute x positions of SNPs in plot
        diffs = pos[1:] - pos[:-1]
        p = 0
        self.plot_pos = [p]
        cc = chroms[0]
        for d, c in zip(diffs, chroms[1:]):
            if c == cc:
                p += d
            else:
                p += self.CHROM_DIFF
                cc = c
            self.plot_pos.append(p)
        self.plot_pos = np.array(self.plot_pos)

    def plot(self, title, wsaf, output_path=None):
        """Plot WSAF along the genome for a sample"""

        fig, ax = plt.subplots(1, 1, figsize=(10, 4))

        ax.scatter(x=self.plot_pos, y=wsaf, c=self.chrom_cols, s=2)

        ax.set_title(title, loc="left")
        ax.set_ylabel("WSAF")
        ax.set_xlabel("Genome Position\n(approx)")

        ax.yaxis.set_major_locator(plt.MultipleLocator(0.1))
        ax.yaxis.set_minor_locator(plt.MultipleLocator(0.05))

        ax.set_axisbelow(True)
        ax.grid(ls="dotted", axis="y", alpha=0.5)

        ax.set_xlim(
            (self.plot_pos.min() - self.EDGE_BUF, self.plot_pos.max() + self.EDGE_BUF)
        )

        if output_path is not None:
            fig.savefig(output_path, dpi=200, bbox_inches="tight", pad_inches=0.5)
            plt.close()
