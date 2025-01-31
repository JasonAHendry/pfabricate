import numpy as np
import pandas as pd
from typing import List, Tuple, Dict
from itertools import combinations
from pfabricate.mixed.chromosomes import ChromosomeFactory
from dataclasses import fields
from pfabricate.util.ibd import (
     IBDSegment,
     convert_ibd_segments_to_dataframe,
     IBDSummaryStatistics
)
from .meiosis import (
    partition_strains,
    MeiosisEngine
)


class MalariaInfectionMaker:
    """
    Create a malaria infection from a set of input_haplotypes

    """
    
    # Parameters
    MAX_K = 4
    MAX_M = 4
    MIN_SITES = 100

    def __init__(self,
                 input_haplotypes: np.ndarray,
                 chroms: np.ndarray,
                 pos: np.ndarray,
                 B: int,
                 M: int,
                 props: np.ndarray, # should I allow this to be optional?
                 input_names: List[str] = None
                 ) -> None:
        """
        Initialise key variables and sanity check inputs, with the
        COI (K) being set by the number of input_haplotypes provided

        """
    
        # Parse inputs
        self.input_haplotypes = input_haplotypes
        self.chroms, self.pos = chroms, pos
        self.K, self.n_sites = input_haplotypes.shape
        self.B = B  # Number of infectious bites
        self.M = M  # Number of rounds of meiosis for haplotypes in bites
        self.props = np.copy(props) # avoid sorting original
        self.props.sort()
        self.Keff = 1 / (self.props**2).sum()
        self.input_names = input_names

        # This is dumb and repetitive /w MeiosisEngine, but needed for IBD calculations
        self.chromosomes = ChromosomeFactory(chroms, pos).create_chromosomes()
        self.genome_length = sum([c.l_kbp*10**3 for c in self.chromosomes])

        # Checks
        self._sanity_checks()

        # Infection contents
        self.sample_id = None
        self.infection_haplotypes = None
        self.ibd_segments = {(i, j): [] for i, j in combinations(range(self.K), 2)}

    def _sanity_checks(self):
        """
        Perform a set of (non-exhaustive) sanity checks on the input data
        """
        if self.K > self.MAX_K:
            raise ValueError(f"Number of haplotypes exceeds maximum of {self.MAX_K}.")
        if self.n_sites <= self.MIN_SITES:
            raise ValueError(f"Provided haplotypes have only {self.n_sites} segregating sites,"
                             "minimum is set to {self.MIN_SITES}.")
        if self.B > self.K:
            raise ValueError("Cannot have more infectious bites than available haplotypes.")
        if self.B < 1:
            raise ValueError("Must have at least one infectious bite to create infection.")
        if not 0 <= self.MAX_M <= 4:
            raise ValueError(f"Number of meioses exceeds maximum of {self.MAX_M}.")
        if self.props.shape[0] != self.K:
            raise ValueError(f"Number of proportions {self.props.shape} must match COI={self.K}.")
        if len(self.input_names) != self.K:
            raise ValueError(f"Number of input strain names {len(self.input_names)} must match COI={self.K}.")

    def assign_sample_id(self, id_prefix: str, index: int) -> None:
        """Assign an ID to the infection"""
        self.sample_id = f"{id_prefix}{index:04d}-K{self.K:02d}-B{self.B:02d}-M{self.M:02d}"

    def _update_ibd_segments(self, new_ibd: Dict[Tuple[int, int],List[IBDSegment]]):
        """
        Update the IBD segments dictionary with an additional set of segments

        TODO
        * Should implement a check to merge IBD segments that are overlapping
        
        """
        
        for strain_pair, ibd_segments in new_ibd.items():
            self.ibd_segments[strain_pair].extend(ibd_segments)

    def generate(self) -> None:
        """
        Generate the infection by simulating individual bites and rounds of meiosis

        The key result of this function is the generation of `self.infection_haplotypes`,
        which are binary numpy arrays representing the haplotype of each strain in the 
        infection; and populating `self.ibd_segments`, which contains a list of IBD
        segments for each pair of strains in the infection

        """

        if self.B == self.K: # Each strain was transferred separately, no meiosis.
            self.infection_haplotypes = self.input_haplotypes
            return
        
        # Partition strains into infectious bites
        bites = partition_strains(self.K, self.B) # e.g. [[0], [1,2]]

        # Simulate bites
        transmitted = [] # Store haplotypes that get transmitted, across the bites
        for bite in bites:
            # If one strain in the bite, no meiosis
            if len(bite) == 1:
                transmitted.append(self.input_haplotypes[bite])
                continue
            
            # Simulate meiosis
            meiosis_engine = MeiosisEngine(
                    haplotypes=self.input_haplotypes[bite],
                    chroms=self.chroms,
                    pos=self.pos,
                    haplotype_ixs=bite
            )
            meiosis_engine.run(
                n_rounds=self.M, 
                n_select=len(bite), 
                force_outbred=True  # same in as out
            )
            meiosis_engine.detect_generated_ibd()

            # Store new IBD segments and haplotypes
            self._update_ibd_segments(meiosis_engine.ibd_segments)
            transmitted.append(meiosis_engine.get_progeny_haplotypes())
        
        # Concatenate final infection haplotypes across all bites
        self.infection_haplotypes = np.vstack(transmitted)

    def get_infection_haplotypes(self) -> np.ndarray:
        return self.infection_haplotypes
    
    def get_ibd_segments(self) -> pd.DataFrame:
        """
        Combine information about IBD segments across all
        pairwise strain comparisons into a dataframe
        """
        dfs = [pd.DataFrame([], columns=[f.name for f in fields(IBDSegment)])] # handles K=1
        for strains, segments in self.ibd_segments.items():
            df = convert_ibd_segments_to_dataframe(segments)

            # if df is empty, below just adds column names; no values
            df.insert(0, "sample_id", self.sample_id) 
            df.insert(1, "strain1", strains[0])
            df.insert(2, "strain2", strains[1])

            dfs.append(df)
        return pd.concat(dfs)

    def get_ibd_pairwise(self) -> pd.DataFrame:
        dts = []
        for (s1, s2), segments in self.ibd_segments.items():
            dt = {
                "sample_id": self.sample_id,
                "strain1": s1,
                "strain2": s2,
                "prop_strain1": self.props[s1],
                "prop_strain2": self.props[s2]
            }
            dt.update(
                IBDSummaryStatistics.from_segments(segments, self.genome_length).as_dict()
            )
            dts.append(dt)
        return pd.DataFrame(dts)
    
    def get_sample_summary(self, **annot_kwargs) -> Dict:
        """Return a sample-level summary of the infection"""
        
        # Prepare sample names and props into string
        # TODO: Can be done in init
        samps_str = None
        if self.input_names is not None:
            samps_str = ";".join(self.input_names)
        props_str = ";".join([str(p) for p in self.props])

        # Get all the other attributes I want to add
        samp_dt = {
            k: getattr(self, k) for k in ["sample_id", "K", "Keff", "B", "M"]
        }
        samp_dt["props"] = props_str
        samp_dt["samps"] = samps_str

        # Add optional annotations
        for annot, value in annot_kwargs.items():
            samp_dt[annot] = value

        # Compute IBD summary statistics
        all_segements = []
        for strains, segments in self.ibd_segments.items():
            all_segements.extend(segments)
        ibd_dt = IBDSummaryStatistics.from_segments(
            all_segements,
            self.genome_length*len(self.ibd_segments)
        ).as_dict()
        samp_dt.update(ibd_dt)

        return samp_dt



        

