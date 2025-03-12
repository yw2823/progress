from os import path
import numpy as np
import pandas as pd
from pyranges import read_bed
from pyranges import PyRanges as pr
from scipy.sparse import csr_matrix, save_npz


def read_bed_4col(bedfile):
    bed = pr(read_bed(bedfile, as_df=True).rename({'Name':'Score'}, axis=1).query('Score>0'), int64=True)
    return bed

class atac_seq(object):
    def __init__(self):
        self.peak_bed = self.read_atac()

        self.accessibility = self.peak_bed.as_df().iloc[:, 3].values
        self.promoter_atac = self.get_promoter_atac()


    #bed file remove 0 intervals
    def extract_peak(self, slop=100, target_length=2000):
        """Extracts sequences for ATAC peaks."""
        genome = Genome(assembly=self.assembly, fasta_file=f'{self.assembly}.fa')
        bed = GenomicRegionCollection(genome, self.peak_bed.df)
        seq = bed.collect_sequence(upstream=slop, downstream=slop, target_length=target_length)
        return seq



    def peak_motif(self, with_cutoff=False, save_as_bed=False):
        """
        Reads a peak motif bed file and returns a dataframe with the motifs in each peak.
        
        Args:
        - with_cutoff (bool): If True, applies a threshold cutoff for motif scores.
        - save_as_bed (bool): If True, saves the DataFrame as a BED file.

        Returns:
        - DataFrame containing peaks with motif scores.
        """

        # Define file paths
        bed_output_file = f"{self.sample}.atac.motif.output.bed"

        # Read peak regions
        peaks = self.peak_bed.as_df().reset_index()

        # Read peak motif BED file
        peak_motif = pd.read_csv(
            f"{self.sample}.peak_motif.bed", sep='\t', header=None,
            names=['chr', 'start', 'end', 'motif', 'score']
        ).drop_duplicates()

        # Pivot the table so that each motif gets its own column
        peak_motif = peak_motif.pivot(index=['chr', 'start', 'end'], 
                                        columns='motif', 
                                        values='score').reset_index().fillna(0)

        # Merge with peak regions
        peak_motif = pd.merge(
            peaks, peak_motif,
            left_on=['Chromosome', 'Start', 'End'], 
            right_on=['chr', 'start', 'end'], 
            how='left'
        )

        # Drop duplicate or redundant columns
        peak_motif.drop(['chr', 'start', 'end'], axis=1, inplace=True)
        peak_motif.set_index('index', inplace=True)

        # Apply cutoff threshold if required
        if with_cutoff:
            cutoff_values = [self.get_motif_cutoff(motif, peak_motif) for motif in peak_motif.columns[3:]]
        else:
            cutoff_values = [0 for _ in peak_motif.columns[3:]]

        # Apply thresholding
        peak_motif.iloc[:, 3:] = (peak_motif.iloc[:, 3:].values > np.array(cutoff_values)) * peak_motif.iloc[:, 3:].values

        # Optionally save as a BED file
        if save_as_bed:
            bed_columns = ['Chromosome', 'Start', 'End'] + list(peak_motif.columns[3:])
            peak_motif[bed_columns].to_csv(bed_output_file, sep='\t', index=False, header=False)
            print(f"Saved peak motif data as BED: {bed_output_file}")

        return peak_motif




from pysam import tabix_compress, tabix_index

def df_to_bedgraph(df, output_file):
    '''
    Four columns in df: chr, start, end, score
    '''
    df.to_csv(output_file, sep='\t', header=False, index=False)

    zip_file = f"{output_file}.gz"
    pysam.tabix_compress(output_file, zip_file)
    pysam.tabix_index(zip_file, preset="bed")
