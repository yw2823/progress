
import os
import pyBigWig as bw
from pyranges import PyRanges as pr, read_bed


#compute average signal over a region from bigwig files
class bigwig:
    """
    """

    def __init__(self, bigwig_file_path: str, regions: pr):
        self.bigwig_file_path = bigwig_file_path
        self.regions = regions
        self.bw = self._open_bigwig()
        self.average_values = self._compute_average_values()

    def _open_bigwig(self):
        """Opens the BigWig file and returns a handle."""
        if not os.path.exists(self.bigwig_file_path):
            raise FileNotFoundError(f"BigWig file not found: {self.bigwig_file_path}")
        return bw.open(self.bigwig_file_path)

    def _compute_average_values(self):
        """
        Computes the average signal values over each region.
        Returns a list of average values.
        """

        def compute_average(region):
            chrom, start, end = region.Chromosome[0], region.Start[0], region.End[0]
            value = self.bw.stats(chrom, start, end)[0]
            return 0 if value is None else value

        return list(map(compute_average, self.regions.df.itertuples(index=False)))

    def close(self):
        """Closes the BigWig file."""
        self.bw.close()


# %%

if __name__ == "__main__":
    #genomic_regions = read_bed("")
    #test = bigwig("", atac_regions)
    
    #print(test.average_values)

    # Close BigWig file
    #test.close()