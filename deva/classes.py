from typing import List, Optional, Iterable
from pandas import DataFrame
import numpy as np
from deva.constants import col_seps, col_outlier_suffix, col_not_outlier_suffix


def list_to_file(lis: Iterable, filename: str):
    """Takes an iterable and a file path and writes a value per line from the iterable into the new
    file.

    Args:
        lis: Iterable to write to file
        filename: Filename to write to.

    Returns:
        None

    """

    with open(filename, "w") as fh:
        for x in lis:
            fh.write("%s\n" % x)


class OutlierTable:
    """Output of calling outliers. """

    def __init__(
        self,
        df: DataFrame,
        updown: str,
        iqrs: Optional[float],
        samples: Optional[list],
        frac_table: Optional[DataFrame],
    ):
        """Instantiate an OutlierTable

        Args:
            df: DataFrame with outlier and non-outlier columns, and genes/sites as rows.
            updown: Whether the outliers are above or below the median. Options are "up" or "down"
            iqrs: The IQR threshold used to call outliers.
            samples: The samples included in the analysis to define median and IQR.
            frac_table: DataFrame with samples as columns and genes/sites as rows indicating
            what fraction of sites per sample were called as outliers. Useful for visualization.
        """

        self.df = df
        self.up_or_down = updown
        self.iqrs = iqrs
        self.samples = samples
        self.frac_table = frac_table

    def make_frac_table(self):
        """Constructs the fraction table from the outliers table

        Returns: A DataFrame with one column per sample, with the fraction of outliers per row
        per sample. This table is useful for visualization but not statistics.

        """

        df = self.df
        cols_outliers = [x + col_seps + col_outlier_suffix for x in self.samples]
        cols_notOutliers = [x + col_seps + col_not_outlier_suffix for x in self.samples]
        df = df.fillna(0)
        num_total_psites = df[cols_notOutliers].values + df[cols_outliers].values
        with np.errstate(divide="ignore", invalid="ignore"):
            frac_table = df[cols_outliers].values / num_total_psites

        self.frac_table = DataFrame(frac_table, index=df.index, columns=self.samples)

        return self.frac_table


class qValues:
    """Output from comparing groups using outliers. """

    def __init__(self, df: DataFrame, comps: list, frac_filter: Optional[float]):
        """Instantiates a qValues object.

        Args:
            df: DataFrame with genes/sites as rows and comparison_group as columns.
            comps: List of comparisons used to populate table.
            frac_filter: What fraction of samples in group of interest were required to have
            an outliers for any given row to be considered for analysis.
        """

        self.df = df
        self.comps = comps
        self.frac_filter = frac_filter

    def write_gene_lists(
        self,
        fdr_cut_off: float = 0.01,
        output_prefix: str = "outliers",
        comparisons: Optional[List] = None,
    ):
        """ Writes significant gene list files for every column in a qvalue table

        Args:
            fdr_cut_off: FDR threshold for significance
            output_prefix: Output prefix for files
            comparisons: which subset of qvalue columns to write gene lists for. Default will
            write for all columns

        Returns:

        """

        if comparisons is None:
            comparisons = self.df.columns
        else:
            comparisons = [
                col
                for col in self.df.columns
                if (col.startswith(tuple(comparisons))) and (col in self.comps)
            ]

        for comp in comparisons:
            sig_genes = list(self.df.loc[(self.df[comp] < fdr_cut_off), :].index)
            list_to_file(
                sig_genes, gene_list_file_name % (output_prefix, comp, fdr_cut_off)
            )
