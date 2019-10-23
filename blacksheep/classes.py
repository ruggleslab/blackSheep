from typing import List, Optional, Iterable
import logging
import pandas as pd
from pandas import DataFrame
import numpy as np
from blacksheep._constants import col_seps, col_outlier_suffix, col_not_outlier_suffix, gene_list_file_name


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


def make_frac_table(df, samples):
    """Constructs the fraction table from the outliers table

    Returns: A DataFrame with one column per sample, with the fraction of outliers per row
    per sample. This table is useful for visualization but not statistics.

    """
    df = df.copy()
    cols_outliers = [x + col_seps + col_outlier_suffix for x in samples]
    cols_notOutliers = [x + col_seps + col_not_outlier_suffix for x in samples]
    df = df.fillna(0)
    num_total_psites = df[cols_notOutliers].values + df[cols_outliers].values
    with np.errstate(divide="ignore", invalid="ignore"):
        frac_table = df[cols_outliers].values / num_total_psites

    return DataFrame(frac_table, index=df.index, columns=samples)

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
        if frac_table is not None:
            self.frac_table = frac_table
        else:
            self.frac_table = make_frac_table(df, samples)



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

        Returns: None

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


    def make_signed_logqs(self) -> DataFrame:
        """Create a DataFrame with signed log10 qvalues for each comparison. E.g. group1 qvalues
        will be positive, and group 2 qvalues will be negative. Assignment of positive group is
        based on order in qvalues, could be helpful to negate some columns in output depending on
        group of interest.

        Returns: DataFrame with signed qvalues.

        """
        if not (self.comps is None):
            self.comps = [i.split('_', 1)[1].rsplit('_', 1)[0] for i in self.df.columns]
            self.comps = sorted(list(set(self.comps)))

        signed_qs = pd.DataFrame()
        for comp in self.comps:
            cols = [
                col for col in self.df.columns if col.split('_', 1)[1].rsplit('_', 1)[0] == comp
            ]

            if len(cols) > 2 or len(cols) == 0:
                logging.warning("Excluding %s, %s columns are associated with %s, need 1 or 2 "
                                "columns. Annotation value probably has an _ in it. "
                                %(comp, len(cols), comp))
                continue

            if len(cols) == 1:
                temp = self.df[cols[0]]
            elif len(cols) == 2:
                temp = pd.DataFrame(
                    -np.log10(self.df[cols[1]]).subtract(np.log10(self.df[cols[0]]), fill_value=0),
                    columns=['%s_%s'%(comp, cols[1].rsplit('_', 1)[1])]
                )
            signed_qs = pd.concat([signed_qs, temp], join='outer', axis=1, sort=False)

        return signed_qs
