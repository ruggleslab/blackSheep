from typing import List, Optional
from pandas import DataFrame
import numpy as np
from .constants import *
import blacksheep as bsh


class OutlierTable:
    def __init__(
        self,
        df: DataFrame,
        updown: str,
        iqrs: Optional[float],
        samples: list,
        frac_table: Optional[DataFrame],
    ):
        self.df = df
        self.up_or_down = updown
        self.iqrs = iqrs
        self.samples = samples
        self.frac_table = frac_table

    def make_frac_table(self):
        """
        Constructs the fraction table from the outliers table
        :return: A table with one column per sample with the fraction of sites in each row
        that are outliers per each sample. This table is useful for visualization but not statistics.
        """
        df = self.df
        cols_outliers = [x + col_seps + col_outlier_suffix for x in self.samples]
        cols_notOutliers = [x + col_seps + col_not_outlier_suffix for x in self.samples]
        df = df.fillna(0)
        num_total_psites = df[cols_notOutliers].values + df[cols_outliers].values
        with np.errstate(divide="ignore", invalid="ignore"):
            frac_table = df[cols_outliers].values / num_total_psites

        self.frac_table = DataFrame(
            frac_table, index=df.index, columns=self.samples
        )

        return self.frac_table


class qValues:
    def __init__(self, df: DataFrame, comps: list, frac_filter: Optional[float]):
        self.df = df
        self.comps = comps
        self.frac_filter = frac_filter

    def write_gene_lists(
        self,
        fdr_cut_off: float = 0.01,
        output_prefix: str = "outliers",
        comparisons: Optional[List] = None,
    ):
        """
        Writes signficant gene list files for every column of a qvalue table
        :param fdr_cut_off: FDR threshold for signficance
        :param output_prefix: Output prefix for files
        :param comparisons: which subset of qvalue columns to write gene lists for. Default will
        write for all columns
        :return: None
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
            bsh.parsers.list_to_file(sig_genes, gene_list_file_name % (output_prefix, comp,
                                                                   fdr_cut_off))
