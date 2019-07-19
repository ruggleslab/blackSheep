from typing import List, Optional
from pandas import DataFrame


def listToFile(lis, file_name):
    with open(file_name, "w") as fh:
        for x in lis:
            fh.write("%s\n" % x)


class OutlierTable:
    def __init__(
        self,
        df: DataFrame,
        updown: str,
        iqrs: float,
        samples: List,
        frac_table: Optional[DataFrame],
    ):
        self.df = df
        self.up_or_down = updown
        self.iqrs = iqrs
        self.samples = samples
        self.frac_table = frac_table


class qValues:
    def __init__(self, df: DataFrame, comps: List, frac_filter: Optional[float]):
        self.df = df
        self.comps = comps
        self.frac_filter = frac_filter

    def write_gene_lists(
        self,
        fdr_cut_off: float = 0.01,
        prefix: str = "outliers",
        comparisons: Optional[List] = None,
    ):

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
            listToFile(
                sig_genes, "%s.%s.sig_genes.fdr%s.txt" % (prefix, comp, fdr_cut_off)
            )
