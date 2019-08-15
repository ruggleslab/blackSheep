<!-- blacksheep documentation master file, created by
sphinx-quickstart on Thu Aug 15 14:49:18 2019.
You can adapt this file completely to your liking, but it should at least
contain the root `toctree` directive. -->
# Welcome to blacksheep’s documentation!

# DEVA (Differential Extreme Value Analysis)


#### deva.make_outliers_table(df: pandas.core.frame.DataFrame, iqrs: float = 1.5, up_or_down: str = 'up', aggregate: bool = True, save_outlier_table: bool = False, save_frac_table: bool = False, output_prefix: str = 'outliers', ind_sep: str = '-')
Converts a DataFrame of values into an OutlierTable object, which includes a DataFrame
of outlier and non-outlier count values.


* **Parameters**

    * **df** – Input DataFrame with samples as columns and sites/genes as columns.

    * **iqrs** – The number of inter-quartile ranges (IQRs) above or below the median to consider a         value as an outlier.

    * **up_or_down** – Whether to call up or down outliers. Up is above the median; down         is below the median. Options “up” or “down”.

    * **aggregate** – Whether to sum outliers across a grouping (e.g. gene-level) than individual         sites. For instance if columns indicate phosphosites on proteins, with the format         “RAG2-S365”, output will show counts of outliers per protein (e.g. RAG2) rather than on         individual sites (e.g. RAG2-S365).

    * **save_outlier_table** – Whether to write a file with the outlier count table.

    * **save_frac_table** – Whether to write a file with the outlier fraction table.

    * **output_prefix** – If files are written, a prefix for the files.

    * **ind_sep** – The separator used in sites, for instance, to separate a gene and site.         If just using genes (i.e. no separator), or not aggregating this parameter has no effect.


Returns: outliers

    Returns an OutlierTable object, with outlier and non-outlier counts and metadata
    about how the outliers were called.


#### deva.compare_groups_outliers(outliers: deva.classes.OutlierTable, annotations: pandas.core.frame.DataFrame, frac_filter: Optional[float] = 0.3, save_qvalues: bool = False, output_prefix: str = 'outliers', save_comparison_summaries: bool = False)
Takes an OutlierTable object and a sample annotation DataFrame and performs comparisons for
any column in annotations with exactly 2 groups. For each group identified in the annotations
DataFrame, this function will calculate the q-values of enrichment of outliers for each row in
each group.


* **Parameters**

    * **outliers** – An OutlierTable, with a DataFrame of outlier and non-outlier counts,         as well as parameters for how outliers were calculated.

    * **annotations** – A DataFrame with samples as rows and annotations as columns. Each         column must contain exactly 2 different categories, not counting missing values. Columns         without 2 options will be ignored.

    * **frac_filter** – The fraction of samples in the group of interest that must         have an outlier value to be considered in the comparison. Float between 0 and 1 or None.

    * **save_qvalues** – Whether to write a file with a table of qvalues.

    * **output_prefix** – If files are written, a prefix for the files.

    * **save_comparison_summaries** – Whether to write a file for each annotation column with the         counts in the fisher table, pvalues and q values per row.


Returns: qvals

    A qValues object, which includes a DataFrame of q-values for each comparison,         as well as some metadata about how the comparisons were performed.


#### deva.run_outliers(df: pandas.core.frame.DataFrame, annotations: pandas.core.frame.DataFrame, iqrs: float = 1.5, up_or_down: str = 'up', aggregate: bool = True, save_outlier_table: bool = False, save_frac_table: bool = False, frac_filter: Optional[float] = 0.3, save_qvalues: bool = False, output_prefix: str = 'outliers', ind_sep: str = '-', save_comparison_summaries: bool = False)
Takes a DataFrame of values and returns OutlierTable and qValues objects. This command runs
the whole outliers pipeline. The DataFrame in the OutlierTable object can be used to run more
comparisons in future. The qValues object can be used for visualization, or writing
significant gene lists.


* **Parameters**

    * **df** – Input DataFrame with samples as columns and sites/genes as rows.

    * **annotations** – A DataFrame with samples as rows and annotations as columns. Each         column must contain exactly 2 different values, not counting missing         values. Other columns will be ignored.

    * **iqrs** – The number of interquartile ranges (IQRs) above or below the median to consider a         value as an outlier.

    * **up_or_down** – Whether to call up or down outliers. Up is above the median; down         is below the median. Options “up” or “down”.

    * **aggregate** – Whether to sum outliers across a grouping (e.g. gene-level) than individual         sites. For instance if columns indicate phosphosites on proteins, with the format         “RAG2-S365”, output will show counts of outliers per protein (e.g. RAG2) rather than on         individual sites (e.g. RAG2-S365).

    * **save_outlier_table** – Whether to write a file with the outlier count table.

    * **save_frac_table** – Whether to write a file of the fraction of outliers.

    * **frac_filter** – The fraction of samples in the group of interest that must         have an outlier value to be considered in the comparison. Float between 0 and 1 or None.

    * **save_qvalues** – Whether to write a file of qvalues.

    * **output_prefix** – If files are written, a prefix for the files.

    * **ind_sep** – The separator used in the columns, for instance, to separate a gene and site.         If just using genes (i.e. no separator), or not aggregating this parameter         has no effect.

    * **save_comparison_summaries** – Whether to write a table for each comparison with the         counts in the fisher table, pvalues and qvalues per row.


Returns: outliers, qvals

    Returns an OutlierTable object and qValues object.


#### deva.plot_heatmap(annotations: pandas.core.frame.DataFrame, qvals: pandas.core.frame.DataFrame, col_of_interest: str, vis_table: pandas.core.frame.DataFrame, fdr: float = 0.05, red_or_blue: str = 'red', output_prefix: str = 'outliers', colors: Optional[str] = None, savefig: bool = False)
Plots a heatmap of significantly enriched values for a given comparison.


* **Parameters**

    * **annotations** – Annotations DataFrame, samples as rows, annotations as columns

    * **qvals** – qvalues DataFrame with genes/sites as rows and comparisons as columns

    * **col_of_interest** – Which column from qvalues should be used to find signficant genes

    * **vis_table** – Table to be visualized in heatmap. Index values should correspond to the         annotation df index, column names should correspond to qvals df index

    * **fdr** – FDR threshold to for significance

    * **red_or_blue** – Whether heatmap should be in red or blue color scale

    * **output_prefix** – If saving files, output prefix

    * **colors** – File to find color map for annotation header

    * **savefig** – Whether to save the plot to a pdf


Returns: [annot_ax, vals_ax, cbar_ax, leg_ax]

    List of matplotlib axs, can be further customized before saving. In order the axes         contain: annotation header, the heatmap, the color bar, and the legend.


#### deva.binarize_annotations(df: pandas.core.frame.DataFrame)
Takes an annotation DataFrame, checks each column for the number of possible values,
and adjusts based on that. If the column has 0 or 1 options, it is dropped. Cols with 2
possible values are retained as-is. Cols with more than 2 values are expanded. For each
value in that column, a new column is created with val and not_val options.


* **Parameters**

    **df** – Annotations DataFrame.


Returns: new_df

    Refactored annotations DataFrame.


#### deva.read_in_values(path: str)
Figures out sep and parsing file into dataframe.


* **Parameters**

    **path** – File path


Returns: df

    DataFrame from table in file


#### deva.read_in_outliers(path: str, updown: str, iqrs: float)
Parses a file into an OutlierTable object.


* **Parameters**

    * **path** – File path

    * **updown** – Whether the outliers represent up or down outliers

    * **iqrs** – How many IQRs were used to define an outlier


Returns: outliers

    OutlierTable object

# Indices and tables

* Index

* Module Index

* Search Page
