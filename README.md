# BlackSheep
##### A tool for differential extreme-value analysis

### Installation
With pip
```bash
pip install blacksheep-outliers
```
With conda
```bash
conda install -c bioconda blacksheep-outliers
```

### Requirements
pandas  
numpy  
scipy  
scikit-learn  


### Usage
##### In python
```python
import deva

# Read in data
values_file = '' #insert values file here
annotations_file = '' #insert annotations file here
values = deva.read_in_values(values_file)
annotations = deva.read_in_values(annotations_file)

# Binarize annotation columns
annotations = deva.binarize_annotations(annotations)

# Run outliers comparative analysis
outliers, qvalues = deva.run_outliers(
    values, annotations, 
    save_outlier_table=True,
    save_qvalues=True, 
    save_comparison_summaries=True
)

# Pull out results
qvalues_table = qvalues.df
vis_table = outliers.frac_table

# Make heatmaps for significant genes
for col in annotations.columns:
    axs = deva.plot_heatmap(annotations, qvalues_table, col, vis_table, savefig=True)

```

##### Command line interface
```bash
usage: deva outliers_table [-h] [--output_prefix OUTPUT_PREFIX] [--iqrs IQRS]
                           [--up_or_down {up,down}] [--ind_sep IND_SEP]
                           [--do_not_aggregate] [--write_frac_table]
                           values

Takes a table of values and converts to a table of outlier counts.

positional arguments:
  values                File path to input values. Columns must be samples,
                        genes must be sites or genes. Only .tsv and .csv
                        accepted.

optional arguments:
  -h, --help            show this help message and exit
  --output_prefix OUTPUT_PREFIX
                        Output prefix for writing files. Default outliers.
  --iqrs IQRS           Number of interquartile ranges (IQRs) above or below
                        the median to consider a value an outlier. Default is
                        1.5 IQRs.
  --up_or_down {up,down}
                        Whether to look for up or down outliers. Choices are
                        up or down. Default up.
  --ind_sep IND_SEP     If site labels have a parent molecule (e.g. a gene
                        name such as ATM) and a site identifier (e.g. S365)
                        this is the delimiter between the two elements.
                        Default is -
  --do_not_aggregate    Use flag if you do not want to sum outliers based on
                        site prefixes.
  --write_frac_table    Use flag if you want to write a table with fraction of
                        values per site, per sample that are outliers. Will
                        not be written by default. Useful for visualization.


usage: deva binarize [-h] [--output_prefix OUTPUT_PREFIX] annotations

Takes an annotation table where some columns may have more than 2 possible
values (not including empty/null values) and outputs an annotation table with
only two values per annotation. Propagates null values.

positional arguments:
  annotations           Annotation table with samples as rows and annotation
                        labels as columns.

optional arguments:
  -h, --help            show this help message and exit
  --output_prefix OUTPUT_PREFIX
                        Output prefix for writing files. Default outliers.


usage: deva compare_groups [-h] [--output_prefix OUTPUT_PREFIX]
                           [--frac_filter FRAC_FILTER]
                           [--write_comparison_summaries] [--iqrs IQRS]
                           [--up_or_down {up,down}] [--write_gene_list]
                           [--make_heatmaps] [--fdr FDR]
                           [--red_or_blue {red,blue}]
                           [--annotation_colors ANNOTATION_COLORS]
                           outliers_table annotations

Takes an annotation table and outlier count table (output of outliers_table)
and outputs qvalues from a statistical test that looks for enrichment of
outlier values in each group in the annotation table. For each value in each
comparison, the qvalue table will have 1 column, if there are any genes in
that comparison.

positional arguments:
  outliers_table        Table of outlier counts (output of outliers_table).
                        Must be .tsv or .csv file, with outlier and non-
                        outlier counts as columns, and genes/sites as rows.
  annotations           Table of annotations. Must be .csv or .tsv. Samples as
                        rows and comparisons as columns. Comparisons must have
                        only unique values (not including missing values). If
                        there are more options than that, you can use binarize
                        to prepare the table.

optional arguments:
  -h, --help            show this help message and exit
  --output_prefix OUTPUT_PREFIX
                        Output prefix for writing files. Default outliers.
  --frac_filter FRAC_FILTER
                        The minimum fraction of samples per group that must
                        have an outlier in a gene toconsider that gene in the
                        analysis. This is used to prevent a high number of
                        outlier values in 1 sample from driving a low qvalue.
                        Default 0.3
  --write_comparison_summaries
                        Use flag to write a separate file for each column in
                        the annotations table, with outlier counts in each
                        group, p-values and q-values in each group.
  --iqrs IQRS           Number of IQRs used to define outliers in the input
                        count table. Optional.
  --up_or_down {up,down}
                        Whether input outlier table represents up or down
                        outliers. Needed for output file labels. Default up
  --write_gene_list     Use flag to write a list of significantly enriched
                        genes for each value in each comparison. If used, need
                        an fdr threshold as well.
  --make_heatmaps       Use flag to draw a heatmap of signficantly enriched
                        genes for each value in each comparison. If used, need
                        an fdr threshold as well.
  --fdr FDR             FDR cut off to use for signficantly enriched gene
                        lists and heatmaps. Default 0.05
  --red_or_blue {red,blue}
                        If --make_heatmaps is called, color of values to draw
                        on heatmap. Default red.
  --annotation_colors ANNOTATION_COLORS
                        File with color map to use for annotation header if
                        --make_heatmaps is used. Must have a 'value color'
                        format for each value in annotations. Any value not
                        represented will be assigned a new color.

usage: deva visualize [-h] [--output_prefix OUTPUT_PREFIX]
                      [--annotations_to_show ANNOTATIONS_TO_SHOW [ANNOTATIONS_TO_SHOW ...]]
                      [--fdr FDR] [--red_or_blue {red,blue}]
                      [--annotation_colors ANNOTATION_COLORS]
                      [--write_gene_list]
                      comparison_qvalues annotations visualization_table
                      comparison_of_interest

Used to make custom heatmaps from significant genes.

positional arguments:
  comparison_qvalues    Table of qvalues, output from compare_groups. Must be
                        .csv or .tsv. Has genes/sites as rows and comparison
                        values as columns.
  annotations           Table of annotations used to generate qvalues.
  visualization_table   Values to visualize in heatmap. Samples as columns and
                        genes/sites as rows. Using outlier fraction table is
                        recommended, but original values can also be used if
                        no aggregation was used.
  comparison_of_interest
                        Name of column in qvalues table from which to
                        visualize significant genes.

optional arguments:
  -h, --help            show this help message and exit
  --output_prefix OUTPUT_PREFIX
                        Output prefix for writing files. Default outliers.
  --annotations_to_show ANNOTATIONS_TO_SHOW [ANNOTATIONS_TO_SHOW ...]
                        Names of columns from the annotation table to show in
                        the header of the heatmap. Default is all columns.
  --fdr FDR             FDR threshold to use to select genes to visualize.
                        Default 0.05
  --red_or_blue {red,blue}
                        Color of values to draw on heatmap. Default red.
  --annotation_colors ANNOTATION_COLORS
                        File with color map to use for annotation header. Must
                        have a line with 'value color' format for each value
                        in annotations. Any value not represented will be
                        assigned a new color.
  --write_gene_list     Use flag to write a list of significantly enriched
                        genes for each value in each comparison.


usage: deva outliers [-h] [--output_prefix OUTPUT_PREFIX] [--iqrs IQRS]
                     [--up_or_down {up,down}] [--do_not_aggregate]
                     [--write_outlier_table] [--write_frac_table]
                     [--ind_sep IND_SEP] [--frac_filter FRAC_FILTER]
                     [--write_comparison_summaries] [--fdr FDR]
                     [--write_gene_list] [--make_heatmaps]
                     [--red_or_blue {red,blue}]
                     [--annotation_colors ANNOTATION_COLORS]
                     values annotations

Runs whole outliers pipeline. Has options to output every possible output.

positional arguments:
  values                File path to input values. Samples are columns and
                        genes/sites are rows. Only .tsv and .csv accepted.
  annotations           File path to annotation values. Rows are sample names,
                        header is different annotations. e.g. mutation status.

optional arguments:
  -h, --help            show this help message and exit
  --output_prefix OUTPUT_PREFIX
                        Output prefix for writing files. Default outliers.
  --iqrs IQRS           Number of inter-quartile ranges (IQRs) above or below
                        the median to consider a value an outlier. Default is
                        1.5.
  --up_or_down {up,down}
                        Whether to look for up or down outliers. Choices are
                        up or down. Default up.
  --do_not_aggregate    Use flag if you do not want to sum outliers based on
                        site prefixes.
  --write_outlier_table
                        Use flag to write a table of outlier counts.
  --write_frac_table    Use flag if you want to write a table with fraction of
                        values per site per sample that are outliers. Useful
                        for custom visualization.
  --ind_sep IND_SEP     If site labels have a parent molecule (e.g. a gene
                        name such as ATM) and a site identifier (e.g. S365)
                        this is the delimiter between the two elements.
                        Default is -
  --frac_filter FRAC_FILTER
                        The minimum fraction of samples per group that must
                        have an outlier in a gene toconsider that gene in the
                        analysis. This is used to prevent a high number of
                        outlier values in 1 sample from driving a low qvalue.
                        Default 0.3
  --write_comparison_summaries
                        Use flag to write a separate file for each column in
                        the annotations table, with outlier counts in each
                        group, p-values and q-values in each group.
  --fdr FDR             FDR threshold to use to select genes to visualize.
                        Default 0.05
  --write_gene_list     Use flag to write a list of significantly enriched
                        genes for each value in each comparison.
  --make_heatmaps       Use flag to draw a heatmap of significantly enriched
                        genes for each value in each comparison. If used, need
                        an fdr threshold as well.
  --red_or_blue {red,blue}
                        Color of values to draw on heatmap. Default red.
  --annotation_colors ANNOTATION_COLORS
                        File with color map to use for annotation header. Must
                        have a line with 'value color' format for each value
                        in annotations. Any value not represented will be
                        assigned a new color.

``` 