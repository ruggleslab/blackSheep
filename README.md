# BlackSheep
##### A tool for differential extreme-value analysis

### Installation
With pip
```bash
pip install blksheep
```
With conda
```bash
conda install -c bioconda blksheep
```

### Requirements (automatically taken care of with pip and conda)
pandas  
numpy  
matplotlib  
seaborn  
scipy  
scikit-learn  
statsmodels  

### Documentation
https://blacksheep.readthedocs.io/en/latest/index.html

### Usage
##### In python
```python
import blacksheep

# Read in data
values_file = '' #insert values file here
annotations_file = '' #insert annotations file here
values = deva.read_in_values(values_file)
annotations = deva.read_in_values(annotations_file)

# Binarize annotation columns
annotations = blacksheep.binarize_annotations(annotations)

# Run outliers comparative analysis
outliers, qvalues = blacksheep.deva(
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
    axs = blacksheep.plot_heatmap(annotations, qvalues_table, col, vis_table, savefig=True)

# Normalize values
phospho = blacksheep.read_in_values('') #Fill in file here
protein = blacksheep.read_in_values('') #Fill in file here


```

##### Command line interface
*Example*
```bash
blacksheep binarize annotations.tsv --output_prefix annotations_test
blacksheep deva values.csv annotations_test.binarized.tsv --output_prefix test \
--write_outlier_table --write_comparison_summaries --write_gene_list \
--make_heatmaps
```

*Full help*  
Just make the outliers table:
```bash
usage: blacksheep outliers_table [-h] [--output_prefix OUTPUT_PREFIX] [--iqrs IQRS]
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
```

Binarize the columns in an annotations table.   
**Warning: do not include non-categorical columns, or columns you don't want binarized. You'll
 end up with a huge un-wieldly table. **
```bash
usage: blacksheep binarize [-h] [--output_prefix OUTPUT_PREFIX] annotations

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
```

Compare all the groups described in columns of an annotation table using outlier counts
```bash
usage: blacksheep compare_groups [-h] [--output_prefix OUTPUT_PREFIX]
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
```

Make heatmaps visualizing enriched genes for each group in an annotation table
```bash
usage: blacksheep visualize [-h] [--output_prefix OUTPUT_PREFIX]
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
```

Run the whole pipeline: call outliers, perform comparisons on all groups in an annotation table
, optionally make heatmaps for each group.
```bash
usage: blacksheep deva [-h] [--output_prefix OUTPUT_PREFIX] [--iqrs IQRS]
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

For finding the value differences that cannot be explained by a different data level. For example
, this can be used to find out variation due to differential phosphorylation (phospho as
 target_values) not due to protein abundance variation (protein as normalizer_values).   
**Warning: Row IDs between the two tables must match**  
```bash
usage: blacksheep normalize [-h] [--ind_sep IND_SEP] [--output_prefix OUTPUT_PREFIX]
                      target_values normalizer_values

Takes a target table and a normalizer table, and returns a normalized target
table. Builds a regularized linear model for each line in the target table
using the matching row ID in the normalizer table, and finds the residuals of
that model for each value. for example, this could be used to normalize
phospho-peptide data by protein abundance data; resulting values will reflect
only abundance differences due to phosphorylation changes, not peptide
abundances. Another use could be normalizing RNA by CNA.

positional arguments:
  target_values         Table of values to be normalized. Sites/genes as rows,
                        samples as columns. Row identifiers must be unique.
  normalizer_values     Table of values to use for normalization. Sites/genes
                        as rows, samples as columns. Row identifiers must be
                        unique, and must match the pre-ind_sep part of the
                        target values identifiers.

optional arguments:
  -h, --help            show this help message and exit
  --ind_sep IND_SEP     Separator used in index if target is site specific.
                        Row IDs before ind_sep in the target must match the
                        row IDs in normalizer_values. If row IDs already
                        match, leave blank.
  --output_prefix OUTPUT_PREFIX
                        Prefix for output file. Suffix will be
                        '.normalized.tsv'

```

For a more thorough vignette, refer to our [supplementary notebooks](https://github.com/ruggleslab/blacksheep_supp)
