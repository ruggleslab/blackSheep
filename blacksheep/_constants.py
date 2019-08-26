# Used in parsers
binarized_col_name = "%s_%s"  # % (col, val)
outgroup_val = "not-%s"  # % val


# Used primarily in outlierTable
row_iqr_name = "row_iqr"
row_median_name = "row_median"
row_upper_bound_name = "row_medPlus"
row_lower_bound_name = "row_medMinus"
col_seps = "_"
col_not_outlier_suffix = "notOutliers"
col_outlier_suffix = "outliers"
agg_col = "gene"


# Used primarily in comparisons
outlier_count_lab = "Outliers"
not_outlier_count_lab = "NotOutlier"
fisherp_col = "fisherp"
fisherfdr_col = "fisherFDR"
mult_hypoth_method = "fdr_bh"


# Used in outliers
fdr_col_label = "fisherFDR_%s_%s"  # % (comp, group_label)
general_group_label_0 = "0"
general_group_label_1 = "1"
comp_group_suffix = "_%s_%s"
general_fisher_p = "fisherp"
specific_fisher_p = "fisherp_%s_%s"


# Used in visualization
default_palette = "Set2"
cbar_label = "Fraction\nOutliers"
plot_title = "Outliers in %s"


# File name outputs
gene_list_file_name = "%s.%s.sig_genes.fdr%s.txt"
frac_table_file_name = "%s.%s.fraction_table.tsv"
outlier_table_file_name = "%s.%s.count_table.tsv"
ind_comparison_file_name = "%s.%s.%s.qvalues.tsv"
qvalues_file_name = "%s.%s.qvalues.tsv"
figure_file_name = "%s.%s.fdr%s.heatmap.pdf"
parameters_file_name = "%s.parameters.txt"
