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
deva 
``` 