# Gene Platform Utils API Documentation

## Modules

### transformation

This module provides functions for variance stabilizing transformation (VST) and other data transformations.

#### Functions

- `vst(counts: np.ndarray, min_mu: float = 0.5, min_disp: float = 1e-8, max_disp: float = 10.0, beta_tol: float = 1e-8, subset_n_genes: Optional[int] = 1200, rand_seed: Optional[int] = 123) -> NDArray`
  
  Performs variance stabilizing transformation on the input counts matrix.

- `log2_1p(counts: np.ndarray) -> np.ndarray`
  
  Applies log2(x + 1) transformation.

- `ln_1p(counts: np.ndarray) -> np.ndarray`
  
  Applies natural log(x + 1) transformation.

- `log10_1p(counts: np.ndarray) -> np.ndarray`
  
  Applies log10(x + 1) transformation.

### plot_de

This module provides functions for creating differential expression plots.

#### Functions

- `create_volcano_plot(data: np.ndarray, row_names: list[str], column_names: list[str], lfc_thresh: float, pval_thresh: float, cohort_name: str) -> str`
  
  Creates an interactive volcano plot using Plotly.

- `create_mean_difference_plot(data: np.ndarray, row_names: list[str], column_names: list[str], fdr: float, cohort_name: str) -> str`
  
  Creates an interactive mean difference plot using Plotly.

### plot_eda

This module provides functions for creating exploratory data analysis plots.

#### Functions

- `calculate_correlation(counts: np.ndarray) -> np.ndarray`
  
  Calculates the sample-to-sample correlation matrix.

- `create_correlation_heatmap(counts: np.ndarray, sample_ids: list[str]) -> str`
  
  Creates an interactive heatmap of sample-to-sample correlation using Plotly.

- `create_pca_plot(counts, sample_ids)`
  
  Creates a 3D PCA plot.

- `create_tsne_plot(counts, sample_ids)`
  
  Creates a 3D t-SNE plot.

### plot_gsea

This module provides functions for creating gene set enrichment analysis plots.

#### Functions

- `gene_concept_network_plot(gsea_res: pd.DataFrame, de_res: pd.DataFrame, genes_df: pd.DataFrame, color_metric: str = "log2FoldChange", pvalue_threshold: float = 0.05, layout_seed: int = 0, color_seed: int = 0) -> str`
  
  Creates an interactive gene concept network plot using Plotly.

### generate_html_3d_embedding

This module provides functions for generating HTML content for 3D embeddings.

#### Functions

- `generate_scattergl_html(title1: str, title2: str, data_csv_content: str, metadata_csv_content: str) -> str`
  
  Generates HTML content for a 3D scatter plot using Plotly.

## Dependencies

The package has the following dependencies:

- networkx
- numpy
- plotly
- scipy
- scikit-learn

These dependencies are specified in the `pyproject.toml` and `requirements.txt` files.


```1:15:pyproject.toml
[build-system]
requires = ["setuptools"]
build-backend = "setuptools.build_meta"

[project]
name = "gene-platform-utils"
version = "0.0.1"
description = "Python helpers for the Gene Platform (run with pyodide)"
dependencies = ["networkx", "numpy", "plotly", "scipy", "scikit-learn"]
requires-python = ">= 3.9"

classifiers = [
    "License :: OSI Approved :: MIT License",
]

```

```1:6:requirements.txt
networkx
numpy
plotly
scipy
scikit-learn

```
