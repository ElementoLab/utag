# Unsupervised discovery of tissue architechture with graphs (UTAG)
<!-- 
[![Zenodo badge](https://zenodo.org/badge/doi/___doi1___.svg)](https://doi.org/___doi1___)
[![Biorxiv badge](https://zenodo.org/badge/doi/__doi1___.svg)](https://doi.org/__doi1___) ⬅️ read the preprint here
 -->
 
This package implements a microanatomical structure segmentation method and outlines possible downstream structural analysis for multiplexed histology images such as Imaging Mass Cytometry (IMC), CODEX, Multiplexed Ion Beam Imaging by Time Of Flight (MIBI-TOF), Cyclic Immunofluorescence (CyCIF), and etc in ```python```. 

## Getting Started
### Install using pip
```
pip install utag
```
[comment]: <> (Needs to be uploaded to pypi)

### Install by cloning repository and running setup.py
```
git clone https://github.com/ElementoLab/utag.git
python3 setup.py install
```
[comment]: <> (Need to check whether this works)

#### Requirements

- Python 3.7+ (was run on 3.8.2)
- Python packages:
  - numpy
  - pandas
  - scipy
  - anndata
  - scanpy
  - squidpy
  - tifffile
  - scikit-image
  - tensorflow
  - matplotlib
  - tqdm
  - parmap


## Basic Usage Principles

To run the method on a single slide:
```python
from utag import utag
utag_results = utag(
    adata,
    slide_key=None,
    max_dist=20,
    normalization_mode='l1_norm',
    apply_clustering=True,
    clustering_method = 'leiden', 
    resolutions = [0.05, 0.1, 0.3]
)
```

To run the method on multiple slides in batch mode:
```python
from utag import utag
utag_results = utag(
    adata,
    slide_key="roi",
    max_dist=20,
    normalization_mode='l1_norm',
    apply_clustering=True,
    clustering_method = 'leiden', 
    resolutions = [0.05, 0.1, 0.3]
)
```

To vizually inspect the results of the method:
```python
import scanpy as sc
for roi in utag_results.obs['roi'].unique():
    result = utag_results[utag_results.obs['roi'] == roi].copy()
    sc.pl.spatial(result, color = 'UTAG Label_leiden_0.1', spot_size = 10)
```


## Key Parameters
| Input Parameter | Description |
| ---------- |----------|
| `adata` | (`anndata.AnnData`) n_cells x n_features. `AnnData` of cells with spatial coordinates stored in `adata.obsm['spatial']` as `numpy.ndarray`. |
| `max_dist` | (`float`, default = 20.0) Threshold euclidean distance to determine whether a pair of cell is adjacent in graph structure. Recommended values are between 10 to 100 depending on magnification. |
| `slide_key` | (`str`, optional, default = 'Slide') Key required for running UTAG across multiple images. Unique image identifiers should be placed under `adata.obs`. Use `None` to run UTAG on a single slide. |
| `save_key` | (`str`, default = 'UTAG Label') Key to be added to adata object holding the UTAG clusters. Depending on the values of `clustering_method` and `resolutions`, the final keys will be of the form: {save_key}\_{method}\_{resolution}". |
| `normalization_mode` |  (`str`, default = 'l1_norm') Method to normalize adjacency matrix. 'l1_norm' will behave as mean-aggregation during message passing. Default is 'l1_norm'. Any other value will not perform normalization, leading to a sum-aggregation. |
| `apply_clustering` |  (bool, default = True) Whether to cluster the message passed matrix. |
| `clustering_method` |  (Sequence[str], default = ['leiden', 'parc']) Which clustering method(s) to use for clustering of the message passed matrix. |
| `resolutions` |  (Sequence[float], default = [0.05, 0.1, 0.3, 1.0]) Resolutions the methods in `clustering_method` should be run at. |



For more detailed usage of the package and downstream analysis, please refer to [IMC Healthy Lung.ipynb](https://github.com/ElementoLab/utag/blob/main/documentation/IMC%20Healthy%20Lung.ipynb) in the documentation folder.