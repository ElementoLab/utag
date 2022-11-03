# Unsupervised discovery of tissue architechture with graphs (UTAG)
[![Zenodo badge](https://zenodo.org/badge/doi/10.1038/s41592-022-01657-2.svg)](https://doi.org/10.1038/s41592-022-01657-2) ⬅️ read the published article here <br>
[![Biorxiv badge](https://zenodo.org/badge/doi/10.1101/2022.03.15.484534.svg)](https://doi.org/10.1101/2022.03.15.484534) ⬅️ read the preprint here <br>
[![Zenodo badge](https://zenodo.org/badge/doi/10.5281/zenodo.6376767.svg)](https://doi.org/10.5281/zenodo.6376767) ⬅️ Preprocessed Multiplexed Image Data and UTAG results <br>
 
 
This package implements segmentation of multiplexed imaging data into microanatomical domains.
Multiplexed imaging data types are typically imaging mass cytometry (IMC), co-detection by indexing (CODEX), multiplexed ion beam imaging by time of flight (MIBI-TOF), cyclic immunofluorescence (CyCIF), and others.
The package also provides functions for the downstream analysis of the detected micro-anatomical structure.


## Getting Started

### Install from github

```bash
pip install git+https://github.com/ElementoLab/utag.git@main
```
Installation should take less than 10 seconds.

#### Requirements
There are no specific hardware requirements.

Software requirements:
- UTAG has been tested on Mac OS, Linux, and Windows [WSL](https://docs.microsoft.com/en-us/windows/wsl/about) (Windows Subsystem for Linux).
- Python 3.7+ (tested on 3.8.2)
- Python packages (automatically installed by `pip`):
  - numpy
  - pandas
  - anndata
  - scanpy
  - parc
  - squidpy
  - scipy
  - matplotlib
  - tqdm
  - networkx
  - parmap
  - scikit-learn
  - setuptools_scm (may require manual installation for MacOS through pip)

Specific versions of Python packages have been pinned to the [setup.cfg](setup.cfg) file.

## Basic Usage Principles

The UTAG process can be run with a single function call `utag.utag`.
The input is a [AnnData](https://anndata.readthedocs.io/) object which should have the position of cells (typically centroids) in the `spatial` slot of `adata.obsm`.
The function will output domain classes for each cell stored in the `obs` slot of the returned AnnData object.


### Running an example/demo dataset

Please refer to the [notebook directory](documentation/), and to the notebook on [running UTAG on healthy lung data](documentation/IMC%20Healthy%20Lung.ipynb) for a reproducible example.
All data and respective results used for analysis can be downloaded from [![Zenodo badge](https://zenodo.org/badge/doi/10.5281/zenodo.6376767.svg)](https://doi.org/10.5281/zenodo.6376767).

All data could alternatively be downloaded through command line:
```bash
pip install zenodo_get
zenodo_get -d 10.5281/zenodo.6376767 -o data
```

### Running on your data

To run the method on multiple images/slides in batch mode:
```python
from utag import utag

# Use Scanpy to get a h5ad file with provided data
import scanpy as sc
adata = sc.read(
    'data/healthy_lung_adata.h5ad',
    backup_url='https://zenodo.org/record/6376767/files/healthy_lung_adata.h5ad?download=1')

# Run UTAG on provided data
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
UTAG should take around \~2 min on a local machine for the batch mode on the data.

To run the method on a single image, pass `None` to the slide_key argument:
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

To visually inspect the results of the method:
```python
import scanpy as sc
for roi in utag_results.obs['roi'].unique():
    result = utag_results[utag_results.obs['roi'] == roi].copy()
    sc.pl.spatial(result, color = 'UTAG Label_leiden_0.1', spot_size = 10)
```

## User Guide on UTAG (Hyperparameters to test and tune)


Although UTAG greatly reduces manual labor involved in segmentation of microanatomical domains across, successful application of UTAG depends on three key user inputs. First is the `max_dist` parameter which defines the threshold distance. Second is the clustering resolution (under list of `resolutions`) to determine the coarsity of the clustering. Last is user interpretation of the resulting clusters to identify the structure.

We intentionally leave the optimization of max_dist open to users to maximize the applicability of UTAG to unseen datasets. This is because this parameter is tightly related with the resolution or magnification of the data under use. In our manuscript, we apply the method on IMC data and optical imaging-based CyCIF, which have different per unit area pixel densities. In the case of IMC, we suggest that a well working max_dist is between 10 and 20 as 1 pixel exactly maps to 1 micrometer. With an imaging-based technique like CyCIF, the optimal distance can vary with magnification, focal lengths, distance to tissue, and other factors, which make it hard to suggest a one-fits-all rule. Also there might be nuanced differences for the exact tissue of interest that may vary across specimens under examination.

We believe that the optimal clustering resolution is a hyperparameter that should be explored to suit their biological question of interest. For such reasons, we provide a list of resolutions as default to be explored by the users. A general rule here is that increasing the resolution parameter will return more refined substructures, while decreasing it will return coarser, more broad structures. We also recommend users to use a higher resolution parameter when screening for a rare structure, as a higher resolution will capture more structures, and vice versa. In our benchmarking, we saw that with the exception of extreme hyperparameter values, UTAG’s performance was fairly robust across various clustering resolutions (Extended Data Figure S3).

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
| `parallel` | Whether the message passing part of the method should be parallelized. This is done using the `parmap` package and the `multiprocessing` module from the standard library. |

For more detailed usage of the package and downstream analysis, please refer to [IMC Healthy Lung.ipynb](documentation/IMC%20Healthy%20Lung.ipynb) in the documentation folder.


## Running UTAG on R

To make UTAG available to R users, we port the python code to using the `reticulate` package. The code was tested under after installing `UTAG` natively for python 3.8.10 and under conda environment with python 3.10.4, on Ubuntu 20.04.3 LTS and R 4.2.0.

Nonetheless, we highly recommend that users use our package in python for more involved analysis as the package has been developed and tested more thoroughly in python.

```R
install.packages('reticulate')
library(reticulate)

# grab python interpreter
use_python('/usr/bin/python3') # in case UTAG is installed on native python
# use_condaenv('utag') # in case UTAG is installed using conda
# use_virtualevn('utag') # in case UTAG is installed under virtualenv


# import necessary python packages
utag <- import('utag')
scanpy <- import('scanpy')

# read anndata with cell expressions and and locations 
adata <- scanpy$read(
    'data/healthy_lung_adata.h5ad',
    backup_url='https://zenodo.org/record/6376767/files/healthy_lung_adata.h5ad?download=1'
)

# show general content of the data
print(adata)

# run UTAG
utag_results <- utag$utag(
    adata,
    slide_key = "roi",
    max_dist = 20,
    normalization_mode = 'l1_norm',
    apply_clustering = TRUE,
    clustering_method = 'leiden', 
    resolutions = c(0.05, 0.1, 0.3)
)

# show content of the data that now includes UTAG results for various resolutions
print(utag_results)
```
Also available as a [R markdown file in](documentation/Running%20in%20R.Rmd).

## Development

We are happy to receive community contributions to UTAG through pull requests on Github.

Please run tests after re-installing the package:
```bash
pytest --pyargs utag
```
or before by:
```bash
pytest utag/tests/utag_test.py
```
