# Unsupervised discovery of tissue architechture with graphs (GATDU)
<!-- 
[![Zenodo badge](https://zenodo.org/badge/doi/___doi1___.svg)](https://doi.org/___doi1___)
[![Biorxiv badge](https://zenodo.org/badge/doi/__doi1___.svg)](https://doi.org/__doi1___) ⬅️ read the preprint here
 -->

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

#### Install

```bash
#git clone https://github.com/afrendeiro/gatdu.git
git clone -b dev https://github.com/afrendeiro/gatdu.git
cd gatdu
pip install --editable .
```

#### Usage Principles

GATDU repository input requires an AnnData `adata` with preprocessed expression profiles in `adata.X` and spatial coordinates stored in `adata.obsm['spatial']`.
For more detailed usage of the package, please refer to [IMC Image Segmentation.ipynb](https://github.com/afrendeiro/gatdu/blob/dev/gatdu/documentation/IMC%20Image%20Segmentation.ipynb) in the documentation folder.
