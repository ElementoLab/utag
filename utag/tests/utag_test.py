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
