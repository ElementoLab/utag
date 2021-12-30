import scanpy as sc
import squidpy as sq
import os

from utag.types import Path, Array, AnnData
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from tqdm import tqdm
import warnings
from scipy.sparse import (spdiags, SparseEfficiencyWarning, csc_matrix,
    csr_matrix, isspmatrix, dok_matrix, lil_matrix, bsr_matrix)

import anndata
warnings.simplefilter('ignore',SparseEfficiencyWarning)

def UTAG(
    adata: AnnData,
    channels_to_use: Array = None,
    slide_key: str = 'Slide', # None in case only one slide
    save_key: str = 'UTAG Label',
    filter_by_variance: bool = False,
    max_dist: float = 20, # Recommended value of 20 ~ 50 depending on magnification
    normalization_mode: str = 'l1_norm',
    parc_pcs: int = 10,
    apply_umap: bool = False,
    umap_gamma: float = 10.0,
    apply_clustering: bool = True,
    clustering_method: Array = ['leiden', 'parc'],
    resolutions: Array = [0.05, 0.1, 0.3, 1.0]
) -> AnnData:
    '''
    Parameters:
        adata: AnnData - AnnData with spatial key to build graph
        max_dist: float - Maximum distance to cut edges within a graph
        slide_key(opt): str - batch key that needs to be provided in case where the input AnnData comes from more than one image
    '''
    
    ad = adata.copy()
    
    if channels_to_use:
        ad = ad[:,channels_to_use]
        
    if filter_by_variance:
        ad = low_variance_filter(ad)
        
    if isinstance(clustering_method, list):
        clustering_method = [m.upper() for m in clustering_method]
    elif isinstance(clustering_method, str):
        clustering_method = [clustering_method.upper()]
    else:
        print('Invalid Clustering Method. Clustering Method Should Either be a string or a list')
        return
    assert('LEIDEN' in clustering_method or 'PARC' in clustering_method)
        
    print('Applying UTAG Algorithm...')
    if slide_key:
        ad_ = ad
        ad_list = []
        for slide in tqdm(ad_.obs[slide_key].unique()):
            ad = ad_[ad_.obs[slide_key] == slide]
            ad_copy = ad.copy()

            sq.gr.spatial_neighbors(ad_copy, radius = max_dist, coord_type = 'generic', set_diag = True)
            
            ad_copy = custom_message_passing(ad_copy, mode = normalization_mode)
            ad_list.append(ad_copy)
        ad_result = anndata.concat(ad_list)
    else:
        sq.gr.spatial_neighbors(ad, radius = max_dist, coord_type = 'generic', set_diag = True)
        ad_result = custom_message_passing(ad, mode = normalization_mode)
    
    if apply_clustering:
        sc.tl.pca(ad_result, n_comps = parc_pcs)
        sc.pp.neighbors(ad_result, n_pcs = parc_pcs)

        if apply_umap:
            print('Running UMAP on Input Dataset...')
            sc.tl.umap(ad_result, gamma = umap_gamma)
        
        for resolution in tqdm(resolutions):
            
            res_key1 = save_key + '_leiden_' + str(resolution)
            res_key2 = save_key + '_parc_' + str(resolution)
            if 'LEIDEN' in clustering_method:
                print(f'Applying Leiden Clustering at Resolution: {resolution}...')
                sc.tl.leiden(ad_result, resolution = resolution, key_added = res_key1)
            
            if 'PARC' in clustering_method:
                print(f'Applying PARC Clustering at Resolution: {resolution}...')
                Parc1 = parc.PARC(ad_result.obsm['X_pca'], neighbor_graph=ad_result.obsp["connectivities"], resolution_parameter = resolution, random_seed = 1, small_pop = 1000)
                Parc1.run_PARC()
                ad_result.obs[res_key2] = pd.Categorical(Parc1.labels)
                ad_result.obs[res_key2] = ad_result.obs[res_key2].astype('category')
            
    return ad_result 

def custom_message_passing(
    adata: AnnData,
    mode: str = 'l1_norm'
) -> AnnData:
    
    from scipy.linalg import sqrtm
    if mode == 'l1_norm':
        A = adata.obsp['spatial_connectivities']
        A_mod = A + np.eye(A.shape[0])
        
        from sklearn.preprocessing import normalize
        affinity = normalize(A_mod, axis=1, norm='l1')
    else:
        # Plain A_mod multiplication
        A = adata.obsp['spatial_connectivities']
        affinity = A
        
    adata.X = affinity @ adata.X
    return adata

def low_variance_filter(
    adata: AnnData
) -> AnnData:
    return adata[:,adata.var['std']  > adata.var['std'].median()]

def evaluate_performance(
    adata: AnnData,
    batch_key: str = 'Slide',
    truth_key: str = 'DOM_argmax',
    pred_key: str = 'cluster',
    method: str = 'rand'
) -> Array:
    assert(method in ['rand', 'homogeneity'])
    from sklearn.metrics import rand_score, homogeneity_score
           
    score_list = []
    for key in adata.obs[batch_key].unique():
        batch = adata[adata.obs[batch_key] == key]
        if method == 'rand':
            score = rand_score(batch.obs[truth_key], batch.obs[pred_key])
        elif method == 'homogeneity':
           score = homogeneity_score(batch.obs[truth_key], batch.obs[pred_key])
        score_list.append(score)
    return score_list