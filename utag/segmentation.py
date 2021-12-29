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


def segment_slide(
    adata: AnnData,
    channels_to_use: Array = None,
    slide_key: str = 'Slide', # None in case only one slide
    batch_key: str = 'sample',
    save_key: str = 'UTAG Label',
    filter_by_variance:bool = False,
    normalize_adjacency_matrix: bool = False,
    max_dist: float = 1e20,
    kernel: str = None,
    sigma: float = 1.0,
    gamma: float = 10.0,
    n_pcs: int = 10,
    apply_umap: bool = False,
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
            
            ad_copy = segment_graph(ad_copy, sigma = sigma, normalize_adjacency_matrix = normalize_adjacency_matrix, kernel = kernel)
            ad_list.append(ad_copy)
        ad_result = anndata.concat(ad_list)
    else:
        sq.gr.spatial_neighbors(ad, radius = max_dist, coord_type = 'generic', set_diag = True)
        ad_result = segment_graph(ad, sigma = sigma, kernel = kernel)
    
    if apply_clustering:
        import bbknn, parc
        sc.tl.pca(ad_result, n_comps = n_pcs)
        sc.pp.neighbors(ad_result, n_pcs = n_pcs)
        #bbknn.bbknn(ad_result, batch_key = batch_key, n_pcs = n_pcs)
        if apply_umap:
            print('Running UMAP on Input Dataset...')
            sc.tl.umap(ad_result, gamma = gamma)
        
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

def segment_graph(
    adata: AnnData,
    kernel: str = 'rbf',
    sigma: float = 1.0,
    normalize_adjacency_matrix: bool = False
) -> AnnData:
    
    from scipy.linalg import sqrtm
    
    # Radial Basis Function
    if kernel == 'rbf':
        #fully_connected = adata.copy()
        # radius supposed to be very large number
        #sq.gr.spatial_neighbors(fully_connected, radius = 1e100, coord_type = 'generic')

        A = adata.obsp['spatial_distances']
        A_mod = A
        affinity = np.exp(-(A_mod/(2 * sigma**2)))
        
        #A_mod = A + np.eye(A.shape[0])
        #affinity = affinity * A_mod
        connectivity = adata.obsp['spatial_connectivities'] + np.eye(adata.obsp['spatial_connectivities'].shape[0])
        affinity = affinity * connectivity
    elif kernel == 'inverse_quadratic':
        #fully_connected = adata.copy()
        # radius supposed to be very large number
        #sq.gr.spatial_neighbors(fully_connected, radius = 1e100, coord_type = 'generic')

        A = adata.obsp['spatial_distances']
        A_mod = A
        affinity = 1/(1+sigma ** 2 * A_mod)
        
        connectivity = adata.obsp['spatial_connectivities'] + np.eye(adata.obsp['spatial_connectivities'].shape[0])
        affinity = affinity * connectivity
    elif normalize_adjacency_matrix:
        # D for A_mod:
        A = adata.obsp['spatial_connectivities']
        A_mod = A + np.eye(A.shape[0])
        
        '''
        D_mod = np.zeros_like(A_mod)
        np.fill_diagonal(D_mod, A_mod.sum(axis=1).flatten())

        # Inverse square root of D:
        D_mod_invroot = np.linalg.inv(sqrtm(D_mod))
        
        # technically not affinity
        affinity = D_mod_invroot @ A_mod @ D_mod_invroot
        '''
        #affinity = A_mod / A_mod.sum(axis=1)[:,np.newaxis]
        from sklearn.preprocessing import normalize
        affinity = normalize(A_mod, axis=1, norm='l1')
        
    else:
        # Plain A_mod multiplication
        A = adata.obsp['spatial_connectivities']
        affinity = A
        
    adata.X = affinity @ adata.X
    return adata


def segment_graph_backup(
    adata: AnnData,
    sigma: float = 10.0,
    n_pcs: int = 10,
    resolution: float = 0.5,
    kernel: str = 'rbf',
    apply_clustering: bool = True
) -> AnnData:
    
    #A = adata.obsp['spatial_distances']
    A = adata.obsp['spatial_distances']
    A_mod = A + np.eye(A.shape[0])
    
    # Radial Basis Function
    if kernel == 'rbf':
        affinity = np.exp(-((sigma)**2 * A_mod))
    if kernel == 'inverse_quadratic':
        affinity = 1/(1+sigma ** 2 * A_mod)
    
    adata.X = affinity @ adata.X
    
    if apply_clustering:
        sc.pp.neighbors(adata, n_pcs = n_pcs)
        sc.tl.umap(adata)
        sc.tl.leiden(adata, resolution = resolution)
        freq = adata.obs.groupby('leiden').count()
        
        freq = freq / freq.sum()
        low_freq = freq.loc[freq['area'] < 0.01,:].index
        low_freq = {key: 'Empty' for key in low_freq}
        adata.obs['leiden'] = adata.obs['leiden'].replace(low_freq)

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