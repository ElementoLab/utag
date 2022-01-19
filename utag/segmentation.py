import typing as tp
import warnings
import os

import scanpy as sc
import squidpy as sq
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from tqdm import tqdm
import anndata
import parmap

from utag.types import Path, Array, AnnData


def utag(
    adata: AnnData,
    channels_to_use: tp.Sequence[str] = None,
    slide_key: tp.Optional[str] = "Slide",
    save_key: str = "UTAG Label",
    filter_by_variance: bool = False,
    max_dist: float = 20.0,
    normalization_mode: str = "l1_norm",
    pca_kwargs: tp.Dict[str, tp.Any] = dict(n_comps=10),
    apply_umap: bool = False,
    umap_kwargs: tp.Dict[str, tp.Any] = dict(),
    apply_clustering: bool = True,
    clustering_method: tp.Sequence[str] = ["leiden", "parc"],
    resolutions: tp.Sequence[float] = [0.05, 0.1, 0.3, 1.0],
    parallel: bool = True
) -> AnnData:
    """
    Discover tissue architechture in single-cell imaging data
    by combining phenotypes and positional information of cells.

    Parameters
    ----------
    adata: AnnData
        AnnData object with spatial positioning of cells in obsm 'spatial' slot.
    channels_to_use: Optional[Sequence[str]]
        An optional sequence of strings used to subset variables to use.
        Default (None) is to use all variables.
    max_dist: float
        Maximum distance to cut edges within a graph.
        Should be adjusted depending on resolution of images.
        For imaging mass cytometry, where resolution is 1um, 20 often gives good results.
        Default is 20.
    slide_key: {str, None}
        Key of adata.obs containing information on the batch structure of the data.
        In general, for image data this will often be a variable indicating the image
        so image-specific effects are removed from data.
        Default is "Slide".
    save_key: str
        Key to be added to adata object holding the UTAG clusters.
        Depending on the values of `clustering_method` and `resolutions`,
        the final keys will be of the form: {save_key}_{method}_{resolution}".
        Default is "UTAG Label".
    filter_by_variance: bool
        Whether to filter vairiables by variance.
        Default is False, which keeps all variables.
    max_dist: float
        Recommended values are between 20 to 50 depending on magnification.
        Default is 20.
    normalization_mode: str
        Method to normalize adjacency matrix.
        Default is "l1_norm", any other value will not use normalization.
    pca_kwargs: Dict[str, Any]
        Keyword arguments to be passed to scanpy.pp.pca for dimensionality reduction after message passing.
        Default is to pass n_comps=10, which uses 10 Principal Components.
    apply_umap: bool
        Whether to build a UMAP representation after message passing.
        Default is False.
    umap_kwargs: Dict[str, Any]
        Keyword arguments to be passed to scanpy.tl.umap for dimensionality reduction after message passing.
        Default is 10.0.
    apply_clustering: bool
        Whether to cluster the message passed matrix.
        Default is True.
    clustering_method: Sequence[str]
        Which clustering method(s) to use for clustering of the message passed matrix.
        Default is ["leiden", "parc"].
    resolutions: Sequence[float]
        What resolutions should the methods in `clustering_method` be run at.
        Default is [0.05, 0.1, 0.3, 1.0].
    parallel: bool
        Whether to run message passing part of algorithm in parallel.
        Will accelerate the process but consume more memory.
        Default is True.

    Returns
    -------
    adata: AnnData
        AnnData object with UTAG domain predictions for each cell in adata.obs, column `save_key`.
    """
    ad = adata.copy()

    if channels_to_use:
        ad = ad[:, channels_to_use]

    if filter_by_variance:
        ad = low_variance_filter(ad)

    if isinstance(clustering_method, list):
        clustering_method = [m.upper() for m in clustering_method]
    elif isinstance(clustering_method, str):
        clustering_method = [clustering_method.upper()]
    else:
        print(
            "Invalid Clustering Method. Clustering Method Should Either be a string or a list"
        )
        return
    assert "LEIDEN" in clustering_method or "PARC" in clustering_method

    if "PARC" in clustering_method:
        from parc import PARC  # early fail if not available

    print("Applying UTAG Algorithm...")
    if slide_key:

        if not parallel:
            ad_ = ad
            ad_list = []
            for slide in tqdm(ad_.obs[slide_key].unique()):
                ad = ad_[ad_.obs[slide_key] == slide]
                ad_copy = ad.copy()

                sq.gr.spatial_neighbors(
                    ad_copy, radius=max_dist, coord_type="generic", set_diag=True
                )

                ad_copy = custom_message_passing(ad_copy, mode=normalization_mode)
                ad_list.append(ad_copy)
        else:
            # Parallel mode will consume more memory
            ads = [ad[ad.obs[slide_key] == slide].copy() for slide in ad.obs[slide_key].unique()]
            ad_list = parmap.map(
                _parallel_message_pass,
                ads,
                radius=max_dist, coord_type="generic", set_diag=True,
                mode=normalization_mode,
                pm_pbar=True)
        ad_result = anndata.concat(ad_list)
    else:
        sq.gr.spatial_neighbors(ad, radius=max_dist, coord_type="generic", set_diag=True)
        ad_result = custom_message_passing(ad, mode=normalization_mode)

    if apply_clustering:
        sc.tl.pca(ad_result, **pca_kwargs)
        sc.pp.neighbors(ad_result)

        if apply_umap:
            print("Running UMAP on Input Dataset...")
            sc.tl.umap(ad_result, **umap_kwargs)

        for resolution in tqdm(resolutions):

            res_key1 = save_key + "_leiden_" + str(resolution)
            res_key2 = save_key + "_parc_" + str(resolution)
            if "LEIDEN" in clustering_method:
                print(f"Applying Leiden Clustering at Resolution: {resolution}...")
                sc.tl.leiden(ad_result, resolution=resolution, key_added=res_key1)

            if "PARC" in clustering_method:
                from parc import PARC

                print(f"Applying PARC Clustering at Resolution: {resolution}...")
                Parc1 = PARC(
                    ad_result.obsm["X_pca"],
                    neighbor_graph=ad_result.obsp["connectivities"],
                    resolution_parameter=resolution,
                    random_seed=1,
                    small_pop=1000,
                )
                Parc1.run_PARC()
                ad_result.obs[res_key2] = pd.Categorical(Parc1.labels)
                ad_result.obs[res_key2] = ad_result.obs[res_key2].astype("category")

    return ad_result


def _parallel_message_pass(
    ad: AnnData,
    radius: int,
    coord_type: str,
    set_diag: bool,
    mode: str,
):
    sq.gr.spatial_neighbors(
        ad, radius=radius, coord_type=coord_type, set_diag=set_diag
    )
    ad = custom_message_passing(ad, mode=mode)
    return ad


def custom_message_passing(adata: AnnData, mode: str = "l1_norm") -> AnnData:

    from scipy.linalg import sqrtm

    if mode == "l1_norm":
        A = adata.obsp["spatial_connectivities"]
        A_mod = A + np.eye(A.shape[0])

        from sklearn.preprocessing import normalize

        affinity = normalize(A_mod, axis=1, norm="l1")
    else:
        # Plain A_mod multiplication
        A = adata.obsp["spatial_connectivities"]
        affinity = A

    adata.X = affinity @ adata.X
    return adata


def low_variance_filter(adata: AnnData) -> AnnData:
    return adata[:, adata.var["std"] > adata.var["std"].median()]


def evaluate_performance(
    adata: AnnData,
    batch_key: str = "Slide",
    truth_key: str = "DOM_argmax",
    pred_key: str = "cluster",
    method: str = "rand",
) -> Array:
    assert method in ["rand", "homogeneity"]
    from sklearn.metrics import rand_score, homogeneity_score

    score_list = []
    for key in adata.obs[batch_key].unique():
        batch = adata[adata.obs[batch_key] == key]
        if method == "rand":
            score = rand_score(batch.obs[truth_key], batch.obs[pred_key])
        elif method == "homogeneity":
            score = homogeneity_score(batch.obs[truth_key], batch.obs[pred_key])
        score_list.append(score)
    return score_list
