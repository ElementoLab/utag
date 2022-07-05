import typing as tp

import pytest
import numpy as np
import scanpy as sc
from anndata import AnnData

from utag import utag


kwargs = dict(
    slide_key="roi",
    max_dist=20,
    normalization_mode="l1_norm",
    apply_clustering=True)


@pytest.fixture
def adata() -> AnnData:
    return sc.read(
        "data/healthy_lung_adata.h5ad",
        backup_url="https://zenodo.org/record/6376767/files/healthy_lung_adata.h5ad?download=1",
    )


def check(
    utag_results: AnnData,
    probabilities: list[float],
    n: int,
    clustering: tp.Sequence[str],
) -> None:
    assert utag_results.obs.columns.str.contains("UTAG Label").any()
    for cluster in clustering:
        for prob in probabilities:
            col = f"UTAG Label_{cluster}_{prob}_probabilities"
            assert col in utag_results.obsm
            assert utag_results.obsm[col].shape[0] == n
            assert np.allclose(utag_results.obsm[col].sum(1), [1] * n)


def test_subsample_serial(adata: AnnData) -> None:
    n = 10_000
    clustering = ["leiden", "parc"]
    utag_results = utag(
        adata[:n],
        **kwargs
        clustering_method=clustering,
        parc_kwargs=dict(small_pop=10),
        resolutions=[0.3],
        parallel=False,
    )
    check(utag_results, [0.3], n, clustering)


def test_subsample_parallel(adata: AnnData) -> None:
    n = 10_000
    clustering = ["leiden", "parc"]
    utag_results = utag(
        adata[:n],
        **kwargs
        clustering_method=clustering,
        parc_kwargs=dict(small_pop=10),
        resolutions=[0.3],
        parallel=True,
    )
    check(utag_results, [0.3], n, clustering)


def test_full(adata: AnnData) -> None:
    probabilities = [0.05, 0.1, 0.3]
    n = adata.shape[0]
    utag_results = utag(
        adata,
        **kwargs
        clustering_method="leiden",
        resolutions=probabilities,
    )
    check(utag_results, probabilities, adata.shape[0], ["leiden"])
