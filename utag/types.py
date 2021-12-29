#!/usr/bin/env python

"""
Specific data types used for type annotations in the package.
"""

from __future__ import annotations
import os
import typing as tp
import pathlib

import matplotlib
from matplotlib.figure import Figure as _Figure
import numpy
import pandas
import skimage
import skimage.future
import networkx
import anndata

__all__ = [
    "Shape",
    "Shape2D",
    "Shape3D",
    "Shape4D",
    "Array",
    "Graph",
    "DataFrame",
    "Figure",
    "Axis",
    "Path",
    "AnnData",
]


class Path(pathlib.Path):
    """
    A pathlib.Path child class that allows concatenation with strings
    by overloading the addition operator.

    In addition, it implements the ``startswith`` and ``endswith`` methods
    just like in the base :obj:`str` type.

    The ``replace_`` implementation is meant to be an implementation closer
    to the :obj:`str` type.

    Iterating over a directory with ``iterdir`` that does not exists
    will return an empty iterator instead of throwing an error.

    Creating a directory with ``mkdir`` allows existing directory and
    creates parents by default.
    """

    _flavour = (
        pathlib._windows_flavour  # type: ignore[attr-defined]  # pylint: disable=W0212
        if os.name == "nt"
        else pathlib._posix_flavour  # type: ignore[attr-defined]  # pylint: disable=W0212
    )

    def __add__(self, string: str) -> Path:
        return Path(str(self) + string)

    def startswith(self, string: str) -> bool:
        return str(self).startswith(string)

    def endswith(self, string: str) -> bool:
        return str(self).endswith(string)

    def replace_(self, patt: str, repl: str) -> Path:
        return Path(str(self).replace(patt, repl))

    def iterdir(self) -> tp.Generator:
        if self.exists():
            yield from [Path(x) for x in pathlib.Path(str(self)).iterdir()]
        yield from []

    def mkdir(self, mode=0o777, parents: bool = True, exist_ok: bool = True) -> Path:
        super().mkdir(mode=mode, parents=parents, exist_ok=exist_ok)
        return self


Shape = tp.Tuple[int, ...]
Shape2D = tp.Tuple[int, int]
Shape3D = tp.Tuple[int, int, int]
Shape4D = tp.Tuple[int, int, int, int]

Array = tp.Union[numpy.ndarray]
Graph = tp.Union[networkx.Graph, skimage.future.graph.RAG]

# MultiIndexSeries = tp.Union[pandas.Series]
DataFrame = tp.Union[pandas.DataFrame]
AnnData = tp.Union[anndata.AnnData]

Figure = tp.Union[_Figure]
Axis = tp.Union[matplotlib.axis.Axis]

# Patch = tp.Union[matplotlib.patches.Patch]
# ColorMap = tp.Union[matplotlib.colors.LinearSegmentedColormap]
