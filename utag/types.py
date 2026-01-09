#!/usr/bin/env python

"""
Specific data types used for type annotations in the package.
"""

import typing as tp
from pathlib import Path


import numpy
import pandas
import anndata
import networkx
import matplotlib
from matplotlib.figure import Figure as _Figure


__all__ = [
    "Path",
    "Array",
    "Graph",
    "DataFrame",
    "Figure",
    "Axis",
    "AnnData",
]


Array = tp.Union[numpy.ndarray]
Graph = tp.Union[networkx.Graph]

DataFrame = tp.Union[pandas.DataFrame]
AnnData = tp.Union[anndata.AnnData]

Figure = tp.Union[_Figure]
Axis = tp.Union[matplotlib.axis.Axis]
