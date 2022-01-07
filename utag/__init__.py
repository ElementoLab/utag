import warnings

from scipy.sparse import SparseEfficiencyWarning

from utag.segmentation import utag

warnings.simplefilter("ignore", FutureWarning)
warnings.simplefilter("ignore", SparseEfficiencyWarning)
