import pandas
import anndata
import numpy as np
from utag import utag

clinical_data = pandas.read_csv("clinical_data.tsv", sep='\t')
sc = pandas.read_csv("SingleCells.csv")
sc = sc.rename(columns={'metabric_id': 'sample', 'ImageNumber': 'slide', 'ObjectNumber': 'obj_id'})

np.random.seed(0)
slides = sc.slide.unique()[np.random.random_sample(len(sc.slide.unique())) > 0.9]

sc = sc[sc.slide.isin(slides)]

clinical_data = clinical_data.rename(columns={'Sample ID': 'sample'})
x_cols = ['Histone H3', 'SMA', 'CK5', 'CD38', 'HLA-DR', 'CK8-18', 'CD15', 'FSP1',
          'CD163', 'ICOS', 'OX40', 'CD68', 'HER2 (3B5)', 'CD3', 'Podoplanin',
          'CD11c', 'PD-1', 'GITR', 'CD16', 'HER2 (D8F12)', 'CD45RA', 'B2M',
          'CD45RO', 'FOXP3', 'CD20', 'ER', 'CD8', 'CD57', 'Ki-67', 'PDGFRB',
          'Caveolin-1', 'CD4', 'CD31-vWF', 'CXCL12', 'HLA-ABC', 'panCK',
          'c-Caspase3', 'DNA1', 'DNA2',
          'AreaShape_Area']

X = sc[x_cols]

coords = sc[['Location_Center_X', 'Location_Center_Y']].to_numpy()

clinical_data = clinical_data[clinical_data['sample'].isin(sc['sample'].unique())]

just_names = pandas.DataFrame(sc[['sample', 'slide']])

obs = just_names.merge(clinical_data, on='sample', how='left')

obs = obs.reset_index()
del obs['index']
X = X.reset_index()
del X['index']

adata = anndata.AnnData(X=X, obs=obs, obsm={'spatial': coords})

utag_results = utag(
    adata,
    slide_key="slide",
    max_dist=20,
    normalization_mode='l1_norm',
    apply_clustering=True,
    clustering_method='leiden',
    resolutions=[0.05],
    parallel=False
)
