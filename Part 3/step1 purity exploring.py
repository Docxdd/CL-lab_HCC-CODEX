# Load necessary libraries
import sys
import os
import anndata as ad
import pandas as pd
import scanpy as sc
import seaborn as sns; sns.set(color_codes=True)
import matplotlib.pyplot as plt

# Import Scimap
import scimap as sm

#leidenalg

# Set the working directory
os.chdir (r"E:\myPythonProject\scimapHCCcodex20240308")






# Load data
rawdata = pd.read_csv (r'E:\11. CODEX\TrainingData\20220308\AllCN_r50_CN40_0331\All_CN_TP_0424.csv') 
rawdata.columns
temp = list(rawdata.columns)
len(temp)
for index, value in enumerate(temp):
  print(index, " ::  ", value)



select_cols_data=rawdata.columns[10:45]
data=rawdata[select_cols_data]

select_cols_meta = rawdata.columns[list(range(1, 9)) + list(range(45, 58))]
meta = rawdata[select_cols_meta]


# combine the data and metadata file to generate the AnnData object
adata = ad.AnnData (data)
adata.obs = meta

# Print adata to check for it's content
adata
adata.obs # prints the meta data
adata.X # prints the counts table
adata.var # prints the channel or marker names


# Summary of the phenotyping
#
adata.obs.columns
adata.obs['celltype'].value_counts()
adata.obs['Class'].value_counts()
adata.obs.rename(columns={'Class':"imageid"}, inplace=True)
adata.obs.rename(columns={'celltype':"phenotype"}, inplace=True)

##reindex

adata.obs_names = rawdata['cellid'].apply(lambda x: 'cellid' + str(int(x)))  ####需要把index改成字符串，这边spatial_cluster才跑的出来
adata.raw = adata
### Save AnnData object to h5ad file
adata.write('adata_original_0308.h5ad')



########reload dat
adata = ad.read('adata_original_0308.h5ad')
adata.obs_keys()

# Analyze spatial aggregation purity with celltype


adata = sm.tl.spatial_aggregate(adata, x_coordinate= "XMin",
                                y_coordinate= "YMin",
                                phenotype='phenotype', method='radius', radius=50, purity=90,
                                imageid='imageid', subset=None, label='spatial_aggregate_radius_purity90')


adata = sm.tl.spatial_aggregate(adata, x_coordinate= "XMin",
                                y_coordinate= "YMin",
                                phenotype='phenotype', method='radius', radius=50, purity=80,
                                imageid='imageid', subset=None, label='spatial_aggregate_radius_purity80')


adata = sm.tl.spatial_aggregate(adata, x_coordinate= "XMin",
                                y_coordinate= "YMin",
                                phenotype='phenotype', method='radius', radius=50, purity=70,
                                imageid='imageid', subset=None, label='spatial_aggregate_radius_purity70')

adata = sm.tl.spatial_aggregate(adata, x_coordinate= "XMin",
                                y_coordinate= "YMin",
                                phenotype='phenotype', method='radius', radius=50, purity=60,
                                imageid='imageid', subset=None, label='spatial_aggregate_radius_purity60')


adata = sm.tl.spatial_aggregate(adata, x_coordinate= "XMin",
                                y_coordinate= "YMin",
                                phenotype='phenotype', method='radius', radius=50, purity=50,
                                imageid='imageid', subset=None, label='spatial_aggregate_radius_purity50')


adata = sm.tl.spatial_aggregate(adata, x_coordinate= "XMin",
                                y_coordinate= "YMin",
                                phenotype='phenotype', method='radius', radius=50, purity=40,
                                imageid='imageid', subset=None, label='spatial_aggregate_radius_purity40')

adata = sm.tl.spatial_aggregate(adata, x_coordinate= "XMin",
                                y_coordinate= "YMin",
                                phenotype='phenotype', method='radius', radius=50, purity=30,
                                imageid='imageid', subset=None, label='spatial_aggregate_radius_purity30')

adata = sm.tl.spatial_aggregate(adata, x_coordinate= "XMin",
                                y_coordinate= "YMin",
                                phenotype='phenotype', method='radius', radius=50, purity=20,
                                imageid='imageid', subset=None, label='spatial_aggregate_radius_purity20')

adata = sm.tl.spatial_aggregate(adata, x_coordinate= "XMin",
                                y_coordinate= "YMin",
                                phenotype='phenotype', method='radius', radius=50, purity=10,
                                imageid='imageid', subset=None, label='spatial_aggregate_radius_purity10')



totalnumber = adata.obs['phenotype'].value_counts()
nonsig90 = adata[adata.obs.spatial_aggregate_radius_purity90 == "non-significant"].obs['phenotype'].value_counts()
nonsig80 = adata[adata.obs.spatial_aggregate_radius_purity80 == "non-significant"].obs['phenotype'].value_counts()
nonsig70 = adata[adata.obs.spatial_aggregate_radius_purity70 == "non-significant"].obs['phenotype'].value_counts()
nonsig60 = adata[adata.obs.spatial_aggregate_radius_purity60 == "non-significant"].obs['phenotype'].value_counts()
nonsig50 = adata[adata.obs.spatial_aggregate_radius_purity50 == "non-significant"].obs['phenotype'].value_counts()
nonsig40 = adata[adata.obs.spatial_aggregate_radius_purity40 == "non-significant"].obs['phenotype'].value_counts()
nonsig30 = adata[adata.obs.spatial_aggregate_radius_purity30 == "non-significant"].obs['phenotype'].value_counts()
nonsig20 = adata[adata.obs.spatial_aggregate_radius_purity20 == "non-significant"].obs['phenotype'].value_counts()
nonsig10 = adata[adata.obs.spatial_aggregate_radius_purity10 == "non-significant"].obs['phenotype'].value_counts()


totalnumber.name = 'totalnumber'
nonsig90.name = 'nonsig90'
nonsig80.name = 'nonsig80'
nonsig70.name = 'nonsig70'
nonsig60.name = 'nonsig60'
nonsig50.name = 'nonsig50'
nonsig40.name = 'nonsig40'
nonsig30.name = 'nonsig30'
nonsig20.name = 'nonsig20'
nonsig10.name = 'nonsig10'


dfPurity = pd.DataFrame([totalnumber, nonsig90,nonsig80, nonsig70, nonsig60, nonsig50, nonsig40, nonsig30, nonsig20, nonsig10])

dfPurity.to_csv("dfPurity.csv", index=False)
