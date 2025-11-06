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


########reload dat
adata = ad.read('adata_original_0308.h5ad')
adata.obs_keys()




# Analyze spatial aggregation purity with celltype, purity 50

adata = sm.tl.spatial_aggregate(adata, x_coordinate= "XMin",
                                y_coordinate= "YMin",
                                phenotype='phenotype', method='radius', radius=50, purity=50,
                                imageid='imageid', subset=None, label='spatial_aggregate_radius_purity50')

adata.obs.spatial_aggregate_radius_purity50.value_counts()

adata_nonsig = adata [adata.obs.spatial_aggregate_radius_purity50 == 'non-significant']                 

adata_nonsig.obs.columns
adata_nonsig.obs.Allsubtypes.value_counts()
adata.obs.Allsubtypes.value_counts()


###先purity 50往下看看CN——cell type的情况

sm.tl.spatial_count(adata, 
                    x_coordinate= "XMin",
                    y_coordinate= "YMin",
                    z_coordinate=None, 
                    phenotype='phenotype', 
                    method='radius', radius=50,
                    imageid='imageid', subset=None, 
                    label='spatial_count_CT')

print(adata.uns['spatial_count_CT'])


adata_nonsig = adata [adata.obs.spatial_aggregate_radius_purity50 == 'non-significant']       

cellid = adata_nonsig.obs.index

filter_uns_SCC = adata.uns['spatial_count_CT']

filter_uns_SCC = filter_uns_SCC.loc[cellid]

adata_nonsig.uns['spatial_count_CT'] = filter_uns_SCC

adata_nonsig = sm.tl.spatial_cluster(adata_nonsig, df_name='spatial_count_CT', method='kmeans', k=10, label='cluster_kmeans')  ##nan
adata_nonsig.obs['cluster_kmeans'].value_counts()
print(adata_nonsig.obs['cluster_kmeans'])

#adata_nonsig.write('adata_nonsig_purity50_CN_CT_radius50.h5ad')



df_CN = adata_nonsig.uns['spatial_count_CT']

combinedDf = pd.concat([df_CN, adata_nonsig.obs], axis=1)
combinedDf.columns
combinedDf.cluster_kmeans.value_counts()

combinedDf.to_csv('adata_nonsig_purity50_CN_CT_radius50.csv', index=False)