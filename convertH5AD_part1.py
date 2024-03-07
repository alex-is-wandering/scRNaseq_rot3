#in Python:
import scipy.io
import scanpy as sc

adata=sc.read("/workdir/data/crosstissue_data/human/gut/gutcellatlas/Full_obj_log_counts_soupx_v2.h5ad") #"FutureWarning: `anndata.read` is deprecated, use `anndata.read_h5ad` instead. `ad.read` will be removed in mid 2024"

scipy.io.mmwrite("/workdir/data/crosstissue_data/human/gut/gutcellatlas/gutatlas_processed_v.mtx",adata.X) #keyerror: 'winsorized' - doesn't show up as an obs or var; make it? substitute for something else? if the data hasn't been processed yet (raw) is that a part of it?
adata.var.to_csv("/workdir/data/crosstissue_data/human/gut/gutcellatlas/gutatlas_processed_var.csv")
adata.obs.to_csv("/workdir/data/crosstissue_data/human/gut/gutcellatlas/gutatlas_processed_obs.csv")

#doesn't work to get these files, come back to later maybe
#import pandas as pd
#import numpy as np
#uns_df = pd.DataFrame(adata.uns)
#uns_dict = adata.uns
#array_lengths = {key: len(value) for key, value in uns_dict.items() if isinstance(value, (list, np.ndarray))}
#if len(set(array_lengths.values())) > 1: raise ValueError("Arrays in 'uns' must be of the same length.") #yep
#uns_df.to_csv("/workdir/data/crosstissue_data/human/gut/gutcellatlas/gutatlas_processed_uns.csv",index=True)
#adata.uns.to_csv("/workdir/data/crosstissue_data/human/gut/gutcellatlas/gutatlas_processed_uns.csv")
#for key, matrix in adata.obsm.items(): pd.DataFrame(matrix).to_csv(f"/workdir/data/crosstissue_data/human/gut/gutcellatlas/gutatlas_processed_obsm_{key}.csv",index=True)
#adata.obsm.to_csv("/workdir/data/crosstissue_data/human/gut/gutcellatlas/gutatlas_processed_obsm_{key}.csv")


#for R:
#library("Seurat")
#library("SeuratObject")
#library("Matrix")
#obj=readMM("gutatlas_processed_v.mtx")
#obj=t(obj)
#obs=read.csv("gutatlas_processed_obs.csv")
#var=read.csv("gutatlas_processed_var.csv")
#uns=read.csv("gutatlas_processed_uns.csv")
#obsm=read.csv("gutatlas_processed_obsm.csv")
#rownames(obj)=var$featurekey
#colnames(obj)=obs$barcodes
#gutatlas_processed=CreateSeuratObject(obj,project="gutatlas_processed")

