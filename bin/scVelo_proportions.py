#!/usr/bin/env python

import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
import scvelo as scv
import matplotlib
import scipy
import argparse
from pyroe import load_fry
import os
import sys


andata_list = sys.argv[1:-1]

print(andata_list)

#parser = argparse.ArgumentParser(description='Do velocity analysis')
#parser.add_argument('alevin_fry_quant', help = 'Alevin output directory')
#parser.add_argument('out_file', help = 'output file name for plot')
#args = parser.parse_args() 

#alevin_out=args.alevin_fry_quant
out = sys.argv[-1]


# Run some checks in the Alevin output

# if not os.path.isdir(alevin_out):
#     print("{} is not a directory".format( alevin_out ))
#     sys.exit(1)

# Read mtx from alevin_fry_quant 
adata_appended = []
for andata_path in andata_list:
    andata = load_fry(andata_path, output_format = "velocity")
    adata_appended.append(andata)
adata = ad.concat(adata_appended)

scv.pl.proportions(adata, save=out)
#basic filtering

sc.pp.filter_cells(adata, min_genes=750)
sc.pp.filter_genes(adata, min_cells=3)

sc.pp.normalize_total(adata, target_sum = 1000000)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata)

sc.pp.scale(adata)


# make sure gene names are unique


# get embeddings
sc.tl.pca(adata)
sc.pp.neighbors(adata)
sc.tl.tsne(adata)
sc.tl.umap(adata, n_components = 2)

# housekeeping
matplotlib.use('AGG')
scv.settings.set_figure_params('scvelo')

# get the proportion of spliced and unspliced count
scv.utils.show_proportions(adata)

# filter cells and genes, then normalize expression values
scv.pp.filter_and_normalize(adata, min_shared_counts=20, n_top_genes=2000,enforce=True)

# scVelo pipeline
scv.pp.moments(adata, n_pcs=30, n_neighbors=15)
scv.tl.recover_dynamics(adata, n_jobs = 11)
scv.tl.velocity(adata, mode = 'dynamical')
scv.tl.velocity_graph(adata)
#adata.write('p.h5ad', compression='gzip')
scv.pl.velocity_embedding_stream(adata, basis='X_umap', save=out)

