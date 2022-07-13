#!/usr/bin/env python

import numpy as np
import pandas as pd
import scanpy as sc
import anndata
import scvelo as scv
import matplotlib
import scipy
import argparse
from pyroe import load_fry
import os
import sys

parser = argparse.ArgumentParser(description='Do velocity analysis')
parser.add_argument('alevin_fry_quant', help = 'Alevin output directory')
parser.add_argument('out_file', help = 'output file name for plot')
args = parser.parse_args() 

alevin_out=args.alevin_fry_quant
out = args.out_file


# Run some checks in the Alevin output

if not os.path.isdir(alevin_out):
    print("{} is not a directory".format( alevin_out ))
    sys.exit(1)

# Read mtx from alevin_fry_quant 

adata = load_fry(alevin_out, output_format = "velocity")


scv.pl.proportions(adata, save=out)
#basic filtering

sc.pp.filter_cells(adata, min_genes=750)
sc.pp.filter_genes(adata, min_cells=3)


# make sure gene names are unique
adata.var_names_make_unique()

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

