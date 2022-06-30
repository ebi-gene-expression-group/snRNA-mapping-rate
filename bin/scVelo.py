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



#basic filtering

sc.pp.filter_cells(adata, min_genes=750)
sc.pp.filter_genes(adata, min_cells=3)

scv.pl.proportions(adata)