#!/usr/bin/env python

# Read the results of Alevin and write outputs to a .mtx file readable by tools
# expecting 10X outputs. Adapted from 
# https://github.com/ebi-gene-expression-group/scxa-droplet-quantification-workflow/blob/develop/bin/alevinMtxTo10x.py
# which is adapted from 
# https://github.com/k3yavi/vpolo/blob/master/vpolo/alevin/parser.py 

from __future__ import print_function
from collections import defaultdict
from struct import Struct
import pandas as pd
import gzip
import sys
import os
from scipy.io import mmread,mmwrite
from scipy.sparse import *
from shutil import copyfile
import pathlib
import numpy as np
import argparse
import pyroe

parser = argparse.ArgumentParser(description='Convert Alevin outputs to 10X .mtx.')
parser.add_argument('alevin_fry_quant', help = 'Alevin output directory')
parser.add_argument('mtx_out', help = 'Output directory for converted results')
parser.add_argument('--cell_prefix', dest='cell_prefix', default='', help = 'Prefix to apply to cell barcodes')
parser.add_argument('mode', default='scRNA', help='scRNA or snRNA')
args = parser.parse_args() 

alevin_out=args.alevin_fry_quant
mtx_out=args.mtx_out
cell_prefix=args.cell_prefix
mode=args.mode

# Run some checks in the Alevin output

if not os.path.isdir(alevin_out):
    print("{} is not a directory".format( alevin_out ))
    sys.exit(1)

# Read mtx from alevin_fry_quant 
ad = pyroe.load_fry(alevin_out, output_format=mode)

#pd.DataFrame(ad.var.index).to_csv(os.path.join(destination, "genes.tsv" ),   sep = "\t", index_col = False)
#pd.DataFrame(ad.obs.index).to_csv(os.path.join(destination, "barcodes.tsv"), sep = "\t", index_col = False)
# ad.obs.to_csv(os.path.join(destination, "metadata.tsv"), sep = "\t", index_col = True
#scipy.io.mmwrite(os.path.join(destination, "matrix.mtx"), ad.X)

cb_names = [cell_prefix + s for s in ad.obs_names]
gene_names = ad.var_names
umi_counts = ad.X
    
# Write outputs to a .mtx file readable by tools expecting 10X outputs.
# Barcodes file works as-is, genes need to be two-column, duplicating the
# identifiers. Matrix itself needs to have genes by row, so we transpose. 

pathlib.Path(mtx_out).mkdir(parents=True, exist_ok=True)
mmwrite('%s/matrix.mtx' % mtx_out, umi_counts.transpose(), field='real') 

genes_frame = pd.DataFrame(zip(gene_names, gene_names))
genes_frame.to_csv(path_or_buf='%s/genes.tsv' % mtx_out, index=False, sep="\t", header = False)

with open('%s/barcodes.tsv' % mtx_out, 'w') as f:
    f.write("\n".join(cb_names))    
    f.write("\n")    
