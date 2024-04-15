import os
os.environ["OMP_NUM_THREADS"] = "32"
from graph_tool.all import *
import pandas as pd
import numpy as np
import scipy as sp
from sklearn.covariance import LedoitWolf, OAS
import matplotlib.pyplot as plt
import matplotlib.cm as mpl
import seaborn as sns
import statsmodels.api as sm
from multipy.fdr import qvalue
from multipy.fdr import lsu

import dill
import sys

#  read command line argument for current tissue variable
current_tissue = sys.argv[1]
print('Current tissue is ' + current_tissue)

g_path = '../../SBM/snakemake/cache/trimmed_graph/fdr-1e-3/layered/'
tissues = ['head', 'body']
conditions = ['hs', 'ctrl']

print('Loading graphs...')
graphs = {f'{tissue}':load_graph(g_path + f'{tissue}.xml.gz') for tissue in tissues}
for t in tissues:
    remove_parallel_edges(graphs[t])

print('Reading gene expression data...')
data_dir = '/Genomics/argo/users/damelo/projects/HS-Expression-GxE/SBM/rawData/layered/' + current_tissue
input_list = []
for file in os.listdir(data_dir):
    if file.endswith(".tsv"):
        input_list.append(file)
input_names = list(map(lambda p: p[:p.rfind('.')], input_list))
gene_expr = []
for file in input_list:
    gene_expr_raw = pd.read_table(os.path.join(data_dir, file))
    gene_expr.append(gene_expr_raw.T)

gene_expr_concat = pd.concat([gene_expr[0], gene_expr[1]], 
                   axis=0, 
                   keys=input_names,
                   names=['source']).reset_index(level=[0])

# Extract source column
source = gene_expr_concat['source']
ge = gene_expr_concat.drop(columns=['source'])

print('Normalizing gene expression data...')
z_scores = (ge - ge.mean()) / np.sqrt(ge.var())

# Model matrix
X = pd.get_dummies(source, drop_first=True)
X = np.column_stack((np.ones(X.shape[0]), X))

g = graphs[current_tissue]

genes = g.vp.genes

clip_a = g.new_ep("double", 0)
clip_b = g.new_ep("double", 0)
clip_p = g.new_ep("double", 0)
clip_shift = g.new_ep("double", 0)

print('Computing edge coefficients...')
for e in g.edges():
    v_i, v_j = e.source(), e.target()
    gene_i = genes[v_i]
    gene_j = genes[v_j]
    y = z_scores[gene_i] * z_scores[gene_j]
    mod = sm.OLS(y, X)
    fii = mod.fit()
    intercept = fii.summary2().tables[1].iloc[0,0]
    coef = fii.summary2().tables[1].iloc[1,0]
    pvalue = fii.summary2().tables[1].iloc[1,3]
    zero_shift = abs(intercept + coef) - abs(intercept)

    clip_a[e] = intercept
    clip_b[e] = coef
    clip_p[e] = pvalue
    clip_shift[e] = zero_shift

print('Creating edge properties...')
g.edge_properties["clip_a"] = clip_a
g.edge_properties["clip_b"] = clip_b
g.edge_properties["clip_p"] = clip_p
g.edge_properties["clip_shift"] = clip_shift

print('Saving graph...')
g.save('clip_g_' + current_tissue +  '.xml.gz')

print('Done!')