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

# Loading blocks...
def load_blocks(blocks):
    with open (blocks, "rb") as fh:
        bs = dill.load(fh)[0:6]
    return bs

def filterByFDR(g, level, pval, keepOnlyMain=False):
    # Filtering edges
    pvals = np.array(g.edge_properties[pval].a)

    fdr_ep = g.new_ep("bool", True)
    fdr_ep.a = lsu(pvals, q=level)

    tv = GraphView(g, efilt=fdr_ep)

    # Keeping largest component
    if keepOnlyMain:
        comp, hist = label_components(tv)
        main_component = tv.new_vp("bool", (comp.a == np.where(hist == max(hist))[0][0]))
        tv.vertex_properties["main_component"] = main_component
        tv.set_vertex_filter(main_component)
    return tv

#  read command line argument for current tissue variable
current_tissue = sys.argv[1]
print('Current tissue is ' + current_tissue)

g_path = '../../SBM/snakemake-layer/cache/trimmed_graph/fdr-1e-3/layered/'
tissues = ['head', 'body']
conditions = ['hs', 'ctrl']

print('Loading graphs...')
graphs = {f'{tissue}':load_graph(g_path + f'{tissue}.xml.gz') for tissue in tissues}
b_path = '../../SBM/snakemake-layer/cache/MCMC/blocks/fdr-1e-3/layered/'

print('Loading blocks...')
blocks = {f'{tissue}':load_blocks(b_path + f'{tissue}.dill') for tissue in tissues}

labels = [f'{tissue}-{condition}' for tissue in tissues for condition in conditions]

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
gene_expr_concat
# remove first column
ge = gene_expr_concat.drop(columns=['source'])

print('Normalizing gene expression data...')
z_scores = (ge - ge.mean()) / np.sqrt(ge.var())

# Extract source column
source = gene_expr_concat['source']


# Model matrix
X = pd.get_dummies(source, drop_first=True)
X = np.column_stack((np.ones(X.shape[0]), X))

g = graphs[current_tissue]

genes = g.vp.genes

clip_c = g.new_ep("double", 0)
clip_p = g.new_ep("double", 0)

print('Computing edge coefficients...')
for e in g.edges():
    v_i, v_j = e.source(), e.target()
    gene_i = genes[v_i]
    gene_j = genes[v_j]
    y = z_scores[gene_i] * z_scores[gene_j]
    mod = sm.OLS(y, X)
    fii = mod.fit()
    coef = fii.summary2().tables[1].iloc[1,0]
    pvalue = fii.summary2().tables[1].iloc[1,3]

    clip_c[e] = coef
    clip_p[e] = pvalue

print('Creating edge properties...')
g.edge_properties["clip_c"] = clip_c
g.edge_properties["clip_p"] = clip_p

print('Saving graph...')
g.save('clip_g_head.xml.gz')

print('Done!')