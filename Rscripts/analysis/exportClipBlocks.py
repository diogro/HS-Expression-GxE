import os
os.environ["OMP_NUM_THREADS"] = '16' # export OMP_NUM_THREADS=4
os.environ["OPENBLAS_NUM_THREADS"] = '16' # export OPENBLAS_NUM_THREADS=4
os.environ["NUMEXPR_NUM_THREADS"] = '16' # export NUMEXPR_NUM_THREADS=6

import dill

from graph_tool.all import *
import numpy as np
import pandas as pd
import scipy as sp
from sklearn.covariance import LedoitWolf, OAS
import matplotlib.pyplot as plt
import matplotlib.cm as mpl
import seaborn as sns
import statsmodels.api as sm
from multipy.fdr import qvalue
from multipy.fdr import lsu
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

def filterBySign(g, ep, positive = True, keepOnlyMain = False):
    # Filtering edges
    corr = g.edge_properties[ep]
    sign = g.new_ep("bool", True)
    if positive:
        sign.a = np.array(corr.a > 0)
    else:
        sign.a = np.array(corr.a < 0)    

    tv = GraphView(g, efilt=sign)

    # Keeping largest component
    if keepOnlyMain:
        comp, hist = label_components(tv)
        main_component = tv.new_vp("bool", (comp.a == np.where(hist == max(hist))[0][0]))
        tv.vertex_properties["main_component"] = main_component
        tv.set_vertex_filter(main_component)
    return tv
    
def first_occurrence(x):
    _, idx = np.unique(x, return_index=True)
    return idx

def labelVertices(state):
    g = state.g
    g.vp.level_0 = g.new_vertex_property("double", state.get_bs()[0])
    first = first_occurrence(np.array([g.vp.level_0.a]))
    state.g.vp.labels = g.new_vp("string", [str(int(x)) if i in first else "" for i, x in enumerate(g.vp.level_0.a)])
    return state

def ensure_dir(file_path):
    directory = os.path.dirname(file_path)
    if not os.path.exists(directory):
        os.makedirs(directory)
    if not os.path.exists(file_path):
        os.makedirs(file_path)

def create_nestedBlock_df(g, corr, state):
    genes = g.vertex_properties["genes"]
    nested_block_df = pd.DataFrame(columns=('Gene', "Degree", "Average", 'Sum', 'B1', "B2", "B3", "B4", "B5", "B6"))
    for v in g.vertex_index:
        line = [genes[v]]
        line.append(g.get_total_degrees([v])[0])
        line.append(np.mean(g.get_all_edges(v, [corr] )[:,2]))
        line.append(np.sum(g.get_all_edges(v, [corr] )[:,2]))
        [line.append(i) for i in get_group(v, state)]
        nested_block_df.loc[v] = line
    nested_block_df = nested_block_df[nested_block_df.Degree > 0]
    # Sort by Degree, Sum and Average
    nested_block_df = nested_block_df.sort_values(by=['Degree', 'Sum', 'Average'], ascending=False)
    return nested_block_df

def get_group(x, state):
    levels = state.get_levels()
    n_levels = len(levels)
    r = np.zeros(n_levels)
    r[0] = levels[0].get_blocks()[x]
    for i in range(1, n_levels):
        r[i] = levels[i].get_blocks()[r[i-1]]
    r = r.astype(int)
    return r

def makeClipGraph (current_tissue, clip_fdr):

    # output path
    output_path = "../../tmp/clip/" + "fdr-" + str(clip_fdr) + "/" + current_tissue + "-nEdges_" + str(n_edges) 
    # create output directory
    os.makedirs(os.path.dirname(output_path), exist_ok=True)

    g_clip = clip_graphs[current_tissue]
    gNeg= filterByFDR(g_clip, clip_fdr, 'clip_p')
    gNeg = Graph(gNeg, prune=True)
    gNeg = filterBySign(gNeg, 'clip_shift', positive = False)
    gNeg = Graph(gNeg, prune=True)
    gNeg.ep.clip_shift.a = np.abs(gNeg.ep.clip_shift.a)

    bs = blocks[current_tissue]

    s_Neg = NestedBlockState(gNeg, bs=bs,
                             state_args=dict(recs=[gNeg.ep.clip_shift],
                                             rec_types=["real-normal"]))

    g_clip = clip_graphs[current_tissue]
    gPos = filterByFDR(g_clip, clip_fdr, 'clip_p')
    gPos = Graph(gPos, prune=True)
    gPos = filterBySign(gPos, 'clip_shift', positive = True)
    gPos = Graph(gPos, prune=True)

    bs = blocks[current_tissue]
    s_Pos = NestedBlockState(gPos, bs=bs,
                            state_args=dict(recs=[gPos.ep.clip_shift],
                                            rec_types=["real-normal"]))

    return {"decohere": s_Neg, "integrate": s_Pos}


if __name__ == '__main__':

    tissues = ['head', 'body']
    conditions = ['hs', 'ctrl']

    b_path = '../../SBM/snakemake/cache/MCMC/blocks/fdr-1e-3/layered/'
    print("Reading blocks...")
    blocks = {f'{tissue}':load_blocks(b_path + f'{tissue}.dill') for tissue in tissues}

    # Read clip graphs
    clip_path = '../../cache/'
    print("Reading clip graphs...")
    clip_graphs = {f'{tissue}':load_graph(clip_path + f'clip_g_{tissue}.xml.gz') for tissue in tissues}

    # Make clip graphs for various FDRs
    print("Making clip block states for various FDRs...")
    fdrs = [1e-1, 1e-2, 1e-3, 1e-4]
    n_edges = 50000
    clip_state = {}
    for tissue in tissues:
        print(f"Tissue: {tissue}")
        clip_state[tissue] = {}
        for fdr in fdrs:
            print(f"Making clip graphs for FDR {fdr}")
            l = tissue + '-ctrl'
            clip_state[tissue][fdr] = makeClipGraph(tissue, fdr)

    print("Creating nested block DataFrame...")
    block_dict = {}
    for t in tissues:
        decohere_df = create_nestedBlock_df(clip_state[t][1e-2]['decohere'].g, 
                                            clip_state[t][1e-2]['decohere'].g.ep.clip_shift, 
                                            clip_state[t][1e-2]['decohere'])
        integrate_df = create_nestedBlock_df(clip_state[t][1e-2]['integrate'].g, 
                                            clip_state[t][1e-2]['integrate'].g.ep.clip_shift, 
                                            clip_state[t][1e-2]['integrate'])
        merge = pd.merge(decohere_df, integrate_df, on='Gene', suffixes=('_decohere', '_integrate'))
        clip_Df = pd.concat([decohere_df, integrate_df],
                        axis=0, 
                        keys=['decohere', 'integrate'],
                        names=['direction']).reset_index(level=[0])
        # Write df to file
        clip_Df.to_csv("../../cache/" + f'clip_fdr-1e2_{t}.csv', sep="\t")
        block_dict[t] = clip_Df

    print("Creating block summary...")
    directions = ['decohere', 'integrate']
    for t in tissues:
        print("Creating block DataFrame...")
        block_df = block_dict[t]

        print("Calculating block sizes...")
        blocks = [list(set(block_df[b])) for b in block_df.filter(like='B', axis=1)]
        block_sizes = [len(b) for b in blocks]
        block_sizes = -np.sort(-np.array(list(set(block_sizes))))
        block_sizes = [x for x in block_sizes if x >= 2]

        print("Creating gene lists...")
        n_levels = len(block_sizes)
        output_df = pd.DataFrame(columns=('Nested_Level', 'Block', 'Path', 'Direction', 'N_genes', 'Internal_degree', 'Average_internal_degree', 'Total_degree', 'Average_total_degree' , 'Assortativity'))
        l = 0
        for d in directions:
            print(f"Processing {t} - {d}")
            state = clip_state[t][1e-2][d]
            g = state.g
            corr = g.ep.clip_shift


            for i in range(n_levels):
                print("At level: " + str(i+1))
                bl = blocks[i]
                for b in bl:

                    line = [i+1]
                    line.append(b)


                    df = block_df[block_df['B' + str(i+1)]==b]
                    genes = df["Gene"]
                               
                    block_path = '-'.join([str(num) for num in list(df.filter(like='B', axis=1).iloc[0,range(i, n_levels)])]) 
                    line.append(block_path)

                    line.append(d)

                    N_genes = genes.shape[0]
                    line.append(N_genes)

                    # Weighted
                    ers = adjacency(state.levels[i].bg, weight=state.levels[i].mrs)
                    B = len(bl)
                    E = ers.sum()
                    q_r = (B/E) * (ers[b,b] - (ers[b,:].sum()**2/E))
                    if np.abs(q_r) > 1:
                        sys.error("q_r is larger than one.")

                    line.append(ers[b,b])
                    line.append(ers[b,b] / N_genes)
                    line.append(ers[b,:].sum())
                    line.append(ers[b,:].sum() / N_genes)
                    line.append(q_r)

                    output_df.loc[l] = line
                    l = l + 1

        print(f"Outputing block summary for {t}...")
        output_df.to_csv("../../cache/" + f'blockSummary_clip_fdr-1e2_{t}.csv', sep="\t")

        print("Done!")
