import urllib
import StringIO
import gzip
import pandas as pd
import scipy.stats
import numpy as np

import sys
import io
sys.path.append("/PATH/TO/l1ktools/python")
import cmap.io.gct as gct
import cmap.io.plategrp as grp

gctx_path = '/PATH/TO/lincs/modzs.gctx'


# function that takes gzipped github tsv files and creates dataframe
def url_to_df(path):
    url = urllib.urlopen(path)
    url_f = StringIO.StringIO(url.read())
    g = gzip.GzipFile(fileobj=url_f)
    return pd.read_table(g)

def extract_from_gctx(path, probes, signatures):
    """Returns a DataFrame with probes as rows and signatures as columns."""
    gct_object = gct.GCT(path)
    gct_object.read_gctx_matrix(cid = signatures, rid = probes)
    return pd.DataFrame(gct_object.matrix, index=probe_ids, columns=sig_ids)

def probes_to_genes(df, probe_to_gene):
    """Converts probe level dataframe to gene level dataframe"""
    get_gene = lambda probe: probe_to_gene.get(probe)
    grouped = df.groupby(by=get_gene, axis=0)
    gene_df = grouped.mean()
    return gene_df

def get_consensus_signatures(df, pert_to_sigs):
    """pert_to_sigs is a dictionary of context_id to a list of signatures."""
    consensuses = dict()
    for pert, sigs in pert_to_sigs.items():
        consensuses[pert] = get_consensus_signature(df[sigs])
    return pd.DataFrame(consensuses)

def get_consensus_signature(df):
    """TODO"""
    weights = weight_signature(df)
    consensus = df.apply(stouffer, axis=1, weights=weights)
    return consensus

def stouffer(z_scores, weights):
    assert len(z_scores) == len(weights)
    z_scores = np.array(z_scores)
    weights = np.array(weights)
    return np.sum(z_scores * weights) / np.sqrt(np.sum(weights ** 2))

def weight_signature(df, min_cor = 0.05):
    """
    """
    if len(df.columns) == 1:
        return np.array([1])
    
    if len(df.columns) == 2:
        return np.array([0.5, 0.5])
    
    rho, p = scipy.stats.spearmanr(df, axis=0)
    mean_cor = (np.sum(rho, axis=0) - 1) / (len(rho) - 1)
    weights = np.maximum(mean_cor, min_cor)
    weights /= np.sum(weights)
    return weights

""" TESTS
probes = ['218075_at', '218434_s_at', '202852_s_at', '201511_at']
sigs = ['AML001_CD34_24H:BRD-A19037878:1.11111', 'AML001_CD34_24H:BRD-A19037878:3.33333', 'AML001_CD34_24H:BRD-A19500257:0.37037']
df = extract_from_gctx(gctx_path, probes, sigs)
df

probes_to_genes(df, {'218075_at': 1, '218434_s_at': 2, '202852_s_at': 1, '201511_at': 2})

weight_signature(df)

get_consensus(df)

context_to_sigs = dict()
for sig in df.columns:
    context = sig.split(':')[1]
    context_to_sigs.setdefault(context, []).append(sig)
signatures_to_context(df, context_to_sigs)
"""
