from . import __path__, __version__
version = __version__

from typing import Optional
import os
import anndata
import numpy as np
import shutil

import anndata2ri
from rpy2.robjects import r
anndata2ri.activate()
import rpy2.robjects as ro
from anndata2ri import py2rpy
from rpy2.rinterface_lib import callbacks
callbacks._WRITECONSOLE_EXCEPTION_LOG = "    %s"
from rpy2.rinterface_lib.callbacks import logger as rpy2_logger
import logging
rpy2_logger.setLevel(logging.ERROR)

from rpy2.robjects.packages import importr

if not shutil.which("R"):
    raise Exception(
        "R installation is necessary for converting anndata.\
        \nPlease install R and try again"
    )
   

try:
    pagoda2 = importr("pagoda2")

except Exception as e:
    raise Exception(
        'R package "pagoda2" is necessary for converting anndata!'
    )
    


Rfun = os.path.join(__path__[0], "to_p2w.R")

rsource=ro.r["source"]
rsource(Rfun)
to_p2w = ro.globalenv['to_p2w']


help = 'Generates a pagoda2 web object from a processed anndata object.'

import argparse
from argparse import RawTextHelpFormatter
parser = argparse.ArgumentParser(description=help,formatter_class=RawTextHelpFormatter)
parser.add_argument("--adata","-a", help="Path to anndata file.",required=True)
parser.add_argument("--layer","-l", help="layer from anndata (default: None).", default=None)
parser.add_argument("--use_rep","-r", help="basis of 2D embedding to use (default: X_umap).", default="X_umap")
parser.add_argument("--clustering","-c", help="cluster key to include (default: leiden).",default="leiden")
parser.add_argument("--key_to_include","-k", help="key(s) to include (for more than one, add space) (default: empty).",
                    nargs='+',default=[])
parser.add_argument("--title","-t", help="title of the pagoda2 web object (default: p2w).", default="p2w")
parser.add_argument("--filename","-f", help="filename of pagoda2 web (default: p2w.bin).", default="p2w.bin")


def _dict_to_rlist(dd):
    pairs = []
    for k,v in dd.items():
        if isinstance(v, str):
            pairs.append("{}=\"{}\"".format(k, v))
        else:
            pairs.append("{}={}".format(k, v))
    cmd = "list({})".format(",".join(pairs))
    cmd=cmd.replace("[","c(")
    cmd=cmd.replace("]",")")
    lst = ro.r(cmd)
    return lst


def pagoda2web(adata,
               layer = None,
               use_rep: str = "X_umap",
               clustering: str = "leiden",
               key_to_include: list = [],
               title: str = "p2w",
               filename:str = "p2w.bin"):
    
    layer_name="X" if layer is None else layer
    
    if clustering not in adata.obs:
        raise Exception(f"{clustering} is not present in .obs, available are: {', '.join(adata.obs.columns)}")
        
    if use_rep not in adata.obsm:
        raise Exception(f"{use_rep} is not present in .obsm, available are: {', '.join(list(adata.obsm.keys()))}")
    
    print("Converting to pagoda2 web object from anndata: layer %s, %s embedding, %s obs keys" %(layer_name,
                                                                                                 use_rep,
                                                                                                 key_to_include))
    
    adata = adata.copy()
    
    clusters=clustering
    try:
        adata.obs[clustering].cat.categories.astype(int)
    except TypeError:
        adata.obs[clustering+"_noname"]=adata.obs[clustering].cat.rename_categories(range(len(adata.obs[clustering].cat.categories)))
        clusters = clustering+"_noname"
    
    if adata.obs[clusters].cat.categories.astype(int).min()==0:
        cl = adata.obs[clusters].cat.categories.astype(int)+1
    elif adata.obs[clusters].cat.categories.astype(int).min()==1:
        cl = adata.obs[clusters].cat.categories.astype(int)
    else:
        raise Exception("clustering does not start with zero or one!")
    adata.obs[clusters]=adata.obs[clusters].cat.rename_categories(cl.astype(str))
    
    ks = list(adata.uns.keys())
    for k in ks:
        if (k.find("colors"))==-1:
            del adata.uns[k]
        else:
            kcol = adata.uns[k]
            adata.uns[k] = kcol.tolist() if type(kcol) == np.ndarray else kcol
    
    for k, v in adata.uns.items():
        adata.uns[k] = v.tolist() if type(v) == np.ndarray else v
    
    
    palettes = {}
    for key in key_to_include:
        palettes[key] = adata.uns[key+"_colors"]
        
    if layer is not None:
        adata.X = adata.layers[layer].copy()
    
    adata.layers=None
    print("    converting to SingleCellExperiment")
    adata = py2rpy(adata)
        
    palettes = _dict_to_rlist(palettes)
    
    if use_rep=="X_umap":
        use_rep="UMAP"
    if use_rep=="X_tsne":
        use_rep="TSNE"
    elif use_rep=="X_pca":
        use_rep="PCA"
    
    to_p2w(adata,use_rep,clustering,palettes,title,filename)
    
    print("    done! Saved as %s" %filename)
    

def main():
    args = parser.parse_args()
    adata=anndata.read_h5ad(args.adata)
    pagoda2web(adata,args.layer,args.use_rep,args.clustering,args.key_to_include,args.title,args.filename)
