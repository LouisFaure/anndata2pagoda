from . import __path__, __version__
import os
import anndata
version = __version__
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
parser.add_argument("--use_rep","-r", help="basis of 2D embedding to use (default: UMAP).", default="UMAP")
parser.add_argument("--key_to_include","-k", help="key to include (for more than one, add space) (default: leiden).",
                    nargs='+',default=["leiden"])
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
             layer=None,
             use_rep="UMAP",
             key_to_include=["leiden"],
             title="p2w",
             filename="p2w.bin"):
    
    layer_name="X" if layer is None else layer
    print("Converting to pagoda2 web object from anndata: layer %s, %s embedding, %s obs keys" %(layer_name,
                                                                                                 use_rep,
                                                                                                 key_to_include))
    
    adata = adata.copy()
    
    try:
        adata.obs.leiden.cat.categories.astype(int)
        leiden = "leiden"
    except TypeError:
        adata.obs["leiden_noname"]=adata.obs.leiden.cat.rename_categories(range(len(adata.obs.leiden.cat.categories)))
        leiden = "leiden_noname"
    
    adata.obs[leiden]=adata.obs[leiden].cat.rename_categories((adata.obs[leiden].cat.categories.astype(int)+1).astype(str))

    
    palettes = {}
    for key in key_to_include:
        palettes[key] = adata.uns[key+"_colors"].tolist()
        
    if layer is not None:
        adata.X = adata.layers[layer].copy()
    
    adata.layers=None
    print("    converting to SingleCellExperiment")
    adata = py2rpy(adata)
        
    palettes = _dict_to_rlist(palettes)
    
    to_p2w(adata,use_rep,palettes,title,filename)
    
    print("    done! Saved as %s" %filename)
    

def main():
    args = parser.parse_args()
    adata=anndata.read_h5ad(args.adata)
    pagoda2web(adata,args.layer,args.use_rep,args.key_to_include,args.title,args.filename)
