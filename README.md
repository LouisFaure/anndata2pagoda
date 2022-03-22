# anndata2pagoda

Small pagoda2 web app generator from anndata object.


## Installation

```
mamba create -n anndata2pagoda -c bioconda rpy2=3.4.2 python=3.8 bioconductor-s4vectors bioconductor-singlecellexperiment r-base64enc r-stringr
conda activate anndata2pagoda
Rscript -e 'install.packages("pagoda2",repo="https://cran.wu.ac.at/")'
pip install git+https://github.com/LouisFaure/anndata2pagoda.git
```

## Usage

### In python
```python
import scanpy as sc
from anndata2pagoda import pagoda2web

adata = sc.read("adata.h5ad")
pagoda2web(adata)
```

### In command line
```bash
anndata2pagoda -a adata.h5ad
```
