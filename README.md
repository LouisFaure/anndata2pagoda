# anndata2pagoda

Small pagoda2 web app generator from anndata object.


## Installation

```
pip install git+https://github.com/LouisFaure/anndata2pagoda.git
```

## Usage

### In python
```python
import scanpy as sc
sc.read("adata.h5ad")
from anndata2pagoda import pagoda2web
pagoda2web(adata)
```

### In command line
```bash
anndata2pagoda -a adata.h5ad
```