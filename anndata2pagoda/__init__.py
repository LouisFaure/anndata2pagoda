try:
    from importlib_metadata import version
except:
    from importlib.metadata import version
__version__ = version(__name__)
del version
from .anndata2pagoda import pagoda2web
