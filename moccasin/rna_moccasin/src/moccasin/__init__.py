try:
    from ._version import version as __version__
except (ModuleNotFoundError, ImportError):
    __version__ = "0.26.0unknown"
except Exception:
    raise
