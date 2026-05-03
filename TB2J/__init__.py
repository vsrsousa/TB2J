import importlib.metadata

# Apply seekpath patch early for k-path generation using SeekPath
try:
    from . import seekpath_patch  # type: ignore
except Exception:
    pass

__version__ = importlib.metadata.version("TB2J")
