from .equations import equations

try:
    from ._version import version as __version__
except ModuleNotFoundError:
    __version__ = ""