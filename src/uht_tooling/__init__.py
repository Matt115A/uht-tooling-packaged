"""Core package initialization for uht_tooling."""

from importlib.metadata import PackageNotFoundError, version

__all__ = ["__version__"]

try:
    __version__ = version("uht-tooling")
except PackageNotFoundError:
    __version__ = "0.0.0"
