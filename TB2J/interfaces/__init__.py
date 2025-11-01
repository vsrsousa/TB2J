# Lazy imports to avoid requiring all optional dependencies
# This allows wann2J.py to work without HamiltonIO (needed by abacus/siesta/lawaf)

from .manager import Manager
from .wannier90_interface import WannierManager, gen_exchange
from .paoflow_interface import PAOflowManager, gen_exchange_paoflow

__all__ = [
    "Manager",
    "gen_exchange_siesta",
    "WannierManager",
    "gen_exchange",
    "gen_exchange_abacus",
    "PAOflowManager",
    "gen_exchange_paoflow",
]

def __getattr__(name):
    """Lazy import for optional interfaces that require HamiltonIO."""
    if name == "gen_exchange_abacus":
        from .abacus import gen_exchange_abacus
        return gen_exchange_abacus
    elif name == "gen_exchange_siesta":
        from .siesta_interface import gen_exchange_siesta
        return gen_exchange_siesta
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
