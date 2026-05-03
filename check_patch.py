#!/usr/bin/env python
"""
Check if patch was applied to Cell.bandpath
"""
from ase.cell import Cell

# Check the actual function
print(f"Cell.bandpath: {Cell.bandpath}")
print(f"Is wrapper? {'cell_bandpath_wrapper' in str(Cell.bandpath)}")

# Try to import patch and check
import TB2J.seekpath_patch as sp
print(f"\nAfter import TB2J.seekpath_patch:")
print(f"Cell.bandpath: {Cell.bandpath}")
print(f"Is wrapper? {'cell_bandpath_wrapper' in str(Cell.bandpath)}")

# Check what the function actually is
import inspect
sig = inspect.signature(Cell.bandpath)
print(f"\nSignature: {sig}")
