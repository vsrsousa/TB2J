#!/usr/bin/env python
"""
Test different ways to call ASE.bandpath
"""
import numpy as np
from ase.io import read
from ase.cell import Cell
from ase.dft.kpoints import bandpath as ase_bandpath

# Read Gd structure
atoms = read('/Users/vinicius/Documents/cifs/Gd.cif')
cell = atoms.get_cell()

print("=" * 70)
print("Testing different ASE.bandpath calls:")
print("=" * 70)

# Test 1: with empty list
print("\n1. ase_bandpath([], cell, npoints=50):")
try:
    result = ase_bandpath([], cell, npoints=50)
    xs, Xs, knames = result.get_linear_kpoint_axis()
    print(f"   Path: {result.path}")
    print(f"   knames: {knames}")
    print(f"   kpts shape: {result.kpts.shape}")
except Exception as e:
    print(f"   ERROR: {e}")

# Test 2: with None
print("\n2. ase_bandpath(None, cell, npoints=50):")
try:
    result = ase_bandpath(None, cell, npoints=50)
    xs, Xs, knames = result.get_linear_kpoint_axis()
    print(f"   Path: {result.path}")
    print(f"   knames: {knames}")
    print(f"   kpts shape: {result.kpts.shape}")
except Exception as e:
    print(f"   ERROR: {e}")

# Test 3: with string 'GXM'
print("\n3. ase_bandpath('GXM', cell, npoints=50):")
try:
    result = ase_bandpath('GXM', cell, npoints=50)
    xs, Xs, knames = result.get_linear_kpoint_axis()
    print(f"   Path: {result.path}")
    print(f"   knames: {knames}")
    print(f"   kpts shape: {result.kpts.shape}")
except Exception as e:
    print(f"   ERROR: {e}")

# Test 4: using Cell.bandpath() with no args
print("\n4. Cell(cell).bandpath(npoints=50):")
try:
    result = Cell(cell).bandpath(npoints=50)
    xs, Xs, knames = result.get_linear_kpoint_axis()
    print(f"   Path: {result.path if hasattr(result, 'path') else 'N/A'}")
    print(f"   knames: {knames}")
    print(f"   kpts shape: {result.kpts.shape}")
except Exception as e:
    print(f"   ERROR: {e}")

# Test 5: Cell.bandpath() with string path
print("\n5. Cell(cell).bandpath('GXM', npoints=50):")
try:
    result = Cell(cell).bandpath('GXM', npoints=50)
    xs, Xs, knames = result.get_linear_kpoint_axis()
    print(f"   Path: {result.path if hasattr(result, 'path') else 'N/A'}")
    print(f"   knames: {knames}")
    print(f"   kpts shape: {result.kpts.shape}")
except Exception as e:
    print(f"   ERROR: {e}")

# Test 6: Check what bandpath signature is
print("\n6. ASE.bandpath signature:")
import inspect
sig = inspect.signature(ase_bandpath)
print(f"   {sig}")

# Test 7: Try with empty string
print("\n7. ase_bandpath('', cell, npoints=50):")
try:
    result = ase_bandpath('', cell, npoints=50)
    xs, Xs, knames = result.get_linear_kpoint_axis()
    print(f"   Path: {result.path}")
    print(f"   knames: {knames}")
    print(f"   kpts shape: {result.kpts.shape}")
except Exception as e:
    print(f"   ERROR: {e}")
