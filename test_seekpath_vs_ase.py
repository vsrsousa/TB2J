#!/usr/bin/env python
"""
Compare raw SeekPath vs ASE bandpath (without TB2J patch)
"""
import numpy as np
from ase.io import read
from ase.cell import Cell
from ase.dft.kpoints import bandpath as ase_bandpath
import seekpath

# Read Gd structure from CIF
atoms = read('/Users/vinicius/Documents/cifs/Gd.cif')
cell = atoms.get_cell()
positions_frac = atoms.get_scaled_positions()
positions_cart = atoms.get_positions()
atomic_numbers = atoms.get_atomic_numbers()

print("=" * 70)
print("Structure: Gd (HCP)")
print("=" * 70)
print(f"Cell:\n{cell}")
print(f"Atomic numbers: {atomic_numbers}")
print(f"Fractional positions:\n{positions_frac}")

# ============================================================================
print("\n" + "=" * 70)
print("ASE.bandpath (original, no patch):")
print("=" * 70)

try:
    result_ase = ase_bandpath(None, cell, npoints=50)
    xs_ase, Xs_ase, knames_ase = result_ase.get_linear_kpoint_axis()
    print(f"Special points: {list(result_ase.special_points.keys())}")
    print(f"Path: {result_ase.path}")
    print(f"knames: {knames_ase}")
    print(f"Xs: {Xs_ase}")
    print(f"kpts shape: {result_ase.kpts.shape}")
except Exception as e:
    print(f"ERROR: {e}")
    import traceback
    traceback.print_exc()

# ============================================================================
print("\n" + "=" * 70)
print("SeekPath.get_path (raw, no ASE):")
print("=" * 70)

try:
    result_seekpath = seekpath.get_path(
        (cell.tolist(), positions_frac.tolist(), atomic_numbers.tolist())
    )
    print(f"Special points: {list(result_seekpath['point_coords'].keys())}")
    print(f"Path segments: {result_seekpath['path']}")
    
    # Extract coordinates
    print(f"\nPoint coordinates:")
    for name, coord in result_seekpath['point_coords'].items():
        print(f"  {name}: {coord}")
    
except Exception as e:
    print(f"ERROR: {e}")
    import traceback
    traceback.print_exc()

# ============================================================================
print("\n" + "=" * 70)
print("SeekPath.get_path (with Cartesian positions):")
print("=" * 70)

try:
    result_seekpath_cart = seekpath.get_path(
        (cell.tolist(), positions_cart.tolist(), atomic_numbers.tolist())
    )
    print(f"Special points: {list(result_seekpath_cart['point_coords'].keys())}")
    print(f"Path segments: {result_seekpath_cart['path']}")
    
except Exception as e:
    print(f"ERROR: {e}")
