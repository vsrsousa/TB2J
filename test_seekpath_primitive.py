#!/usr/bin/env python
"""
Test SeekPath with primitive cell (to match ASE comparison)
"""
import numpy as np
from ase.io import read
from ase.dft.kpoints import bandpath as ase_bandpath
import seekpath

# Read Gd structure
atoms = read('/Users/vinicius/Documents/cifs/Gd.cif')
cell = atoms.get_cell()
positions_frac = atoms.get_scaled_positions()
atomic_numbers = atoms.get_atomic_numbers()

print("=" * 70)
print("ORIGINAL Gd structure (from CIF - is already PRIMITIVE):")
print("=" * 70)
print(f"Cell:\n{cell}")
print(f"Scaled positions:\n{positions_frac}")
print(f"Atomic numbers: {atomic_numbers}")

# ============================================================================
print("\n" + "=" * 70)
print("ASE.bandpath (with original/primitive cell):")
print("=" * 70)

result_ase = ase_bandpath(None, cell, npoints=50)
xs, Xs, knames = result_ase.get_linear_kpoint_axis()
print(f"Path: {result_ase.path}")
print(f"knames: {knames}")
print(f"kpts shape: {result_ase.kpts.shape}")
print(f"Special points: {list(result_ase.special_points.keys())}")

# ============================================================================
print("\n" + "=" * 70)
print("SeekPath.get_path output:")
print("=" * 70)

result_seekpath = seekpath.get_path(
    (cell.tolist(), positions_frac.tolist(), atomic_numbers.tolist())
)

print(f"\nPrimitive cell (from SeekPath):")
prim_cell = np.array(result_seekpath['primitive_lattice'])
print(prim_cell)

print(f"\nPrimitive positions:")
prim_positions = np.array(result_seekpath['primitive_positions'])
print(prim_positions)

print(f"\nPrimitive types: {result_seekpath['primitive_types']}")

# ============================================================================
print("\n" + "=" * 70)
print("SeekPath on PRIMITIVE cell:")
print("=" * 70)

result_seekpath_prim = seekpath.get_path(
    (prim_cell.tolist(), prim_positions.tolist(), result_seekpath['primitive_types'])
)

print(f"Path segments: {result_seekpath_prim['path']}")
print(f"Special points: {list(result_seekpath_prim['point_coords'].keys())}")
print(f"Bravais lattice: {result_seekpath_prim['bravais_lattice']}")

# ============================================================================
print("\n" + "=" * 70)
print("COMPARISON:")
print("=" * 70)

print("\nASE (primitive cell):")
print(f"  Path: {result_ase.path}")
print(f"  Points: {list(result_ase.special_points.keys())}")

print("\nSeekPath (original/CIF cell):")
print(f"  Path: {result_seekpath['path']}")
print(f"  Points: {list(result_seekpath['point_coords'].keys())}")

print("\nSeekPath (primitive cell - extracted from above):")
print(f"  Path: {result_seekpath_prim['path']}")
print(f"  Points: {list(result_seekpath_prim['point_coords'].keys())}")

print("\nCells used:")
print(f"  ASE cell:\n{result_ase.cell}")
print(f"  SeekPath primitive:\n{prim_cell}")
print(f"  Match: {np.allclose(result_ase.cell, prim_cell)}")
