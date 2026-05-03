#!/usr/bin/env python
"""
Test if spglib standardization helps ASE.bandpath generate paths
"""
import numpy as np
from ase.io import read
from ase.cell import Cell
from ase.dft.kpoints import bandpath as ase_bandpath
import seekpath

try:
    import spglib
except ImportError:
    print("ERROR: spglib not installed")
    exit(1)

# Read Gd structure from CIF
atoms = read('/Users/vinicius/Documents/cifs/Gd.cif')
cell = atoms.get_cell()
positions_frac = atoms.get_scaled_positions()
atomic_numbers = atoms.get_atomic_numbers()

print("=" * 70)
print("ORIGINAL Gd structure (from CIF):")
print("=" * 70)
print(f"Cell:\n{cell}")
print(f"Atomic numbers: {atomic_numbers}")
print(f"Fractional positions:\n{positions_frac}")

# ============================================================================
print("\n" + "=" * 70)
print("ASE.bandpath on ORIGINAL cell:")
print("=" * 70)

try:
    result_ase_orig = ase_bandpath(None, cell, npoints=50)
    xs_ase, Xs_ase, knames_ase = result_ase_orig.get_linear_kpoint_axis()
    print(f"Special points: {list(result_ase_orig.special_points.keys())}")
    print(f"Path: {result_ase_orig.path}")
    print(f"knames: {knames_ase}")
    print(f"kpts shape: {result_ase_orig.kpts.shape}")
except Exception as e:
    print(f"ERROR: {e}")

# ============================================================================
print("\n" + "=" * 70)
print("Standardizing with SPGLIB:")
print("=" * 70)

# Use spglib to standardize the cell
dataset = spglib.get_symmetry_dataset((cell, positions_frac, atomic_numbers))
print(f"Space group number: {dataset['number']}")
print(f"Space group symbol: {dataset['international']}")
print(f"Point group: {dataset['pointgroup']}")

# Get standardized cell
std_cell = np.array(dataset['std_lattice'])
std_positions = np.array(dataset['std_positions'])
std_atomic_numbers = np.array(dataset['std_types'])

print(f"\nStandardized cell:\n{std_cell}")
print(f"Standardized positions:\n{std_positions}")
print(f"Standardized atomic numbers: {std_atomic_numbers}")

# ============================================================================
print("\n" + "=" * 70)
print("ASE.bandpath on STANDARDIZED cell:")
print("=" * 70)

try:
    result_ase_std = ase_bandpath(None, std_cell, npoints=50)
    xs_ase_std, Xs_ase_std, knames_ase_std = result_ase_std.get_linear_kpoint_axis()
    print(f"Special points: {list(result_ase_std.special_points.keys())}")
    print(f"Path: {result_ase_std.path}")
    print(f"knames: {knames_ase_std}")
    print(f"kpts shape: {result_ase_std.kpts.shape}")
except Exception as e:
    print(f"ERROR: {e}")

# ============================================================================
print("\n" + "=" * 70)
print("SeekPath on STANDARDIZED structure:")
print("=" * 70)

try:
    result_seekpath_std = seekpath.get_path(
        (std_cell.tolist(), std_positions.tolist(), std_atomic_numbers.tolist())
    )
    print(f"Special points: {list(result_seekpath_std['point_coords'].keys())}")
    print(f"Path segments: {result_seekpath_std['path']}")
    print(f"Bravais lattice: {result_seekpath_std['bravais_lattice']}")
except Exception as e:
    print(f"ERROR: {e}")
