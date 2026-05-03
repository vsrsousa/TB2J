#!/usr/bin/env python
"""
Test ASE bandpath generation CORRECTLY for GdNi4Si
Try different approaches and parameters
"""
import numpy as np
from ase.io import read
from ase.dft.kpoints import BandPath, bandpath
from ase.cell import Cell
import seekpath

# Read structure
atoms = read('/Users/vinicius/Documents/cifs/GdNi4Si.cif')
cell = atoms.get_cell()
positions_frac = atoms.get_scaled_positions()
atomic_numbers = atoms.get_atomic_numbers()

print("=" * 80)
print("GdNi4Si: TESTING ASE BANDPATH CORRECTLY")
print("=" * 80)

print(f"\nStructure: orthorhombic (Cmmm, space group 65)")
print(f"Cell: {cell.cellpar()}")

# ========== TEST 1: BandPath with path=None ==========
print("\n### TEST 1: BandPath(path=None, cell=cell) ###")
try:
    bp1 = BandPath(path=None, cell=cell)
    print(f"Path: '{bp1.path}'")
    print(f"Special points: {bp1.special_points}")
    print(f"K-points shape: {bp1.kpts.shape}")
except Exception as e:
    print(f"Error: {e}")

# ========== TEST 2: BandPath with density parameter ==========
print("\n### TEST 2: BandPath(path=None, cell=cell, density=10) ###")
try:
    bp2 = BandPath(path=None, cell=cell, density=10)
    print(f"Path: '{bp2.path}'")
    print(f"Special points: {bp2.special_points}")
    print(f"K-points shape: {bp2.kpts.shape}")
except Exception as e:
    print(f"Error: {e}")

# ========== TEST 3: bandpath function (module level) ==========
print("\n### TEST 3: bandpath(path=None, cell=cell) ###")
try:
    bp3 = bandpath(None, cell)
    print(f"Path: '{bp3.path}'")
    print(f"Special points: {bp3.special_points}")
    print(f"K-points shape: {bp3.kpts.shape}")
except Exception as e:
    print(f"Error: {e}")

# ========== TEST 4: Cell.bandpath() method ==========
print("\n### TEST 4: Cell(cell).bandpath() ###")
try:
    bp4 = Cell(cell).bandpath()
    print(f"Path: '{bp4.path}'")
    print(f"Special points: {bp4.special_points}")
    print(f"K-points shape: {bp4.kpts.shape}")
except Exception as e:
    print(f"Error: {e}")

# ========== TEST 5: Check if ASE can find path with spglib ==========
print("\n### TEST 5: ASE with spglib standardization ###")
try:
    from ase.build import bulk
    from ase.spacegroup import get_spacegroup
    
    sg = get_spacegroup(atoms)
    print(f"Space group detected: {sg}")
    
    # Try bandpath again
    bp5 = BandPath(path=None, cell=cell)
    print(f"Path after spacegroup check: '{bp5.path}'")
    print(f"K-points shape: {bp5.kpts.shape}")
except Exception as e:
    print(f"Error: {e}")

# ========== TEST 6: Try with manual hexagonal path ==========
print("\n### TEST 6: BandPath with explicit path string ###")
try:
    # Try common paths for orthorhombic
    for path_str in ['GMK', 'GMKX', 'GMRX', 'GAMMA', 'G', '']:
        try:
            bp6 = BandPath(path=path_str, cell=cell)
            if bp6.kpts.shape[0] > 0:
                print(f"Path '{path_str}': Found {bp6.kpts.shape[0]} k-points")
                print(f"  Special points: {bp6.special_points}")
                break
        except:
            pass
except Exception as e:
    print(f"Error: {e}")

# ========== SeekPath for comparison ==========
print("\n" + "=" * 80)
print("### COMPARISON: SeekPath (ground truth) ###")
print("=" * 80)

structure = (cell.tolist(), positions_frac.tolist(), atomic_numbers.tolist())
result = seekpath.get_path(structure)

print(f"\nSpace group: {result['spacegroup_international']} ({result['spacegroup_number']})")
print(f"Bravais lattice: {result['bravais_lattice']}")
print(f"Special points: {len(result['point_coords'])}")
print(f"Path segments: {len(result['path'])}")

print("\nSeekPath special points:")
for name in sorted(result['point_coords'].keys()):
    print(f"  {name}")

print("\nSeekPath path:")
for i, (start, end) in enumerate(result['path'], 1):
    print(f"  {i:2d}. {start} → {end}")
