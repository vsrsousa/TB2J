#!/usr/bin/env python
"""
Generate magnon band with SeekPath patch using real TB2J data
"""
import sys
import pickle
import warnings
import numpy as np

warnings.filterwarnings("ignore")

work_dir = "/Users/vinicius/Documents/gdni4al/TB2J_results_Gd_rcut10_kpt9"

print("=" * 90)
print("GENERATING MAGNON BAND WITH SEEKPATH PATCH")
print("=" * 90)

# Load pickle
print(f"\n1. Loading TB2J data...")
with open(f"{work_dir}/TB2J.pickle", "rb") as f:
    data = pickle.load(f)

atoms = data['atoms']
cell = atoms.get_cell()
pos = atoms.get_scaled_positions()
atomic_numbers = atoms.get_atomic_numbers()

print(f"   ✓ Loaded structure:")
print(f"     Cell: {cell.shape}")
print(f"     Atoms: {len(atoms)} ({atoms.get_chemical_symbols()})")

# Import patch BEFORE testing
from TB2J.seekpath_patch import set_structure_context
from ase.cell import Cell

# Test 1: WITHOUT patch (ASE baseline)
print(f"\n2. Generating k-path WITHOUT context (ASE)...")
bp_ase = Cell(cell).bandpath(npoints=100)
print(f"   Result: {len(bp_ase.special_points)} special points, {bp_ase.kpts.shape[0]} k-points")
print(f"   Special points: {list(bp_ase.special_points.keys())}")

# Test 2: WITH patch (SeekPath)
print(f"\n3. Generating k-path WITH context (SeekPath)...")
with set_structure_context(cell, pos, atomic_numbers):
    bp_seekpath = Cell(cell).bandpath(npoints=100)

print(f"   Result: {len(bp_seekpath.special_points)} special points, {bp_seekpath.kpts.shape[0]} k-points")
print(f"   Special points: {sorted(bp_seekpath.special_points.keys())}")

# Compare
print(f"\n" + "=" * 90)
print("COMPARISON")
print("=" * 90)
diff = len(bp_seekpath.special_points) - len(bp_ase.special_points)
print(f"\nASE:      {len(bp_ase.special_points):2d} special points")
print(f"SeekPath: {len(bp_seekpath.special_points):2d} special points")
print(f"Gain:     {diff:+2d} special points ({100*diff/max(len(bp_ase.special_points), 1):+.0f}%)")

if diff > 0:
    print(f"\n✓ SEEKPATH IS WORKING!")
    print(f"\n  Extra points discovered by SeekPath:")
    extra = set(bp_seekpath.special_points.keys()) - set(bp_ase.special_points.keys())
    for name in sorted(extra):
        coord = bp_seekpath.special_points[name]
        print(f"    {name:10s}: {coord}")
else:
    print(f"\n⚠ WARNING: SeekPath did not find more points than ASE")
    print(f"  Possible reasons:")
    print(f"    - Context not being used")
    print(f"    - SeekPath fallback to dummy atom")
    print(f"    - Issue with patch application")

print(f"\n" + "=" * 90)
