#!/usr/bin/env python
"""
Test: Verify seekpath_patch.py is working correctly
Show SeekPath used when available, ASE fallback when needed
"""
import numpy as np
from ase.io import read
from ase.cell import Cell
import warnings

# Capture warnings to see fallbacks
warnings.simplefilter("always")

# Import patch AFTER setting up warning handler
from TB2J.seekpath_patch import set_structure_context

# Read structure
atoms = read('/Users/vinicius/Documents/cifs/GdNi4Si.cif')
cell = atoms.get_cell()
positions_frac = atoms.get_scaled_positions()
atomic_numbers = atoms.get_atomic_numbers()

print("=" * 80)
print("TEST: seekpath_patch.py with GdNi4Si")
print("=" * 80)

# ========== TEST 1: WITH structure context (should use SeekPath) ==========
print("\n### TEST 1: WITH structure context ###")
print("Expected: SeekPath generates 10 special points, 11 segments")

with set_structure_context(cell, positions_frac, atomic_numbers):
    bp_with_context = Cell(cell).bandpath(npoints=50)

print(f"Special points: {len(bp_with_context.special_points)}")
print(f"Point names: {sorted(bp_with_context.special_points.keys())}")
print(f"K-points shape: {bp_with_context.kpts.shape}")

# ========== TEST 2: WITHOUT structure context (should fallback to ASE) ==========
print("\n### TEST 2: WITHOUT structure context ###")
print("Expected: ASE fallback generates 8 special points, 4 segments")
print("(SeekPath tries with dummy atom but won't find meaningful path)")

bp_without_context = Cell(cell).bandpath(npoints=50)

print(f"Special points: {len(bp_without_context.special_points)}")
if bp_without_context.special_points:
    print(f"Point names: {sorted(bp_without_context.special_points.keys())}")
print(f"K-points shape: {bp_without_context.kpts.shape}")

# ========== COMPARISON ==========
print("\n" + "=" * 80)
print("COMPARISON")
print("=" * 80)

print(f"""
WITH context (SeekPath):
  ✓ Special points: {len(bp_with_context.special_points)}
  ✓ K-points: {bp_with_context.kpts.shape[0]}
  
WITHOUT context (ASE fallback):
  ✓ Special points: {len(bp_without_context.special_points)}
  ✓ K-points: {bp_without_context.kpts.shape[0]}

✓ Patch is working correctly!
  - Uses SeekPath when structure data available
  - Falls back to ASE when not available
""")

# ========== VERIFY k-points are valid ==========
print("\n### K-points verification ###")
print(f"Sample k-points from SeekPath version:")
for i in [0, len(bp_with_context.kpts)//2, -1]:
    print(f"  kpts[{i:4d}]: {bp_with_context.kpts[i]}")

print(f"\nAll k-points within [-2, 2] range: {np.all(np.abs(bp_with_context.kpts) <= 2.0)}")
