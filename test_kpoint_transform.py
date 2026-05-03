#!/usr/bin/env python
"""
Test: Convert k-points from primitive (seekpath) back to original hexagonal cell
"""
import numpy as np
from ase.io import read
import seekpath

# Read structure
atoms = read('/Users/vinicius/Documents/cifs/Gd.cif')
cell = atoms.get_cell()
positions_frac = atoms.get_scaled_positions()
atomic_numbers = atoms.get_atomic_numbers()

print("=" * 80)
print("ORIGINAL HEXAGONAL CELL")
print("=" * 80)
print(f"Cell:\n{cell}")

# Get seekpath result
structure = (cell.tolist(), positions_frac.tolist(), atomic_numbers.tolist())
result = seekpath.get_path(structure)

print(f"\nPrimitive transformation matrix (original → primitive):")
P = np.array(result['primitive_transformation_matrix'])
print(P)

print(f"\nInverse matrix (primitive → original):")
P_inv = np.array(result['inverse_primitive_transformation_matrix'])
print(P_inv)

# ============================================================================
# KEY: k-point transformation in RECIPROCAL space
# If real-space: r_prim = P @ r_original
# Then reciprocal-space: k_original = (P^-T) @ k_prim = (P^T)^-1 @ k_prim
# Or equivalently: k_original = k_prim @ (P^-1)^T
# ============================================================================

print("\n" + "=" * 80)
print("SPECIAL POINTS TRANSFORMATION")
print("=" * 80)

print(f"\nSeekPath special points (in PRIMITIVE k-space):")
for name, coord in result['point_coords'].items():
    print(f"  {name}: {coord}")

print(f"\nConverted to ORIGINAL HEXAGONAL k-space:")
print(f"Using: k_original = k_prim @ (P_inv)^T")

P_inv_T = P_inv.T  # Transpose of inverse
print(f"\n(P_inv)^T =\n{P_inv_T}")

for name, coord in result['point_coords'].items():
    k_prim = np.array(coord)
    # Transform from primitive to original reciprocal space
    k_orig = k_prim @ P_inv_T
    print(f"  {name}: {coord} → {k_orig}")

# ============================================================================
print("\n" + "=" * 80)
print("PATH VERIFICATION")
print("=" * 80)

print(f"\nSeekPath path segments:")
for seg in result['path']:
    k1 = np.array(result['point_coords'][seg[0]])
    k2 = np.array(result['point_coords'][seg[1]])
    
    k1_orig = k1 @ P_inv_T
    k2_orig = k2 @ P_inv_T
    
    print(f"  {seg[0]} → {seg[1]}")
    print(f"    Primitive: {k1} → {k2}")
    print(f"    Hexagonal: {k1_orig} → {k2_orig}")
