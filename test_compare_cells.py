#!/usr/bin/env python
"""
Compare primitive cells and scaled positions for ASE vs SeekPath
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
print("ORIGINAL Gd structure (from CIF):")
print("=" * 70)
print(f"Cell:\n{cell}")
print(f"Scaled positions:\n{positions_frac}")
print(f"Atomic numbers: {atomic_numbers}")

# ============================================================================
print("\n" + "=" * 70)
print("ASE.bandpath output:")
print("=" * 70)

result_ase = ase_bandpath(None, cell, npoints=50)
print(f"Path: {result_ase.path}")
print(f"Special points keys: {list(result_ase.special_points.keys())}")

# ASE doesn't directly return primitive cell, but let's check if it stores it
if hasattr(result_ase, 'cell'):
    print(f"Cell in bandpath object:\n{result_ase.cell}")

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

print(f"\nPrimitive positions (fractional):")
prim_positions = np.array(result_seekpath['primitive_positions'])
print(prim_positions)

print(f"\nPrimitive atomic types:")
prim_types = np.array(result_seekpath['primitive_types'])
print(prim_types)

print(f"\nConventional cell (from SeekPath):")
conv_cell = np.array(result_seekpath['conv_lattice'])
print(conv_cell)

print(f"\nConventional positions (fractional):")
conv_positions = np.array(result_seekpath['conv_positions'])
print(conv_positions)

print(f"\nConventional atomic types:")
conv_types = np.array(result_seekpath['conv_types'])
print(conv_types)

print(f"\nBravais lattice: {result_seekpath['bravais_lattice']}")
print(f"Primitive transformation matrix P:")
P = np.array(result_seekpath['primitive_transformation_matrix'])
print(P)

print(f"\nInverse P matrix:")
invP = np.array(result_seekpath['inverse_primitive_transformation_matrix'])
print(invP)

# Verify P relationship: conv_cell = prim_cell @ P
print(f"\n\nVerification:")
print(f"prim_cell @ P =\n{prim_cell @ P}")
print(f"conv_cell =\n{conv_cell}")
print(f"Match: {np.allclose(prim_cell @ P, conv_cell)}")

# ============================================================================
print("\n" + "=" * 70)
print("Comparison:")
print("=" * 70)
print(f"\nOriginal cell shape: {cell.shape}")
print(f"Original positions shape: {positions_frac.shape}")
print(f"Original atoms: {atomic_numbers}")

print(f"\nPrimitive cell shape: {prim_cell.shape}")
print(f"Primitive positions shape: {prim_positions.shape}")
print(f"Primitive atoms: {len(prim_types)}")

print(f"\nConventional cell shape: {conv_cell.shape}")
print(f"Conventional positions shape: {conv_positions.shape}")
print(f"Conventional atoms: {len(conv_types)}")

print(f"\nVolume ratio (original / primitive): {result_seekpath['volume_original_wrt_prim']}")
print(f"Volume ratio (original / conventional): {result_seekpath['volume_original_wrt_conv']}")
