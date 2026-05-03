#!/usr/bin/env python
"""
Show EXACTLY what SeekPath receives and returns
"""
import numpy as np
from ase.io import read
import seekpath

# Read Gd structure
atoms = read('/Users/vinicius/Documents/cifs/Gd.cif')
cell = atoms.get_cell()
positions_frac = atoms.get_scaled_positions()
atomic_numbers = atoms.get_atomic_numbers()

print("=" * 70)
print("INPUT to seekpath.get_path():")
print("=" * 70)
print(f"\nCell (input):")
print(cell)
print(f"type: {type(cell)}")

print(f"\nScaled positions (input):")
print(positions_frac)
print(f"type: {type(positions_frac)}")

print(f"\nAtomic numbers (input):")
print(atomic_numbers)
print(f"type: {type(atomic_numbers)}")

# Call seekpath with EXACT input
structure_input = (cell.tolist(), positions_frac.tolist(), atomic_numbers.tolist())
print(f"\nStructure tuple passed:")
print(f"cell: {structure_input[0]}")
print(f"positions: {structure_input[1]}")
print(f"numbers: {structure_input[2]}")

# ============================================================================
result = seekpath.get_path(structure_input)

print("\n" + "=" * 70)
print("OUTPUT from seekpath.get_path():")
print("=" * 70)

print(f"\nAll keys in result:")
print(result.keys())

print(f"\nBravais lattice: {result['bravais_lattice']}")
print(f"Bravais lattice extended: {result['bravais_lattice_extended']}")

print(f"\nPrimitive lattice:")
print(result['primitive_lattice'])

print(f"\nPrimitive positions:")
for i, pos in enumerate(result['primitive_positions']):
    print(f"  {i}: {pos}")

print(f"\nPrimitive types: {result['primitive_types']}")

print(f"\n### PATH ###")
print(f"Path segments: {result['path']}")
print(f"\nSpecial points:")
for name, coord in result['point_coords'].items():
    print(f"  {name}: {coord}")

print(f"\n### TRANSFORMATIONS ###")
print(f"Primitive transformation matrix:")
print(result['primitive_transformation_matrix'])
