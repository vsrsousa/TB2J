#!/usr/bin/env python
"""
Test GdNi4Si structure with SeekPath k-path generation
"""
import numpy as np
from ase.io import read
from ase.cell import Cell

# Import patch to activate it
from TB2J.seekpath_patch import set_structure_context

# Read structure
atoms = read('/Users/vinicius/Documents/cifs/GdNi4Si.cif')
cell = atoms.get_cell()
positions_frac = atoms.get_scaled_positions()
atomic_numbers = atoms.get_atomic_numbers()

print("=" * 80)
print("GdNi4Si STRUCTURE ANALYSIS")
print("=" * 80)

print(f"\nCrystal structure: Cmmm (orthorhombic)")
print(f"Cell parameters:")
print(f"  a = {cell[0, 0]:.4f} Å")
print(f"  b = {cell[1, 1]:.4f} Å")
print(f"  c = {cell[2, 2]:.4f} Å")

print(f"\nCell vectors:")
print(cell)

print(f"\nAtomic positions (fractional):")
for i, (symbol, pos) in enumerate(zip(atoms.get_chemical_symbols(), positions_frac)):
    print(f"  {i+1}. {symbol}: {pos}")

print(f"\nAtomic numbers: {atomic_numbers}")

# Generate bandpath with context
print("\n" + "=" * 80)
print("SeekPath k-PATH GENERATION")
print("=" * 80)

with set_structure_context(cell, positions_frac, atomic_numbers):
    bp = Cell(cell).bandpath(npoints=50)

print(f"\nSpace group: Cmmm (65)")
print(f"Bravais lattice: oC (orthorhombic centered)")

print(f"\nSpecial points in reciprocal space:")
for name, coord in sorted(bp.special_points.items()):
    print(f"  {name:15s}: {coord}")

print(f"\nPath segments:")
knames = getattr(bp, '_knames', [])
if knames:
    segments = []
    for i in range(len(knames) - 1):
        if knames[i] != knames[i + 1]:
            seg = f"{knames[i]} → {knames[i + 1]}"
            if seg not in segments:
                segments.append(seg)
    
    for i, seg in enumerate(segments, 1):
        print(f"  {i:2d}. {seg}")

print(f"\nTotal path segments: {len(segments)}")
print(f"K-points array shape: {bp.kpts.shape}")

# Plot coordinates
x, Xs, knames_plot = bp.get_linear_kpoint_axis()
print(f"\nBand plot axis:")
print(f"  x range: {x[0]:.4f} to {x[-1]:.4f}")
print(f"  Special point positions: {[f'{x:.3f}' for x in Xs]}")
print(f"  Labels: {knames_plot}")

# Show some k-points
print(f"\nSample k-points:")
print(f"  First: {bp.kpts[0]}")
print(f"  Middle: {bp.kpts[len(bp.kpts)//2]}")
print(f"  Last: {bp.kpts[-1]}")
