#!/usr/bin/env python
"""
Test new structure from user input
"""
import numpy as np
from ase import Atoms
from ase.cell import Cell
import warnings

warnings.filterwarnings("ignore")

# Parse structure from user input
atoms_data = [
    ('Gd', np.array([-0.0000000000, 0.0000000000, 0.0000000000])),
    ('Ni', np.array([0.6666670000, 0.3333330000, 0.0000000000])),
    ('Ni', np.array([0.3333330000, 0.6666670000, 0.0000000000])),
    ('Ni', np.array([-0.0000000000, 0.5000000000, 0.5000000000])),
    ('Ni', np.array([0.5000000000, 0.0000000000, 0.5000000000])),
    ('Al', np.array([0.5000000000, 0.5000000000, 0.5000000000])),
]

cell_vectors = np.array([
    [2.48000000000000, -4.29548600000000, 0.00000000000000],
    [2.48000000000000, 4.29548600000000, 0.00000000000000],
    [0.00000000000000, 0.00000000000000, 4.03700000000000],
])

# Build ASE Atoms object
symbols = [sym for sym, _ in atoms_data]
positions = np.array([pos for _, pos in atoms_data])

atoms = Atoms(
    symbols=symbols,
    scaled_positions=positions,
    cell=cell_vectors,
    pbc=True
)

cell = atoms.get_cell()
positions_frac = atoms.get_scaled_positions()
atomic_numbers = atoms.get_atomic_numbers()

print("=" * 80)
print("TESTING NEW STRUCTURE: GdNi2AlNi2 (Hexagonal)")
print("=" * 80)

print(f"\nStructure composition:")
for i, (sym, pos) in enumerate(atoms_data, 1):
    print(f"  {i}. {sym:3s}: {pos}")

print(f"\nCell parameters:")
cellpars = atoms.get_cell().cellpar()
print(f"  a = {cellpars[0]:.4f} Å")
print(f"  b = {cellpars[1]:.4f} Å")
print(f"  c = {cellpars[2]:.4f} Å")
print(f"  α = {cellpars[3]:.1f}°, β = {cellpars[4]:.1f}°, γ = {cellpars[5]:.1f}°")

print(f"\nCell vectors:")
print(cell)

# Test ASE
print("\n" + "=" * 80)
print("### ASE bandpath ###")
bp_ase = Cell(cell).bandpath(npoints=50)
print(f"Special points: {len(bp_ase.special_points)}")
if bp_ase.special_points:
    for name in sorted(bp_ase.special_points.keys()):
        print(f"  {name}")
print(f"K-points: {bp_ase.kpts.shape}")

# Test SeekPath
print("\n" + "=" * 80)
print("### SeekPath ###")
from TB2J.seekpath_patch import set_structure_context

with set_structure_context(cell, positions_frac, atomic_numbers):
    bp_seekpath = Cell(cell).bandpath(npoints=50)

print(f"Special points: {len(bp_seekpath.special_points)}")
for name in sorted(bp_seekpath.special_points.keys()):
    coord = bp_seekpath.special_points[name]
    print(f"  {name:10s}: {coord}")
print(f"K-points: {bp_seekpath.kpts.shape}")

knames = getattr(bp_seekpath, '_knames', [])
if knames:
    # Show unique segments
    segments = []
    for i in range(len(knames) - 1):
        seg = f"{knames[i]} → {knames[i + 1]}"
        if seg not in segments:
            segments.append(seg)
    print(f"\nPath segments ({len(segments)}):")
    for i, seg in enumerate(segments, 1):
        print(f"  {i:2d}. {seg}")

# Comparison
print("\n" + "=" * 80)
print("COMPARISON")
print("=" * 80)
print(f"""
ASE:      {len(bp_ase.special_points)} special points, {bp_ase.kpts.shape[0]} k-points
SeekPath: {len(bp_seekpath.special_points)} special points, {bp_seekpath.kpts.shape[0]} k-points

Improvement: {len(bp_seekpath.special_points) - len(bp_ase.special_points):+d} special points
""")
