#!/usr/bin/env python
"""
Debug: Check if set_structure_context is working in plot_magnon_band
"""
import sys
import pickle
from TB2J.spinham.hamiltonian import SpinHamiltonian
from ase.cell import Cell

work_dir = "/Users/vinicius/Documents/gdni4al/TB2J_results_Gd_rcut10_kpt9"

# Load pickle
with open(f"{work_dir}/TB2J.pickle", "rb") as f:
    data = pickle.load(f)

atoms = data['atoms']
cell = atoms.get_cell()
pos = atoms.get_scaled_positions()
atomic_numbers = atoms.get_atomic_numbers()

print("=" * 80)
print("DEBUG: set_structure_context in plot_magnon_band")
print("=" * 80)

# Test 1: Without context (simulating old code)
print("\n1. WITHOUT set_structure_context:")
bp1 = Cell(cell).bandpath(npoints=50)
print(f"   Result: {len(bp1.special_points)} special points")
print(f"   Points: {sorted(bp1.special_points.keys())}")

# Test 2: With context (what plot_magnon_band does)
print("\n2. WITH set_structure_context:")
from TB2J.seekpath_patch import set_structure_context

context_manager = set_structure_context(
    cell, pos, 
    [int(z) for z in atomic_numbers]
)

if context_manager:
    context_manager.__enter__()

try:
    bp2 = Cell(cell).bandpath(npoints=50)
    print(f"   Result: {len(bp2.special_points)} special points")
    print(f"   Points: {sorted(bp2.special_points.keys())}")
finally:
    if context_manager:
        context_manager.__exit__(None, None, None)

print("\n" + "=" * 80)
if len(bp2.special_points) > len(bp1.special_points):
    print(f"✓ SEEKPATH WORKING: +{len(bp2.special_points) - len(bp1.special_points)} points")
else:
    print(f"✗ SEEKPATH NOT WORKING: Same number of points")
