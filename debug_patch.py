#!/usr/bin/env python
"""
Debug: Verify if patch is actually being applied
"""
import sys
import warnings

# Capture all warnings
warnings.simplefilter("always")

# Check what's happening at import time
print("=" * 80)
print("DEBUG: Checking seekpath_patch activation")
print("=" * 80)

print(f"\n1. Checking if seekpath is available...")
try:
    import seekpath
    print(f"   ✓ seekpath found: {seekpath.__version__ if hasattr(seekpath, '__version__') else 'version unknown'}")
except Exception as e:
    print(f"   ✗ seekpath NOT available: {e}")

print(f"\n2. Checking if ASE is available...")
try:
    import ase
    from ase.dft import kpoints as ase_kpoints
    print(f"   ✓ ASE found: {ase.__version__}")
except Exception as e:
    print(f"   ✗ ASE NOT available: {e}")

print(f"\n3. Checking original ASE bandpath function...")
print(f"   Original bandpath: {ase_kpoints.bandpath}")

print(f"\n4. Importing TB2J (should apply patch)...")
import TB2J
print(f"   ✓ TB2J imported")

print(f"\n5. Checking if bandpath was patched...")
print(f"   Patched bandpath: {ase_kpoints.bandpath}")

# Check if it's been replaced
if "bandpath_wrapper" in str(ase_kpoints.bandpath):
    print(f"   ✓ PATCH IS ACTIVE (wrapper function detected)")
else:
    print(f"   ✗ PATCH NOT ACTIVE (still original function)")

print(f"\n6. Checking Cell.bandpath...")
from ase.cell import Cell
print(f"   Cell.bandpath: {Cell.bandpath}")

if "cell_bandpath_wrapper" in str(Cell.bandpath):
    print(f"   ✓ PATCH IS ACTIVE (wrapper function detected)")
else:
    print(f"   ✗ PATCH NOT ACTIVE (still original function)")

# Now test actual behavior
print("\n" + "=" * 80)
print("DEBUG: Testing actual patch behavior")
print("=" * 80)

import numpy as np
from ase import Atoms

# Simple hexagonal cell
cell_hex = np.array([
    [2.48, -4.295486, 0.0],
    [2.48, 4.295486, 0.0],
    [0.0, 0.0, 4.037],
])

atoms_data = [
    ('Gd', np.array([-0.0, 0.0, 0.0])),
    ('Ni', np.array([0.666667, 0.333333, 0.0])),
    ('Ni', np.array([0.333333, 0.666667, 0.0])),
    ('Ni', np.array([-0.0, 0.5, 0.5])),
    ('Ni', np.array([0.5, 0.0, 0.5])),
    ('Al', np.array([0.5, 0.5, 0.5])),
]

symbols = [sym for sym, _ in atoms_data]
positions = np.array([pos for _, pos in atoms_data])

atoms = Atoms(symbols=symbols, scaled_positions=positions, cell=cell_hex, pbc=True)
cell = atoms.get_cell()
pos = atoms.get_scaled_positions()
nums = atoms.get_atomic_numbers()

from TB2J.seekpath_patch import set_structure_context

print(f"\nTest 1: WITHOUT structure context (should fallback to ASE)...")
bp_no_context = Cell(cell).bandpath(npoints=50)
print(f"  Result: {len(bp_no_context.special_points)} special points")

print(f"\nTest 2: WITH structure context (should use SeekPath)...")
with set_structure_context(cell, pos, nums):
    bp_with_context = Cell(cell).bandpath(npoints=50)
print(f"  Result: {len(bp_with_context.special_points)} special points")

print(f"\nTest 3: Direct check - is _structure_context being used?...")
from TB2J import seekpath_patch
print(f"  _structure_context = {seekpath_patch._structure_context}")

# Test setting context
print(f"\nTest 4: Setting context manually and checking...")
with set_structure_context(cell, pos, nums):
    print(f"  Inside context: _structure_context = {seekpath_patch._structure_context}")

print(f"  Outside context: _structure_context = {seekpath_patch._structure_context}")

# Try direct function call
print(f"\n" + "=" * 80)
print("DEBUG: Testing _seekpath_bandpath_from_cell directly")
print("=" * 80)

seekpath_patch._structure_context['cell'] = cell
seekpath_patch._structure_context['positions'] = pos
seekpath_patch._structure_context['atomic_numbers'] = nums

print(f"\nCalling _seekpath_bandpath_from_cell with context set...")
try:
    bp_direct = seekpath_patch._seekpath_bandpath_from_cell(cell, npoints=50)
    print(f"  ✓ Result: {len(bp_direct.special_points)} special points")
    for name in sorted(bp_direct.special_points.keys()):
        print(f"    {name}")
except Exception as e:
    print(f"  ✗ Error: {e}")
    import traceback
    traceback.print_exc()
