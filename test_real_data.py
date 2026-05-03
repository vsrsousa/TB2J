#!/usr/bin/env python
"""
Test seekpath_patch with real TB2J data from Multibinit folder
Load TB2J.pickle and generate magnon band structure
"""
import sys
import pickle
import warnings

warnings.filterwarnings("ignore")

# Set up path
work_dir = "/Users/vinicius/Documents/gdni4al/TB2J_results_Gd_rcut10_kpt9"
sys.path.insert(0, work_dir)

print("=" * 90)
print("TEST: SeekPath Patch with Real TB2J Data")
print("=" * 90)

# Load TB2J pickle
print(f"\n1. Loading TB2J data from {work_dir}/TB2J.pickle...")
try:
    with open(f"{work_dir}/TB2J.pickle", "rb") as f:
        tb2j = pickle.load(f)
    print(f"   ✓ Loaded successfully")
except Exception as e:
    print(f"   ✗ Error: {e}")
    sys.exit(1)

# Check structure
print(f"\n2. Structure information:")
print(f"   Cell shape: {tb2j.cell.shape}")
print(f"   Positions shape: {tb2j.pos.shape}")
print(f"   Atomic numbers: {tb2j.zion if hasattr(tb2j, 'zion') else 'N/A'}")

# Get spin hamiltonian
print(f"\n3. Getting SpinHamiltonian...")
try:
    # The TB2J object might have different structure
    # Let's try to get the magnetic structure
    if hasattr(tb2j, 'cell'):
        print(f"   TB2J object has cell attribute")
    if hasattr(tb2j, 'spinham'):
        spinham = tb2j.spinham
        print(f"   ✓ Found spinham")
    else:
        print(f"   ✗ No spinham attribute")
        # Try to list available attributes
        attrs = [a for a in dir(tb2j) if not a.startswith('_')]
        print(f"   Available attributes: {attrs[:10]}")
except Exception as e:
    print(f"   ✗ Error: {e}")

# Try to generate magnon band directly
print(f"\n4. Testing magnon band generation...")
try:
    from TB2J.seekpath_patch import set_structure_context
    from ase.cell import Cell
    import numpy as np
    
    # Get structure data
    cell = tb2j.cell
    pos = tb2j.pos
    zion = [int(z) for z in tb2j.zion] if hasattr(tb2j, 'zion') else None
    
    print(f"   Structure data:")
    print(f"     cell: {cell.shape}")
    print(f"     pos: {pos.shape}")
    print(f"     zion: {zion}")
    
    # Test WITHOUT context (ASE)
    print(f"\n   a) WITHOUT context (ASE behavior):")
    bp_ase = Cell(cell).bandpath(npoints=50)
    print(f"      ASE: {len(bp_ase.special_points)} special points, {bp_ase.kpts.shape[0]} k-points")
    
    # Test WITH context (SeekPath)
    print(f"\n   b) WITH context (SeekPath):")
    with set_structure_context(cell, pos, zion):
        bp_seekpath = Cell(cell).bandpath(npoints=50)
    print(f"      SeekPath: {len(bp_seekpath.special_points)} special points, {bp_seekpath.kpts.shape[0]} k-points")
    
    # Show improvement
    diff = len(bp_seekpath.special_points) - len(bp_ase.special_points)
    print(f"\n   Improvement: {diff:+d} special points ({100*diff/len(bp_ase.special_points):+.0f}%)")
    
    if len(bp_seekpath.special_points) > len(bp_ase.special_points):
        print(f"\n   ✓ SEEKPATH IS WORKING - Found more special points!")
        print(f"\n   SeekPath special points:")
        for name in sorted(bp_seekpath.special_points.keys()):
            print(f"     {name}")
    else:
        print(f"\n   ⚠ Same number of points - check if context is being used")
        
except Exception as e:
    print(f"   ✗ Error: {e}")
    import traceback
    traceback.print_exc()

print("\n" + "=" * 90)
