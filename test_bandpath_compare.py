#!/usr/bin/env python
"""Compare seekpath vs ASE bandpath output for debugging."""

import numpy as np
from ase.build import bulk
from ase.cell import Cell
from ase.dft.kpoints import bandpath as ase_bandpath

# Disable seekpath patch temporarily to get ASE original
import sys
import importlib

# Test with simple structure (1 atom)
print("=" * 70)
print("Comparing SeekPath vs ASE BandPath")
print("=" * 70)
print()

# Single atom in cubic cell
atoms = bulk('Fe', 'bcc', a=2.87)  # Body-centered cubic
print(f"Test structure: {atoms.get_chemical_formula()}")
print(f"Cell: {atoms.get_cell()}")
print()

# ===== ASE Original =====
print("ASE Original (ase.dft.kpoints.bandpath):")
try:
    ase_result = ase_bandpath([], atoms.get_cell(), npoints=50)
    print(f"  Type: {type(ase_result).__name__}")
    print(f"  kpts shape: {np.array(ase_result[0]).shape if isinstance(ase_result, tuple) else np.array(ase_result.kpts).shape}")
    if hasattr(ase_result, 'special_points'):
        print(f"  special_points: {ase_result.special_points}")
    if hasattr(ase_result, 'get_linear_kpoint_axis'):
        x, Xs, knames = ase_result.get_linear_kpoint_axis()
        print(f"  knames: {knames}")
        print(f"  Xs positions: {Xs}")
except Exception as e:
    print(f"  Error: {e}")
    import traceback
    traceback.print_exc()

print()

# ===== SeekPath (via patch) =====
print("SeekPath (via monkey-patch):")
try:
    import TB2J  # This applies the patch
    import ase.dft.kpoints as kp_patched
    from ase.cell import Cell
    
    seekpath_result = kp_patched.bandpath([], atoms.get_cell(), npoints=50)
    print(f"  Type: {type(seekpath_result).__name__}")
    print(f"  kpts shape: {np.array(seekpath_result.kpts).shape if hasattr(seekpath_result, 'kpts') else 'N/A'}")
    if hasattr(seekpath_result, 'special_points'):
        print(f"  special_points: {seekpath_result.special_points}")
    if hasattr(seekpath_result, 'get_linear_kpoint_axis'):
        x, Xs, knames = seekpath_result.get_linear_kpoint_axis()
        print(f"  knames: {knames}")
        print(f"  Xs positions: {Xs}")
        print(f"  x (full): {x[:10]}... (first 10)")
except Exception as e:
    print(f"  Error: {e}")
    import traceback
    traceback.print_exc()

print()

# ===== Test group_band_path =====
print("Testing group_band_path with both:")
from TB2J.spinham.plot import group_band_path

try:
    # Reload to get original ASE
    from importlib import reload
    import ase.dft.kpoints as kp_orig
    reload(kp_orig)
    ase_bp = ase_bandpath([], atoms.get_cell(), npoints=50)
    
    print(f"\n  ASE original via group_band_path:")
    xlist_ase, kptlist_ase, Xs_ase, knames_ase = group_band_path(ase_bp)
    print(f"    xlist segments: {len(xlist_ase)}")
    print(f"    kptlist shapes: {[k.shape for k in kptlist_ase]}")
    print(f"    knames: {knames_ase}")
except Exception as e:
    print(f"    Error: {e}")
    import traceback
    traceback.print_exc()

try:
    # Now with seekpath patch
    import TB2J
    seekpath_bp = kp_patched.bandpath([], atoms.get_cell(), npoints=50)
    
    print(f"\n  SeekPath via group_band_path:")
    xlist_sk, kptlist_sk, Xs_sk, knames_sk = group_band_path(seekpath_bp)
    print(f"    xlist segments: {len(xlist_sk)}")
    print(f"    kptlist shapes: {[k.shape for k in kptlist_sk]}")
    print(f"    knames: {knames_sk}")
except Exception as e:
    print(f"    Error: {e}")
    import traceback
    traceback.print_exc()

print()
print("=" * 70)
