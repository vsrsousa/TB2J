#!/usr/bin/env python
import pickle
import numpy as np

work_dir = "/Users/vinicius/Documents/gdni4al/TB2J_results_Gd_rcut10_kpt9"

with open(f"{work_dir}/TB2J.pickle", "rb") as f:
    data = pickle.load(f)

print("Pickle content keys:", list(data.keys()))
print("\nAtoms info:")
atoms = data['atoms']
print(f"  Type: {type(atoms)}")
print(f"  Number of atoms: {len(atoms)}")
print(f"  Atomic numbers: {atoms.get_atomic_numbers()}")
print(f"  Symbols: {atoms.get_chemical_symbols()}")

print("\nSpinHamiltonian info (if present):")
if 'SpinHamiltonian' in data:
    sh = data['SpinHamiltonian']
    print(f"  natom: {sh.natom if hasattr(sh, 'natom') else 'N/A'}")
    print(f"  zion: {sh.zion if hasattr(sh, 'zion') else 'N/A'}")
    print(f"  pos shape: {sh.pos.shape if hasattr(sh, 'pos') else 'N/A'}")
else:
    print("  SpinHamiltonian not in pickle directly")
