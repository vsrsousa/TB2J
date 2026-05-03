#!/usr/bin/env python
"""Debug compute_charge_and_magnetic_moments"""

import sys
sys.path.insert(0, '/Users/vinicius/opt/projects/TB2J')

import os
os.chdir('/Users/vinicius/opt/data/TB2J_examples/Siesta/bccFe/DFT')

import numpy as np
from TB2J.interfaces.siesta_interface import gen_exchange_siesta
import TB2J.exchangeCL2 as excl2_mod

print("=" * 70)
print("DEBUG: compute_charge_and_magnetic_moments details")
print("=" * 70)

original_compute = excl2_mod.ExchangeCL2.compute_charge_and_magnetic_moments

def debug_compute_charges(self):
    print("\n[compute_charge_and_magnetic_moments]")
    
    if not hasattr(self, "G_diagonal_up") or not self.G_diagonal_up:
        print("  No Green's function diagonals!")
        return

    self.charges = np.zeros(len(self.atoms))
    self.spinat = np.zeros((len(self.atoms), 3))

    for iatom in range(len(self.atoms)):
        if not self.G_diagonal_up[iatom] or not self.G_diagonal_dn[iatom]:
            continue
        
        print(f"\n  Atom {iatom}:")
        
        G_up_diags = np.array(self.G_diagonal_up[iatom])
        G_dn_diags = np.array(self.G_diagonal_dn[iatom])
        
        print(f"    Shape: {G_up_diags.shape} (n_energies, n_orbitals)")
        print(f"    G_up_diags[0,0]: {G_up_diags[0,0]:.6f}")
        print(f"    G_up_diags dtype: {G_up_diags.dtype}")
        
        # Contour integration
        result_up = self.contour.integrate_values(G_up_diags)
        result_dn = self.contour.integrate_values(G_dn_diags)
        
        print(f"    contour.integrate_values(G_up): {result_up:.6f}")
        print(f"    imag part: {np.imag(result_up):.6f}")
        print(f"    -imag/π: {-np.imag(result_up)/np.pi:.6f}")
        
        integrated_up = -np.imag(result_up) / np.pi
        integrated_dn = -np.imag(result_dn) / np.pi
        
        print(f"    sum: {np.sum(integrated_up) + np.sum(integrated_dn):.6f}")
        
        self.charges[iatom] = np.sum(integrated_up) + np.sum(integrated_dn)
        self.spinat[iatom, 2] = np.sum(integrated_up) - np.sum(integrated_dn)

excl2_mod.ExchangeCL2.compute_charge_and_magnetic_moments = debug_compute_charges

gen_exchange_siesta(
    fdf_fname='siesta.fdf',
    kmesh=[9, 9, 9],
    emin=-10,
    magnetic_elements=['Fe'],
)

print("\n✅ Done")
