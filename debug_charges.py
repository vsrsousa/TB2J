#!/usr/bin/env python
"""Debug script to trace charge/magmom calculation step-by-step."""

import sys
import numpy as np
sys.path.insert(0, '/Users/vinicius/opt/projects/TB2J')

# Set up paths
import os
os.chdir('/Users/vinicius/opt/data/TB2J_examples/Siesta/bccFe/DFT')

# Import
from HamiltonIO.siesta import SislParser
from TB2J.exchangeCL2 import ExchangeCL2

print("=" * 70)
print("DEBUG: Charge/Magmom Calculation Trace (seekpath-clean branch)")
print("=" * 70)

try:
    # Parse DFT output
    parser = SislParser(fdf_fname='siesta.fdf', read_H_soc=False)
    is_colinear = parser.read_spin() == 2
    H_up = parser.get_model(ispin=0)
    H_dn = parser.get_model(ispin=1)
    print(f"\n✓ Parsed Hamiltonian: {H_up.nbasis} basis functions, colinear={is_colinear}")
    
    # Create exchange calculator
    kmesh = [9, 9, 9]
    emin = -10
    
    exchange = ExchangeCL2(
        efermi=-2.5,  # rough estimate
        kmesh=kmesh,
        emin=emin,
        output_path='TB2J_results_debug',
    )
    
    # Set models for colinear case
    exchange.set_tbmodels((H_up, H_dn))
    
    print(f"\n1. Green's function objects created:")
    print(f"   - Gup.nbasis: {exchange.Gup.nbasis}")
    print(f"   - Gdn.nbasis: {exchange.Gdn.nbasis}")
    print(f"   - Gup.nkpts: {exchange.Gup.nkpts}")
    
    # Check kweights
    print(f"\n2. K-point information:")
    print(f"   - Number of k-points: {exchange.Gup.nkpts}")
    print(f"   - K-weights shape: {np.array(exchange.Gup.kweights).shape}")
    print(f"   - K-weights sum: {sum(exchange.Gup.kweights):.6f}")
    print(f"   - K-weight[0]: {exchange.Gup.kweights[0]:.6f}")
    
    # Check density matrices
    print(f"\n3. Density matrix traces:")
    print(f"   - rho_up shape: {exchange.rho_up.shape}")
    print(f"   - rho_dn shape: {exchange.rho_dn.shape}")
    rho_up_trace = np.trace(exchange.rho_up)
    rho_dn_trace = np.trace(exchange.rho_dn)
    print(f"   - trace(rho_up): {rho_up_trace:.6f}")
    print(f"   - trace(rho_dn): {rho_dn_trace:.6f}")
    print(f"   - Total electrons: {rho_up_trace + rho_dn_trace:.6f}")
    
    # Extract charges
    exchange.get_rho_atom()
    
    print(f"\n4. Orbital mapping:")
    print(f"   - Number of atoms: {len(exchange.atoms)}")
    print(f"   - Atom types: {exchange.atoms}")
    for iatom in exchange.orb_dict:
        iorb = exchange.iorb(iatom)
        print(f"   - Atom {iatom}: orbitals {iorb} (count: {len(iorb)})")
    
    print(f"\n5. Per-atom charges and moments:")
    for iatom in exchange.orb_dict:
        iorb = exchange.iorb(iatom)
        tup = np.real(np.trace(exchange.rho_up[np.ix_(iorb, iorb)]))
        tdn = np.real(np.trace(exchange.rho_dn[np.ix_(iorb, iorb)]))
        charge = tup + tdn
        magmom = tup - tdn
        
        print(f"   Atom {iatom}:")
        print(f"      - Orbitals: {list(iorb)}")
        print(f"      - trace(rho_up[orbs]): {tup:.6f}")
        print(f"      - trace(rho_dn[orbs]): {tdn:.6f}")
        print(f"      - Charge: {charge:.6f} (expect ~8)")
        print(f"      - Magmom(z): {magmom:.6f} (expect ~2.3)")
    
    print(f"\n6. Final charges and moments (from exchange object):")
    print(f"   - Charges: {exchange.charges}")
    print(f"   - Spinat[:,2]: {exchange.spinat[:, 2]}")
    
    print("\n✅ Completed successfully!")
    
except Exception as e:
    print(f"\n❌ Error: {e}")
    import traceback
    traceback.print_exc()
