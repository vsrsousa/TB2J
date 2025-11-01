"""
Example: Calculate exchange interactions from PAOFLOW Hamiltonian

This example demonstrates how to use TB2J with PAOFLOW-generated Hamiltonians
for a collinear magnetic system (e.g., bcc Fe).

Prerequisites:
1. PAOFLOW calculation completed with hamiltonian.dat_0 and hamiltonian.dat_1
2. POSCAR or structure file
3. TB2J installed

Usage:
    python run_paoflow2J.py
"""

import os
import sys

# Check if required files exist
if not os.path.exists('hamiltonian.dat_0'):
    print("Error: hamiltonian.dat_0 not found")
    print("Please run PAOFLOW first and use paoflow.write_Hamiltonian()")
    sys.exit(1)

if not os.path.exists('hamiltonian.dat_1'):
    print("Error: hamiltonian.dat_1 not found")
    print("For collinear calculations, both spin channels are required")
    sys.exit(1)

if not os.path.exists('POSCAR'):
    print("Error: POSCAR not found")
    print("Please provide a structure file (POSCAR, cif, xyz, etc.)")
    sys.exit(1)

# Import TB2J function
from TB2J.manager import gen_exchange_paoflow

# Configuration
config = {
    'hr_fname': 'hamiltonian.dat_0',        # Spin-up Hamiltonian from PAOFLOW
    'hr_dn_fname': 'hamiltonian.dat_1',     # Spin-down Hamiltonian from PAOFLOW
    'atoms_fname': 'POSCAR',                # Structure file
    'colinear': True,                        # Collinear calculation
    'magnetic_elements': ['Fe'],             # Magnetic elements (modify as needed)
    'kmesh': [7, 7, 7],                     # k-point mesh
    'efermi': 0.0,                          # Fermi energy (eV) - PAOFLOW shifts to 0
    'emin': -10.0,                          # Min energy below E_F (eV)
    'emax': 0.5,                            # Max energy above E_F (eV)
    'nz': 100,                              # Energy integration steps
    'output_path': 'TB2J_results',          # Output directory
    'description': 'Exchange calculation from PAOFLOW Hamiltonian',
}

# Optional: Uncomment to enable parallel processing
# config['np'] = 4  # Use 4 CPU cores

# Optional: Uncomment to limit interaction range
# config['Rcut'] = 10.0  # Cutoff distance in Angstroms

print("=" * 60)
print("TB2J Exchange Calculation from PAOFLOW Hamiltonian")
print("=" * 60)
print(f"\nConfiguration:")
for key, value in config.items():
    print(f"  {key:20s}: {value}")
print()

# Run the calculation
try:
    gen_exchange_paoflow(**config)
    print("\n" + "=" * 60)
    print("Calculation completed successfully!")
    print("=" * 60)
    print(f"\nResults are in: {config['output_path']}/")
    print("\nTo plot magnon band structure:")
    print(f"  cd {config['output_path']}/Multibinit")
    print("  TB2J_magnon.py --figfname magnon.png")
    
except Exception as e:
    print("\n" + "=" * 60)
    print("Error during calculation:")
    print("=" * 60)
    print(f"\n{e}")
    import traceback
    traceback.print_exc()
    sys.exit(1)
