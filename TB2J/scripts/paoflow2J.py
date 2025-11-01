#!/usr/bin/env python3
"""
Command-line interface for calculating exchange parameters from PAOflow output.

PAOflow is an alternative to Wannier90 for generating tight-binding Hamiltonians
from DFT calculations using projections of atomic orbitals.
"""

import argparse
import sys

from TB2J.exchange_params import add_exchange_args_to_parser, parser_argument_to_dict
from TB2J.versioninfo import print_license


def run_paoflow2J():
    """
    Main function for paoflow2J command-line tool.
    
    Calculates magnetic exchange parameters using magnetic force theorem
    from PAOflow tight-binding Hamiltonians.
    """
    print_license()
    parser = argparse.ArgumentParser(
        description="paoflow2J: Calculate exchange parameters J from PAOflow Hamiltonians using magnetic force theorem. "
                    "PAOflow uses an unusual file naming convention where spin components are appended AFTER the .dat extension "
                    "(e.g., hamiltonian.dat_0 for spin-up, hamiltonian.dat_1 for spin-down).",
        epilog="Examples:\n"
               "  # Basic usage with PAOflow output in current directory:\n"
               "  paoflow2J.py --efermi 18.35 --elements Fe --posfile POSCAR\n\n"
               "  # With custom Hamiltonian file locations:\n"
               "  paoflow2J.py --efermi 18.35 --elements Fe --posfile POSCAR \\\n"
               "    --prefix_up output/hamiltonian.dat_0 --prefix_down output/hamiltonian.dat_1\n\n"
               "  # Non-collinear calculation:\n"
               "  paoflow2J.py --efermi 18.35 --elements Fe --posfile POSCAR --spinor \\\n"
               "    --prefix_spinor hamiltonian.dat\n",
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    
    # PAOflow-specific arguments
    parser.add_argument(
        "--path",
        help="Path to directory containing PAOflow Hamiltonian files (.dat files)",
        default="./",
        type=str,
    )
    parser.add_argument(
        "--posfile",
        help="Name of the structure file (POSCAR, CIF, etc.). REQUIRED for PAOflow since "
             "Hamiltonian files don't contain atomic positions",
        default=None,
        type=str,
    )
    parser.add_argument(
        "--prefix_spinor",
        help="Hamiltonian file for non-collinear/SOC calculation. PAOflow names files as 'hamiltonian.dat'. "
             "Provide the complete filename or just the prefix (e.g., 'hamiltonian' or 'hamiltonian.dat')",
        default="hamiltonian.dat",
        type=str,
    )
    parser.add_argument(
        "--prefix_up",
        help="Hamiltonian file for spin-up channel. PAOflow appends spin index AFTER .dat extension, "
             "e.g., 'hamiltonian.dat_0'. Provide complete filename like 'output/hamiltonian.dat_0' or "
             "'hamiltonian.dat_0' (default: 'hamiltonian.dat_0')",
        default="hamiltonian.dat_0",
        type=str,
    )
    parser.add_argument(
        "--prefix_down",
        help="Hamiltonian file for spin-down channel. PAOflow appends spin index AFTER .dat extension, "
             "e.g., 'hamiltonian.dat_1'. Provide complete filename like 'output/hamiltonian.dat_1' or "
             "'hamiltonian.dat_1' (default: 'hamiltonian.dat_1')",
        default="hamiltonian.dat_1",
        type=str,
    )
    
    # Add common TB2J exchange calculation arguments
    add_exchange_args_to_parser(parser)
    
    args = parser.parse_args()
    
    # Validate required arguments
    if args.efermi is None:
        print("ERROR: Please provide the Fermi energy using --efermi <value>")
        print("Example: paoflow2J.py --efermi 6.15 --elements Mn Fe")
        sys.exit(1)
    
    if args.elements is None and args.index_magnetic_atoms is None:
        print("ERROR: Please specify magnetic atoms using --elements or --index_magnetic_atoms")
        print("Example: paoflow2J.py --efermi 6.15 --elements Mn Fe")
        print("   or:   paoflow2J.py --efermi 6.15 --index_magnetic_atoms 1 2 3")
        sys.exit(1)
    
    # Import here to avoid issues with other interfaces
    from TB2J.interfaces.paoflow_interface import gen_exchange_paoflow
    
    # Convert argument namespace to dictionary for Manager
    kwargs = parser_argument_to_dict(args)
    
    # Adjust index_magnetic_atoms to 0-based indexing if provided
    if args.index_magnetic_atoms is not None:
        kwargs['index_magnetic_atoms'] = [i - 1 for i in args.index_magnetic_atoms]
    
    print("\n" + "="*70)
    print("PAOflow to TB2J Interface")
    print("="*70)
    print(f"Reading PAOflow data from: {args.path}")
    if args.spinor:
        print(f"Non-collinear calculation with prefix: {args.prefix_spinor}")
    else:
        print(f"Collinear calculation with prefixes: {args.prefix_up}, {args.prefix_down}")
    print("="*70 + "\n")
    
    # Create manager and run calculation
    try:
        manager = gen_exchange_paoflow(
            path=args.path,
            prefix_up=args.prefix_up,
            prefix_dn=args.prefix_down,
            prefix_SOC=args.prefix_spinor,
            colinear=not args.spinor,
            posfile=args.posfile,
            **kwargs,
        )
        print("\n" + "="*70)
        print("Calculation completed successfully!")
        print("="*70)
        return manager
    except Exception as e:
        print("\n" + "="*70)
        print("ERROR: Calculation failed!")
        print("="*70)
        print(f"Error message: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    run_paoflow2J()
