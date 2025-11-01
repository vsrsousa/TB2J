#!/usr/bin/env python3
"""
paoflow2J: Script to calculate magnetic exchange parameters from PAOFLOW Hamiltonians.

PAOFLOW outputs tight-binding Hamiltonians in real space that are compatible
with Wannier90's hr.dat format. This script reads those Hamiltonians and
calculates exchange interactions using the magnetic force theorem.
"""

from TB2J.manager import gen_exchange_paoflow
from TB2J.versioninfo import print_license
import sys
import argparse
import os


def run_paoflow2J():
    print_license()
    parser = argparse.ArgumentParser(
        description="paoflow2J: Calculate exchange parameter J from PAOFLOW Hamiltonian using magnetic force theorem"
    )
    
    parser.add_argument(
        "--path",
        help="Path to directory containing PAOFLOW output files. Default: current directory",
        default="./",
        type=str,
    )
    
    parser.add_argument(
        "--hr_up",
        help="Filename of spin-up Hamiltonian. For collinear: 'hamiltonian.dat_0'. "
             "For non-collinear: 'hamiltonian.dat'. Default: hamiltonian.dat_0",
        default="hamiltonian.dat_0",
        type=str,
    )
    
    parser.add_argument(
        "--poscar",
        help="Filename of structure file (POSCAR, cif, xyz, etc.) readable by ASE. Default: POSCAR",
        default="POSCAR",
        type=str,
    )
    
    parser.add_argument(
        "--hr_dn",
        help="Filename of spin-down Hamiltonian for collinear calculations. "
             "Typically 'hamiltonian.dat_1'. Default: hamiltonian.dat_1",
        default="hamiltonian.dat_1",
        type=str,
    )
    
    parser.add_argument(
        "--positions_fname",
        help="Path to file with orbital positions (optional). Format: 3 columns (x, y, z) in Cartesian coordinates.",
        default=None,
        type=str,
    )
    
    parser.add_argument(
        "--elements",
        help="Magnetic elements to be considered in Heisenberg model (e.g., Fe Ni)",
        default=None,
        type=str,
        nargs="*",
    )
    
    parser.add_argument(
        "--non_colinear",
        help="Use non-collinear spin mode. By default, collinear mode is used.",
        action="store_true",
        default=False,
    )
    
    parser.add_argument(
        "--rcut",
        help="Cutoff of spin pair distance. Default: calculate all commensurate R points to the k mesh.",
        default=None,
        type=float,
    )
    
    parser.add_argument(
        "--efermi",
        help="Fermi energy in eV. Default: 0.0",
        default=0.0,
        type=float,
    )
    
    parser.add_argument(
        "--kmesh",
        help="k-mesh in the format kx ky kz. Default: 5 5 5",
        type=int,
        nargs="*",
        default=[5, 5, 5],
    )
    
    parser.add_argument(
        "--emin",
        help="Energy minimum below Fermi level (eV). Default: -14.0",
        type=float,
        default=-14.0,
    )
    
    parser.add_argument(
        "--emax",
        help="Energy maximum above Fermi level (eV). Default: 0.0",
        type=float,
        default=0.0,
    )
    
    parser.add_argument(
        "--nz",
        help="Number of steps for semicircle contour integration. Default: 100",
        default=100,
        type=int,
    )
    
    parser.add_argument(
        "--cutoff",
        help="Minimum J amplitude to write (eV). Default: 1e-5",
        default=1e-5,
        type=float,
    )
    
    parser.add_argument(
        "--exclude_orbs",
        help="Indices of orbitals to exclude from magnetic sites (counting from 0)",
        default=[],
        type=int,
        nargs="+",
    )
    
    parser.add_argument(
        "--np",
        help="Number of CPU cores for parallel processing. Default: 1",
        default=1,
        type=int,
    )
    
    parser.add_argument(
        "--use_cache",
        help="Use disk cache for wavefunctions and Hamiltonian to reduce memory usage",
        action="store_true",
        default=False,
    )
    
    parser.add_argument(
        "--description",
        help="Description of the calculation (e.g., XC functional, U values, magnetic state)",
        type=str,
        default="Calculated with TB2J from PAOFLOW Hamiltonian.",
    )
    
    parser.add_argument(
        "--orb_decomposition",
        default=False,
        action="store_true",
        help="Perform orbital decomposition in non-collinear mode",
    )
    
    parser.add_argument(
        "--output_path",
        help="Output directory path. Default: TB2J_results",
        type=str,
        default="TB2J_results",
    )
    
    parser.add_argument(
        "--write_dm",
        help="Write density matrix",
        action="store_true",
        default=False,
    )
    
    args = parser.parse_args()
    
    # Validate required arguments
    if args.elements is None:
        print("ERROR: Please specify magnetic elements using --elements (e.g., --elements Fe Ni)")
        sys.exit(1)
    
    # Construct full paths to files
    hr_up_path = os.path.join(args.path, args.hr_up)
    hr_dn_path = os.path.join(args.path, args.hr_dn)
    poscar_path = os.path.join(args.path, args.poscar)
    positions_path = os.path.join(args.path, args.positions_fname) if args.positions_fname else None
    
    # Call the main function
    gen_exchange_paoflow(
        hr_up=hr_up_path,
        poscar=poscar_path,
        colinear=not args.non_colinear,
        hr_dn=hr_dn_path,
        positions_fname=positions_path,
        magnetic_elements=args.elements,
        kmesh=args.kmesh,
        emin=args.emin,
        emax=args.emax,
        nz=args.nz,
        exclude_orbs=args.exclude_orbs,
        Rcut=args.rcut,
        np=args.np,
        efermi=args.efermi,
        use_cache=args.use_cache,
        output_path=args.output_path,
        orb_decomposition=args.orb_decomposition,
        write_density_matrix=args.write_dm,
        description=args.description,
    )


if __name__ == "__main__":
    run_paoflow2J()
