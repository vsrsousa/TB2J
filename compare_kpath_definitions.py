#!/usr/bin/env python
"""
Compare ASE vs SeekPath k-path generation definitions
"""
import numpy as np
from ase.io import read
from ase.dft.kpoints import BandPath
import seekpath

# Read structure
atoms = read('/Users/vinicius/Documents/cifs/Gd.cif')
cell = atoms.get_cell()
positions_frac = atoms.get_scaled_positions()
atomic_numbers = atoms.get_atomic_numbers()

print("=" * 80)
print("ASE k-path DEFINITION")
print("=" * 80)

# Get ASE bandpath
bp_ase = BandPath(path='GMK', cell=cell)
print(f"\nASE.dft.kpoints.BandPath documentation:")
print(f"- Uses simple high-symmetry points lookup")
print(f"- Path argument: string like 'GMK' or list of points")
print(f"- Looks up special points from standard tables (Setyawan & Curtarolo 2010)")
print(f"- Special points database: ase/data/spacegroups/")
print(f"- For HCP (Gd): uses standard hexagonal special points")

print(f"\nASE special_points for HCP:")
print(bp_ase.special_points)

print(f"\nASE bandpath attributes:")
print(f"  .path: {bp_ase.path}")
print(f"  .kpts: shape {bp_ase.kpts.shape}")

# Get auto path from ASE
bp_auto = BandPath(path=None, cell=cell)
print(f"\nASE auto-path (path=None):")
print(f"  .path: {bp_auto.path}")
print(f"  .special_points: {bp_auto.special_points}")

print("\n" + "=" * 80)
print("SeekPath k-path DEFINITION")
print("=" * 80)

structure = (cell.tolist(), positions_frac.tolist(), atomic_numbers.tolist())
result = seekpath.get_path(structure)

print(f"\nSeekPath algorithm: HPKOT (Hinuma, Pizzi, Kumagai, Tanaka, Otani)")
print(f"Paper: Physical Review B 100, 104302 (2019)")
print(f"https://journals.aps.org/prb/abstract/10.1103/PhysRevB.100.104302")

print(f"\nSeekPath process:")
print(f"1. Standardizes input structure to conventional cell")
print(f"2. Detects space group and Bravais lattice (here: {result['bravais_lattice']})")
print(f"3. Uses crystal symmetry to determine high-symmetry directions")
print(f"4. Generates k-path respecting all symmetries")
print(f"5. Returns primitive cell coordinates")

print(f"\nSeekPath space group info:")
print(f"  Space group number: {result['spacegroup_number']}")
print(f"  Space group (International): {result['spacegroup_international']}")
print(f"  Bravais lattice: {result['bravais_lattice']}")
print(f"  Bravais lattice extended: {result['bravais_lattice_extended']}")

print(f"\nSeekPath special points (in primitive cell):")
for name, coord in result['point_coords'].items():
    print(f"  {name}: {coord}")

print(f"\nSeekPath path (segments):")
for seg in result['path']:
    print(f"  {seg[0]} → {seg[1]}")

print("\n" + "=" * 80)
print("KEY DIFFERENCES")
print("=" * 80)

print(f"""
ASE:
  • Uses pre-computed special points tables (Setyawan & Curtarolo 2010)
  • Hardcoded paths for each structure type
  • Does NOT detect space group automatically
  • Returns fewer high-symmetry points (typically 4-6 for HCP)
  • Simple, fast, but limited in completeness
  • Example Gd HCP: GMKGALHA,LM,KH (6 points, 4 segments)

SeekPath (HPKOT):
  • Automatically determines space group from atomic structure
  • Uses rigorous symmetry analysis (group theory)
  • Adapts path based on crystal's actual symmetries
  • Finds more complete set of high-symmetry points (typically 8-12 for HCP)
  • Can discover symmetries ASE's tables miss
  • Example Gd HCP: 10 points, 11 segments (GAMMA, Y, T, Z, S, R, SIGMA_0, C_0, A_0, E_0)

INPUT SPACE:
  ASE: Original cell as-is
  SeekPath: Transforms to conventional/standardized cell for analysis

OUTPUT SPACE:
  ASE: Returns k-points in original reciprocal space
  SeekPath: Returns k-points in primitive reciprocal space (requires transformation back if needed)
""")

print("\n" + "=" * 80)
print("TRANSFORMATION MATRIX (SeekPath input → primitive)")
print("=" * 80)
print(f"\nPrimitive transformation matrix P:")
print(result['primitive_transformation_matrix'])
print(f"\nInverse (to convert back to original):")
print(result['inverse_primitive_transformation_matrix'])
