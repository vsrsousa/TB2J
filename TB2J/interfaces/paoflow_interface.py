"""
Interface between PAOflow and TB2J.

PAOflow is an alternative to Wannier90 for constructing tight-binding
Hamiltonians from DFT calculations using projections of atomic orbitals (PAO).

PAOflow outputs Hamiltonians in Wannier90-compatible format (_hr.dat files),
which can be read using TB2J's Wannier90 interface. This module provides
a convenient wrapper for PAOflow users.

Note: PAOflow doesn't generate _centres.xyz files, so orbital positions are
assigned based on atomic positions.
"""

import os
import numpy as np
from ase.io import read

from TB2J.myTB import MyTB
from TB2J.utils import auto_assign_basis_name
from TB2J.wannier import parse_ham

from .manager import Manager


class PAOflowManager(Manager):
    """
    Manager class for interfacing PAOflow Hamiltonians with TB2J.
    
    PAOflow uses `write_Hamiltonian()` to output _hr.dat files in Wannier90
    format. This class handles the case where PAOflow doesn't generate
    _centres.xyz files (orbital positions are assigned from atomic positions).
    
    Note: Since PAOflow uses atomic orbital projections, orbitals are centered
    on atoms. TB2J will automatically assign orbital positions.
    """
    
    def __init__(
        self,
        path,
        prefix_up="hamiltonian_0",
        prefix_dn="hamiltonian_1", 
        prefix_SOC="hamiltonian",
        colinear=True,
        posfile=None,
        **kwargs,
    ):
        """
        Initialize PAOflowManager.
        
        Parameters
        ----------
        path : str
            Directory containing PAOflow output files
        prefix_up : str
            Prefix for spin-up _hr.dat file (default: "hamiltonian_0" for collinear)
        prefix_dn : str
            Prefix for spin-down _hr.dat file (default: "hamiltonian_1" for collinear)
        prefix_SOC : str
            Prefix for non-collinear/SOC _hr.dat file (default: "hamiltonian")
        colinear : bool
            Whether calculation is collinear (True) or non-collinear (False)
        posfile : str, optional
            Name of structure file (POSCAR, CIF, etc.) - REQUIRED for PAOflow
        **kwargs : dict
            Additional arguments passed to Manager
        """
        # Prepare atoms
        atoms = self.prepare_atoms(path, posfile, prefix_up, prefix_SOC, colinear)
        output_path = kwargs.get("output_path", "TB2J_results")
        
        # Prepare models and basis
        if colinear:
            tbmodels, basis = self.prepare_model_colinear(
                path, prefix_up, prefix_dn, atoms, output_path=output_path
            )
        else:
            tbmodels, basis = self.prepare_model_ncl(
                path, prefix_SOC, atoms, output_path=output_path
            )
        
        description = self.description(path, prefix_up, prefix_dn, prefix_SOC, colinear)
        kwargs["description"] = description
        
        super().__init__(atoms, tbmodels, basis, colinear=colinear, **kwargs)
    
    def prepare_atoms(self, path, posfile, prefix_up, prefix_SOC, colinear=True):
        """
        Prepare atomic structure.
        
        For PAOflow, a structure file is required since _hr.dat files don't
        contain atomic positions.
        """
        if posfile is None:
            raise ValueError(
                "PAOflow requires a structure file. Please provide --posfile argument.\n"
                "PAOflow's _hr.dat files don't contain atomic positions, unlike Wannier90 .win files."
            )
        
        fname = os.path.join(path, posfile)
        print(f"Reading atomic structure from file {fname}.")
        atoms = read(fname)
        return atoms
    
    def prepare_model_colinear(self, path, prefix_up, prefix_dn, atoms, output_path):
        """
        Prepare tight-binding models for collinear spin calculation.
        
        PAOflow doesn't generate _centres.xyz files, so we create orbital
        positions based on the number of orbitals and atomic positions.
        
        Note: PAOflow names Hamiltonian files as prefix.dat_0 and prefix.dat_1
        (with spin component after the .dat extension).
        """
        print("Reading PAOflow hamiltonian: spin up.")
        # PAOflow uses unusual naming: hamiltonian.dat_0, hamiltonian.dat_1
        # So don't append .dat if prefix already looks like a complete filename
        if not (prefix_up.endswith('.dat_0') or prefix_up.endswith('.dat_1') or prefix_up.endswith('.dat')):
            hr_fname_up = os.path.join(path, prefix_up + ".dat")
        else:
            hr_fname_up = os.path.join(path, prefix_up)
        nbasis_up, data_up, R_degens_up = parse_ham(fname=hr_fname_up)
        
        print("Reading PAOflow hamiltonian: spin down.")
        if not (prefix_dn.endswith('.dat_0') or prefix_dn.endswith('.dat_1') or prefix_dn.endswith('.dat')):
            hr_fname_dn = os.path.join(path, prefix_dn + ".dat")
        else:
            hr_fname_dn = os.path.join(path, prefix_dn)
        nbasis_dn, data_dn, R_degens_dn = parse_ham(fname=hr_fname_dn)
        
        # Create orbital positions based on atomic positions
        # PAOflow uses atomic orbital projections, so orbitals are centered on atoms
        positions = self._create_orbital_positions(nbasis_up, atoms)
        
        cell = atoms.get_cell()
        xred = cell.scaled_positions(positions)
        
        # Create TB models
        tbmodel_up = MyTB(
            nbasis=nbasis_up,
            data=data_up,
            R_degens=R_degens_up,
            positions=xred
        )
        tbmodel_up.set_atoms(atoms)
        
        tbmodel_dn = MyTB(
            nbasis=nbasis_dn,
            data=data_dn,
            R_degens=R_degens_dn,
            positions=xred
        )
        tbmodel_dn.set_atoms(atoms)
        
        # Assign basis names
        basis, _ = auto_assign_basis_name(
            xred,
            atoms,
            write_basis_file=os.path.join(output_path, "assigned_basis.txt"),
        )
        
        tbmodels = (tbmodel_up, tbmodel_dn)
        return tbmodels, basis
    
    def prepare_model_ncl(self, path, prefix_SOC, atoms, output_path):
        """
        Prepare tight-binding model for non-collinear spin calculation.
        
        Note: PAOflow may name files as prefix.dat or with other extensions.
        """
        print("Reading PAOflow hamiltonian: non-collinear spin.")
        # PAOflow uses unusual naming convention - don't append .dat if already present
        if not prefix_SOC.endswith('.dat'):
            hr_fname = os.path.join(path, prefix_SOC + ".dat")
        else:
            hr_fname = os.path.join(path, prefix_SOC)
        nbasis, data, R_degens = parse_ham(fname=hr_fname)
        
        # Create orbital positions based on atomic positions
        positions = self._create_orbital_positions(nbasis, atoms)
        
        cell = atoms.get_cell()
        xred = cell.scaled_positions(positions)
        
        # Create TB model
        tbmodel = MyTB(
            nbasis=nbasis,
            data=data,
            R_degens=R_degens,
            positions=xred,
            nspin=2  # Non-collinear has spinor structure
        )
        tbmodel.set_atoms(atoms)
        
        # Assign basis names
        basis, _ = auto_assign_basis_name(
            xred,
            atoms,
            write_basis_file=os.path.join(output_path, "assigned_basis.txt"),
        )
        
        return tbmodel, basis
    
    def _create_orbital_positions(self, nbasis, atoms):
        """
        Create orbital positions based on atomic positions.
        
        PAOflow uses atomic orbital projections, so orbitals should be
        centered on atoms. This method distributes orbitals evenly among atoms.
        
        Parameters
        ----------
        nbasis : int
            Number of basis functions (orbitals)
        atoms : Atoms
            ASE Atoms object with atomic positions
            
        Returns
        -------
        positions : ndarray
            Cartesian coordinates of orbital centers (shape: nbasis x 3)
        """
        natoms = len(atoms)
        atomic_pos = atoms.get_positions()
        
        # Distribute orbitals among atoms
        norb_per_atom = nbasis // natoms
        remainder = nbasis % natoms
        
        positions = []
        for i in range(natoms):
            # Each atom gets norb_per_atom orbitals
            # First 'remainder' atoms get one extra orbital
            norb = norb_per_atom + (1 if i < remainder else 0)
            for _ in range(norb):
                positions.append(atomic_pos[i])
        
        return np.array(positions)
    
    def description(self, path, prefix_up, prefix_dn, prefix_SOC, colinear=True):
        """
        Generate description text for the calculation.
        """
        if colinear:
            description = f""" Input from collinear PAOflow data.
 Tight binding data from {path}.
 PAOflow output files: {prefix_up}.dat and {prefix_dn}.dat (Wannier90 _hr.dat format).
 PAOflow generates Hamiltonians using projections of atomic orbitals.
 Note: Orbital positions assigned based on atomic positions (PAOflow doesn't generate _centres.xyz).
Warning: Please check if the noise level of the PAOflow Hamiltonian is much smaller than the exchange values.
\n"""
        else:
            description = f""" Input from non-collinear PAOflow data.
Tight binding data from {path}.
PAOflow output file: {prefix_SOC}.dat (Wannier90 _hr.dat format).
PAOflow generates Hamiltonians using projections of atomic orbitals.
Note: Orbital positions assigned based on atomic positions (PAOflow doesn't generate _centres.xyz).
Warning: Please check if the noise level of the PAOflow Hamiltonian is much smaller than the exchange values.
The DMI component parallel to the spin orientation, the Jani which has the component of that orientation should be disregarded
e.g. if the spins are along z, the xz, yz, zz, zx, zy components and the z component of DMI.
If you need these components, try to do three calculations with spin along x, y, z, or use structure with z rotated to x, y and z. 
And then use TB2J_merge.py to get the full set of parameters.\n"""
        
        return description


# Convenience function for backward compatibility
gen_exchange_paoflow = PAOflowManager
