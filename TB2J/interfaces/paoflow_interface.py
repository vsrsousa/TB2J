"""
Interface between PAOflow and TB2J.

PAOflow is an alternative to Wannier90 for constructing tight-binding
Hamiltonians from DFT calculations using projections of atomic orbitals (PAO).
"""

import os
from ase.io import read

from TB2J.myTB import MyTB
from TB2J.utils import auto_assign_basis_name

from .manager import Manager
from .paoflow_parser import parse_paoflow_files


class PAOflowManager(Manager):
    """
    Manager class for interfacing PAOflow Hamiltonians with TB2J.
    
    This class handles both collinear and non-collinear spin calculations
    from PAOflow output files.
    """
    
    def __init__(
        self,
        path,
        prefix_up,
        prefix_dn,
        prefix_SOC,
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
            Prefix for spin-up channel files (collinear case)
        prefix_dn : str
            Prefix for spin-down channel files (collinear case)
        prefix_SOC : str
            Prefix for non-collinear/SOC files
        colinear : bool
            Whether calculation is collinear (True) or non-collinear (False)
        posfile : str, optional
            Name of structure file (POSCAR, etc.). If None, structure is read
            from PAOflow HDF5 file
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
        
        First tries to read from posfile if provided, otherwise reads from
        PAOflow HDF5 file.
        """
        if posfile is not None:
            fname = os.path.join(path, posfile)
            try:
                print(f"Reading atomic structure from file {fname}.")
                atoms = read(fname)
                return atoms
            except Exception as e:
                print(f"Warning: Could not read atomic structure from {fname}: {e}")
        
        # Read from PAOflow HDF5 file
        print("Reading atomic structure from PAOflow HDF5 file.")
        if colinear:
            prefix = prefix_up
        else:
            prefix = prefix_SOC
        
        _, _, _, atoms, _ = parse_paoflow_files(path, prefix)
        return atoms
    
    def prepare_model_colinear(self, path, prefix_up, prefix_dn, atoms, output_path):
        """
        Prepare tight-binding models for collinear spin calculation.
        
        Parameters
        ----------
        path : str
            Directory containing PAOflow files
        prefix_up : str
            Prefix for spin-up files
        prefix_dn : str
            Prefix for spin-down files
        atoms : Atoms
            ASE Atoms object
        output_path : str
            Output directory for basis assignment file
            
        Returns
        -------
        tbmodels : tuple
            Tuple of (tbmodel_up, tbmodel_dn)
        basis : list
            List of basis function labels
        """
        print("Reading PAOflow hamiltonian: spin up.")
        data_up, R_degens_up, positions_up, _, nbasis_up = parse_paoflow_files(
            path, prefix_up
        )
        tbmodel_up = MyTB(
            nbasis=nbasis_up,
            data=data_up,
            R_degens=R_degens_up,
            positions=positions_up
        )
        tbmodel_up.set_atoms(atoms)
        
        print("Reading PAOflow hamiltonian: spin down.")
        data_dn, R_degens_dn, positions_dn, _, nbasis_dn = parse_paoflow_files(
            path, prefix_dn
        )
        tbmodel_dn = MyTB(
            nbasis=nbasis_dn,
            data=data_dn,
            R_degens=R_degens_dn,
            positions=positions_dn
        )
        tbmodel_dn.set_atoms(atoms)
        
        # Assign basis names
        basis, _ = auto_assign_basis_name(
            positions_up,
            atoms,
            write_basis_file=os.path.join(output_path, "assigned_basis.txt"),
        )
        
        tbmodels = (tbmodel_up, tbmodel_dn)
        return tbmodels, basis
    
    def prepare_model_ncl(self, path, prefix_SOC, atoms, output_path):
        """
        Prepare tight-binding model for non-collinear spin calculation.
        
        Parameters
        ----------
        path : str
            Directory containing PAOflow files
        prefix_SOC : str
            Prefix for non-collinear/SOC files
        atoms : Atoms
            ASE Atoms object
        output_path : str
            Output directory for basis assignment file
            
        Returns
        -------
        tbmodel : MyTB
            Tight-binding model
        basis : list
            List of basis function labels
        """
        print("Reading PAOflow hamiltonian: non-collinear spin.")
        data, R_degens, positions, _, nbasis = parse_paoflow_files(path, prefix_SOC)
        
        tbmodel = MyTB(
            nbasis=nbasis,
            data=data,
            R_degens=R_degens,
            positions=positions,
            nspin=2  # Non-collinear has spinor (2-component) structure
        )
        tbmodel.set_atoms(atoms)
        
        # Assign basis names
        basis, _ = auto_assign_basis_name(
            positions,
            atoms,
            write_basis_file=os.path.join(output_path, "assigned_basis.txt"),
        )
        
        return tbmodel, basis
    
    def description(self, path, prefix_up, prefix_dn, prefix_SOC, colinear=True):
        """
        Generate description text for the calculation.
        """
        if colinear:
            description = f""" Input from collinear PAOflow data.
 Tight binding data from {path}.
 Prefix of PAOflow HDF5 files: {prefix_up} and {prefix_dn}.
Warning: Please check if the noise level of the PAOflow Hamiltonian is much smaller than the exchange values.
\n"""
        else:
            description = f""" Input from non-collinear PAOflow data.
Tight binding data from {path}.
Prefix of PAOflow HDF5 files: {prefix_SOC}.
Warning: Please check if the noise level of the PAOflow Hamiltonian is much smaller than the exchange values.
The DMI component parallel to the spin orientation, the Jani which has the component of that orientation should be disregarded
e.g. if the spins are along z, the xz, yz, zz, zx, zy components and the z component of DMI.
If you need these components, try to do three calculations with spin along x, y, z, or use structure with z rotated to x, y and z. 
And then use TB2J_merge.py to get the full set of parameters.\n"""
        
        return description


# Convenience function for backward compatibility
gen_exchange_paoflow = PAOflowManager
