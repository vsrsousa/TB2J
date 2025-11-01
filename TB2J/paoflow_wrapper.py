"""
Wrapper for reading PAOFLOW Hamiltonian output for TB2J calculations.

PAOFLOW outputs tight-binding Hamiltonians in real space (HRs) in a format
compatible with Wannier90's hr.dat format. This wrapper facilitates reading
PAOFLOW data into TB2J.
"""

import os
import numpy as np
from ase.atoms import Atoms
from TB2J.myTB import MyTB
from TB2J.wannier import parse_ham
from TB2J.utils import auto_assign_basis_name


class PAOFLOWWrapper:
    """
    Wrapper class for reading PAOFLOW Hamiltonian output.
    
    PAOFLOW can write Hamiltonians in Wannier90-compatible format using
    the write_Hamiltonian() method. This wrapper reads those files and
    provides them in a format compatible with TB2J's exchange calculations.
    """
    
    @staticmethod
    def read_paoflow_hr(hr_fname, atoms, positions_fname=None):
        """
        Read PAOFLOW Hamiltonian from hr.dat format file.
        
        Parameters
        ----------
        hr_fname : str
            Path to the Hamiltonian file in Wannier90 hr.dat format
            (typically 'hamiltonian.dat' from PAOFLOW's write_Hamiltonian)
        atoms : ase.Atoms
            ASE Atoms object with structure information
        positions_fname : str, optional
            Path to file containing orbital positions. If not provided,
            positions will be auto-assigned based on atomic positions.
            
        Returns
        -------
        MyTB
            TB2J tight-binding model object
        """
        # Read Hamiltonian using the existing Wannier90 parser
        nbasis, data, R_degens = parse_ham(fname=hr_fname)
        
        # Handle orbital positions
        if positions_fname is not None and os.path.exists(positions_fname):
            # Read positions from file if provided
            positions = np.loadtxt(positions_fname)
            if positions.shape[0] != nbasis:
                raise ValueError(
                    f"Number of positions ({positions.shape[0]}) does not match "
                    f"number of orbitals ({nbasis})"
                )
            cell = atoms.get_cell()
            xred = cell.scaled_positions(positions)
        else:
            # Auto-assign positions based on atomic structure
            # Create positions array with one orbital per atom (can be refined)
            natoms = len(atoms)
            norb_per_atom = nbasis // natoms
            
            if nbasis % natoms != 0:
                # If orbitals don't divide evenly, create approximate positions
                xred = np.zeros((nbasis, 3))
                for i in range(nbasis):
                    atom_idx = i % natoms
                    xred[i] = atoms.get_scaled_positions()[atom_idx]
            else:
                # Replicate atomic positions for each orbital
                xred = np.repeat(atoms.get_scaled_positions(), norb_per_atom, axis=0)
        
        # Create TB2J model
        ind, positions = auto_assign_basis_name(xred, atoms)
        m = MyTB(nbasis=nbasis, data=data, positions=xred, R_degens=R_degens)
        nm = m.shift_position(positions)
        nm.set_atoms(atoms)
        
        return nm
    
    @staticmethod
    def read_paoflow_collinear(hr_up_fname, hr_dn_fname, atoms, positions_fname=None):
        """
        Read PAOFLOW collinear spin Hamiltonians.
        
        For collinear spin calculations, PAOFLOW writes separate files for
        spin up and spin down channels (hamiltonian.dat_0 and hamiltonian.dat_1).
        
        Parameters
        ----------
        hr_up_fname : str
            Path to spin-up Hamiltonian file
        hr_dn_fname : str
            Path to spin-down Hamiltonian file
        atoms : ase.Atoms
            ASE Atoms object with structure information
        positions_fname : str, optional
            Path to file containing orbital positions
            
        Returns
        -------
        tuple of MyTB
            (tbmodel_up, tbmodel_dn) - TB2J tight-binding model objects for
            spin up and spin down channels
        """
        tbmodel_up = PAOFLOWWrapper.read_paoflow_hr(
            hr_up_fname, atoms, positions_fname
        )
        tbmodel_dn = PAOFLOWWrapper.read_paoflow_hr(
            hr_dn_fname, atoms, positions_fname
        )
        
        return tbmodel_up, tbmodel_dn
