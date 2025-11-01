"""
Parser for PAOflow Hamiltonian files.

PAOflow outputs tight-binding Hamiltonians in HDF5 format.
This module provides functions to read these files and convert them
to TB2J's internal format.
"""

import numpy as np
import h5py
from ase import Atoms
from ase.geometry import Cell


def parse_paoflow_hdf5(fname):
    """
    Parse PAOflow HDF5 file containing Hamiltonian data.
    
    Parameters
    ----------
    fname : str
        Path to PAOflow HDF5 file
        
    Returns
    -------
    data : dict
        Dictionary mapping R vectors (tuples) to Hamiltonian matrices
    R_degens : dict
        Degeneracy of each R vector
    positions : ndarray
        Orbital positions in reduced coordinates
    atoms : Atoms
        ASE Atoms object containing structure information
    nbasis : int
        Number of basis functions
    """
    with h5py.File(fname, 'r') as f:
        # Read Hamiltonian data
        # PAOflow typically stores H(R) in format:
        # - 'hamiltonian' or 'H_R': Hamiltonian matrices
        # - 'R_vectors' or 'Rlist': R lattice vectors
        # - 'orbital_positions' or 'wann_centers': orbital centers
        # - 'cell' or 'lattice': lattice vectors
        # - 'positions': atomic positions
        # - 'numbers': atomic numbers
        
        # Try different possible dataset names
        if 'hamiltonian' in f:
            H_R = f['hamiltonian'][:]
        elif 'H_R' in f:
            H_R = f['H_R'][:]
        else:
            raise KeyError("Could not find Hamiltonian data in HDF5 file. "
                          "Expected 'hamiltonian' or 'H_R' dataset.")
        
        if 'R_vectors' in f:
            R_vectors = f['R_vectors'][:]
        elif 'Rlist' in f:
            R_vectors = f['Rlist'][:]
        else:
            raise KeyError("Could not find R vectors in HDF5 file. "
                          "Expected 'R_vectors' or 'Rlist' dataset.")
        
        if 'orbital_positions' in f:
            orbital_pos = f['orbital_positions'][:]
        elif 'wann_centers' in f:
            orbital_pos = f['wann_centers'][:]
        else:
            raise KeyError("Could not find orbital positions in HDF5 file. "
                          "Expected 'orbital_positions' or 'wann_centers' dataset.")
        
        if 'cell' in f:
            cell = f['cell'][:]
        elif 'lattice' in f:
            cell = f['lattice'][:]
        else:
            raise KeyError("Could not find cell information in HDF5 file. "
                          "Expected 'cell' or 'lattice' dataset.")
        
        if 'positions' in f:
            atomic_pos = f['positions'][:]
        else:
            raise KeyError("Could not find atomic positions in HDF5 file. "
                          "Expected 'positions' dataset.")
        
        if 'numbers' in f:
            atomic_numbers = f['numbers'][:]
        elif 'atomic_numbers' in f:
            atomic_numbers = f['atomic_numbers'][:]
        else:
            raise KeyError("Could not find atomic numbers in HDF5 file. "
                          "Expected 'numbers' or 'atomic_numbers' dataset.")
        
        # Optional: degeneracy weights
        if 'R_degeneracy' in f:
            R_degens_array = f['R_degeneracy'][:]
        elif 'weights' in f:
            R_degens_array = f['weights'][:]
        else:
            # Default to 1 for all R vectors
            R_degens_array = np.ones(len(R_vectors))
    
    # Convert to TB2J format
    nbasis = H_R.shape[1]  # Assuming H_R has shape (nR, nbasis, nbasis)
    
    # Create data dictionary mapping R tuples to matrices
    data = {}
    R_degens = {}
    for iR, R in enumerate(R_vectors):
        R_tuple = tuple(R.astype(int))
        data[R_tuple] = H_R[iR]
        R_degens[R_tuple] = R_degens_array[iR]
    
    # Create ASE Atoms object
    atoms = Atoms(numbers=atomic_numbers, 
                  positions=atomic_pos, 
                  cell=cell, 
                  pbc=True)
    
    # Convert orbital positions to reduced coordinates
    cell_obj = Cell(cell)
    positions = cell_obj.scaled_positions(orbital_pos)
    
    return data, R_degens, positions, atoms, nbasis


def parse_paoflow_files(path, prefix):
    """
    Parse PAOflow output files for a given prefix.
    
    Parameters
    ----------
    path : str
        Directory containing PAOflow output files
    prefix : str
        Prefix for PAOflow files (e.g., 'paoflow', 'paoflow_up', 'paoflow_dn')
        
    Returns
    -------
    data : dict
        Hamiltonian data
    R_degens : dict
        R vector degeneracies
    positions : ndarray
        Orbital positions
    atoms : Atoms
        Atomic structure
    nbasis : int
        Number of basis functions
    """
    import os
    
    # Try different possible file extensions
    possible_files = [
        f"{prefix}.hdf5",
        f"{prefix}.h5",
        f"{prefix}_hr.hdf5",
        f"{prefix}_hr.h5",
    ]
    
    fname = None
    for possible_file in possible_files:
        test_path = os.path.join(path, possible_file)
        if os.path.exists(test_path):
            fname = test_path
            break
    
    if fname is None:
        raise FileNotFoundError(
            f"Could not find PAOflow HDF5 file with prefix '{prefix}' in {path}. "
            f"Tried: {', '.join(possible_files)}"
        )
    
    return parse_paoflow_hdf5(fname)
