"""
Interface between PAOflow and TB2J.

PAOflow is an alternative to Wannier90 for constructing tight-binding
Hamiltonians from DFT calculations using projections of atomic orbitals (PAO).

PAOflow outputs Hamiltonians in Wannier90-compatible format (_hr.dat files),
which can be read using TB2J's Wannier90 interface. This module provides
a convenient wrapper for PAOflow users.
"""

from .wannier90_interface import WannierManager


class PAOflowManager(WannierManager):
    """
    Manager class for interfacing PAOflow Hamiltonians with TB2J.
    
    PAOflow uses `write_Hamiltonian()` to output _hr.dat files in Wannier90
    format. This class is a thin wrapper around WannierManager that provides
    PAOflow-specific documentation and defaults.
    
    Note: PAOflow generates Wannier90-compatible files, so you can also use
    the standard `wann2J.py` command with PAOflow output.
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
            Name of structure file (POSCAR, CIF, etc.)
        **kwargs : dict
            Additional arguments passed to WannierManager
        """
        # Call parent WannierManager with PAOflow-specific settings
        super().__init__(
            path=path,
            prefix_up=prefix_up,
            prefix_dn=prefix_dn,
            prefix_SOC=prefix_SOC,
            colinear=colinear,
            posfile=posfile,
            wannier_type="paoflow",
            **kwargs,
        )
    
    def description(self, path, prefix_up, prefix_dn, prefix_SOC, colinear=True):
        """
        Generate description text for the calculation.
        """
        if colinear:
            description = f""" Input from collinear PAOflow data.
 Tight binding data from {path}.
 PAOflow output files: {prefix_up}.dat and {prefix_dn}.dat (Wannier90 _hr.dat format).
 PAOflow generates Hamiltonians using projections of atomic orbitals.
Warning: Please check if the noise level of the PAOflow Hamiltonian is much smaller than the exchange values.
\n"""
        else:
            description = f""" Input from non-collinear PAOflow data.
Tight binding data from {path}.
PAOflow output file: {prefix_SOC}.dat (Wannier90 _hr.dat format).
PAOflow generates Hamiltonians using projections of atomic orbitals.
Warning: Please check if the noise level of the PAOflow Hamiltonian is much smaller than the exchange values.
The DMI component parallel to the spin orientation, the Jani which has the component of that orientation should be disregarded
e.g. if the spins are along z, the xz, yz, zz, zx, zy components and the z component of DMI.
If you need these components, try to do three calculations with spin along x, y, z, or use structure with z rotated to x, y and z. 
And then use TB2J_merge.py to get the full set of parameters.\n"""
        
        return description


# Convenience function for backward compatibility
gen_exchange_paoflow = PAOflowManager
