Use TB2J with PAOflow
=====================

This tutorial explains how to use TB2J with PAOflow to calculate magnetic exchange parameters.

What is PAOflow?
----------------

PAOflow (Projections of Atomic Orbitals flow) is an alternative to Wannier90 for constructing 
tight-binding Hamiltonians from DFT calculations. It uses projections of atomic orbitals (PAO) 
to generate localized basis functions and outputs the Hamiltonian in HDF5 format.

Key features of PAOflow:

* Generates tight-binding Hamiltonians from DFT calculations
* Uses HDF5 format for efficient data storage
* Supports both collinear and non-collinear spin calculations
* Alternative workflow to Wannier90

Installation
------------

To use PAOflow with TB2J, you need to install the optional PAOflow dependency:

.. code-block:: bash

   pip install TB2J[paoflow]

This will install ``h5py``, which is required to read PAOflow's HDF5 output files.

PAOflow Output Files
--------------------

PAOflow typically generates the following files that TB2J can read:

* ``paoflow_up.hdf5`` or ``paoflow_up.h5`` - Hamiltonian for spin-up channel (collinear)
* ``paoflow_dn.hdf5`` or ``paoflow_dn.h5`` - Hamiltonian for spin-down channel (collinear)
* ``paoflow_soc.hdf5`` or ``paoflow_soc.h5`` - Hamiltonian for non-collinear/SOC calculation

The HDF5 files should contain:

* ``hamiltonian`` or ``H_R`` - Hamiltonian matrices H(R) for each R vector
* ``R_vectors`` or ``Rlist`` - Real-space lattice vectors R
* ``orbital_positions`` or ``wann_centers`` - Positions of orbital centers
* ``cell`` or ``lattice`` - Lattice vectors of the unit cell
* ``positions`` - Atomic positions
* ``numbers`` or ``atomic_numbers`` - Atomic numbers
* ``R_degeneracy`` or ``weights`` (optional) - Degeneracy weights for R vectors

Running TB2J with PAOflow
--------------------------

Collinear Calculation
~~~~~~~~~~~~~~~~~~~~~

For a collinear spin-polarized calculation, you need two PAOflow output files 
(one for spin-up and one for spin-down):

.. code-block:: bash

   paoflow2J.py --path /path/to/paoflow/output \
                --prefix_up paoflow_up \
                --prefix_down paoflow_dn \
                --efermi 6.15 \
                --kmesh 5 5 5 \
                --elements Mn Fe \
                --emin -10.0 \
                --emax 0.0

Parameters:

* ``--path`` - Directory containing PAOflow HDF5 files
* ``--prefix_up`` - Prefix for spin-up HDF5 file (default: ``paoflow_up``)
* ``--prefix_down`` - Prefix for spin-down HDF5 file (default: ``paoflow_dn``)
* ``--efermi`` - Fermi energy in eV (required)
* ``--kmesh`` - k-point mesh for integration (default: 5 5 5)
* ``--elements`` - Magnetic elements (e.g., Mn Fe Ni)
* ``--emin`` - Lower energy limit relative to Fermi energy (in eV)
* ``--emax`` - Upper energy limit relative to Fermi energy (in eV, typically 0.0)

Non-collinear Calculation
~~~~~~~~~~~~~~~~~~~~~~~~~~

For non-collinear or spin-orbit coupling calculations:

.. code-block:: bash

   paoflow2J.py --path /path/to/paoflow/output \
                --spinor \
                --prefix_spinor paoflow_soc \
                --efermi 6.15 \
                --kmesh 5 5 5 \
                --elements Mn Fe \
                --emin -10.0 \
                --emax 0.0

The ``--spinor`` flag indicates a non-collinear calculation.

Structure File
~~~~~~~~~~~~~~

By default, TB2J reads the atomic structure from the PAOflow HDF5 file. 
However, you can optionally provide a separate structure file:

.. code-block:: bash

   paoflow2J.py --posfile POSCAR --path /path/to/paoflow/output ...

The ``--posfile`` parameter accepts any format supported by ASE (POSCAR, CIF, XYZ, etc.).

Complete Example
----------------

Here's a complete example for calculating exchange parameters in SrMnO\ :sub:`3`:

Step 1: Run PAOflow
~~~~~~~~~~~~~~~~~~~

First, generate the tight-binding Hamiltonian using PAOflow (refer to PAOflow documentation):

.. code-block:: python

   # Example PAOflow script (pseudocode)
   from paoflow import PAOflow
   
   pao = PAOflow(calculation_type='scf',
                 inputfile='scf.in',
                 output_format='hdf5')
   pao.projectability()
   pao.pao_hamiltonian()
   pao.write_hamiltonian('paoflow_up.hdf5')  # For spin-up
   pao.write_hamiltonian('paoflow_dn.hdf5')  # For spin-down

Step 2: Run TB2J
~~~~~~~~~~~~~~~~

Calculate the exchange parameters:

.. code-block:: bash

   paoflow2J.py --path ./ \
                --prefix_up paoflow_up \
                --prefix_down paoflow_dn \
                --efermi 6.15 \
                --kmesh 4 4 4 \
                --elements Mn \
                --emin -10.0 \
                --emax 0.0 \
                --rcut 10.0

This will create a ``TB2J_results`` directory containing:

.. code-block:: text

   TB2J_results/
   ├── exchange.txt          # Human-readable exchange parameters
   ├── assigned_basis.txt    # Basis function assignments
   ├── Multibinit/          # Input files for Multibinit
   ├── TomASD/              # Input files for ASD
   └── Vampire/             # Input files for Vampire

Output Files
------------

The main output file is ``exchange.txt``, which contains:

* Isotropic exchange parameters (J)
* Anisotropic exchange (Jani)
* Dzyaloshinskii-Moriya interaction (DMI)
* Distances and coordination information

Tips and Troubleshooting
-------------------------

1. **Check Hamiltonian Quality**: Ensure that the PAOflow Hamiltonian is well-converged 
   and the noise level is much smaller than the exchange energies.

2. **Orbital Localization**: Like Wannier90, the quality of results depends on having 
   well-localized orbitals. Check the spreads in PAOflow output.

3. **Energy Window**: Choose ``--emin`` and ``--emax`` to include all relevant orbitals 
   (typically d or f orbitals for magnetic atoms and p orbitals for ligands).

4. **k-mesh Convergence**: Test convergence with respect to the k-point mesh using 
   ``--kmesh``. A denser mesh may be needed for complex materials.

5. **File Format**: TB2J looks for files with extensions ``.hdf5`` or ``.h5``. Ensure 
   your PAOflow outputs use one of these extensions.

6. **Missing Dependencies**: If you get an import error for ``h5py``, install it:

   .. code-block:: bash

      pip install h5py

Comparison with Wannier90
--------------------------

PAOflow vs Wannier90:

+------------------------+------------------------+------------------------+
| Feature                | Wannier90              | PAOflow                |
+========================+========================+========================+
| File format            | Text files (_hr.dat)   | HDF5 (.h5, .hdf5)      |
+------------------------+------------------------+------------------------+
| Orbital construction   | Maximally localized WFs| Projected atomic orbs  |
+------------------------+------------------------+------------------------+
| Interface in TB2J      | ``wann2J.py``          | ``paoflow2J.py``       |
+------------------------+------------------------+------------------------+
| Typical use case       | General tight-binding  | Alternative workflow   |
+------------------------+------------------------+------------------------+

Both methods should give similar results if the Hamiltonians are well-converged.

Advanced Options
----------------

Index-based Magnetic Atom Selection
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Instead of specifying elements, you can select magnetic atoms by index:

.. code-block:: bash

   paoflow2J.py --index_magnetic_atoms 1 3 5 ...

(Note: indices are 1-based in the command-line interface)

Custom Output Path
~~~~~~~~~~~~~~~~~~

Specify a custom output directory:

.. code-block:: bash

   paoflow2J.py --output_path my_results ...

Orbital Decomposition
~~~~~~~~~~~~~~~~~~~~~

Analyze orbital contributions to exchange:

.. code-block:: bash

   paoflow2J.py --orb_decomposition ...

For more advanced options, see:

.. code-block:: bash

   paoflow2J.py --help

Further Reading
---------------

* PAOflow documentation: https://github.com/marcobn/PAOflow
* TB2J documentation: https://tb2j.readthedocs.io
* Wannier90 comparison: :doc:`wannier`
* Exchange parameter reference: :doc:`parameters`

References
----------

If you use PAOflow with TB2J, please cite:

* TB2J paper: [TB2J citation]
* PAOflow paper: Marco Buongiorno Nardelli, "PAOflow: A utility for constructing 
  Hamiltonians from projections of atomic orbitals", Comp. Mat. Sci. (2020)
