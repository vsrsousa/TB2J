Use TB2J with PAOflow
=====================

This tutorial explains how to use TB2J with PAOflow to calculate magnetic exchange parameters.

What is PAOflow?
----------------

PAOflow (Projections of Atomic Orbitals flow) is an alternative to Wannier90 for constructing 
tight-binding Hamiltonians from DFT calculations. It uses projections of atomic orbitals (PAO) 
to generate localized basis functions.

Key features of PAOflow:

* Generates tight-binding Hamiltonians from DFT calculations (QE, VASP, etc.)
* Uses projections of atomic orbitals instead of maximally localized Wannier functions
* Outputs Hamiltonians in **Wannier90-compatible format** (_hr.dat files)
* Supports both collinear and non-collinear spin calculations
* Alternative workflow to Wannier90 with different localization approach

PAOflow Output Format
---------------------

PAOflow outputs Hamiltonians in Wannier90-compatible ``_hr.dat`` format using the 
``write_Hamiltonian()`` method. This means:

* **You can use the standard** ``wann2J.py`` **command with PAOflow output!**
* The ``paoflow2J.py`` command is a convenience wrapper with PAOflow-specific defaults

For collinear calculations, PAOflow generates:

* ``hamiltonian_0.dat`` - Hamiltonian for spin-up channel (Wannier90 _hr.dat format)
* ``hamiltonian_1.dat`` - Hamiltonian for spin-down channel (Wannier90 _hr.dat format)

For non-collinear/SOC calculations:

* ``hamiltonian.dat`` - Hamiltonian with spin-orbit coupling (Wannier90 _hr.dat format)

The _hr.dat format contains:

* Line 1: Header comment
* Line 2: Number of Wannier functions
* Line 3: Number of R-lattice vectors
* Lines 4+: Degeneracy weights (all 1's in PAOflow output)
* Data lines: ``Rx Ry Rz orb1 orb2 H_real H_imag`` for each matrix element

Running TB2J with PAOflow
--------------------------

Method 1: Using paoflow2J.py (Recommended)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The ``paoflow2J.py`` command provides PAOflow-specific defaults:

.. code-block:: bash

   paoflow2J.py --path /path/to/paoflow/output \
                --efermi 6.15 \
                --kmesh 5 5 5 \
                --elements Mn Fe \
                --emin -10.0 \
                --emax 0.0

This automatically looks for ``hamiltonian_0.dat`` and ``hamiltonian_1.dat`` for 
collinear calculations.

Method 2: Using wann2J.py (Also Works)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Since PAOflow outputs Wannier90-compatible files, you can use:

.. code-block:: bash

   wann2J.py --path /path/to/paoflow/output \
             --prefix_up hamiltonian_0 \
             --prefix_down hamiltonian_1 \
             --efermi 6.15 \
             --kmesh 5 5 5 \
             --elements Mn Fe \
             --emin -10.0 \
             --emax 0.0

Parameters:

* ``--path`` - Directory containing PAOflow output files
* ``--prefix_up`` - Prefix for spin-up _hr.dat file (default for paoflow2J.py: ``hamiltonian_0``)
* ``--prefix_down`` - Prefix for spin-down _hr.dat file (default for paoflow2J.py: ``hamiltonian_1``)
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
                --efermi 6.15 \
                --kmesh 5 5 5 \
                --elements Mn Fe \
                --emin -10.0 \
                --emax 0.0

The ``--spinor`` flag indicates a non-collinear calculation. This automatically uses 
``--prefix_spinor hamiltonian`` (the default).

Structure File
~~~~~~~~~~~~~~

You need to provide a structure file containing atomic positions. PAOflow's _hr.dat 
files don't include atomic structure information:

.. code-block:: bash

   paoflow2J.py --posfile POSCAR --path /path/to/paoflow/output ...

The ``--posfile`` parameter accepts any format supported by ASE (POSCAR, CIF, XYZ, etc.). 
You can typically use the same structure file that was input to your DFT calculation.

Complete Example
----------------

Here's a complete workflow for calculating exchange parameters in a magnetic material:

Step 1: Run PAOflow
~~~~~~~~~~~~~~~~~~~

First, generate the tight-binding Hamiltonian using PAOflow:

.. code-block:: python

   from PAOFLOW import PAOFLOW
   
   # Initialize PAOflow with your DFT output
   paoflow = PAOFLOW.PAOFLOW(savedir='yourcode.save', outputdir='./output')
   
   # Read atomic projections from DFT output
   paoflow.read_atomic_proj_QE()  # Or read_atomic_proj_VASP() for VASP
   
   # Set projectability threshold
   paoflow.projectability(pthr=0.95)
   
   # Build PAO Hamiltonian in real space
   paoflow.pao_hamiltonian()
   
   # Write Hamiltonian in Wannier90 format
   # For collinear: writes hamiltonian_0.dat and hamiltonian_1.dat
   # For non-collinear: writes hamiltonian.dat
   paoflow.write_Hamiltonian()
   
   paoflow.finish_execution()

Step 2: Run TB2J
~~~~~~~~~~~~~~~~

Calculate the exchange parameters:

.. code-block:: bash

   paoflow2J.py --path ./output \
                --posfile POSCAR \
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
   well-localized orbitals. Check the projectability in PAOflow output.

3. **Energy Window**: Choose ``--emin`` and ``--emax`` to include all relevant orbitals 
   (typically d or f orbitals for magnetic atoms and p orbitals for ligands).

4. **k-mesh Convergence**: Test convergence with respect to the k-point mesh using 
   ``--kmesh``. A denser mesh may be needed for complex materials.

5. **File Names**: By default, paoflow2J.py looks for ``hamiltonian_0.dat`` and 
   ``hamiltonian_1.dat`` (collinear) or ``hamiltonian.dat`` (non-collinear). If PAOflow 
   writes different names, use ``--prefix_up``, ``--prefix_down``, or ``--prefix_spinor``.

6. **Structure File**: Always provide a structure file with ``--posfile`` since PAOflow's 
   _hr.dat files don't contain atomic positions.

Comparison with Wannier90
--------------------------

PAOflow vs Wannier90:

+------------------------+--------------------------------+---------------------------+
| Feature                | Wannier90                      | PAOflow                   |
+========================+================================+===========================+
| Output format          | Text files (_hr.dat)           | Same (_hr.dat)            |
+------------------------+--------------------------------+---------------------------+
| Orbital construction   | Maximally localized WFs        | Projected atomic orbitals |
+------------------------+--------------------------------+---------------------------+
| Interface in TB2J      | ``wann2J.py``                  | ``paoflow2J.py`` or       |
|                        |                                | ``wann2J.py``             |
+------------------------+--------------------------------+---------------------------+
| Structure info         | Included in .win file          | Separate file needed      |
+------------------------+--------------------------------+---------------------------+

**Key Point**: PAOflow outputs are **Wannier90-compatible**, so both ``wann2J.py`` and 
``paoflow2J.py`` work with PAOflow output. The main difference is the orbital construction 
method (MLWFs vs. PAOs), which may affect localization quality.

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

* PAOflow GitHub: https://github.com/marcobn/PAOflow (original)
* PAOflow fork: https://github.com/vsrsousa/PAOFLOW (with TB2J examples)
* TB2J documentation: https://tb2j.readthedocs.io
* Wannier90 tutorial: :doc:`wannier`
* Exchange parameter reference: :doc:`parameters`

References
----------

If you use PAOflow with TB2J, please cite:

* **TB2J**: He Xu et al., "TB2J: A Package for Computing Magnetic Interaction Parameters", 
  Comp. Phys. Commun. (2020)

* **PAOflow**: F.T. Cerasoli et al., "Advanced modeling of materials with PAOFLOW 2.0: 
  New features and software design", Comp. Mat. Sci. 200, 110828 (2021)

* **PAOflow original**: M. Buongiorno Nardelli et al., "PAOFLOW: A utility to construct 
  and operate on ab initio Hamiltonians from the Projections of electronic wavefunctions 
  on Atomic Orbital bases", Comp. Mat. Sci. 143, 462 (2018)
