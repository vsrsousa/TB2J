# PAOFLOW Integration for TB2J

This directory contains an example of how to use TB2J with PAOFLOW-generated Hamiltonians.

## Overview

PAOFLOW (Projections of DFT wavefunctions on Atomic Orbital bases) can output tight-binding Hamiltonians in real space that are compatible with Wannier90's format. TB2J can read these Hamiltonians to calculate magnetic exchange interactions.

## Workflow

### 1. Generate Hamiltonian with PAOFLOW

First, run your PAOFLOW calculation and output the Hamiltonian:

```python
from PAOFLOW import PAOFLOW

# Initialize PAOFLOW with your settings
paoflow = PAOFLOW(savedir='./nscf/', outputdir='./output/', verbose=True, dft="VASP")

# Configure basis and compute Hamiltonian
basis_config = {'Fe': ['3S', '3P', '3D', '4S']}
paoflow.projections(basispath='./BASIS/', configuration=basis_config)
paoflow.pao_hamiltonian()

# Write Hamiltonian in Wannier90-compatible format
paoflow.write_Hamiltonian(fname='hamiltonian.dat')
```

**For collinear spin calculations**, PAOFLOW will create:
- `hamiltonian.dat_0` (spin up)
- `hamiltonian.dat_1` (spin down)

**For non-collinear calculations**, PAOFLOW creates:
- `hamiltonian.dat` (spinor)

### 2. Prepare Structure File

Ensure you have an atomic structure file readable by ASE (Atomic Simulation Environment):
- POSCAR (VASP format)
- CIF file
- XYZ file
- Any format supported by ASE

### 3. Run TB2J Calculation

#### For Non-Collinear Calculations:

```bash
paoflow2J.py \
    --hr_fname hamiltonian.dat \
    --atoms_fname POSCAR \
    --elements Fe Ni \
    --efermi 0.0 \
    --kmesh 7 7 7 \
    --emin -10.0 \
    --emax 0.5 \
    --output_path TB2J_results
```

#### For Collinear Calculations:

```bash
paoflow2J.py \
    --hr_fname hamiltonian.dat_0 \
    --hr_dn_fname hamiltonian.dat_1 \
    --atoms_fname POSCAR \
    --elements Fe \
    --colinear \
    --efermi 0.0 \
    --kmesh 7 7 7 \
    --emin -10.0 \
    --emax 0.5 \
    --output_path TB2J_results
```

## Command-Line Options

| Option | Description | Default |
|--------|-------------|---------|
| `--hr_fname` | Path to Hamiltonian file | hamiltonian.dat |
| `--atoms_fname` | Path to structure file | POSCAR |
| `--hr_dn_fname` | Spin-down Hamiltonian (collinear only) | None |
| `--elements` | Magnetic elements (e.g., Fe Ni) | Required |
| `--colinear` | Use collinear mode | False |
| `--efermi` | Fermi energy (eV) | 0.0 |
| `--kmesh` | k-point mesh (kx ky kz) | 5 5 5 |
| `--emin` | Min energy below E_F (eV) | -14.0 |
| `--emax` | Max energy above E_F (eV) | 0.0 |
| `--nz` | Energy integration steps | 100 |
| `--rcut` | Distance cutoff for interactions | Auto |
| `--np` | Number of parallel processes | 1 |
| `--output_path` | Output directory | TB2J_results |

## Output Files

TB2J creates several output files in the specified output directory:

- `exchange.xml` - Exchange parameters in XML format
- `exchange.txt` - Human-readable exchange parameters
- `Multibinit/` - Files for MULTIBINIT spin dynamics
- `assigned_basis.txt` - Automatically assigned basis information

## Using TB2J Results

After the calculation, you can:

1. **Plot magnon band structure**:
```bash
cd TB2J_results/Multibinit
TB2J_magnon.py --figfname magnon.png
```

2. **Analyze exchange interactions**:
Read the `exchange.txt` or `exchange.xml` files to examine:
- Isotropic exchange (J)
- Dzyaloshinskii-Moriya interaction (DMI)
- Anisotropic exchange

3. **Use in spin dynamics**:
The output is compatible with MULTIBINIT for Monte Carlo and spin dynamics simulations.

## Tips

1. **Fermi Energy**: PAOFLOW typically shifts the Fermi level to zero. Use `--efermi 0.0` unless you have a specific reason to change it.

2. **k-mesh**: Use a denser k-mesh (e.g., 7×7×7 or higher) for more accurate results. The mesh should be consistent with your PAOFLOW calculation.

3. **Energy Window**: Adjust `--emin` and `--emax` to include the relevant electronic states around the Fermi level. The default window (-14 to 0 eV) is reasonable for many systems.

4. **Magnetic Elements**: Only specify the elements that carry magnetic moments. Other elements will be treated as non-magnetic.

5. **Orbital Positions**: If needed, you can provide orbital positions via `--positions_fname` (3 columns: x, y, z in Cartesian coordinates).

## Python API

You can also use the Python API directly:

```python
from TB2J.manager import gen_exchange_paoflow

gen_exchange_paoflow(
    hr_fname='hamiltonian.dat_0',
    hr_dn_fname='hamiltonian.dat_1',
    atoms_fname='POSCAR',
    colinear=True,
    magnetic_elements=['Fe'],
    kmesh=[7, 7, 7],
    efermi=0.0,
    emin=-10.0,
    emax=0.5,
    output_path='TB2J_results'
)
```

## References

For PAOFLOW:
- Cerasoli et al., Comp. Mat. Sci. 200, 110828 (2021)
- Buongiorno Nardelli et al., Comp. Mat. Sci. 143, 462 (2018)

For TB2J:
- He et al., Comp. Phys. Commun. 264, 107938 (2021)

## Support

For issues related to:
- PAOFLOW output: https://github.com/marcobn/PAOFLOW
- TB2J calculations: https://github.com/mailhexu/TB2J
