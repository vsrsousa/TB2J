# PAOFLOW Example for TB2J

This directory contains an example of using TB2J with PAOFLOW-generated Hamiltonians.

## Files

- `run_paoflow2J.py` - Python script showing how to use the Python API
- `README.md` - This file

## Quick Start

### Step 1: Generate PAOFLOW Hamiltonian

First, run PAOFLOW to generate the Hamiltonian files. Here's a minimal example:

```python
from PAOFLOW import PAOFLOW

# Initialize PAOFLOW
paoflow = PAOFLOW(
    savedir='./nscf/',
    outputdir='./output/',
    verbose=True,
    dft="VASP"  # or "QE" for Quantum ESPRESSO
)

# Configure atomic orbital basis
basis_config = {'Fe': ['3S', '3P', '3D', '4S']}
paoflow.projections(basispath='../../BASIS/', configuration=basis_config)

# Build PAO Hamiltonian
paoflow.pao_hamiltonian()

# Write Hamiltonian in Wannier90-compatible format
paoflow.write_Hamiltonian(fname='hamiltonian.dat')
```

For **collinear calculations** (nspin=2), this will create:
- `hamiltonian.dat_0` (spin up)
- `hamiltonian.dat_1` (spin down)

For **non-collinear calculations** (nspin=1 with SOC or nspin=4), this creates:
- `hamiltonian.dat` (spinor)

### Step 2: Prepare Structure File

Make sure you have a structure file (e.g., POSCAR) in the same directory.

### Step 3: Run TB2J

#### Using Python API:

```bash
python run_paoflow2J.py
```

#### Using Command Line:

For collinear calculations (default):
```bash
paoflow2J.py \
    --hr_up hamiltonian.dat_0 \
    --hr_dn hamiltonian.dat_1 \
    --poscar POSCAR \
    --elements Fe \
    --efermi 0.0 \
    --kmesh 7 7 7 \
    --output_path TB2J_results
```

For non-collinear calculations:
```bash
paoflow2J.py \
    --hr_up hamiltonian.dat \
    --poscar POSCAR \
    --elements Fe \
    --non_colinear \
    --efermi 0.0 \
    --kmesh 7 7 7 \
    --output_path TB2J_results
```

## Expected Output

After successful execution, you'll find these files in `TB2J_results/`:

- `exchange.txt` - Exchange parameters in text format
- `exchange.xml` - Exchange parameters in XML format
- `Multibinit/` - Input files for MULTIBINIT
- `assigned_basis.txt` - Orbital assignments

## Visualizing Results

To plot the magnon dispersion:

```bash
cd TB2J_results/Multibinit
TB2J_magnon.py --figfname magnon.png
```

## Notes

1. **Fermi Energy**: PAOFLOW typically shifts the Fermi level to zero in the output, so use `--efermi 0.0`

2. **k-mesh**: Should be consistent with or denser than your PAOFLOW calculation

3. **Energy Window**: Adjust `--emin` and `--emax` to include relevant states near the Fermi level

4. **Magnetic Elements**: Only specify elements that carry magnetic moments

## Troubleshooting

**Issue**: "hamiltonian.dat_0 not found"
- **Solution**: Make sure PAOFLOW's `write_Hamiltonian()` was called and completed successfully

**Issue**: "For collinear calculations, hr_dn must be provided"
- **Solution**: For nspin=2 calculations, both spin channels are required. Use `--hr_up hamiltonian.dat_0 --hr_dn hamiltonian.dat_1`

**Issue**: Orbital positions mismatch
- **Solution**: Provide orbital positions explicitly using `--positions_fname` if auto-assignment fails

## For More Information

See the main documentation: `docs/paoflow_integration.md`
