# SACLA Light-data TR-SFX Preprocessing Pipeline

This repository runs the SACLA phytochrome **light-data** preprocessing pipeline in Python.  It follows the same practical layout as the dark-data preprocessing code, but the light-data workflow keeps the light-specific fields such as `delay`, `runID`, `eventID`, `OSF`, `relB`, and `DRL`.

The pipeline can be run in two modes:

* Remember to change the usrid in your file/path


1. **Full light-data preprocessing**: run `C1 -> C13` from `.params`, `.stream`, and delay files.
2. **Fast test from an existing C7 MATLAB/Python packed file**: run only `C8 -> C13`.

Use the existing C7 file and run from `C8` avoids re-parsing the large `.stream` file.

---

## 1. Repository layout

Recommended project layout, example:

```text
Data-cxfel/phytochrome2026/dataPhyto_reduction_light/
├── driver_light.py
├── environment.yml
├── Functions_python/
│   ├── C1_read_partialator_params_light.py
│   ├── ...
│   ├── C13_pack_phyto_light_lpf_drl_scl_bst.py
│   └── sparse_file_management.py
├── inputs/
│   ├── config_light_from_c7.yaml
│   ├── config_light_upto1ps.yaml
│   ├── delays_original.mat
│   ├── upto1ps.params
│   ├── upto1ps.stream
│   └── DATA_upto1ps/
│       ├── upto1ps_AVG_of_All.hkl
│       └── dataPhyto_upto1ps_int_sortdelay_nS135041_nBrg62530.mat
├── jobs/
│   └── tr_sfx_light_preprocess_job01.sbatch
├── logs/
│   ├── xxx.err
│   └── xxx.out
├── outputs/
│   └── upto1ps_test01/
└── results_metric/        # optional analysis/comparison scripts
```

The important input files should be placed like this:

```text
inputs/DATA_upto1ps/dataPhyto_upto1ps_int_sortdelay_nS135041_nBrg62530.mat
inputs/delays_original.mat
inputs/upto1ps.params
inputs/upto1ps.stream
```


---

## 2. What each step does

| Step | Main purpose | Main input | Main output |
|---|---|---|---|
| C1 | Read partialator parameters | `upto1ps.params` | `runID`, `eventID`, `OSF`, `relB` |
| C2 | Read stream metadata | `upto1ps.stream` | `runID`, `eventID`, `nRefl`, `DRL` |
| C3 | Extract indexed reflections | `upto1ps.stream` | per-snapshot `h k l I` files |
| C4 | Find max Miller index range | C3 hklI files | `h_max`, `k_max`, `l_max` |
| C5 | Compute reflection redundancy | C3 hklI files | redundancy array |
| C6 | Pack light-data metadata | C1, C2, `delays_original.mat` | `allInfo` with delay |
| C7 | Pack final sparse matrices | C3, C5, C6 | sorted sparse `T`, `M` |
| C8 | Uniformize delay distribution | C7 file | reduced light-data HDF5 |
| C9 | Apply LPF correction | C8 file | LPF-corrected HDF5 |
| C10 | Generate DRL mask, optional | C9 file | `maskDRL` |
| C11 | Apply or pass DRL mask | C9/C10 | LPF+DRL HDF5 |
| C12 | Apply scaling by OSF/relB | C11 | scaled HDF5 and average `.hkl` |
| C13 | Boost by multiplicity/counts | C12 | final boosted HDF5 |

By default, C11 matches the current MATLAB test behavior:

```matlab
maskDRL = ones(nS,nBrg)
```

This means Bragg reflections beyond the diffraction-resolution limit are **not removed** unless `make_drl_mask` and `use_drl_mask` are enabled in the config.

---

## 3. Set up the conda environment

### 3.1 Create the environment from `environment.yml`

From the project root:

```bash
cd /home/uwm/usrid/Data-cxfel/phytochrome2026/dataPhyto_reduction_light
conda env create -f environment.yml
```

Then activate it:

```bash
conda activate phytochrome_conda_env
```

### 3.2 Update an existing environment

If `phytochrome_conda_env` already exists, update it with:

```bash
conda env update -n phytochrome_conda_env -f environment.yml --prune
conda activate phytochrome_conda_env
```

### 3.3 Minimal `environment.yml` example

If you need to recreate the file, a minimal version is:

```yaml
name: phytochrome_conda_env
channels:
  - defaults
dependencies:
  - python=3.14
  - scipy
  - h5py
  - numpy
  - pyxtal
  - git
  - matplotlib
  - jupyterlab
  - ipykernel
```

Check the environment:

```bash
which python
python -c "import numpy, scipy, h5py, yaml; print('environment OK')"
```

---

## 4. Prepare input files

Create the input folders (existed if you clone from GitHub):

```bash
cd /home/uwm/usrid/Data-cxfel/phytochrome2026/dataPhyto_reduction_light
mkdir -p inputs/DATA_upto1ps logs outputs
```

Put the files in the following locations:

```text
inputs/DATA_upto1ps/dataPhyto_upto1ps_int_sortdelay_nS135041_nBrg62530.mat
inputs/delays_original.mat
inputs/upto1ps.params
inputs/upto1ps.stream
```

You can check them with:

```bash
ls -lh inputs/
ls -lh inputs/DATA_upto1ps/
```

Expected examples:

```text
inputs/upto1ps.params
inputs/upto1ps.stream
inputs/delays_original.mat
inputs/DATA_upto1ps/dataPhyto_upto1ps_int_sortdelay_nS135041_nBrg62530.mat
```

---

## 5. Configure the run

There are two typical config files in `inputs/`.

### 5.1 Full C1 to C13 run

Use:

```text
inputs/config_light_upto1ps.yaml
```

Example:

```yaml
project_root: /home/uwm/usrid/Data-cxfel/phytochrome2026/dataPhyto_reduction_light
functions_dir: /home/uwm/usrid/Data-cxfel/phytochrome2026/dataPhyto_reduction_light/Functions_python

start_from: C1

data_name: upto1ps
stream: /home/uwm/usrid/Data-cxfel/phytochrome2026/dataPhyto_reduction_light/inputs/upto1ps.stream
params: /home/uwm/usrid/Data-cxfel/phytochrome2026/dataPhyto_reduction_light/inputs/upto1ps.params
delays: /home/uwm/usrid/Data-cxfel/phytochrome2026/dataPhyto_reduction_light/inputs/delays_original.mat
output: /home/uwm/usrid/Data-cxfel/phytochrome2026/dataPhyto_reduction_light/outputs/upto1ps_test01

a: 54.22
b: 115.78
c: 117.08
wavelength: 1.77

ss: 155
tmin: -84
tmax: 550
seed: 123

remove_negative_pixels: false
make_drl_mask: false
use_drl_mask: false
generate_hkl_avg: true
```

This runs:

```text
C1 -> C2 -> C3 -> C4 -> C5 -> C6 -> C7 -> C8 -> C9 -> C11 -> C12 -> C13
```

### 5.2 Fast test from existing C7 file

Use:

```text
inputs/config_light_from_c7.yaml
```

Example:

```yaml
project_root: /home/uwm/usrid/Data-cxfel/phytochrome2026/dataPhyto_reduction_light
functions_dir: /home/uwm/usrid/Data-cxfel/phytochrome2026/dataPhyto_reduction_light/Functions_python

start_from: C8

data_name: upto1ps
c7_file: /home/uwm/usrid/Data-cxfel/phytochrome2026/dataPhyto_reduction_light/inputs/DATA_upto1ps/dataPhyto_upto1ps_int_sortdelay_nS135041_nBrg62530.mat
output: /home/uwm/usrid/Data-cxfel/phytochrome2026/dataPhyto_reduction_light/outputs/upto1ps_from_C7_test01

a: 54.22
b: 115.78
c: 117.08
wavelength: 1.77

ss: 155
tmin: -84
tmax: 550
seed: 123

remove_negative_pixels: false
make_drl_mask: false
use_drl_mask: false
generate_hkl_avg: true
```

This starts from:

```text
C8 -> C9 -> C11 -> C12 -> C13
```


---

## 6. Run with the driver

The driver reads the YAML config, calls the C1-C13 functions, and prints step-by-step progress and memory information.

### 6.1 Fast local/interactice test from C7

```bash
cd /home/uwm/usrid/Data-cxfel/phytochrome2026/dataPhyto_reduction_light
conda activate phytochrome_conda_env
python ./driver_light.py --config ./inputs/config_light_from_c7.yaml
```

This is the recommended first test.

### 6.2 Full local/interactice run

```bash
cd /home/uwm/usrid/Data-cxfel/phytochrome2026/dataPhyto_reduction_light
conda activate phytochrome_conda_env
python ./driver_light.py --config ./inputs/config_light_upto1ps.yaml
```

The full run reads the large `.stream` file, so it should usually be submitted as a SLURM job.

---

## 7. Run with SLURM sbatch

The job file is:

```text
jobs/tr_sfx_light_preprocess_job01.sbatch
```

Submit it with:

```bash
cd /home/uwm/usrid/Data-cxfel/phytochrome2026/dataPhyto_reduction_light
sbatch jobs/tr_sfx_light_preprocess_job01.sbatch
```

Check job status:

```bash
squeue -u $USER
```

Check logs:

```bash
ls -lh logs/
tail -f logs/tr_sfx_light_job01-<JOBID>.out
```

After the job finishes, check SLURM resource usage:

```bash
sacct -j <JOBID> --format=JobID,State,Elapsed,ReqMem,MaxRSS,ExitCode
```

### 7.1 Switching between full run and C7 test in sbatch

Inside `jobs/tr_sfx_light_preprocess_job01.sbatch`, the full run line is:

```bash
python ./driver_light.py --config ./inputs/config_light_upto1ps.yaml
```

For a faster first test from the existing C7 file, comment that line and use:

```bash
python ./driver_light.py --config ./inputs/config_light_from_c7.yaml
```

### 7.2 Typical sbatch settings

For large light-data processing, start with something like:

```bash
#SBATCH --time=12:00:00
#SBATCH --mem=120G
#SBATCH --ntasks=2
#SBATCH --cpus-per-task=2
#SBATCH --partition=bioxfel
```

If C7 or C8 runs out of memory, increase `--mem`. If the stream parsing or sparse matrix packing takes too long, increase `--time`.

---

## 8. Output files

Outputs are written to the folder specified in the config, for example:

```text
outputs/upto1ps_from_C7_test01/
outputs/upto1ps_test01/
```

Important outputs include:

```text
dataPhyto_upto1ps_C8_unifdelay.hdf5
dataPhyto_upto1ps_C9_LPF.hdf5
dataPhyto_upto1ps_C11_LPF_DRL.hdf5
dataPhyto_upto1ps_C12_LPF_DRL_SCL.hdf5
dataPhyto_upto1ps_C13_LPF_DRL_SCL_BST.hdf5
```

If `generate_hkl_avg: true`, C12 also writes an average `.hkl` file for result-metric analysis.

The HDF5 files store vectors such as:

```text
/delay
/OSF
/relB
/DRL
/runID
/eventID
/miller_h
/miller_k
/miller_l
```

Sparse matrices are stored in SciPy CSC-style HDF5 groups:

```text
/T/data
/T/indices
/T/indptr
/T/shape

/M/data
/M/indices
/M/indptr
/M/shape
```

---

## 9. Understanding driver memory logs

The driver prints logs like:

```text
[START] C8 generate_uniform_delay | PID=12345 | current=0.50 GB | process_peak=0.50 GB
[MONITOR] C8 generate_uniform_delay | elapsed=300s | current=72.80 GB | step_peak=78.40 GB | process_peak=78.40 GB
[ END ] C8 generate_uniform_delay | current=21.30 GB | step_peak=101.20 GB | process_peak=101.20 GB
```

Meanings:

| Field | Meaning |
|---|---|
| `current` | Current RSS memory of the Python process |
| `step_peak` | Highest sampled RSS during the current step |
| `process_peak` | Highest RSS since the driver started |

Use this to identify whether C7, C8, or another step is the memory bottleneck.

---

## 10. Result-metric analysis and MATLAB/Python comparison

The safest way to validate the Python pipeline is to compare it step by step against the MATLAB output.

### 10.1 Quick HDF5 shape check

Run:

```bash
python - <<'PY'
import h5py
from pathlib import Path

files = [
    "outputs/upto1ps_from_C7_test01/dataPhyto_upto1ps_C8_unifdelay.hdf5",
    "outputs/upto1ps_from_C7_test01/dataPhyto_upto1ps_C9_LPF.hdf5",
    "outputs/upto1ps_from_C7_test01/dataPhyto_upto1ps_C11_LPF_DRL.hdf5",
    "outputs/upto1ps_from_C7_test01/dataPhyto_upto1ps_C12_LPF_DRL_SCL.hdf5",
    "outputs/upto1ps_from_C7_test01/dataPhyto_upto1ps_C13_LPF_DRL_SCL_BST.hdf5",
]

for name in files:
    p = Path(name)
    if not p.exists():
        print(f"missing: {p}")
        continue
    with h5py.File(p, "r") as f:
        print("\n", p)
        for key in ["delay", "OSF", "relB", "DRL", "runID", "eventID", "miller_h"]:
            if key in f:
                print(f"  {key}: {f[key].shape}")
        for key in ["T", "M"]:
            if key in f:
                shape = tuple(f[key]["shape"][:])
                nnz = len(f[key]["data"])
                print(f"  {key}: shape={shape}, nnz={nnz}")
PY
```

For the C8 uniform-delay output, check:

```text
nS after C8 should match the MATLAB uniform-delay result, if the same seed and sampling rule are used.
each delay bin should have at most ss snapshots.
```

### 10.2 Check delay uniformization

```bash
python - <<'PY'
import h5py
import numpy as np

path = "outputs/upto1ps_from_C7_test01/dataPhyto_upto1ps_C8_unifdelay.hdf5"
ss = 155

with h5py.File(path, "r") as f:
    delay = f["delay"][:]

u, c = np.unique(delay, return_counts=True)
print("nS =", delay.size)
print("unique delays =", u.size)
print("max count per delay =", c.max())
print("all bins <= ss:", bool(np.all(c <= ss)))
print("delay range =", delay.min(), delay.max())
PY
```

### 10.3 Compare average `.hkl` files with `results_metric`

If your repository has a `results_metric/` folder, use it to compare the MATLAB and Python `.hkl` average outputs.

Example expected files:

```text
MATLAB output:
  path/to/matlab/dataPhyto_upto1ps_C12_LPF_DRL_SCL_AVG.hkl

Python output:
  outputs/upto1ps_from_C7_test01/dataPhyto_upto1ps_C12_LPF_DRL_SCL_AVG.hkl
```

Run shell/q-bin intensity analysis, for example:

```bash
python results_metric/shell_intensity.py \
  outputs/upto1ps_from_C7_test01/dataPhyto_upto1ps_C12_LPF_DRL_SCL_AVG.hkl
```

Run the same command for the MATLAB `.hkl` file, then compare the generated `.txt` and `.png` files.

The useful comparison points are:

```text
1. Same number of common reflections.
2. Similar q-shell average intensity trend.
3. Similar per-shell counts.
4. High correlation between MATLAB and Python intensities for common h/k/l reflections.
5. No unexpected sign flip or q-shell-dependent scaling error.
```

### 10.4 Direct common-reflection correlation check

If you do not yet have a dedicated comparison script, this quick check compares two `.hkl` files by common `h k l`:

```bash
python - <<'PY'
from pathlib import Path
import numpy as np
import pandas as pd

matlab_hkl = Path("path/to/matlab/dataPhyto_upto1ps_C12_LPF_DRL_SCL_AVG.hkl")
python_hkl = Path("outputs/upto1ps_from_C7_test01/dataPhyto_upto1ps_C12_LPF_DRL_SCL_AVG.hkl")

# Adjust column names if your .hkl file has a different number/order of columns.
cols = ["h", "k", "l", "I", "sigma"]

def read_hkl(path):
    df = pd.read_csv(path, delim_whitespace=True, comment="#", header=None)
    df = df.iloc[:, :len(cols)]
    df.columns = cols
    return df

m = read_hkl(matlab_hkl)
p = read_hkl(python_hkl)

merged = m.merge(p, on=["h", "k", "l"], suffixes=("_matlab", "_python"))
print("MATLAB reflections:", len(m))
print("Python reflections:", len(p))
print("Common reflections:", len(merged))

x = merged["I_matlab"].to_numpy(float)
y = merged["I_python"].to_numpy(float)
mask = np.isfinite(x) & np.isfinite(y)
x = x[mask]
y = y[mask]

corr = np.corrcoef(x, y)[0, 1]
mae = np.mean(np.abs(x - y))
rel = np.mean(np.abs(x - y) / np.maximum(np.abs(x), 1e-12))

print("Pearson correlation:", corr)
print("MAE:", mae)
print("mean relative error:", rel)
PY
```

If the correlation is unexpectedly low, debug in this order:

```text
C1: compare runID/eventID/OSF/relB counts.
C2: compare runID/eventID/nRefl/DRL counts.
C6: compare allInfo row count and delay matching.
C7: compare nS, nBrg, miller_h/k/l order, and sparse matrix nnz.
C8: compare nS after uniform delay and delay-bin counts.
C9-C13: compare a small number of common sparse entries or final .hkl intensity values.
```

---

## 11. Common problems

### Problem: `ModuleNotFoundError`

Check that `functions_dir` in the config points to the correct folder:

```yaml
functions_dir: /home/uwm/usrid/Data-cxfel/phytochrome2026/dataPhyto_reduction_light/Functions_python
```

### Problem: input file not found

Check absolute paths in the config:

```bash
ls -lh inputs/upto1ps.params
ls -lh inputs/upto1ps.stream
ls -lh inputs/delays_original.mat
ls -lh inputs/DATA_upto1ps/dataPhyto_upto1ps_int_sortdelay_nS135041_nBrg62530.mat
```

### Problem: job killed by memory limit

Check the end of the `.err` file and `sacct`:

```bash
tail -n 50 logs/*.err
sacct -j <JOBID> --format=JobID,State,Elapsed,ReqMem,MaxRSS,ExitCode
```

If C8 reaches close to the requested memory, increase:

```bash
#SBATCH --mem=160G
```

### Problem: C11 result differs from MATLAB or other data resourses

Make sure the DRL mask settings match the MATLAB behavior. For current MATLAB-style pass-through behavior:

```yaml
make_drl_mask: false
use_drl_mask: false
```

For real DRL masking:

```yaml
make_drl_mask: true
use_drl_mask: true
```

---


