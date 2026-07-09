# Light-data driver template

This folder follows the same idea as the dark-data preprocessing folder:

```text
project/
  Functions_python/      # C1-C13 light-data functions
  inputs/                # config files and raw inputs
  jobs/                  # SLURM sbatch files
  logs/                  # SLURM logs
  outputs/               # pipeline outputs
  driver_light.py        # YAML-driven driver with memory monitor
```

## 1. Fast first test from existing C7 file

If you already have

```text
inputs/DATA_upto1ps/dataPhyto_upto1ps_int_sortdelay_nS135041_nBrg62530.mat
```

run only C8-C13 first:

```bash
cd /home/uwm/gengh/Data-cxfel/phytochrome2026/dataPhyto_reduction_light
python ./driver_light.py --config ./inputs/config_light_from_c7.yaml
```

This is the recommended first check because it avoids re-parsing the 4.5 GB stream file.

## 2. Full C1-C13 run

After the C8-C13 test works:

```bash
cd /home/uwm/gengh/Data-cxfel/phytochrome2026/dataPhyto_reduction_light
python ./driver_light.py --config ./inputs/config_light_upto1ps.yaml
```

or submit to SLURM:

```bash
sbatch jobs/tr_sfx_light_preprocess_job01.sbatch
```

## 3. Important config keys

```yaml
start_from: C1       # C1 for full pipeline, C8 for existing C7 file
stream: ...          # needed only for start_from: C1
params: ...          # needed only for start_from: C1
delays: ...          # needed only for start_from: C1
c7_file: ...         # needed only for start_from: C8
output: ...          # output folder

a: 54.22
b: 115.78
c: 117.08
wavelength: 1.77

ss: 155
tmin: -84
tmax: 550
seed: 123

make_drl_mask: false
use_drl_mask: false
```

By default C11 matches the MATLAB test behavior:

```matlab
maskDRL = ones(nS,nBrg)
```

so Bragg reflections beyond DRL are **not** removed unless both `make_drl_mask` and `use_drl_mask` are enabled.

## 4. Notes

- `C2_read_stream_files_light.py` and `C3_grab_hklI_sacla2021_light.py` read the stream file chunk by chunk instead of loading the whole 4.5 GB stream into memory.
- `driver_light.py` prints `[START]`, `[MONITOR]`, `[END]`, and `[FAIL]` messages similar to the dark-data driver.
- C8 can read either native Python HDF5 sparse matrices or MATLAB v7.3 sparse groups.
- For debugging, use `config_light_from_c7.yaml` first.
