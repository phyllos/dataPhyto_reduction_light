# Core light-data Python rewrite

This folder is a minimal SACLA light-data Python translation in the same spirit as the dark-data Python `Functions` folder and the original light-data MATLAB `Functions` folder.

It intentionally keeps only core functions/scripts and avoids the larger full-pipeline driver.

## Main difference from dark data

Dark data can be packed mainly by snapshot/serial number. Light data must keep these extra columns aligned:

```text
runID, eventID/tag, nBrg, DRL, delay, OSF, relB
```

Therefore C1/C2/C3/C6/C7 are light-specific. C4/C5 are close to the dark-data logic.

## Files

```text
C1_read_partialator_params_light.py       params -> nRun,nEvent,OSF,relB
C2_read_stream_files_light.py             stream -> nRun,nEvent,nRefl,DRL
C3_grab_hklI_sacla2021_light.py           stream -> per-snapshot hklI DAT files
C4_find_hkl_max_light.py                  hklI files -> hmax,kmax,lmax
C5_find_hkl_redundancy_light.py           hklI files -> redundancy
C6_pack_phyto_allInfo_light.py            C1+C2+delays -> allInfo
C7_pack_phyto_final_light.py              allInfo+hklI+redundancy -> sparse T/M sorted by delay
C8_generate_uniform_delay_sacla2021_hdf5.py  uniform delay subset
C9_pack_phyto_light_lpf.py                LPF correction
C10_generate_maskDRL_light.py             optional real DRL mask
C11_pack_phyto_light_lpf_drl.py           apply DRL mask, or default mask=ones
C12_pack_phyto_light_lpf_drl_scl.py       scaling by OSF/relB
C13_pack_phyto_light_lpf_drl_scl_bst.py   boosting by reflection multiplicity
sparse_file_management.py                 small HDF5/sparse read-write helper
```

## Minimal C8-C13 run from your existing MATLAB C7 file

You can start from the uploaded MATLAB v7.3 file:

```bash
python C8_generate_uniform_delay_sacla2021_hdf5.py \
  dataPhyto_upto1ps_int_sortdelay_nS135041_nBrg62530.mat \
  --ss 155 --tmin -84 --tmax 550 --seed 123

python C9_pack_phyto_light_lpf.py \
  dataPhyto_upto1ps_int_sortdelay_unifdelay_nS68084_nBrg62530.hdf5

python C11_pack_phyto_light_lpf_drl.py \
  dataPhyto_upto1ps_int_sortdelay_unifdelay_nS68084_nBrg62530_LPF.hdf5

python C12_pack_phyto_light_lpf_drl_scl.py \
  dataPhyto_upto1ps_int_sortdelay_unifdelay_nS68084_nBrg62530_LPF_DRL.hdf5

python C13_pack_phyto_light_lpf_drl_scl_bst.py \
  dataPhyto_upto1ps_int_sortdelay_unifdelay_nS68084_nBrg62530_LPF_DRL_SCL.hdf5
```

C11 default matches your current MATLAB behavior:

```matlab
maskDRL = ones(nS,nBrg)
```

To apply a real DRL mask, first run C10 and then pass `--mask-file` to C11.

## Notes

- Outputs are native HDF5 `.hdf5`, not MATLAB `.mat` files.
- Sparse matrices are saved in scipy CSC format: `data`, `indices`, `indptr`, `shape`.
- C8 can read both MATLAB v7.3 sparse groups and the native HDF5 sparse format used here.
- `delays_original.mat` for `upto1ps` uses `delay_1ps` and `tag_1ps`; C6 handles this automatically.
