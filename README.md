# SBND TPC calibration macros

## Run gain calibration loop

- 1. Producing MPV dE/dx vs dQ/dx 2D distributions using calibration ntuples

Modified box model\

```
$ root -l -b -q "src/run_recom_loop_mb.C(<run_number>)"
```

Ellipsoidal box model\

```
$ root -l -b -q "src/run_recom_loop_emb.C(<run_number>)"
```

For data, replace <run_number> with run number. A sample list file should present in data directory.
For MC, use "0" for <run_number>.