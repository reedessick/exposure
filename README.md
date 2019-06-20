# exposure

This is a simple library that's meant to compute the exposure of a network of GW detectors. In particular, this is decomposed into 

  - estimating the PSD within (join) lock segments as a function of time via massive parallelization
    - `compute_psd` (parallelized via `condor-compute_psd`)

  - estimation of horizon for a particular detector given a particular PSD
    - `compute_horizon` (parallelized via `condor-compute_horizon`)

  - estimation of network sensitivies given a set of detectors and their associated horizons
    - `compute_network_sensitivity` (parallelized via `condor-compute_network_sensitivity`)

  - estimation of overall exposure by integrating network sensitivities over time
    - `compute_exposure`

The work flow will likely run in this order. We note that there is further parallelization available by writing a single DAG that combines `compute_psd`, `compute_horizon`, and `compute_network_sensitivity` together with appropriate parent/child relationships. We could add `compute_exposure` to that as well, but it will require all jobs to finish prior and therefore might as well just be a postscript or run by hand afterwards.

We also provide some basic functionality to estimate diurnal cycles.

  - compute a stacked histogram (with some user-specified periodicity)
    - `diurnal_hist`

  - compute the dft of lock segments and run a peak finding algorithm to identify periodic components
    - `diuranl_dft`

  - compute a time-frequency representation of lock segments to look for variabilitoy of periodic elements
    - `diurnal_timefreq`

## Dependencies

  - [skymap_statistics](https://github.com/reedessick/skymap_statistics)
  - numpy
  - matplotlib
  - healpy
  - lalsuite (lalframe.frread, gpstime.tconvert)
  - corner.py

**NOTE**, the library is not Python3 compatible (it uses deprecated functionality only available in Python2)
