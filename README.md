# exposure

This is a simple library that's meant to compute the exposure of a network of GW detectors. In particular, this is decomposed into 

  - estimating the PSD within (join) lock segments as a function of time via massive parallelization
    - compute_psd (parallelized via condor-compute_psd)

  - estimation of range for a particular detector given a particular PSD
    - compute_range (parallelized via condor-compute_range)

  - estimation of network sensitivies given a set of detectors and their associated ranges
    - compute_network_sensitivity (parallelized via condor-compute_network_sensitivity)

  - estimation of overall exposure by integrating network sensitivities over time
    - `compute_exposure`

The work flow will likely run in this order. We note that there is further parallelization available by writing a single DAG that combines compute_psd, compute_range, and compute_network_sensitivity together with appropriate parent/child relationships. We could add compute_exposure to that as well, but it will require all jobs to finish prior and therefore might as well just be a postscript or run by hand afterwards.
