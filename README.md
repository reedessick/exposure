# exposure

This is a simple library that's meant to compute the exposure of a network of GW detectors. In particular, this is decomposed into 

  - estimating the PSD within (join) lock segments as a function of time via massive parallelization
  - estimation of range for a particular detector given a particular PSD
  - estimation of network sensitivies given a set of detectors and their associated ranges
  - estimation of overall exposure by integrating network sensitivities over time


