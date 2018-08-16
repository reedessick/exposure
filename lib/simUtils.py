__doc__ = "a module housing useful routines for simulation"
__author__ = "reed.essick@ligo.org"

#-------------------------------------------------

import numpy as np
import healpy as hp

### non-standard libraries
import triangulate

#-------------------------------------------------

DEFAULT_SIZE = 100
DEFAULT_NSIDE = 64

DEFAULT_SNR_THRESHOLD = 8.
DEFAULT_FLOW = 32. ### Hz
DEFAULT_FHIGH = 1024. ### Hz

DEFAULT_TIME_DELAY_ERROR = 1.e-3 ### sec

#-------------------------------------------------

def draw_times(start, stop, size=DEFAULT_SIZE):
    """
    draw size random times within [start, stop)
    """
    ans = start + (stop-start)*np.random.rand(size)
    return ans[ans.argsort()]

#-------------------------------------------------

def horizon2geo_antenna(detector, horizon, theta, phi, psi):
    return horizons2geo_antenna([(detector, horizon)], theta, phi, psi)[0]

def horizons2geo_antenna(horizons, theta, phi, psi):
    """
    compute antenna patterns in geographic coordinates cooresponding to detectors
    horizons = [(detector, horizon), (detector, horizon), ...]
    """
    ### set up locations
    assert (len(theta)==len(phi)) and (len(theta)==len(psi)), 'all position vectors must have the same dimension'

    ### sum over detectors
    network = np.zeros_like(theta, dtype='float')
    for detector, horizon in horizons:
        fp, fx = detector.antenna_patterns(theta, phi, psi)
        network += (np.abs(fp)**2 + np.abs(fx)**2)*horizon**2

    ### record the sensitivity linear in ~F/ASD ~ range
    network = network**0.5

    return network

def detectors2horizons(detectors, snr_threshold=DEFAULT_SNR_THRESHOLD, flow=DEFAULT_FLOW, fhigh=DEFAULT_FHIGH):
    """
    iterates over detectors and delegates to detector2horizon
    """
    return [(detector, detector2horizon(detector, snr_threshold=snr_threshold, flow=flow, fhigh=fhigh)) for detector in detectors]
        
def detector2horizon(detector, snr_threshold=DEFAULT_SNR_THRESHOLD, flow=DEFAULT_FLOW, fhigh=DEFAULT_FHIGH):
    """
    a janky way of estimating the detector horizon.
    THIS COULD BE GREATLY IMPROVED
    """
    freq = detector.psd.get_freqs()
    psd = detector.psd.get_psd()
    truth = (flow<=freq)*(freq<=fhigh)

    h_of_f = 1e-23 * (freq[truth]/100.)**(7./6) ### FIXME use lalsim or something to get units correct...

    ### compute the optimal SNR for this strain
    snr = (4.*np.trapz(h_of_f**2./psd[truth], x=freq[truth]))**0.5 ### approximate the integral

    ### compute the horizon from this
    horizon = snr/snr_threshold ### this is the scaling, again making some broad approximations

    return horizon

#-------------------------------------------------

def draw_positions(theta, phi, network, size=DEFAULT_SIZE):
    order = network.argsort()
    args = np.random.randint(0, high=len(theta), size=size)
    return zip(theta[args], phi[args])

def geo_position2geo_skymap(
        (source_theta, source_phi),
        (map_theta, map_phi),
        detectors,
        network,
        time_delay_error=DEFAULT_TIME_DELAY_ERROR,
    ):
    """
    compute modulated triangulation rings given the detectors and network (antenna pattern)
    """
    logskymap = np.zeros_like(map_theta, dtype=float)

    ### add in triangulation rings
    source_n = hp.ang2vec(source_theta, source_phi)
    map_n = hp.ang2vec(map_theta, map_phi)

    prefact = 0.5/time_delay_error**2
    for ind, d1 in enumerate(detectors):
        r1 = d1.dr
        for d2 in detectors[ind+1:]:
            dr = r1-d2.dr
            source_dt = np.sum(source_n*dr)
            map_dt = np.sum(map_n*dr, axis=1)

            logskymap -= prefact*(source_dt-map_dt)**2 ### assume Gaussian in time-delay errors

    logskymap += np.log(network) ### multiply by antenna patterns to modulate triangulation rings

    ### return a normalized skymap
    logskymap -= np.max(logskymap)
    return np.exp(logskymap)/np.sum(np.exp(logskymap))

#-------------------------------------------------

def simulate_geo_skymaps(
        start,
        stop,
        detectors,
        size=DEFAULT_SIZE,
        nside=DEFAULT_NSIDE,
        snr_threshold=DEFAULT_SNR_THRESHOLD,
        flow=DEFAULT_FLOW,
        fhigh=DEFAULT_FHIGH,
        time_delay_error=DEFAULT_TIME_DELAY_ERROR,
    ):
    """
    the whole kit-and-kaboodle
    """
    ### figure out network antenna pattern in geographic coordinages
    npix = hp.nside2npix(nside)
    theta, phi = hp.pix2ang(nside, np.arange(npix))
    psi = np.zeros_like(theta, dtype='int')

    horizons = detectors2horizons(detectors, snr_threshold=snr_threshold, flow=flow, fhigh=fhigh) ### compute horizons
    network = horizons2geo_antenna(horizons, theta, phi, psi) ### compute network sensitivity
    network = network**3 ### weight by volume, not range

    ### iterate over simulated times, generating a skymap for each, rotate it to celestial coords
    return [(time, geo_position2geo_skymap(position, (theta, phi), detectors, network, time_delay_error=time_delay_error)) for time, position in zip(draw_times(start, stop, size=size), draw_positions(theta, phi, network, size=size))]

def simulate_cel_skymaps(*args, **kwargs):
    """
    the whole kit-and-kaboodle. Delegates to simulate_geo_skymaps and then rotates the results into celestial coordinates
    """
    return [(time, triangulate.rotateMapE2C(skymap, time)) for time, skymap in simulate_geo_skymaps(*args, **kwargs)]
