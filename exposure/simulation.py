__doc__ = "a module housing useful routines for simulation"
__author__ = "reed.essick@ligo.org"

#-------------------------------------------------

import numpy as np
import healpy as hp

### non-standard libraries
from skymap_statistics import triangulate

#-------------------------------------------------

PI_2 = np.pi*0.5

DEFAULT_SIZE = 100
DEFAULT_NSIDE = 64

DEFAULT_SNR_THRESHOLD = 8.
DEFAULT_FLOW = 32. ### Hz
DEFAULT_FHIGH = 1024. ### Hz

DEFAULT_TIME_ERROR = 1.e-4 ### sec
DEFAULT_SEGLEN = 600. ### sec

#-------------------------------------------------

def time2sec_ns(time):
    sec = int(time)
    ns = int(round((time-sec)*1e9, 0))
    return sec, ns

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
        time_error=DEFAULT_TIME_ERROR,
    ):
    """
    compute modulated triangulation rings given the detectors and network (antenna pattern)
    """
    ### add in triangulation rings
    source_n = hp.ang2vec(source_theta, source_phi)
    map_n = hp.ang2vec(map_theta, map_phi)

    ### generate time-of-arrivals based on Gaussian distribs
    ### compute the likelihood after marginalizing over the time-of-arrival at geocenter
    ### logL = -0.5*sum( (t_i - (tg + n*r_i))**2/error_i**2 )  --> and then marginalize over tg
    ### marg logL = -0.5*sum((t_i-n*r_i)**2/error_i**2) + 0.5*(sum((t_i-n*r_i)/error_i**2))**2/sum(error_i**2)

    aprefact = 1./time_error**2
    a = np.zeros_like(map_theta, dtype=float)

    bprefact = 0.5/time_error**2
    b = np.zeros_like(map_theta, dtype=float)

    for detector in detectors:
        time = np.sum(source_n*detector.dr)+np.random.randn()*time_error
        map_time = np.sum(map_n*detector.dr, axis=1)

        a += aprefact*(time - map_time)
        b -= bprefact*(time - map_time)**2

    ### construct skymap as the sum of its parts
    logskymap = +0.5*(a**2)/(len(detectors)*aprefact) + b

    ### multiply by antenna patterns to modulate triangulation rings
    logskymap += np.log(network)

    ### return a normalized skymap
    logskymap -= np.max(logskymap)
    return np.exp(logskymap)/np.sum(np.exp(logskymap))

#-------------------------------------------------

def simulate_geo_exposure(
        start,
        stop,
        detectors,
        nside=DEFAULT_NSIDE,
        snr_threshold=DEFAULT_SNR_THRESHOLD,
        flow=DEFAULT_FLOW,
        fhigh=DEFAULT_FHIGH,
    ):
    """
    simulate exposure in geographic coordinates
    Note, we ignore start, stop when actually calculculating exposure because
        the antenna patterns are simulated as stationary in this frame (stationary PSDs)
    """
    npix = hp.nside2npix(nside)
    theta, phi = hp.pix2ang(nside, np.arange(npix))
    psi = np.zeros_like(theta, dtype='int')

    horizons = detectors2horizons(detectors, snr_threshold=snr_threshold, flow=flow, fhigh=fhigh)
    network = horizons2geo_antenna(horizons, theta, phi, psi)
    network = network **3 ### weight by volume, not range
    network *= (stop-start) ### weight by time

    start_sec, start_ns = time2sec_ns(start)
    stop_sec, stop_ns = time2sec_ns(stop)

    head = {
        'gps_start_sec':start_sec,
        'gps_start_nsec':start_ns,
        'gps_stop_sec':stop_sec,
        'gps_stop_nsec':stop_ns,
        'ifos':[detector.name for detector in detectors],
    }
    return head, network

def simulate_cel_exposure(
        start,
        stop,
        detectors,
        nside=DEFAULT_NSIDE,
        snr_threshold=DEFAULT_SNR_THRESHOLD,
        flow=DEFAULT_FLOW,
        fhigh=DEFAULT_FHIGH,
        seglen=DEFAULT_SEGLEN
    ):
    """
    simulate the exposure in celestial coordinates
    more expensive than simulate_geo_exposure because we re-compute antenna patterns at each time-step

    NOTE: this does NOT account for leap-seconds
    """
    npix = hp.nside2npix(nside)
    theta, ra = hp.pix2ang(nside, np.arange(npix))
    psi = np.zeros_like(theta, dtype='int')

    horizons = detectors2horizons(detectors, snr_threshold=snr_threshold, flow=flow, fhigh=fhigh)

    exposure = np.zeros(npix, dtype=float)
    while start < stop:
        dt = min(stop, start+seglen) - start
        t = 0.5*dt + start

        network = horizons2geo_antenna(horizons, theta, triangulate.rotateRAC2E(ra, t), psi)
        network = network**3 ### weigh by volume
        network *= dt ### weigh by time

        exposure += network
        start += dt

    start_sec, start_ns = time2sec_ns(start)
    stop_sec, stop_ns = time2sec_ns(stop)

    head = {
        'gps_start_sec':start_sec,
        'gps_start_nsec':start_ns,
        'gps_stop_sec':stop_sec,
        'gps_stop_nsec':stop_ns,
        'ifos':[detector.name for detector in detectors],
    }
    return head, exposure

def simulate_geo_skymaps(
        start,
        stop,
        detectors,
        size=DEFAULT_SIZE,
        nside=DEFAULT_NSIDE,
        snr_threshold=DEFAULT_SNR_THRESHOLD,
        flow=DEFAULT_FLOW,
        fhigh=DEFAULT_FHIGH,
        time_error=DEFAULT_TIME_ERROR,
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

    ifos = [detector.name for detector in detectors]

    ### iterate over simulated times, generating a skymap for each, rotate it to celestial coords
    skymaps = []
    for time, (_theta, _phi) in zip(draw_times(start, stop, size=size), draw_positions(theta, phi, network, size=size)):
        sec, ns = time2sec_ns(time)
        head = {
            'gps_sec':sec,
            'gps_nsec':ns,
            'theta':_theta,
            'phi':_phi,
            'dec':PI_2 - _theta,
            'ra':triangulate.rotateRAE2C(_phi, time),
            'ifos':ifos,
        }
        skymaps.append( (head, geo_position2geo_skymap((_theta, _phi), (theta, phi), detectors, network, time_error=time_error)) )

    return skymaps

def simulate_cel_skymaps(
        start,
        stop,
        detectors,
        size=DEFAULT_SIZE,
        nside=DEFAULT_NSIDE,
        snr_threshold=DEFAULT_SNR_THRESHOLD,
        flow=DEFAULT_FLOW,
        fhigh=DEFAULT_FHIGH,
        time_error=DEFAULT_TIME_ERROR,
    ):
    """
    the whole kit-and-kaboodle. Delegates to simulate_geo_skymaps and then rotates the results into celestial coordinates
    more expensive than geo_skymaps because we have to re-compute network for each event...

    NOTE: this does NOT account for leap-seconds
    """
    ### figure out network antenna pattern in geographic coordinages
    npix = hp.nside2npix(nside)
    theta, ra = hp.pix2ang(nside, np.arange(npix))
    psi = np.zeros_like(theta, dtype='int')

    horizons = detectors2horizons(detectors, snr_threshold=snr_threshold, flow=flow, fhigh=fhigh) ### compute horizons

    ifos = [detector.name for detector in detectors]

    skymaps = []
    for time in draw_times(start, stop, size=size):
        phi = triangulate.rotateRAC2E(ra, time)
        network = horizons2geo_antenna(horizons, theta, phi, psi) ### compute network sensitivity
        network = network**3 ### weight by volume, not range

        _theta, _phi = draw_positions(theta, phi, network, size=1)[0]
        sec, ns = time2sec_ns(time)
        head = {
            'gps_sec':sec,
            'gps_nsec':ns,
            'theta':_theta,
            'phi':_phi,
            'dec':PI_2 - _theta,
            'ra':triangulate.rotateRAE2C(_phi, time),
            'ifos':ifos,
        }
        skymaps.append( (head, geo_position2geo_skymap((_theta, _phi), (theta, phi), detectors, network, time_error=time_error)) )

    return skymaps
