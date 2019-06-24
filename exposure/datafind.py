__doc__ = ""
__author__ = "Reed Essick (reed.essick@ligo.org)"

#-------------------------------------------------

import glob
import time

import numpy as np

import subprocess as sp

### non-standard libraries
try:
    from lalframe import frread
except ImportError:
    frread = None

try:
    from lal.gpstime import tconvert as lal_tconvert
except ImportError:
    lal_tconvert = None

from . import utils

#-------------------------------------------------

DEFAULT_STRIDE = 60
DEFAULT_DELAY = 0
DEFAULT_MAX_LATENCY = 2*DEFAULT_STRIDE
DEFAULT_CADENCE = 1e-3

#-------------------------------------------------

def tconvert(arg='now'):
    if lal_tconvert is None:
        raise ImportError('could not import lal.gpstime.tconvert!')
    return float(lal_tconvert(arg))

def latency(target, delay=DEFAULT_DELAY):
    return tconvert('now')-delay - target

def wait(target, delay=DEFAULT_DELAY, logger=None):
    w = max(0, -latency(target, delay=delay))
    if logger is not None:
        logger.info('sleeping for %.3f sec'%w)
    sleep(w)

def sleep(dt):
    if dt > 0:
        time.sleep(dt)

#-------------------------------------------------

### SEGMENT DISCOVERY

query_cmd = "ligolw_segment_query_dqsegdb -q -t https://segments.ligo.org -a %(flag)s -s %(gpsstart)d -e %(gpsstop)d"
print_cmd = "ligolw_print -c start_time -c end_time -t segment".split()
def query_flag(flag, start, stride, verbose=False):
    ### FIXME: use the SegDb Python API directly instead of this hack...
    cmd = query_cmd%{
        'flag':flag,
        'gpsstart':start,
        'gpsstop':start+stride,
    }
    if verbose:
        print( cmd )
    segs = sp.Popen(cmd.split(), stdout=sp.PIPE).communicate()[0]
    return [[int(_) for _ in seg.split(',')] for seg in sp.Popen(print_cmd, stdin=sp.PIPE, stdout=sp.PIPE).communicate(segs)[0].strip('\n').split('\n') if seg]

def include_flags(segments, flags, start, stride, verbose=False):
    for flag in flags:
        segments = utils.andsegments(segments, query_flag(flag, start, stride, verbose=verbose))

    return segments

def exclude_flags(segments, flags, start, stride, verbose=False):
    if not segments:
        return segments

    gpsstart = segments[0][0]
    gpsstop = segments[-1][1]
    for flag in flags:
        segments = utils.andsegments(segments, invsegments(gpsstart, gpsstop, query_flag(flag, start, stride, verbose=verbose)))

    return segments

#-------------------------------------------------

### DATA EXTRACTION AND MANIPULATION

def vec_from_frames(frames, channel, start, stop, verbose=False):
    """
    returns a numpy array of the data inculded in frames between start and stop
    CURRENTLY ASSUME CONTIGUOUS DATA, but we should check this
    """
    if frread is None:
        raise ImportError('could not import lalframe.frread!')
    vecs = []
    dt = 0
    for frame, frame_start, frame_dur in frames:
        if verbose:
            print( frame )
        s = max(frame_start, start)
        d = min(frame_start+frame_dur,stop) - s
        if d > 0:
            out = frread.read_timeseries(frame, channel, start=s, duration=d)
            vecs.append( out.data.data )
            dt = out.deltaT
    vec = np.concatenate(vecs)
    return vec, dt

def extract_scisegs(frames, channel, bitmask, start, stride, verbose=False):
    """
    extract scisegs from channel in frames using bitmask
    """
    if not frames: ### empty list, so no segments
        return []

    ### extract vectors and build segments
    segset = []
    vec, dt = vec_from_frames(frames, channel, start, start+stride, verbose=verbose)
    n = len(vec)

    ### build time vector        add starting time
    t = np.arange(0, dt*n, dt) + start

    ### determine whether state acceptable
    ### add "False" buffers to set up the computation of start and end time
    state = np.concatenate( ([False], (vec >> bitmask) & 1, [False])) ### bitwise operation

    ### determine beginning of segments
    ###      i=False      i+1 = True  strip the trailing buffer
    b = ( (1-state[:-1])*(state[1:]) )[:-1].astype(bool)
    b = t[b] ### select out times

    ### determine end of segments
    ###     i=True     i+1=False      strip the leading buffer
    e = ( (state[:-1])*(1-state[1:]) )[1:].astype(bool)
    e = t[e] + dt ### select out times
                  ### extra dt moves these markers to the end of segments

    ### stitch together start and end times, append to global list
    segset = list( np.transpose( np.array( [b, e] ) ) )

    if not segset: ### empty list
        return []

    ### clean up segs!
    return andsegments(mergesegments(segset), [(start, start+stride)])

#------------------------

### FRAME DISCOVERY

gw_data_find_cmd = "gw_data_find -o %(observatory)s --type %(frametype)s --url file -s %(gpsstart)d -e %(gpsstop)d"
def gw_data_find(ifo, ldr_type, start, stride, verbose=False):
    """a routine to automate calls to gw_data_find to get frame locations, starts, and durations
    """
    cmd = gw_data_find_cmd%{
        'observatory':ifo,
        'frametype':ldr_type,
        'gpsstart':start,
        'gpsend':start+stride,
    }
    if verbose:
        print( cmd )
    frames = []
    for frame in [frame for frame in sp.Popen(cmd.split(), stdout=sp.PIPE).communicate()[0].replace('file://localhost','').strip().split('\n') if frame]:
        if verbose:
            print( frame )
        s, d = utils.extract_start_dur(frame)
        frames.append( (frame, s, d) )
    return frames

shm_glob_tmp = "%s/%s1/%s-%s-*-*.gwf"
def shm_data_find(ifo, ldr_type, start, stride, directory='.', verbose=False):
    """a routine to automate discovery of frames within /dev/shm
    """
    end = start+stride
    frames = []
    for frame in sorted(glob.glob(shm_glob_tmp%(directory, ifo, ifo, ldr_type))):
          s, d = utils.extract_start_dur(frame, suffix=".gwf")
          if (s <= end) and (s+d > start): ### there is some overlap!
                frames.append( (frame, s, d) )

    return frames

def coverage(frames, start, stride):
    """
    determines the how much of [start, start+stride] is covered by these frames
    assumes non-overlapping frames!
    """
    ### generate segments from frame names
    segs = [(s, s+d) for _, s, d in sorted(frames)]

    ### check whether segments overlap with desired time range
    covered = 1.0*stride

    end = start + stride
    for s, e in segs:
        if (s < end) and (start < e): ### at least some overlap
            covered -= min(e, end) - max(s, start) ### subtract the overlap

        if covered <= 0:
            break

    return 1 - covered/stride ### return fraction of coverage

def extract_scisegs(frames, channel, bitmask, start, stride):
    """
    extract scisegs from channel in frames using bitmask
    """
    if not frames: ### empty list, so no segments
        return []

    ### extract vectors and build segments
    segset = []
    for frame in frames:
        ### extract the vector from the frame
        vect, s, ds, dt, xunit, yunit = Fr.frgetvect1d(frame, channel)
        n = len(vect)

        ### build time vector        add starting time
        t = np.arange(0, dt*n, dt) + s+ds

        ### determine whether state acceptable
        ### add "False" buffers to set up the computation of start and end time
#       state = np.concatenate( ([False], vect == bitmask, [False])) ### exact integer match
        state = np.concatenate( ([False], (vect >> bitmask) & 1, [False])) ### bitwise operation

        ### determine beginning of segments
        ###      i=False      i+1 = True  strip the trailing buffer
        b = ( (1-state[:-1])*(state[1:]) )[:-1].astype(bool)
        b = t[b] ### select out times

        ### determine end of segments
        ###     i=True     i+1=False      strip the leading buffer
        e = ( (state[:-1])*(1-state[1:]) )[1:].astype(bool)
        e = t[e] + dt ### select out times
                              ### extra dt moves these markers to the end of segments

        ### stitch together start and end times, append to global list
        segset += list( np.transpose( np.array( [b, e] ) ) )

    if not segset: ### empty list
        return []

    ### clean up segs!
    segs = []
    seg1 = segset[0]
    for seg2 in segset[1:]:
        if seg1[1] == seg2[0]:
            seg1[1] = seg2[1] ### join the segments
        else:
            ### check segment for sanity
            append_seg = utils.check_seg( seg1, window=(start, start+stride) )
            if append_seg:
                segs.append( append_seg )
            seg1 = seg2
    append_seg = utils.check_seg( seg1, window=(start, start+stride) )
    if append_seg:
        segs.append( append_seg )

    ### return final list of lists!
    return segs
