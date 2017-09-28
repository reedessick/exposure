__description__ = "general utilities"
__author__ = "Reed Essick (reed.essick@ligo.org)"

#-------------------------------------------------

import numpy as np
from pylal import Fr

#-------------------------------------------------

def extract_start_dur(path, suffix='.gwf'):
    return [int(_) for _ in path[:-len(suffix)].split('-')[-2:]]

def vec_from_frames(frames, channel, start, stop, verbose=False):
        """
        returns a numpy array of the data inculded in frames between start and stop
        CURRENTLY ASSUME CONTIGUOUS DATA, but we should check this

        meant to be used with files_from_cache
        """
        vecs = []
        dt = 0
        for frame, frame_start, frame_dur in frames:
                if verbose: 
                    print( frame )
                s = max(frame_start, start)
                d = min(frame_start+frame_dur,stop) - s
                vec, gpstart, offset, dt, _, _ = Fr.frgetvect1d(frame, channel, start=s, span=d)
                vecs.append( vec )
        vec = np.concatenate(vecs)
        return vec, dt

def andsegments(list1, list2):
    """
    computes the intersection of 2 lists of segments, returning a new list
    taken from laldetchar.idq.event.andtwosegmentlists
    """
    newsegments = []
    index = 0
    for seg in list1:
        while index > 0 and list2[index][0] > seg[0]:
            index -= 1
        while index < len(list2) and list2[index][1] <= seg[0]:
            index += 1
        while index < len(list2) and list2[index][0] < seg[1]:
            newsegments.append([max(seg[0], list2[index][0]), min(seg[1], list2[index][1])])
            index += 1
        if index > 0:
            index -= 1
    return newsegments

def invsegments(start, stop, segs):
    """
    takes the inverse of a list of segments between start and stop
    assumes that all segments are already between start and stop
    """
    newsegments = []
    for s, e in segs:
        if start < s <= stop:
            newsegments.append( [start, s] )
        else:
            start = e
    if e < stop:
        newsegments.append( [e, stop] )
    return newsegments
