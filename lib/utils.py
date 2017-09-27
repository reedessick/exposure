__description__ = "general utilities"
__author__ = "Reed Essick (reed.essick@ligo.org)"

#-------------------------------------------------

import numpy as np
from pylal import Fr

#-------------------------------------------------

def vec_from_frames(frames, start, stop, verbose=False):
        """
        returns a numpy array of the data inculded in frames between start and stop
        CURRENTLY ASSUME CONTIGUOUS DATA, but we should check this

        meant to be used with files_from_cache
        """
        vecs = []
        dt = 0
        for frame, strt, dur in frames:
                if verbose: 
                    print( frame )
                s = max(strt, start)
                d = min(start+dur,stop) - s
                vec, gpstart, offset, dt, _, _ = Fr.frgetvect1d(frame, ifo_chan, start=s, span=d)
                vecs.append( vec )
        vec = np.concatenate(vecs)
        return vec, dt

