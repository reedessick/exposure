__description__ = "general utilities"
__author__ = "Reed Essick (reed.essick@ligo.org)"

#-------------------------------------------------

import gzip
import numpy as np

#-------------------------------------------------

### define template strings here for easy reference between multiple scripts!
compute_horizon_sub = '''\
universe = %(universe)s
executable = %(exe)s
arguments = "$(path) --output-dir $(outdir) $(tag) --flow %(flow)d --fhigh %(fhigh)d --distance %(distance)f --snr-thr %(snr_thr)f -v"
getenv = true
accounting_group = %(accounting_group)s
accounting_group_user = %(accounting_group_user)s
log    = $(outdir)/condor-compute_horizon%(filetag)s_$(Cluster)-$(Process).log
error  = $(outdir)/condor-compute_horizon%(filetag)s_$(Cluster)-$(Process).err 
output = $(outdir)/condor-compute_horizon%(filetag)s_$(Cluster)-$(Process).out
notification = never
queue 1'''

compute_horizon_dag = '''\
JOB    compute_horizon_%(jobid)s %(sub)s
VARS   compute_horizon_%(jobid)s path="%(path)s" tag="%(tag)s" outdir="%(outdir)s"
RETRY  compute_horizon_%(jobid)s %(retry)d
'''

compute_network_sensitivity_sub = '''\
universe = %(universe)s
executable = %(exe)s
arguments = "$(gps) $(ifo_horizons) --nside %(nside)d --output-dir $(outdir) --tag $(tag) -V"
getenv = true
accounting_group = %(accounting_group)s
accounting_group_user = %(accounting_group_user)s
log    = $(outdir)/condor-compute_network_sensitivity%(filetag)s_$(Cluster)-$(Process).log
error  = $(outdir)/condor-compute_network_sensitivity%(filetag)s_$(Cluster)-$(Process).err 
output = $(outdir)/condor-compute_network_sensitivity%(filetag)s_$(Cluster)-$(Process).out
notification = never
queue 1'''

compute_network_sensitivity_dag = '''\
JOB    compute_network_sensitivity_%(jobid)s %(sub)s
VARS   compute_network_sensitivity_%(jobid)s gps="%(gps)d" tag="%(tag)s" outdir="%(outdir)s" ifo_horizons="%(ifo_horizons)s"
RETRY  compute_network_sensitivity_%(jobid)s %(retry)d
'''

compute_psd_sub = '''\
universe = %(universe)s
executable = %(exe)s
arguments = "%(channel)s %(frametype)s $(gpsstart) $(gpsstop) --seglen %(seglen)d --overlap %(overlap)d --tukey-alpha %(tukey_alpha)f --output-dir $(outdir) %(tag)s -V"
getenv = true
accounting_group = %(accounting_group)s
accounting_group_user = %(accounting_group_user)s
log    = $(outdir)/condor-compute_psd%(filetag)s_$(Cluster)-$(Process).log
error  = $(outdir)/condor-compute_psd%(filetag)s_$(Cluster)-$(Process).err 
output = $(outdir)/condor-compute_psd%(filetag)s_$(Cluster)-$(Process).out
notification = never
queue 1'''

compute_psd_dag = '''\
JOB    compute_psd_%(jobid)s %(sub)s
VARS   compute_psd_%(jobid)s gpsstart="%(gpsstart)d" gpsstop="%(gpsstop)d" outdir="%(outdir)s"
RETRY  compute_psd_%(jobid)s %(retry)d
'''

plot_psd_sub = '''\
universe = %(universe)s
executable = %(exe)s
arguments = "$(psd) --f-min %(fmin)d --f-max %(fmax)d --y-min %(ymin)e --y-max %(ymax)e --output-dir %(outdir)s %(tag)s -v"
getenv = true
accounting_group = %(accounting_group)s
accounting_group_user = %(accounting_group_user)s
log    = %(outdir)s/condor-compute_psd%(filetag)s_$(Cluster)-$(Process).log
error  = %(outdir)s/condor-compute_psd%(filetag)s_$(Cluster)-$(Process).err 
output = %(outdir)s/condor-compute_psd%(filetag)s_$(Cluster)-$(Process).out
notification = never
queue 1'''

plot_psd_dag = '''\
JOB    plot_psd_%(jobid)s %(sub)s
VARS   plot_psd_%(jobid)s psd="%(path)s"
RETRY  plot_psd_%(jobid)s %(retry)d
'''

compute_exposure_sub = '''\
universe = %(universe)s
executable = %(exe)s
arguments = "$(FITS) $(normalize) --nside %(nside)d --index %(index)d --output-dir %(outdir)s %(tag)s"
getenv = true
accounting_group = %(accounting_group)s
accounting_group_user = %(accounting_group_user)s
log    = %(outdir)s/condor-compute_exposure%(filetag)s_$(Cluster)-$(Process).log
error  = %(outdir)s/condor-compute_exposure%(filetag)s_$(Cluster)-$(Process).err 
output = %(outdir)s/condor-compute_exposure%(filetag)s_$(Cluster)-$(Process).out
notification = never
queue 1'''

compute_exposure_dag = '''\
JOB    compute_exposure_%(jobid)s %(sub)s
VARS   compute_exposure_%(jobid)s FITS="%(FITS)s" normalize="%(normalize)s"
RETRY  compute_exposure_%(jobid)s %(retry)d
'''

plot_maps_sub = '''\
universe = %(universe)s
executable = %(exe)s
arguments = "$(base),$(path) --coord %(coord)s --grid --projection %(projection)s $(tag) -v"
getenv = true
accounting_group = %(accounting_group)s
accounting_group_user = %(accounting_group_user)s
log    = $(outdir)/condor-plot_maps%(filetag)s_$(Cluster)-$(Process).log
error  = $(outdir)/condor-plot_maps%(filetag)s_$(Cluster)-$(Process).err 
output = $(outdir)/condor-plot_maps%(filetag)s_$(Cluster)-$(Process).out
notification = never
queue 1'''

plot_maps_dag = '''\
JOB    compute_horizon_%(jobid)s %(sub)s
VARS   compute_horizon_%(jobid)s path="%(path)s" tag="%(tag)s"
RETRY  compute_horizon_%(jobid)s %(retry)d
'''

#-------------------------------------------------

MOD_STRIDE = 100000

def segs_path(output_dir, tag, gpsstart, gpsdur):
    return "%s/seg%s-%d-%d.txt.gz"%(gps2moddir(directory, gpsstart), tag, (int(gpsstart)/MOD_STRIDE)*MOD_STRIDE, MOD_STRIDE)

def psd_path(output_dir, tag, gpsstart, gpsdur):
    return "%s/psd%s-%d-%d.txt.gz"%(output_dir, tag, gpsstart, gpsdur)

def horizon_path(output_dir, tag, gpstart, gpsdur):
    return "%s/horizon%s-%d-%d.txt.gz"%(output_dir, tag, gpstart, gpsdur)

def sensitivity_path(output_dir, tag, gps, gzip=False):
    ans = "%s/network_sensitivity%s-%d.fits"%(output_dir, tag, gps)
    if gzip:
        ans = ans+".gz"
    return ans

def gps2moddir(directory, start):
    return "%s/%d/%d-%d/"%(directory, int(start)/MOD_STRIDE)

def gps2dir(directory, start, dur):
    return "%s/%d-%d/"%(gps2moddir(directory, start), start, dur)

#-------------------------------------------------

def extract_start_dur(path, suffix='.gwf'):
    return [int(_) for _ in path[:-len(suffix)].split('-')[-2:]]

def report_psd(path, freqs, psd):
    """writes the PSD to disk
    """
    np.savetxt(psdpath, np.array(zip(freqs, psd)), delimiter=' ')

def report_segs(path, new_segs):
    """reads in the segments contained in path and appends the current seg of segs to them
    assumes you are only adding segments in time order, so it will either start a new segment or merge the last segment in the file
    """
    ### read in and merge segs if necessary
    segs = np.loadtxt(path)
    if len(segs):
        if (segs[-1][1]==new_segs[0][0]):
            segs[-1][1] = segs.pop(0)[1] ### remove the first element from segs

        segs = list(segs)+new_segs

    else:
        segs = new_segs

    ### write the result back out
    np.savetxt(path, segs)

#-------------------------------------------------

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
    if not segs:
        return [[start, stop]]

    newsegments = []
    for s, e in segs:
        if start < s <= stop:
            newsegments.append( [start, s] )
        start = e
    if e < stop:
        newsegments.append( [e, stop] )

    return newsegments

def mergesegments(segments):
    """
    assumes segments are in the correct order and non-overlapping.
    simply merges any touching segments into one bigger segemnt
    """
    if len(segments)<2:
        return segments

    segs = []
    s, e = segments[0]
    for S, E in segments[1:]:
        if e==S:
            e=E
        else:
            segs.append( [s,e] )
            s = S
            e = E
    segs.append( [s, e] )
    return segs

def checksegments(seg, window=None):
    """
    ensures segment makes sense
    if provided, segment is sliced so that we only keep the part that interesects with window=[start, end]
    """
    s, e = seg 
    if window:
        start, end = window
        if s < start:
            s = start # truncate the start of seg
        elif not (s < end):
            return False # no overlap between current segment and window
        if e > end:
            e = end # truncate the end of seg
        elif not (e > start):
            return False # no overlap between currnet segment and window

    if s < e:
        return [s,e]
    elif s > e:
        raise ValueError("something is very wrong with segment generation... seg[1]=%.3f < seg[0]=%.3f"%tuple(seg))
    else:
        return False