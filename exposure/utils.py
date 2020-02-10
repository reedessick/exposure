__description__ = "general utilities"
__author__ = "Reed Essick (reed.essick@ligo.org)"

#-------------------------------------------------

import os
import shutil
import gzip

import h5py
import numpy as np

#-------------------------------------------------

DEFAULT_ACCOUNTING_GROUP = 'ligo.sim.o3.cbc.bayesianpopulations.parametric' 
DEFAULT_RETRY = 0

#-------------------------------------------------

### define template strings here for easy reference between multiple scripts!
compute_horizon_sub = '''\
universe = %(universe)s
executable = %(exe)s
arguments = "$(path) --output-dir $(outdir) $(tag) --flow %(flow)d --fhigh %(fhigh)d --distance %(distance)f --snr-thr %(snr_thr)f -v"
getenv = true
accounting_group = %(accounting_group)s
accounting_group_user = %(accounting_group_user)s
log    = $(outdir)/log/condor-compute_horizon%(filetag)s_$(Cluster)-$(Process).log
error  = $(outdir)/log/condor-compute_horizon%(filetag)s_$(Cluster)-$(Process).err 
output = $(outdir)/log/condor-compute_horizon%(filetag)s_$(Cluster)-$(Process).out
notification = never
queue 1'''

compute_horizon_jobid = "compute_horizon_%s"
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
log    = $(outdir)/log/condor-compute_network_sensitivity%(filetag)s_$(Cluster)-$(Process).log
error  = $(outdir)/log/condor-compute_network_sensitivity%(filetag)s_$(Cluster)-$(Process).err 
output = $(outdir)/log/condor-compute_network_sensitivity%(filetag)s_$(Cluster)-$(Process).out
notification = never
queue 1'''

compute_network_sensitivity_jobid = "compute_network_sensitivity_%s"
compute_network_sensitivity_dag = '''\
JOB    compute_network_sensitivity_%(jobid)s %(sub)s
VARS   compute_network_sensitivity_%(jobid)s gps="%(gps)d" tag="%(tag)s" outdir="%(outdir)s" ifo_horizons="%(ifo_horizons)s"
RETRY  compute_network_sensitivity_%(jobid)s %(retry)d
'''

compute_psd_sub = '''\
universe = %(universe)s
executable = %(exe)s
arguments = "%(channel)s %(frametype)s $(gpsstart) $(gpsstop) --seglen %(seglen)d --overlap %(overlap)d --tukey-alpha %(tukey_alpha)f --suffix %(suffix)s --output-dir $(outdir) %(tag)s -V"
getenv = true
accounting_group = %(accounting_group)s
accounting_group_user = %(accounting_group_user)s
log    = $(outdir)/log/condor-compute_psd%(filetag)s_$(Cluster)-$(Process).log
error  = $(outdir)/log/condor-compute_psd%(filetag)s_$(Cluster)-$(Process).err 
output = $(outdir)/log/condor-compute_psd%(filetag)s_$(Cluster)-$(Process).out
notification = never
queue 1'''

compute_psd_jobid = "compute_psd_%s"
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
log    = %(outdir)s/log/condor-compute_psd%(filetag)s_$(Cluster)-$(Process).log
error  = %(outdir)s/log/condor-compute_psd%(filetag)s_$(Cluster)-$(Process).err 
output = %(outdir)s/log/condor-compute_psd%(filetag)s_$(Cluster)-$(Process).out
notification = never
queue 1'''

plot_psd_jobid = "plot_psd_%s"
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
log    = %(outdir)s/log/condor-compute_exposure%(filetag)s_$(Cluster)-$(Process).log
error  = %(outdir)s/log/condor-compute_exposure%(filetag)s_$(Cluster)-$(Process).err 
output = %(outdir)s/log/condor-compute_exposure%(filetag)s_$(Cluster)-$(Process).out
notification = never
queue 1'''

compute_exposure_jobid = " compute_exposure_%s"
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
log    = $(outdir)/log/condor-plot_maps%(filetag)s_$(Cluster)-$(Process).log
error  = $(outdir)/log/condor-plot_maps%(filetag)s_$(Cluster)-$(Process).err 
output = $(outdir)/log/condor-plot_maps%(filetag)s_$(Cluster)-$(Process).out
notification = never
queue 1'''

plot_maps_jobid = "compute_horizon_%s"
plot_maps_dag = '''\
JOB    compute_horizon_%(jobid)s %(sub)s
VARS   compute_horizon_%(jobid)s path="%(path)s" tag="%(tag)s"
RETRY  compute_horizon_%(jobid)s %(retry)d
'''

monte_carlo_vt_sub = '''\
universe = %(universe)s
executable = %(exe)s
arguments = "--verbose %(population)s $(gpsstart) $(gpsstop) --snr-threshold %(snr_threshold).6f --min-num-samples %(min_num_samples)d --fractional-systematic-error %(error).6f --psd-suffix %(psd_suffix)s --samples-suffix %(samples_suffix)s $(sources) --output-dir $(outdir) %(tag)s"
getenv = true
accounting_group = %(accounting_group)s
accounting_group_user = %(accounting_group_user)s
log    = $(outdir)/log/condor-monte-carlo-vt%(filetag)s_$(Cluster)-$(Process).log
error  = $(outdir)/log/condor-monte-carlo-vt%(filetag)s_$(Cluster)-$(Process).err 
output = $(outdir)/log/condor-monte-carlo-vt%(filetag)s_$(Cluster)-$(Process).out
notification = never
queue 1'''

monte_carlo_vt_jobid = "monte_carlo_vt_%s"
monte_carlo_vt_dag = '''\
JOB    monte_carlo_vt_%(jobid)s %(sub)s
VARS   monte_carlo_vt_%(jobid)s gpsstart="%(gpsstart)d" gpsstop="%(gpsstop)d" outdir="%(outdir)s" sources="%(sources)s"
RETRY  monte_carlo_vt_%(jobid)s %(retry)d
'''

compute_psd_cdf_tag = '$(day)%(filetag)s'
compute_psd_cdf_sub = '''\
universe = %(universe)s
executable = %(exe)s
arguments = "$(psds) --output-dir %(outdir)s --tag $(day)%(filetag)s --num-points %(num_points)d --verbose"
getenv = true
accounting_group = %(accounting_group)s
accounting_group_user = %(accounting_group_user)s
log    = %(outdir)s/log/condor-compute-psd-cdf%(filetag)s_$(Cluster)-$(Process).log
error  = %(outdir)s/log/condor-compute-psd-cdf%(filetag)s_$(Cluster)-$(Process).err 
output = %(outdir)s/log/condor-compute-psd-cdf%(filetag)s_$(Cluster)-$(Process).out
notification = never
queue 1'''

compute_psd_cdf_jobid = 'compute_psd_cdf_%s'
compute_psd_cdf_dag = '''\
JOB %(job)s %(sub)s
VARS %(job)s psds="%(psds)s" day="%(day)s"
RETRY %(job)s %(retry)d
'''

combine_psd_cdf_sub = '''\
universe = %(universe)s
executable = %(exe)s
arguments = "$(psdcdfs) --output-dir %(outdir)s --tag %(tag)s --num-points %(num_points)d --verbose"
getenv = true
accounting_group = %(accounting_group)s
accounting_group_user = %(accounting_group_user)s
log    = %(outdir)s/log/condor-combine-psd-cdf%(filetag)s_$(Cluster)-$(Process).log
error  = %(outdir)s/log/condor-combine-psd-cdf%(filetag)s_$(Cluster)-$(Process).err 
output = %(outdir)s/log/condor-combine-psd-cdf%(filetag)s_$(Cluster)-$(Process).out
notification = never
queue 1'''

combine_psd_cdf_jobid = 'combine_psd_cdf_%s'
combine_psd_cdf_dag = '''
JOB %(job)s %(sub)s
VARS %(job)s psdcdfs="%(psdcdfs)s"
RETRY %(job)s %(retry)d
'''

#-------------------------------------------------

MOD_STRIDE = 100000

def segs_path(output_dir, tag, gpsstart, gpsdur):
    return "%s/seg%s-%d-%d.txt.gz"%(gps2moddir(output_dir, gpsstart), tag, (int(gpsstart)/MOD_STRIDE)*MOD_STRIDE, MOD_STRIDE)

def psd_path(output_dir, tag, gpsstart, gpsdur, suffix='csv.gz'):
    return "%s/%s"%(gps2dir(output_dir, gpsstart, gpsdur), psd_basename(tag, gpsstart, gpsdur, suffix=suffix))

def psd_basename(tag, gpsstart, gpsdur, suffix='csv.gz'):
    return "psd%s-%s-%s.%s"%(tag, gpsstart, gpsdur, suffix)

def psd_cdf_path(output_dir, tag):
    return "%s/psdcdf%s.hdf5"%(output_dir, tag)

def samples_path(output_dir, tag, gpsstart, gpsdur, suffix='csv.gz'):
    return "%s/samples%s-%d-%d.%s"%(gps2dir(output_dir, gpsstart, gpsdur), tag, gpsstart, gpsdur, suffix)

def horizon_path(output_dir, tag, gpstart, gpsdur):
    return "%s/horizon%s-%d-%d.txt.gz"%(output_dir, tag, gpstart, gpsdur)

def sensitivity_path(output_dir, tag, gps, gzip=False):
    ans = "%s/network_sensitivity%s-%d.fits"%(output_dir, tag, gps)
    if gzip:
        ans = ans+".gz"
    return ans

def gps2moddir(directory, start):
    return "%s/%d/"%(directory, int(start)/MOD_STRIDE)

def gps2dir(directory, start, dur):
    return "%s/%d-%d/"%(gps2moddir(directory, start), start, dur)

#-------------------------------------------------

def extract_start_dur(path, suffix='.gwf'):
    return [int(_) for _ in path[:-len(suffix)].split('-')[-2:]]

#------------------------
### logic for reading/writing PSDs

def report_psd(path, freqs, psd, tmpdir='/tmp'):
    """writes the PSD to disk
    """
    ### save to a temporary file
    tmp = os.path.join(tmpdir, '.'+os.path.basename(path))
    if path.endswith('csv') or path.endswith('csv.gz'):
        _write_psd_csv(tmp, freqs, psd)
    elif path.endswith('hdf5') or path.endswith('h5'):
        _write_psd_hdf5(tmp, freqs, psd)
    else:
        raise ValueError('suffix not recognized for PSD path %s'%path)

    ### move to the final location
    shutil.move(tmp, path)

def _write_psd_csv(path, freqs, psd):
    np.savetxt(path, np.array(zip(freqs, psd)), comments='', delimiter=',', header='frequency,psd')

def _write_psd_hdf5(path, freqs, psd, name='samples'):
    with h5py.File(path, 'w') as obj:
        obj.create_dataset(name, data=np.array(zip(freqs, psd), dtype=[('frequency',float), ('psd',float)]))

def retrieve_psd(path):
    """read the PSD from disk
    return freqs, psd
    """
    if path.endswith('csv') or path.endswith('csv.gz'):
        return _read_psd_csv(path)
    elif path.endswith('hdf5') or path.endswith('h5'):
        return _read_psd_hdf5(tmp, freqs, psd)
    else:
        raise ValueError('suffix not recognized for PSD path %s'%path)

def _read_psd_csv(path):
    ans = np.genfromtxt(path, delimiter=',', names=True)
    return ans['frequency'], ans['psd']

def _read_psd_hdf5(path, name='samples'):
    with h5py.File(path, 'r') as obj:
        freq = obj[name]['frequency']
        psd = obj[name]['psd']
    return freq, psd

#------------------------
### logic for reading/writing PSD CDFs

def report_psd_cdf(path, freqs, vals, cdfs, Npsd):
    """interact with hdf5 file format for marginal CDFs for a set of PSDs.
Expect vals and cdfs to have the shape (Nfrq, Npts)
    """
    with h5py.File(path, 'w') as obj:
        group = obj.create_group('PSD_CDF')
        group.attrs.create('num_psds', Npsd)
        group.create_dataset('frequencies', data=freqs, dtype=float)

        Nfrq, Npts = np.shape(vals)
        assert len(freqs)==Nfrq

        data = np.empty((Nfrq, 2, Npts), dtype=float)
        data[:,0,:] = vals
        data[:,1,:] = cdfs

        group.create_dataset('CDFs', data=data, dtype=float)

def retrieve_psd_cdf(path):
    """interact with hdf5 file format for marginal CDFs for a set of PSDs"""
    with h5py.File(path, 'r') as obj:
        group = obj['PSD_CDF']
        Npsd = group.attrs['num_psds']
        freqs = group['frequencies'][...]
        data = group['CDFs'][...]

    vals = data[:,0,:]
    cdfs = data[:,1,:]

    return freqs, vals, cdfs, Npsd

#------------------------
### logic for reading/wirting segments

def report_segs(path, new_segs):
    """reads in the segments contained in path and appends the current seg of segs to them
    assumes you are only adding segments in time order, so it will either start a new segment or merge the last segment in the file
    """
    if len(new_segs): ### there must be something to report...
        new_segs = [(a,b) for a, b in new_segs] ### make a copy so we don't mess up a shared reference

        ### read in and merge segs if necessary
        if os.path.exists(path):
            segs = retrieve_segs(path)
            if len(segs):
                if np.ndim(segs)==1: ### a single segment
                    segs = [segs]

                while len(segs) and (segs[-1][0] > new_segs[0][1]): ### remove existing segs that are redundant with where new_segs starts
                    segs.pop(-1)

                if len(segs):
                    if (segs[-1][1] >= new_segs[0][0]): ### these things overlap
                        segs[-1][1] = new_segs.pop(0)[1] ### remove the first element from segs

                segs = segs+list(new_segs)

            else: ### either empty or a single se
                segs = new_segs

        else:
            segs = new_segs

        ### write the result back out to a temporary location
        tmp = os.path.join(os.path.dirname(path), '.'+os.path.basename(path))
        np.savetxt(tmp, segs, fmt='%d')
        ### move to final location
        shutil.move(tmp, path)

def retrieve_segs(path):
    """retrieves segments from disk
    """
    return list(np.loadtxt(path))

#------------------------
### logic for reading/writing horizon estimates

def report_horizon(path, horizon, flow, fhigh, psd_path):
    """report the horizon to disk
    """
    with open(path, 'w') as f:
        f.write('%.9e\n%d\n%d\n%s'%(horizon, flow, fhigh, psd_path))

def retrieve_horizon(path):
    """retrieve the horizon from disk
    """
    with open(path, 'r') as f:
        horizon = float(f.readline().strip())
        flow = float(f.readline().strip())
        fhigh = float(f.readlin().strip())
        psd_path = f.readline().strip()
    return horizon, flow, fhigh, psd_path

#------------------------
### logic for reading/writing sensitivities

def report_sensitivity(path, sensitivity, *extra_header):
    basepath = path.strip('.gz')
    hp.write_map(
        basepath,
        network,
        coord='C',
        extra_header=extra_header,
    )

    if path.endswith('.gz'):
        if opts.verbose:
            print( "gzipping : %s -> %s"%(basepath, path) )
        with open(basepath, 'r') as f:
            with gzip.open(path, 'w') as gzf:
                gzf.write(f.read())

#-------------------------------------------------

def livetime(segs):
    return np.sum(e-s for s, e in segs)

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
