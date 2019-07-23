"""a module that implements Monte-Carlo logic and supporting utilties when defining a population
"""
__author__ = "Reed Essick (reed.essick@gmail.com)"

#-------------------------------------------------

import os
import glob

import numpy as np

from . import utils
from gw_event_gen import eventgen

try:
    from ConfigParser import ConfigParser
except ImportError: ### Python2 vs Python3
    from configparser import ConfigParser

#-------------------------------------------------

DEFAULT_MIN_NUM_SAMPLES = 100
DEFAULT_ERROR = np.infty

#-------------------------------------------------

### basic config parsing

def subclasses(klass):
    ans = dict()
    for obj in klass.__subclasses__():
        ans[obj.__name__] = obj
        ans.update(subclasses(obj))
    return ans

def parse(section, config):
    ans = dict()
    for k in config.options(section):
        try: ### maybe it's an integer
            ans[k] = config.getint(section, k)
            continue
        except ValueError:
            pass

        try: ### maybe a float
            ans[k] = config.getfloat(section, k)
            continue
        except ValueError:
            pass

        try: ### maybe a bool
            ans[k] = config.getbool(section, k)
            continue
        except:
            pass

        ### just grab it as a string
        ans[k] = config.get(section, k)

    return ans

### utilities for setting up Monte-Carlo sampling objects

def path2generator(path):
    """a function that reads in an INI file and instantiates the relevant EventGenerator object
    """
    config = ConfigParser()
    config.read(path)

    samplingdistributions = subclasses(eventgen.SamplingDistribution)
    attributetransformations = subclasses(eventgen.AttributeTransformation)

    ### iterate through sections, instantiating a generator for each one
    generators = []
    timedistribution = None
    for name in config.get('general', 'generators').split(): ### iterate through generators
        gen = samplingdistributions[name](**parse(name, config))
        generators.append(gen)

        if isinstance(gen, eventgen.TimeDistribution):
            timedistribution = gen

    if timedistribution is None:
        timedistribution = eventgen.UniformEventTime(t0=0, dur=10)
        generators.append(timedistribution)

    ### add attribute transformations
    transforms = []
    if config.has_option('general', 'transforms'):
        for name in config.get('general', 'transforms').split(): ### iterate through transforms
            trans = attributetransformations[name]() ### FIXME: some of these take in arguments...
            transforms.append(trans)

    return eventgen.MonteCarloIntegrator(generators=generators, transforms=transforms), timedistribution

#-------------------------------------------------

### I/O and downsampling utilities

def glob_samples(rootdir, tag, start, stop, downsample=1):
    """find all samples within [start, stop) contained within rootdir with tag. return a subset, keeping one out of every `downsample` samples
    """
    # find the integer top-level directory labels for the relevant time range
    intstart = int(start)/utils.MOD_STRIDE
    intstop = int(stop)/utils.MOD_STRIDE

    segs = [[start, stop]]

    events = []
    prob = 1./downsample ### probability of keeping any individual sample
    ### iterate over top-level directories
    for directory in glob.glob('%s/*'%rootdir):
        inttime = int(os.path.basename(directory))
        if (intstart<=inttime) and (inttime<=intstop): ### there is some overlap
            ### iterate over sub-directories
            for subdirectory in glob.glob('%s/*-*'%directory):
                substart, subdur = [int(_) for _ in subdirectory.split('-')]
                intersection = utils.andsegments(segs, [[substart, substart+subdur]]) ### figure out whether this is relevant data
                if intersection: ### there is some overlap
                    path = utils.samples_path(rootdir, tag, substart, subdur) ### predict path
                    if os.path.exists(path):
                        events = eventgen.csv2events(path) ### load samples
                        rand = np.random.random(len(events))<=prob ### boolean array, true for events that we keep

                        ### throw away anything that is not included in the relevant window
                        substart, substop = intersection[0]
                        for event, keep in zip(events, rand):
                            if keep and (substart<=event.time) and (event.time<=substop):
                                events.append(event)

    ### return the final list
    return events
