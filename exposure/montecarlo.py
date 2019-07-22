"""a module that implements Monte-Carlo logic and supporting utilties when defining a population
"""
__author__ = "Reed Essick (reed.essick@gmail.com)"

#-------------------------------------------------

import numpy as np

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

#-------------------------------------------------

### utilities for setting up Monte-Carlo sampling objects

def path2generator(path):
    """a function that reads in an INI file and instantiates the relevant EventGenerator object
    """
    config = ConfigParser()
    config.read(path)

    ### iterate through sections, instantiating a generator for each one
    generator = eventgen.EventGenerator()

    ### keep track of which generators we instantiate
    ### used to make sure we have everything we need and that we don't repeat anything
    gentypes = set()

    # figure out what is declared
    types = eventgen.SamplingDistribution.__subclasses__() ### get the different types of distributions
    samplingdistributions = subclasses(eventgen.SamplingDistribution)       ### get instantiators for all distributions
    attributetransformations = subclasses(eventgen.AttributeTransformation)

    # add generators
    timedistribution = None
    for name in config.get('general', 'generators').split(): ### iterate through generators
        for parent in samplingdistributions[name].__bases__:
            if parent in gentypes: ### check to make sure we haven't repeated anything
                raise RuntimeError('cannot specify more than one instance of %s!'%samplingdistributions[name].__name__)
            gentypes.add(parent)

        gen = samplingdistributions[name](**parse(name, config))
        generator.append_generator(gen)

        if isinstance(gen, eventgen.TimeDistribution):
            timedistribution = gen

    if timedistribution is None:
        timedistribution = eventgen.UniformEventTime(t0=0, dur=10)
        generator.append_generator(timedistribution)
        for parent in eventgen.UniformEventTime.__bases__:
            gentypes.add(parent)

    ### figure out whether we have enough different types of generators
    ### includes declaring required transformations
    # check for things that *must* be specified
    for gentype in [eventgen.TimeDistribution, eventgen.MassDistribution, eventgen.SpinDistribution, eventgen.OrientationDistribution, eventgen.EccentricityDistribution]:
        if gentype not in gentypes:
            raise RuntimeError('must specify an instance of %s!'%gentype.__name__)

    # check for things that may be specified in several ways
    required_transforms = [eventgen.SourceMass2DetectorMass.__name__] ### always require this transformation
    if eventgen.DistanceDistribution in gentypes:
        if eventgen.RedshiftDistribution in gentypes:
            raise RuntimeError('cannot specify both DistanceDistribution and RedshiftDistribution!')
        required_transforms.append(eventgen.LuminosityDistance2Redshift.__name__) ### must include transformation from distance to redshift
    elif eventgen.RedshiftDistribution in gentypes:
        required_transforms.append(eventgen.Redshift2LuminosityDistance.__name__) ### must include transformation from redshift to distance
    else:
        raise RuntimeError('must specify either DistanceDistribution or RedshiftDistribution!')

    # add attribute transformations
    if config.has_option('general', 'transforms'):
        for name in config.get('general', 'transforms').split(): ### iterate through transforms
            if name in required_transforms:
                required_transforms.remove(name)
            trans = attributetransformations[name]() ### FIXME: some of these take in arguments...
            generator.append_transform(trans)

    if required_transforms:
        raise RuntimeError('required AttributeTransformation (%s) not specified!'%(', '.join(required_transforms)))

    return generator, timedistribution

#-------------------------------------------------

### Monte-Carlo logic

def update_max_redshift(generator, network):
    """update the maximum redshift or distance defined in generator to avoid sampling too distant sources

    **NOTE**: right now, we just skip this and leave the generator untouched...this will not scale well if the detector sensitivities change significantly!
"""
    pass
    '''
find maximum distance for detectable things
    -> difficult when we use a non-central chi2 distrib to measure p(det)...
    -> maybe just require cdf(det) to be less than some small fraction --> control the systematic error from this truncation to an acceptable level?

Maybe we just Monte-Carlo at a fixed distance and determine the maximum SNR observed? With that, we can scale the maximum distance based on that fixed distance and the assumption that everything with SNR below some value has negligible probability of being detected
'''

def update_montecarlo_counters(generator, event, e1, e2, de1, de2, Nparams):
    """update the montecarlo counters for our smart termination condition
    """
    # extract basic params
    pdet = event.pdet

    # compute the effective pdf and jacobian
    pdf = 1.
    jac = np.empty(Nparams, dtype=float)
    i = 0
    for gen in generator._generators:
        variates = [getattr(event, lbl) for lbl in gen._variates]

        pdf *= gen.pdf(*variates)
        jac[i:i+len(gen._params)] = gen.jacobian(*variates)

    # update the counters
    e1 += pdet
    e2 += pdet**2
    de1 += pdet*jac/pdf
    de2 += 2*pdet**2*jac/pdf

    return e1, e2, de1, de2

def montecarlo(generator, network, min_num_samples=DEFAULT_MIN_NUM_SAMPLES, error=DEFAULT_ERROR):
    """generate samples from generator until we reach the termination conditions specified by min_num_samples and error
    """
    ### generate the minimum number of samples required
    for event in generator.generate_events(n_itr=np.infty, n_event=min_num_samples):
        pass

    # figure out whether we want to even attempt the "smart termination condition"
    if error!=np.infty:

        ### set up the further iteration
        Nparams = sum(len(gen._params) for gen in generator._generators)
        e1 = 0
        e2 = 0
        de1 = np.zeros(Nparams, dtype=float)
        de2 = np.zeros(Nparams, dtype=float)
        for event in generator._events:
            e1, e2, de1, de2 = update_montecarlo_counters(generator, event, e1, e2, de1, de2, Nparams)

        for event in generator.generate_events(n_itr=np.infty, n_event=np.infty, generator=True): ### go forever until we are told to stop
            if np.all(2*e1**2 * np.abs(de1) * error >= 3*np.abs(2*e2*de1 - e1*de2)):
                break ### all smart termination conditions are satisfied, so we break the broader iteration

            ### update monte-carlo sums
            e1, e2, de1, de2 = update_montecarlo_counters(generator, event, e1, e2, de1, de2, Nparams)

    ### format the samples into a numpy structured array compatible with writing to disk
    attrs = [(attr, 'S50' if attr=='approximant' else 'float') for attr in sorted(generator.attributes)+sorted(eventgen.Event._waveattrs)]
    return np.array([tuple(getattr(event, attr) for attr, _ in attrs) for event in generator._events], dtype=attrs)
