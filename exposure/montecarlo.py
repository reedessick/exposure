"""a module that implements Monte-Carlo logic and supporting utilties when defining a population
"""
__author__ = "Reed Essick (reed.essick@gmail.com)"

#-------------------------------------------------

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
    for k in config.options():
        try: ### maybe it's an integer
            ans[k] = config.getint(section, k)
            continue
        except ValueError:
            pass

        try: ### maybe a float
            ans[k] = config.getint(section, k)
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
    for name in config.get('general', 'generators').split(): ### iterate through generators
        for parent in samplingdistributions[name].__bases__:
            if parent in gentypes: ### check to make sure we haven't repeated anything
                raise RuntimeError('cannot specify more than one instance of %s!'%samplingdistributions[name].__name__)
            gentypes.add(parent)

        gen = samplingdistributions[name](**parse(name, config))
        generator.append_generator(gen)

    ### figure out whether we have enough different types of generators
    ### includes declaring required transformations
    # check for things that *must* be specified
    for gentype in [eventgen.TimeDistribution, eventgen.MassDistribution, eventgen.SpinDistribution, eventgen.OrientationDistribution]:
        if gentype not in gentypes:
            raise RuntimeError('must specify an instance of %s!'%gentype.__name__)

    # check for things that may be specified in several ways
    required_transforms = []
    if eventgen.DistanceDistribution in gentypes:
        if eventgen.RedshiftDistribution in gentypes:
            raise RuntimeError('cannot specify both DistanceDistribution and RedshiftDistribution!')
    elif eventgen.RedshiftDistribution in gentypes:
        required_transforms.append(eventgen.Redshift2LuminosityDistance) ### must include transformation from redshift to distance
    else:
        raise RuntimeError('must specify either DistanceDistribution or RedshiftDistribution!')

    # add attribute transformations
    for name in config.get('general', 'transforms').split(): ### iterate through transforms
        if name in required_transforms:
            required_transforms.remove(name)
        trans = attributetransformations[name]() ### FIXME: some of these take in arguments...
        generator.append_transform(trans)

    if required_transforms:
        raise RuntimeError('required AttributeTransformation (%s) not specified!'%(', '.join(required_transforms)))

    return generator

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

def update_montecarlo_counters(event, e1, e2, de1, de2, Nparams):
    """update the montecarlo counters for our smart termination condition
    """
    raise NotImplementedError('return e1, e2, de1, de2')

def montecarlo(generator, network, min_num_samples=DEFAULT_MIN_NUM_SAMPLES, error=DEFAULT_ERROR):
    """generate samples from generator until we reach the termination conditions specified by min_num_samples and error
    """
    ### generate the minimum number of samples required
    generator.generate_events(n_iter=np.infty, n_event=min_num_samples, generator=False) ### just store these quietly internally

    # figure out whether we want to even attempt the "smart termination condition"
    if error!=np.infty:

        ### set up the further iteration
        Nparams = sum(len(gen._params) for gen in generator._generators)
        e1 = 0
        e2 = 0
        de1 = np.zeros(Nparams, dtype=float)
        de2 = np.zeros(Nparams, dtype=float)
        for event in generator._events:
            e1, e2, de1, de2 = update_montecarlo_counters(event, e1, e2, de1, de2, Nparams)

        for event in generator.generate_events(n_iter=np.infty, n_event=np.infty, generator=True): ### go forever until we are told to stop
            for i in range(Nparams):
                raise NotImplementedError('check smart termination condition. If this is not satisfied, we break')
            else:
                break ### all smart termination conditions are satisfied, so we break the broader iteration

            ### update monte-carlo sums
            e1, e2, de1, de2 = update_montecarlo_counters(event, e1, e2, de1, de2, Nparams)

    ### format the samples into a numpy structured array compatible with writing to disk
    attrs = sorted(generator.attributes)
    return np.array([[getattr(event, attr) for attr in attrs] for event in generator._events], dtype=[(attr, float) for attr in attrs])