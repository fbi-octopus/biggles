"""Simulate a Biggles dataset."""

from __future__ import print_function

import argparse
import json
import logging
import sys
import warnings
from os import path

import numpy as np
from _filedescriptor import filedescriptor
from indexing import indexed_data
import datetime
from indexing import Flag
from hdf5util import write_hdf5_sim
import h5py

logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')
log = logging.getLogger()


CHARS = np.array(list('abcdefghijklmnopqrstuvwxyz0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ'))

# Default options for simulation parameters
SIMULATION = {
    'options': {
        'time_span': [0, 128],
        'velocity_sigma': 0.01, # orignal: velocity state error - both, new: velocity state error - directionality
        'position_sigma': 0.25, # orignal: position state error - both, new: position state error - random
        'initial_velocity_sigma': 0.1,  # orignal: velocity state - random, new: speed state - directionality
        'crop': True,
        'demote': False,
        'constraint_region_radius': 1.0,
        'constraint_region_density': 0.0,
        'directionality_proportion': 0.0,
        'directionality_velocity': [0.5, 1.0], # original: velocity state - directionality, new: not used
        'directionality_position_sigma': 0.0, # original: not used, new: position state error - directionality
        'directionality_speed': 1.0, # original: not used, new: mean intial state speed - directionality
        'creation_events' : 'exponential', # affects 'rates'; valid values: 'exponantial', 'uniform'
        'zigzag' : 0.1,
        'split_probability' : 0.1,
    },

    'parameters': {
        'births_per_frame': 1.5,
        'survival_probability': 0.95,
        'observation_covariance': [
          [ 0.01, 0.00 ],
          [ 0.00, 0.01 ],
        ],
        'clutter_per_frame': 2.0,
        'observation_probability': 0.8,
    },

    # If not None, this is a Nx2 array of centres specifying constraining regions. Each row is the (x,y) position of the
    # region centre. The region radious is the constraint_region_radius option.
    'constraining_regions': None,
    'original' : False,
    'two_mode_tracks' : False,
}
def setup_parser(parser):
    parser.add_argument('--string', metavar='STRING', type=str,
            help='a string that will be included into the output file name: "*.simulation.<string>.input.json"')
    parser.add_argument('--config', '-f', metavar='FILE', type=str, nargs='?',
            help='Read options from a file (overridden by command-line options)')
    parser.add_argument('--no-crop', action='store_true', default=False,
            help='do not crop generated values to simulation domain')
    parser.add_argument('--demote', action='store_true', default=False,
            help='demote all observations from the tracks into the clutter')
    parser.add_argument('--original', action='store_true', default=False,
            help='Use a the original simulation approach')
    parser.add_argument('--to-file', action='store_true', default=False,
            help='Write output to file instead to "stdout"')
    parser.add_argument('--number', '-n', type=int, default = 1,
            help='number of ground truth states to generate with these settings')
    parser.add_argument('--realisations', '-r', type=int, default = 1,
            help='number of a observed realisations to simulate per ground truth state with these settings')
    parser.add_argument('--two-mode-tracks', action='store_true', default=False,
            help='all tracks will have a directional and a random part')
    parser.add_argument('output', metavar='FILE', type=str, nargs='?',
                        help = ('overwriting output filename'))
    parser.add_argument('--dir', type=str, default='.',
                        help = ('target folder (must exist)'))

    model_group = parser.add_argument_group('model parameters')
    model_group.add_argument('--birth-rate', metavar='FLOAT', type=float,
            help='mean new tracks being born per-100x100-pixel-per-frame')
    model_group.add_argument('--clutter-rate', metavar='FLOAT', type=float,
            help='mean spurious clutter observations per-100x100-pixel-per-frame')
    model_group.add_argument('--survival-probability', metavar='FLOAT', type=float,
            help='probability of a track surviving for one time step')
    model_group.add_argument('--mean-track-length', metavar='FLOAT', type=float,
            help='alternative to "survival probability"')
    model_group.add_argument('--observation-probability', metavar='FLOAT', type=float,
            help='probability of a track generating an observation in a frame')
    model_group.add_argument('--observation-error', metavar='FLOAT', type=float,
            help='standard deviation of the observation localisation error')

    sim_group = parser.add_argument_group('simulation parameters')
    sim_group.add_argument('--width', metavar='INT', type=int, default = 128,
            help='width of image in pixels')
    sim_group.add_argument('--height', metavar='INT', type=int, default = 128,
            help='height of image in pixels')
    sim_group.add_argument('--margin', metavar='INT', type=int, default = 10,
            help='Margin of image in pixels. The margin will be added to the image size, '
            'but not birth events will happen there. However clutter events will happen. '
            'This shall minimise tracks moving outside the FoV.')
    sim_group.add_argument('--duration', metavar='INT', type=int,
            help='duration of simulation (number of frames)')
    sim_group.add_argument('--start-time', metavar='INT', type=int,
            help='first time stamp in simulation')
    sim_group.add_argument('--position-sigma', metavar='FLOAT', type=float,
            help='standard deviation of position random walk (diffusion)')
    sim_group.add_argument('--velocity-sigma', metavar='FLOAT', type=float,
            help='standard deviation of velocity random walk')
    sim_group.add_argument('--initial-velocity-sigma', metavar='FLOAT', type=float,
            help='standard deviation of initial velocity')
    sim_group.add_argument('--seed', metavar='INT', type=int,
            help='seed to PRNG')
    sim_group.add_argument('--directionality-proportion', metavar='FLOAT', type=float,
            help='proportion of tracks forced to be directional')
    sim_group.add_argument('--directionality-position-sigma', metavar='FLOAT', type=float,
            help='position component of the system error for diffusion')
    sim_group.add_argument('--directionality-speed', metavar='FLOAT', type=float,
            help='length of the velocity vector for diffusion (random direction)')
    sim_group.add_argument('--directionality-velocity', metavar='X:Y', type=str,
            help='the velocity vector (uniform direction)')
    sim_group.add_argument('--zigzag', metavar='FLOAT', type=float,
            help='the probability that the track will change curse (random direction)')
    sim_group.add_argument('--split-probability', metavar='FLOAT', type=float,
            help='the probability that the track will split')
    sim_group.add_argument('--creation-events', metavar='DISTRIBUTION', type=str,
            default='uniform', choices=['exponential', 'uniform'],
            help='The probability distribution for birth/clutter rate.')

    constraint_group = parser.add_argument_group('constrained motion parameters')
    constraint_group.add_argument('--constraint-region-radius', metavar='FLOAT', type=float,
            help='radius of constraining regions in pixels')
    constraint_group.add_argument('--constraint-region-density', metavar='FLOAT', type=float,
            help='relative density of constraining regions in range [0, 1)')

def regions_containing(state):
    """
    Return an array specifying if the passed state is within any of the corresponding SIMULATION['constraining_regions'] regions. If
    there are no regions, return None.
    """
    if SIMULATION['constraining_regions'] is None:
        return None

    dxs = SIMULATION['constraining_regions'][:,0] - state[0]
    dys = SIMULATION['constraining_regions'][:,1] - state[2]
    radii = np.sqrt(dxs*dxs + dys*dys)
    return radii < SIMULATION['options']['constraint_region_radius']

def evolve_track(track, state, Q, num):
    A = np.array([[1, 1, 0, 0], [0, 1, 0, 0], [0, 0, 1, 1], [0, 0, 0, 1]], dtype=np.float)
    for _ in xrange(num):
        track.append(state.tolist())
        # Evolve state
        state = A.dot(state) + np.random.multivariate_normal(np.zeros_like(state), Q, 1)[0]
    return track, state

def generate_directed_motion(state):
    isotropic = not SIMULATION['original']
    vv = SIMULATION['options']['velocity_sigma'] * SIMULATION['options']['velocity_sigma']
    v_sigma = SIMULATION['options']['initial_velocity_sigma']
    if isotropic:
        phi = np.random.random()*2.0*np.pi
        r = SIMULATION['options']['directionality_speed'] + v_sigma * np.random.randn()
        state[1] = r * np.cos(phi)
        state[3] = r * np.sin(phi)
    else:
        state[1] = SIMULATION['options']['directionality_velocity'][0] + v_sigma * np.random.randn()
        state[3] = SIMULATION['options']['directionality_velocity'][1] + v_sigma * np.random.randn()
    pv = SIMULATION['options']['directionality_position_sigma'] * SIMULATION['options']['directionality_position_sigma']
    Q = np.diag([pv, vv, pv, vv])
    return state, Q

def generate_brownian_motion(state):
    state[1], state[3] = 0.0, 0.0
    pv = SIMULATION['options']['position_sigma'] * SIMULATION['options']['position_sigma']
    Q = np.diag([pv, 0.0, pv, 0.0])
    return state, Q

def sample_first_postion():
    w, h = SIMULATION['options']['image_size']
    margin = float(SIMULATION['options']['image_margin'])
    w0 = w - 2*margin
    h0 = h - 2*margin
    return w0 * np.random.random() + margin, h0 * np.random.random() + margin

def generate_2mode_track(start_time):
    ''' Generate a single molecule's track, with a part of random motion and a part of directed motion
    '''
    x0, y0 = sample_first_postion()
    state = np.array([ x0, 0.0, y0, 0.0 ])
    end_time = int(SIMULATION['options']['time_span'][1])

    track_length = int(np.random.exponential(scale= 1/(1-SIMULATION['parameters']['survival_probability'])))

    len_direct = int(track_length*SIMULATION['options']['directionality_proportion'])
    len_brown = track_length - len_direct

    direct_first = bool(np.random.randint(2))
    this_time = int(start_time)
    track = []
    if direct_first:
        state, Q = generate_directed_motion(state)
        track, state = evolve_track(track, state, Q, min(end_time, this_time + len_direct) - this_time)
        this_time = min(end_time, this_time + len_direct)
        state, Q = generate_brownian_motion(state)
        track, state = evolve_track(track, state, Q, min(end_time, this_time + len_brown) - this_time)
    else:
        state, Q = generate_brownian_motion(state)
        track, state = evolve_track(track, state, Q, min(end_time, this_time + len_brown) - this_time)
        this_time = min(end_time, this_time + len_brown)
        state, Q = generate_directed_motion(state)
        track, state = evolve_track(track, state, Q, min(end_time, this_time + len_direct) - this_time)
    try:
        assert(len(track) == min(end_time, start_time + track_length) - start_time)
    except:
        print('real length %d' % len(track))
        print('preceived length %d' % (min(end_time, start_time + track_length) - start_time))
        print('generated track length %d' % track_length)
        print('---')
        print('start time %d' % start_time)
        print('end_time %d' % end_time)
        print('---')
        if direct_first:
            print('len direct %d' % len_direct)
            print('---')
            print('this_time %d' % this_time)
            print('---')
            print('len brown %d' % len_brown)
        else:
            print('len brown %d' % len_brown)
            print('---')
            print('this_time %d' % this_time)
            print('---')
            print('len direct %d' % len_direct)
        print('---')
        raise
    return (start_time, track)

def zigzag_state(state, r):
    phi = np.arctan2(state[1], state[3]) + 2.0*np.arctan(np.random.randn()/2.0)
    state[1] = r * np.cos(phi)
    state[3] = r * np.sin(phi)

def extend_track(start_time, state, A, Q, r, zigzag):
    t = start_time
    track_states = []
    end_time = SIMULATION['options']['time_span'][1]
    p_survive = SIMULATION['parameters']['survival_probability']
    p_split = SIMULATION['options']['split_probability']
    while t < end_time:
        # Do we survive to this frame?
        if len(track_states) != 0 and np.random.random() >= p_survive:
            break

        # Save state
        track_states.append(state.tolist())

        # Evolve state
        new_state = A.dot(state) + np.random.multivariate_normal(np.zeros(4), Q, 1)[0]
        t = t + 1

        # Do we need to worry about constrained motion?
        if SIMULATION['constraining_regions'] is not None:
            # Determine which regions the previous and new state are in
            regions = regions_containing(state)
            new_regions = regions_containing(new_state)

            # Reject any motion which transitions from within state to outside of state i.e. those where there is at
            # least one region which previous is within and next is without.
            if np.any(np.logical_and(regions, np.logical_not(new_regions))):
                continue

        if zigzag and np.random.random() < SIMULATION['options']['zigzag']:
            zigzag_state(new_state, r)

        # Accept evolved state
        state = new_state

    branches = [(start_time, track_states)]
    if len(track_states) > 1 and p_split > np.random.random():
        rel_split_time = np.random.randint(len(track_states)-1)
        split_state = track_states[rel_split_time]
        split_time = start_time + rel_split_time
        new_state = A.dot(split_state) + np.random.multivariate_normal(np.zeros(4), Q, 1)[0]
        zigzag_state(new_state, r)
        branches.extend(extend_track(split_time+1, new_state, A, Q, r, zigzag))

    return branches

def generate_track(start_time):
    r"""
    Generate a single molecule's track.

    This uses the following motion model. At time t, the state. :math:`s_t` of the molecule is given by the vector

        [x, x', y, y']

    where (x,y) are the position and (x',y') are the instantaneous velocity. The state evolves according to

    .. math::
        s_{t+1} = A s_t + W_t,

        A = \left[
        \begin{array}{cccc}1& 1& 0& 0 \\
        0& 1& 0& 0 \\
        0& 0& 1& 1 \\
        0& 0& 0& 1 \end{array}\right]

    where :math:`W_t` is some Gaussian process with a known covariance matrix Q.
    This is hard-wired to be a diagonal matrix with
    diagonal elements set from the --{position,velocity}-sigma options.

    The track is returned as a pair giving start frame and a list of states,
    """

    A = np.array([[1, 1, 0, 0], [0, 1, 0, 0], [0, 0, 1, 1], [0, 0, 0, 1]], dtype=np.float)
    pv = SIMULATION['options']['position_sigma'] * SIMULATION['options']['position_sigma']
    vv = SIMULATION['options']['velocity_sigma'] * SIMULATION['options']['velocity_sigma']
    Q = np.diag([pv, vv, pv, vv])

    # Compute the initial state
    v_sigma = SIMULATION['options']['initial_velocity_sigma']
    x0, y0 = sample_first_postion()

    zigzag = False
    r = 0.0

    if SIMULATION['original']:

        state = np.array([ x0, v_sigma * np.random.randn(), y0, v_sigma * np.random.randn() ])

        # Should this track be directional?
        if np.random.random() < SIMULATION['options']['directionality_proportion']:
            state[1] = SIMULATION['options']['directionality_velocity'][0]
            state[3] = SIMULATION['options']['directionality_velocity'][1]

    else:
        # state: position -> random, velocity -> zero
        # system error: position -> position_sigma, velocity -> zero
        state = np.array([ x0, 0.0, y0, 0.0 ])
        Q[1, 1]=0.0
        Q[3, 3]=0.0

        # Should this track be directional?
        if np.random.random() < SIMULATION['options']['directionality_proportion']:
            # state: position -> random, velocity -> Gaussian(directionality_speed, initial_velocity_sigma)
            # system error: position -> directionality_postion_sigma, velocity -> velocity_sigma
            phi = np.random.random()*2.0*np.pi
            r = SIMULATION['options']['directionality_speed'] + v_sigma * np.random.randn()
            state[1] = r * np.cos(phi)
            state[3] = r * np.sin(phi)
            pv = SIMULATION['options']['directionality_position_sigma'] * SIMULATION['options']['directionality_position_sigma']
            Q = np.diag([pv, vv, pv, vv])
            zigzag = SIMULATION['options']['zigzag'] > 0

    # Evolve over time
    branches = extend_track(start_time, state, A, Q, r, zigzag)

    #print("Q={0}".format(Q))

    return branches

def expo(x, l):
    return l*np.exp(-x*l)

def cut_track(track, first_ts):
    ''' cutting off the excess from a track.
    Everything that comes before 'first_ts' will be removed. Returns the track or None,
    if everything was cut
    '''
    if track[0] >= first_ts:
        return track
    diff = first_ts - track[0]
    if diff >= len(track[1]):
        return None
    return [ first_ts, track[1][diff:]]

def generate_tracks():
    """
    Generate a ground truth set of tracks. Main function
    """
    start_time, end_time = SIMULATION['options']['time_span']
    w, h = SIMULATION['options']['image_size']
    creation_type = SIMULATION['options']['creation_events']
    tracks = []
    birth_rate = SIMULATION['parameters']['births_per_frame'] *(w*h/10000.0)
    log.info("birth rate {0} events/frame".format(birth_rate))
    if creation_type == 'uniform':
        ''' TODO we need to start before inital time '''
        ''' ln(4)*(mean track length) is 3rd quartille (Q3) of the track length distribution
            starting simulation at (start_time - Q3) to have a smooth transition
        '''
        mean_track_length = 1./(1.-SIMULATION['parameters']['survival_probability'])
        for frame in xrange(start_time - int(np.ceil(np.log(4) * mean_track_length)), end_time):
            n_tracks = np.random.poisson(birth_rate)
            for track_idx in xrange(n_tracks):
                if SIMULATION['two_mode_tracks']:
                    tracks.append(generate_2mode_track(frame))
                else:
                    tracks.extend(generate_track(frame))
        ''' cutting of the excess '''
        tracks = [ cut_track(t, start_time) for t in tracks ]
        tracks = [ t for t in tracks if not t is None ]
    elif creation_type == 'exponential':
        ''' choose a lambda such that the number of events at end_time
            is 5% of the number of events at start_time
        '''
        lambda_ = -np.log(0.05)/(end_time-start_time)
        total_events = (end_time-start_time +1) * birth_rate
        for frame in xrange(start_time, end_time):
            # How many tracks this frame?
            # Birth rate is in events per frame and 100x100 pixels
            n_tracks = np.random.poisson(total_events*expo(frame-start_time, lambda_))
            for track_idx in xrange(n_tracks):
                if SIMULATION['two_mode_tracks']:
                    tracks.append(generate_2mode_track(frame))
                else:
                    tracks.extend(generate_track(frame))
    else:
        log.error('unkown type of event creation "{0}"'.format(creation_type))

    return tracks

def observation_model(tracks):
    r"""
    Generate a set of observations from ground truth tracks.

    For each track, apply the observation model and return the set of observations.
    The observation, :math:`o_t`, of state :math:`s_t`  is given by:

    .. math::
        o_t = B s_t + V_t

        B = \left[ \begin{array}{cccc} 1& 0& 0& 0 \\
             0& 0& 1& 0 \end{array}\right]

    and :math:`V_t` is a Gaussian process with covariance R. R is a parameter to the simulation.
    """
    R = SIMULATION['parameters']['observation_covariance']
    p_observe = SIMULATION['parameters']['observation_probability']
    w, h = SIMULATION['options']['image_size']


    observed_tracks = []
    ground_truth = []
    unobserved = []
    for start_time, states in tracks:
        observations = []
        for state_idx, state in enumerate(states):
            # Skip observation with probability (1 - p_observe)
            if np.random.random() >= p_observe:
                continue
            observation = np.array([state[0], state[2]]) + np.random.multivariate_normal(
                np.zeros(2), np.eye(2)*R, 1)[0]

            # Clip observations to viewable area
            if SIMULATION['options']['crop']:
                if observation[0] < 0 or observation[0] >= w or observation[1] < 0 or observation[1] >= h:
                    continue

            observations.append((observation[0], observation[1], state_idx + start_time))

        if len(observations) > 1:
            observed_tracks.append({
                'observations': observations,
                'time_span': [observations[0][2], observations[-1][2]+1],
            })
            ground_truth.append([start_time, states])
        else:
            unobserved.append([start_time, states])

    return observed_tracks, ground_truth, unobserved

def generate_clutter():
    """
    Generate the false observations.
    """
    clutter = []

    w, h = SIMULATION['options']['image_size']
    start_time, end_time = SIMULATION['options']['time_span']
    clutter_rate = SIMULATION['parameters']['clutter_per_frame'] *(w*h/10000.0)
    log.info("clutter rate {0} events/frame".format(clutter_rate))
    creation_type = SIMULATION['options']['creation_events']
    lambda_ = -np.log(0.05)/(end_time-start_time)
    total_events = (end_time-start_time +1) * clutter_rate
    if creation_type == 'exponential':
        for frame_idx in xrange(start_time, end_time):
            n_clutter = np.random.poisson(total_events*expo(frame_idx-start_time, lambda_))
            for clutter_idx in xrange(n_clutter):
                clutter.append((
                    np.random.random() * w,
                    np.random.random() * h,
                    frame_idx
                ))
    elif creation_type == 'uniform':
        for frame_idx in xrange(start_time, end_time):
            n_clutter = np.random.poisson(clutter_rate)
            for clutter_idx in xrange(n_clutter):
                clutter.append((
                np.random.random() * w,
                np.random.random() * h,
                frame_idx
                ))
    else:
        log.error('unkown type of event creation "{0}"'.format(creation_type))

    return clutter

def random_str(num):
    return ''.join(CHARS[ map(int, np.random.random(num)*len(CHARS)) ])

def generate_constraining_regions():
    """
    Populate SIMULATION['constraining_regions'] with a field of constraining regions.
    """

    w, h = SIMULATION['options']['image_size']
    r = SIMULATION['options']['constraint_region_radius']
    region_density = SIMULATION['options']['constraint_region_density']

    domain_area = float(w*h)
    region_area = np.pi * r * r
    n_regions = int(np.ceil(region_density * domain_area / region_area))
    if n_regions == 0:
        return

    log.info('Generating {0} constraining regions'.format(n_regions))

    constraining_regions = np.random.random((n_regions, 2))
    constraining_regions[:,0] *= w
    constraining_regions[:,1] *= h
    SIMULATION['constraining_regions'] = constraining_regions

def adjust_birt_rate_wrt_margin():
    ''' | adjust the birth rate such that:
        | input BR is the normed BR that achieves a intended density in core FoV
        | adjusted BR is the BR that maintains the intended density in core FoV
        | but has an extended area. Hence the adjusted BR is lower
    '''
    br_intended = SIMULATION['parameters']['births_per_frame']
    w, h = SIMULATION['options']['image_size']
    margin = SIMULATION['options']['image_margin']
    core_fov = float((w - 2*margin) * (w - 2*margin))
    extended_fov = float(w*h)
    SIMULATION['parameters']['births_per_frame'] = br_intended *core_fov/extended_fov

def assign_simulation_parameters(args):
    # Assign model parameters
    if args.birth_rate is not None:
        SIMULATION['parameters']['births_per_frame'] = args.birth_rate
    if args.clutter_rate is not None:
        SIMULATION['parameters']['clutter_per_frame'] = args.clutter_rate
    if args.survival_probability is not None:
        SIMULATION['parameters']['survival_probability'] = args.survival_probability
    elif args.mean_track_length is not None:
        SIMULATION['parameters']['survival_probability'] = 1.0 - 1.0/args.mean_track_length
    if args.observation_probability is not None:
        SIMULATION['parameters']['observation_probability'] = args.observation_probability

    if args.observation_error is not None:
        oe = args.observation_error
        SIMULATION['parameters']['observation_covariance'] = [[oe*oe, 0], [0, oe*oe]]

def assign_simulation_options(args):
    # Assign simulation options
    if args.start_time is not None:
        duration = SIMULATION['options']['time_span'][1] - SIMULATION['options']['time_span'][0]
        SIMULATION['options']['time_span'] = [args.start_time, args.start_time + duration]
    if args.duration is not None:
        SIMULATION['options']['time_span'][1] = SIMULATION['options']['time_span'][0] + args.duration
    if args.position_sigma is not None:
        SIMULATION['options']['position_sigma'] = args.position_sigma
    if args.velocity_sigma is not None:
        SIMULATION['options']['velocity_sigma'] = args.velocity_sigma
    if args.initial_velocity_sigma is not None:
        SIMULATION['options']['initial_velocity_sigma'] = args.initial_velocity_sigma
    if args.no_crop is not None:
        SIMULATION['options']['crop'] = not args.no_crop
    if args.demote is not None:
        SIMULATION['options']['demote'] = args.demote
    if args.constraint_region_radius is not None:
        SIMULATION['options']['constraint_region_radius'] = args.constraint_region_radius
    if args.constraint_region_density is not None:
        SIMULATION['options']['constraint_region_density'] = args.constraint_region_density
    if args.directionality_proportion is not None:
        SIMULATION['options']['directionality_proportion'] = args.directionality_proportion
    if args.directionality_speed is not None:
        SIMULATION['options']['directionality_speed'] = args.directionality_speed
    if args.directionality_velocity is not None:
        SIMULATION['options']['directionality_velocity'] = map(float, args.directionality_velocity.split(":"))
        assert len(SIMULATION['options']['directionality_velocity']) == 2
    if args.directionality_position_sigma is not None:
        SIMULATION['options']['directionality_position_sigma'] = args.directionality_position_sigma
    if args.zigzag is not None:
        SIMULATION['options']['zigzag'] = args.zigzag
    if args.split_probability is not None:
        SIMULATION['options']['split_probability'] = args.split_probability
    if args.creation_events is not None:
        SIMULATION['options']['creation_events'] = args.creation_events
    #log.info(args.creation_events)
    #log.info(args.zigzag)
    assert( 0 <= SIMULATION['options']['zigzag'] <= 1 )
    SIMULATION['original'] = args.original
    # the flag is set and the proportion of dirctional tracks is neither 0 nor 1
    SIMULATION['two_mode_tracks'] = args.two_mode_tracks and (args.directionality_proportion % 1 != 0)

    SIMULATION['options']['image_size'] = [ args.width + 2*args.margin, args.height + 2*args.margin ]
    SIMULATION['options']['image_margin'] = args.margin

def parse_config_file(args):
    ''' parse the configuration file '''
    import ConfigParser
    import StringIO

    # A dirty hack to allow section-less config
    # http://stackoverflow.com/questions/2885190/using-pythons-configparser-to-read-a-file-without-section-name
    ini_str = '[root]\n' + open(args.config, 'r').read()
    ini_fp = StringIO.StringIO(ini_str)
    config = ConfigParser.RawConfigParser()
    config.readfp(ini_fp)

    # Copy configuration into arguments
    for key, value in config.items('root'):
        # replace '-' in key with '_' for checking args object
        key = key.replace('-', '_')

        if not hasattr(args, key):
            log.warning('Configuration has unknown parameter "{0}"'.format(key))

        ## don't replace arguments which have been specified on the command line
        #if key in args and getattr(args, key) is None:
        if True:
            try:
                setattr(args, key, int(value))
            except ValueError:
                try:
                    setattr(args, key, float(value))
                except ValueError:
                    setattr(args, key, value)

def write_output(output, args):
    ''' write a single simulated data set '''
    write_hdf5 = False
    if args.output is not None:
        outfname=path.join(args.dir, args.output)
        if write_hdf5:
            with open(outfname, 'w') as fh:
                write_hdf5_sim(output, fh)
        else:
            with open(outfname, 'w') as fh:
                json.dump(output, fh)
    elif args.to_file:
        ext = ""
        if args.string is not None:
            ext += ".{0}".format(args.string)
        if args.demote:
            ext += ".input"
        froot = path.join(args.dir, output['metadata']['project'])
        if write_hdf5:
            ext += ".hdf5"
            with h5py.File(froot + ext, 'w') as fh:
                write_hdf5_sim(output, fh)
        else:
            ext += ".json"
            with open(froot + ext, 'w') as fh:
                json.dump(output, fh)
        outfname = froot + ext
    else:
        output_file = sys.stdout
        outfname = ''
        if write_hdf5:
            write_hdf5_sim(output, output_file)
        else:
            json.dump(output, output_file)

    print(outfname)

def generate_realisation(args, ground_truth_tracks, gt_str):
    ''' generate a single realisation, given the ground truth and its identfier '''
    # Generate clutter
    clutter = generate_clutter()

    # Observe tracks
    observed_tracks, ground_truth_tracks, unobserved = observation_model(ground_truth_tracks)

    observations, observed_tracks, clutter = indexed_data(observed_tracks, clutter)

    root_tracks = observed_tracks[:]
    root_clutter = clutter[:]

    # Do we demote all the observed tracks?
    if args.demote:
        for track in root_tracks:
            for obs in track['observations']:
                root_clutter.append(obs)
        root_tracks = []

    # Sort the clutter by ascending time stamp
    #clutter = sorted(clutter, key=lambda obs: obs[2])

    # HACK: if constraining_regions is non-None, we need to convert it
    if SIMULATION['constraining_regions'] is not None:
        SIMULATION['constraining_regions'] = SIMULATION['constraining_regions'].tolist()

    file_des = filedescriptor()
    #rstr = ''.join(np.random.choice(np.array(list(chars)), 12))
    #proj = datetime.datetime.utcnow().strftime("%Y%m%d_00%H_%M%S_") + file_des["uuid4"] + ".simulation"
    rstr = gt_str + random_str(6)
    proj = datetime.datetime.utcnow().strftime("%Y%m%d_00%H_%M%S_") + rstr + ".simulation"


    output = {
        'flags' : [Flag.INDEXED],
        'clutter': root_clutter,
        'tracks': root_tracks,
        'ground_truth': ground_truth_tracks,
        'unobserved': unobserved,
        'observations': observations,
        'simulation': SIMULATION,
        'metadata': {'project': proj, 'channel': 0},
        'fileheader': file_des,
        'initial_partition' : { 'tracks' : observed_tracks, 'clutter': clutter}
    }

    write_output(output, args)

def simulate_data_set(args):
    ''' simulate a single data set - potentially with several realisations'''
    # Generate constraining regions
    if args.constraint_region_density > 0:
        generate_constraining_regions()

    # Generate ground truth
    ground_truth_tracks = generate_tracks()
    gt_str = random_str(6)

    for i in range(args.realisations):
        generate_realisation(args, ground_truth_tracks, gt_str)

def main(args):
    warnings.filterwarnings('once')
    # Parse config file if present
    if not path.isdir(args.dir):
        print('"%s" is not a directory' % args.dir)
    if args.config is not None:
        parse_config_file(args)
    assign_simulation_parameters(args)
    assign_simulation_options(args)
    adjust_birt_rate_wrt_margin()

    # Set random seed
    np.random.seed(args.seed)

    for i in range(args.number):
        simulate_data_set(args)

    return 0

if __name__ == '__main__':
    # Parse the command line options
    parser = argparse.ArgumentParser(description=globals()['__doc__'],
                                     formatter_class = argparse.ArgumentDefaultsHelpFormatter)
    setup_parser(parser)
    args = parser.parse_args()

    # Run the main  program
    sys.exit(main(args))

# vim:tw=120
