"""Track partitions from JSON input"""

from __future__ import print_function

import argparse
import json
import logging
import warnings
import os
import sys
#import subprocess
#import tempfile
from _utilities import count2str, float2str, sec2str, fifo, UnsafeOpenError, safe_open
from _utilities import four2str, number_parser, get_random_string, auto_corr_time
from _util_output import print_final_state
from _move_stats import get_acceptance_rate
from _sample_recorder import SampleRecorder, get_output_filename, get_tracking_dir
from _track_utils import init_parameters, convert_initial_partition
from _track_utils import Tracker, TickTimer, get_max_ged
from convergence import ControlData
from indexing import Flag, has_flag, indexed_partition
from hdf5util import write_hdf5_observations, write_hdf5_indexed_partition, hdf5tracks3
from random_partition import max_partition, rev_partition

import time
import datetime
import h5py
import signal
import traceback

logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')
logFormatter = logging.Formatter('%(levelname)s: %(message)s')
log = logging.getLogger()

#consoleHandler = logging.StreamHandler()
#consoleHandler.setFormatter(logFormatter)
#log.addHandler(consoleHandler)

from numpy import mean
from biggles import Stepper
from biggles import IndexedPartition, log_pdf, sample_model_parameters_given_partition

BIGGLESVERSION = '013'

class SignalReceived(Exception):
    def __init__(self, signum, stack):
        self.signum = signum
        self.stack = stack
    def __str__(self):
        return repr(self.signum)

def signal_handler(signum, stack):
    raise SignalReceived(signum, stack)

def setup_parser(parser):
    parser.add_argument('input', metavar='FILE', type=str, nargs='?',
            help='Read partition from FILE. If omitted read from standard input.')
    parser.add_argument('--tick-interval', '-i', metavar='INTERVAL', type=float, default=0.1,
            help='Request tracker state every INTERVAL seconds.')
    parser.add_argument('--unsafe', action='store_true', default=False,
            help='Overwrite existing files. If omitted do not overwrite existing files.')
    parser.add_argument('--norecord', action='store_true', default=False,
            help='''Don't write anything to file''')
    parser.add_argument('--sim-pdf', action='store_true', default=False,
            help='''estimate the best posterior of the input simulations''')
    parser.add_argument('--max', action='store_true', default=False,
            help='''Start with a max track partition (clutter will be minimised before starting)''')
    parser.add_argument('--report-act', action='store_true', default=False,
            help='''report the autocorrelation times''')
    parser.add_argument('--comment', type = str,
            help = '''Have you anything to say?''')
    parser.add_argument('--dual', action='store_true', default=False,
            help='''Start with two inital partitions.''')
    parser.add_argument('-a', "--obs-cov-lag", metavar='LAG', type=int, default = 20,
                        help = "Sample the observation error covariance every LAG samples")
    parser.add_argument('-q', '--process-noise', type = float, nargs = 2,
                        default = [0.1, 0.01], metavar = "Q[i, i]",
                        help = "Kalman filter process noise Q. Q[0, 0] = Q[2, 2] localisation noise, "
                        "Q[1, 1] = Q[3, 3] velocity noise. Q[i, j] = 0 iff i != j. "
                        "Values are the square root of multiples of the light speed."
                        )
    parser.add_argument("--min-samples", type = int, default = 100000,
            help = '''Auto burn-in ONLY. The minimum number of samples to be taken.
            The burn-in convergence
            test ist not optimal yet. This value prevents a too short burn-in period
            for a small number (ca. < 500) of observations.''')
    parser.add_argument("--start", type=str, nargs = 2, default = ["mini", "maxi"],
        choices = ["input", "mini", "maxi", "rev", "split"],
        help = """Initialisation of the two chains for the "--dual" option.
            "input" - the input partition;
            "mini" - the clutter-only partition;
            "maxi" - a random partition that tries to assign as many as possible
                      observations to tracks;
            "rev" - as "maxi" but starting from the last time stamp;
            "split" - uses half the tracks of "maxi". If both values are "split"
            then the second partition will used the other half from the first.""" )

    logging_group = parser.add_argument_group('logging', 'Logging to disk.')
    logging_group.add_argument('--log', '-l', metavar='FILE', type=str, nargs='?', const='track_log.txt', default=None,
            help='Write a copy of the log to FILE. If FILE is omitted, use "track_log.txt"')
    logging_group.add_argument('--log-states', '-s', metavar='FILE', type=str, nargs='?', const='states.txt', default=None,
            help='Append JSON formatted states, one per line, to FILE on each tick. IF FILE is omitted, use "states.txt"')
    logging_group.add_argument('-o', '--output', metavar='FILE', type=str,
            help='''Write output to FILE.
            This will be the new standard output
            ''')
    logging_group.add_argument('-r', '--records', metavar='NUMBER', type=int, default = 2000,
                               help = 'number of samples to record')
    logging_group.add_argument('--run-time', action='store_true', default=False,
                               help = 'record the runtime')
    logging_group.add_argument('--log-best', action='store_true', default=False,
                               help = 'write the best partition')
    logging_group.add_argument('--log-ged', action='store_true', default=False,
                               help = 'write the ged-time response')
    logging_group.add_argument('--log-all', action='store_true', default=False,
                               help = 'write many scores to file')

    terminate_group = parser.add_argument_group('termination',
        'Control termination of tracking. If no options are specified then tracking will continue indefinitely.'
    )
    terminate_group.add_argument('--duration', '-t', metavar='SECONDS', type=int, default=None,
            help='Stop tracking after DURATION seconds.')
    terminate_group.add_argument('--samples', '-n', metavar='NUMBER', type=str, default=None,
            help='Stop tracking after NUMBER samples have been drawn.')
    terminate_group.add_argument('--mean-pdf', '-p', metavar='NUMBER', type=float, default=None,
            help='Stop tracking after the mean log PDF is larger than -NUMBER.')
    terminate_group.add_argument('--convergence', '-c', type=str,
                                 choices=['ged', 'anova', 'gelman', 'ano_ged', 'gel_ged'],
                                 default = 'gel_ged',
                                 help = "method to determine convergence")

    terminate_group.add_argument('--fifo-length', type=int, default=15, metavar = 'FACTOR'
        'length of the chains for convergence test is FACTOR x1000')


def print_r(parameters, printer):
    printer("R={0}".format(parameters.observation_error_covariance))

def input_from_hdf5(fh):
    '''reads input data from HDF5 file'''
    return {
        'observations': fh['observations'][:].tolist(),
        'clutter' : fh['clutter'][:].tolist(),
        'tracks' : hdf5tracks3(fh['tracks']),
        'flags' : [Flag.INDEXED]}

def print_record_act(recorder, logger):
    br_act = auto_corr_time(recorder.records['param_samples']['births_per_frame'])
    cr_act = auto_corr_time(recorder.records['param_samples']['clutter_per_frame'])
    op_act = auto_corr_time(recorder.records['param_samples']['observation_probability'])
    sp_act = auto_corr_time(recorder.records['param_samples']['survival_probability'])
    for key, val in [ ('BR', br_act),  ('CR', cr_act), ('OP', op_act), ('SP', sp_act) ]:
        logger("Auto correlation time %s: %.2f" % (key, val))

def print_header(total_obs, state, have_sims = False):
    ged_hstr = " GED"
    if have_sims:
        ged_hstr = "d-GT"
    headstr = "{0:>12s} {1:>7s} {2:>8s} {3:>6s} {4:>6s} {5:>6s} {6:>8s} {7:>7s} {8:>4s} {9:>7s} {10:4s}".format(
        'Time', 'Samples', 'best PDF', 'Reject', 'Identy', 'Accept',
        'Rate', 'Clutter', '#Tr', '+/-PDF', ged_hstr
    )
    params = json.loads(state.best_model_parameters.to_json())
    log.info('-' * len(headstr))
    volume = state.best_partition.volume()
    fact = 10000.0/volume
    out_list1 = ['observations: {0}'.format(total_obs)]
    out_list2 = []
    out_list2.append('b = {0:.2f}'.format(fact*params['births_per_frame']))
    out_list2.append('c = {0:.2f}'.format(fact*params['clutter_per_frame']))
    out_list2.append('o = {0:.2f}'.format(params['observation_probability']))
    out_list2.append('s = {0:.2f} ({1:.2f})'.format(params['survival_probability'],
                    1.0/(1.0-params['survival_probability'])))
    oc = params['observation_covariance']
    out_list1.append('R = [[{0:.3f}, {1:.3f}], [{2:.3f}, {3:.3f}]]'.format(
        oc[0][0], oc[0][1], oc[1][0], oc[1][1]))
    out_list1.append(' %.1f pix, %d fr.' % (volume, state.best_partition.duration))
    acc_rate = get_acceptance_rate(state.moves_accepted, state.moves_identity, state.moves_rejected)
    ar_list = []
    for sign in sorted(acc_rate.keys()):
        val = acc_rate[sign]
        if val is None:
            ar_list.append("{0}  --- ".format(sign[0]))
        else:
            ar_list.append("{0}={1:4.1f}%".format(sign[0], 100.0*val))
    log.info(', '.join(out_list1))
    log.info(', '.join(out_list2))
    log.info(', '.join(ar_list))
    log.info('-' * len(headstr))
    log.info(headstr)
    log.info('-' * len(headstr))

def write_partition_indices(partition, pool, args, fh):
    '''writes the indices of the partition to file'''
    reference_file = os.path.realpath(args.input)
    index_part = IndexedPartition()
    index_part.from_partition(partition)
    part_dict = {'reference' : reference_file}
    part_dict['tracks'] = index_part.tracks()
    part_dict['clutter'] = index_part.clutter()
    part_dict['observations'] = [ (o.x, o.y, int(o.t)) for o in pool]
    part_dict['flags'] = [Flag.INDEXED]
    call_paras = dict(vars(args))
    del(call_paras['func'])
    part_dict['call parameters'] = call_paras
    json.dump(part_dict, fh)

def write_input_partition(partition, args, fh):
    '''writes a partition  suitable for a input file'''
    part_dict  = json.loads(partition.to_json())
    call_paras = dict(vars(args))
    del(call_paras['func'])
    part_dict['call parameters'] = call_paras
    json.dump(part_dict, fh)

def make_maxi_partition(observations, pool):
    maxi_tracks, maxi_clutter = max_partition(observations)
    maxi_partition = IndexedPartition()
    maxi_partition.set_tracks(maxi_tracks)
    maxi_partition.set_clutter(maxi_clutter)
    return maxi_partition.to_partition(pool)

def make_split_partition(observations, pool):
    maxi_tracks, maxi_clutter = max_partition(observations)
    i = 0
    s_tracks = {0 : [], 1 : []}
    s_clutter = {0 : list(maxi_clutter), 1 : list(maxi_clutter)}
    for track in maxi_tracks:
        s_tracks[i].append(track)
        i = (i + 1) % 2
        s_clutter[i].extend(track['observations'])
    half1 = IndexedPartition()
    half1.set_tracks(s_tracks[0])
    half1.set_clutter(s_clutter[0])
    half2 = IndexedPartition()
    half2.set_tracks(s_tracks[1])
    half2.set_clutter(s_clutter[1])
    return half1.to_partition(pool), half2.to_partition(pool)

def make_rev_partition(observations, pool):
    r_tracks, r_clutter = rev_partition(observations)
    r_partition = IndexedPartition()
    r_partition.set_tracks(r_tracks)
    r_partition.set_clutter(r_clutter)
    return r_partition.to_partition(pool)

def make_mini_partition(observations, pool):
    mini_partition = IndexedPartition()
    mini_partition.set_tracks([])
    mini_partition.set_clutter(range(len(observations)))
    return mini_partition.to_partition(pool)

def simulation_posterior(initial_partition, args):

    if has_flag(initial_partition, Flag.INDEXED):
        pass
    else:
        initial_partition = indexed_partition(initial_partition)

    parameters = init_parameters()
    parameters.initQ(args.process_noise*2)

    have_sims, input_partition, sims_partition, pool = convert_initial_partition(initial_partition)
    if not have_sims:
        log.info("NO simulation found.")
        return

    posts = []
    for i in range(10000):
        sample_model_parameters_given_partition(sims_partition, parameters)
        posts.append(log_pdf(sims_partition, parameters))
    log.info("Best PDF = %f" % max(posts))

def track_dual(initial_partition, args):
    '''does tracking with dual burn-in'''

    out_format = 'json'

    if has_flag(initial_partition, Flag.INDEXED):
        pass
    else:
        initial_partition = indexed_partition(initial_partition)

    if not args.norecord:
        directory  = get_tracking_dir(args.input)
        output_file = open(
            get_output_filename(args.input, directory,
                                "tracked." + initial_partition['rand_id']), 'w')
        recorder = SampleRecorder(initial_partition)

    initial_model = init_parameters()
    initial_model.initQ(args.process_noise*2)

    have_sims, input_partition, sims_partition, pool = convert_initial_partition(initial_partition)

    dancers = []
    start_partitions = []
    split_part = None
    singletons = make_mini_partition(initial_partition['observations'], pool)
    for part_type in args.start:
        if part_type == "input":
           dancers.append(Tracker(Stepper(initial_model, input_partition)))
           start_partitions.append(input_partition)
        elif part_type == "maxi":
            start_partition = make_maxi_partition(initial_partition['observations'], pool)
            dancers.append(Tracker(Stepper(initial_model, start_partition)))
            start_partitions.append(start_partition)
        elif part_type == "rev":
            start_partition = make_rev_partition(initial_partition['observations'], pool)
            dancers.append(Tracker(Stepper(initial_model, start_partition)))
            start_partitions.append(start_partition)
        elif part_type == "mini":
            dancers.append(Tracker(Stepper(initial_model, singletons)))
            start_partitions.append(singletons)
        elif part_type == "split":
            if split_part is None:
                half1, split_part = make_split_partition(initial_partition['observations'], pool)
                dancers.append(Tracker(Stepper(initial_model, half1)))
                start_partitions.append(half1)
            else:
                dancers.append(Tracker(Stepper(initial_model, split_part)))
                start_partitions.append(split_part)

    start_time = time.time()
    min_sample_count = args.min_samples

    for dancer in dancers:
        dancer.sample_observation_error(args.obs_cov_lag)

    log.info('convergence method = %s' % (args.convergence,))
    ctrl = ControlData(start_time, len(pool), args.convergence, args.fifo_length)
    ctrl.print_explaination(log.info)
    ctrl.print_header_explain(log.info)
    ctrl.set_logger(log.info)
    ctrl.report_ac_time(args.report_act)
    ctrl.set_start_partitions(start_partitions)
    if args.log_all:
        fname = get_output_filename(args.input, '.', 'all.' + initial_partition['rand_id'])
        fname = fname.replace('.json', '.csv')
        ctrl.write_all_log(fname)
    max_ged = get_max_ged(initial_partition['observations'])
    ctrl.set_max_ged(max_ged)
    if have_sims:
        log.info("simulation found.")
        ctrl.set_simulation(sims_partition)
    timer = TickTimer(lag = 2)

    loop_size = 25

    do_record = False
    init_ged = dancers[0].state().current_partition.ged(dancers[1].state().current_partition)

    log.info("Max GED = {1}, Initial distance = {0:5.1f}%".format(
        100.0*init_ged/float(max_ged), max_ged))

    try:
        for i in range(ctrl.buffer_size()):
            ctrl.states = [dancer.run(1) for dancer in dancers]
            ctrl.set_sampling_rates([dancer.sample_rate for dancer in dancers])

        while True:
            ctrl.states = [dancer.run(loop_size) for dancer in dancers]
            ctrl.set_sampling_rates([dancer.sample_rate for dancer in dancers])
            tick_time = time.time() - start_time
            if timer.passed(tick_time):
                ctrl.try_output(tick_time)

                if ctrl.converged() and ctrl.sample_count >= min_sample_count:
                    end_time = time.time() - start_time
                    log.info('Tracker has converged after %s' % sec2str(end_time))
                    do_record = not args.norecord
                    if args.run_time:
                        fh = open(get_output_filename(args.input, '.', 'time', insert_time = False), 'a')
                        fh.write("{0:.1f}\n".format(end_time))
                    break

            if args.duration and args.duration < tick_time:
                log.info('Stopping since duration has exceeded %s seconds.' % (args.duration,))
                break

            if args.samples and number_parser(args.samples, int) <= ctrl.sample_count:
                log.info('Stopping since sample count has exceeded %s.' % (args.samples,))
                break
        if not do_record:
            return 1
        dancer = dancers[0]
        if (dancers[1].state().best_log_pdf > dancers[0].state().best_log_pdf):
            dancer = dancers[1]
        recorder.set_burnin_length(dancer.state().sample_count)
        dancer.sample_observation_error(1)
        # this is not necessary
        # dancer.run(5000) # do some sampling with new 'sample_observation_error' setting
        num_partiton_links = list()
        fact = 8.0
        num_steps = max(1, int(fact*0.25/dancer.state().recent_acceptance_rate))
        log.info("Recent acceptance rate = %.2f%%" % (100.0*dancer.state().recent_acceptance_rate))
        log.info("Recording every %d sample (factor = %.2f)" % (num_steps, fact))
        accepted = 0.
        last_partition = dancer.state().current_partition
        for i in range(args.records):
            recorder.record_state(dancer.run(num_steps))
            if last_partition.ged(dancer.state().current_partition) > 0:
                accepted += 1.
            num_partiton_links.append(singletons.ged(dancer.state().current_partition))
            last_partition = dancer.state().current_partition
            #if have_sims:
            #    pass
        recorder.write_results(args, output_file)
        log.info("ACT number of links: %.2f" % auto_corr_time(num_partiton_links))
        print_record_act(recorder, log.info)
        log.info("Partition change rate during recording: %.2f%%" % (100.*accepted/args.records))
    except KeyboardInterrupt:
        log.warning('Terminating tracking due to keyboard interrupt.')
        return 1
    except SignalReceived as sig:
        for line in traceback.format_stack():
            log.warning(line.strip())
        log.warning('Terminating tracking due to user signal {0}'.format(sig))
        return 1
    #except Exception as ex:
        #for line in traceback.format_stack():
        #    log.warning(line.strip())
        #log.warning('Terminating tracking due to exception "{0}"'.format(ex))
        #raise ex
    finally:
        ctrl.print_criterion()
        for chain, dancer in enumerate(dancers):
            log.info("Chain-%d:" % chain)
            print_final_state(log.info, dancer.stepper.tracking_state)
        log.info("Run time %s" % sec2str(time.time()-start_time))
        if args.log_ged:
            fname = get_output_filename(args.input, '.', 'ged_time', insert_time = True)
            ctrl.write_ged_over_time(open(fname.replace('.json', '.csv'), 'w'))
        if args.log_best:
            for chain, dancer in enumerate(dancers):
                fname = get_output_filename(args.input, '.', 'best%d.%s' % (
                        chain, initial_partition['rand_id']), insert_time = True)
                if out_format == 'json':
                    with open(fname, 'w') as fh:
                        #write_input_partition(dancer.stepper.tracking_state.best_partition, args, fh)
                        write_partition_indices(
                            dancer.stepper.tracking_state.best_partition, pool, args, fh)
                if out_format == 'hdf5':
                    index_part = IndexedPartition()
                    index_part.from_partition(dancer.stepper.tracking_state.best_partition)
                    with h5py.File(fname.replace('.json', '.hdf5')) as fh:
                        write_hdf5_observations([ (o.x, o.y, o.t) for o in pool] ,fh)
                        write_hdf5_indexed_partition(index_part, fh)


def record_samples(last_state, recorder, num_records):
    '''Given a burned-in sample, record samples'''
    log.info('start recording')
    num_steps = max(1, int(.25/last_state.recent_acceptance_rate))
    log.info('* num_steps per record = {0}'.format(num_steps))
    dancer = Tracker(Stepper(last_state.current_model_parameters, last_state.current_partition))
    #log.info(last_state.current_model_parameters.process_noise_covariance)
    dancer.sample_observation_error(1)
    for i in range(num_records):
        #if i == 0 or i == 1:
        #   log.info(dancer.state().current_model_parameters.process_noise_covariance)
        recorder.record_state(dancer.run(num_steps))

def track_stepper(initial_partition, args):
    '''Does tracking without burn-in control (uses Stepper)
    Stops after given time or given number of samples or keyboard interrupt.
    An eventual recording follows after the stop signal
    '''
    out_format = 'json'

    if has_flag(initial_partition, Flag.INDEXED):
        pass
    else:
        log.info('convert to indexed partition')
        initial_partition = indexed_partition(initial_partition)

    total_obs =  len(initial_partition['observations'])

    if args.max:
        maxi_tracks, maxi_clutter = max_partition(initial_partition['observations'])

    pdfs_len = 200
    current_pdfs = fifo(pdfs_len)

    if not args.norecord:
        '''open the file asap'''
        directory  = get_tracking_dir(args.input)
        output_file = open(get_output_filename(args.input, directory,
                                               "tracked." + initial_partition['rand_id']), 'w')
        recorder = SampleRecorder(initial_partition)

    initial_model = init_parameters()
    initial_model.initQ(args.process_noise*2)

    have_sims, input_partition, sims_partition, pool = convert_initial_partition(initial_partition)
    singletons = make_mini_partition(initial_partition['observations'], pool)


    if have_sims:
        ged_list = fifo(200)
        if not args.norecord:
            recorder.set_sims_partition(sims_partition)

    # Create a new tracking server
    dancer = None
    if args.max:
        maxi_partition = IndexedPartition()
        maxi_partition.set_tracks(maxi_tracks)
        maxi_partition.set_clutter(maxi_clutter)
        dancer = Stepper(initial_model, maxi_partition.to_partition(pool))
    else:
        #dancer = Stepper(initial_model, input_partition)
        dancer = Stepper(initial_model, singletons)

    print_r(dancer.tracking_state.current_model_parameters, log.info)

    state_file = None
    if args.log_states is not None:
        state_file = safe_open(args.log_states, 'w', args.unsafe, log.error)
        log.info('Writing states one-per-line to: %s' % (args.log_states,))

    start_time = time.time()
    #tick = 0
    #tick_header_interval = 10
    #dancer.step(1000)
    state = dancer.tracking_state
    stored_log_pdf = log_pdf(state.current_partition, initial_model)
    last_best_partition = state.current_partition
    header_counter = 0
    header_lag = 15
    seconds_passed = 0
    one_hour = 3600
    count_chunk = 100000 # for reporting
    samples_passed = 0 # for reporting
    try:
        while True:
            dancer.step(1000)
            state = dancer.tracking_state
            current_pdfs.append(state.current_log_pdf)
            tick_time = time.time() - start_time
            if have_sims:
                ged_list.append(state.current_partition.ged(sims_partition))

            if state.best_log_pdf > stored_log_pdf or (tick_time - seconds_passed) > one_hour or \
                    (state.sample_count-samples_passed)>count_chunk:
                while (tick_time - seconds_passed) > one_hour:
                    seconds_passed += one_hour
                while (state.sample_count-samples_passed)>count_chunk:
                    samples_passed += count_chunk
                if header_counter % header_lag == 0:
                    print_header(total_obs, state, have_sims)
                header_counter = (header_counter + 1) % header_lag
                hist = list(state.move_histogram)
                accepted = float(sum(hist[1:-1]))
                total = float(sum(hist))
                assert(total == state.sample_count)
                total /= 100
                improvement = 100.0*(stored_log_pdf - state.best_log_pdf)/stored_log_pdf
                ged = last_best_partition.ged(state.best_partition)
                if have_sims:
                    ged = min(ged_list)
                    #ged = state.best_partition.ged(sims_partition)
                ged_str = four2str(ged)
                last_best_partition = state.best_partition
                log.info('%12s %7s %8s %6.2f %6.2f %6.2f %8.1f %6.1f%% %4d %6.2f%% %4s' % (
                    sec2str(int(tick_time)), count2str(state.sample_count),
                    float2str(state.best_log_pdf),
                    float(hist[0]) / total, float(hist[-1])/total, accepted/total,
                    state.sample_count/float(tick_time),
                    100.0*float(len(state.best_partition.clutter))/float(total_obs),
                    len(state.best_partition.tracks),
                    improvement, ged_str))
                stored_log_pdf = state.best_log_pdf

            # Log state if required
            if state_file is not None:
                # Add time field to state
                state_d = {
                        'sample_count': state.sample_count,
                        'log_pdf': state.current_log_pdf,
                        'partition': json.loads(state.current_partition.to_json()),
                        'model_parameters': json.loads(state.current_model_parameters.to_json()),
                        'move_type': int(state.current_move_type),
                        'acceptance_rate': state.acceptance_rate,
                        'move_histogram': list(state.move_histogram),
                }
                json.dump(state_d, state_file)
                state_file.write('\n')

            if args.duration and args.duration < tick_time:
                log.info('Stopping since duration has exceeded %s seconds.' % (args.duration,))
                break

            if args.samples and number_parser(args.samples, int) < state.sample_count:
                log.info('Stopping since sample count has exceeded %s.' % (args.samples,))
                break

            if args.mean_pdf and len(current_pdfs) == pdfs_len and mean(current_pdfs) > -args.mean_pdf:
                log.info('Stopping since mean log PDF exceeded %.1f' % (-args.mean_pdf,))
                break

            #tick+=1
    except KeyboardInterrupt:
        log.warning('Terminating tracking due to keyboard interrupt.')
        return 1
    finally:
        state = dancer.tracking_state
        # del dancer
        if not args.norecord:
            recorder.set_burnin_length(state.sample_count)
            record_samples(state, recorder, args.records)
            if out_format == 'json':
                recorder.write_results(args, output_file)
            if out_format == 'hdf5':
                hdf5name = get_output_filename(args.input, 'tracking')
                with h5py.File(hdf5name.replace('.json', '.hdf5')) as fh:
                    recorder.write_hdf5_results(args, fh)
        print_final_state(log.info, state)
        if args.log_best:
            fname = get_output_filename(
                args.input, '.', 'best.' + initial_partition['rand_id'], insert_time = True)
            if out_format == 'json':
                with open(fname, 'w') as fh:
                    #write_input_partition(state.best_partition, args, fh)
                    write_partition_indices(state.best_partition, pool, args, fh)
            if out_format == 'hdf5':
                index_part = IndexedPartition()
                index_part.from_partition(state.best_partition)
                with h5py.File(fname.replace('.json', '.hdf5')) as fh:
                    write_hdf5_observations([ (o.x, o.y, o.t) for o in pool] ,fh)
                    write_hdf5_indexed_partition(index_part, fh)
        log.info("Run time %s" % sec2str(time.time()-start_time))

    ## Stop the tracker.
    #tracker.stop()

    return 0

def main(args):

    signal.signal(signal.SIGUSR1, signal_handler)
    signal.signal(signal.SIGUSR2, signal_handler)
    warnings.simplefilter("error") # turn warnings into errors
    rand_id = get_random_string(4)
    if args.log and  args.input is not None:
        #logfile = get_output_filename(args.input, 'logs', None, insert_time = True )
        logfile = get_output_filename(args.input, 'logs', rand_id, insert_time = True )
        logfile = logfile.replace('.json', '.log')
        handler = logging.StreamHandler(safe_open(logfile, 'w', args.unsafe, log.error))
        handler.setFormatter(logging.Formatter('%(levelname)s: %(message)s'))
        log.info('Writing log to: %s' % (logfile,))
        log.addHandler(handler)

    log.info('Call options:')
    for key, val in vars(args).items():
        if key != 'func':
            log.info("{0}: {1}".format(key, val))
    log.info('--')

    # Load the input partition
    if args.input is not None:
        if args.input.endswith('.json'):
            with open(args.input, 'r') as input_file:
                initial_partition = json.load(input_file)
        elif args.input.endswith('.hdf5'):
            with h5py.File(args.input, 'r') as input_file:
                initial_partition = input_from_hdf5(input_file)
    else:
        initial_partition = json.load(sys.stdin)
    initial_partition['rand_id'] = rand_id
    if args.sim_pdf:
        simulation_posterior(initial_partition, args)
    elif args.dual:
        track_dual(initial_partition, args)
    else:
        track_stepper(initial_partition, args)
    log.info(datetime.datetime.now().strftime("%c"))

if __name__ == '__main__':
    try:
        # Parse the command line options
        parser = argparse.ArgumentParser(description=globals()['__doc__'])
        setup_parser(parser)
        args = parser.parse_args()

        # Run the main  program
        sys.exit(main(args))
    except UnsafeOpenError:
        sys.exit(1)
