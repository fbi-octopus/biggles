# -*- coding: utf-8 -*-
"""
some output stuff for biggles

Created on Tue Apr 26 12:20:24 2016

@author: vcn81216
"""

from _move_stats import print_move_stats

def print_final_state(out, state):
    '''does a final print of a trackers state'''
    best_part = state.best_partition
    out("Sample count = %d" % state.sample_count)
    out("Accteptance rate = %6.2f%%" % (state.acceptance_rate*100.0))
    out("Best PDF = %.2f" % state.best_log_pdf)
    out("Best partition:")
    out("   number of tracks = %d" % len(best_part.tracks))
    num_clutter = len(best_part.clutter)
    out("   number of clutter = %d, rate = %6.2f%%" % (
        num_clutter, 100.0 * float(num_clutter)/best_part.observation_count()))
    out("Best parameters:")
    factor = (100.0*100.0)/best_part.volume()
    best_para = state.best_model_parameters
    out("   b = %.4f" % (best_para.mean_new_tracks_per_frame*factor))
    out("   c = %.4f" % (best_para.mean_false_observations_per_frame*factor))
    out("   o = %.4f" % (best_para.generate_observation_probability))
    surv_prob = best_para.frame_to_frame_survival_probability
    out("   s = %.4f (%.2f)" % (surv_prob, 1./(1.-surv_prob)))
    obs_cov = best_para.observation_error_covariance
    out("   R = [%7.4f, %7.4f]" % (obs_cov[0][0], obs_cov[0][1]))
    out("       [%7.4f, %7.4f]" % (obs_cov[1][0], obs_cov[1][1]))
    pc = best_para.process_noise_covariance
    out("   Q = [{0:7.4f}, {1:7.4f}, {2:7.4f}, {3:7.4f}]".format(*pc[0]))
    out("       [{0:7.4f}, {1:7.4f}, {2:7.4f}, {3:7.4f}]".format(*pc[1]))
    out("       [{0:7.4f}, {1:7.4f}, {2:7.4f}, {3:7.4f}]".format(*pc[2]))
    out("       [{0:7.4f}, {1:7.4f}, {2:7.4f}, {3:7.4f}]".format(*pc[3]))
    print_move_stats(out, state.moves_accepted,
                     state.moves_identity, state.moves_rejected)

def print_parameters(out, parameters, volume):
    '''Prints parameters'''
    factor = (100.0*100.0)/volume
    out("   b = %.4f" % (parameters.mean_new_tracks_per_frame*factor))
    out("   c = %.4f" % (parameters.mean_false_observations_per_frame*factor))
    out("   o = %.4f" % (parameters.generate_observation_probability))
    surv_prob = parameters.frame_to_frame_survival_probability
    out("   s = %.4f (%.2f)" % (surv_prob, 1./(1.-surv_prob)))
    obs_cov = parameters.observation_error_covariance
    out("   R = [%7.4f, %7.4f]" % (obs_cov[0][0], obs_cov[0][1]))
    out("       [%7.4f, %7.4f]" % (obs_cov[1][0], obs_cov[1][1]))
    pc = parameters.process_noise_covariance
    out("   Q = [{0:7.4f}, {1:7.4f}, {2:7.4f}, {3:7.4f}]".format(*pc[0]))
    out("       [{0:7.4f}, {1:7.4f}, {2:7.4f}, {3:7.4f}]".format(*pc[1]))
    out("       [{0:7.4f}, {1:7.4f}, {2:7.4f}, {3:7.4f}]".format(*pc[2]))
    out("       [{0:7.4f}, {1:7.4f}, {2:7.4f}, {3:7.4f}]".format(*pc[3]))

def print_partition_summary(out, partition):
    '''Prints partition summary'''
    out("Number of tracks = %d" % len(partition.tracks))
    num_clutter = len(partition.clutter)
    out("Number of clutter = %d, rate = %6.2f%%" % (
        num_clutter, 100.0 * float(num_clutter)/partition.observation_count()))
    out("Volume = %d" % partition.volume())
