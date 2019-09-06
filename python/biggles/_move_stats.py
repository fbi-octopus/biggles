# -*- coding: utf-8 -*-
"""
Created on Wed Mar  9 15:42:44 2016

@author: vcn81216
"""

from biggles import MoveType, move_name, move_sign

def print_move_hist(out, counts):
    '''prints simple move statistics'''
    total = sum(counts)
    out("Number of moves = {0}".format(total))
    
    for idx, count in enumerate(counts):
        name = move_name(MoveType(idx))[:3]
        out("{0} => {1:6.2f}%".format(name, 100.0*float(count)/float(total)))
        
def print_move_stats(out, accepted, identity, rejected):
    '''print the statistics of moves'''
    count = len(accepted)
    off_set = MoveType.BIRTH
    out("{0:^10s} {1:8s} {2:8s} {3:8s} {4:>10s}".format(
        "move", "rejected", "identity", "accepted", "total"
    ))
    for i in range(count):
        acc = accepted[i]
        idy = identity[i]
        rej = rejected[i]
        total = acc + idy + rej
        name = move_name(MoveType(i+off_set))
        if total == 0:
            out("{0:10s} {1:8s} {2:8s} {3:8s} {4:10d}".format(
                name, "N/A", "N/A", "N/A", total
            ))
        else:
            out("{0:10s} {1:7.1f}% {2:7.1f}% {3:7.1f}% {4:10d}".format(
                name, 100.0*rej/total, 100.0*idy/total, 100.0*acc/total, total
            ))
            
def get_acceptance_rate(accepted, identity, rejected):
    '''returns a dict of acceptance rate per move'''            
    count = len(accepted)
    off_set = MoveType.BIRTH
    acc_rate = dict()
    for i in range(count):
        acc = accepted[i]
        idy = identity[i]
        rej = rejected[i]
        total = acc + idy + rej
        name = move_sign(MoveType(i+off_set))
        if total == 0:
            acc_rate[name] = None
        else:
            acc_rate[name] = float(acc)/float(total)
    return acc_rate

def print_sampling_summary(out, state, tick_time):
    '''print a sampling summary'''
    fields = ('Time (min)', 'Samples', 'best PDF') 
    field2 = ('Rejected', 'Identity', 'Accepted', 'Rate')
    out('-' * (17 * len(fields) - 1) + '-' + '-' * (9 * len(field2) - 1))
    out((' '.join(['%16s'] * len(fields))) % fields + ' ' + (' '.join(['%8s'] * len(field2))) % field2)
    out('-' * (17 * len(fields) - 1) + '-' + '-' * (9 * len(field2) - 1))
    hist = list(state.move_histogram)
    accepted = float(sum(hist[1:-1]))
    total = float(sum(hist))
    total /= 100
    out('%16.2f %16d %16.2f %8.2f %8.2f %8.2f %8.1f' % (
                    tick_time/60.0, state.sample_count, state.best_log_pdf,
                    float(hist[0]) / total, float(hist[-1])/total, accepted/total, 
                    state.sample_count/float(tick_time)))    
