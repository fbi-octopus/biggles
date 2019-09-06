# -*- coding: utf-8 -*-
"""
Biggles two chain convergence control

Created on Wed Jun 22 14:26:31 2016
@author: vcn81216
"""
from __future__ import print_function
from _utilities import count2str, float2str, sec2str, fifo, four2str, percentstr
#from _utilities import autocorr
from _utilities import auto_corr_time, autocorr, sample_pairs
from numpy import mean, std, array, sqrt
from scipy.stats import f_oneway
import csv
from vector_partition import PV2 as PartitionVector
from biggles import IndexedPartition

def sample_cross_part_ged(parts1, parts2):
    '''samples len(part1) cross GED between part1 and part2 (without repetition)'''
    assert(len(parts1) == len(parts2))
    num = len(parts1)
    return array([parts1[i].ged(parts2[j]) for i, j in sample_pairs(num, num)])

def get_cross_part_ged(parts1, parts2):
    result = list()
    for p1 in parts1:
        for p2 in parts2:
            result.append(p1.ged(p2))
    return array(result)

def get_inner_part_ged(part_list):
    num = len(part_list)
    result = list()
    for i in range(num-1):
        for j in range(i+1, num):
            result.append(part_list[i].ged(part_list[j]))
    assert(len(result) == num*(num-1)/2)
    return array(result)

def get_mean_displ_ged(part_list, dt):
    ''' mean displacement ged '''
    num = len(part_list)
    assert (dt < num/2)
    assert(isinstance(dt, int))
    result = list()
    for i in range(num-dt):
        result.append(part_list[i].ged(part_list[i+dt]))
    return mean(result)

def get_traveled_ged(part_list):
    num = len(part_list)
    total = 0.0
    for i in range(num-1):
        total += part_list[i].ged(part_list[i+1])
    crow = float(part_list[0].ged(part_list[-1]))
    return crow/total

class ControlData(object):
    ''' print info and test for convergence '''
    def __init__(self, start_time, total_obs, method, fifo_fact = 15):
        self._fifo_size = 1000 * fifo_fact
        self.header_counter = 0
        self.header_lag = 30
        self.stored_log_pdf = None
        self.last_sample_count = None
        self.start_time = start_time
        self.total_obs = total_obs
        self.last_output_time = 0
        self.last_log_time = 0
        #print("start time {0}".format(start_time))
        self.log_pdfs  = [fifo(self._fifo_size), fifo(self._fifo_size)]
        self.sam_rates = [fifo(self._fifo_size), fifo(self._fifo_size)]
        self.clu_rates = [fifo(self._fifo_size), fifo(self._fifo_size)]
        self.b_para = [fifo(self._fifo_size), fifo(self._fifo_size)]
        self.c_para = [fifo(self._fifo_size), fifo(self._fifo_size)]
        self.o_para = [fifo(self._fifo_size), fifo(self._fifo_size)]
        self.s_para = [fifo(self._fifo_size), fifo(self._fifo_size)]
        self.ged_pipe = fifo(self._fifo_size)
        self.logger = print
        self.method = method
        self.sim_part = None
        self.ged_over_time = list()
        self.best_parts = [None, None]
        self._start_partitions = [None, None]
        self.report_act = False
        self.report_ged = True
        num_parts = 1500
        self.part_samples = [ fifo(num_parts), fifo(num_parts)]
        self._anova_pval = 0.05
        self._gelman_thresh = 1.1
        self._ged_factor = 1.0
        self._log_writer = None

    def _similiar_data(self, data0, data1, fact = 1.0):
        '''tests how similar two data sets are'''
        return abs(mean(data0) - mean(data1)) <= fact * min(std(data0), std(data0))

    def report_ac_time(self, false_or_true):
        self.report_act = bool(false_or_true)

    def set_logger(self, logger):
        self.logger = logger

    def set_simulation(self, partition):
        '''sets the realisation of the ground truth for comparison purposes'''
        self.sim_part = partition

    def buffer_size(self):
        '''the size of the FIFOs used by the object'''
        return self._fifo_size

    def print_criterion(self):
        '''prints the stats for the convergence criterion'''
        self.logger("ANOVA convergence: {0}. GED convergence {1}.".format(
            self.anova_based_test(), self.ged_based_test()  ))
        self.logger("PDF: {0:.1f} +/- {1:.1f}; {2:.1f} +/- {3:.1f}".format(
            mean(self.log_pdfs[0]), std(self.log_pdfs[0]),
            mean(self.log_pdfs[1]), std(self.log_pdfs[1]),
            ))
        self.logger("clutter rate: {0:.1f}% +/- {1:.1f}%; {2:.1f}% +/- {3:.1f}%".format(
            mean(self.clu_rates[0])*100.0, std(self.clu_rates[0])*100.0,
            mean(self.clu_rates[1])*100.0, std(self.clu_rates[1])*100.0,
            ))
        self.logger("sample rate: {0:.1f} +/- {1:.1f}; {2:.1f} +/- {3:.1f}".format(
            mean(self.sam_rates[0]), std(self.sam_rates[0]),
            mean(self.sam_rates[1]), std(self.sam_rates[1]),
            ))
        min_ged = min(self.ged_pipe)
        avg_ged = mean(self.ged_pipe)
        max_ged = max(self.ged_pipe)
        self.logger("GED (parallel pipe): min %d (%.1f%%), mean %.1f (%.1f%%), max %d (%.1f%%)" % (
                min_ged, 100.0 * float(min_ged) / self._max_ged,
                avg_ged, 100.0 * float(avg_ged) / self._max_ged,
                max_ged, 100.0 * float(max_ged) / self._max_ged,
            ))
        #inner_ged1, inner_ged2 = [ mean(get_inner_part_ged(p)) for p in self.part_samples]
        inner_ged1, inner_ged2 = [ mean(sample_cross_part_ged(p, p)) for p in self.part_samples ]
        #cross_ged = mean(get_cross_part_ged(self.part_samples[0], self.part_samples[1]))
        cross_ged = mean(sample_cross_part_ged(self.part_samples[0], self.part_samples[1]))
        self.logger("GED Chain-1: %.1f (%.1f%%), GED Chain-2: %.1f (%.1f%%), Cross GED %.1f (%.1f%%)" % (
                inner_ged1, 100.0 * float(inner_ged1) / self._max_ged,
                inner_ged2, 100.0 * float(inner_ged2) / self._max_ged,
                cross_ged , 100.0 * float(cross_ged ) / self._max_ged,
            ))
        dist_best = self._states[0].best_partition.ged(self._states[1].best_partition)
        self.logger("GED best partitions: %.1f (%.1f%%)" % (
            dist_best, 100.0 * dist_best / self._max_ged))
        if self.sim_part is not None:
            sim_ged1, sim_ged2 = [ get_cross_part_ged(p, [self.sim_part]) for p in self.part_samples ]
            self.logger("GED to simulated partition:")
            sim2best0 = self._states[0].best_partition.ged(self.sim_part)
            sim2best1 = self._states[1].best_partition.ged(self.sim_part)
            self.logger("    Chain-0: best %d (%.1f%%), min %d (%.1f%%), avg %.1f (%.1f%%)" % (
                sim2best0, 100.0 * sim2best0 / self._max_ged,
                min(sim_ged1), 100.0 * min(sim_ged1) / self._max_ged,
                mean(sim_ged1), 100.0 * mean(sim_ged1) / self._max_ged,
            ))
            self.logger("    Chain-1: best %d (%.1f%%), min %d (%.1f%%), avg %.1f (%.1f%%)" % (
                sim2best1, 100.0 * sim2best1 / self._max_ged,
                min(sim_ged2), 100.0 * min(sim_ged2) / self._max_ged,
                mean(sim_ged2), 100.0 * mean(sim_ged2) / self._max_ged,
            ))
        ac0 = auto_corr_time(self.log_pdfs[0])
        ac1 = auto_corr_time(self.log_pdfs[1])
        self.logger("Auto correlation time PDF: %.2f and %.2f" % (ac0, ac1))
        intro = "Auto correlation (lag 1, 5, 10, 20)"
        self.logger(intro + ": %.2f, %.2f, %.2f, %.2f, and %.2f, %.2f, %.2f, %.2f" % (
            autocorr(list(self.log_pdfs[0]),  1),
            autocorr(list(self.log_pdfs[0]),  5),
            autocorr(list(self.log_pdfs[0]), 10),
            autocorr(list(self.log_pdfs[0]), 20),
            autocorr(list(self.log_pdfs[1]),  1),
            autocorr(list(self.log_pdfs[1]),  5),
            autocorr(list(self.log_pdfs[1]), 10),
            autocorr(list(self.log_pdfs[1]), 20),
            ))
        ac0 = auto_corr_time(self.b_para[0])
        ac1 = auto_corr_time(self.b_para[1])
        self.logger("Auto correlation time para B: %.2f and %.2f" % (ac0, ac1))
        ac0 = auto_corr_time(self.c_para[0])
        ac1 = auto_corr_time(self.c_para[1])
        self.logger("Auto correlation time para C: %.2f and %.2f" % (ac0, ac1))
        ac0 = auto_corr_time(self.o_para[0])
        ac1 = auto_corr_time(self.o_para[1])
        self.logger("Auto correlation time para O: %.2f and %.2f" % (ac0, ac1))
        ac0 = auto_corr_time(self.s_para[0])
        ac1 = auto_corr_time(self.s_para[1])
        self.logger("Auto correlation time para S: %.2f and %.2f" % (ac0, ac1))
        ac0 = auto_corr_time(self.clu_rates[0])
        ac1 = auto_corr_time(self.clu_rates[1])
        self.logger("Auto correlation time clutter rate: %.2f and %.2f" % (ac0, ac1))
        #
        return
        vp1 = []
        vp2 = []
        for part in self.part_samples[0]:
            indexed_partiton = IndexedPartition()
            indexed_partiton.from_partition(part)
            vp = PartitionVector(part.observation_count())
            vp.from_indexed_partition(indexed_partiton)
            vp1.append(vp)
        for part in self.part_samples[1]:
            indexed_partiton = IndexedPartition()
            indexed_partiton.from_partition(part)
            vp = PartitionVector(part.observation_count())
            vp.from_indexed_partition(indexed_partiton)
            vp2.append(vp)
        print('calc vp1act ...')
        vp1act = auto_corr_time(vp1)
        print('calc vp2act ...')
        vp2act = auto_corr_time(vp2)
        self.logger("Auto correlation time partitions: %.2f and %.2f" % (vp1act, vp2act))


    @property
    def states(self):
        return self._states

    @states.setter
    def states(self, new_states):
        self._states = new_states
        self.best_log_pdf = max([st.best_log_pdf for st in self._states])
        self.sample_count = min([st.sample_count for st in self._states])
        for i in [0, 1]:
            self.log_pdfs[i].append(self.states[i].current_log_pdf)
            clutter_count = float(len(self.states[i].current_partition.clutter))
            self.clu_rates[i].append(clutter_count/float(self.total_obs))
            self.part_samples[i].append(self.states[i].current_partition)
            paras = self.states[i].current_model_parameters
            self.b_para[i].append(paras.mean_new_tracks_per_frame)
            self.c_para[i].append(paras.mean_false_observations_per_frame)
            self.o_para[i].append(paras.generate_observation_probability)
            self.s_para[i].append(paras.frame_to_frame_survival_probability)
        if self.stored_log_pdf is None:
            self.stored_log_pdf = [st.best_log_pdf for st in self._states]
        self.ged_pipe.append(self.states[0].current_partition.ged(self.states[1].current_partition))

    def set_sampling_rates(self, sr):
        self.sampling_rates = sr
        for i in [0, 1]:
            self.sam_rates[i].append(self.sampling_rates[i])
    def set_max_ged(self, max_ged):
        self._max_ged = float(max_ged)

    def set_start_partitions(self, start_partitions):
        assert(len(start_partitions) == 2)
        self._start_partitions = start_partitions

    def try_output(self, tick_time):
        time_passed = tick_time - self.last_output_time
        log_time = tick_time - self.last_log_time
        if (time_passed > 20 or self.best_log_pdf > min(self.stored_log_pdf)) and self.sample_count > 0:
            # if tick_time - self.last_output_time < 1: return
            if self.header_counter % self.header_lag == 0:
                self.print_header()
            self.header_counter = (self.header_counter + 1) % self.header_lag
            self.print_state(tick_time)
            self.last_output_time = tick_time
            pdf_diff = mean(array(self.log_pdfs[0]) - array(self.log_pdfs[1]))
            cr_diff = mean(array(self.clu_rates[0]) - array(self.clu_rates[1]))
            self.ged_over_time.append(
                (tick_time, mean(self.ged_pipe)/self._max_ged, pdf_diff, cr_diff)
            )
            #print("{0} {1}".format(tick_time, log_time))
        if (log_time > 20) and not self._log_writer is None:
            self.write_log_line()
            self.last_log_time = tick_time


    def print_explaination(self, out):
        out('''Convergence criteria:''')
        out('''1) Gelman-Rubin:''')
        out('''   the Gelman-Rubin statistic for the parameters''')
        out('''   birth rate, clutter rate, observation probability, survival probability''')
        out('''   is less than {0}'''.format(self._gelman_thresh))
        out('''2) ANOVA:''')
        out('''   the p-value for the one-side ANOVA test for the observed clutter rate and the log-PDF''')
        out('''   is larger than {0}'''.format(self._anova_pval))
        out('''3) GED:''')
        out('''   the sum of the mean inner-chain GED is larger than''')
        out('''   {0} times the mean cross-chain GED.'''.format(self._ged_factor))
        out('''   Shown is (2*mean cross-chain GED) / (sum mean inner-chain GED)''')

    def print_header_explain(self, out):
        out('''Legend:''')
        out('''"Sample"       total number of samples drawn so far; sum of chains.''')
        out('''"max PDF"      max. aposteriory log density found so far in chains.''')
        out('''"Accept.R"     recent acceptance rate of the partition sampler for each chain.''')
        out('''"Clutter rate" clutter rate of the currently best partition sample for each chain.''')
        out('''"#Tr"          Gelman-Rubin for the number of tracks.''')
        out('''"S/s"          sample rate; mean of chains.''')
        out('''"GED"          (2*mean cross-chain GED) / (sum mean inner-chain GED)''')
        out('''"PDF"          Gelman-Rubin for the posterior PDF.''')
        out('''"GR b"         Gelman-Rubin for the birth rate.''')
        out('''"GR c"         Gelman-Rubin for the clutter rate.''')
        out('''"GR o"         Gelman-Rubin for the observation probability.''')
        out('''"GR s"         Gelman-Rubin for the survival probability.''')
        if  self.report_ged:
            out('''"Travelled"    The GED between the first and the last partition sample''')
            out('''               relative to the sum of GEDs from one partition to the next''')
            out('''               for each chain pipe.''')
            out('''               Distance as the crow flies relative to the actual travelled distance.''')


    def print_header(self):
        headstr = "{0:6d} obs {1:>6s} {2:7s} {3:^9s} {4:^13s} {5:^4s} {6:^4s} {7:4s}".format(
            self.total_obs, 'Sample', 'max PDF', 'Accept.R',
            'Clutter rate', '#Tr', 'S/s', " GED"
        )
        headstr += " %4s %4s %4s %4s %4s" % ("PDF", "GR b", "GR c", "GR o", "GR s")
        if self.report_act:
            headstr += " %9s %9s %9s %9s %9s %9s" % ("AC-PDF", "AC-Clu-ra", "AC-b-para",
                                                     "AC-c-para", "AC-o-para", "AC-s-para")
        if self.report_ged:
            headstr += " {0:9s}".format("Travelled")
        #intstr = " {0:>9s} {1:>9s}".format("crossable", "transible")
        intstr = ""
        self.logger('-' * len(headstr) + "-" * len(intstr))
        self.logger(headstr + intstr)
        self.logger('-' * len(headstr) + "-" * len(intstr))

    def burn_in_pdfs(self):
        '''returns the log PDFs of both chains as one list'''
        result = list(self.log_pdfs[0])
        for log_pdf in self.log_pdfs[1:]:
            result += list(log_pdf)
        return result

    def ged_based_test(self):
        #return mean(self.ged_pipe) <= sum([ mean(get_inner_part_ged(p)) for p in self.part_samples])

        #return mean(get_cross_part_ged(self.part_samples[0], self.part_samples[1])) <= sum(
        #    [ mean(get_inner_part_ged(p)) for p in self.part_samples])
        fact = self._ged_factor
        p0, p1 = self.part_samples
        return fact*mean(sample_cross_part_ged(p0, p1)) <= sum(
            [ mean(sample_cross_part_ged(p0, p0)) + mean(sample_cross_part_ged(p1, p1)) ])
    def anova_based_test(self):
        pval = self._anova_pval
        return f_oneway(*self.clu_rates)[1] > pval and f_oneway(*self.log_pdfs)[1] > pval

    def gelman_rubin_test(self):
        criterion = self._gelman_thresh
        return all([gelman_rubin(*self.b_para) < criterion,
                    gelman_rubin(*self.c_para) < criterion,
                    gelman_rubin(*self.o_para) < criterion,
                    gelman_rubin(*self.s_para) < criterion])

    def converged(self):
        '''tests if the two trackers have convered'''
        if self.method == 'ged':
            return self.ged_based_test()
        if self.method == 'gelman':
            return self.gelman_rubin_test()
        if self.method == 'anova':
            return self.anova_based_test()
        if self.method == 'ano_ged':
            return self.anova_based_test() and self.ged_based_test()
        if self.method == 'gel_ged':
            return self.gelman_rubin_test() and self.ged_based_test()
        raise('unknown convergence method: "{0}"'.format(self.method))

    def print_state(self, tick_time):
        data = {}
        for i, state in enumerate(self.states):
            record = {}
            #hist = list(state.move_histogram)
            #accepted = float(sum(hist[1:-1]))
            #total = float(sum(hist))
            record['sample_count'] = state.sample_count
            record['best_log_pdf'] = state.best_log_pdf
            #record['acceptance_rate'] = accepted/total*100.0
            record['acceptance_rate'] = state.recent_acceptance_rate
            record['clutter_rate'] = 100.0*float(len(state.best_partition.clutter))/float(self.total_obs)
            record['num_tracks'] = len(state.best_partition.tracks)
            record['sample_rate'] = self.sampling_rates[i]
            record['improvement'] = (self.stored_log_pdf[i] - state.best_log_pdf)/self.stored_log_pdf[i]
            record['intern-merge'] = state.internals[0]
            record['intern-cross'] = state.internals[1]
            record['intern-trans'] = state.internals[2]
            record['intern-extend'] = state.internals[3]
            self.stored_log_pdf[i] = state.best_log_pdf
            data[i] = record

        #ged = self.states[0].best_partition.ged(self.states[1].best_partition)
        #ged_std_str = percentstr(std(self.ged_pipe)/self._max_ged) + "%"
        #ged_str = four2str(ged)
        #cross_ged = get_cross_part_ged(self.part_samples[0], self.part_samples[1])
        cross_ged = sample_cross_part_ged(self.part_samples[0], self.part_samples[1])
        #ged_str = percentstr(mean(self.ged_pipe)/self._max_ged) + "%"
        #ged_str = percentstr(mean(cross_ged)/self._max_ged) + "%"
        #ged_str = percentstr(float(ged)/self._max_ged) + "%"
        #
        #geds1 = get_inner_part_ged(self.part_samples[0])
        #geds2 = get_inner_part_ged(self.part_samples[1])
        geds1 = sample_cross_part_ged(self.part_samples[0], self.part_samples[0])
        geds2 = sample_cross_part_ged(self.part_samples[1], self.part_samples[1])
        #geds1_str = percentstr(mean(geds1)/self._max_ged) + "%"
        #geds2_str = percentstr(mean(geds2)/self._max_ged) + "%"
        #
        ged_diff = 2.0 * mean(cross_ged) / (mean(geds1) + mean(geds2))
        #
        #mean_tr0 = mean([ len(part.tracks) for part in self.part_samples[0]])
        #mean_tr1 = mean([ len(part.tracks) for part in self.part_samples[1]])
        gr_tracks = gelman_rubin([ len(part.tracks) for part in self.part_samples[0]],
                                 [ len(part.tracks) for part in self.part_samples[1]])
        #
        fmtstr = '%10s %6s %7s %3s%% %3s%% %5.1f%% %5.1f%% %4s %4s %4s %4s %4s %4s %4s %4s'
        mainpart = fmtstr % (
            sec2str(int(tick_time)), count2str(data[0]['sample_count'] + data[1]['sample_count']),
            float2str(max(data[0]['best_log_pdf'], data[1]['best_log_pdf'])),
            percentstr(data[0]['acceptance_rate']),
            percentstr(data[1]['acceptance_rate']),
            data[0]['clutter_rate'], data[1]['clutter_rate'],
            #data[0]['num_tracks'], data[1]['num_tracks'],
            four2str(gr_tracks),
            four2str((data[0]['sample_rate'] + data[1]['sample_rate'])/2.0),
            four2str(ged_diff),
            four2str(gelman_rubin(*self.log_pdfs)),
            four2str(gelman_rubin(*self.b_para)),
            four2str(gelman_rubin(*self.c_para)),
            four2str(gelman_rubin(*self.o_para)),
            four2str(gelman_rubin(*self.s_para))
            )
        supplement = ''
        if self.report_act:
            ac0 = four2str(auto_corr_time(self.log_pdfs[0]))
            ac1 = four2str(auto_corr_time(self.log_pdfs[1]))
            supplement += '|%4s %4s' % (ac0, ac1)
            ac0 = four2str(auto_corr_time(self.clu_rates[0]))
            ac1 = four2str(auto_corr_time(self.clu_rates[1]))
            supplement += '|%4s %4s' % (ac0, ac1)
            ac0 = four2str(auto_corr_time(self.b_para[0]))
            ac1 = four2str(auto_corr_time(self.b_para[1]))
            supplement += '|%4s %4s' % (ac0, ac1)
            ac0 = four2str(auto_corr_time(self.c_para[0]))
            ac1 = four2str(auto_corr_time(self.c_para[1]))
            supplement += '|%4s %4s' % (ac0, ac1)
            ac0 = four2str(auto_corr_time(self.o_para[0]))
            ac1 = four2str(auto_corr_time(self.o_para[1]))
            supplement += '|%4s %4s' % (ac0, ac1)
            ac0 = four2str(auto_corr_time(self.s_para[0]))
            ac1 = four2str(auto_corr_time(self.s_para[1]))
            supplement += '|%4s %4s' % (ac0, ac1)
        ged_suppl = ''
        if self.report_ged:
            #traveled0 = percentstr(get_traveled_ged(self.part_samples[0]))
            #traveled1 = percentstr(get_traveled_ged(self.part_samples[1]))
            traveled0 = percentstr(self._start_partitions[0].ged(self.states[0].current_partition)/self._max_ged)
            traveled1 = percentstr(self._start_partitions[1].ged(self.states[1].current_partition)/self._max_ged)
            ged_suppl += ' {0:3s}% {1:3s}%'.format(traveled0, traveled1)
        sim_suppl = ''
        if self.sim_part is not None:
            sim_ged1, sim_ged2 = [ get_cross_part_ged(p, [self.sim_part]) for p in self.part_samples ]
            sim_suppl += ' {0:3s}% {1:3s}%'.format(
                percentstr(min(sim_ged1) / self._max_ged),
                percentstr(min(sim_ged2) / self._max_ged)
            )
        self.logger(mainpart + supplement + ged_suppl + sim_suppl)

            #+ " %4d %4d %4d %4d" % (
            #    data[0]['intern-cross'], data[1]['intern-cross'],
            #    data[0]['intern-trans'], data[1]['intern-trans'],
            #))
    def write_ged_over_time(self, fh):
        writer = csv.writer(fh)
        writer.writerows(self.ged_over_time)

    def write_log_line(self):
        ''' write a single log -line'''
        cross_ged = sample_cross_part_ged(self.part_samples[0], self.part_samples[1])
        geds1 = sample_cross_part_ged(self.part_samples[0], self.part_samples[0])
        geds2 = sample_cross_part_ged(self.part_samples[1], self.part_samples[1])
        ged_diff = 2.0 * mean(cross_ged) / (mean(geds1) + mean(geds2))
        line = [ self.states[0].sample_count,
                self.best_log_pdf, self.log_pdfs[0][-1], self.log_pdfs[1][-1],
                len(self.part_samples[0][-1].tracks), len(self.part_samples[1][-1].tracks),
                self.states[0].recent_acceptance_rate, self.states[1].recent_acceptance_rate,
                gelman_rubin(*self.b_para), gelman_rubin(*self.c_para),
                gelman_rubin(*self.o_para), gelman_rubin(*self.s_para),
                #get_traveled_ged(self.part_samples[0]), get_traveled_ged(self.part_samples[1]),
                self._start_partitions[0].ged(self.states[0].current_partition)/self._max_ged,
                self._start_partitions[1].ged(self.states[1].current_partition)/self._max_ged,
                ged_diff,
                self.b_para[0][-1], self.b_para[1][-1],
                self.c_para[0][-1], self.c_para[1][-1],
                self.o_para[0][-1], self.o_para[1][-1],
                self.s_para[0][-1], self.s_para[1][-1],
                get_det(self.states[0].current_model_parameters.observation_error_covariance),
                get_det(self.states[1].current_model_parameters.observation_error_covariance),
                f_oneway(*self.clu_rates)[1], f_oneway(*self.log_pdfs)[1]
        ]
        self._log_writer.writerow(line)

    def write_all_log(self, fname):
        header = ['sample', 'best_pdf', 'pdf1', 'pdf2', 'tracks1', 'tracks2', 'ar1', 'ar2', 'grb',
                  'grc', 'gro' , 'grs', 'travel1', 'travel2', 'ged',
                  'b1', 'b2', 'c1', 'c2', 'o1', 'o2', 's1', 's2', 'dr1', 'dr2',
                  'anova_cr', 'anova_pdf']
        self._fh = open(fname, 'w')
        self._log_writer = csv.writer(self._fh)
        self._log_writer.writerow(header)

def sample_variance(data):
    n = 0.0
    mean = 0.0
    M2 = 0.0

    for x in data:
        n += 1.0
        delta = x - mean
        mean += delta/n
        delta2 = x - mean
        M2 += delta*delta2

    if n < 2:
        return float('nan')
    else:
        return M2 / (n - 1)

def get_det(mat2x2):
    return mat2x2[0][0]*mat2x2[1][1]-mat2x2[0][1]*mat2x2[1][0]

def gelman_rubin(*args):
    ''' Calculate the Gelman Rubin statistic. Each argument is a chain of a parameter.
        Each chain is assumed to have the same length.
    '''
    n = float(len(args[0]))
    m = float(len(args))
    means = []
    variances = []
    for chain in args:
        means.append(mean(chain))
        variances.append(sample_variance(chain))
    B = sample_variance(means) # times n skipped
    W = mean(variances)
    if W == 0:
        return 999
    mar_post_var = (n-1)/n*W + B # divided by n skipped
    R = mar_post_var/W
    return sqrt((m+1)/m*R - (n-1)/(m*n))
