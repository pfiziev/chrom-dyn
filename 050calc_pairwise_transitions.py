#!/usr/bin/python
import copy
import json
import os

import pprint
import sys
import math
from utils import *
__author__ = 'pf'
import itertools

WINDOW_LENGTH = 10000


def weight(source_start, source_end, target_start, target_end):
    return (min(source_end, target_end) - max(source_start, target_start) + 1)/float(target_end - target_start + 1)

def distance(source_start, source_end, target_start, target_end):
    """ return 0 if the two regions overlap.
        return 1 + distance/BIN_SIZE, otherwise.
    """
    return 0 if partial_overlap(source_start, source_end, target_start, target_end) else 1 + (min(abs(target_start - source_end - 1), abs(source_start - target_end - 1)))/BIN_SIZE



if __name__ == '__main__':
    if len(sys.argv) != 2:
        error('usage: %s output_dir' % __file__)
    output_dir = sys.argv[1]
    segments = {}
    time_points = []
    transitions = json.load(open(os.path.join(output_dir, 'transitions_010.json')))

    _transitions = copy.deepcopy(transitions)


    segments = {}
    time_points = []

    states = set()
    for fname in sorted(os.listdir(output_dir)):
        if 'segments' not in fname: continue

        segments[fname] = {}
        time_points.append(fname)

        for line in open(os.path.join(output_dir, fname)):
            chrom, start, end, state = line.strip().split()
#            if chrom != 'chr22': continue

            if chrom not in segments[fname]:
                segments[fname][chrom] = []

            segments[fname][chrom].append((int(start), int(end) - 1, state))
            states.add(state)

        elapsed(fname)
    states = list(sorted(states, key = lambda s: int(s[1:])))
    nstates = len(states)
    ntimepoints = len(time_points)

    state_i = lambda s: states.index(s)
    elapsed('reading states')

    # initialize the distance weights with 1 - bin_no/bins_per_win
    bins_per_win = 1 + WINDOW_LENGTH/(2*BIN_SIZE)
    distance_weights = [[[1 - float(d)/bins_per_win for d in xrange(bins_per_win)]
                                for state in xrange(nstates)]
                                  for t in xrange(ntimepoints - 1)]



    START, END, STATE = 0, 1, 2

    # calculate the state occurrences at each time point
    occurrences = matrix(ntimepoints, nstates)
    for t in xrange(ntimepoints):
        for chrom in segments[time_points[t]]:
            for start, end, state in segments[time_points[t]][chrom]:
                occurrences[t][state_i(state)] += end - start + 1

    genome_length = sum(occurrences[0])

    elapsed('calculating state frequences')

    # calculate joint_occurrences and transition probabilities
    joint_occurrences = [None] * (ntimepoints - 1)
    for t in xrange(ntimepoints - 1):
        delta = 1
        tpoint1 = segments[time_points[t]]
        tpoint2 = segments[time_points[t+1]]

        print 'time points: ', t, t + 1
        while delta > 0.0001:
            trans = transitions[t]
            dweights = distance_weights[t]

            diff_comparisons = 0
            total_comparisons = 0
#            new_dweights = [[1]*len(dweights[0]) for state in xrange(nstates)]

            j_occurs = matrix(nstates, nstates)

            for chrom in tpoint1:

                tp1_segments = tpoint1[chrom]
                tp2_segments = tpoint2[chrom]

                chrom_length = tp1_segments[-1][END]

                t1_window_start_index = 0
                t1_window_end_index = 0
                t1_overlap_index = 0

                for start, end, state in tp2_segments:
                    state_length = end - start + 1
                    state_index = state_i(state)

#                    middle = start + (end - start)/2
#
#                    while t1_overlap_index < len(tp1_segments) and not (tp1_segments[t1_overlap_index][START] <= middle <= tp1_segments[t1_overlap_index][END]):
#                        t1_overlap_index += 1
#
                    while t1_window_start_index < len(tp1_segments) and not (tp1_segments[t1_window_start_index][START] <= start <= tp1_segments[t1_window_start_index][END]):
                        t1_window_start_index += 1

                    # find the window boundaries
#                    t1_window_start_index = t1_overlap_index
                    t1_window_end_index = t1_window_start_index

#                    while t1_window_start_index > 0 and \
#                          not (tp1_segments[t1_window_start_index][START] <= max(middle - WINDOW_LENGTH/2, 0) <= tp1_segments[t1_window_start_index][END]):
#                        t1_window_start_index -= 1
#
#                    while t1_window_end_index < len(tp1_segments) and \
#                          not (tp1_segments[t1_window_end_index][START] <= min(middle + WINDOW_LENGTH/2, chrom_length) <= tp1_segments[t1_window_end_index][END]):
#                        t1_window_end_index += 1

                    while t1_window_end_index < len(tp1_segments) and \
                          not (tp1_segments[t1_window_end_index][START] <= end <= tp1_segments[t1_window_end_index][END]):
                        t1_window_end_index += 1

                    window = tp1_segments[t1_window_start_index : t1_window_end_index + 1]

                    # calculate the most likely source
#                    mle_score, mle_source = max((weight(source_start, source_end, start, end) * trans[state_i(source)][state_index],
#                                                 source) for source_start, source_end, source in window)

                    mle_distance, mle_source = max(((weight(source_start, source_end, start, end), source)
                                                        for source_start, source_end, source in window),
                                                    key = lambda (w, s):  w * trans[state_i(s)][state_index])

                    mle_distance1, mle_source1 = max(((1, source)
                                                        for source_start, source_end, source in window),
                                                    key = lambda (w, s):  w * trans[state_i(s)][state_index])


                    if mle_source1 != mle_source:
                        diff_comparisons += 1
                    total_comparisons += 1
                    j_occurs[state_i(mle_source)][state_index] += state_length
#                    new_dweights[state_index][mle_distance] += 1

            print total_comparisons, diff_comparisons, 100*float(diff_comparisons)/total_comparisons
            new_trans = matrix(nstates, nstates)
            for i in xrange(nstates):
#                total_w = float(sum(new_dweights[i]))
#                for j in xrange(bins_per_win):
#                    new_dweights[i][j] /= total_w

                total_i = float(sum(j_occurs[i]))
                for j in xrange(nstates):
                    new_trans[i][j] = j_occurs[i][j] / total_i

            # calculate transitions change
            delta = math.sqrt(sum((new - old)**2
                                    for new_r, old_r in itertools.izip(new_trans, trans)
                                        for new, old in itertools.izip(new_r, old_r)) / (len(trans)**2))

#            delta_distances = math.sqrt(sum((new - old)**2
#                                            for new_r, old_r in itertools.izip(new_dweights, dweights)
#                                                for new, old in itertools.izip(new_r, old_r)) / (nstates * bins_per_win))

            print delta   #, delta_distances
            transitions[t] = new_trans
            joint_occurrences[t] = j_occurs

#            distance_weights[t] = new_dweights


    enrichments = [matrix(nstates, nstates) for t in xrange(ntimepoints - 1)]
    for t in xrange(ntimepoints - 1):
        for s1 in xrange(nstates):
            for s2 in xrange(nstates):
                enrichments[t][s1][s2] = joint_occurrences[t][s1][s2]/(occurrences[t][s1]*occurrences[t+1][s2]/float(genome_length))

    json.dump(transitions, open(os.path.join(output_dir, 'em_transitions_050.json'), 'w'), indent = 1)
    json.dump(enrichments, open(os.path.join(output_dir, 'em_enrichments_050.json'), 'w'), indent = 1)
    elapsed('output files: ' + os.path.join(output_dir, 'em_transitions_050.json') + ', ' + os.path.join(output_dir, 'em_enrichments_050.json'))


"""
# PLOT ENRICHMENTS

import json
import math
trans = json.load(open('/home/pf/hoffman2/chromHMM/adipogenesis/human_control_in_binarization/output/em_enrichments_050.json'))
from matplotlib.pylab import *
close_to_zero = 10**-5
def samplemat(d):
    orig_max_enr = max(e for r in d for e in r)
    orig_min_enr = min(e or close_to_zero for r in d for e in r)
    range_end = min(math.log(orig_max_enr,2),  -math.log(orig_min_enr,2))
    print range_end, math.log(orig_max_enr,2),  -math.log(orig_min_enr,2)
    aa = zeros((len(d), len(d)))
    for i in range(len(d)):
        for j in range(len(d[i])):
            aa[i,j] = min(max(math.log(d[i][j] or close_to_zero, 2), -range_end), range_end)
    return aa

for d in trans:
    matshow(samplemat(d), cmap = cm.seismic)
"""