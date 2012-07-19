#!/usr/bin/python
import json
import os

import pprint
import sys
from utils import *
__author__ = 'pf'



if __name__ == '__main__':
    if len(sys.argv) != 2:
        error('usage: %s output_dir' % __file__)
    output_dir = sys.argv[1]
    segments = {}
    time_points = []
    states = set()
    for fname in sorted(os.listdir(output_dir)):
        if 'segments' not in fname: continue

        segments[fname] = {}
        time_points.append(fname)
        with open(os.path.join(output_dir, fname)) as f:
            for line in f:
                chrom, start, end, state = line.strip().split()
                if chrom not in segments[fname]: segments[fname][chrom] = []
                segments[fname][chrom].append((int(start), int(end), state))
                states.add(state)
        elapsed(fname)
    states = list(sorted(states, key = lambda s: int(s[1:])))
    nstates = len(states)
    state_i = lambda s: states.index(s)
    elapsed('reading states')

    total_bases = 0

    transitions = [matrix(nstates, nstates) for t in xrange(len(time_points) - 1)]
    enrichments = [matrix(nstates, nstates) for t in xrange(len(time_points) - 1)]

    START, END, STATE = 0, 1, 2
    for t in xrange(len(time_points) - 1):
        tpoint1 = segments[time_points[t]]
        tpoint2 = segments[time_points[t+1]]
        trans = transitions[t]
        enr = enrichments[t]

        for chrom in tpoint1:
            chrom_len = tpoint1[chrom][-1][END]
            s1_i = 0
            s2_i = 0
            for pos in xrange(0, chrom_len, BIN_SIZE):

                if pos > tpoint1[chrom][s1_i][END]: s1_i += 1
                if pos > tpoint2[chrom][s2_i][END]: s2_i += 1

                s1 = state_i(tpoint1[chrom][s1_i][STATE])
                s2 = state_i(tpoint2[chrom][s2_i][STATE])

                trans[s1][s2] += 1



        # sum over all rows for each state j (the total number of bins that have state j at time t2)
        total_j = [sum(trans[k][j] for k in xrange(nstates)) for j in xrange(nstates)]
        total_bases = sum(total_j)
        for i in xrange(nstates):
            total_i = float(sum(trans[i][j] for j in xrange(nstates)))
            for j in xrange(nstates):
                enr[i][j] = trans[i][j] / ((total_i*total_j[j])/total_bases)
                trans[i][j] /= total_i
#                trans[i][j] /= float(total_i*total_j[j])/total_bases

        elapsed('t%d -> t%d' % (t, t+1))

    json.dump(enrichments, open(os.path.join(output_dir, 'enrichments_010.json'), 'w'), indent = 1)
    json.dump(transitions, open(os.path.join(output_dir, 'transitions_010.json'), 'w'), indent = 1)
    elapsed('output files: ' + os.path.join(output_dir, 'enrichments_010.json')+', '+ os.path.join(output_dir, 'transitions_010.json'))

"""
## PLOTS ENRICHMENT IN LOG SCALE

import json
import math
trans = json.load(open('/home/pf/hoffman2/chromHMM/adipogenesis/human_control_in_binarization/output/enrichments_010.json'))
from matplotlib.pylab import *

def samplemat(d):
    orig_max_enr = max(e for r in d for e in r)
    orig_min_enr = min(e for r in d for e in r if e > 0)
    range_end = min(math.log(orig_max_enr,2),  -math.log(orig_min_enr,2))
    print range_end
    aa = zeros((len(d), len(d)))
    for i in range(len(d)):
        for j in range(len(d[i])):
            aa[i,j] = min(max(math.log(d[i][j] or 10**-5, 2), -range_end), range_end)
    return aa

for d in trans:
    matshow(samplemat(d), cmap = cm.seismic)

"""


"""
## PLOTS TRANSITION PROBABILITIES

import json
import math
trans = json.load(open('/home/pf/hoffman2/chromHMM/adipogenesis/human_control_in_binarization/output/transitions_010.json'))
from matplotlib.pylab import *

def samplemat(d):
    aa = zeros((len(d), len(d)))
    for i in range(len(d)):
        for j in range(len(d[i])):
            aa[i,j] = d[i][j]
    return aa

for d in trans:
    matshow(samplemat(d), cmap = cm.Blues)
"""