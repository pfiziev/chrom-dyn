import json
import os
import re
import sys
import math
import operator
from utils import *
import numpy as np
import cPickle as pickle

__author__ = 'pf'


def get_n_states(fname):
    with open(fname) as f:
        _ = f.readline()
        states = f.readline().strip().split('\t')
        return len(states)

def count_bins(fname):
    bins = 0
    with open(fname) as f:
        _ = f.readline()
        _ = f.readline().strip().split('\t')
        for l in f:
            bins += 1
    elapsed('counting bins in %s: %d' % (fname, bins))
    return bins

def init_tree(n_time_points, states):
    tree = {}
    for state in states:
        subtree = tree
        for _ in xrange(n_time_points):
            subtree[state] = {}
            subtree = subtree[state]
    return tree



def iterpath(tree):
    stack = [([key], tree[key]) for key in tree]
    while len(stack) > 0:
        path, subtree = stack.pop()
        keys = [key for key in subtree if isinstance(key, int)]
        if not keys:
            yield path
        else:
            for key in keys:
                stack.append((path + [key], subtree[key]))



def calc_initial_best_paths(tree, posteriors):
    scores = dict((chrom,
                   np.zeros(len(posteriors[chrom])))
                        for chrom in posteriors)

    for chrom in posteriors:
        for bin_no, bin_probs in enumerate(posteriors[chrom]):
            scores[chrom][bin_no] = max(calc_path_score(path, bin_probs) for path in iterpath(tree))


    return scores

def calc_path_score(path, bin_probs):
#    score = 0
#    for tp, state in enumerate(path):
#        prob = bin_probs[tp][state]
#        if not prob:
#            return float('-inf')
#        score += math.log(prob)
#    return score
#    return sum(math.log(bin_probs[tp][state], 2) for tp, state in enumerate(path))
    score = 1
    for tp, state in enumerate(path):
        score *= bin_probs[tp][state]
    return score
#    return reduce(operator.mul, (bin_probs[tp][state] for tp, state in enumerate(path)))



def iter_new_paths(tree, states, max_depth):
    stack = [([key], tree[key]) for key in tree]
    while len(stack) > 0:
        path, subtree = stack.pop()
        keys = [key for key in subtree if isinstance(key, int)]
        if keys:
            for extension in [[state]*(max_depth - len(path)) for state in states if state not in subtree]:
                yield path + extension

            for key in keys:
                stack.append((path + [key], subtree[key]))


def add_path(tree, path):
    subtree = tree
    for state in path:
        if state not in subtree:
            subtree[state] = {}
        subtree = subtree[state]

def update_scores(scores, posteriors, best_path):
    # update the best scores according the new best_path
    for chrom in scores:
        for bin_no, bin_score in enumerate(scores[chrom]):
            scores[chrom][bin_no] = max(calc_path_score(best_path, posteriors[chrom][bin_no]),   # get the maximum between the score along the current path of that bin
                                        bin_score)                                               # and the previous best score for that bin

def optimize_tree(tree, posteriors, scores, n_paths, states, n_timepoints):
    for i in xrange(n_paths):
        best_path = max(iter_new_paths(tree, states, n_timepoints),
                        key = lambda path: sum(max(calc_path_score(path, bin_probs), # get the maximum between the score along the current path of that bin
                                                    scores[chrom][bin_no])           # and the previous best score for that bin
                                                        for chrom in posteriors
                                                            for bin_no, bin_probs in enumerate(posteriors[chrom])))
        update_scores(scores, posteriors, best_path)
        add_path(tree, best_path)
        elapsed('%d best path: %s' % (i , str(tuple(map(lambda s: 'E%d' % (s+1), best_path)))))



def get_n_timepoints(posterior_dir):
    # returns the number of time points
    tpoints = set()
    for fname in sorted(os.listdir(posterior_dir),
        key = lambda fn: re.search(r't(\d+)_\d+_(\w+)_posterior', fn).groups()[::-1]):
        fname = os.path.join(posterior_dir, fname)


        tpoint, _ = re.search(r't(\d+)_\d+_(\w+)_posterior', fname).groups()
        tpoints.add(tpoint)
    return len(tpoints)

if __name__ == '__main__':
    if len(sys.argv) != 2:
        error('usage: %s output_dir' % __file__)
    output_dir = sys.argv[1]
    posterior_dir = os.path.join(output_dir, 'POSTERIOR')

    with open(os.path.join(output_dir, '037_bins_that_change_state.data')) as cbf:
        changing_bins = pickle.load(cbf)
    elapsed('reading changing bins')

    posteriors = {}
    c_chrom = None
    n_states = None
    states = None

    n_timepoints = get_n_timepoints(posterior_dir)

    for fname in sorted(os.listdir(posterior_dir),
                        key = lambda fn: re.search(r't(\d+)_\d+_(\w+)_posterior', fn).groups()[::-1]):
        fname = os.path.join(posterior_dir, fname)

        if n_states is None:
            n_states = get_n_states(fname)

        _, chrom = re.search(r't(\d+)_\d+_(\w+)_posterior', fname).groups()
        if chrom not in ['chr22']:
            continue

        print fname

        if chrom != c_chrom:
            c_chrom = chrom

#            posteriors[c_chrom] = np.zeros((100, n_timepoints, n_states))
            posteriors[c_chrom] = np.zeros((len(changing_bins[chrom]), n_timepoints, n_states))
            tpoint = -1

        tpoint += 1

        with open(fname) as f:
            _ = f.readline()
            states = range(len(f.readline().strip().split('\t')))
            bin_no = 0
            bin_real_no = 0
            for line in f:
                if bin_real_no in changing_bins[chrom]:
                    bin_probs = [float(p) for p in line.strip().split()]
                    for s_i in xrange(n_states):
                        posteriors[c_chrom][bin_no][tpoint][s_i] = bin_probs[s_i] or 10**-10

                    bin_no += 1
                bin_real_no += 1
#                if bin_no == 100:
#                    break
        elapsed('chrom:%s\tbins:%d' % (c_chrom, bin_no))

    elapsed('reading posterior probabilities')


    tree = init_tree(n_timepoints, states)
#    scores = []
    scores = calc_initial_best_paths(tree, posteriors)
#    elapsed('calc_initial_scores')
    n_paths = 10
    optimize_tree(tree, posteriors, scores, n_paths, states, n_timepoints)
    elapsed('optimize_tree')
    pickle.dump(tree, open(os.path.join(output_dir, '040_optimal_tree.data'), 'w'))

#    from guppy import *
#    h = hpy()
#    print h.heap()
