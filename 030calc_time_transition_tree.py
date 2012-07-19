import json
import os
from pprint import pformat
import re

__author__ = 'pf'
from utils import *

def add_to_tree(tree, probs):

    for bin_no in xrange(len(probs[0])):
        c_tree = tree

        for t in xrange(len(probs)):
            state, posterior = probs[t][bin_no]

            if state not in c_tree:
                c_tree[state] = {'_posterior' : 0, '_total' : 0}

            c_tree = c_tree[state]
            c_tree['_posterior'] += posterior
            c_tree['_total'] += 1



def calc_average_posterior(tree):
    for s in tree:
        if s[0] != '_':
            tree[s]['_average_posterior'] = tree[s]['_posterior'] / tree[s]['_total']
            calc_average_posterior(tree[s])



if __name__ == '__main__':
    if len(sys.argv) != 2:
        error('usage: %s output_dir' % __file__)
    output_dir = sys.argv[1]
    posterior_dir = os.path.join(output_dir, 'POSTERIOR')

    tree = {}
    c_chrom = None
    for fname in sorted(os.listdir(posterior_dir),
                        key = lambda fn: re.search(r't(\d+)_\d+_(\w+)_posterior', fn).groups()[::-1]):
        print fname
        tpoint, chrom = re.search(r't(\d+)_\d+_(\w+)_posterior', fname).groups()

        if chrom != c_chrom:
            if c_chrom is not None:
                add_to_tree(tree, probs)
                elapsed(c_chrom)
            probs = []
            c_chrom = chrom



        probs.append([])

        with open(os.path.join(posterior_dir, fname)) as f:
            header = f.readline()
            states = f.readline().strip().split('\t')
#            total = 0
            for line in f:
#                total += 1
#                if total == 1000: break

                max_p = 0
                max_i = 0
                for state_i, prob in enumerate(map(float, line.strip().split())):
                    if max_p <= prob:
                        max_p = prob
                        max_i = state_i

                probs[-1].append((states[max_i], max_p))



    add_to_tree(tree, probs)
    elapsed(c_chrom)
    calc_average_posterior(tree)
    json.dump(tree, open(os.path.join(output_dir, '030_tree.json'), 'w'), indent = 1)
    #print pformat(tree)
    elapsed('done')