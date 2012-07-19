from itertools import izip
import cPickle as pickle
import os
import re
import sys
from utils import *

__author__ = 'pf'

if __name__ == '__main__':
    if len(sys.argv) != 2:
        error('usage: %s output_dir' % __file__)
    output_dir = sys.argv[1]
    statesbyline_dir = os.path.join(output_dir, 'STATEBYLINE')

    files = {}
    for fname in sorted(os.listdir(statesbyline_dir)):
        chrom = re.search(r'\w+_\d+_(\w+)_statebyline', fname).group(1)
        if chrom not in files:
            files[chrom] = []
        files[chrom].append(os.path.join(statesbyline_dir, fname))

    bins = {}
    total = 0
    change = 0
    for chrom in sorted(files):
        inputf = map(open, files[chrom])

        # read out the two header lines
        map(lambda f: f.readline(), inputf)
        map(lambda f: f.readline(), inputf)

        bins[chrom] = {}
        bin_no = 0
        for l1, l2, l3, l4 in izip(*inputf):
            if len(set([l1, l2, l3, l4])) != 1:
                bins[chrom][bin_no] = tuple(map(int, [l1, l2, l3, l4]))
            bin_no += 1

        total += bin_no
        change += len(bins[chrom])
        elapsed(chrom+'\ttotal: %d\tchange:%d' % (bin_no, len(bins[chrom])))


        map(lambda f: f.close(), inputf)
    print "Total:%d\tChange:%d\tFraction:%.2lf%%" % (total, change, 100*float(change)/total)
    pickle.dump(bins, open(os.path.join(output_dir, '037_bins_that_change_state.data'), 'w'))