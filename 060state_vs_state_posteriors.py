from itertools import izip, chain
import os
import re
import sys
from utils import error, elapsed

if __name__ == '__main__':
    if len(sys.argv) != 2:
        error('usage: %s output_dir' % __file__)
    output_dir = sys.argv[1]
    posterior_dir = os.path.join(output_dir, 'POSTERIOR')

    states = [2, 3]

    timepoints = [1, 2]

    posteriors = {}
    tp_fnames = [[],[]]
    for fname in sorted(os.listdir(posterior_dir),
                        key = lambda fn: re.search(r't(\d+)_\d+_(\w+)_posterior', fn).groups()[::-1]):


        #print fname
        tpoint, chrom = re.search(r't(\d+)_\d+_(\w+)_posterior', fname).groups()
        if chrom != 'chr1':
            continue
        tpoint = int(tpoint)

        if tpoint not in timepoints:
            continue

        tp_fnames[timepoints.index(tpoint)].append(os.path.join(posterior_dir,fname))



    for tp1, tp2 in chain.from_iterable(izip(open(tp_fname1), open(tp_fname2)) for tp_fname1, tp_fname2 in izip(tp_fnames[0], tp_fnames[1])):

        if tp1.startswith('hASC'):
            print 'scanning', tp1.strip(), tp2.strip()
            continue
        if tp1[0] == 'E':
            continue

        values1 = map(float, tp1.strip().split())
        values2 = map(float, tp2.strip().split())
        max_i1, max_p1 = max(enumerate(values1), key = lambda (a,b): b)
        max_i2, max_p2 = max(enumerate(values2), key = lambda (a,b): b)

        if max_i1 in states and max_i2 in states:
            key = (max_i1, max_i2)
            if key not in posteriors:
                posteriors[key] = [[],[]]
            posteriors[key][0].append(max_p1)
            posteriors[key][1].append(max_p2)


    import numpy as np
    import matplotlib
    import matplotlib.pyplot as plt

    for key in posteriors:
        fname = os.path.join(output_dir, 'E%d_vs_E%d.png' % (key[0] + 1, key[1] + 1))
        print 'saving', fname
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.plot(posteriors[key][0], posteriors[key][1], '.')
        ax.set_title(fname)
        plt.savefig(fname)




#        if max_i in states:
#            for state_id in states:
#                probs[state_id].append(values[state_id])
#
#        probs[-1].append((states[max_i], max_p))
