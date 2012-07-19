from itertools import izip
import os
import sys
from utils import *

__author__ = 'pf'




if __name__ == '__main__':
    if len(sys.argv) != 2:
        error('usage: %s output_dir' % __file__)
    output_dir = sys.argv[1]



    statebyline_dir = os.path.join(output_dir, 'STATEBYLINE')

    chrom = 'chr22'
    path = ('E8', 'E7', 'E8', 'E8')

    f1, f2, f3, f4 = sorted([open(os.path.join(statebyline_dir, f)) for f in sorted(os.listdir(statebyline_dir)) if '_'+chrom+'_' in f])

    in_path = False
    bin_no = -3
    for l1, l2, l3 , l4 in izip(f1, f2, f3, f4):
        bin_no += 1
        if bin_no < 0:
            # skip the two header lines in the beginning
            continue
        if path == tuple('E'+l.strip() for l in [l1, l2, l3, l4]):
            if not in_path:
                start_pos = bin_no*200
            in_path = True
        else:
            if in_path:
                print '%d-%d' % (start_pos, bin_no*200)
            in_path = False
        if bin_no > 100000:
            break



