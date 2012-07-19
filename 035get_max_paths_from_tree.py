import json
import os
import itertools

__author__ = 'pf'

from utils import *

if __name__ == '__main__':
    if len(sys.argv) != 2:
        error('usage: %s output_dir' % __file__)
    output_dir = sys.argv[1]

    tree = json.load(open(os.path.join(output_dir, '030_tree.json')))

    def path_score(path):
        t = tree
        score = 0
        for state in path:
            t = t.get(state, None)
            if t is None:
                return 0
            score += t['_average_posterior']
        return score/len(path)


    print sorted([(path_score(path), path) for path in itertools.product(['E%d' % i for i in xrange(1,11)], repeat=4)], reverse = True)[:10]