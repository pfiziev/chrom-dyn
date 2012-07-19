import os
from utils import *

if __name__ == '__main__':
    if len(sys.argv) != 2:
        error('usage: %s output_dir' % __file__)
    output_dir = sys.argv[1]
    posterior_dir = os.path.join(output_dir, 'POSTERIOR')

    for fname in sorted(os.listdir(posterior_dir)):
        total = 0
        noisy = 0
        with open(os.path.join(posterior_dir, fname)) as f:
            header = f.readline()
            states = sorted(int(s[1:]) for s in f.readline().strip().split('\t'))
            for line in f:
                total += 1
                probs = map(float, line.strip().split())
                p1, p2 = sorted(probs, reverse = True)[:2]
                if p1 <= 0.5 or p2 >= 0.4:
                    noisy += 1
        print '%s\ttotal:%d\tnoisy:%d\t%.2f%%' % (fname, total, noisy, (100*float(noisy)/total))
#        elapsed(fname)

