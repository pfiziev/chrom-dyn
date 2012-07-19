import datetime
import sys

def error(msg):
    print >> sys.stderr, 'ERROR: %s' % msg
    exit(1)


global_stime = datetime.datetime.now()
def elapsed(msg = None):
    print "[%s]" % msg if msg is not None else "+", "Last:" , datetime.datetime.now() - elapsed.stime, '\tTotal:', datetime.datetime.now() - global_stime

    elapsed.stime = datetime.datetime.now()

elapsed.stime = datetime.datetime.now()


def open_log(fname):
    open_log.logfile = open(fname, 'w', 1)

def logm(message):
    open_log.logfile.write("[ %s ] %s\n" % (datetime.datetime.now().strftime('%Y-%m-%d T%H:%M:%S'), message))

def close_log():
    open_log.logfile.close()



BIN_SIZE = 200

def matrix(d1, d2):
    return [[0 for j in xrange(d2)] for i in xrange(d1)]



def partial_overlap(s1, e1, s2, e2):
    return s1 <= e2 and s2 <= e1 #and min(e1, e2) - max(s1, s2) >= 10
