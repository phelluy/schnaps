#!/usr/bin/python

# A script to do timing tests to compare CPU and GPU implementations.

import time
from subprocess import * # for popen, running processes
import csv
import getopt
import sys # to check for file existence, etc.
import copy

# Create an empty tsv file with comments
def create_csv(filename, comments):
    f = open(filename, 'w')
    f.write(comments)
    f.close()

# Add the row (consisting of data) to the file with name filename.
def append_to_csv(filename, data):
    fd = open(filename, 'a')
    a = csv.writer(fd, delimiter='\t');
    a.writerows([data])
    fd.close()

def main(argv):
    filename = "time.csv"
    gpu = False
    nmax = 64
    nmin = 1
    nplat = 0
    ndef = 0
    do_append = False

    usage = "./vtime.py" \
            " -a: append to output instead of overwriting it\n" \
            " -g<0 or 1>: use GPU?\n" \
            " -P<int>: Platform number\n" \
            " -D<int>: Device\n" \
            " -n<int>: max number of subcells in each direction\n" \
            " -m<int>: min number of subcells in each direction\n" \
            " -f<filename>: output filename\n" 
    try:
        opts, args = getopt.getopt(argv,"ag:f:n:m:P:D:")
    except getopt.GetoptError:
        print usage
    for opt, arg in opts:
        if opt in ("-a"):
            do_append = True
        if opt in ("-g"):
            gpu = (arg == "True" or arg == "true" or arg == "1")
        if opt in ("-f"):
            filename = arg
        if opt in ("-P"):
            nplat = int(arg)
        if opt in ("-D"):
            ndev = int(arg)
        if opt in ("-n"):
            nmax = int(arg)
        if opt in ("-m"):
            nmin = int(arg)

    cmd0 = []
    cmd0.append("./testmanyv")
    cmd0.append("-t 1")
    if(gpu):
        cmd0.append("-g1")
        cmd0.append("-P" + str(nplat))
        cmd0.append("-D" + str(ndev))
    else:
        cmd0.append("-g0")

    if not do_append:
        create_csv(filename, "#n\ttime(s)\tcommand: "+str(cmd0) +"\n")

    n = nmin
    while n <= nmax:
        cmd = copy.deepcopy(cmd0)
        cmd.append("-x" + str(n))
        cmd.append("-y" + str(n))

        ts = time.time()
        p = Popen(cmd, stdout = PIPE, stderr = PIPE)
        p.wait()
        te = time.time()
        elapsed_time = te - ts
        print "n:" + str(n) + ", elapsed time: " + str(elapsed_time)
        append_to_csv(filename, [n, elapsed_time])

        n *= 2

if __name__ == "__main__":
    main(sys.argv[1:])
