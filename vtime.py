#!/usr/bin/python

# A script to do timing tests to compare CPU and GPU implementations.

import time
from subprocess import * # for popen, running processes
import csv
import getopt
import sys # to check for file existence, etc.

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
    
    usage = "./vtime.py" \
            " -g<0 or 1>: use GPU?\n" \
            " -n<int>: max number of subcells in each direction\n" \
            " -f<filename>: output filename\n" 
    try:
        opts, args = getopt.getopt(argv,"g:f:n:")
    except getopt.GetoptError:
        print usage
    for opt, arg in opts:
        if opt in ("-g"):
            gpu = (arg == "True" or arg == "true" or arg == "1")
        if opt in ("-f"):
            filename = arg
        if opt in ("-n"):
            nmax = int(arg)

    cmd0 = []
    cmd0.append("./testmanyv")
    cmd0.append("-t 1")
    if(gpu):
        cmd0.append("-g1")
    else:
        cmd0.append("-g0")

    create_csv("time.csv", "#n\ttime(s)\tcommand: "+str(cmd0) +"\n")


    n = 1
    while n <= nmax:
        cmd = cmd0
        cmd.append("-x"+str(n))
        cmd.append("-y"+str(n))

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
