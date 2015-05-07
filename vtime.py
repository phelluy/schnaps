#!/usr/bin/python

# A script to do timing tests to compare CPU and GPU implementations.

import time
from subprocess import * # for popen, running processes
import csv

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

create_csv("cputime.csv", "#n\ttime(s)")
create_csv("gputime.csv", "#n\ttime(s)")

cmd0 = []
cmd0.append("./testmanyv")
cmd0.append("-t 1")
nmax = 64

n = 1
while n <= nmax:

    cmd = cmd0
    cmd.append("-x"+str(n))
    cmd.append("-y"+str(n))
    cmd.append("-g0")

    ts = time.time()
    p = Popen(cmd, stdout = PIPE, stderr = PIPE)
    p.wait()
    te = time.time()
    elapsed_time = te - ts
    print "cpu: n:" + str(n) + ", elapsed time: " + str(elapsed_time)
    append_to_csv("cputime.csv", [n, elapsed_time])

    cmd = cmd0
    cmd.append("-x"+str(n))
    cmd.append("-y"+str(n))
    cmd.append("-g1")

    ts = time.time()
    p = Popen(cmd, stdout = PIPE, stderr = PIPE)
    p.wait()
    te = time.time()
    elapsed_time = te - ts
    print "gpu: n:" + str(n) + ", elapsed time: " + str(elapsed_time)
    append_to_csv("gputime.csv", [n, elapsed_time])

    n *= 2
