#!/usr/bin/python -u
import os
import sys
import csv
import numpy as np
from subprocess import * # for popen, running processes
import re # regexp package

# usage: speed_test <output filename>

filename = sys.argv[1]
progname = "schnaps"

f = open(filename, 'wb') # erase file
f.write("#deg\tnraf\tdof\texectime\tperRK2time\n")
f.close()

# Return the line which appears just after the search string
def lineafter(searchstring, output):
    outlines = out.split('\n')
    dataline = ""
    itline = 0
    while itline < len(outlines):
        #print itline
        line = outlines[itline]
        #print line
        if re.search(searchstring, line) is not None:
            #print "\t"+str(outlines[itline])
            #print "\t"+str(outlines[itline + 1])
            dataline = outlines[itline + 1]
            itline = len(outlines)
        itline += 1

        if not dataline == "":
            #print dataline
            return dataline
    return ""

deg = 1
degmax = 4
while(deg <= degmax):

    if(deg == 1):
        nrafmax = 4
    if(deg == 2):
        nrafmax = 4
    if(deg == 3):
        nrafmax = 4
    if(deg == 4):
        nrafmax = 4

    tmax = 2.0
    
    m = 8
    dim = 3
    nraf = 1
    dof = (deg + 1) * (deg + 1) * (deg + 1) * nraf * nraf * nraf * m
    
    devmem = 3e9 # 3GB on the Tahiti of gpu3
    nbuffers = 3 # RK2 has wn, wnp, and dtw

    while(3 * dof * 8 < devmem):
        dof = (deg + 1) * (deg + 1) * (deg + 1) * nraf * nraf * nraf * m
        print "deg: " + str(deg)
        print "nraf: " + str(nraf)
        print "\tdof: " + str(dof)
        
        cfl = 0.25

        cmd = []

        cmd.append("./" +  progname)
        cmd.append("-n" + str(dim))
        cmd.append("-m" + str(m))
        cmd.append("-fMaxwell3DNumFluxClean_uncentered")
        cmd.append("-bMaxwell3DBoundaryFlux_uncentered")
        cmd.append("-iMaxwell3DInitData")
        cmd.append("-IMaxwell3DImposedData")
        cmd.append("-T" + str(0.1)) # tmax
        cmd.append("-w" + str(0))   # do not write output do disk
        
        cmd.append("-d" + str(deg))
        cmd.append("-r" + str(nraf))

        cmd.append("-P" + str(0))
        cmd.append("-D" + str(0))
        
        print "\t" + str(cmd)
        print "\t" + " ".join(cmd)

        p = Popen(cmd, stdout = PIPE, stderr = PIPE)
        p.wait()
        prc = p.returncode
        out, err = p.communicate()
            
        if (prc == 0):
            exectime = lineafter("Total RK time", out)
            print "\texectime:", exectime
            perRK2time = lineafter("Total RK time per time-step", out)
            print "\tperRK2time:", perRK2time

            f = open(filename, 'a') # now we append
            datawriter = csv.writer(f, delimiter = '\t')
            datawriter.writerow([deg, nraf, dof, exectime, perRK2time])
            f.close()

        else: 
            print "cout:"
            print out
            print "cerr:"
            print err

        nraf *= 2
    deg += 1        

