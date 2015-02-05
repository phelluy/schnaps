#!/usr/bin/python -u
import os
import sys
import csv
import numpy as np
from subprocess import * # for popen, running processes
import re # regexp package

# usage: ./errortest.py <output filename>

filename = sys.argv[1]
progname = "testmanyv"

degmax = 3
nrafmax = 4

f = open(filename, 'wb') # erase file
f.write("#deg\tnraf\tdof\ttime(s)\titmax\ttimepRK(s)\terr\tdt\tdx\n");
f.close()

def lineafter(searchstring, output):
    #print "Searching for " + searchstring
    #print output

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

nvx = 30
nvy = 30

deg = 3
while(deg <= degmax):

    if(deg == 1):
        nrafmax = 32
    if(deg == 2):
        nrafmax = 16
    if(deg == 3):
        nrafmax = 32
    if(deg == 4):
        nrafmax = 8

    tmax = 0.4
    
    nraf = 4
    dof = 4 * (deg + 1)^2 * nraf^2 * nvx * nvy
    while(dof < 1e7):
        dof = 4 * (deg + 1) * (deg + 1) * nraf * nraf * nvx * nvy


        print "deg: " + str(deg)
        print "nraf: " + str(nraf)
        print "dof: " + str(dof)
        
        
        dt = 0.01
        cfl = 0.25

        cmd = []
        cmd.append("./" +  progname)
        cmd.append("-d " + str(deg))
        cmd.append("-n " + str(nraf))
        cmd.append("-t " + str(tmax))
#        tmax /= 2
        cmd.append("-C2")
        cmd.append("-g0")
        cmd.append("-X" + str(nvx))
        cmd.append("-Y" + str(nvy))

        print "\t" + str(cmd)
        p = Popen(cmd, stdout = PIPE, stderr = PIPE)
        p.wait()
        prc = p.returncode
        out, err = p.communicate()
            
        if (prc == 0):
            error = lineafter("L2", out)
            DOF = lineafter("DOF", out)
            deltax = lineafter("deltax", out)
            deltat = lineafter("deltat", out)
            exectime = lineafter("executiontime", out)
            perRK2time = lineafter("perRK2time", out)
            itermax = lineafter("itermax", out)

            f = open(filename, 'a') # now we append
            datawriter = csv.writer(f, delimiter = '\t')
            datawriter.writerow([deg, nraf, DOF, exectime, itermax, perRK2time, error, deltat, deltax])
            f.close()

        else: 
            print "cout:"
            print out
            print "cerr:"
            print err


        nraf *= 2
    deg += 1        

