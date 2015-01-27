#!/usr/bin/python -u
import os
import sys
import csv
from subprocess import * # for popen, running processes
import re # regexp package

# usage: ./errortest.py <output filename>

filename = sys.argv[1]
progname = "testmanyv"

degmax = 4
nrafmax = 64

f = open(filename, 'wb') # erase file
f.write("#deg\tnraf\tDOF\terror\n");
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
            print "\t"+str(outlines[itline])
            print "\t"+str(outlines[itline + 1])
            dataline = outlines[itline + 1]
            itline = len(outlines)
        itline += 1

        if not dataline == "":
            #print dataline
            return dataline
    return ""
    
nraf = 1
while(nraf <= nrafmax):

    deg = 1
    while(deg <= degmax):

        print "deg: " + str(deg)
        print "nraf: " + str(nraf)

        dt = 0.01
        cfl = 1.0

        cmd = []
        cmd.append("./" +  progname)
        cmd.append("-d " + str(deg))
        cmd.append("-n " + str(nraf))
        cmd.append("-t " + str(0.4))
        cmd.append("-C")
        cmd.append("-X30")
        cmd.append("-Y30")
        cmd.append("-s" + str(dt))

        L2error = []
        DOF = 0

        maxtests = 4
        i = 0
        while(i < maxtests):
            #print command0 + str(cfl) + command1
            print "\tdt: " + str(dt)
                        
            cmd[len(cmd) - 1] = "-s" + str(dt)

            print cmd
            p = Popen(cmd, stdout = PIPE, stderr = PIPE)
            p.wait() # sets the return code
            prc = p.returncode
            out, err = p.communicate() # capture output
            
            if (prc == 0): # did the process succeed?
                L2error.append(lineafter("L2", out))
                DOF = lineafter("DOF", out)

            cfl *= 0.5
            if(i > 0):
                if(L2error[i] == L2error[i-1]):
                    break;

            dt *= 0.5
            i += 1

        error = L2error[len(L2error) - 1]
        print "\terror: " + str(error)

        if(True or error < 1.0):
            print "append!"
            f = open(filename, 'a') # now we append
            datawriter = csv.writer(f, delimiter = '\t')
            datawriter.writerow([deg, nraf, DOF, error])
            f.close()

        deg += 1
    nraf *= 2
        

