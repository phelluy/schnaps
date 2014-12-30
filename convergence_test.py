#!/usr/bin/python -u
import os
import sys
import csv

# usage: ./errortest.py <output filename>

filename = sys.argv[1]
progname = "testmanyv"

degmax = 4
nrafmax = 64

f = open(filename, 'wb') # erase file
f.close()

nraf = 1
while(nraf <= nrafmax):

    deg = 1
    while(deg <= degmax):

        print "deg: " + str(deg)
        print "nraf: " + str(nraf)

        cfl = 1.0

        command0 = " ./" +  progname + " -c "
        command1 = " -d " + str(deg) + " -n " + str(nraf) + " -t 0.5 | grep L2 -A 1 | tail -n 1"

        L2error = []

        maxtests = 4
        i = 0
        while(i < maxtests):
            #print command0 + str(cfl) + command1
            print "\tcfl: " + str(cfl)
            L2error.append(float(os.popen(command0 + str(cfl) + command1).read()))
            cfl *= 0.5
            if(i > 0):
                if(L2error[i] == L2error[i-1]):
                    break;
                    
            i += 1

        error = L2error[len(L2error) - 1]
        print "\t" + str(error)

        if(error < 1.0):
            f = open(filename, 'a') # now we append
            datawriter = csv.writer(f, delimiter = '\t')
            datawriter.writerow([deg, nraf, error])
            f.close()

        deg += 1
    nraf *= 2
        

