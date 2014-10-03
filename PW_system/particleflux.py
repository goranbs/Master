# -*- coding: utf-8 -*-
"""
Created on Thu Oct  2 10:19:05 2014

@author: goran

Python script that calculates the particle flux for some particles in a
set of statefiles for a particle system.
Reads a set of .txt -LAMMPS output files with format:

id type xs ys zs

The script should for som particle types check if the particle has crossed
some border in a particlular direction, and count the number of particles
that passes either way from one timestep to another (e.g from one statefile 
to the next)
"""

import sys
import os
import numpy as np

def readfile(filename, readlines):
    
    ofile = open(filename, 'r')
    for i in range(readlines):
        line = ofile.readline()

    aid = []
    atype = []
    xs = []
    ys = []
    zs = []
    
    for line in ofile:
        ID,TYPE,XS,YS,ZS = line.split()
        aid.append(int(ID)),atype.append(int(TYPE)),xs.append(float(XS)),ys.append(float(YS)),zs.append(float(ZS))
       
    readfilearray = [aid,atype,xs,ys,zs]
    return readfilearray
    

def compare(time1,time2,types,x):
    """
    Takes two timesteps of the system, where time1,time2 contains
    information about the system at timestep t and t+1.
    The number of particles crossing the boarder in a given direction has
    to be recorded and retured from the function.
    """
    
    lentime1 = len(time1[0]); lentime2 = len(time2[0])
    ltime1 = len(time1); ltime2 = len(time2)
    
    if (lentime1 != lentime2):
        print "ERROR. The lenght of statefiles are different. lentime"
        sys.exit()
    if(ltime1 != ltime2):
        print "ERROR! The number of indices of the state arrays are different!"
        print "len(time1)=%g len(time2)=%g" % (ltime1,ltime2)

    ntypes = 0    # number of particles we track
    nsearch = 0   # number of times that we have to search for the particle
    pos_cross = 0; neg_cross = 0; no_cross = 0
    for i in range(lentime1):
        if (time1[1][i] in types):
            ntypes += 1   # count the number of particles we track. 
            #x1 = time1[2] # get the x-positions at time t
            x1 = 0.001
            if (time2[0][i] == time1[0][i]):
                # find if atom has crossed line:
                x2 = time1[2][i]
                if (x1 < x and x < x2):
                    pos_cross += 1 # particle has crossed the line in pos x dir
                elif (x1 > x and x > x2):
                    neg_cross += 1 # particle has crossed the line in neg x dir
                else:
                    no_cross += 1 # the particle did not cross the line
            
            else:
                index = None
                nsearch += 1
#               print "the particle index is not at the same line! Searching..."
                for j in range(lentime2):
                    if (time2[0][j] == time1[0][i]):
                        index = j
#                        print "found it at line %g, index %g == %g" % (j,time2[0][j],time1[0][i])
                if (index is not None):
                    # find if atom has crossed line:
                    x2 = time2[2][j]
                    if (x1 < x and x < x2):
                        pos_cross += 1 # particle has crossed the line in pos x dir
                    elif (x1 > x and x > x2):
                        neg_cross += 1 # particle has crossed the line in neg x dir
                    else:
                        no_cross += 1 # the particle did not cross the line
                        
            
            
    #print "particles: %g == %g ?" % (ntypes, (no_cross+pos_cross+neg_cross))
    #print "number of times have to search for the correct atom: %g " % nsearch
    crossings = pos_cross - neg_cross
    ncrossings = [crossings, pos_cross, neg_cross, no_cross, ntypes, nsearch]
    return ncrossings
    
def gothroughfiles(path):
    filenames = []
    numbers = []
    for name in os.listdir(path):
        if (name[-4:] == ".txt"):
            filenames.append(name)
            numbers.append(int(name[-10:-4]))

    time = sorted(numbers)
    sortedfilenames = sorted(filenames, key = lambda x: x[:-4])
    
    return sortedfilenames,time
    
    
def plotting(compares,time,types):

    import matplotlib.pyplot as plt
    plt.rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
    ## for Palatino and other serif fonts use:
    #rc('font',**{'family':'serif','serif':['Palatino']})
    plt.rc('text', usetex=True)
    #print " len(time)=%g \n len(compares)=%g " % (len(time),len(compares))
    
    mean_cross = np.mean(compares[:,0])
    std_cross = np.std(compares[:,0])
    mean_pos = np.mean(compares[:,1])
    mean_neg = np.mean(compares[:,2])
    mean_search = np.mean(compares[:,5])
    std_search = np.std(compares[:,5])
    mean_nocross = np.mean(compares[:,3])
    std_nocross = np.std(compares[:,3])
    
    print "-------------------------"
    print "particle types:"
    for i in range(len(types)):
        print "                %g" % types[i]
    
    plt.figure()
    plt.plot(time[1:],compares[:,0],'b-')
    plt.hold(True)
    plt.plot(time[1:],compares[:,1],'g--')
    plt.plot(time[1:],compares[:,2],'y--')
    plt.hold(False)
    plt.title(r'Number of crossings as function of time. $ N_{cross} = $ %.f $ \pm $ %.f' % (mean_cross, std_cross))
    plt.xlabel('Time [fsec]')
    plt.ylabel('crossings')
    plt.legend([r'$n(t)$',r'$n_{pos}(t)$',r'$n_{neg}(t)$'],loc='upper left')
    
    plt.figure()
    plt.plot(time[1:],compares[:,3],'r-*')
    plt.title(r'Number of particles that did not cross. $ N_{nocross} = $ %.f $ \pm $ %.f' % (mean_nocross, std_nocross))
    plt.xlabel('Time [fsec]')
    plt.ylabel('Particles')
    plt.legend([r'$n(t)$'],loc='lower left')
    
    plt.figure()
    plt.plot(time[1:],compares[:,5],'r-*')
    plt.title('Number of times we had to search for the comparison particle.$ N_{search} = $ %.f $ \pm $ %.f' % (mean_search, std_search))
    plt.xlabel('Timestep file')
    plt.ylabel('Number of searches')
    plt.legend([r'$N(file)$'],loc='upper left')
    plt.show(True)
   
def testtwofiles():
    path = "/home/goran/Goran/PythonScripting/statefiles/"
    filename1 = "dump.PW_flow.200700.txt" 
    filename2 = "dump.PW_flow.200800.txt"
    file1 = os.path.join(path, filename1)
    file2 = os.path.join(path, filename2)
    
    types = [4] # just track the oxygen atoms in water
    x = 0.5    
    readlines = 9
    
    print "running testtwofiles ..."
    n = compare(readfile(file1,readlines),readfile(file2,readlines),types,x)
    print n

def testplot():
    # run the function plotting for test.

    time = np.linspace(0,10,100)
    compares = []
    types = [4]
    for i in range(len(time)-1):
        j = i+1
        pos = 0.1*j
        neg = (j - 0.1*j)
        compares.append([pos - neg ,pos, neg,j-(pos-neg),j,(0.1*j*j - 0.2*j*(j-0.2))])
    compares = np.array(compares)
    plotting(compares,time,types)

def main():
    
    path = "/home/goran/Goran/PythonScripting/statefiles/"

    filenames,time= gothroughfiles(path)
    compares = []
    types = [4]    # only looking for oxygen in water
    x = 0.5       # passing through boarder at 0.5*(xhi - xlo)
    readlines = 9  # number of lines in beginning of files that should be neglected 
    iterations = (len(filenames) - 1)

    run = "test"
    #run = "run"
    if (run == "test"):
        testplot()
        testtwofiles()
        
    else: 
        for i in range(iterations):
            filename1 = filenames[i]
            filename2 = filenames[i+1]
            file1 = os.path.join(path,filename1)
            file2 = os.path.join(path,filename2)
            print " (%g / %g) Comparing files %s and %s " % (i+1,iterations,filename1,filename2)
            n = compare(readfile(file1,readlines),readfile(file2,readlines),types,x)
            compares.append(n)

        compares = np.array(compares)
        plotting(compares,time,types)
        print "tot crossings, pos, neg, no_cross, tot particles, n_searches"
        print compares

if (__name__ == "__main__"):
    main()
    