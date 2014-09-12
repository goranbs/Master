# -*- coding: utf-8 -*-
"""
Created on Fri Sep 12 14:13:43 2014

@author: goran

This python script should read in an system-output file for a given atom
structure. Read the atoms positions, and find out if atoms overlap.
If some atoms overlap, then the system is buildt wrong, and it is nice to know
atoms that we should be looking at!
The file should look like:

N_atoms
timestep
Atom-ID   x1   y1   z1
Atom-ID   x2   y2   z2
   .      .    .    .
   .      .    .    .
   .      .    .    .
"""

import numpy as np

import sys, getopt

def main(argv):
   inputfile = ''
   try:
      opts, args = getopt.getopt(argv,"h")
   except getopt.GetoptError:
      print 'ERROR!!!'
      print 'Useage: python find_overlaooing_atoms.py <inputfile> '
      print '-----------------------------------------------------'
      print 'Where <inputfile> is the full path of the inputfile'
      print '-----------------------------------------------------'
      print ' The file is asumed to look like: '
      print ' N_atoms \n timestep \n Atom-ID   x1   y1   z1 \n Atom-ID   x2   y2   z2 \n    .      .    .    . \n    .      .    .    . \n    .      .    .    . \n '
      sys.exit(2)
   for opt, arg in opts:
      if opt == '-h':
         print 'Useage: python find_overlapping_atoms.py <inputfile>'
         print '-----------------------------------------------------'
         print 'Where <inputfile> is the full path of the inputfile'
         print '-----------------------------------------------------'
         print ' The file is asumed to look like: '
         print ' N_atoms \n timestep \n Atom-ID   x1   y1   z1 \n Atom-ID   x2   y2   z2 \n    .      .    .    . \n    .      .    .    . \n    .      .    .    . \n '
      
         sys.exit()
   if (len(args) == 1):
        inputfile = args[0]
   else:
       print 'ERROR!'
       print 'Useage: python find_overlapping_atoms.py <inputfile>'
       print '-----------------------------------------------------'
       print 'Where <inputfile> is the full path of the inputfile'
       print '-----------------------------------------------------'
       print ' The file is asumed to look like: '
       print ' N_atoms \n timestep \n Atom-ID   x1   y1   z1 \n Atom-ID   x2   y2   z2 \n    .      .    .    . \n    .      .    .    . \n    .      .    .    . \n '
      
       sys.exit()
        
   print 'Input file is located at: \n ', inputfile , '-----------------------------------------------------------'
   
   read_inputfile(inputfile)


def read_inputfile(inputfile):
    try:
        # is the file possible to locate/open?
        infile = open(inputfile, "r") 
    except (OSError, IOError) as e:
        print 'IOError!'
        print 'input file or outputfile or both is not readable/reachable'
    
    line1 = infile.readline()
    line2 = infile.readline()

    atoms = []
    linenumber = 0
    for line in infile:
        l = line.split()
        atom = []
        for ii in range(1,len(l)):
            atom.append(float(l[ii]))
    
        atoms.append(atom)
        linenumber += 1
    
    i = 0
    overlapping = []    
    for atom_i in atoms:
        # go through atoms to find if some atoms overlap
        j = 0
        for atom_j in atoms:
            if (i != j and i > j):
                # compare position of atom_i with atom_j
                test = "OVERLAP"
                for index in range(3):
                    if (atom_i[index] != atom_j[index]):
                        test = "OK"
                if (test != "OK"):
                    overlap = [i+1,j+1,atom_i, atom_j]
                    overlapping.append(overlap)                     
                            
            j += 1
        i += 1
                            
    if len(overlapping) > 0 :
        print "------------------------------------------------------------------------"    
        print "--------------------------------OVERLAPPS-------------------------------"
        for overlap in overlapping:
            print " atom %g [%.3f, %.3f, %.3f] overlaps with \n atom %g [%.3f, %.3f, %.3f] " % (overlap[0], overlap[2][0], overlap[2][1], overlap[2][2], overlap[1], overlap[3][0], overlap[3][1], overlap[3][2])
            print "------------------------------------------------------------------------"
        
    else:
        print "No overlapping atoms found" 



if __name__ == "__main__":
   main(sys.argv[1:])