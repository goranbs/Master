# -*- coding: utf-8 -*-
"""
Created on Tue Nov 11 14:56:47 2014

@author: goran

Thetrahedral structure parameter is a way to describe how tetrahedral the structure
of the water is, and therefore a way to describe how close the water is to be in
the aggregate state of ice.

Qk = 1 - 3/8 (3,i)sum((4,j=i+1)sum(cos(theta_{ikj}) + 1/3)^2)

where theta_{ikj} is the angle i-k-j (k is the vertex angle), and the two sums goes
over the 6 possible angles theta, between the main point and its four nearest neighbours.

# We have to run through all the oxygen atoms in the water, and for each, we have to find
  its four nearest neighbours.
  
# Then we might calculate the order parameter for this molecule.
"""

def writeLAMMPSfile(Natoms,matrix,nearest,k):
    filename = 'testfile_LAMMPS.txt'
    toptext  = 'ITEM: TIMESTEP \n2500 \nITEM: NUMBER OF ATOMS\n%g\n' % Natoms
    nexttext = 'ITEM: BOX BOUNDS pp pp pp\n'
    xlohi    = '-3.32356 34.3236 \n'
    ylohi    = '-3.32356 34.3236 \n'
    zlohi    = '-3.32356 34.3236 \n'
    lotext   = 'ITEM: ATOMS id type xs yx zs charge\n'
    ofile = open(filename,'w')    
        
    ofile.write(toptext)
    ofile.write(nexttext)
    ofile.write(xlohi)
    ofile.write(ylohi)
    ofile.write(zlohi)
    ofile.write(lotext)
    for i in range(Natoms):
        neighbour = 0
        if (i in nearest):
            neighbour = 1
        if (i == k):
            neighbour = 2
        line = '%g %g %g %g %g %g \n' % (matrix['id'][i],matrix['type'][i],matrix['xs'][i],matrix['ys'][i],matrix['zs'][i],neighbour)
        ofile.write(line)
    
    
    ofile.close()
    print "done :-)"

def fourNearest(Natoms,matrix,Type):
    '''
    Go through matrix and find the four nearest atoms of type Type. For every atom 
    of type Type, calculate the tetrehedral order parameter...
    Find the four nearest atoms of a 
    '''
    Noxygen = 0
    for k in range(Natoms):             # find number of water oxygen molecules there are...
        atomtype = matrix['type'][k]    # test if it is the correct atom type
        if atomtype == Type:
            distancetoeverybody = []    # list that contain k's distance to atom i
            atomindex_i = []
            Noxygen += 1
            k_id = matrix['id'][k]
            x = matrix['xs'][k]
            y = matrix['ys'][k]
            z = matrix['zs'][k]
            #r = np.sqrt(x*x + y*y + z*z)
            for i in range(Natoms):
                atomtype_i = matrix['type'][i]
                if (atomtype_i == Type):
                    i_id = matrix['id'][i]
                    if (k_id != i_id): 
                        xi = matrix['xs'][i]
                        yi = matrix['ys'][i]
                        zi = matrix['zs'][i]
                        #ri = np.sqrt(xi*xi + yi*yi + zi*zi)
                        
                        distance = np.sqrt((x-xi)**2 + (y-yi)**2 + (z-zi)**2 )
                        distancetoeverybody.append(distance)
                        atomindex_i.append(i)
                        
                        
            fourNearest = np.zeros((4,1))
            min1 = np.argmin(distancetoeverybody) # index of the smallest value in distancetoeverybody           
            dist1 = distancetoeverybody.pop(min1) # remove this value, and store it in dist1
            fourNearest[0] = min1
            
            min2 = np.argmin(distancetoeverybody)
            dist2 = distancetoeverybody.pop(min2)
            fourNearest[1] = min2
            
            min3 = np.argmin(distancetoeverybody)
            dist3 = distancetoeverybody.pop(min3)
            fourNearest[2] = min3
            
            min4 = np.argmin(distancetoeverybody)
            dist4 = distancetoeverybody.pop(min4)
            fourNearest[3] = min4
            
            if (k == 2970):
                writeLAMMPSfile(Natoms,matrix,fourNearest,k)
            
 


def main():
    '''
    Thetrahedral Order Parameter calculation for a statefile produced with LAMMPS.
    TetrahedralParameter is dependent on os,numpy and displacement (made by goran)
    '''
    
    path = "/home/goran/lammps-28Jun14/examples/water/water_from_298K/dump"    
    filename = "water.210000.txt"
    Type = 2  # oxygen
    
    filetoread = os.path.join(path,filename)
    print filetoread
    
    t, Natoms, system_size, matrix, readstructure, entries,types = readfile(filetoread)   
    fourNearest(Natoms,matrix,Type)
    
    
    
    
    
if __name__ == "__main__":
    import os
    import numpy as np
    from displacement import readfile
    main()