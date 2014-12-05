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
    lotext   = 'ITEM: ATOMS id type xs ys zs charge\n'
    ofile = open(filename,'w')    
        
    ofile.write(toptext)
    ofile.write(nexttext)
    ofile.write(xlohi)
    ofile.write(ylohi)
    ofile.write(zlohi)
    ofile.write(lotext)
    print "################################"
    print "line of nearest atoms to atom at line k=%g" % k
    print nearest
    print "################################"
    for i in range(Natoms):
        neighbour = 0
        if (i in nearest):
            neighbour = 1
            print matrix['id'][i]
        if (i == k):
            neighbour = 2
            print matrix['id'][i],k
        line = '%g %g %g %g %g %g \n' % (matrix['id'][i],matrix['type'][i],matrix['xs'][i],matrix['ys'][i],matrix['zs'][i],neighbour)
        ofile.write(line)
    
    ofile.close()
    print "write to file done :-)"


def fourNearest(Natoms,matrix,Type,output_filename,end_title_of_plot,showplot=True,saveplot=False):
    '''
    Go through matrix and find the four nearest atoms of type Type. For every atom 
    of type Type, calculate the tetrehedral order parameter...
    Find the four nearest atoms of a 
    '''
    Noxygen = 0   # number of oxygen atoms in system. (# atoms of type Type)
    Qk = 0        # average Tetrahedral order parameter 
    Qk_sq = 0     # avarage of squared Tetrahedral order parameter
    counter = 0   # number of estimates on the Tetrahedral order parameter
    Qk_ = []      # tetrahedral order parameter for oxygen atom k.
        
    for k in range(Natoms):             # leap over all atoms in system
        atomtype = matrix['type'][k]    
        if atomtype == Type:            # test if it is the correct atom type Type
            distancetoeverybody = []    # list that contain k's distance to atom i
            atomindex_i = []            # holds line index of atoms that is of type Type, and i != k.
            Noxygen += 1
            k_id = matrix['id'][k]      # atom id number of atom on line k
            x = matrix['xs'][k]         # positions
            y = matrix['ys'][k]
            z = matrix['zs'][k]
            #r = np.sqrt(x*x + y*y + z*z)
            for i in range(Natoms):
                atomtype_i = matrix['type'][i]
                if (atomtype_i == Type):
                    i_id = matrix['id'][i]      # id of atom i
                    if (k_id != i_id):          # cannot take the distance between an atom and itself!
                    #if (k != i):               # would be the same as the above statement
                        xi = matrix['xs'][i]
                        yi = matrix['ys'][i]
                        zi = matrix['zs'][i]
                        #ri = np.sqrt(xi*xi + yi*yi + zi*zi)
                        distance = np.sqrt((x-xi)**2 + (y-yi)**2 + (z-zi)**2 )
                        distancetoeverybody.append(distance) # append distance between atom k and i
                        atomindex_i.append(i)                # and what number is atom i
            
            ##################################################################
            fourNearest = []
            min1 = np.argmin(distancetoeverybody)  # index of the smallest value in distancetoeverybody 
            #dist1 = distancetoeverybody[min1]     # watch over this value, if we need it...
            distancetoeverybody[min1] = 299793458  # set the value in the list to something absurdly large
            fourNearest.append(atomindex_i[min1])               # index of nearest neighbour in distance to everybody
            
            min2 = np.argmin(distancetoeverybody)
            #dist2 = distancetoeverybody[min2]            
            distancetoeverybody[min2] = 299793458            
            fourNearest.append(atomindex_i[min2])
            
            min3 = np.argmin(distancetoeverybody)
            #dist3 = distancetoeverybody[min3]
            distancetoeverybody[min3] = 299793458
            fourNearest.append(atomindex_i[min3])
            
            min4 = np.argmin(distancetoeverybody)
            #dist4 = distancetoeverybody[min4]
            fourNearest.append(atomindex_i[min4])
                 
            ###################################################################
            # --------  Calculate tetrahedral order parameter  -------------- #
            
            COS = 0
            avg_theta = 0
            for i in range(3):
                # vector of particle j:
                xi = matrix['xs'][fourNearest[i]] - x
                yi = matrix['ys'][fourNearest[i]] - y
                zi = matrix['zs'][fourNearest[i]] - z
                Ri = np.sqrt(xi*xi + yi*yi + zi*zi)
                j = i+1
                for l in range(j,4,1):
                    # print i,l
                    # vector of particle ll
                    xj = matrix['xs'][fourNearest[l]] - x
                    yj = matrix['ys'][fourNearest[l]] - y
                    zj = matrix['zs'][fourNearest[l]] - z
                    Rj = np.sqrt(xj*xj + yj*yj + zj*zj)
                    hosl = (xj*xi + yj*yi + zj*zi)
                    hyp = Rj*Ri
                    #print hosl, hyp
                    theta_ikj = np.arccos(hosl/hyp)
                    avg_theta += theta_ikj
                    COS += ( np.cos(theta_ikj) + 1/3.0 )**2
                    #print "i=%g j=%g :: COS = %g" % (i,j,3*COS/8.0)
            
            #print (3/8.0)*COS
            avg_theta = avg_theta/6.0
            #print avg_theta, avg_theta*(180/3.14159)  # average internal angles
            Q = 1 - (3/8.0)*COS  # The calculated Tetrahedral order param for this k
            Qk += Q              # Sum of order parameters in system
            Qk_sq += Q*Q         # square sum of order parameters in system 
            Qk_.append(Q)        # calculated order parameters
            counter += 1         # number of estimates
            
            '''
            ###################################################################
            # Test if we really find the closest atoms!
            if ((k_id == 2455)):
                # tested: 1687, 1912, 2455(lies at the edge; does not take care of periodic bounary conditions!)  
                writeLAMMPSfile(Natoms,matrix,fourNearest,k)
            ###################################################################
            '''
            
    Qk = Qk/counter                  # mean Qk ; <Qk> (expectationvalue)
    Qk_sq = Qk_sq/counter            # <(Qk)^2>
    std_Qk = np.sqrt(Qk_sq - Qk*Qk)  # standard deviation
    print "#################################################"
    print "# Qk = %g \pm %g " % (Qk,std_Qk)
    print "#################################################"
    
    minQk = min(Qk_)
    maxQk = 1.0
    
    nbins = 100
    Qk_values = np.linspace(minQk,maxQk,nbins)
    dQk = Qk_values[1] - Qk_values[0]             # should be the same as (max-min)/(nbins)
    #print dQk, ((maxQk - minQk)/nbins)
    # Now in what bin should we place the different Qk values?
    Qk_bins = np.zeros((nbins,1))
    hist, bin_edges = np.histogram(Qk_,nbins,density=True)
    dx = bin_edges[1] - bin_edges[0]
    
    
    rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
    ## for Palatino and other serif fonts use:
    #rc('font',**{'family':'serif','serif':['Palatino']})
    rc('text', usetex=True) 
    
    center = np.linspace(minQk+dx,maxQk-dx, nbins)
    plt.hold(True)
    line = np.linspace(0,max(hist),3)
    Qk_avg = np.linspace(Qk,Qk,3)
    plt.plot(center,hist,'b-',linewidth='3')
    plt.plot(Qk_avg,line,'y--',linewidth='5')
    plt.hist(Qk_,nbins,normed=True,color='red')
    Title = r'Tetrahedral order parameter $ Q_k = %.2f \pm %.3f $ %s' % (Qk,std_Qk,end_title_of_plot)
    plt.title(Title)
    plt.xlabel(r'$Q_k$'), plt.ylabel('Relative occurence $P(Q_k) [\%]$')
    plt.hold(False)
    if (saveplot):
        fig_name = output_filename
        plt.savefig(fig_name,format='png')

    plt.show(showplot)
    plt.close('all')
    

def main():
    '''
    Thetrahedral Order Parameter calculation for a statefile produced with LAMMPS.
    The state file has to be consisting of water only. The program will search for 
    atoms of type = Type only, and calculate the orderparameter according to the four
    nearest atoms of type = Type.
    TetrahedralParameter is dependent on os,numpy and displacement (made by goran)
    '''
    
    #path = "/home/goran/lammps-28Jun14/examples/water/water_from_298K/dump" 
    #filename = "water.210000.txt"
    #path = "/home/goran/lammps-28Jun14/examples/water_portlandite_system/water/statefiles"
    #filename = "dump.bulkwater_nvt.2300000.txt"
    path = "/home/goran/lammps-28Jun14/examples/water_portlandite_system/water/dump_pure_water"
    filename = "dump.water.0.txt"
    Type = 1  # oxygen
    showplot = True
    saveplot = True
    output_filename = 'TetrahedralOrderParam_initial_stat.png'
    end_title_of_plot   = 'for initial state water'
    
    filetoread = os.path.join(path,filename)
    #print filetoread
    
    t, Natoms, system_size, matrix, readstructure, entries,types = readfile(filetoread)   
    fourNearest(Natoms,matrix,Type,output_filename,end_title_of_plot,showplot,saveplot)
    
    
    

if __name__ == "__main__":
    import os
    import numpy as np
    from displacement import readfile
    import matplotlib.pyplot as plt
    from matplotlib import rc
    main()