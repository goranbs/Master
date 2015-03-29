# -*- coding: utf-8 -*-
"""
Created on Mon Mar 23 18:21:10 2015

@author: goran
"""

def readfile(fullpath):
    # Header of file looks like:
    '''
        ITEM: TIMESTEP
        350000
        ITEM: NUMBER OF ATOMS
        111
        ITEM: BOX BOUNDS xy xz yz pp pp pp
        -5.0153 10.0306 -5.0153
        0 8.6868 0
        0 34 0
        ITEM: ATOMS id type x y z vx vy vz fx fy fz 
    '''
    
    try:
        ofile = open(fullpath,'r') 
    except:
        print "Error! Cannot open file. path is incorrect!"
        print fullpath
    
    ATOMS = False
    item = None
    atoms = None
    lastline = None
    while (ATOMS is not True):
        aline = ofile.readline()
        words = aline.split()
        try:
            item = words[0]
        except:
            print "Error! The file seems to have a wrong format!"

        if (item == "ITEM:"):        
            try:
                atoms = words[1]
            except:
                continue
        
            if (atoms == "ATOMS"):
                lastline = words[2:]
                ATOMS = True
    
            

    Noutputs = len(lastline)             # number of output variables
    dictionary = {}
    for i in range(Noutputs):
        dictionary[lastline[i]] = []     # add dict list to item
    
        
    for line in ofile:
        values = line.split()
        try:
            float(values[0])
        except:
            #print "#--------------------------LAST LINE-----------------------------#"
            #print values
            break
        
        for i in range(Noutputs):
            dictionary[lastline[i]].append(float(values[i]))
        
        
    return dictionary, lastline
    
    
#def plotshit(dictionary, dictionary_keywords):
 
def calc_total_force(dictionary,dictionary_keywords,Title, ids):
    
    ID = dictionary["id"]
    forcex = dictionary["fx"]
    forcey = dictionary["fy"]
    forcez = dictionary["fz"]
    lforcex = len(forcex)

    maxforcex = max(forcex)
    maxforcey = max(forcey)
    maxforcez = max(forcez)
    maxforce = np.sqrt(maxforcex**2 + maxforcey**2 + maxforcez**2)
    threshold = 0.5
    if maxforce < threshold:
        print "\n\n\n\n\n\n\n\n\n\n\n\n !!!!!!!!!!!!!!!!!!!!!!\n maxforce < %.3f " % threshold
        maxforce = 1.0
    Nbins = 100
    forces = np.linspace(0,maxforce,Nbins)
    df = forces[1]- forces[0]
    Nforces = np.linspace(0,0,Nbins)
    allforces = []
    Fx = 0; Fy = 0; Fz = 0; F = 0; totforce_CO2 = 0; nparticles = 0
    for i in range(lforcex):
        fx = forcex[i]
        fy = forcey[i]
        fz = forcez[i]
        Fx += fx
        Fy += fy
        Fz += fz
        F_ = np.sqrt(fx**2 + fy**2 + fz**2)
        allforces.append(F_)
        F += F_
        Bin = np.floor(F_/df)
        if (Bin > Nbins):
            print "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n Bin > Nbins : (%f / %f)\n" % (Bin,Nbins)
            Bin = Nbins-1
        Nforces[Bin] += 1
        if (ID[i] in ids):
            nparticles += 1
            totforce_CO2 += F_

    # Units: REAL   
    #force = Kcal/mole-Angstrom
    kcalmol2joule = 6.9477E-21  # convertionfactor from [kcal/mol] -> [Joule]
    kcalmol2eV = 0.043          # convertionfactor from [kcal/mol] -> [eV]
    Fj = F*kcalmol2joule/lforcex
    FeV = kcalmol2eV*F/lforcex
    totforce_CO2_eV = kcalmol2eV*totforce_CO2

    for i in range(Nbins):
        Nforces[i] = float(Nforces[i])/lforcex
        forces[i] = forces[i]*kcalmol2eV
        #Nforces[i] = float(Nforces[i])
    
    plt.figure()
    plt.plot(forces,Nforces)
    plt.title(Title); plt.xlabel('Force [eV/Angstrom]'); plt.ylabel('possibility [\%]')
    
    # the histogram of the data with histtype='step'
    #plt.figure()
    #n, bins, patches = plt.hist(allforces, 50, normed=1, histtype='stepfilled')
    #plt.setp(patches, 'facecolor', 'g', 'alpha', 0.75)
    #plt.title(Title); plt.xlabel('Force [eV/Angstrom]'); plt.ylabel('possibility [\%]')
    plt.show()
    # add a line showing the expected distribution
    #y = plt.normpdf( bins, mu, sigma)
    #plt.plot(bins, y, 'k--', linewidth=1.5)
        
    print "#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#"
    print " %s " % Title
    print "#=========================================#"
    #print "#-------  Total forces in system  -------#"
    #print "# Fx = %.16f [N]" % Fx
    #print "# Fy = %.16f [N]" % Fy
    #print "# Fz = %.16f [N]" % Fz
    print "#------  Total force per particle   ------#"
    print "# F = %g [eV/Angstrom]" % (FeV)
    print "#-----------------------------------------#"
    print "# F_CO2 = %.3f [ev/Angstrom] (total force)" % totforce_CO2_eV
    print "# F_CO2 = %.3f [ev/Angstrom] (total force per particle)" % (totforce_CO2_eV/nparticles)
    print "# Number of particles counted: %g" % nparticles
    
def main():

    units = "REAL"    
    
    # Calculate total forces in the system given an input file
    #path1 = "/home/goran/lammps-28Jun14/theproject/alpha-quartz/reaxff/Alexandru_sio2_files/carbondioxide_SiO2/CO3_1/run1_0K_nvt_periodic"
    #filename1 = "dump.SiO2_periodic.350000.txt"
    path1 = "/home/goran/lammps-28Jun14/theproject/alpha-quartz/reaxff/Alexandru_sio2_files/carbondioxide_SiO2/CO3_1/minimize"
    filename1 = "dump.quartz_n_carbs.370000.txt"
    title = "CO3_1"
    ids = [109,110,111]
    id1 = ids
    fullpath1 = os.path.join(path1,filename1)
    dic1, types1 = readfile(fullpath1)
    calc_total_force(dic1,types1,title,ids)
    print "Number 1: %s done" % title

    #path1 = "/home/goran/lammps-28Jun14/theproject/alpha-quartz/reaxff/Alexandru_sio2_files/carbondioxide_SiO2/CO3_1_f/run2_0K_nvt_periodic"
    #filename1 = "dump.SiO2_CO2_periodic.350000.txt"
    path1 = "/home/goran/lammps-28Jun14/theproject/alpha-quartz/reaxff/Alexandru_sio2_files/carbondioxide_SiO2/CO3_1_f/minimize"
    filename1 = "dump.quartz_n_carbs.20000.txt"
    title = "CO3_1_f"
    ids = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,113,114,115,116]
    id2 = ids
    fullpath1 = os.path.join(path1,filename1)
    dic1, types1 = readfile(fullpath1)
    calc_total_force(dic1,types1,title,ids)
    print "Number 2: %s done" % title
    
    
    #path1 = "/home/goran/lammps-28Jun14/theproject/alpha-quartz/reaxff/Alexandru_sio2_files/carbondioxide_SiO2/CO3_2_f/run2_0K_nvt_periodic"
    #filename1 = "dump.SiO2_CO2_periodic.350000.txt"
    path1 = "/home/goran/lammps-28Jun14/theproject/alpha-quartz/reaxff/Alexandru_sio2_files/carbondioxide_SiO2/CO3_2_f/minimize"
    #filename1 = "dump.quartz_n_carbs.390000.txt"
    filename1 = "dump.quartz_n_carbs.425000.txt"
    title = "CO3_2_f"
    ids = [1,2,3,4,5,6,7,8,13,14,15,16,17,18,19,20,109,110,111,112,129,130,131,132]
    id3 = ids
    fullpath1 = os.path.join(path1,filename1)
    dic1, types1 = readfile(fullpath1)
    calc_total_force(dic1,types1, title,ids)
    print "Number 3: %s done" % title
    
    #path1 = "/home/goran/lammps-28Jun14/theproject/alpha-quartz/reaxff/Alexandru_sio2_files/carbondioxide_SiO2/SiO4_C_O/run1_0K_nvt_periodic"
    #filename1 = "dump.SiO2_periodic.350000.txt"
    path1 = "/home/goran/lammps-28Jun14/theproject/alpha-quartz/reaxff/Alexandru_sio2_files/carbondioxide_SiO2/SiO4_C_O/minimize"
    filename1 = "dump.quartz_n_carbs.370000.txt"
    title = "SiO4_C_O"
    ids = [109,110,111]
    id4 = ids
    fullpath1 = os.path.join(path1,filename1)
    dic1, types1 = readfile(fullpath1)
    calc_total_force(dic1,types1,title,ids)
    print "Number 4: %s done" % title
    
    
    path1 = "/home/goran/lammps-28Jun14/theproject/alpha-quartz/reaxff/Alexandru_sio2_files/slab/minimize"
    filename1 = "dump.quartz_n_carbs.350000.txt"
    title = "Slab"
    ids = [1]
    fullpath1 = os.path.join(path1,filename1)
    dic1, types1 = readfile(fullpath1)
    calc_total_force(dic1,types1,title,ids)
    print "Number 5: %s done" % title
    
    #print "#############################################"
    #print len(id2), len(id3), len(id1), len(id4)
        
    #plotshit(dic,types)    
    
    #print dic
    #¤print types
    
    ###################### EnD of MaiN ############################
    
if (__name__ == "__main__"):
    import os
    import matplotlib.pyplot as plt
    import numpy as np
    main()
    