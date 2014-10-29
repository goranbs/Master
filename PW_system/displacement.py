# -*- coding: utf-8 -*-
"""
Created on Fri Oct 24 09:57:09 2014

@author: goran

Python sctipt that should take a set of LAMMPS output files and find the 
displacement of the particles in the system as a function of distance to
a known surface.
The script should take the first file and locate the initial positions of the 
interesting particles (provided by user) and store them according to a 1D
system volume grid.
The ID's of the particles are then found in the next statefile, and the displacement
of the particles are calculated as a funciton of timestep.

"""

def writetofile(filename,toptext,contents,column1=None,column2=None,column3=None,column4=None):
    '''
    write information to output file "filename". 
    '''
    generate = True
    if (type(toptext)!= str):
        print "ERROR! type(toptext) = %s " % type(toptext)
        print "when type(toptext) should be <str>"
        print "Outputfile is not generated"
        generate = False
    if (type(contents)!= str):
        print "ERROR! type(contents) = %s " % type(contents)
        print "when type(contents) should be <str>"
        print " outputfile is not created"
        generate = False
    if (column1 == None):
        print "failed to give datas in column1"
        print "outputfile is not created"
        generate = False
        
    if (generate == True):
        write_to_file = open(filename,'w')
        toptext = '# ' + toptext + '\n'
        contents = '# ' + contents + '\n'
        write_to_file.write(toptext)
        write_to_file.write(contents)
        kk = 1
        a = column1
        b = c = d = []
        string = 0 #
        for j in range(len(column1)):
            b.append(string)
            c.append(string)
            d.append(string)
        if column2 is not None: 
            kk +=1
            b = column2
        if column3 is not None: 
            kk +=1
            c = column3
        if column4 is not None: 
            kk +=1
            d = column4

            
        for i in range(len(column1)):
            #print i, len(column1), len(column2)
            line = "%g   %g   %g   %g \n" % (a[i], b[i], c[i], d[i])
            write_to_file.write(line)
    
        
        write_to_file.close()
        print "file %s is done." % filename
        
        
def gothroughfiles(path,arg=None):
    '''
    Go through stystem state files in "path", if they are .txt files, the first 
    character in the filename which is a number will be acknowledged as the 
    timestep of the system state file, and the number will be read as:
    somefilename<somenumber>.txt , when <somenumber> is the time that will be
    read. 
    Returns:
            sortedfilenames, time
    '''
    
    filenames = []
    numbers = []
    for name in os.listdir(path):
        if arg is not None:
            if (arg in name):
                if (name[-4:] == ".txt"):
                    char = 'string'
                    start = 0
                    for character in name:
                        # find the first character in the filename which is a number
                        try:
                            char = int(character)
                        except:
                            pass
                        if (type(char) == int):
                            break
                        start += 1
                
                    filenames.append(name)
                    numbers.append(int(name[start:-4]))
        else:
            if (name[-4:] == ".txt"):
                    char = 'string'
                    start = 0
                    for character in name:
                        # find the first character in the filename which is a number
                        try:
                            char = int(character)
                        except:
                            pass
                        if (type(char) == int):
                            break
                        start += 1
                
                    filenames.append(name)
                    numbers.append(int(name[start:-4]))
                    
    time = sorted(numbers)
    sortedfilenames = sorted(filenames, key = lambda x: x[:-4])
    
    return sortedfilenames,time

def readfile(filename):
    '''
    Returns timestep, number of atoms, system sizes and
    a matrix that contains information about the atoms\n
    Returns:\n
    t, Natoms, system_size, matrix, readstructure, entries, types
    '''
    infile = open(filename,'r')
    infile.readline()
    t = int(infile.readline()) # second line is the timestep
    infile.readline()
    Natoms = int(infile.readline()) # fourth line is number of atoms
    infile.readline()
    line1 = infile.readline() # system size x-axis
    line2 = infile.readline()
    line3 = infile.readline()
    
    xlo, xhi = line1.split()
    ylo, yhi = line2.split() # system size y-axis
    zlo, zhi = line3.split() # system size z-axis
    
    system_size = [float(xlo),float(xhi),float(ylo),float(yhi),float(zlo),float(zhi)]
    #structure of file: 
    
    line1 = infile.readline()
    columns = line1.split()
    # the two first words are ITEM: ATOMS
    c_id=c_mol=c_type=c_xs=c_ys=c_zs=c_vx=c_vy=c_vz=c_fx=c_fy=c_fz= -1
    ID=MOL=TYPE=XS=YS=ZS=VX=VY=VZ=FX=FY=FZ= None
    col = []
    entries = []
    it = -2
    for column in columns:

        if column == "id":
            c_id = it
            ID = []
            entries.append('id')
            col.append(c_id)
        if column == "mol":
            c_mol = it
            MOL = []
            entries.append('mol')
            col.append(c_mol)
        if column == "type":
            c_type = it
            TYPE = []
            entries.append('type')
            col.append(c_type)
        if column == "xs":
            c_xs = it
            XS = []
            entries.append('xs')
            col.append(c_xs)
        if column == "ys":
            c_ys = it
            YS = []
            entries.append('ys')
            col.append(c_ys)
        if column == "zs":
            c_zs= it
            ZS = []
            entries.append('zs')
            col.append(c_zs)
        if column == "vx":
            c_vx = it
            VX = []
            entries.append('vx')
            col.append(c_vx)
        if column == "vy":
            c_vy = it
            VY = []
            entries.append('vy')
            col.append(c_vy)
        if column == "vz":
            c_vz = it
            VZ = []
            entries.append('vz')
            col.append(c_vz)
        if column == "fx":
            c_fx = it
            FX = []
            entries.append('fx')
            col.append(c_fx)
        if column == "fy":
            c_fy = it
            FY = []
            entries.append('fy')
            col.append(c_fy)
        if column == "fz":
            c_fz = it
            FZ = []
            entries.append('fz')
            col.append(c_fz)
            
        it += 1
    
    readstructure = []
    for index in col:
        if (index >= 0):
            readstructure.append(index)
  
    matrix = {'id':ID,'mol':MOL,'type':TYPE,'xs':XS,'ys':YS,'zs':ZS,'vx':VX,'vy':VY,'vz':VZ,'fx':FX,'fy':FY,'fz':FZ}
    types = [] # list of atom types

    indexes = len(readstructure)
    for i in range(Natoms):
        line = infile.readline()
        values = line.split()
        for index in readstructure:
            matrix[entries[index]].append(values[index])
    
    for i in range(Natoms):
        for j in range(indexes):
            if matrix[entries[j]] is not None:
                if (entries[j] in ['id', 'mol', 'type']): # if entries is either id, mol or type:
                    matrix[entries[j]][i] = int(matrix[entries[j]][i])   # ints
                    if (entries[j] == 'type'):
                        if (matrix[entries[j]][i] not in types):
                            types.append(matrix[entries[j]][i]) # add atom type to list of atom types
                else:
                    matrix[entries[j]][i] = float(matrix[entries[j]][i]) # floats
    '''
    for j in range(len(entries)):
        print matrix[entries[j]]
    '''
    return t,Natoms,system_size,matrix,readstructure, entries,types
    
    
def initialize(Natoms,matrix0,types,Nx,Nz,xmin,xmax,zmin,zmax):
    '''
    initialize(Natoms, matrix = {}, types = [], Nx=#bins in x dir. Nz=#bins in z dir, xmin E [0,1], xmax E [0,1])\n
    Where xmin < xmax\n
    Takes the initial statefile information and bins the atoms according to their position in space.
    It returns an array of bins that contain the atoms located at the bins position in the system.\n
    Returns two numpy arrays:\n
    binsx, binsz\n
    which are three dimentional matrices, where an atom j in bin i is accessed by:
    binsx[i][j] = [<id>, <xs>, <ys>, <zs>]
    '''

    binsx = []
    binsz = []
    natoms_x = np.zeros((Nx,1))
    natoms_z = np.zeros((Nx,1))
    boxsize_x = (xmax-xmin)/float(Nx)
    boxsize_z = (zmax-zmin)/float(Nz)
    for Bin in range(Nx): # fill with empty bins
        binsx.append([])
    for Bin in range(Nz):
        binsz.append([])
    
    for i in range(Natoms):
        if (matrix0['type'][i] in types):                     # only atoms in type are traced
            x = matrix0['xs'][i]
            if ((xmin < x) and (x < xmax)):                   # tracking only atoms inside nanopore
                index_z = int(round(matrix0['zs'][i]*(Nz-1))) # index of z bin in system
                index_x = int(round(x*(Nx-1)))                # index of x bin in system
                atom = [matrix0['id'][i],x,matrix0['ys'][i],matrix0['zs'][i]]
                natoms_x[index_x] += 1
                natoms_z[index_z] += 1
                binsx[index_x].append(atom)
                binsz[index_z].append(atom)

    print "##################################################################"
    print natoms_x
    print "##################################################################"
    print natoms_z
    print "##################################################################"
    return np.array(binsx),np.array(binsz)
    

def displacement(t, initial_binsx, initial_binsz, matrix, Natoms, types, Nx, Nz, xmin, xmax, zmin,zmax):
    '''
    displacementprofile(time, initial_binsx, initial_binsz, matrix, Natoms, types, Nx, Nz, xmin, xmax)\n 
    Takes the system state at the starting time t0, and a system state at time t.
    It creates bins in the x-y plane (along z-axis) and in the y-z plane (along x-axis) where it places the 
    atoms according to their location in space.
    Every atoms displacement from t1 to t2 is calculated only if it is still located in the bin it was in at
    time t1. If not. The atom gets a new initial position in the new bin, and it will therefore not contribute
    to any displacement.
    '''

    binsx = []
    binsz = []
    boxsize_x = (xmax-xmin)/float(Nx)
    boxsize_z = (zmax-zmin)/float(Nz)
    
    #still_in_binx = []
    #still_in_binz = []
    #atoms_to_removex = []
    #atoms_to_removez = []
    
    for Bin in range(Nx): # fill with empty bins
        binsx.append([])
        #still_in_binx.append([])
        #atoms_to_removex.append([])
    for Bin in range(Nz):
        binsz.append([])
        #still_in_binz.append([])
        #atoms_to_removez.append([])
    
    for i in range(Natoms):
        if (matrix['type'][i] in types):                      # only atoms in type are traced
            x = matrix['xs'][i]
            if ((xmin < x) and (x < xmax)):                   # tracking only atoms inside nanopore
                index_z = int(round(matrix['zs'][i]*(Nz-1)))  # index of z bin in system
                index_x = int(round(x*(Nx-1)))                # index of x bin in system
                atom = [matrix['id'][i],x,matrix['ys'][i],matrix['zs'][i]]

                binsx[index_x].append(atom)
                binsz[index_z].append(atom)

    # Now the atoms in the next timestep is added to bins. We must further go
    # through the bins, find the atoms that are still in the bin and calculate
    # the displacement.
    
    msd_x = np.zeros((Nx,1))
    msd_z = np.zeros((Nz,1))
    contributing_x = 0
    contributing_z = 0
    for j in range(Nx):
        for k in range(len(initial_binsx[j])): # all atoms in initial bin j
            a = initial_binsx[j][k][0]         # atom indexes in bin at last timestep
            count = 0
            for l in range(len(binsx[j])):     # all atoms in new bin j
                b = binsx[j][l][0]             # atom indexes in bin at this timestep.
                if (a == b):                   # If it is the same atom: 
                    dx = binsx[j][l][1] - initial_binsx[j][k][1]
                    dy = binsx[j][l][2] - initial_binsx[j][k][2]
                    dz = binsx[j][l][3] - initial_binsx[j][k][3]
                    msd_x[j] += dx*dx + dz*dz + dy*dy
                    #still_in_binx[j].append(b)
                    count += 1 # number of atoms that contribute to the msd
            contributing_x += count        
                    
            if (count != 0):
                msd_x[j] = msd_x[j]/count # mean square displacement in bin j.

    contributing_x = contributing_x/Nx # average number of contributing particles
    
    for j in range(Nz):
        for k in range(len(initial_binsz[j])): # all atoms in initial bin j
            a = initial_binsz[j][k][0]
            count = 0
            for l in range(len(binsz[j])):     # all atoms in new bin j
                if (a == binsz[j][l][0]):      # If it is the same atom: 
                    dx = binsz[j][l][1] - initial_binsz[j][k][1]
                    dy = binsz[j][l][2] - initial_binsz[j][k][2]
                    dz = binsz[j][l][3] - initial_binsz[j][k][3]
                    msd_z[j] += dx*dx + dz*dz + dy*dy
                    #still_in_binz[j].append(b)
                    count += 1 # number of atoms that contribute to the msd
            contributing_z += count
            if (count != 0):
                msd_z[j] = msd_z[j]/count # mean square displacement in bin j.
    

    print "avg # of contributing particles = %g , %g" % (contributing_x, contributing_z/Nz)
    print "----------------------------------------------"
    return initial_binsx, initial_binsz, msd_x, msd_z  # return the updated initial bins
        

def main():
    
    path = "/home/goran/lammps-28Jun14/examples/water_portlandite_system/npt_run_and_energyminimization/statefiles"
    
    arg = 'nve'  # use only filnames containing nve!
    filenames,time = gothroughfiles(path,arg) # get filenames in location "path"
    
    name = os.path.join(path,filenames[0])

    t,Natoms,system_size,matrix,readstructure,entries,types = readfile(name)

    lorock = 7.2   # height of lower portlandite part
    hirock = 36.7  # height of upper portlandite part
    xmin = 0.25
    xmax = 0.75
    zmin = 0.0
    zmax = 1.0    
    Nx = Nz = 15
    Types = [4]    # oxygen
    binsx, binsz = initialize(Natoms,matrix,Types,Nx,Nz,xmin,xmax,zmin,zmax)
    
    msdx = []
    msdz = []
    time = []
    conversionfactor = 10**(-4.0)
    kk = 0
    nfiles = len(filenames) -1
    for name in filenames[2:]:
        kk += 1
        print "Processing %s . file nr %g / %g" % (name,kk,nfiles)
        filename = os.path.join(path,name)
        t,Natoms,system_size,matrix,readstructure,entries,types = readfile(filename)
        binsx,binsz,msd_x,msd_z = displacement(t,binsx,binsz,matrix,Natoms,Types,Nx,Nz,xmin,xmax,zmin,zmax)
        msdx.append(msd_x), msdz.append(msd_z), time.append(t)
    
    xlo = system_size[0]
    xhi = system_size[1]
    zlo = system_size[4]
    zhi = system_size[5]
    lx = (xhi-xlo)          # length of system in x dir.
    xstart = lx*xmin        # starting pos of portlandite in x dir
    xstop = lx*xmax         # ending pos in x dir
    dx = (xstop-xstart)/Nx 
    dz = (zhi-zlo)/Nz
    binposx = np.zeros((Nx,1))
    binposz = np.zeros((Nz,1))
    msdx_time = np.zeros((Nx,len(time)))
    msdz_time = np.zeros((Nz,len(time)))
    Dx_time = np.zeros((Nx,len(time)))
    Dz_time = np.zeros((Nz,len(time)))
    lower_rock = np.linspace(lorock,lorock,Nx)
    upper_rock = np.linspace(hirock,hirock,Nx)
    bins = np.linspace(1,Nx,Nx)  # for now
    
    # --------- Averaging the displacement in the bins ---------#    
    mean_sdx = np.zeros((Nx,1))
    mean_sdz = np.zeros((Nx,1))
    length_msdx = len(msdx)
    length_msdz = len(msdz)    
    for i in range(length_msdx):
        time[i] = time[i]*conversionfactor
        for j in range(Nx):
            mean_sdx[j] += msdx[i][j]
            msdx_time[j][i] = msdx[i][j]
            if (time[i] > 0):
                Dx_time[j][i] = mean_sdx[j]/(6*time[i])
                
    for i in range(length_msdz):
        for j in range(Nz):
            mean_sdz[j] += msdz[i][j]
            msdz_time[j][i] = msdz[i][j]
            if (time[i] > 0):
                Dz_time[j][i] = mean_sdz[j]/(6*time[i])
    
    for j in range(Nx):
        mean_sdx[j] = mean_sdx[j]/length_msdx
        binposx[j] = dx*(j+0.5) + xstart
    for j in range(Nz):
        mean_sdz[j] = mean_sdz[j]/length_msdz
        binposz[j] = dz*(j+0.5)
    
    time = np.array(time)

    #---------------------------------- plotting ---------------------------#
    fig = plt.figure()
    plt.hold(True)
    #plt.plot(binposx,mean_sdx,'b-*')
    for i in range(length_msdx):
#        plt.plot(binposx,msdx[i])
        plt.plot(bins,msdx[i])
    plt.hold(False)
    plt.title('Mean square displacement of particles as function of\ndistance into the gap/nano-tube.')
    plt.xlabel('x [Angstrom]'),plt.ylabel('msd [A^2]')

    fig2 = plt.figure()
    plt.hold(True)
    #plt.plot(binposz,mean_sdz,'b-*')
    #lspace = np.linspace(min(mean_sdz),max(mean_sdz),Nx)
    #plt.plot(lower_rock,lspace,'k--')
    #plt.plot(upper_rock,lspace,'k--')
    for i in range(length_msdz):
#        plt.plot(binposz,msdz[i])
        plt.plot(bins,msdz[i])
    plt.hold(False)
    plt.title('Mean square displacement of particles as function of\ndistance from the walls')
    plt.xlabel('z [Angstrom]'),plt.ylabel('msd [A^2]')

    fig3 = plt.figure()
    plt.hold(True)
    legends = []
    for j in range(int(round(Nx/2.0)+1)):
        legend = "bin %g" % (j+1)
        plt.plot(time[:],msdx_time[j][:],'-')
        legends.append(legend)
    plt.hold(False)
    plt.title('Displacement as funciton of time for x')
    plt.xlabel('time [ps]'),plt.ylabel('msd [A^2]')
    plt.legend(legends,loc='upper left')    
    
    fig4 = plt.figure()
    plt.hold(True)
    legends = []
    for j in range(int(round(Nz/2.0)+1)):
        legend = "bin %g" % (j+1)   
        plt.plot(time[:],msdz_time[j][:],'-')
        legends.append(legend)
    plt.hold(False)
    plt.title('Displacement as funciton of time for z')
    plt.xlabel('time [ps]'),plt.ylabel('msd [A^2]')
    plt.legend(legends,loc='upper left')
    
    fig5 = plt.figure()
    plt.hold(True)    
    for j in range(int(round(Nx/2.0)+1)):
        plt.plot(time,Dx_time[j][:],'-')
    plt.hold(False)
    plt.title('Diffution for bins along the x-axis.\nThat is, as a funciton debth into the nanopore')
    plt.xlabel('time [ps]'),plt.ylabel('D [A^2/s]')
    plt.legend(legends,loc='upper left')
    
    fig5 = plt.figure()
    plt.hold(True)    
    for j in range(int(round(Nz/2.0)+1)):
        plt.plot(time,Dz_time[j][:],'-')
    plt.hold(False)
    plt.title('Diffusion for bins in z-direction.\nThat is, parallell to the rock surface')
    plt.xlabel('time [ps]'),plt.ylabel('D [A^2/s]')
    plt.legend(legends,loc='upper left')
    
    plt.show(True)

        

if (__name__ == "__main__"):
    import numpy as np
    import os
    import matplotlib.pyplot as plt
    main()