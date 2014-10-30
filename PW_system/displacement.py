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
    
    
def initialize(Natoms,matrix0,types,Nz,xmin,xmax,z0,z1):
    '''
    initialize(Natoms, matrix = {}, types = [], Nz=#bins in z dir, xmin E [0,1], xmax E [0,1]), z0,z1 is the z position of the surfaces\n
    Where xmin < xmax\n
    Takes the initial statefile information and bins the atoms according to their position in space.
    It returns an array of bins that contain the atoms located at the bins position in the system.\n
    Returns two numpy arrays:\n
    binsx, binsz\n
    which are three dimentional matrices, where an atom j in bin i is accessed by:
    binsz[i][j] = [<id>, <xs>, <ys>, <zs>]
    '''
    binsz = []
    natoms_z = np.zeros((Nz,1))
    d_port = (z1-z0)                       # crack opening
    half_d_port = d_port/2.0               # half crack opening
    boxsize_z = (half_d_port)/float(Nz)    # boxsize for boxes from surface to center of crack.
    center_of_port = half_d_port + z0      # center of system

    for Bin in range(Nz):
        binsz.append([])
    
    for i in range(Natoms):
        if (matrix0['type'][i] in types):                     # only atoms in type are traced
            x = matrix0['xs'][i]
            if ((xmin < x) and (x < xmax)):                   # tracking only atoms inside nanopore
                z = matrix0['zs'][i]
                if (z < center_of_port):
                    dz = (z-z0)  # distance from lower
                else:
                    dz = (z1-z)  # distance from upper
                if (dz < 0):
                    index_z = 0
                    print "Shit! index_z=%g, dz=%g, z=%g : center_of_port=%g, boxsize=%g" % (index_z,dz,z,center_of_port,boxsize_z)
                else:
                    index_z = int(np.floor(dz/boxsize_z))  # bin index for the i'th atom
                if (index_z >= Nz):
                    print "Shit! index_z=%g, dz=%g, z=%g : center_of_port=%g, boxsize=%g" % (index_z,dz,z,center_of_port,boxsize_z)
                    index_z = Nz-1

                atom = [matrix0['id'][i],x,matrix0['ys'][i],z]
                natoms_z[index_z] += 1 # count the number of atoms that are appended in each
                binsz[index_z].append(atom)

    print "##################################################################"
    print natoms_z
    print "##################################################################"
    return np.array(binsz)
    

def displacement(t, initial_binsz, matrix, Natoms, types, Nz, xmin, xmax, z0,z1):
    '''
    displacementprofile(time, initial_binsx, initial_binsz, matrix, Natoms, types, Nx, Nz, xmin, xmax)\n 
    Takes the system state at the starting time t0, and a system state at time t.
    It creates bins in the x-y plane (along z-axis) and in the y-z plane (along x-axis) where it places the 
    atoms according to their location in space.
    Every atoms displacement from t1 to t2 is calculated only if it is still located in the bin it was in at
    time t1. If not. The atom gets a new initial position in the new bin, and it will therefore not contribute
    to any displacement.
    '''

    binsz = []
    d_port = (z1-z0)                       # crack opening
    half_d_port = d_port/2.0               # half crack opening
    boxsize_z = (half_d_port)/float(Nz)    # boxsize for boxes from surface to center of crack.
    center_of_port = half_d_port + z0

    for Bin in range(Nz):
        binsz.append([])
    
    for i in range(Natoms):
        if (matrix['type'][i] in types):                      # only atoms in type are traced
            x = matrix['xs'][i]
            if ((xmin < x) and (x < xmax)):                   # tracking only atoms inside nanopore
                z = matrix['zs'][i]
                if (z < center_of_port):
                    dz = (z-z0)      # distance from closest crack surface
                else:
                    dz = (z1-z)      # distance from closest crack surface
                if (dz < 0):
                    index_z = 0
                else:
                    index_z = int(np.floor(dz/boxsize_z))  # bin index for atom
                if (index_z >= Nz):
                    #print "Shit! index_z=%g, dz=%g, z=%g : center_of_port=%g, boxsize=%g" % (index_z,dz,z,center_of_port,boxsize_z)
                    index_z = Nz-1
                    
                atom = [matrix['id'][i],x,matrix['ys'][i],z]

                binsz[index_z].append(atom)

    # Now the atoms in the next timestep is added to bins. We must further go
    # through the bins, find the atoms that are still in the bin and calculate
    # the displacement.
    
    msd_z = np.zeros((Nz,1))
    contributing_z = 0
    
    for j in range(Nz):
        for k in range(len(initial_binsz[j])): # all atoms in initial bin j
            a = initial_binsz[j][k][0]         # atom_ID
            count = 0
            for l in range(len(binsz[j])):     # all atoms in new bin j
                if (a == binsz[j][l][0]):      # If the atom is still in bin j: 
                    dx = binsz[j][l][1] - initial_binsz[j][k][1]
                    dy = binsz[j][l][2] - initial_binsz[j][k][2]
                    dz = binsz[j][l][3] - initial_binsz[j][k][3]
                    msd_z[j] += dx*dx + dz*dz + dy*dy
                    count += 1 # number of atoms that contribute to the msd
            contributing_z += count
            if (count != 0):
                msd_z[j] = msd_z[j]/count # mean square displacement in bin j.
    

    print "avg # of contributing particles = %g" % (contributing_z/Nz)
    print "----------------------------------------------"
    return initial_binsz, msd_z  # return the updated initial bins
        

def main():
    
    path = "/home/goran/lammps-28Jun14/examples/water_portlandite_system/npt_run_and_energyminimization/statefiles"
    
    arg = 'nve'  # use only filnames containing nve!
    filenames,time = gothroughfiles(path,arg) # get filenames in location "path"
    
    name = os.path.join(path,filenames[0])

    t,Natoms,system_size,matrix,readstructure,entries,types = readfile(name)

    portlandite_thickness = 11.3  # Angstrom
    lorock = 7.0         # height of lower portlandite part
    hirock = 37.8        # height of upper portlandite part
    zlo = system_size[4]
    zhi = system_size[5]
    zsize = (zhi-zlo)
    #print zlo, zhi,zsize
    z0 = portlandite_thickness/zsize
    z1 = (zsize - portlandite_thickness)/zsize
    xmin = 0.28
    xmax = 0.72    
    Nz = 5
    d_port = (z1-z0)                       # crack opening
    half_d_port = d_port/2.0               # half crack opening
    boxsize_z = (half_d_port)/float(Nz)    # boxsize for boxes from surface to center of crack.
    center_of_port = half_d_port + z0
    Types = [4]          # oxygen
    binsz = initialize(Natoms,matrix,Types,Nz,xmin,xmax,z0,z1)
    
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
        binsz,msd_z = displacement(t,binsz,matrix,Natoms,Types,Nz,xmin,xmax,z0,z1)
        msdz.append(msd_z), time.append(t)
    
    
    
    # ----------- Plotting n'shit ----------------------------
    
    save = True    # save figures
    Formats = ['png','jpg','jpeg']
    fig1_name = 'msd_evolution_binnr_' + arg +'.png'
    fig2_name = 'msd_evolution_dist_' + arg +'.png'
    fig3_name = 'msd_time_' + arg +'.png'
    fig4_name = 'Diffusion_bins' + arg +'.png'
    
    ntimesteps = len(msdz)  
    
    fig = plt.figure()
    fig1_name = 'msd_evolution_binnr_' + arg
    bins = np.linspace(1,Nz,Nz) # bins
    dist = np.linspace(zsize*(boxsize_z/2.0),zsize*(Nz*boxsize_z/2.0),Nz)  # midpoints of boxes/bins
    plt.hold(True)
    for i in range(ntimesteps):
        plt.plot(bins,msdz[i],'-*')
    plt.hold(False)
    plt.xlabel('bin [number]'),plt.ylabel('msd [A^2/ps]')
    plt.title('Mean square displacement evolution\nAs a function of  distance from the rock surface')
    if (save == True):
        for Format in Formats:
            plt.savefig(fig1_name,format=Format)
    
    fig2 = plt.figure()
    plt.hold(True)
    for i in range(ntimesteps):
        plt.plot(dist,msdz[i],'-*')
    plt.hold(False)
    plt.xlabel('z [Angstrom]'),plt.ylabel('msd [A^2/ps]')
    plt.title('Mean square displacement evolution\nAs a function of  distance from the rock surface')
    if (save == True):
        for Format in Formats:
            plt.savefig(fig2_name,format=Format)
        
    
    bin_msdz = np.zeros((Nz,ntimesteps))
    Dz = np.zeros((Nz,ntimesteps))
    Dz_si = np.zeros((Nz,ntimesteps))
    Aaps_to_si_units = 10**(-8.0)
    tid = np.zeros((ntimesteps,1))
    for i in range(ntimesteps):
        for j in range(Nz):
            value = msdz[i][j]
            bin_msdz[j][i] = value
            tid[i] = time[i]*conversionfactor
            Dz[j][i] = value/(6.0*tid[i])
            Dz_si[j][i] = Dz[j][i]*Aaps_to_si_units 
            

    legends = []
    Title1 = 'Mean square displacement for different distances\nfrom the surface wall of the portlandite'
    Title2 = 'Diffusion coefficient for different distances\nfrom the surface wall of the portlandite. [Aa^2/ps] = %g [m^2/s]' % Aaps_to_si_units
    xlabel = 'time [ps]'
    ylabel1 = 'msd []'
    ylabel2 = 'D [A^2/ps]'
    location = 'upper left'
    
    linestyle = ['-','--','-.',':','-','--','-.',':','-','--','-.',':','-','--','-.',':','-','--','-.',':','-','--','-.',':','-','--','-.',':','-','--','-.',':','-','--','-.',':']
    markers = ['*','o','d','x','v','^','<','>','1','2','3','4','8','s','p','+','*','o','d','x','v','^','<','>','1','2','3','4','8','s','p','+']
    
    fig3 = plt.figure()
    plt.hold(True)
    for j in range(Nz):
        string = linestyle[j] + markers[j]
        legend = 'bin %g' % j
        legends.append(legend)
        plt.plot(tid,bin_msdz[j],string)
    plt.hold(False)
    plt.title(Title1)
    plt.xlabel(xlabel), plt.ylabel(ylabel1), plt.legend(legends,loc=location)
    if (save == True):
        for Format in Formats:
            plt.savefig(fig3_name,format=Format)
        
    fig4 = plt.figure()
    plt.hold(True)
    for j in range(Nz):
        string = linestyle[j] + markers[j]        
        plt.plot(tid,Dz[j],string)
    plt.hold(False)
    plt.title(Title2)
    plt.xlabel(xlabel), plt.ylabel(ylabel2), plt.legend(legends,loc=location)
    if (save == True):
        for Format in Formats:
            plt.savefig(fig4_name,format=Format)
        
    
    #-- Calculate the diffusion coefficient for the bins, and plot them as function of dist from surface--
    straight_line = 1 # degree of polynomial
    
    p,v = np.polyfit(tid,bin_msdz[0],straight_line)
    # p = ndarray, polynomial coefficients 
    plt.show(True)

    '''
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
    '''
        

if (__name__ == "__main__"):
    import numpy as np
    import os
    import matplotlib.pyplot as plt
    main()