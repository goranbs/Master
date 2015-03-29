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

def get_savepath(currentpath,pathgiven,truevalues):
    #print "path is from the submit-folder. That is: %s \n" % currentpath
    saveindatabase = raw_input("Do you want to save the plots in the directory you gave on the cmd (yes/no)? ")
    if (saveindatabase in truevalues):
        fullpathtosavedir = os.path.abspath(os.path.join(currentpath,pathgiven))
        fullpathtosavedir = os.path.dirname(fullpathtosavedir)
        print "save in %s" % fullpathtosavedir
    else:    
        print "Where do you want to save the plots?"
        savepath = raw_input("Give path: ")
        if (savepath[0] == "~"):
            home = os.path.expanduser("~")
            savepath = home + savepath[1:]
            fullpathtosavedir = os.path.abspath(savepath)
        else: 
            fullpathtosavedir = os.path.abspath(os.path.join(currentpath,savepath))
            print "path chosen: %s " % fullpathtosavedir
    
        if not os.path.isdir(fullpathtosavedir): # if the path is unique, do you wish to create it?
            createfolder = raw_input("not a valid path. Do you wish to create a folder here? (yes/no) ")
            if (createfolder in ["Y","y","Yes","yes","true","True"]): 
                print "creating folder...."
                os.makedirs(fullpathtosavedir) # use default mode
            else:
                print "Ok then..."
                fullpathtosavedir = get_savepath(currentpath)
        else: 
            print "The path '%s'  is an existing path," % fullpathtosavedir
            savehere = raw_input("Do you wish to save plots here? (yes/no) ")
            if not (savehere in ["Y","y","Yes","yes","true","True","jepp","Jepp"]):
                print "Ok then..."
                fullpathtosavedir = get_savepath(currentpath)
    
    #print "save in: %s " % fullpathtosavedir
    return fullpathtosavedir

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
        writefile = open(filename,'w')
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
            writefile.write(line)
    
        
        writefile.close()
        print "file %s is done." % filename
        
        
def gothroughfiles(path,arg=None):
    '''
    Go through stystem state files in "path", if they are .txt files, the first 
    character in the filename which is a number will be acknowledged as the 
    timestep of the system state file, and the number will be read as:
    somefilename<somenumber>.txt , then <somenumber> is the time that will be
    read.\n
    Returns:\n
            sortedfilenames, time
    '''
    import os
    
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
                    #print start
                    try:
                        anumber = int(name[start:-4])                
                        numbers.append(anumber)
                    except:
                        pass
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
    boxbounds = infile.readline()
    nvals = boxbounds.split()
    line1 = infile.readline() # system size x-axis
    line2 = infile.readline()
    line3 = infile.readline()
    if (len(nvals) == 9):
        systembounds = 'triclinic'
        #print systembounds
        xlo,xhi,xy = line1.split()
        ylo,yhi,xz = line2.split()
        zlo,zhi,yz = line3.split()
        
    elif (len(nvals) == 6):
        systembounds = 'orthogonal'
        #print systembounds
        xlo, xhi = line1.split() # system size y-axis
        ylo, yhi = line2.split() # system size y-axis
        zlo, zhi = line3.split() # system size z-axis
    else:
        print "Error! triclinic or orthogonal?"
        print boxbounds, len(nvals)
    


    
    system_size = [float(xlo),float(xhi),float(ylo),float(yhi),float(zlo),float(zhi)]
    #structure of file: 
    
    line1 = infile.readline()
    columns = line1.split()
    # the two first words are ITEM: ATOMS
    c_id=c_mol=c_type=c_xs=c_ys=c_zs=c_vx=c_vy=c_vz=c_fx=c_fy=c_fz=c_x=c_y=c_z= -1
    ID=MOL=TYPE=XS=YS=ZS=VX=VY=VZ=FX=FY=FZ=X=Y=Z= None
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
        if column == "x":
            c_x = it
            X = []
            entries.append('x')
            col.append(c_x)
        if column == "y":
            c_y = it
            Y = []
            entries.append('y')
            col.append(c_y)
        if column == "z":
            c_z = it
            Z = []
            entries.append('z')
            col.append(c_z)
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
  
    matrix = {'id':ID,'mol':MOL,'type':TYPE,'xs':XS,'ys':YS,'zs':ZS,'vx':VX,'vy':VY,'vz':VZ,'fx':FX,'fy':FY,'fz':FZ,'x':X,'y':Y,'z':Z}
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
    initialize(Natoms, matrix = {}, types = [], Nz=#bins in z dir, xmin E [0,1], xmax E [0,1]), y0,y1 is the y position of the surfaces\n
    Where xmin < xmax\n
    Takes the initial statefile information and bins the atoms according to their position in space.
    It returns an array of bins that contain the atoms located at the bins position in the system.\n
    Returns two numpy arrays:\n
    binsx, binsy\n
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
                z = matrix0['ys'][i]
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

                atom = [matrix0['id'][i],matrix0['xs'][i],matrix0['ys'][i],matrix0['zs'][i]]
                natoms_z[index_z] += 1 # count the number of atoms that are appended in each
                binsz[index_z].append(atom)

    print "##################### Initial distribution #######################"
    for i in range(Nz):
        print "Binnr: %g #atoms= %g" % ((i+1),natoms_z[i])
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
    
    msd = np.zeros((Nz,1)); msd_z = np.zeros((Nz,1)); msd_r = np.zeros((Nz,1))
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
                    msd[j] += dx*dx + dz*dz + dy*dy               # total disp
                    msd_z[j] += dz*dz                             # disp in z direction
                    msd_r[j] += dx*dx + dy*dy                     # disp in the xy-plane
                    count += 1 # number of atoms that contribute to the msd
            contributing_z += count
            if (count != 0):
                msd[j] = msd[j]/count      # msd in bin j.
                msd_z[j] = msd_z[j]/count
                msd_r[j] = msd_r[j]/count
    

    print "avg # of contributing particles = %g" % (contributing_z/Nz)
    print "----------------------------------------------"
    return initial_binsz, msd, msd_z, msd_r
        

def diffusion(ofile,path,date,Formats=None,arg=None,showplots=True,saveplots=True):
    '''
    Goes through LAMMPS output files stored in "path", and selects the contains "arg". If "arg" is\n
    not given, then all files in the directory given by "path" that are .txt files are read.\n
    The value of "date" should be given to specify namecontent of the output plots
    '''    
    rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
    ## for Palatino and other serif fonts use:
    #rc('font',**{'family':'serif','serif':['Palatino']})
    rc('text', usetex=True) 
    
    save = saveplots
    if (Formats is None):
        Formats = ['png']
    if (arg is not None):
        name = arg + '_' + date
    else:
        name = date
    
    fig1_name = 'msd_evolution_binnr_' + name                # msd bin
    fig2_name = 'msd_evolution_dist_' + name                 # msd dist
    fig3_name = 'msdz_binnr_' + name                         # msdz bin
    fig4_name = 'msdr_binnr' + name                          # msdr bin
    fig5_name = 'msd_time_' + name                           # msd as function of time
    fig6_name = 'Diffusion_time_' + name                     # Diffusion development
    fig7_name = 'Diffusion_dist_' + name                     # Diffusion constant
    
    xlabel = 'time [ps]'
    ylabel1 = 'msd []'
    ylabel2 = r'D $[A^2/ps]$'
    location = 'upper left'
    
    linestyle = ['-','--','-.',':','-','--','-.',':','-','--','-.',':','-','--','-.',':','-','--','-.',':','-','--','-.',':','-','--','-.',':','-','--','-.',':','-','--','-.',':']
    markers = ['*','o','d','x','v','^','<','>','1','2','3','4','8','s','p','+','*','o','d','x','v','^','<','>','1','2','3','4','8','s','p','+']
    
    
    filenames,time = gothroughfiles(path,arg) # get filenames in location "path"
    
    name = os.path.join(path,filenames[0])

    t,Natoms,system_size,matrix,readstructure,entries,types = readfile(name)

    portlandite_thickness = 21.0  # Angstrom
    ylo = system_size[2]          # system origin y
    yhi = system_size[3]          # system end y
    ysize = (yhi-ylo)             # length of system in y dir.
    y0 = portlandite_thickness/ysize             # relative thickness of the portlandite
    y1 = (ysize - portlandite_thickness)/ysize   # approximate relative starting point of upper portlandite layer
    xmin = 0.0                                   # we only want to look at the diffusion inside the pore
    xmax = 1.0                                   # we only want to look at the diffusion inside the pore
    Nz = 5                                       # number of boxes!
    d_port = (y1-y0)                             # width of pore opening
    half_d_port = d_port/2.0                     # half crack opening
    boxsize_y = (half_d_port)/float(Nz)          # boxsize for boxes from surface to center of crack.
    center_of_port = half_d_port + y0            # approximate center
    Types = [4]                                  # oxygen. only trace oxygen!
    binsy = initialize(Natoms,matrix,Types,Nz,xmin,xmax,y0,y1)
    
    msd = []; msdz = []; msdr = []; time = []
    conversionfactor = 10**(-4.0)
    kk = 0
    nfiles = len(filenames) -1
    for name in filenames[2:]:
        kk += 1
        print "Processing %s . file nr %g / %g" % (name,kk,nfiles)
        filename = os.path.join(path,name)
        t,Natoms,system_size,matrix,readstructure,entries,types = readfile(filename)
        binsy,msd_tot,msd_z,msd_r = displacement(t,binsy,matrix,Natoms,Types,Nz,xmin,xmax,y0,y1)
        msd.append(msd_tot), msdz.append(msd_z), msdr.append(msd_r), time.append(t)
    
    zsize = ysize
    boxsize_z = boxsize_y
    
    # ----------- Plotting n'shit ----------------------------    
    
    ntimesteps = len(msd)  
    legends = []

    #------------------------------------- FIGURE 1 -----------------------------------------------
    # MSD evolution in time as a funtion of distance from the surface.

    plt.figure()
    Title = 'Mean square displacement evolution in time\nas a function of distance from the rock surface'
    Xlabel = 'bin [number]'
    Ylabel = 'msd $[A^2/ps]$'
    linestyles = []
    bins = np.linspace(1,Nz,Nz) # bins
    plt.hold(True)
    for i in range(ntimesteps):
        linestyle = '-*'
        linestyles.append(linestyle)
        plt.plot(bins,msd[i],linestyle)
    plt.hold(False)
    plt.xlabel(Xlabel),plt.ylabel(r'msd $[A^2/ps]$')
    plt.title(Title)
    if (save == True):
        for Format in Formats:
            fig_name = fig1_name + '.' + Format
            plt.savefig(fig_name,format=Format)
    
    write_to_file(ofile,Title,Xlabel,Ylabel,fig1_name,legends,linestyles,[bins],msd)
    
    #------------------------------------- FIGURE 2 -----------------------------------------------
    # MSD evolution in time as a funtion of distance from the surface.
    strt = zsize*(boxsize_z/2.0)
    end = Nz*zsize*boxsize_z - strt
    dist = np.linspace(strt,end,Nz)  # midpoints of boxes/bins
    
    plt.figure()
    Title = 'Mean square displacement evolution in time \nas a function of distance from the rock surface'
    Xlabel = 'z [Angstrom]'
    Ylabel = 'msd $[A^2/ps]$'
    linestyles = []
    plt.hold(True)
    for i in range(ntimesteps):
        linestyle = '-*'
        linestyles.append(linestyle)
        plt.plot(dist,msd[i],linestyle)
    plt.hold(False)
    plt.xlabel(Xlabel),plt.ylabel(r'msd $[A^2/ps]$')
    plt.title(Title)
    if (save == True):
        for Format in Formats:
            fig2_name = fig2_name + '.' + Format
            plt.savefig(fig2_name,format=Format)

    write_to_file(ofile,Title,Xlabel,Ylabel,fig2_name,legends,linestyles,[dist],msd)

    #------------------------------------- FIGURE 3 -----------------------------------------------
    # MSD using the z component. As function of bins to the surface for the different times.
    
    plt.figure()
    Title = 'Mean square displacement normal to the surface\nas a function of distance from the surface'
    Xlabel = 'bin [number]'
    Ylabel = 'msd $[Aa^2/ps]$'
    plt.hold(True)
    for i in range(ntimesteps):
        plt.plot(bins,msdz[i],'-*')
    plt.hold(False)
    plt.xlabel(Xlabel),plt.ylabel(r'msd $[Aa^2/ps]$')
    plt.title(Title)
    if (save == True):
        for Format in Formats:
            fig3_name = fig3_name + '.' + Format
            plt.savefig(fig3_name,format=Format)

    write_to_file(ofile,Title,Xlabel,Ylabel,fig3_name,legends,linestyles,[bins],msdz)

    #------------------------------------- FIGURE 4 -----------------------------------------------
    # MSD using the x,y component. As function of bins to the surface for the different times.            
    plt.figure()
    Title = 'Mean square displacement parallell to the surface\nas a function of distance from the surface'
    Xlabel = 'bin [number]'
    Ylabel = 'msd $[Aa^2/ps]$'
    plt.hold(True)
    for i in range(ntimesteps):
        plt.plot(bins,msdr[i],'-*')
    plt.hold(False)
    plt.xlabel(Xlabel),plt.ylabel(r'msd $[Aa^2/ps]$')
    plt.title(Title)
    if (save == True):
        for Format in Formats:
            fig_name = fig4_name + '.' + Format
            plt.savefig(fig_name,format=Format)
        
    write_to_file(ofile,Title,Xlabel,Ylabel,fig4_name,legends,linestyles,[bins],msdr)
    
    bin_msd = np.zeros((Nz,ntimesteps));bin_msdz = np.zeros((Nz,ntimesteps));bin_msdr = np.zeros((Nz,ntimesteps))

    Dz = np.zeros((Nz,ntimesteps))      # aangstrom^2/picoseconds
    Dz_si = np.zeros((Nz,ntimesteps))    # SI units
    Aaps_to_si_units = 10**(-8.0)
    tid = []
    for i in range(ntimesteps):
        tid.append(time[i]*conversionfactor)        
        for j in range(Nz):
            value = msd[i][j]
            value_z = msdz[i][j]
            value_r = msdr[i][j]
            bin_msd[j][i] = value
            bin_msdz[j][i] = value_z
            bin_msdr[j][i] = value_r
            Dz[j][i] = value/(6.0*tid[i])
            Dz_si[j][i] = Dz[j][i]*Aaps_to_si_units 
    
    
    #------------------------------------- FIGURE 5 -----------------------------------------------
    # MSD as function of time. Also, we are calculating the diffusion coefficients
    Title1 = 'Mean square displacement for different distances from the\n surface wall of the portlandite'
    Title2 = r'Diffusion coefficient for different distances from the surface\n $[Aa^2/ps] = %g [m^2/s]$' % Aaps_to_si_units
   
    degree = 1 # degree of polynomial
    D = []; D_line = []; start = int(len(tid)/4.0)
    Dz_ = []; Dr_ = []
    plt.figure()
    plt.hold(True)
    legends = []

    for j in range(Nz):
        line = "--"
        p = np.polyfit(tid[start:],bin_msd[j][start:],degree)
        f = np.polyval(p,tid[start:])
        D.append(p[0])
        D_line.append(f)

        pz = np.polyfit(tid[start:],bin_msdz[j][start:],degree)
        #fz = np.polyval(p2,tid[start:])

        pr = np.polyfit(tid[start:],bin_msdr[j][start:],degree)
        #fr = np.polyval(p2,tid[start:])
          
        Dz_.append(pz[0])
        Dr_.append(pr[0])
        
        string = linestyle[j] + markers[j]
        legend = r'$bin %g$' % (j+1)
        legend2 = r"$D_{%g}=%g$" % ((j+1),D[j])
        legends.append(legend)
        legends.append(legend2)
        plt.plot(tid,bin_msd[j],string)
        plt.plot(tid[start:],f,line)
    plt.hold(False)
    plt.title(Title1)
    plt.xlabel(xlabel), plt.ylabel(ylabel1), plt.legend(legends,loc=location)
    if (save == True):
        for Format in Formats:
            fig_name = fig5_name + '.' + Format
            plt.savefig(fig_name,format=Format)

    #------------------------------------- FIGURE 6 -----------------------------------------------
    # Diffution coefficient as a function of time for the different bins.
 
    plt.figure()
    plt.hold(True)
    legends = []
    linestyles = []
    for j in range(Nz):
        linestyle = "-" + markers[j]
        linestyles.append(linestyle)
        plt.plot(tid,Dz[j],linestyle)
        legend1 = r'$D_{j=%g}=%g$' % (j+1,D[j])
        legends.append(legend1)
    plt.hold(False)
    plt.title(Title2)
    plt.xlabel(xlabel), plt.ylabel(ylabel2), plt.legend(legends,loc=location)
    if (save == True):
        for Format in Formats:
            fig_name = fig6_name + '.' + Format            
            plt.savefig(fig_name,format=Format)

    write_to_file(ofile,Title2,xlabel,ylabel2,fig6_name,legends,linestyles,[tid],[Dz])

    #------------------------------------- FIGURE 7 -----------------------------------------------
    # Approximation of the diffusion coefficient as a function of distance to the surface.

    plt.figure()
    Title = 'Diffusion constant as function of distance from the surface'
    Xlabel = 'd $[Aa]$'
    legends = ['D(d)','Dz(d)','Dr(d)']
    linestyles = ['--*','--o','--x']
    plt.plot(dist,D,linestyles[0])
    plt.hold(True)
    plt.plot(dist,Dz_,linestyles[1])
    plt.plot(dist,Dr_,linestyles[2])
    plt.hold(False)
    plt.title(Title)
    plt.xlabel(r'd $[Aa]$'), plt.ylabel(ylabel2), plt.legend(legends,loc='upper left')
    if (save == True):
        for Format in Formats:
            fig_name = fig7_name + '.' + Format
            plt.savefig(fig_name, format=Format)
            

    write_to_file(ofile,Title,Xlabel,Ylabel,fig7_name,legends,linestyles,[dist],[D,Dz_,Dr_])
    
    plt.figure()
    figurename = fig7_name + '_funcOfBins.'
    plt.plot(bins,D,'--*')
    plt.hold(True)
    plt.plot(bins,Dz_,'--d')
    plt.plot(bins,Dr_,'--x')
    plt.hold(False)
    plt.title('Diffusion constant as function of bins from the surface')
    plt.xlabel('bin'), plt.ylabel(ylabel2), plt.legend(['D(d)','Dz(d)','Dr(d)'],loc='upper left')
    if (save == True):
        for Format in Formats:
            fig_name =  + Format
            plt.savefig(fig_name, format=Format)
     
    write_to_file(ofile,Title,Xlabel,Ylabel,figurename,legends,linestyles,[dist],[D,Dz_,Dr_])

    plt.show(showplots)
    plt.close('all')


def printhelp(progname):
    '''
    Print run information about this to the screen
    '''
    print "#########################################"
    print "# INFO about %s \n" % progname
    print "# >>python %s [path-to-dir] [showplots=yes/no] [saveplots=yes/no] \n#"
    print "########################################"

def main(argv):


    #path = "/home/goran/lammps-28Jun14/examples/Abel_runs/evaporate/co2_system_npt"
    
    pathtodir = argv[1]
    
    truevalues = ["Y","y","Yes","yes","true","True"]
    falsevalues = ["N","n","No","no","false","False"]
    
    if (argv[1] in ["-h", "h", "Help", "help", "HELP", "--help"]):
        printhelp(argv[0])
        sys.exit(1)
    
    if ((len(argv) < 2) or (len(argv) > 4)):
        sys.stderr.write("Useage:python %s <path-to-plotfile> <showplot=yes/no/None> <saveplot=yes/no/None>" % (argv[0]))
        return 1
    if not os.path.exists(pathtodir):
        sys.stderr.write("Error: Database %s was not found!" % (argv[1]))
        return 1
    if ((len(argv) == 3) or (len(argv) == 4)):
        showplot = argv[2]
        if (showplot in truevalues):
            showplot = True
            print "Show plots = True"
        if (showplot in falsevalues):
            showplot = False
            print "Show plots = False"
        if (len(argv) == 3):
            print "Commandline argument for saveplot not given. saveplot = False"
    if (len(argv) == 4):
        saveplot = argv[3]
        if (saveplot in truevalues):
            saveplot = True
            print "Save plots = True"
            pathtosavedir = get_savepath(os.getcwd(),pathtodir,truevalues)
        if (saveplot in ["N","n","No","no","false","False"]):
            saveplot = False
            print "Save plots = False"
    
    date = str(datetime.date.today())
    showplots = True
    saveplots = False
    Formats = ['png'] #Formats = ['png','jpg','jpeg']
    outputfilename = os.path.basename(pathtodir) + '-' + date + '.dat'
    ofile = open(outputfilename,'w')
    
    print "Do you wish to give keyword for files in database? None for no kwyword:"
    arg = raw_input('Give keyword: ')
    
    
    diffusion(ofile,pathtodir,date,Formats,arg,showplots,saveplots)

    ofile.close()
         

if (__name__ == "__main__"):
    import numpy as np
    import sys, os
    import matplotlib.pyplot as plt
    from matplotlib import rc
    import datetime
    from statistics_from_dumpfiles import write_to_file
    sys.exit(main(sys.argv))