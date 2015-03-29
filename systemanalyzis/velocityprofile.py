# -*- coding: utf-8 -*-
"""
Created on Wed Oct 15 10:21:40 2014

@author: goran

Takes one statefile from LAMMPS output. 
Read initial line and identify index, molecule, type, position, velocity
and forces acting on the particle.
Extract information from file.

The system size is in a scaled range, so that x,y,z E [0,1]

Tasks:

1. Bin particles in central region of the tube in z-direction, average the
   velocities of the particles and create a velocity profile. This will
   indicate the structure of the flow inside the tube.

2. Estimate the density of the sytem in regions

3. Bin particles in left and right region of the system (water reservoairs)
   Find a way to calculate the pressure in these regions.

4. Later work, do these calculations for multiple files, and generate a 
    time evolution of the flow, density and the pressure.
"""

def readfile(filename):
    '''
    Returns timestep, number of atoms, system sizes and
    a matrix that contains information about the atoms
    Returns:
    t,Natoms,system_size,matrix,readstructure, entries,types
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


def velocityprofile(time,Natoms,system_size,matrix,entries,save_fig=None,showplot=False):
    '''
    velocityprofile should divide the central part of the system into bins,
    spanning across the shaft/tube in its z-direction. Find what atoms that
    belong to each bin, and average their velocity.
    returns:
    nothing, it returns nothing! only plots a velocityprofile for the given timestep
    '''

    dx = 0.25     # remember that the system is scaled to be 1 wide, long and high.
    x = 0.5       # middle of the system
    oxygen = 4    # atom type
    nbins = 20
    bins = []

    for i in range(nbins):
        abin = []
        bins.append(abin)
    
    vx = np.zeros(nbins)
    vy = np.zeros(nbins)
    vz = np.zeros(nbins)
    natoms_in_bins = np.zeros(nbins)
    
    for i in range(Natoms):
        if (matrix['type'][i] == oxygen): # !!! We get all the oxygen atoms, only water!
            if ((matrix['xs'][i] < (x + dx)) and (matrix['xs'][i] > (x - dx))):
                # then the values lie in the desired x-range
                indx = int(round(matrix['zs'][i]*(nbins-1)))
                #atom = Atom(matrix['id'])
                #atom.set_all(matrix['mol'],matrix['type'],matrix['xs'],matrix['ys'],matrix['zs'],matrix['vx'],matrix['vy'],matrix['vz'],matrix['fx'],matrix['fy'],matrix['fz'])
                #bins[indx].append(atom)
                natoms_in_bins[indx] += 1
                vx[indx] += matrix['vx'][i]
                vy[indx] += matrix['vy'][i]
                vz[indx] += matrix['vz'][i]
    
    #print natoms_in_bins
    for i in range(nbins):
        if (natoms_in_bins[i] != 0):
            vx[i] = vx[i]/natoms_in_bins[i]
            vy[i] = vy[i]/natoms_in_bins[i]
            vz[i] = vz[i]/natoms_in_bins[i]
    
    zlo = system_size[4]
    zhi = system_size[5]
    zax = np.linspace(zlo,zhi,nbins)

    # now plotting:

    factor = 10**5  # conversion factor from [Aangstrom/fsec] to [m/s]
    fig1 = plt.figure()
    plt.plot(vx*factor,zax,'b-')
    plt.hold(True)
    plt.plot(vz*factor,zax,'r--')
    plt.plot(vy*factor,zax,'y--')
    plt.hold(False)
    plt.title('velocity along the z-axis at timestep t=%g' % time)
    plt.xlabel('v [m/s]')
    plt.ylabel('z [Aangstrom]')
    plt.legend(['v_x(z)','v_z(z)','v_y(z)'],loc='lower left')
    plt.show(showplot)
    
    if save_fig is not None:
        if save_fig is not False:
            filename = 'velocityprofile_' + str(time) + '.png'
            plt.savefig(filename,format='png')
            print " saving plot 1: Velocity plot"
    plt.close(fig1)
        
    
    fig2 = plt.figure()
    plt.plot(natoms_in_bins,np.linspace(0,nbins,nbins),'b--o')
    plt.title('Number of oxygen atoms counted in bin along z-axis, \nand then hopefully the number of water molecules. t=%g ps' % time)
    plt.xlabel('number of oxygen atoms')
    plt.ylabel('bin')
    plt.legend(['oxygen atoms'],loc='upper right')
    plt.show(showplot)
    
    if save_fig is not None:
        if save_fig is not False:
            filename = 'oxygendistribution_' + str(time) + '.png'
            plt.savefig(filename,format='png')
            print " saving plot 2: oxygen distribution"
    plt.close(fig2)
    

def density(time,Natoms,system_size,matrix,entries,save_fig=None,showplot=False):
    '''
    Calculate the density of the water in bulk and inside the tube.
    Density is only mass per unit volume, so we can again count the number 
    of oxygen atoms present in some volum, multiply with the mass of a
    water molecule, and divide by the volume.
    returns:
    it returns nothing, just plots the density profile at the given timestep
    '''
    oxygen = 4       # index type of oxygen in water
    x = 0.5          # Center of system 
    dx_tube = 0.2    # Create box in middle of system
    dx_bulk = 0.15   #  Bulk box
    
    nbins = 15
    natoms_tube_dist = np.zeros(nbins)
    natoms_bulk = 0
    
    for i in range(Natoms):
        if (matrix['type'][i] == oxygen): # only count oxygen atoms
            
            if (matrix['xs'][i] < (x + dx_tube) and matrix['xs'][i] > (x - dx_tube)):
                # divide into layers in z-direction to find a density distribution
                indx = int(round(matrix['zs'][i]*(nbins-1)))
                natoms_tube_dist[indx] += 1
            
            if ((matrix['xs'][i] < dx_bulk) or (matrix['xs'][i] > (1 - dx_bulk))):
                # The particle is in lower or upper region of the system
                # e.g is in bulk.
                natoms_bulk += 1
    
    #print natoms_tube_dist
    #print natoms_bulk
#     need to calculate the volumes of the boxes!
    
    factor = 1 #10**(-7)    # Conversion factor from Aangstrom to centimeter 
    factor1 = 10**(-21)      # conversion Aangstrom**3 to cm**3
    #size_x = (system_size[1] - system_size[0])*factor   # [A]
    size_y = (system_size[3] - system_size[2])*factor    # [A] system length y-dir
    size_z = (system_size[5] - system_size[4])*factor    # [A] system length z-dir
    dz = (size_z/nbins)*factor                           # [A] bin hight
    dx_tube = dx_tube*factor                             # [A]
    dx_bulk = dx_bulk*factor                             # [A]
    
    volume_bulk = size_y*size_z*(2*dx_bulk)*factor1              # [cm**3]
    volume_tube = size_y*(2*dx_tube)*dz*factor1                  # [cm**3]
    
    Na = 6.0221413*10**(23)            # Avogadro constant [molecules/mol]
    M_oxygen = 15.9994                 # [u] = [g/mol]
    M_hydrogen = 1.008                 # [u] = [g/mol]
    M_water = M_oxygen*(2*M_hydrogen)  # Molar mass of water [g/mol]
    #mass_water_molecule = M_water/Na   # Mass of one water moecule [g]     
    mass_water_molecule = 2.99153*10**(-23)
    
    density_water_bulk = mass_water_molecule*natoms_bulk/volume_bulk # [g/m**3]  
    
    zax = np.linspace(system_size[4], system_size[5],nbins)
    density_water_tube = np.zeros(nbins)
    f = mass_water_molecule/volume_tube
                      
    average = 0
    n_in_average = 0                 
    for i in range(nbins):
        if (natoms_tube_dist[i] != 0):
            density_water_tube[i] = f*natoms_tube_dist[i]
            average += density_water_tube[i]
            n_in_average += 1
    
    average = average/n_in_average
    avg_array = np.ones(len(density_water_tube))
    for i in range(len(density_water_tube)):
        avg_array[i] = average
    
    # plotting of density porfile in tube:
    fig3 = plt.figure()
    plt.plot(density_water_tube,zax,'bo--')
    plt.hold(True)
    plt.plot(avg_array,zax,'r--')
    plt.title('Density profile for water between the plates.\nAt time %g ps' % time)
    plt.xlabel('Density [g/cm**3]')
    plt.ylabel('z [Aangsrom]')
    plt.legend(['rho_{water}','average'],loc='lower right') 
    plt.show(showplot)

    if save_fig is not None:
        if save_fig is not False:
            filename = 'densityprofile_water_ ' + str(time) + '.png'
            plt.savefig(filename, format='png')
            print " saving plot 3: density profile of water"
    plt.close(fig3)
    
    print "----------------------------------------------"    
    print "Density of water in bulk, rho=%.3f [g/cm**3]" % density_water_bulk
    print "Number of molecules counted: %g " % natoms_bulk 
    print "Mass of one water molecule = %g " % mass_water_molecule
    print "----------------------------------------------"
    print density_water_tube
    
def density_system_CO2(time,Natoms,system_size,matrix,entries,masses,save_fig=None,showplot=False):
    '''
    density_system_CO2 prints the density of a CO2 system.
    Use readfile(filename) to get system information before using this function.
    returns:
                density
    '''
    Na = 6.0221413*10**(23)                   # Avogadro constant [molecules/mol]
    conversionfactor = 10**(-24)               # [A**3 -> cm**3]
    
    size_x = (system_size[1] - system_size[0]) # [A] system length x-dir
    size_y = (system_size[3] - system_size[2]) # [A] system length y-dir
    size_z = (system_size[5] - system_size[4]) # [A] system length z-dir
    vol_Aa = size_x*size_y*size_z              # [A**3]
    #vol_cm = vol_Aa*conversionfactor           # [cm**3]        
    mass_CO2 = masses['C'] + 2*masses['O']     # Mass CO2 [au] = [g/mol]    
    
    # count number of C particles in the system:
    carbon = 1
    nmolecules = 0
    for i in range(Natoms):
        if (matrix['type'][i] == carbon):
            nmolecules += 1
            
    conv = 0.60221413                         # Na/conversionfactor [particles/mol]
    tot_mass_CO2 = nmolecules*mass_CO2        # [particles*g/mol]
    density = (tot_mass_CO2/vol_Aa)*conv      # [g/cm**3]
    #dens_kg_m3 = density*10**(3.0)            # [kg/m**3]
    
    '''
    #print size_x, size_y, size_z
    print "-------- Density of system ----------"
    print " # molecules = %g " % nmolecules
    print " mass CO2    = %g [g/mol]" % mass_CO2
    print " volume      = %g [cm**3]" % vol_cm 
    print " rho         = %g [g/cm**3]" % density
    print " rho         = %g [kg/m**3]" % dens_kg_m3
    print " number density = %g [particles/Aa] " % (nmolecules/vol_Aa)
    print "-------------------------------------"
    '''
    return density
    

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
    if (arg is not None):
        for name in os.listdir(path):
            print name
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
        for name in os.listdir(path):
            if (name[-4:] == ".txt"):
                    filenames.append(name)
                    
    time = sorted(numbers)
    sortedfilenames = sorted(filenames, key = lambda x: x[:-4])
    
    return sortedfilenames,time
    
def get_velocityprofile(t,Natoms,system_size,matrix,entries,types):
    '''
    get_velocityprofile takes a timevalue t, system_size which is a list containing
    the system size [xlo xhi ylo yhi zlo zhi], a library matrix element that contains
    information about the atoms in the system. entries holds known entries in matrix
    types contain known types of atoms in the system.
    get_velocityprofile should return two arrays that desctibes the velocityprofile
    of the system at a given timestep t.
    returns:
            vx,vy,vz,nbins
    '''

    dx = 0.2     # remember that the system is scaled to be 1 wide, long and high.
    x = 0.5       # middle of the system
    oxygen = 4    # atom type that we are interested in knowing the velocity distribution about.
    nbins = 15    # number of bins that we divide the system hight into
    
    vx = np.zeros((nbins,1))
    vy = np.zeros((nbins,1))
    vz = np.zeros((nbins,1))
    natoms_in_bins = np.zeros(nbins)
    
    for i in range(Natoms):
        if (matrix['type'][i] == oxygen): # !!! We get all the oxygen atoms, only water!
            if ((matrix['xs'][i] < (x + dx)) and (matrix['xs'][i] > (x - dx))):
                # then the values lie in the desired x-range
                indx = int(round(matrix['zs'][i]*(nbins-1)))
                natoms_in_bins[indx] += 1
                vx[indx] += matrix['vx'][i] # add velocities of the oxygen atoms into the correct bin
                vy[indx] += matrix['vy'][i]
                vz[indx] += matrix['vz'][i]
                #print matrix['vx'][i]
                if (matrix['vx'][i] == matrix['vy'][i] and  matrix['vx'][i] == matrix['vz'][i]):
                    print "vx = vy something is wrong. go to line 402!!"
    
    for i in range(nbins):
        # average the velocities of the atoms present in the bins:
        if (natoms_in_bins[i] != 0):
            vx[i] = vx[i]/natoms_in_bins[i]
            vy[i] = vy[i]/natoms_in_bins[i]
            vz[i] = vz[i]/natoms_in_bins[i]
    
    # Velocities are given in [Aangstrom/fsec] from LAMMPS, but we leave this to the plotting function for now

    '''
    plt.figure()
    plt.plot(np.linspace(0,nbins,nbins),vx,'r--')
    plt.hold(True)
    plt.plot(np.linspace(0,nbins,nbins),vy,'b--*')
    plt.plot(np.linspace(0,nbins,nbins),vz,'y--d')
    plt.hold(False)

    plt.plot(np.linspace(0,nbins,nbins),vx)
    plt.title('PLOTTIPLOTT')
    plt.show(True)
    '''
    
    return vx,vy,vz,nbins

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
    
def main():
    
    system = "CO2"
    #system = "water"
    #system = "PW"   
    if (system == "CO2"):
        
        savefile = False
        showplot = False   # there is no plotting funtion here!
        
        '''
        generalpath = "/home/goran/lammps-28Jun14/examples/water_portlandite_system"
        path_from_PW_system = "CarbonDioxide/statefiles/"
        '''
        generalpath = "/home/goran/lammps-28Jun14/examples/Abel_runs/carbondioxide"
            
        masses = {'C' : 12.0107, 'O' : 15.9994, 'H' : 1.00794, 'Ca' : 40.078}
        filename = "dump.CO2_res.200000.txt"
        
        dens = []; t = []; tc = 1.0/(10000)
        filenames = gothroughfiles(generalpath,arg='CO2.')  
        print filenames[1]
        for name in filenames[:20]:    
            thisfile = os.path.join(generalpath,name)
            time,Natoms,system_size,matrix,readstructure,entries,types = readfile(thisfile)
            density = density_system_CO2(time,Natoms,system_size,matrix,entries,masses,savefile,showplot)
            dens.append(density), t.append(float(time)*tc)
        
        plt.figure()
        plt.plot(t,dens)
        plt.xlabel('time [ps]'), plt.ylabel('density [g/cm^3]'), plt.legend(['rho(t)'],loc='upper left')
        plt.show(True)
    
    if (system == "PW"):
        generalpath = "/home/goran/lammps-28Jun14/examples/water_portlandite_system"
        path_from_PWsystem = "npt_run_and_energyminimization/statefiles/"
        only_nve = True
    
        path = os.path.join(generalpath,path_from_PWsystem)
        print "--------------------------------------------------------"
        print " PATH where files are found:"
        print path
        print "--------------------------------------------------------"    
        
        filenames, time = gothroughfiles(path)
        
        filestructure = "dump.PW_nve"
        #filestructure = "dump.PW_npt"
        length = len(filestructure)
        if (only_nve == True):
            filenames_to_use = []
            for filename in filenames:
                # if filename contains "nve" in the filename: use these filenames:
                if (filename[0:length] == filestructure):
                    filenames_to_use.append(filename)
                
            filenames = filenames_to_use
    
        savefile_velocityprofile = False
        showplot_velocityprofile = True #False
        savefile_density = False
        showplot_density = True #False
        savefile_sysdens = False
        showplot_sysdens = True
        
        vel_x = []
        vel_y = []
        vel_z = []
        len_filenames = len(filenames)
        times = np.zeros(len_filenames)
        ii = 0
        for filename in filenames:
            print "processing file: %s . Number %g of total %g" % (filename,(ii+1),len_filenames)
            filename = os.path.join(path,filename)
            t,Natoms,system_size, matrix, readstructure, entries, types = readfile(filename)
            vx,vy,vz,nbins = get_velocityprofile(t,Natoms,system_size,matrix,entries,types)
            vel_x.append(vx[:])
            vel_y.append(vy[:])
            vel_z.append(vz[:])
            '''
            plt.plot(np.linspace(0,nbins,nbins),vz[:])
            plt.show(True)
            '''
            times[ii] = t
            ii += 1
    
        # Now we have averaged over time the velocity profile of the system
        # testplotting:
     
        zlo = system_size[4]
        zhi = system_size[5]
        zax = np.linspace(zlo,zhi,nbins)
        len_vel = len(vel_x)
        len_vx = len(vel_x[0])
        
        avg_vel_x = np.zeros((len_vx,1))
        avg_vel_y = np.zeros((len_vx,1))
        avg_vel_z = np.zeros((len_vx,1))
        max_val = 0
        for i in range(len_vel):
            for j in range(len(vel_x[i])):
                velx = vel_x[i][j]
                vely = vel_y[i][j]
                velz = vel_z[i][j]
                avg_vel_x[j] += velx
                avg_vel_y[j] += vely
                avg_vel_z[j] += velz
                
                if (velx > max_val):
                    max_val = velx
                if (vely > max_val):
                    max_val = vely
                if (velz > max_val):
                    max_val = velz
                #print vel_z[i][j]
     
    
        for j in range(len(avg_vel_x)):
            avg_vel_x[j] = avg_vel_x[j]/len_vel
            avg_vel_y[j] = avg_vel_y[j]/len_vel
            avg_vel_z[j] = avg_vel_z[j]/len_vel
    
     
        Title = 'Velocity profile through tube.\nAverage over %g number of timesteps' % len(time)
        name = 'velocityprofile_npt.png'
        legends = ['vx(z)','vy(z)','vz(z)']
        xlabel = 'system width [Aangstrom]'
        ylabel = 'velocity [Aangstrom/fsec]'
        
        fig = plt.figure()        
        plt.plot(zax,avg_vel_x,'b-')
        plt.hold(True)
        plt.plot(zax,avg_vel_y,'r--')
        plt.plot(zax,avg_vel_z,'y--d')
        line = np.linspace(min(avg_vel_z),max_val,len_vx)
        zlo_surf = 7.4   # [Aangstrom]
        zhi_surf = 37.0  # [Aangstrom]
        surf_lo = np.linspace(zlo_surf,zlo_surf,len_vx)
        surf_hi = np.linspace(zhi_surf,zhi_surf,len_vx)
        plt.plot(surf_lo,line,'k--')
        plt.plot(surf_hi,line,'k--')
        legends.append('surface')
        plt.hold(False)
        plt.title(Title)
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
        plt.legend(legends,loc='upper right')
        plt.savefig(name,format='png')
        plt.show(False)
        plt.close(fig)
    
        
        contents = 'zax[Aangstrom] vx(z) vy(z) vz(z)'
        toptext = 'Velocity profile for water-portalndite system. Generated from %g timesteps. from files: %s ' % (len_vel, filestructure)
        writetofile('velocityprofile_nve.txt',toptext,contents,avg_vel_x,avg_vel_y,avg_vel_z)
            
        #----------------- End Of PW system analyze -----------------------------
        


class Atom:
    
    def __init__(self,atom_index):
        self.atom_index = atom_index
        self.pos = [None,None,None]  # holds position of atom
        self.vel = [None,None,None]  # holds velocities
        self.forc = [None,None,None] # holds forces on atom
        self.mol = None              # holds molecule id
        self.TYPE = None               # holds type id
        
    def set_all(self,MOL,TYPE,x,y,z,vx,vy,vz,fx,fy,fz):
        self.TYPE = TYPE
        self.mol = MOL
        self.pos[0] = x
        self.pos[1] = y
        self.pos[2] = z
        self.vel[0] = vx
        self.vel[1] = vy
        self.vel[2] = vz
        self.forc[0] = fx
        self.forc[1] = fy
        self.forc[2] = fz
        
    def set_Type(self,TYPE):
        self.TYPE = TYPE
    
    def set_molID(self,molecule_id):
        self.mol = molecule_id
    
    def set_position(self,x,y,z):
        self.pos[0] = x
        self.pos[1] = y
        self.pos[2] = z
    
    def set_velocity(self,vx,vy,vz):
        self.vel[0] = vx
        self.vel[1] = vy
        self.vel[2] = vz
    
    def set_force(self,fx,fy,fz):
        self.forc[0] = fx
        self.forc[1] = fy
        self.forc[2] = fz
        
    def get_ID(self):
        return self.atom_index
    
    def get_position(self):
        return self.pos
    
    def get_velocity(self):
        return self.vel
        
        

if __name__ == "__main__":
    import matplotlib.pyplot as plt
    import numpy as np
    import os
    main()
       
       
    

