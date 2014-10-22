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
  
    #Natoms = 5  # !!!! TEMPORARY while debugging
    # now it's time to read in the datas:
    #for i in range(Natoms):
    #matrix = [ID, MOL, TYPE, XS, YS, ZS, VX, VY, VZ, FX, FY, FZ]
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
    '''
    plt.figure()
    plt.plot(vx,zax,'b--')
    plt.hold(True)
    plt.plot(vz,zax,'r--')
    plt.plot(vy,zax,'y--')
    plt.hold(False)
    plt.title('velocity along the z-axis at timestep t=%g ps' % time)
    plt.xlabel('v [Aangstrom/fsec]')
    plt.ylabel('z [Aangstrom]')
    plt.legend(['v_x(z)','v_z(z)','v_y(z)'],loc='lower left')
    '''
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
    Na = 6.0221413*10**(23)                   # Avogadro constant [molecules/mol]
    conversionfactor = 10**(-21)               # [A**3 -> cm**3]
    
    size_x = (system_size[1] - system_size[0]) # [A] system length x-dir
    size_y = (system_size[3] - system_size[2]) # [A] system length y-dir
    size_z = (system_size[5] - system_size[4]) # [A] system length z-dir
    vol_Aa = size_x*size_y*size_z              # [A**3]
    vol_cm = vol_Aa*conversionfactor           # [cm**3]        
    mass_CO2 = masses['C'] + 2*masses['O']     # Mass CO2 [au]
    CO2_mass_g = mass_CO2/Na                   # mass of CO2 molecule [g]
    
    # count number of C particles in the system:
    
    carbon = 1
    nmolecules = 0
    for i in range(Natoms):
        if (matrix['type'][i] == carbon):
            nmolecules += 1
            
    tot_mass_CO2 = nmolecules*CO2_mass_g      # [g]
    density = tot_mass_CO2/vol_cm             # [g/cm**3]
    dens_kg_m3 = density*10**(-3)              # [kg/m**3]
    
    #print size_x, size_y, size_z
    print "-------- Density of system ----------"
    print " # molecules = %g " % nmolecules
    print " mass CO2    = %g [g/mol]" % mass_CO2
    print " mass CO2    = %g [g]" % CO2_mass_g 
    print " volume      = %g [cm**3]" % vol_cm 
    print " rho         = %g [g/cm**3]" % density
    print " rho         = %g [kg/m**3]" % dens_kg_m3
    print " number density = %g [particles/Aa] " % (nmolecules/vol_Aa)
    print "-------------------------------------"
    return density
    

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
    
def main():
    generalpath = "/home/goran/lammps-28Jun14/examples/water_portlandite_system"
    path_from_PWsystem = None
    filename = "dump.PW_npt_minimized.20000.txt"
    #filename = "dump.water.191000.txt"
    #filename = "dump.PW_centralfolw.94200.txt"
    #filename = "dump.CO2.200000.txt"
    #filename = "dump.CO2.0.txt"
    
    t,Natoms,system_size, matrix, readstructure, entries, types = readfile(filename)

    savefile_velocityprofile = False
    showplot_velocityprofile = True #False
    savefile_density = False
    showplot_density = True #False
    savefile_sysdens = False
    showplot_sysdens = True
    
    dt = 0.1
    factor = 10**(-3)     # [fsec to ps]
    time = t*dt*factor    # [should get ps]
    
    masses = {'C' : 12.0107, 'O' : 15.9994, 'H' : 1.00794, 'Ca' : 40.078}
        
    velocityprofile(time,Natoms,system_size,matrix,entries,savefile_velocityprofile,showplot_velocityprofile)
    density(time,Natoms,system_size,matrix,entries,savefile_density,showplot_density)

    #density = density_system_CO2(time, Natoms,system_size,matrix,entries,masses,savefile_sysdens,showplot_sysdens)    


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
    main()
       
       
    

