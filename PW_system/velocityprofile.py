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

    indexes = len(readstructure)
    for i in range(Natoms):
        line = infile.readline()
        values = line.split()
        for index in readstructure:
            matrix[entries[index]].append(values[index])
    
    for i in range(Natoms):
        for j in range(indexes):
            if matrix[entries[j]] is not None:
                if (entries[j] in ['id', 'mol', 'type']):
                    matrix[entries[j]][i] = int(matrix[entries[j]][i])
                else:
                    matrix[entries[j]][i] = float(matrix[entries[j]][i])
    '''
    for j in range(len(entries)):
        print matrix[entries[j]]
    '''
    return t,Natoms,system_size,matrix,readstructure, entries


def velocityprofile(t,Natoms,system_size,matrix,entries):
    '''
    velocityprofile should divide the central part of the system into bins,
    spanning across the shaft/tube in its z-direction. Find what atoms that
    belong to each bin, and average their velocity.
    '''
    import numpy as np

    dx = 0.25     # remember that the system is scaled to be 1 wide, long and high.
    x = 0.5       # middle of the system
    oxygen = 4    # atom type
    nbins = 15
    bins = []

    for i in range(nbins):
        abin = []
        bins.append(abin)
    
    vx = np.zeros(nbins)
    vy = np.zeros(nbins)
    vz = np.zeros(nbins)
    natoms_in_bins = np.zeros(nbins)
    
    for i in range(Natoms):
        if (matrix['type'][i] == oxygen): # !!! We get all the oxygen atoms, not only in water!
            if ((matrix['xs'][i] < (x + dx)) and (matrix['xs'][i] > (x - dx))):
                # then the values lie in the desired x-range
                indx = int(round(matrix['ys'][i]*(nbins-1)))
                atom = Atom(matrix['id'])
                atom.set_all(matrix['mol'],matrix['type'],matrix['xs'],matrix['ys'],matrix['zs'],matrix['vx'],matrix['vy'],matrix['vz'],matrix['fx'],matrix['fy'],matrix['fz'])
                bins[indx].append(atom)
                natoms_in_bins[indx] += 1
                vx[indx] += matrix['vx'][i]
                vy[indx] += matrix['vy'][i]
                vz[indx] += matrix['vz'][i]
    
    '''
    atom = Atom(10)
    atom.set_molID('test')
    bins[0].append(atom)
    print bins
    '''
    print natoms_in_bins
    for i in range(nbins):
        vx[i] = vx[i]/natoms_in_bins[i]
        vy[i] = vy[i]/natoms_in_bins[i]
        vz[i] = vz[i]/natoms_in_bins[i]
    
    zlo = system_size[4]
    zhi = system_size[5]
    zax = np.linspace(zlo,zhi,nbins)

    # now plotting:
    import matplotlib.pyplot as plt
    
    plt.figure()
    plt.plot(vx,zax,'b--')
    plt.hold(True)
    plt.plot(vz,zax,'r--')
    plt.plot(vy,zax,'y--')
    plt.hold(False)
    plt.title('velocity along the z-axis')
    plt.xlabel('v [Aangstrom/fsec]')
    plt.ylabel('z [Aangstrom]')
    plt.legend(['v_x(y)','v_z(y)','v_y(y)'],loc='lower left')
    
    factor = 10**5  # conversion factor from [Aangstrom/fsec] to [m/s]
    plt.figure()
    plt.plot(vx*factor,zax,'b--')
    plt.hold(True)
    plt.plot(vz*factor,zax,'r--')
    plt.plot(vy*factor,zax,'y--')
    plt.hold(False)
    plt.title('velocity along the z-axis')
    plt.xlabel('v [m/s]')
    plt.ylabel('z [Aangstrom]')
    plt.legend(['v_x(y)','v_z(y)','v_y(y)'],loc='lower left')
    
    plt.show(True)
    
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
    
    
def main():
    path = "/home/goran/Goran/PythonScripting/statefiles/"
    filename = "dump.PW_npt_minimized.20000.txt"
    
    t,Natoms,system_size, matrix, readstructure, entries = readfile(filename)

    velocityprofile(t,Natoms,system_size,matrix,entries)

if __name__ == "__main__":
    main()
       
       
    

