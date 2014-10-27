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
        
        
def gothroughfiles(path):
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


def displacementprofile(t1, t2, matrix1, matrix2, Natoms, sytem_size, readstructure, entries, types):
    '''
    displacementprofile() shuld take a system state at time t (t1), and a system state at time t+1 (t2).
    It creates bins in the x-y plane (along z-axis) and in the y-z plane (along x-axis) where it places the 
    atoms according to their location in space.
    Every atoms displacement from t1 to t2 is calculated only if it is still located in the bin it was in at
    time t1. If not. The atom gets a new initial position in the new bin, and it will therefore not contribute
    to any displacement.
    '''

    #1 define area of boxes. x E [0.25, 0.75], y E [0,1], z E [0,1]
    # number of boxes in each dim.

    Nx = 20 # Number of boxes in x direction
    Nz = 20 # Number of boxes in z direction
    
    binsx1 = np.zeros((Nx,1)), binsz1 = np.zeros((Nz,1)) # contains Nx and Nz bins that contains a number of atoms in that region.
    binsx2 = np.zeros((Nx,1)), binsz2 = np.zeros((Nz,1)) # contains Nx and Nz bins that contains a number of atoms in that region.

    for Bin in range(Nx): # fill the arrays with empty bins that will contain the atoms.
        binsx1[Bin] = []
        binsx2[Bin] = []
    for Bin in range(Nz):
        binsz1[Bin] = []
        binsz2[Bin] = []
        
    for i in range(Natoms): # run through all the atoms, and use their positions to place them in their bins.
        #atoms[i] = Atom(matrix1['id'][i])
        #atoms[i].set_all(matrix1['mol'][i],matrix1['type'][i],matrix1['xs'][i],matrix1['ys'][i],matrix1['zs'][i],matrix1['vx'][i],matrix1['vy'][i],matrix1['vz'][i],matrix1['fx'][i],matrix1['fy'][i],matrix1['fz'][i])
        if (matrix1['type'][i] in types): # if the atom is one that we are following, go on:
            if (matrix1['xs'] > xmin and matrix21['xs'] < xmax):
                index_z = int(round(matrix1['zs'][i]*(Nz-1)) # hight in system
                index_x = int(round(matrix1['xs'][i]*(Nx-1)) # x pos in system
                binsz1[index_z].append(matrix1['id'][i])     # add atom id to the bin it's located within
                binsx1[index_x].append(matrix1['id'][i])     # add atom id to the bin it's located within
            
        # do the same for the next timestep:
        if (matrix2['type'][i] in types):
            if (matrix2['xs'][i] > xmin and matrix2['xs'][i] < xmax):
                index_z = int(round(matrix2['zs'][i]*(Nz-1))
                index_x = int(round(matrix2['xs'][i]*(Nx-1)) 
                binsz2[index_z].append(matrix2['id'][i])
                binsx2[index_x].append(matrix2['id'][i])  
    
    # Now we have to go through the bins and calculate the msd for every bin.
    
            
        

def main():
    
    path = "/home/goran/lammps-28Jun14/examples/water_portlandite_system/npt_run_and_energyminimization/statefiles"
    
    filenames,time = gothroughfiles(path) # get filenames in location "path"
    start = 300
    stop = -1
    filenames2 = filenames[start:stop]
    length = len(filenames)
    length2 = len(filenames2)
    print "before: %g   after: %g " % (length, length2)
    for name in filenames2:
        t,Natoms,system_size,matrix,readstructure,entries,types = readfile(name)
        
        
        
class Atom:
    
    def __init__(self,atom_index):
        self.atom_index = atom_index
        self.pos = [None,None,None]           # holds position of atom
        self.inititial_pos = [None,None,None] # holds the initial position of the particle
        self.vel = [None,None,None]           # holds velocities
        self.forc = [None,None,None]          # holds forces on atom
        self.mol = None                       # holds molecule id
        self.TYPE = None                      # holds type id
        
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
        self.inititial_pos[0] = x
        self.inititial_pos[1] = y
        self.inititial_pos[2] = z
        
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
    
    def get_square_displacement(self):
        dx = self.pos[0] - self.initial_pos[0]
        dy = self.pos[1] - self.initial_pos[1]
        dz = self.pos[2] - self.initial_pos[2]
        return (dx*dx + dy*dy + dz*dz)
        

if (__name__ == "__main__"):
    import numpy as np
    import os
    import matplotlib.pyplot as plt
    main()