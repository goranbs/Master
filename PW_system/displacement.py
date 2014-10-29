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
    
    
def initialize(Natoms,matrix0,types,Nx,Nz,xmin,xmax):
    '''
    initialize(matrix = {}, types = [], Nx=#bins in x dir. Nz=#bins in z dir, xmin E [0,1], xmax E [0,1])\n
    Where xmin < xmax\n
    Takes the initial statefile information and bins the atoms
    according to their position in space.
    It returns an array of bins that contain the atoms located at the bins
    position in the system.\n
    Returns:\n
    binsx, binsz\n
    which are wo three dimentional matrices, where an atom j in bin i is accessed by:
    binsx[i][j] = [<id>, <xs>, <ys>, <zs>]
    '''


    binsx = []
    binsz = []
    
    for Bin in range(Nx): # fill with empty bins
        binsx.append([])
    for Bin in range(Nz):
        binsz.append([])
    
    for i in range(Natoms):
        if (matrix0['type'][i] in types):
            x = matrix0['xs'][i]
            if ((x < xmin) and (x < xmax)):
                index_z = int(round(matrix0['zs'][i]*(Nx-1))) # index of z bin in system
                index_x = int(round(x*(Nz-1)))                # index of x bin in system
                atom = [matrix0['id'][i],x,matrix0['ys'][i],matrix0['zs'][i]]

                binsx[index_x].append(atom)
                binsz[index_z].append(atom)
    '''    
    for i in range(Nx):
        print "#------------------------------------------------------#"
        print "Containernr: %g " % i
        for j in range(len(binsz[i])):
            print "atom: %g :: x=%g  y=%g  z=%g" % (binsz[i][j][0],binsz[i][j][1],binsz[i][j][2],binsz[i][j][3])

    '''
    return np.array(binsx),np.array(binsz)
    

def displacement(t, initial_binsx, initial_binsz, matrix, Natoms, types, Nx, Nz, xmin, xmax):
    '''
    displacementprofile() shuld take a system state at time t (t1), and a system state at time t+1 (t2).
    It creates bins in the x-y plane (along z-axis) and in the y-z plane (along x-axis) where it places the 
    atoms according to their location in space.
    Every atoms displacement from t1 to t2 is calculated only if it is still located in the bin it was in at
    time t1. If not. The atom gets a new initial position in the new bin, and it will therefore not contribute
    to any displacement.\n
    t=time, osv.
    '''

    binsx = []
    binsz = []
    still_in_binx = []
    still_in_binz = []
    atoms_to_removex = []
    atoms_to_removez = []
    
    for Bin in range(Nx): # fill with empty bins
        binsx.append([])
        still_in_binx.append([])
        atoms_to_removex.append([])
    for Bin in range(Nz):
        binsz.append([])
        still_in_binz.append([])
        atoms_to_removez.append([])
    
    for i in range(Natoms):
        if (matrix['type'][i] in types):
            x = matrix['xs'][i]
            if ((x < xmin) and (x < xmax)):
                index_z = int(round(matrix['zs'][i]*(Nx-1))) # index of z bin in system
                index_x = int(round(x*(Nz-1)))                # index of x bin in system
                atom = [matrix['id'][i],x,matrix['ys'][i],matrix['zs'][i]]

                binsx[index_x].append(atom)
                binsz[index_z].append(atom)

    # Now the atoms in the next timestep is added to bins. We must further go
    # through the bins, find the atoms that are still in the bin and calculate
    # the displacement. If an atom has moved to another bin, then  we have to 
    # delete that atom from that bin, and add it to the new bin that it is
    # located within.

    
    msd_x = np.zeros((Nx,1))
    msd_z = np.zeros((Nz,1))
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
                    still_in_binx[j].append(b)
                    count += 1 # number of atoms that contribute to the msd
                    
                    
            if (count != 0):
                msd_x[j] = msd_x[j]/count # mean square displacement in bin j.

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
                    still_in_binz[j].append(b)
                    count += 1 # number of atoms that contribute to the msd

            if (count != 0):
                msd_z[j] = msd_z[j]/count # mean square displacement in bin j.
    
    # Moving atoms to their correct bin, if they are not already located there!
    # and also make a list of atoms to be removed from bins they were located in before
    for i in range(Nx):
        shit = False
        if (i == (Nx-1)):
            shit = True # we are at the last bin
        for j in range(len(initial_binsx[i])):
            if (initial_binsx[i][j][0] not in still_in_binx[i]):
                atomtofind = initial_binsx[i][j][0]
                foundatom = False
                # then it has moved to another box/bin.
                # 1) find the new bin and append the atom info to this
                # search in the closest bins!
                for m in range(Nx/2): # go through alternating (i+1, i-1, i+2, i-2,...) bins:
                    alpha = 1
                    if (m%2 == 0 and i > 0):
                        alpha = -1
                    if (shit):
                        alpha = -1
                    index = i + alpha
                    for k in range(len(binsx[index])):
                        if (binsx[index][k][0] == atomtofind):
                            x = binsx[index][k][1]
                            index_x = int(round(x*(Nx-1)))
                            np.append(initial_binsx[index_x],binsx[index_x][k])
                            #initial_binsx[index_x].append(binsx[index][k])
                            #add = initial_binsx[i][j][0]
                            # 2) add indexes of atoms to be removed from the original bins
                            atoms_to_removex[i].append(j)
                            foundatom = True
                            '''
                            print "##--------------------------------------------------------------------##"
                            print "##-- Adding/Removing || alpha= %g || i = %g || mod = %g || shit= %g --##" % (alpha, i, m%2, shit)
                            print "Moved   : %g from bin %g to %g" % (add,i,index_x)
                            '''
                        if (foundatom == True):
                            break # if we have found the atom, then break the loop!
                    if (foundatom == True):
                        break     # if we have found the atom, then break the loop!
    
    for i in range(Nz):
        shit = False
        if (i == (Nz-1)):
            shit = True # we are at the last bin
        for j in range(len(initial_binsz[i])):
            if (initial_binsz[i][j][0] not in still_in_binz[i]):
                atomtofind = initial_binsz[i][j][0]
                foundatom = False
                # then it has moved to another box/bin.
                # 1) find the new bin and append the atom info to this
                # search in the closest bins!
                for m in range(Nz/2): # go through alternating (i+1, i-1, i+2, i-2,...) bins:
                    alpha = 1
                    if (m%2 == 0 and i > 0):
                        alpha = -1
                    if (shit):
                        alpha = -1
                    index = i + alpha
                    for k in range(len(binsz[index])):
                        if (binsz[index][k][0] == atomtofind):
                            z = binsz[index][k][3]         # the z-value!
                            index_z = int(round(z*(Nx-1))) # index of z bin in system
                            np.append(initial_binsz[index_z],binsz[index_z][k])
                            #initial_binsz[index_z].append(binsz[index_z][k])
                            #add = initial_binsx[i][j][0]
                            # 2) add indexes of atoms to be removed from the original bins
                            atoms_to_removez[i].append(j)
                            foundatom = True
                            '''
                            print "##--------------------------------------------------------------------##"
                            print "##-- Adding/Removing || alpha= %g || i = %g || mod = %g || shit= %g --##" % (alpha, i, m%2, shit)
                            print "Moved   : %g from bin %g to %g" % (add,i,index_z)
                            '''
                        if (foundatom == True):
                            break # if we have found the atom, then break the loop!
                    if (foundatom == True):
                        break     # if we have found the atom, then break the loop!    
    
    print initial_binsx[1][:]
    print type(initial_binsx)
    print type(initial_binsx[1][1])
    print atoms_to_removex[1]
    for i in range(Nx):
        new_array = np.delete(initial_binsx[i],atoms_to_removex[i],0)
        initial_binsx[i] = new_array
    #print initial_binsx
  
    for i in range(Nz):
        new_array = np.delete(initial_binsz[i],atoms_to_removez[i],0)
        initial_binsz[i] = new_array
        
    print "#------------- Time %g ---------------#" % t
    #print msd_x
    #print "---------------------------------------"
    #print msd_z
    
    return initial_binsx, initial_binsz  # return the updated initial bins
        

def main():
    
    path = "/home/goran/lammps-28Jun14/examples/water_portlandite_system/npt_run_and_energyminimization/statefiles"
    
    filenames,time = gothroughfiles(path) # get filenames in location "path"
    start = 300
    stop = -1
    filenames = filenames[start:stop]
    #length = len(filenames)
    #length2 = len(filenames)
    #print "before: %g   after: %g " % (length, length2)
    
    name = os.path.join(path,filenames[0])

    t,Natoms,system_size,matrix,readstructure,entries,types = readfile(name)

    xmin = 0.25
    xmax = 0.75    
    Nx = Nz = 20
    Types = [4] # oxygen
    binsx, binsz = initialize(Natoms,matrix,Types,Nx,Nz,xmin,xmax)
    '''    
    dispx = [], dispz = [], timearray = np.zeros(len(filenames))

    for i in range(Nx): # for all bins, create a timearray
        dispx.append(timearray)
    for j in range(Nz):
        dispz.append(timearray)
    '''
    for name in filenames[1:]:
        filename = os.path.join(path,name)
        t,Natoms,system_size,matrix,readstructure,entries,types = readfile(filename)
        binsx,binsz = displacement(t,binsx,binsz,matrix,Natoms,types,Nx,Nz,xmin,xmax)
        
        
        
        
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