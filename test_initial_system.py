# -*- coding: utf-8 -*-
"""
Created on Tue Dec  2 13:36:36 2014

@author: goran
"""

def main():
    path = '/home/goran/lammps-28Jun14/examples/Abel_runs/x_moltemplates'
    filename = 'coto_evap_sys_standing.data'
    
    fullname = os.path.join(path,filename)
    #print fullname
    ofile = open(fullname,'r')

    condition = True
    while (condition == True): #read until...
        line = ofile.readline()
        #print line
        if ('Atoms' in line):
            condition == False
            break
        
    
    ofile.readline() # mpty line
    condition = True
    total_charge = 0
    while (condition == True):
        line = ofile.readline()
        values = line.split()
        if (len(values) <= 1):
            condition = False
            break
        ID = int(values[0])
        bond = int(values[1])
        mol = int(values[2])
        try:
            charge = float(values[3])
            # gidder ikke ta resten, siden det bare er charge jeg er interessert i nå :-)
            total_charge += charge
        except:
            continue
        
    
    print total_charge
                
    
if (__name__ == "__main__"):
    import os
    print "Running main..."
    main()