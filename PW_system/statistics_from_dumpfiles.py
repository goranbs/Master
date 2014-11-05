# -*- coding: utf-8 -*-
"""
Created on Wed Nov  5 08:46:09 2014

@author: goran

This program should read and process dump files located in the dump directories
of the master project.
These files contain averaged datas. Filetypes are:

Adump.NDens        -  files that contain system density information 
Adump.*_RDF_       -  files that contain Radial Density information
Adump.*Scalarinfo  -  files that contain scalar values like bondlengths, ROG, Temperature..

This program should search for these types of files, or you may specify the files yourself.
The absolute path of where the files are located has to be given!
Then the program must exrtact the information from the files, process it, find averages and
standard deviations, and finally plot the results.
"""
def gothroughfiles(path,arg=None):
    '''
    Go through stystem state files in "path".
    '''
    
    filenames = []
    
    if arg is not None:
        for name in os.listdir(path):
            if (arg in name):
                filenames.append(name)
    else:
      print "arg shuld be specified, now the function will return all filenames in the given path!"
      for name in os.listdir(path):
          filenames.append(name)
                    
    sortedfilenames = sorted(filenames, key = lambda x: x[:-4])
    
    return sortedfilenames





def main():
    generalpath = "/home/goran/lammps-28Jun14/examples/water_portlandite_system"    

    wp_path = "npt_run_and_energyminimization/dump/"
    arg = "Adump"
    path = os.path.join(generalpath,wp_path)
    filenames = gothroughfiles(path,arg)
    
    '''
    There are different types of files in filename that we are going to do different operation on
    so we should find these types of files, and send them to the correct functions.    
    '''

    dens_files = []
    RDF_files = []
    ScalarInfo_files = []
    Energy_files = []
    
    for name in filenames:
        if ("Dens" in name):
            dens_files.append(name)
        if ("RDF" in name):
            RDF_files.append(name)
        if ("ScalarInfo" in name):
            ScalarInfo_files.append(name)
        if (("Energy" or "energy") in name):
            Energy_files.append(name)
    
    
    
    

if (__name__ == "__main__"):
    import numpy as np
    import matplotlib.pyplot as plt
    import os
    main()


