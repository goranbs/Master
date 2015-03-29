# -*- coding: utf-8 -*-
"""
Created on Wed Mar 18 13:44:54 2015

@author: goran
"""

def read_logfile(path,logfilename):
    
    totalpath = os.path.join(path,logfilename)
    
    ofile = open(totalpath,'r')
    rec = 0
    timestep = None
    variable = None
    word = None
    dt = None
    t_value = None
    run = None

    while (run == None):
        line = ofile.readline()
        words = line.split()
        try:
            word = words[0]
            if (word == "run"):
                run = True
                print run
            dt = words[1]
            t_value = words[3]
            
        except:
            continue
    
        if ((word == 'variable') and (dt == 'dt')):
            print " timestep found!"
            timestep = t_value
        
    
    while (rec != "Memoryusageperprocessor"):
        line = ofile.readline()
        words = line.split() # split line on whitespaces
        #print line
        try:
            rec = words[0] + words[1] + words[2] + words[3]
            variable = words[0]
            dt = words[1]
            t_value = words[3]
        except:
            continue

        if ((variable == 'variable') and (dt == 'dt')):
            timestep = t_value
        
        
    print timestep
    
    line = ofile.readline()
    outputs = line.split()
    Noutputs = len(outputs)     # number of output variables
    #Outputs = []
    dictionary = {}
    for i in range(Noutputs):
        #new_list = []
        #Outputs.append(new_list)        # add empty list
        dictionary[outputs[i]] = []     # add dict list to item
    
        
    for line in ofile:
        values = line.split()
        try:
            float(values[0])
        except:
            #print "#--------------------------LAST LINE-----------------------------#"
            #print values
            break
        
        for i in range(Noutputs):
            #Outputs[i].append(float(values[i]))
            dictionary[outputs[i]].append(float(values[i]))
        
    return dictionary, outputs, float(timestep)
        

def plot_the_shit(t,a,b,Title):

    temp = 'Temp'
    press = 'Press'
    dens = 'systemDe'
    Ekin = 'KinEng'
    Epot = 'PotEng'
    
    Time = []
    convertionfactor = float(t/1000.0)   # 1000fsec in one ps. t=timestep
    ltime = 0
    for time in a['Step']:
        #print time
        value = float(time)*convertionfactor
        Time.append(value)
        ltime += 1


    ######### Does density values exist ? ##########
    try:
        Dens = a[dens]
    except:
        dens = None
    
    ############################################################ Figure 1
    
    Temp = a[temp]
    Press = a[press]
    Ek = a[Ekin]
    Ep = a[Epot]
    start = int(round(time/100))
    mid = int(round((ltime-start)/2.0)) + start
    if dens is not None:
        
        mean_dens = np.mean(Dens[start:])
        std_dens = np.std(Dens[start:])
        mean_temp = np.mean(Temp[start:])
        std_temp = np.std(Temp[start:])
        mean_press = np.mean(Press[start:])
        std_press = np.std(Press[start:])
        
        ml_temp = np.linspace(mean_temp,mean_temp,3)
        ml_press = np.linspace(mean_press,mean_press,3)
        ml_dens = np.linspace(mean_dens,mean_dens,3)
        
        # Three subplots, the axes array is 1-d
        f, axarr = plt.subplots(3, sharex=True)
        axarr[0].plot(Time[start:], Temp[start:], 'r-', label='Temp')
        axarr[0].errorbar([Time[start],Time[mid],Time[-1]],ml_temp,yerr=std_temp,errorevery=1,elinewidth=3,marker='s')
        #plt.errorbar()
        axarr[1].plot(Time[start:], Press[start:], 'y-',label='Press')
        axarr[1].errorbar([Time[start],Time[mid],Time[-1]],ml_press,yerr=std_press,errorevery=1,elinewidth=3,marker='s')
        axarr[2].plot(Time[start:], Dens[start:], 'g-',label='Dens')
        axarr[2].errorbar([Time[start],Time[mid],Time[-1]],ml_dens,yerr=std_dens,errorevery=1,elinewidth=3,marker='s')        

        axarr[0].set_title(Title)
        axarr[0].set_ylabel('temperature [K]')
        axarr[1].set_ylabel('pressure [atm]')
        axarr[2].set_ylabel('density [g/cm^3]')
        axarr[2].set_xlabel('time [ps]')
    
        axarr[0].legend(loc='best', fancybox=True, framealpha=0.5)    
        axarr[1].legend(loc='best', fancybox=True, framealpha=0.5)
        axarr[2].legend(loc='best', fancybox=True, framealpha=0.5)
    
    else:
        # Two subplots, the axes array is 1-d
        f, axarr = plt.subplots(2, sharex=True)
        axarr[0].plot(Time, Temp,label='Temp')
        axarr[1].plot(Time, Press,label='Press')
        
        axarr[0].set_title(Title)
        axarr[0].set_ylabel('temperature [K]')
        axarr[1].set_ylabel('pressure [atm]')
        axarr[1].set_xlabel('time [ps]')
    
        axarr[0].legend(loc='best', fancybox=True, framealpha=0.5)    
        axarr[1].legend(loc='best', fancybox=True, framealpha=0.5)
    
    
    ############################################################# Figure 2
    
    mean_ek = np.mean(Ek[start:])
    std_ek = np.std(Ek[start:])
    mean_ep = np.mean(Ep[start:])
    std_ep = np.std(Ep[start:])
    
    ml_ek = np.linspace(mean_ek,mean_ek,3)
    ml_ep = np.linspace(mean_ep,mean_ep,3)
    
    f2, ax = plt.subplots(2,sharex=True)
    ax[0].plot(Time[start:],Ek[start:],'r-',label='$E_{kin}$')
    ax[0].errorbar([Time[start],Time[mid],Time[-1]],ml_ek,yerr=std_ek,errorevery=1,elinewidth=3,marker='s')
    
    ax[1].plot(Time[start:],Ep[start:],'g-',label='$E_{pot}$')
    ax[1].errorbar([Time[start],Time[mid],Time[-1]],ml_ep,yerr=std_ep,errorevery=1,elinewidth=3,marker='s')

    ax[0].set_title(Title)
    ax[0].set_ylabel('energy [kcal/mol]')
    ax[1].set_ylabel('energy [kcal/mol]')
    ax[1].set_xlabel('time [ps]')
    ax[0].legend(loc='best', fancybox=True, framealpha=0.5)    
    ax[1].legend(loc='best', fancybox=True, framealpha=0.5)

    ##############################################################
    # Print output info

    print "#-----------------------------------------------------#"
    print "# T   = %.2f +- %.3f K " % (mean_temp, std_temp)
    print "# p   = %.2f +- %.3f atm " % (mean_press, std_press)
    print "# rho = %.3f +- %.3f g/cm^3 " % (mean_dens, std_dens)
    print "# Ek  = %.2f +- %.3f kcal/mol " % (mean_ek, std_ek)
    print "# Ep  = %.2f +- %.3f kcal/mol " % (mean_ep, std_ep)
    

    plt.show()
        
    

def main():
    # Read log.lammps and plot temperature, pressure and density in one plot

    ##########################################################################
    # Choose dataset:

    #Title = 'TraPPE CO2'
    #Title = 'EMP2 CO2'
    #Title = 'FPF CO2'
    Title = 'SPC/E water'

    logfilename = "log.lammps"
    timestep = 2.0   # [fsec] timestep used in simulation

    ##########################################################################
    if Title == 'TraPPE CO2':
        path = "/home/goran/lammps-28Jun14/theproject/carbondioxide/TraPPE/run_2015_02_27"         # TraPPE CO2

    if Title == 'EMP2 CO2':
        path = "/home/goran/lammps-28Jun14/theproject/carbondioxide/emp2/run_2015_02_27"           # EMP2 CO2
    
    if Title == 'FPF CO2':
        path = "/home/goran/lammps-28Jun14/theproject/carbondioxide/FPF/run_2015_02_26"            # FPF CO2        
    
    if Title == 'SPC/E water':
        #path = "/home/goran/lammps-28Jun14/theproject/water/run1_300K_500ps_npt" # bulk water system. 500ps simulation at 300K
        path = "/home/goran/lammps-28Jun14/theproject/water/run2_300K_500ps_npt" # bulk water system. 500ps simulation at 300K
        #path = "/home/goran/lammps-28Jun14/examples/Abel_runs/water/pure_H2O/Small_bulk_system"    # Small system of bulk water
    

    a,b,t = read_logfile(path,logfilename)
    

    if t is not None:
        print "#######################################"
        print "# Found time step value dt = %g" % t
        print "#-------------------------------------#"
        print "# We will therefore be using this!    #"
    else:
        print "#######################################"
        print "# Using time step dt = %g" % timestep
        print "#-------------------------------------#"
    
    plot_the_shit(timestep,a,b,Title)
    
    
    # End of main -------------------------------------------------------------    

    
if __name__ == "__main__":
    import os
    import matplotlib.pyplot as plt
    import numpy as np
    
    '''
    font = {'family' : 'normal',
        'size'   : 10}

    plt.rc('font', **font)
    '''
    main()
    