# -*- coding: utf-8 -*-
"""
Created on Tue Feb 10 18:41:05 2015

@author: goran

Read and display graphically ScalarInfo-file

"""



def main():
    
#    path1 = "/home/goran/lammps-28Jun14/examples/Abel_runs/alpha-quartz/nonbond"
#    path2 = "/home/goran/lammps-28Jun14/examples/Abel_runs/alpha-quartz/morse/run-2015-02-12"
#    path3 = "/home/goran/lammps-28Jun14/examples/Abel_runs/alpha-quartz/reax/bulksystem"
#    filename1 = "Adump.ScalarInfo_alpha-quartz"    
#    filename2 = "Adump.ScalarInfo_alpha-quartz"
#    filename3 = "Adump.ScalarInfo_alpha-quartz"
    
    path1 = "/home/goran/lammps-28Jun14/theproject/carbondioxide/FPF/run_2015_02_26"
    path2 = "/home/goran/lammps-28Jun14/theproject/carbondioxide/emp2/run_2015_02_27"
    path3 = "/home/goran/lammps-28Jun14/theproject/carbondioxide/TraPPE/run_2015_02_27"        
    filename1 = "Adump.CO2_FPF_scalarInfo.txt"    
    filename2 = "Adump.CO2_EMPtwo_scalarInfo.txt"
    filename3 = "Adump.CO2_TraPPE_scalarInfo.txt"
    
    nve_timesteps = 10000      # energy equilibtation
    npt_timesteps = 240000     # pressure equilibration
    nvt_timesteps = 0    # nvt run. main run
    
    dt = 2.0
    fsec_to_psec = dt/1000   # conversionfactor from timestep to ps
    
    nve = fsec_to_psec*nve_timesteps
    npt = fsec_to_psec*npt_timesteps + nve
    nvt = fsec_to_psec*nvt_timesteps + npt
    
    nve_line = np.linspace(nve,nve,5)
    npt_line = np.linspace(npt,npt,5)
    nvt_line = np.linspace(nvt,nvt,5)
    
    fullpath1 = os.path.join(path1,filename1)
    fullpath2 = os.path.join(path2,filename2)
    fullpath3 = os.path.join(path3,filename3)
    
    ofile1 = open(fullpath1,'r')
    ofile2 = open(fullpath2,'r')
    ofile3 = open(fullpath3,'r')
    
    # jump two fires lines:
    ofile1.readline(); ofile2.readline(); ofile3.readline()   
    ofile1.readline(); ofile2.readline(); ofile3.readline()    
    
    ###########################################################################
    #                  Density calculation from system volume
    factor = 1/0.6022
    Nmolec_reax = 144
    Nmolecules = 144                 # Check out the number of SiO2 molecules in the systems! Are they equal? 
    m_silicon = 28.0855
    m_oxygen = 15.9994
    mass = Nmolecules*(m_silicon + 2*m_oxygen)
    massreax = Nmolec_reax*(m_silicon + 2*m_oxygen)
    ###########################################################################    
    
    time1 = []; time2 = []; time3 = [] 
    temp1 = []; temp2 = []; temp3 = []
    press1 = []; press2 = []; press3 = []
    vol1 = []; vol2 = []; vol3 = []
    pot1 = []; pot2 = []; pot3 = []
    dens1 = []; dens2 = []; dens3 = []
    convert_pot_eng = True
    
    if ('alpha_quartz' in filename1):
        print "Running for alpha-quartz"
        for line in ofile1:
            t,K,P,V,Ek,Ep = line.split()
            time1.append(int(t)); temp1.append(float(K)); press1.append(float(P)); vol1.append(float(V)); pot1.append(float(Ep))
            dens1.append(factor*(mass/float(V)))
        for line in ofile2:
            t,K,P,V,Ek,Ep = line.split()
            time2.append(int(t)); temp2.append(float(K)); press2.append(float(P)); vol2.append(float(V)); pot2.append(float(Ep))
            dens2.append(factor*(mass/float(V)))
        for line in ofile3:
            t,K,P,V,Ek,Ep = line.split()
            # The reax calc gives pressure in bar. 1bar = 0.987 atm
            time3.append(int(t)); temp3.append(float(K)); press3.append(float(P)*0.987); vol3.append(float(V)); pot3.append(float(Ep))
            dens3.append(factor*(massreax/float(V)))
        
    if ('CO2' in filename1):
        print "Running for CO2"
        convert_pot_eng = False
        for line in ofile1:
            t,K,P,D,V,Ek,Ep = line.split()
            time1.append(int(t)); temp1.append(float(K)); press1.append(float(P)); dens1.append(float(D)); vol1.append(float(V)); pot1.append(float(Ep))
        for line in ofile2:
            t,K,P,D,V,Ek,Ep = line.split()
            time2.append(int(t)); temp2.append(float(K)); press2.append(float(P)); dens2.append(float(D)); vol2.append(float(V)); pot2.append(float(Ep))
        for line in ofile3:
            t,K,P,D,V,Ek,Ep = line.split()
            time3.append(int(t)); temp3.append(float(K)); press3.append(float(P)); dens3.append(float(D)); vol3.append(float(V)); pot3.append(float(Ep))
        
    Time1 = []; Time2 = []; Time3 = []
    index50 = 0; index100 = 0; index250 = 0;ii = 0
    for t in time1:
        thetime = t*fsec_to_psec
        Time1.append(thetime)
        if ((thetime < 51) and (thetime > 49)):
            index50 = ii
        if ((thetime < 101) and (thetime > 99)):
            index100 = ii
        if (thetime > 249):
            index250 = ii
        ii += 1
    
    Mindex50 = 0; Mindex100 = 0; Mindex250 = 0; ii = 0
    for t in time2:
        Time2.append(t*fsec_to_psec)
        if ((thetime < 51) and (thetime > 49)):
            Mindex50 = ii
        if ((thetime < 101) and (thetime > 99)):
            Mindex100 = ii
        if (thetime > 249):
            Mindex250 = ii
        ii += 1  
        
    Rindex50 = 0; Rindex100 = 0; Rindex250 = 0; ii = 0
    for t in time3:
        Time3.append(t*fsec_to_psec)
        if ((thetime < 51) and (thetime > 49)):
            Rindex50 = ii
        if ((thetime < 101) and (thetime > 99)):
            Rindex100 = ii
        if (thetime > 249):
            Rindex250 = ii
        ii += 1
        
    if (convert_pot_eng == True):
        eVtoJoule = 1.60217733*10**(-19)
        Na = 6.022*10**(23)
        Jouletokcal = 1/4184.0    
        #somefactor = (eVtoJoule*Na)*Jouletokcal
        somefactor = 1.60217733*6.022*Jouletokcal*10000
        print somefactor
        for i in range(len(pot3)):
            # convert from eV to kcal/mol
            pot3[i] = somefactor*pot3[i]/(Nmolec_reax*3)
    
    ################################################
    #              Print shit
    if ('alpha-quartz' in filename1):
        print "#####################################################################################"
        print "Potential energies [kcal/mol]:"
        print "         time:      t=0      t=50ps      t=100ps      t=250ps      density [g/cm^3]"
        print "ClayFF              %g    %g    %g    %g       %g" % (pot1[0], pot1[index50],pot1[index100],pot1[index250],dens1[-1])
        print "Morse               %g    %g    %g    %g       %g" % (pot2[0], pot2[Mindex50],pot2[Mindex100],pot2[Mindex250],dens2[-1])
        print "ReaxFF              %g    %g    %g    %g       %g" % (pot3[0], pot3[Rindex50],pot3[Rindex100],pot3[Rindex250],dens3[-1])
        print "#####################################################################################"
    else:
        lend = len(dens1)
        mean_dens1 = np.mean(dens1[int(round(lend/2.0)):-1]); std_dens1 = np.std(dens1[int(round(lend/2.0)):-1])
        mean_dens2 = np.mean(dens2[int(round(lend/2.0)):-1]); std_dens2 = np.std(dens2[int(round(lend/2.0)):-1])
        mean_dens3 = np.mean(dens3[int(round(lend/2.0)):-1]); std_dens3 = np.std(dens3[int(round(lend/2.0)):-1])
        print "#####################################################################################"
        print "Potential energies [kcal/mol]:"
        print "       time:   t=0      t=50ps      t=100ps      t=250ps      density [g/cm^3]"
        print "FPF            %.3f    %.3f    %.3f    %.3f       %.6f pm %.7f" % (pot1[0], pot1[index50],pot1[index100],pot1[index250],mean_dens1,std_dens1)
        print "EMP2           %.3f    %.3f    %.3f    %.3f       %.6f pm %.7f" % (pot2[0], pot2[Mindex50],pot2[Mindex100],pot2[Mindex250],mean_dens2,std_dens2)
        print "TraPPE         %.3f    %.3f    %.3f    %.3f       %.6f pm %.7f" % (pot3[0], pot3[Rindex50],pot3[Rindex100],pot3[Rindex250],mean_dens3,std_dens3)
        print "#####################################################################################"
        
    ################################################
    #              Figure 1  Potential energy
    if ('CO2' in filename1):
        legends = ['FPF','EMP2','TraPPE']
    else:
        legends = ['ClayFF','Morse','ReaxFF']
    rc = ['b-','r-','y-']
    location = 'upper right'
    Title = 'Potential energy'
    Xlabel = 'time [ps]'
    Ylabel = 'Ep [Kcal/mol]'
    plt.figure()
    min_pot1 = min(pot1); max_pot1 = max(pot1)
    min_pot2 = min(pot2); max_pot2 = max(pot2)
    minmin = min([min_pot1,min_pot2])
    maxmax = max([max_pot1,max_pot2])
    yline = np.linspace(minmin,maxmax,5)
    plt.hold(True)
    plt.plot(Time1,pot1,rc[0])
    plt.plot(Time2,pot2,rc[1])
    plt.plot(Time3,pot3,rc[2])
    if (nve_timesteps > 0):
        plt.annotate('nve', xy=(1.0, 1.0), xytext=(nve/2, (maxmax+minmin)/2))
        plt.plot(nve_line,yline,'k--')
        legends.append('nve')
    if (npt_timesteps > 0):
        plt.annotate('npt', xy=(1.0, 1.0), xytext=((npt+nve)/2, (maxmax+minmin)/2))
        plt.plot(npt_line,yline,'k--')
        legends.append('npt')
    if (nve_timesteps > 0):
        plt.annotate('nvt', xy=(1.0, 1.0), xytext=((nvt+npt)/2, (maxmax+minmin)/2))
        plt.plot(nvt_line,yline,'k--')
        legends.append('nvt')
    plt.hold(False)
    plt.title(Title)
    plt.xlabel(Xlabel); plt.ylabel(Ylabel); plt.legend(legends,loc=location)

    ################################################
    #              Figure 2
    Title = 'Density'
    Ylabel = r'$\rho $ [g/cm^3]'
    if ('CO2' in filename1):
        legends = ['FPF','EMP2','TraPPE','expected']
        rho = 0.0017966 # density of CO2 300K 1atm. Isothermal Properties of CO2; NIST webbook.nist.gov
    else:
        legends = ['ClayFF','Morse','ReaxFF','expected']
        rho = 2.65 # density of alpha_quartz 300K 1atm
        
    expected = np.linspace(rho,rho,len(Time1))
    plt.figure()
    plt.hold(True)
    plt.plot(Time1,dens1,rc[0])
    plt.plot(Time2,dens2,rc[1])
    plt.plot(Time3,dens3,rc[2])
    plt.plot(Time1,expected,'k--')
    plt.hold(False)
    plt.title(Title)
    plt.xlabel(Xlabel); plt.ylabel(Ylabel); plt.legend(legends,loc=location)

    ################################################
    #              Figure 3
    Title = 'Temperature'
    Ylabel = 'temperature [K]'
    plt.figure()
    plt.hold(True)
    plt.plot(Time1,temp1,rc[0])
    plt.plot(Time2,temp2,rc[1])
    plt.plot(Time3,temp3,rc[2])
    plt.hold(False)
    plt.title(Title)
    plt.xlabel(Xlabel); plt.ylabel(Ylabel); plt.legend(legends,loc=location)
    
    ################################################
    #              Figure 4
    Title = 'Pressure'
    Ylabel = 'pressure [atm]'
    plt.figure()
    plt.hold(True)
    plt.plot(Time1,press1,rc[0])
    plt.plot(Time2,press2,rc[1])
    plt.plot(Time3,press3,rc[2])
    plt.hold(False)
    plt.title(Title)
    plt.xlabel(Xlabel); plt.ylabel(Ylabel); plt.legend(legends,loc=location)
        
    
        
if (__name__ == "__main__"):
    import matplotlib.pyplot as plt
    import numpy as np
    import os
    main()