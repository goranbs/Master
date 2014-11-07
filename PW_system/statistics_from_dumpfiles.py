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


def plotEnergy(energyfiles,tc):
    '''
    plots the energy of the system as a funciton of time.\n
    Takes LAMMPS outputfiles that contains columns of timestep and energy\n
    tc = time conversion factor
    '''

def plotScalarInfo(path,filename,tc,timeunits,Format,saveplot=False,showplot=True):
    '''
    path is the full path to filename. Where the filename is a scalar info -file.
    ScalarInfoFile contains information like:\n
    Ts = Temperature of system\n
    Tp = Temperature of portlandite\n
    Tw = Temperature of water\n
    p_OH = Intermolecular distance (bond length) of OH in portlandite\n
    w_OH = Intermolecular distance (bond length) of OH in water\n
    p_ROG = Radius of gyration\n
    w_ROG = Radius of gyration\n
    P = global pressure\n
    vx,vy,vz = average velocity in x,y,z-dir
    D = displacement
    '''

    ScalarInfoFile = os.path.join(path,filename)
        
    t = []
    Ts = []; Tp = []; Tw = []; P = []; D = []; p_OH = []; w_OH = []; p_ROG = []; w_ROG = []
    vx = []; vy = []; vz = []; Ekin = []; Epot = []
    
    ofile = open(ScalarInfoFile,'r')
    ofile.readline() # the first line is a header, and should be skipped
    scalarvalues = ofile.readline()
    scalarvalues = scalarvalues.split() # plint on whitespaces
    scalarvalues.pop(0)  # pop the first value (remove a #)
    
    for line in ofile:
        ll = line.split()
        t.append(float(ll[0])*tc)  # convert directly to the disirered timeunits
        index = 0
        for i in scalarvalues:
            if (i == "c_myTemp"):
                Ts.append(float(ll[index]))
            if (i == "c_tempW"):
                Tw.append(float(ll[index]))
            if (i == "c_tempP"):
                Tp.append(float(ll[index]))
            if (i in ("v_blW" , "c_blW")):
                w_OH.append(float(ll[index]))
            if (i in ("v_blP" , "c_blP")):
                p_OH.append(float(ll[index]))
            if (i == "c_Press"):
                P.append(float(ll[index]))
            if (i in ("c_vwx" , "c_vx")):
                vx.append(float(ll[index]))
            if (i in ("c_vwy" , "c_vy")):
                vy.append(float(ll[index]))
            if (i in ("c_vwz" , "c_vz")):
                vz.append(float(ll[index]))
            if (i in ("c_disp" , "c_Disp")):
                D.append(float(ll[index]))
            if ("Ekin" in i):
                Ekin.append(float(ll[index]))
            if ("Epot" in i):
                Epot.append(float(ll[index]))
                
            index += 1
            
    #########################################################################
    #             Averaging and finding standard deviations
    mean_wOH = np.mean(w_OH)
    std_wOH = np.std(w_OH)
    mean_pOH = np.mean(p_OH)
    std_pOH = np.std(p_OH)
    mean_wROG = np.mean(w_ROG)
    std_wROG = np.std(w_ROG)
    mean_pROG = np.mean(p_ROG)
    std_pROG = np.std(p_ROG)
    
    lt = len(t)
    lenvx = len(vx)
    do_not_plot_vel = False
    if ((lenvx > 0) and (len(vx) == len(vy)) and (len(vy) == len(vz))):        
        vw = []
        for i in range(lenvx):
            vw.append(np.sqrt(vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i]))
        
        mean_v = np.mean(vw[int(lt/10.0):-1])   # v
        std_v = np.std(vw[int(lt/10.0):-1])
        mean_vx = np.mean(vx[int(lt/10.0):-1])  # vx
        std_vx = np.std(vx[int(lt/10.0):-1])
        mean_vy = np.mean(vy[int(lt/10.0):-1])  # vy
        std_vy = np.std(vy[int(lt/10.0):-1])
        mean_vz = np.mean(vz[int(lt/10.0):-1])  # vz
        std_vz = np.std(vz[int(lt/10.0):-1])
    else:
        print "###############################################"
        print " Error!\n len(vx)=%g, len(vy) =%g, len(vz)=%g " % (len(vx),len(vy),len(vz))
        print "###############################################"
        do_not_plot_vel = True
   
    mean_P = np.mean(P[int(lt/10.0):-1])
    std_P = np.std(P[int(lt/10.0):-1])
    mean_Ts = np.mean(Ts[int(lt/10.0):-1])
    std_Ts = np.std(Ts[int(lt/10.0):-1])
    mean_Tp = np.mean(Tp[int(lt/10.0):-1])
    std_Tp = np.std(Tp[int(lt/10.0):-1])
    mean_Tw = np.mean(Tw[int(lt/10.0):-1])
    std_Tw = np.std(Tw[int(lt/10.0):-1])
        
    #########################################################################
    #                             Plotting

    from matplotlib import rc
    rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
    ## for Palatino and other serif fonts use:
    #rc('font',**{'family':'serif','serif':['Palatino']})
    rc('text', usetex=True)    
    
    if ((len(Tp) == lt) and (len(Tp) == len(Tw))):
        plt.figure()
        plt.plot(t,Ts,'b-')
        plt.hold(True)
        plt.plot(t,Tp,'y--')
        plt.plot(t,Tw,'r--')
        plt.hold(False)
        plt.title(r'Temperature in the system\newline $T_{system} = $ %g $ \pm $ %g K \newline $T_{portlandite} = $ %g $ \pm $ %.2f K \newline  $T_{water} = $ %g $ \pm $ %g K' % (mean_Ts, std_Ts,mean_Tp, std_Tp,mean_Tw, std_Tw) )
        plt.legend([r'$ T_{system}(t)$ ',r'T_{portlandite}(t)$',r'$T_{water}(t)$'],loc='lower right')
        plt.xlabel('Time [%s]' % timeunits)
        plt.ylabel('Temperature [K]')
        if (saveplot):
            fig_name = 'Temp_%s.png' % filename
            print "saving %s " % fig_name
            plt.savefig(fig_name,format=Format)
            
    else:
        plt.figure()
        plt.plot(t,Ts,'b-')
        plt.title('Temperature in the system. $ T_{system} = $ %g $ \pm $ %g K' % (mean_Ts, std_Ts))
        plt.xlabel('Time [%s]'  % timeunits), plt.ylabel('Temperature [K]'), plt.legend(['T(t)'])
        if (saveplot):
            fig_name = 'Temperature_%s.png' % filename
            print "saving %s " % fig_name
            plt.savefig(fig_name,format=Format)
            
    start = int(lt/10.0)

    
    #print "len(P) = %g " % len(P)
    if (lt == len(P)):
        plt.figure()
        plt.plot(t,P,'b-')
        plt.title(r'Pressure in the system\newline $P = $ %g $ \pm $ %g atm' % (mean_P,std_P))
        plt.xlabel('Time [%s]' % timeunits), plt.ylabel('Pressure [atm]'), plt.legend(['P(t)'],loc='lower right')
        if (saveplot):
            fig_name = 'Press_%s.png' % filename
            print "saving %s " % fig_name
            plt.savefig(fig_name,format=Format)
            
    #print "len(w_OH) = %g " % len(w_OH)        
    if (lt == len(w_OH)):
        plt.figure()
        plt.plot(t,w_OH,'b-')
        plt.title(r'Bond lenght for OH in water\newline $d_{wOH} = $ %g $ \pm $ %g Aa' % (mean_wOH,std_wOH))
        plt.xlabel('Time [%s]' % timeunits), plt.ylabel(r'$d_{wOH}$ [Aa]'), plt.legend([r'$d_{wOH}(t)$'],loc='lower right')
        if (saveplot):
            fig_name = 'Press_%s.png' % filename
            print "saving %s " % fig_name
            plt.savefig(fig_name,format=Format)

    #print "len(p_OH) = %g " % len(p_OH)        
    if (lt == len(p_OH)):
        plt.figure()
        plt.plot(t,p_OH,'b-')
        plt.title(r'Bond lenght for OH in portlandite\newline $d_{pOH} = $ %g $ \pm $ %g Aa' % (mean_pOH,std_pOH))
        plt.xlabel('Time [%s]' % timeunits), plt.ylabel(r'$d_{wOH}$ [Aa]'), plt.legend([r'$d_{pOH}(t)$'],loc='lower right')
        if (saveplot):
            fig_name = 'Press_%s.png' % filename
            print "saving %s " % fig_name
            plt.savefig(fig_name,format=Format)
            
    if (do_not_plot_vel == False):
        plt.figure()
        plt.plot(t[start:-1],vw[start:-1],'b-')
        plt.hold(True)
        plt.plot(t[start:-1],vx[start:-1],'g--')
        plt.plot(t[start:-1],vy[start:-1],'y--')
        plt.plot(t[start:-1],vz[start:-1],'m--')
        plt.hold(False)
        plt.title(r'Velocityprofile \newline $v_{rms} = $ %g $ \pm $ %g Aa/%s \newline $v_x = $%g $ \pm $ %g Aa/%s \newline $v_y = $%g $ \pm $ %g Aa/%s \newline $v_z = $%g $ \pm $ %g Aa/%s' % (mean_v,std_v,timeunits,mean_vx,std_vx,timeunits, mean_vy,std_vy,timeunits, mean_vz,std_vz, timeunits))
        plt.xlabel('Time [%s]' % timeunits)
        plt.ylabel('Velocity [Aangstrom/%s]' % timeunits)
        plt.legend(['v','vx','vy','vz'],loc='lower right')
        if (saveplot):
            fig_name = 'Velocityprofile_%s.png' % filename
            print "saving %s " % fig_name
            plt.savefig(fig_name,format=Format)
    
    if (lt == len(D)):
        plt.figure()
        plt.plot(t,D,'b-')
        plt.title('Mean displacement of the water in the system')
        plt.legend(['D(t)'],loc='upper left')
        plt.xlabel('t [%s]' % timeunits)
        plt.ylabel('Displacement [Aangsrom]')
        if (saveplot):
            fig_name = 'MSD_%s.png' % filename
            print "saving %s " % fig_name
            plt.savefig(fig_name,format=Format)
    
        Coeff_D = np.array(D)/(6*np.array(t))
        mean_D = np.mean(Coeff_D[int(lt/2.0):-1])
        std_D = np.std(Coeff_D[int(lt/2.0):-1])

        tt = np.array(t[int(lt/2.0):-1])
        Diffusion = np.zeros([len(tt),1])
        Diffusion.fill(mean_D) 
        yerror = np.zeros([len(tt),1])
        yerror.fill(std_D)
        
        plt.figure()
        plt.plot(t[int(lt/10.0):-1],Coeff_D[int(lt/10.0):-1],'b-')
        plt.hold(True)
        plt.plot(tt,Diffusion,'r--')
        plt.hold(False)
        plt.title(r'Self diffusion of water in the system.\newline Diffusion constant $D=$ %.2f $ \pm $ %.3f $ Aa^2/%s $' % (mean_D,std_D,timeunits))
        plt.xlabel('Time [%s]' % timeunits)
        plt.ylabel(r'D [${Aangsrom}^2$/%s]' % timeunits)
        if (saveplot):
            fig_name = 'Diffusion_%s.png' % filename
            print "saving %s " % fig_name
            plt.savefig(fig_name,format=Format)
            
    lEkin = len(Ekin)
    if ((lEkin == lt) and (lEkin == len(Epot))):
        Etot = []
        for i in range(lEkin):
            Etot.append((Ekin[i] + Epot[i]))
        
        plt.figure()
        plt.plot(t,Ekin,'r--')
        plt.hold(True)
        plt.plot(t,Epot,'b--')
        plt.plot(t,Etot,'y-')
        plt.title('Energy in the system')
        plt.xlabel('Time [%s]' % timeunits), plt.ylabel('Energy [Kcal/mole]')
        plt.legend(['Ekin','Epot','Etot'],loc='lower right')
        plt.hold(False)
        if (saveplot):
            fig_name = 'Energy_%s.png' % filename
            print "saving %s " % fig_name
            plt.savefig(fig_name,format=Format)
        
    else:
        print "####################################"
        print "len(Ekin) = %g \nlen(Epot) = %g " % (len(Ekin),len(Epot))
    
        
    plt.show(showplot)
    plt.close('all')

def plotRDF(path,filename,tc,timeunits,Format,saveplot,showplot):
    '''
    plot the evolution of the radial distribution function.\n
    plotDRF takes the path to the filename, and the filename as the two firs
    input args. tc is the time conversion factor that converts the timestep
    unit of the input file into the disiered time units.\n
    Choose the format of the output plots as Format = 'png','jpg','jpeg'...\n
    '''
    RDF_file = os.path.join(path,filename)

    ofile = open(RDF_file,'r') # open the filename in read mode
    ofile.readline()           # read unimportant first line of the file
    ofile.readline()           # and second
    header = ofile.readline()  # contains information about the following lines.
    header = header.split()    # split line on whitespaces
    header.pop(0)              # remove the first element (it is a #)

    t = []
    Row = []; BinCoord = []; RDF = []; Coord = []    

    MoreToRead = True
    time = 0
    while (MoreToRead):
        try:
            tr_line = ofile.readline()
            tr_line = (tr_line.split()) # tr_line[0] = TimeStep, tr_line[1] = Number-of-rows 
            t.append(float(tr_line[0])*tc)
            MoreToRead = True
            Row.append([]);BinCoord.append([]);RDF.append([]);Coord.append([])

        except:
            MoreToRead = False
            break
        
        for i in range(int(tr_line[1])): # go through the next tr_line[1] number of lines in ofile
            line = ofile.readline()
            line = line.split()
            Row[time].append( int(line[0]))
            BinCoord[time].append(float(line[1]))
            RDF[time].append(float(line[2]))
            Coord[time].append(float(line[3]))
        time += 1  # number of timesteps that are included (RDF(t))
    ###########################################################################
    # -------------- Plotting the radial distribution function ----------------
    
    plt.figure()
    legends = []
    plt.hold(True)        
    for j in range(len(t)):
        plt.plot(BinCoord[j],RDF[j])
        legend = 't=%.1f %s' % (t[j],timeunits)
        legends.append(legend)
    plt.hold(False)
    plt.title('Radial Distribution Function')
    plt.xlabel('distance [Aa]'), plt.ylabel('average number of particles'), plt.legend(legends,loc='upper left')
    if (saveplot):
        fig_name = 'RDF_%s.png' % filename
        print "saving %s " % fig_name
        plt.savefig(fig_name,format=Format)
    
    plt.figure()
    plt.plot(BinCoord[-1],RDF[-1])
    plt.title('Radial Distribution Function')
    plt.xlabel('distance [Aa]'), plt.ylabel('average number of particles'), plt.legend(['RDF(t=%g)' % t[-1]],loc='upper left')
    if (saveplot):
        fig_name = 'RDF_%g%s_%s.png' % (t[-1],timeunits,filename)
        print "saving %s " % fig_name 
        plt.savefig(fig_name,format=Format)
        
    plt.show(showplot)
    plt.close('all')
    
    
def main():
    #generalpath = "/home/goran/lammps-28Jun14/examples/water_portlandite_system"    

    generalpath = "/home/goran/lammps-28Jun14/examples/"
    wp_path = "Abel_runs/carbondioxide"
    #wp_path = "npt_run_and_energyminimization/dump/"
    arg = "Adump"
    path = os.path.join(generalpath,wp_path)
    filenames = gothroughfiles(path,arg)
    
    '''
    There are different types of files in filename that we are going to do different operation on
    so we should find these types of files, and send them to the correct functions.    
    '''

    #######################################################
    # Sort the filenames in the directory specified by path
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
    

    timeunits='ps'
    ftop = 1000           # [fsec] in [ps]
    tfac = 100            # timesteps in one fsec
    psfactor = ftop*tfac  # convert to [ps] 
    tc = 1.0/psfactor     # time conversionfactor
    
    
    saveplot = True
    showplot = False
    Format = 'png'

    for filename in ScalarInfo_files:
        print filename
        plotScalarInfo(path,filename,tc,timeunits,Format,saveplot,showplot)
    
    
    for filename in RDF_files:
        print filename
        plotRDF(path,filename,tc,timeunits,Format,saveplot,showplot)

if (__name__ == "__main__"):
    import numpy as np
    import matplotlib.pyplot as plt
    import os
    main()


