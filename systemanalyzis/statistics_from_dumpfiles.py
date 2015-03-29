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


def plotScalarInfo(writetothisfile,path,filename,tc,timeunits,Format,saveplot=False,showplot=True):
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

    if ("CO2" in filename):
        system_type = 'CO2'
    else:
        system_type = 'water'
    
    ScalarInfoFile = os.path.join(path,filename)
        
    t = []
    Ts = []; Tp = []; Tw = []; P = []; D = []; p_OH = []; w_OH = []
    vx = []; vy = []; vz = []; Ekin = []; Epot = []; msd = []; Tco2 = []; bondl_CO2 = []
    
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
            if ("myTemp" in i):
                Ts.append(float(ll[index]))
            if ("temperature" in i):
                Ts.append(float(ll[index]))
            if ("tempW" in i):
                Tw.append(float(ll[index]))
            if ("tempP" in i):
                Tp.append(float(ll[index]))
            if ("tempCO2" in i):
                Tco2.append(float(ll[index]))
            if (i in ("v_blW" , "c_blW")):
                w_OH.append(float(ll[index]))
            if (i in ("v_blP" , "c_blP")):
                p_OH.append(float(ll[index]))
            if ("ress" in i):
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
            if ("MSD" in i):
                msd.append(float(ll[index]))
            if ("Bondl_CO2" in i):
                bondl_CO2.append(float(ll[index]))
                
            index += 1
            
    #########################################################################
    #             Averaging and finding standard deviations
    mean_wOH = np.mean(w_OH)
    std_wOH = np.std(w_OH)
    mean_pOH = np.mean(p_OH)
    std_pOH = np.std(p_OH)
    
    lt = len(t)    
    start = int(3*lt/4.0)
    lenvx = len(vx)
    do_not_plot_vel = False
    if ((lenvx > 0) and (len(vx) == len(vy)) and (len(vy) == len(vz))):        
        vw = []
        for i in range(lenvx):
            vw.append(np.sqrt(vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i]))
        
        mean_v = np.mean(vw[start:-1])   # v
        std_v = np.std(vw[start:-1])
        mean_vx = np.mean(vx[start:-1])  # vx
        std_vx = np.std(vx[start:-1])
        mean_vy = np.mean(vy[start:-1])  # vy
        std_vy = np.std(vy[start:-1])
        mean_vz = np.mean(vz[start:-1])  # vz
        std_vz = np.std(vz[start:-1])
    else:
        print "####################################################"
        print " No plots of velocity dist will be generated, cause:"
        print " len(vx)=%g, len(vy) =%g, len(vz)=%g " % (len(vx),len(vy),len(vz))
        print "####################################################"
        do_not_plot_vel = True
   
    mean_P = np.mean(P[start:])
    std_P = np.std(P[start:])
    mean_Ts = np.mean(Ts[start:])
    std_Ts = np.std(Ts[start:])
        
    #########################################################################
    #                             Plotting

    rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
    ## for Palatino and other serif fonts use:
    #rc('font',**{'family':'serif','serif':['Palatino']})
    rc('text', usetex=True)   
    
    lTp = len(Tp); lTw = len(Tw); lTs = len(Ts)
    print "####################################################"
    print " INFO!\n len(t) = %g. len(Ts)=%g, len(Tw)=%g, len(Tp)=%g " % (lt,lTs, lTw, lTp) 
    print "####################################################"

    if ((len(Tp) == lt) and (len(Tp) == len(Tw))):
        mean_Tp = np.mean(Tp[start:])
        std_Tp = np.std(Tp[start:])
        mean_Tw = np.mean(Tw[start:])
        std_Tw = np.std(Tw[start:])
        
        Title = ' Temperature in the system. $T_{sys} = %.1f \pm %.2f K$' % (mean_Ts, std_Ts)
        Xlabel = 'Time [%s]' % timeunits
        Ylabel = 'Temperature [K]'
        legends = ['$T_{system}(t)$ ','T_{portlandite}(t)$','$T_{water}(t)$']
        linestyles = ['b-','y--','r--']
        plt.figure()
        plt.plot(t,Ts,linestyles[0])
        plt.hold(True)
        plt.plot(t,Tp,linestyles[1])
        plt.plot(t,Tw,linestyles[2])
        plt.hold(False)
        plt.title(r'Temperature in the system\newline $T_{system} = %g \pm %g K $\newline $T_{portlandite} = %g \pm %.2f K $\newline  $T_{water} = %g \pm %g K $' % (mean_Ts, std_Ts,mean_Tp, std_Tp,mean_Tw, std_Tw))
        plt.legend([r'$ T_{system}(t)$ ',r'T_{portlandite}(t)$',r'$T_{water}(t)$'],loc='upper right')
        plt.xlabel(Xlabel)
        plt.ylabel(Ylabel)
        
        fig_name = 'Temp_%s' % filename
        write_to_file(writetothisfile,Title,Xlabel,Ylabel,fig_name,legends,linestyles,[t],[Ts,Tp,Tw])
        
        if (saveplot):
            fig_name = fig_name + '.png'
            print "saving %s " % fig_name
            plt.savefig(fig_name,format=Format)
            
    else:
        if (len(Tco2) == lt):
            mean_Tco2 = np.mean(Tco2[start:])
            std_Tco2 = np.std(Tco2[start:])

            Title = 'Temperature in the system. $T_{sys} = %.1f \pm %.2f K$' % (mean_Ts, std_Ts)
            Xlabel = 'Time [%s]' % timeunits
            Ylabel = 'Temperature [K]'
            legends = ['$T_{system}(t)$','$T_{CO_2}(t)$']
            linestyles = ['b-','r--']
            
            plt.figure()
            plt.plot(t,Ts,linestyles[0])
            plt.hold(True)
            plt.plot(t,Tco2,linestyles[1])
            plt.title(r'Temperature in the system.\newline $ T_{system} = %g \pm %g K $\newline $T_{CO_2} = %g \pm %g K $' % (mean_Ts, std_Ts, mean_Tco2, std_Tco2))
            plt.legend(legends,loc='upper right'), plt.xlabel(Xlabel), plt.ylabel(Ylabel)
            
            fig_name = 'Temp_CO2_%s' % filename
            write_to_file(writetothisfile,Title,Xlabel,Ylabel,fig_name,legends,linestyles,[t],[Ts,Tco2])
            
            if (saveplot):
                fig_name += '.png'
                print "saving %s " % fig_name
                plt.savefig(fig_name,format=Format)
        else:
            
            Title = 'Temperature in the system. $T_{system} = %.1f \pm %.2f K $' % (mean_Ts, std_Ts)
            Xlabel = 'Time [%s]'  % timeunits; Ylabel = 'Temperature [K]'; legends = ['T(t)']
            linestyles = ['b-d']
            
            plt.figure()
            plt.plot(t,Ts,linestyles[0])
            plt.title(Title)
            plt.xlabel(Xlabel), plt.ylabel(Ylabel), plt.legend(legends,loc='upper right')
            
            fig_name = 'Temp_%s' % filename
            write_to_file(writetothisfile,Title,Xlabel,Ylabel,fig_name,legends,linestyles,[t],[Ts])
            if (saveplot):
                fig_name += '.png'
                print "saving %s " % fig_name
                plt.savefig(fig_name,format=Format)
            
    
    #print "len(P) = %g " % len(P)
    if (lt == len(P)):
        
        Title = 'Pressure in the system. $P = %.1f \pm %.2f atm $' % (mean_P,std_P)
        Xlabel = 'Time [%s]' % timeunits
        Ylabel = 'Pressure [atm]'
        legends = ['P(t)','mean(P)']
        linestyles = ['b-','r--x']
        
        plt.figure()
        plt.plot(t,P,linestyles[0])
        plt.hold(True)
        tmean = t[start:]
        Pmean = np.linspace(mean_P,mean_P,len(t[start:]))
        plt.plot(tmean,Pmean,linestyles[1])
        plt.hold(False)
        plt.title(r'Pressure in the system\newline $P = %g \pm %g atm $' % (mean_P,std_P))
        plt.xlabel('Time [%s]' % timeunits), plt.ylabel('Pressure [atm]'), plt.legend(['P(t)'],loc='upper right')
  
        fig_name = 'Press_%s' % filename
        write_to_file(writetothisfile,Title,Xlabel,Ylabel,fig_name,legends,linestyles,[t,tmean],[P,Pmean])

        if (saveplot):
            fig_name += '.png'
            print "saving %s " % fig_name
            plt.savefig(fig_name,format=Format)
            
    #print "len(w_OH) = %g " % len(w_OH)        
    if (lt == len(w_OH)):
        
        Title = ' Bond lenght for OH in water\\newline $d_{wOH} = %g \pm %g \AA{} $' % (mean_wOH,std_wOH)
        Xlabel = 'Time [%s]' % timeunits
        Ylabel = '$d_{wOH}$ [\AA{}]'
        legends = ['$d_{wOH}(t)$']
        linestyles = ['b-d']
        
        plt.figure()
        plt.plot(t,w_OH,linestyles)
        plt.title(r'Bond lenght for OH in water\newline $d_{wOH} = %g \pm %g Aa $' % (mean_wOH,std_wOH))
        plt.xlabel('Time [%s]' % timeunits), plt.ylabel(r'$d_{wOH}$ [Aa]'), plt.legend([r'$d_{wOH}(t)$'],loc='lower right')
        fig_name = 'Bondlength_w_%s' % filename
        write_to_file(writetothisfile,Title,Xlabel,Ylabel,fig_name,legends,linestyles,[t],[w_OH])
        
        if (saveplot):
            fig_name = fig_name + '.png'
            print "saving %s " % fig_name
            plt.savefig(fig_name,format=Format)

    #print "len(p_OH) = %g " % len(p_OH)        
    if (lt == len(p_OH)):
        
        Title = ' Bond lenght for OH in portlandite\\newline $d_{pOH} = %g \pm %g \AA{} $' % (mean_pOH,std_pOH)
        Xlabel = 'Time [%s]' % timeunits
        Ylabel = '$d_{pOH}$ [\AA{}]'
        legends = ['$d_{pOH}(t)$']
        linestyles = ['b-d']

        plt.figure()
        plt.plot(t,p_OH,'b-d')
        plt.title(r'Bond lenght for OH in portlandite\newline $d_{pOH} = %g \pm %g Aa $' % (mean_pOH,std_pOH))
        plt.xlabel('Time [%s]' % timeunits), plt.ylabel(r'$d_{pOH}$ [Aa]'), plt.legend([r'$d_{pOH}(t)$'],loc='lower right')
        fig_name = 'Bondlength_p_%s' % filename
        write_to_file(writetothisfile,Title,Xlabel,Ylabel,fig_name,legends,linestyles,[t],[p_OH])

        if (saveplot):
            fig_name = fig_name + '.png'
            print "saving %s " % fig_name
            plt.savefig(fig_name,format=Format)
            
    if (do_not_plot_vel == False):
        Title = 'Velocityprofile\\newline $v_{rms} = %g \pm %g \AA{}/%s $' % (mean_v,std_v,timeunits)
        Xlabel = 'Time [%s]' % timeunits
        Ylabel = 'Velocity [$\AA{}/%s$]' % timeunits
        legends = ['v','vx','vy','vz']
        linestyles = ['b-','g--','y--','m--']
        tt = t[start:-1]
        plt.figure()
        plt.plot(tt,vw[start:-1],'b-')
        plt.hold(True)
        plt.plot(tt,vx[start:-1],'g--')
        plt.plot(tt,vy[start:-1],'y--')
        plt.plot(tt,vz[start:-1],'m--')
        plt.hold(False)
        plt.title(r'Velocityprofile \newline $v_{rms} = %g \pm %g Aa/%s $ \newline $v_x = %g \pm %g Aa/%s $\newline $v_y = %g \pm %g Aa/%s $\newline $v_z = %g \pm %g Aa/%s $' % (mean_v,std_v,timeunits,mean_vx,std_vx,timeunits, mean_vy,std_vy,timeunits, mean_vz,std_vz, timeunits))
        plt.xlabel('Time [%s]' % timeunits)
        plt.ylabel('Velocity [Aangstrom/%s]' % timeunits)
        plt.legend(['v','vx','vy','vz'],loc='lower right')
        fig_name = 'Velocityprofile_%s' % filename
        write_to_file(writetothisfile,Title,Xlabel,Ylabel,fig_name,legends,linestyles,[tt],[vw[start:-1],vx[start:-1],vy[start:-1],vz[start:-1]])
        if (saveplot):
            fig_name = 'Velocityprofile_%s.png' % filename
            print "saving %s " % fig_name
            plt.savefig(fig_name,format=Format)
    
    if (lt == len(msd)):
        start = int(3*lt/4.0)
        tt = t[start:]
        deg = 1
        p,vv = np.polyfit(tt,msd[start:],deg, full=False,cov=True)
        std_D = np.sqrt(vv[0,0])
        f = np.polyval(p,t)
        Title = 'Mean square displacement of %s in the system' % system_type
        Xlabel = 'Time [%s]' % timeunits
        Ylabel = 'Displacement [\AA{}]'
        legends = ['msd(t)','D=%.2f(%.3f)[$\AA{}^2/%s$]' % (p[0]/6.0,std_D,timeunits) ]
        linestyles = ['b-','r--']
        
        plt.figure()
        plt.plot(t,msd,linestyles[0])
        plt.hold(True)
        plt.plot(t,f,linestyles[1])
        plt.title(Title)
        #plt.legend(['msd(t)',r'b + ax = %.2f(%.3f) + %.2f(%.3f)x ' % (p[1],np.sqrt(v[1,1]),p[0],std_D) ],loc='upper left')
        plt.legend(['msd(t)',r'D = %.2f(%.3f) [$Aa^2/%s$]' % (p[0]/6.0,std_D,timeunits) ],loc='upper left')
        plt.xlabel(Xlabel)
        plt.ylabel('Displacement [Aa]')
        fig_name = 'MSD_%s' % filename
        write_to_file(writetothisfile,Title,Xlabel,Ylabel,fig_name,legends,linestyles,[t],[msd,f])
        
        if (saveplot):
            fig_name = fig_name + '.png'
            print "saving %s " % fig_name
            plt.savefig(fig_name,format=Format)
    
        Coeff_D = np.array(msd)/(6*np.array(t))
        mean_D = np.mean(Coeff_D[start:])
        std_D = np.std(Coeff_D[start:])
        
        tt = np.array(t[start:])
        Diffusion = np.linspace(mean_D,mean_D,len(tt))
        
        #yerror = np.linspace(std_D,std_D,len(tt))
        
        plt.figure()
        plt.plot(t,Coeff_D,'b-')
        plt.hold(True)
        plt.plot(tt,Diffusion,'r--')
        plt.hold(False)
        plt.title(r'Self diffusion of %s in the system.\newline Diffusion constant $D= %.2f \pm %.3f Aa^2/%s $' % (system_type,mean_D,std_D,timeunits))
        plt.xlabel(Xlabel)
        plt.ylabel(r'D [${Aa}^2/%s $]' % timeunits)
        Title = 'Self diffusion of %s in the system.\\newline Diffusion constant $D= %.2f \pm %.3f \AA{}^2/%s $' % (system_type,mean_D,std_D,timeunits)
        Ylabel = 'Displacement [\AA{}]'
        legends = ['$D(t)$','$D_{mean}(t)$']
        linestyles = ['b-','r--']
        write_to_file(writetothisfile,Title,Xlabel,Ylabel,fig_name,legends,linestyles,[t,tt],[Coeff_D,Diffusion])
        if (saveplot):
            fig_name = 'Diffusion_%s.png' % filename
            print "saving %s " % fig_name
            plt.savefig(fig_name,format=Format)
            
    lEkin = len(Ekin)
    if ((lEkin == lt) and (lEkin == len(Epot))):
        Etot = []; Etot_rel = []; Epot_rel = []; Ekin_rel = []
        start = int(lEkin/4.0)
        meanEk = np.mean(Ekin[start:])
        meanEp = np.mean(Epot[start:])
        stdEk = np.std(Ekin[start:])
        stdEp = np.std(Epot[start:])
        for i in range(lEkin):
            Etot.append((Ekin[i] + Epot[i]))
            
        meanEt = np.mean(Etot[start:])
        stdEt = np.std(Etot[start:])
        for i in range(lEkin):
            Etot_rel.append(Etot[i]/meanEt)
            Ekin_rel.append(Ekin[i]/meanEt)
            Epot_rel.append(Epot[i]/meanEt)
            
        Title = 'Kinetic energy. $ E_k = %.1f \pm %.2f [Kcal/mol]$' % (meanEk, stdEk)
        Xlabel = 'Time [%s]' % timeunits
        Ylabel = 'Energy $[Kcal/mole]$'
        legends = ['$E_k$']
        linestyles = ['r--']
        plt.figure()
        plt.plot(t,Ekin,linestyles[0])
        plt.title(r'Kinetic energy. $ E_k = %.1f \pm %.2f [Kcal/mol]$' % (meanEk, stdEk))
        plt.xlabel(Xlabel), plt.ylabel(Ylabel)
        plt.legend([r'$E_k$'],loc='upper right')
        fig_name = 'KineticEnergy_%s' % filename
        write_to_file(writetothisfile,Title,Xlabel,Ylabel,fig_name,legends,linestyles,[t],[Ekin])
        if (saveplot):
            fig_name = fig_name + '.png'
            print "saving %s" % fig_name
        
        Title = 'Potential energy. $ E_p = %.1f \pm %.2f [Kcal/mol]$' % (meanEp, stdEp)
        legends = ['$E_p$']
        plt.figure()
        plt.plot(t,Epot,linestyles[0])
        plt.title(r'Potential energy. $ E_p = %g \pm %g [Kcal/mol]$' % (meanEp, stdEp))
        plt.xlabel(Xlabel), plt.ylabel(Ylabel)
        plt.legend([r'$E_p$'],loc='upper right')
        fig_name = 'PotentialEnergy_%s' % filename
        write_to_file(writetothisfile,Title,Xlabel,Ylabel,fig_name,legends,linestyles,[t],[Epot])
        if (saveplot):
            fig_name = fig_name + '.png'
            print "saving %s" % fig_name
            
        Title = 'Energy in the system. $E_{tot} = %.1f \pm %.2f [Kcal/mol]$' % (meanEt,stdEt)
        legends = ['Ekin','Epot','Etot']
        linestyles = ['r--','b--','y-']
        Ylabel = 'Relative Energy'
        
        plt.figure()
        plt.plot(t,Ekin_rel,linestyles[0])
        plt.hold(True)
        plt.plot(t,Epot_rel,linestyles[1])
        plt.plot(t,Etot_rel,linestyles[2])
        plt.title(r'Energy in the system. $E_{tot} = %g \pm %g [Kcal/mol]$' % (meanEt,stdEt))
        plt.xlabel(Xlabel), plt.ylabel(Ylabel)
        plt.legend(legends,loc='center right')
        plt.hold(False)
        fig_name = 'Energy_%s' % filename
        write_to_file(writetothisfile,Title,Xlabel,Ylabel,fig_name,legends,linestyles,[t],[Ekin_rel,Epot_rel,Etot_rel])
        if (saveplot):
            fig_name = fig_name + '.png'
            print "saving %s " % fig_name
            plt.savefig(fig_name,format=Format)
        
    else:
        print "####################################"
        print "len(Ekin) = %g \nlen(Epot) = %g " % (len(Ekin),len(Epot))
    
        
    plt.show(showplot)
    plt.close('all')

def plotRDF(writetofile,path,filename,tc,timeunits,Format,saveplot,showplot):
    '''
    plot the evolution of the radial distribution function.\n
    plotDRF takes the path to the filename, and the filename as the two firs
    input args. tc is the time conversion factor that converts the timestep
    unit of the input file into the disiered time units.\n
    Choose the format of the output plots as Format = 'png','jpg','jpeg'...\n
    '''
    RDF_file = os.path.join(path,filename)
    string = "substance"
    if ((("Water") in filename) or (("water") in filename) or ("h2o" in filename)): string = "water"
    if ((("Port") in filename) or (("port") in filename) or ("Ca" in filename)): string = "portlandite"
    if (('CO2' in filename) or ("co2" in filename)): string = '$CO_2$'
            
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
    
    Title = 'Radial Distribution Function for %s' % string
    Xlabel = 'distance [\AA{}]'
    Ylabel = 'average number of particles'
    linestyles = []
    legends = []
    plt.figure()
    plt.hold(True)        
    for j in range(len(t)):
        plt.plot(BinCoord[j],RDF[j])
        legend = 't=%.1f%s' % (t[j],timeunits)
        legends.append(legend)
    plt.hold(False)
    plt.title(Title)
    plt.xlabel(Xlabel), plt.ylabel(Ylabel), plt.legend(legends,loc='upper left')
    fig_name = 'RDF_%s' % filename
    write_to_file(writetofile,Title,Xlabel,Ylabel,fig_name,legends,linestyles,BinCoord,RDF)
    if (saveplot):
        fig_name = fig_name + '.png'
        print "saving %s " % fig_name
        plt.savefig(fig_name,format=Format)
    
    legends = ['RDF(t=%g)' % t[-1]]
    plt.figure()
    plt.plot(BinCoord[-1],RDF[-1])
    plt.title(Title)
    plt.xlabel(Xlabel), plt.ylabel(Ylabel), plt.legend(legends,loc='upper left')
    fig_name = 'RDF_%g%s_%s' % (t[-1],timeunits,filename)
    write_to_file(writetofile,Title,Xlabel,Ylabel,fig_name,legends,linestyles,[BinCoord[-1]],[RDF[-1]])
    if (saveplot):
        fig_name = fig_name + '.png'
        print "saving %s " % fig_name 
        plt.savefig(fig_name,format=Format)
        
    plt.show(showplot)
    plt.close('all')
    
def plotDensityOfCO2andH2O(writetofile,path,co2_densfiles,water_densfiles,tc,timeunits,Format,saveplot,showplot):
 
    n = 4  # convinient to know the number of sample values in the datafile. Could be given as raw input, but usualy set as 4.

    print "#################################################"
    print " plot density of CO2 and H2O in system"
    for afile in co2_densfiles:
        print "co2_densfile: %s " % afile
        if (len(co2_densfiles) > 1):
            print "Error! more than one co2 density file"
    for afile in water_densfiles:
        print "h2o_densfile: %s " % afile
        if (len(water_densfiles) > 1):
            print "Error! more than one h2o density file"
            
    co2 = os.path.join(path,co2_densfiles[0])
    h2o = os.path.join(path,water_densfiles[0])
    
    co2dens = open(co2,'r'); co2dens.readline(); co2dens.readline(); co2labels = co2dens.readline(); co2labels = co2labels.split(); co2labels.pop(0)
    h2odens = open(h2o,'r'); h2odens.readline(); h2odens.readline(); h2olabels = h2odens.readline(); h2olabels = h2olabels.split(); h2olabels.pop(0)

    print "Labels: ", co2labels
    
    gobacktothispoint = co2dens.tell()  # we want to remember this point in the filestructure
    firstvals = co2dens.readline()
    firstvals = firstvals.split(); N = int(firstvals[1])
    co2dens.seek(gobacktothispoint)    # now go back to this point in the filestructure!
    
    
    t1 = np.zeros((n,1))             # time [ps]
    BIN = np.zeros((n,N))            # bin number [1:N+1]
    Coord = np.zeros((n,N))          # coord [Å]
    Ncount_co2 = np.zeros((n,N))     # number of particles counted in the bin
    Ndens_co2 = np.zeros((n,N))      # number density                         [particles/Å**3]
    Mdens_co2 = np.zeros((n,N))      # mass density                           [(g/mol)/Å**3]
    Ncount_h2o = np.zeros((n,N))     # number of particles counted in the bin
    Ndens_h2o = np.zeros((n,N))      # number density                        [particles/Å**3]
    Mdens_h2o = np.zeros((n,N))      # mass density                          [(g/mol)/Å**3]
    co2vx = np.zeros((n,N))          # velocity in x-dir
    co2vy = np.zeros((n,N))          # velocity in y-dir
    co2vz = np.zeros((n,N))          # velocity in z-dir
    h2ovx = np.zeros((n,N))          # h2o vel x
    h2ovy = np.zeros((n,N))          # h2o vel y
    h2ovz = np.zeros((n,N))          # h2o vel z
    
    for i in range(n):
        co2_line = co2dens.readline()
        h2o_line = h2odens.readline()
        time, nrows = co2_line.split() 
        t1[i,0] = (int(time)*tc)
        nrows = int(nrows)
        kk = 0
        for row in range(nrows):        
            co2_line = co2dens.readline()
            h2o_line = h2odens.readline()
            
            binnr_co2,coord_co2,ncount_co2,ndens_co2,mdens_co2,cvx,cvy,cvz = co2_line.split()
            binnr_h2o,coord_h2o,ncount_h2o,ndens_h2o,mdens_h2o,wvx,wvy,wvz = h2o_line.split()
            
            BIN[i,kk] = (int(binnr_co2))
            Coord[i,kk] = (float(coord_co2))
            Ncount_co2[i,kk] = (float(ncount_co2))
            Ndens_co2[i,kk]  = (float(ndens_co2))
            Mdens_co2[i,kk] = (float(mdens_co2))
            co2vx[i,kk] = float(cvx)
            co2vy[i,kk] = float(cvy)
            co2vz[i,kk] = float(cvz)
            Ncount_h2o[i,kk] = (float(ncount_h2o))
            Ndens_h2o[i,kk]  = (float(ndens_h2o))
            Mdens_h2o[i,kk] = (float(mdens_h2o))
            h2ovx[i,kk] = float(wvx)
            h2ovy[i,kk] = float(wvy)
            h2ovz[i,kk] = float(wvz)
            kk += 1
                        
    plt.figure()
    plt.hold(True)
    w_line = 'b-'
    c_line = 'r-'
    legends = []
    linestyles = []
    yarrays = []
    for i in range(n):
        if (i == 0):
            add = '*'
        if (i == 1):
            add = 'd'
        if (i == 2):
            add = 'x'
        if (i == 3):
            add = 'o'
        
        linestyles.append((w_line + add));linestyles.append((c_line + add))
        yarrays.append(Mdens_co2[i,:]); yarrays.append(Mdens_h2o[i,:])
        plt.plot(Coord[0,:],Mdens_co2[i,:],(c_line+add))
        plt.plot(Coord[0,:],Mdens_h2o[i,:],(w_line+add))
        carbdiox = r'$\rho_{CO2}-$'+ str(t1[i,0]) + timeunits
        water = r'$\rho_{H_2O}-$'+ str(t1[i,0]) + timeunits
        legends.append(carbdiox)       
        legends.append(water)

    Title = 'Density distribution across the system'
    Xlabel = 'distance [\AA{}]'
    Ylabel = 'Density [$g/cm^3$]'

    plt.hold(False)
    plt.title(Title)
    plt.xlabel('z [Aa]')
    plt.ylabel(r'Density $[g/cm^3]$')
    plt.legend(legends,loc='upper right')
    fig_name = "DensityDist_CO2_H2O"
    write_to_file(writetofile,Title,Xlabel,Ylabel,fig_name,legends,linestyles,[Coord[0,:]],yarrays)
    if (saveplot):
            fig_name = fig_name + "." + Format
            plt.savefig(fig_name, format=Format)
    
    plt.show(showplot)
    plt.close('all')
 
def plotDensityOfCO2(writetofile,path,files,tc,timeunits,Format,saveplot,showplot):
    '''
    Takes densityprofile output files generated with LAMMPS to display the 
    density profile through the system.
    
    Here we expect the input files to describe the density profile of CO2,
    the oxygen atoms and the carbon atoms in the system. These datasets
    are given in three different files. files = [carbDens_c,carbDens_o,carbDens_co2].
    We also expect the files to have the same structure and sizes.
    '''
    carbfile = None
    oxfile = None
    for afile in files:
        if (('_c' in afile) and ('_co2' not in afile)):
            carbon = os.path.join(path,afile)  
            carbfile = afile
        if ('_o' in afile):            
            oxygen = os.path.join(path,afile)
            oxfile = afile
        if ('_co2' in afile):
            co2 = os.path.join(path,afile)
            co2file = afile
    
    print "###############################################################"
    if(carbfile is not None): 
        print "carbon file:         %s " % carbfile
        carbdens = open(carbon,'r')

    if(oxfile is not None): 
        print "oxygen file:         %s " % oxfile
        oxdens = open(oxygen,'r')
        
    print "carbon dioxide file: %s " % co2file
    co2dens = open(co2,'r')
    
    if ((oxfile == None) and (carbfile == None )):
        print "ERROR!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
        print "we would in this function also like the density dist. of Carbon and Oxygen!"
        return 1
    
    co2dens.readline(), oxdens.readline(), carbdens.readline()  # comment
    co2dens.readline(), oxdens.readline(), carbdens.readline()  # Timestep Number-of-bins
    CO2labels = co2dens.readline()
    Olabels = oxdens.readline()
    Clabels = carbdens.readline()  # Bin Coord density/number density/mass
    
    CO2labels = CO2labels.split(); CO2labels.pop(0) # split on whitespaces and pop first element
    Olabels = Olabels.split(); Olabels.pop(0) # split on whitespaces and pop first element
    Clabels = Clabels.split(); Clabels.pop(0) # split on whitespaces and pop first element
    
    n = 4  # convinient to know the number of sample values in the datafile
    
    Nlabels = len(CO2labels)    
    if (Nlabels != len(Olabels)):
        print "Error! Number of columns in carbDens_c is %g and carbDens_co2 has %g" % (len(Olabels),n)
    if (Nlabels != len(Clabels)):
        print "Error! Number of columns in carbDens_c is %g and carbDens_co2 has %g" % (len(Clabels),n)        

    gobacktothispoint = co2dens.tell()  # we want to remember this point in the filestructure
    firstvals = co2dens.readline()
    firstvals = firstvals.split(); N = int(firstvals[1])
    co2dens.seek(gobacktothispoint)    # now go back to this point in the filestructure!
    
    #print "Samplevalues = %g. Number of rows = %g. Number of columns = %g" % (n,N, Ncolumns)
    
    t1 = np.zeros((n,1))             # time [ps]
    BIN = np.zeros((n,N))            # bin number [1:N+1]
    Coord = np.zeros((n,N))          # coord [Å]
    Ncount_co2 = np.zeros((n,N))       # water number of particles counted in the bin
    Ndens_co2 = np.zeros((n,N))        # water number density                         [particles/Å**3]
    Mdens_co2 = np.zeros((n,N))        # water mass density                           [(g/mol)/Å**3]
    Ncount_o = np.zeros((n,N))       # oxygen number of particles counted in the bin
    Ndens_o = np.zeros((n,N))        # oxygen number density                        [particles/Å**3]
    Mdens_o = np.zeros((n,N))        # oxygen mass density                          [(g/mol)/Å**3]
    Ncount_c = np.zeros((n,N))       # hydrogen number of particles counted in the bin
    Ndens_c = np.zeros((n,N))        # hydrogen number density                      [particles/Å**3]
    Mdens_c = np.zeros((n,N))        # hydrogen mass density                        [(g/mol)/Å**3]
    co2vx = np.zeros((n,N))             # water velocity in x-dir
    co2vy = np.zeros((n,N))             # water velocity in y-dir
    co2vz = np.zeros((n,N))             # water velocity in z-dir
    
    for i in range(n):
        line_co2 = co2dens.readline()
        line_o = oxdens.readline()
        line_c = carbdens.readline()
        time, nrows = line_co2.split()  # is the timestep for all three files, also same numer of rows/columns
        t1[i,0] = (int(time)*tc)
        nrows = int(nrows)
        kk = 0
        for row in range(nrows):
            line_co2 = co2dens.readline()
            line_o = oxdens.readline()
            line_c = carbdens.readline()
            # spce (water molecule, is out first. we use BIN and Coord only from this split()
            binnumber,coord,ncount_co2,ndens_co2,mdens_co2,vx,vy,vz = line_co2.split()
            BIN[i,kk] = (int(binnumber))
            Coord[i,kk] = (float(coord))
            Ncount_co2[i,kk] = (float(ncount_co2))
            Ndens_co2[i,kk]  = (float(ndens_co2))
            Mdens_co2[i,kk] = (float(mdens_co2))
            co2vx[i,kk] = float(vx)
            co2vy[i,kk] = float(vy)
            co2vz[i,kk] = float(vz)
            # oxygen is up next
            binnumbero,coordo,ncount_o,ndens_o,mdens_o,ovx,ovy,ovz = line_o.split()
            Ncount_o[i,kk] = (float(ncount_o))
            Ndens_o[i,kk]  = (float(ndens_o))
            Mdens_o[i,kk] = (float(mdens_o))
            # last put is hydrogen
            binnumberc,coordc,ncount_c,ndens_c,mdens_c,cvx,cvy,cvz = line_c.split()
            Ncount_c[i,kk] = (float(ncount_c))
            Ndens_c[i,kk]  = (float(ndens_c))
            Mdens_c[i,kk] = (float(mdens_c))
            kk += 1

    m_CO2 = CO2 = m_ox = ox = m_carb = carb = 0
    start = int((len(Mdens_co2[3]))/4.0)
    stop = int(3.0*start)
    counter = 0
    for i in range((stop-start)):
        kk = start + i
        a = Mdens_co2[3,kk]
        b = Mdens_o[3,kk]
        c = Mdens_c[3,kk]
        #print " %.5f    %.5f    %.5f " % (a,b,c)
        m_CO2 += a
        m_ox += b
        m_carb += c
        CO2 += a*a
        ox += b*b
        carb += c*c
        counter += 1
    m_CO2 = m_CO2/counter; m_ox = m_ox/counter; m_carb = m_carb/counter
    std_CO2 = np.sqrt(CO2/counter - m_CO2*m_CO2) 
    std_ox = np.sqrt(ox/counter - m_ox*m_ox)
    std_carb = np.sqrt(carb/counter - m_carb*m_carb)    

    rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
    ## for Palatino and other serif fonts use:
    #rc('font',**{'family':'serif','serif':['Palatino']})
    rc('text', usetex=True) 
    
    print "###############################################################"
    print "# ------------------ Densities across the pore -------------- #"
    print "# density in the pore center. CO2:      %.3f $ \pm $ %.4f  #" % (m_CO2,std_CO2)
    print "# density in the pore center. oxygen:   %.3f $ \pm $ %.4f  #" % (m_ox,std_ox)
    print "# density in the pore center. carbon:   %.3f $ \pm $ %.4f  #" % (m_carb,std_carb)
    print "###############################################################"
    
    # time for some plotting
    
    plt.figure()
    plt.hold(True)
    w_line = 'b-'
    h_line = 'y-'
    o_line = 'r-'
    legends = []
    min_W = 10     # default startvalue
    max_W = 0      # default startvalue 
    
    linestyles = []
    yarrays = []
    for i in range(n):
        if (i == 0):
            add = '*'
        if (i == 1):
            add = '-'
        if (i == 2):
            add = 'x'
        linestyles.append((w_line + add));linestyles.append((o_line + add));linestyles.append((h_line + add))
        yarrays.append(Mdens_co2[i,:]); yarrays.append(Mdens_o[i,:]); yarrays.append(Mdens_c[i,:])
        plt.plot(Coord[0,:],Mdens_co2[i,:],(w_line+add))
        plt.plot(Coord[0,:],Mdens_o[i,:],(o_line+add))
        plt.plot(Coord[0,:],Mdens_c[i,:],(h_line+add))
        water = r'$\rho_{CO2}$'+ str(t1[i,0]) + timeunits
        oxygen = r'$\rho_O-$'+ str(t1[i,0]) + timeunits
        hydrogen = r'$\rho_C-$'+ str(t1[i,0]) + timeunits
        legends.append(water)
        legends.append(oxygen)
        legends.append(hydrogen)
    min_w = min(Mdens_co2[i,:])
    max_w = max(Mdens_co2[i,:])
    if (min_w < min_W):
        min_W = min_w
    if (max_w > max_W):
        max_W = max_w

    Title = 'Density distribution across the pore\\newline Average $CO_2$ density inside the pore $ \rho _{CO_2} = %.2f \pm  %.3f [g/cm^3] $' % (m_CO2,std_CO2)
    Xlabel = 'distance [\AA{}]'
    Ylabel = 'Density [$g/cm^3$]'

    plt.hold(False)
    #plt.axis([0, 60, min_W, max_W])
    plt.title(r'Density distribution across the pore\newline Average density of carbon dioxide inside the pore $\rho _{CO_2} = $ %.2f $ \pm $ %.3f $ [g/cm^3] $' % (m_CO2,std_CO2))
    plt.xlabel('z [Aa]')
    plt.ylabel(r'Density $[g/cm^3]$')
    plt.legend(legends,loc='upper right')
    fig_name = "DensityDist_" + files[3]
    write_to_file(writetofile,Title,Xlabel,Ylabel,fig_name,legends,linestyles,[Coord[0,:]],yarrays)
    if (saveplot):
            fig_name = fig_name + "." + Format
            plt.savefig(fig_name, format=Format)
    
    plt.show(showplot)
    plt.close('all')

def plotDensity(writetofile,path,files,tc,timeunits,Format,saveplot,showplot):
    '''
    Takes densityprofile output files generated with LAMMPS to display the 
    density profile through the system.
    
    Here we expect the input files to describe the density profile of water,
    the oxygen atoms and the hydrogen atoms in the system. These datasets
    are given in three different files. files = [waterDens,oxDens,HyDens].
    We also expect the files to have the same structure and sizes.
    '''
    
    water = os.path.join(path,files[0])
    oxygen = os.path.join(path,files[1])
    hydrogen = os.path.join(path,files[2])
    
    print "###############################################################"
    print "water file:    %s " % files[0]
    print "oxygen file:   %s " % files[1]
    print "hydrogen file: %s " % files[2]

    waterdens = open(water,'r')  # the water file [h2o]    
    oxdens = open(oxygen,'r')  # the oxygen file
    hydens = open(hydrogen,'r')  # the hydrogen file
    
    waterdens.readline(), oxdens.readline(), hydens.readline()  # comment
    waterdens.readline(), oxdens.readline(), hydens.readline()  # Timestep Number-of-bins
    Wlabels = waterdens.readline()
    Olabels = oxdens.readline()
    Hlabels = hydens.readline()  # Bin Coord density/number density/mass
    
    Wlabels = Wlabels.split(); Wlabels.pop(0) # split on whitespaces and pop first element
    Olabels = Olabels.split(); Olabels.pop(0) # split on whitespaces and pop first element
    Hlabels = Hlabels.split(); Hlabels.pop(0) # split on whitespaces and pop first element
    
    n = 4  # convinient to know the number of sample values in the datafile
    
    Nlabels = len(Wlabels)    
    if (Nlabels != len(Olabels)):
        print "Error! Number of columns in Dens_ow is %g and Dens_water has %g" % (len(Olabels),n)
    if (Nlabels != len(Hlabels)):
        print "Error! Number of columns in Dens_hw is %g and Dens_water has %g" % (len(Hlabels),n)        

    gobacktothispoint = waterdens.tell()  # we want to remember this point in the filestructure
    firstvals = waterdens.readline()
    firstvals = firstvals.split(); N = int(firstvals[1])
    waterdens.seek(gobacktothispoint)    # now go back to this point in the filestructure!
    
    #print "Samplevalues = %g. Number of rows = %g. Number of columns = %g" % (n,N, Ncolumns)
    
    t1 = np.zeros((n,1))             # time [ps]
    BIN = np.zeros((n,N))            # bin number [1:N+1]
    Coord = np.zeros((n,N))          # coord [Å]
    Ncount_w = np.zeros((n,N))       # water number of particles counted in the bin
    Ndens_w = np.zeros((n,N))        # water number density                         [particles/Å**3]
    Mdens_w = np.zeros((n,N))        # water mass density                           [(g/mol)/Å**3]
    Ncount_o = np.zeros((n,N))       # oxygen number of particles counted in the bin
    Ndens_o = np.zeros((n,N))        # oxygen number density                        [particles/Å**3]
    Mdens_o = np.zeros((n,N))        # oxygen mass density                          [(g/mol)/Å**3]
    Ncount_h = np.zeros((n,N))       # hydrogen number of particles counted in the bin
    Ndens_h = np.zeros((n,N))        # hydrogen number density                      [particles/Å**3]
    Mdens_h = np.zeros((n,N))        # hydrogen mass density                        [(g/mol)/Å**3]
    wvx = np.zeros((n,N))             # water velocity in x-dir
    wvy = np.zeros((n,N))             # water velocity in y-dir
    wvz = np.zeros((n,N))             # water velocity in z-dir
    
    for i in range(n):
        line_w = waterdens.readline()
        line_o = oxdens.readline()
        line_h = hydens.readline()
        time, nrows = line_w.split()  # is the timestep for all three files, also same numer of rows/columns
        t1[i,0] = (int(time)*tc)
        nrows = int(nrows)
        kk = 0
        for row in range(nrows):
            line_w = waterdens.readline()
            line_o = oxdens.readline()
            line_h = hydens.readline()
            # spce (water molecule, is out first. we use BIN and Coord only from this split()
            binnumber,coord,ncount_w,ndens_w,mdens_w,vx,vy,vz = line_w.split()
            BIN[i,kk] = (int(binnumber))
            Coord[i,kk] = (float(coord))
            Ncount_w[i,kk] = (float(ncount_w))
            Ndens_w[i,kk]  = (float(ndens_w))
            Mdens_w[i,kk] = (float(mdens_w))
            wvx[i,kk] = float(vx)
            wvy[i,kk] = float(vy)
            wvz[i,kk] = float(vz)
            # oxygen is up next
            binnumbero,coordo,ncount_o,ndens_o,mdens_o,ovx,ovy,ovz = line_o.split()
            Ncount_o[i,kk] = (float(ncount_o))
            Ndens_o[i,kk]  = (float(ndens_o))
            Mdens_o[i,kk] = (float(mdens_o))
            # last put is hydrogen
            binnumberh,coordh,ncount_h,ndens_h,mdens_h,hvx,hvy,hvz = line_h.split()
            Ncount_h[i,kk] = (float(ncount_h))
            Ndens_h[i,kk]  = (float(ndens_h))
            Mdens_h[i,kk] = (float(mdens_h))
            kk += 1

    m_wp = wp = m_op = op = m_hp = hp = 0
    start = int((len(Mdens_w[3]))/4.0)
    stop = int(3.0*start)
    counter = 0
    for i in range((stop-start)):
        kk = start + i
        a = Mdens_w[3,kk]
        b = Mdens_o[3,kk]
        c = Mdens_h[3,kk]
        #print " %.5f    %.5f    %.5f " % (a,b,c)
        m_wp += a
        m_op += b
        m_hp += c
        wp += a*a
        op += b*b
        hp += c*c
        counter += 1
    m_wp = m_wp/counter; m_op = m_op/counter; m_hp = m_hp/counter
    std_wp = np.sqrt(wp/counter - m_wp*m_wp) 
    std_op = np.sqrt(op/counter - m_op*m_op)
    std_hp = np.sqrt(hp/counter - m_hp*m_hp)    

    rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
    ## for Palatino and other serif fonts use:
    #rc('font',**{'family':'serif','serif':['Palatino']})
    rc('text', usetex=True) 
    
    print "###############################################################"
    print "# ------------------ Densities across the pore -------------- #"
    print "# density in the pore center. water:    %.3f $ \pm $ %.4f  #" % (m_wp,std_wp)
    print "# density in the pore center. oxygen:   %.3f $ \pm $ %.4f  #" % (m_op,std_op)
    print "# density in the pore center. hydrogen: %.3f $ \pm $ %.4f  #" % (m_hp,std_hp)
    print "###############################################################"
    
    # time for some plotting
    
    plt.figure()
    plt.hold(True)
    w_line = 'b-'
    h_line = 'y-'
    o_line = 'r-'
    legends = []
    min_W = 10     # default startvalue
    max_W = 0      # default startvalue 
    
    linestyles = []
    yarrays = []
    for i in range(n):
        if (i == 0):
            add = '*'
        if (i == 1):
            add = '-'
        if (i == 2):
            add = 'x'
        linestyles.append((w_line + add));linestyles.append((o_line + add));linestyles.append((h_line + add))
        yarrays.append(Mdens_w[i,:]); yarrays.append(Mdens_o[i,:]); yarrays.append(Mdens_h[i,:])
        plt.plot(Coord[0,:],Mdens_w[i,:],(w_line+add))
        plt.plot(Coord[0,:],Mdens_o[i,:],(o_line+add))
        plt.plot(Coord[0,:],Mdens_h[i,:],(h_line+add))
        water = r'$\rho_w-$'+ str(t1[i,0]) + timeunits
        oxygen = r'$\rho_o-$'+ str(t1[i,0]) + timeunits
        hydrogen = r'$\rho_h-$'+ str(t1[i,0]) + timeunits
        legends.append(water)
        legends.append(oxygen)
        legends.append(hydrogen)
    min_w = min(Mdens_w[i,:])
    max_w = max(Mdens_w[i,:])
    if (min_w < min_W):
        min_W = min_w
    if (max_w > max_W):
        max_W = max_w

    Title = 'Density distribution across the pore\\newline Average water  $\rho _w = %.2f \pm %.3f [g/cm^3] $' % (m_wp,std_wp)
    Xlabel = 'distance [\AA{}]'
    Ylabel = 'Density [$g/cm^3$]'

    plt.hold(False)
    #plt.axis([0, 60, min_W, max_W])
    plt.title(r'Density distribution across the pore\newline Average density of water inside the pore $\rho _w = $ %.2f $ \pm $ %.3f $ [g/cm^3] $' % (m_wp,std_wp))
    plt.xlabel('z [Aa]')
    plt.ylabel(r'Density $[g/cm^3]$')
    plt.legend(legends,loc='upper right')
    fig_name = "DensityDist_" + files[0]
    write_to_file(writetofile,Title,Xlabel,Ylabel,fig_name,legends,linestyles,[Coord[0,:]],yarrays)
    if (saveplot):
            fig_name = fig_name + "." + Format
            plt.savefig(fig_name, format=Format)
    
    plt.show(showplot)
    plt.close('all')

def plotEnergy(writetofile, path, files, tc, timeunits, Format, saveplot, showplot):
    '''
    Plots the energy as a function of time. It could be necessary to take multiple files as
    argument. But this is not requirered.\n
    An energyfile contains the timestep in the first column, and then could contain columns
    of kinetic, potential, bond, angle, vdw (van der Waal) or coulumb potential energy.
    '''
    Etot = []; t = []
    Ekin = []; Epot = []; Ebond = []; Eangle = []; Evdw = []; Ecoul = []
    gen = 0
    for afile in files:
        pathtofile = os.path.join(path, afile)
        efile = open(pathtofile,'r')
        efile.readline()           # read the first line
        columns = efile.readline() # read description
        columns = columns.split()            # split on whitespaces
        columns.pop(0)             # remove the hash from the list!
        for line in efile:
            line = line.split() # slit on whitespaces
            index = 0
            for header in columns:
                if (("TimeStep") in header):
                    if (gen == 0):
                        # only use the timesteps from the first file!
                        t.append(float(line[index])*tc)
                if (("kin" or "Kin") in header):
                    Ekin.append(float(line[index]))
                if (("pot" or "Pot") in header):
                    Epot.append(float(line[index]))
                if (("coul" or "Coul")in header):
                    Ecoul.append(float(line(index)))
                if (("bond" or "Bond")in header):
                    Ebond.append(float(line(index)))
                if (("angle" or "Angel")in header):
                    Eangle.append(float(line(index)))
                if (("vdw" or "VDW")in header):
                    Evdw.append(float(line(index)))
                
                index += 1
        gen += 1
    #############################################################
    # Time for some plotting!
    ll = len(t)
    if (ll <= 1):
        print " ################################################ "
        print " # The first file %s contain <= 1 timesteps!" % files[0]
        print " # This will give you trouble!!! :-)            # "
    for i in range(ll):
        Etot.append((Ekin[i] + Epot[i]))

    start = int(ll/3.0)
    meanEt = np.mean(Etot[start:])
    meanEk = np.mean(Ekin[start:])
    meanEp = np.mean(Epot[start:])
    stdEt = np.std(Etot[start:])
    stdEk = np.std(Ekin[start:])
    stdEp = np.std(Epot[start:])
    
    rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
    ## for Palatino and other serif fonts use:
    #rc('font',**{'family':'serif','serif':['Palatino']})
    rc('text', usetex=True) 
    
    plt.figure()
    plt.plot(t,Ekin)
    plt.title(r'Kinetic energy. $ E_k = $ %.1f $ \pm $ %.2f $[Kcal/mol]$' % (meanEk,stdEk))
    plt.xlabel('Time [%s]' % timeunits),plt.ylabel('Energy [Kcal/mole]');plt.legend([r'$E_{k}(t)$'],loc='lower right')
    if (saveplot):    
        fig_name = 'KineticEnergy_' + files[0] + '.' + Format
        plt.savefig(fig_name, format = Format)
    
    plt.figure()
    plt.plot(t,Epot)
    plt.title(r'Potential energy. $ E_p = $ %.1f $ \pm $ %.2f $ [Kcal/mol] $' % (meanEp,stdEp))
    plt.xlabel('Time [%s]' % timeunits),plt.ylabel('Energy [Kcal/mole]');plt.legend([r'$E_{p}(t)$'],loc='upper right')
    if (saveplot):    
        fig_name = 'PotentialEnergy_' + files[0] + '.' + Format
        plt.savefig(fig_name, format = Format)
    
    plt.figure()
    Title = 'The energy of the system as a function of time\\newline $E_{tot} = $ %g $ \pm $ %g $ [Kcal/mol] $' % (meanEt,stdEt)
    Xlabel = 'Time [%s]' % timeunits
    Ylabel = 'Energy [Kcal/mole]'
    legends = ['$E_{tot}(t)$','$E_{k}(t)$','$E_{p}(t)$']
    linestyles = ['b-','r-','g-']
    plt.hold(True)
    plt.plot(t,Etot,linestyles[0])
    plt.plot(t,Ekin,linestyles[1])
    plt.plot(t,Epot,linestyles[2])
    plt.hold(False)
    plt.title(r'The energy of the system as a function of time\newline $E_{tot} = $ %g $ \pm $ %g $ [Kcal/mol] $' % (meanEt,stdEt))
    plt.xlabel(Xlabel),plt.ylabel(Ylabel);plt.legend([r'$E_{tot}(t)$',r'$E_{k}(t)$',r'$E_{p}(t)$'],loc='center right')
    fig_name = 'Energy_' + files[0]
    write_to_file(writetofile,Title,Xlabel,Ylabel,fig_name,legends,linestyles,[t],[Etot,Ekin,Epot])
    if (saveplot):    
        fig_name = fig_name + '.' + Format
        plt.savefig(fig_name, format = Format)
    
    plt.show(showplot)    
    plt.close('all')
    
def plotAutoCorr(writetofile,path,Acorr_files,tc,timeunits,Format,saveplot,showplot):
    '''
    Plot the auto correlation function
    '''
    for name in Acorr_files:
        filename = os.path.join(path,name)
        ofile = open(filename,'r')
        ofile.readline()
        labels = ofile.readline()
        labels = labels.split()
        labels.pop(0)
        
        time = []; acorr = []
        ncorr = len(labels)-1 # number of correlation functions
        for i in range(ncorr):
            acorr.append([])
        
        for line in ofile:
            line = line.split()
            time.append(float(line[0]))
            for j in range(ncorr):
                acorr[j].append(abs(float(line[j+1])))
        
        #######################################################
        #                 PLOTTING                            "
        Title = 'Autocorrelation function'
        Xlabel = 'Time steps'
        Ylabel = 'correlation'
        linestyles = []
        plt.figure()
        plt.hold(True)
        legends = []
        for i in range(ncorr):
            linestyle = 'b-*'
            if (i==1):
                linestyle = 'y-'
            else:
                linestyle = '-'
            linestyles.append(linestyle)
            plt.plot(time,acorr[i],linestyle)
            legend = 'AutoCorr%g' % i
            legends.append(legend)
        plt.hold(False)
        plt.title(Title)
        plt.xlabel(Xlabel), plt.ylabel(Ylabel), plt.legend(legends,loc='upper right')
        fig_name = 'AutoCorrelation_%s' % name
        write_to_file(writetofile,Title,Xlabel,Ylabel,fig_name,legends,linestyles,[time],acorr)
        if (saveplot):
            fig_name = 'AutoCorrelation_%s.png' % name
            plt.savefig(fig_name,format=Format)
            
        #plt.figure()
        #plt.plot(time,acorr[1])
        
        plt.show(showplot)
        plt.close('all')
        
def write_to_file(ofile,title,xlabel,ylabel,plotname,legends,linestyles,xarrays,yarrays):
    '''
    ofile = opened writable file\n
    title = string containing the title of the plot (string)\n
    ylabel, xlabel = y and x label name of the plot (string)\n
    plotname = filename of the plot. No whitespaces (string)\n
    legends = list of legends\n
    linestyles = list of linestyles\n
    xarray = list of x-arrays to plot against the yarrays. One single x-array can be given in x-array list if it is to be
    plotted against multiple yarrays. The xarray still has to be contained within a list!\n
    yarray = list of y-arrays. one for every x-array. If only one x-array is given, the y-arrays will be plotted against these!
    '''
    
    titleline = "title " + title + "\n"
    xlabelline = "xlabel " + xlabel + "\n"
    ylabelline = "ylabel " + ylabel + "\n"
    plotnameline = "plotname " + plotname + "\n"
    legendline = "legends"
    for legend in legends:
        legendline = legendline + " " + legend
    legendline += "\n"
    linestyleline = "linestyles"
    for linestyle in linestyles:
        linestyleline = linestyleline + " " + linestyle
    linestyleline += "\n"
    ofile.write(titleline)
    ofile.write(xlabelline)
    ofile.write(ylabelline)
    ofile.write(plotnameline)
    ofile.write(legendline)
    ofile.write(linestyleline)
    counter = 0
    for xarray in xarrays:
        counter += 1
        xarrayline = "x-array" + str(counter)
        for i in xarray:
            xarrayline = xarrayline + " " + str(i)
        xarrayline += "\n"
        ofile.write(xarrayline)
        
    counter = 0
    for yarray in yarrays:
        counter += 1
        yarrayline = "y-array"+str(counter)
        for i in yarray:
            yarrayline = yarrayline + " " + str(i)
        yarrayline += "\n"
        ofile.write(yarrayline)
        
    if (len(legends) != len(yarrays)):
        print "###################################################################"
        print "#                              Warning!                           #"
        print "# The number of legends no not match the number of yarrays given! #"
        print "# ----------------------------------------------------------------#"
        print "# Nlegends=%g, Nyarrays=%g, for plotname:                         #" % (len(legends),len(yarrays))
        print "#  %s " % (plotname)
        
def save_in_new_dir(base,filename,trueargs):

    oldpath = os.path.join(base,filename)
    print "General path                  : %s " % base 
    newpath = raw_input("Give new path from generalpath: ")
    if (os.path.isdir(newpath)):
        newabspath = (os.path.abspath(os.path.join(base,newpath)))
        print "The given (%s) path is an existing path." % (newabspath)
        savehere = raw_input(" Do you wish to save here(yes/no)? ")
        if (savehere in trueargs):
            print "oldpath = %s" % oldpath
            print "newpath = %s" % os.path.join(newabspath,filename)
        else:
            print "Ok, so so what is your problem?"
            save_in_new_dir(base,filename,trueargs)
    else:
        print "The path given (%s) is not an existing path." % os.path.join(base,newpath)
        createfolder = raw_input("Do you wish to create a folder here? ")
        if (createfolder in trueargs):
            newabspath = os.path.abspath(os.path.join(base,newpath))
            print "newpath = %s" % os.path.join(newabspath,filename)
            print "creating folder...."
            os.makedirs(newabspath) # use default mode
        else:
            print "Ok then. So what is your problem with this?"
            save_in_new_dir(base,filename,trueargs)

    print "########################################"
    newpath = os.path.abspath(os.path.join(newabspath,filename))
    os.rename(oldpath,newpath)    
    print "Alrighty then! Take care now, bye-bye then..."
    
def useage(arg):
    '''
    print usage info for "statistics_from_dumpfiles.py"
    '''
    print "#####################################################################"
    print "# Name of program:    %s   " % arg
    print "# Useage: python %s [showplots] [saveplots]" % arg
    print "#  showplots = (yes/no)\n#  saveplots = (yes/no)\n#"
    print "# Then just follow the instructions. ctrl+c for abortion"
    print "#####################################################################"


def main(argv):
    '''
    ###########################################################################\n
    #                   These values needs to be set                          #\n
    #                                                                         #\n
    #  generalpath = path to local area containing datafile hierarchy         #\n
    #  wp_path     = path from generalpath to data files                      #\n
    #  timeunits   = desiered timeunits (string shown on output figures)      #\n
    #  tc          = time conversion factor from timestep used to timeunits   #\n
    #  saveplot ?  = True/False                                               #\n
    #  showplot ?  = True/False                                               #\n
    #  Format      = 'png'/'jpg'/'jpeg'...                                    #\n
    ###########################################################################\n
    '''
    
    #generalpath = "/home/goran/lammps-28Jun14/examples/water_portlandite_system"    
    HELP = ['--help','-h','help','HELP','info']
    if (argv[1] in HELP):
        useage(argv[0])
        return 1
        
    generalpath = "/home/goran/lammps-28Jun14/examples/Abel_runs"
    
    print "This is the general path from where you may choose where to search for inputfiles:\n %s" % generalpath
    print "where do you wish to retrive files from?"
    wp_path = raw_input("Give path from here: ")
    if (wp_path[0] == '/'):
        wp_path = wp_path[1:]
    arg = raw_input("Give search keyword for keyword in filenames (None for no keyword): " )
    print "Keyword given: %s" %arg
    path = os.path.join(generalpath,wp_path)
    if (arg in ['None','none','no','nope','no keyword']):
        arg = None
    filenames = gothroughfiles(path,arg)
    timestep = raw_input("Give timestep used in simulation (fsec): ")
    try:
        timestep = float(timestep)
    except:
        print "Timestep given is not valid! Type(timestep) = <int> or <float>. You gave:"
        print timestep
        print "exiting..."
        return 1
    
    
    timeunits='ps'                   # desired output units of time
    ftop = 1000                      # [fsec] in [ps]
    psfactor = ftop/float(timestep)  # convert to [ps] 
    tc = 1.0/psfactor                # time conversionfactor
    nargs = len(argv)
    trueargs = ['y','Y','yes','Yes','Yeap','True','true','correct','yes ', 'Yes ', 'y ', 'Y ', 'YES ']
    falseargs = ['n','N','no','No','Nope','False','false','incorrect']
    saveplot = False      # save all plots
    showplot = True      # show all plots
    if (nargs > 1):
        if (argv[1] in trueargs):
            showplot = True
        elif (argv[1] in falseargs):
            showplot = False
        else:
            print "Invalid cmd arg given!! argv[1] = %s" % argv[1]
            print "Showplot is therefore set to defaultvalue. Showplot = True"
    if (nargs > 2):
        if (argv[2] in trueargs):
            saveplot = True
        elif (argv[2] in falseargs):
            saveplot = False
        else:
            print "Invalid cmd arg given!! argv[2] = %s" % argv[2]
            print "Showplot is therefore set to defaultvalue. Saveplot = False"


    Format = 'png'        # output format figures
    yes_plotScalarInfo = False       # Scalar info files present
    yes_plotRDF        = False       # RDF files present
    yes_plotDensity    = False       # Density files present
    yes_plotEnergy     = False       # Energy files present
    yes_plotAutoCorr   = False       # Auto Correlation function
    ###########################################################################
    
    '''
    There are different types of files in filename that we are going to do different operation on
    so we should find these types of files, and send them to the correct functions.    
    '''

    ########################################################
    # Sort the filenames in the directory specified by path:
    dens_files = []
    RDF_files = []
    ScalarInfo_files = []
    Energy_files = []
    Acorr_files = []
    
    for name in filenames:
        if ("Dens" in name):
            dens_files.append(name)
        if ("RDF" in name):
            RDF_files.append(name)
        if ("ScalarInfo" in name):
            ScalarInfo_files.append(name)
        if ((("Energy") in name) or ("energy" in name)):
            Energy_files.append(name)
        if (("acorr") in name or ("Acorr") in name):
            Acorr_files.append(name)

    print "########################################"    
    print "Energy files found:"
    for efile in Energy_files:
        print "    %s" % efile
    print "ScalarInfo files found:"
    for scifile in ScalarInfo_files:
        print "    %s" % scifile
    print "RDF files found:"
    for rdfile in RDF_files:
        print "    %s" % rdfile
    print "Density files found:"
    for densfile in dens_files:
        print "    %s" % densfile
    print "Auto correlation files found:"
    for acorrfile in Acorr_files:
        print "    %s" % acorrfile
    
    print "#######################################"
    print " What filetypes do you wish to process?\n You can give a set of characters for the filetypes you wish to include.\n s = ScalarInfo files\n r = RDF files\n d = density files\n e = Energy files\n a = AutoCorrelation funcction files\n "
    yes_plot = raw_input("Give characters for filetypes: ") 
    if ('s' in yes_plot):
        yes_plotScalarInfo = True 
    if ('r' in yes_plot):
        yes_plotRDF        = True 
    if ('d' in yes_plot):
        yes_plotDensity    = True  
    if ('e' in yes_plot):
        yes_plotEnergy     = True
    if('a' in yes_plot):
        yes_plotAutoCorr   = True 
        
    basefilename = os.path.basename(wp_path)
    today = datetime.date.today()
    writetofile = basefilename + "_plotdate_" + str(today) + ".dat"
    
    ofile = open(writetofile,'w')

    if (yes_plotScalarInfo):
        print "#####################################"
        print "plotScalarInfo..."
        if (len(ScalarInfo_files) < 1):
            print " There were no ScalarInfofiles found at this location: %s" % path
            print " No plots will be generated..."
        else:
            for filename in ScalarInfo_files:
                print filename
                plotScalarInfo(ofile,path,filename,tc,timeunits,Format,saveplot,showplot)
    
    if (yes_plotRDF):
        print "#####################################"
        print "plotRDF..."
        if (len(RDF_files) < 1):
            print " There were no RDF_files found at this location: %s" % path
            print " No plots will be generated..."
        else:
            for filename in RDF_files:
                print filename
                plotRDF(ofile,path,filename,tc,timeunits,Format,saveplot,showplot)


    water_densfiles = []; waterkey = ["water","Water","Wdens","oW","OW", "h2o","H2O"]
    co2_densfiles = []; carbkey = ["carbondioxide","co2","CO2"]
    
    for afile in dens_files:
        for key in waterkey:
            if ((key in afile) and (afile not in water_densfiles)):
                print "water density file: %s " % afile
                water_densfiles.append(afile)
    for afile in dens_files:
        for key in carbkey:
            if ((key in afile) and (afile not in water_densfiles) and (afile not in co2_densfiles)):
                print "CO2 density file: %s " % afile
                co2_densfiles.append(afile)
            
            
    if (yes_plotDensity):
        print "#####################################"
        print "plotDensity..."
        if ((len(co2_densfiles) >= 1) and 'carb' in co2_densfiles[0]):
            # CO2 files
            print "carb in densfile... proceeding..."
            plotDensityOfCO2(ofile,path,co2_densfiles,tc,timeunits,Format,saveplot,showplot)
        elif (len(co2_densfiles) >= 1 ):
            if (("co2" in co2_densfiles[0]) or ("h2o" in water_densfiles[0])):
                print "co2 in co2_densfiles or h2o in water_densfiles... proceeding..."
                plotDensityOfCO2andH2O(ofile,path,co2_densfiles,water_densfiles,tc,timeunits,Format,saveplot,showplot)

        else:
            # water density
            plotDensity(ofile,path,water_densfiles,tc,timeunits,Format,saveplot,showplot)
    
    if (yes_plotEnergy):
        print "#####################################"
        print "plotEnergy..."
        if (len(Energy_files) <= 0):
            print "Error plotting Energy files! No energyfiles found"
        else:
            plotEnergy(ofile,path,Energy_files,tc,timeunits,Format,saveplot,showplot)

    if(yes_plotAutoCorr):
        print "#####################################"
        print "plotAutoCorr..."
        plotAutoCorr(ofile,path,Acorr_files,tc,timeunits,Format,saveplot,showplot)

    print "name of output-plot-file: %s" % writetofile
    print " location of output-plot-file: %s" % os.getcwd()
    ofile.close()

    movefile = raw_input("Do you wish to move the file to a different location (yes/no)?: ")
    if (movefile in trueargs):
        save_in_new_dir(os.getcwd(),writetofile,trueargs)
    else:
        print "Alrighty then. Take care now, bye-bye then :-)"
    return 0

if (__name__ == "__main__"):
    import numpy as np
    import datetime
    import matplotlib.pyplot as plt
    from matplotlib import rc
    import sys, os
    sys.exit(main(sys.argv))

