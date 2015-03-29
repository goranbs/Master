# -*- coding: utf-8 -*-
"""
Created on Wed Jan 14 13:20:37 2015

@author: goran
"""



def displacement(ofile,timestep,Natoms,matrix0,Type,Time,frontname,pathtodatabase,systemsize):
    '''
    timestep,system_size,matrix0 contain the initial information about the system.
    we have to go through sortedfiles to find information about the particles at
    later timesteps.
    '''
    factor = 10**(-8)                                               # conversionfactor from [A^2/ps] -> [m^2/s]

    Lx = systemsize[1] -systemsize[0]    
    Ly = systemsize[3] -systemsize[2]
    Lz = systemsize[5] -systemsize[4]
    halfsystemsize = 0.5*(Lx+Ly+Lz)/3.0   # half of average system size
    
    time = []
    for t in Time:
        time.append(timestep*t/1000.0)   # to picoseconds
        #time.append(t)
        
    MSDX = [];MSDY = [];MSDZ = [];MSD = [];D_local = []
    dt = (time[1]-time[0])
    for T in Time:
        path = frontname + '.%r.txt' % (T)
        print path
        msdx = 0; msdy = 0; msdz = 0; msd = 0
        filename = os.path.abspath(os.path.join(pathtodatabase,path))
        t, Natoms, system_size, matrix, readstructure, entries, types = readfile(filename)
        counter = 0
        for i in range(Natoms):
            if (matrix['type'][i] == Type):
                ID = matrix['id'][i]
                initial_index = None
                for j in range(Natoms):
                    k = matrix0['id'][j]
                    if (k == ID):
                        initial_index = j
                        break
                
                dx = matrix['x'][i] - matrix0['x'][initial_index]
                dy = matrix['y'][i] - matrix0['y'][initial_index]
                dz = matrix['z'][i] - matrix0['z'][initial_index]
                msdx += (dx)**2
                msdy += (dy)**2
                msdz += (dz)**2
                ########################################
                #        MINIMUM IMAGE CONVENTION      #
                if (dx < -halfsystemsize): dx = dx + Lx
                if (dx > halfsystemsize): dx = dx - Lx
                if (dy < -halfsystemsize): dy = dy + Ly
                if (dy > halfsystemsize): dy = dy - Ly
                if (dz < -halfsystemsize): dz = dz + Lz
                if (dz > halfsystemsize): dz = dz - Lz
                ########################################
                msd += dx**2 + dy**2 + dz**2
                counter += 1
    
        D_local.append((msd/(counter*6*T)))
        MSDX.append(msdx/(6*counter));MSDY.append(msdy/(6*counter));MSDZ.append(msdz/(6*counter));MSD.append(msd/(6*counter))
        

    # ballistic_end = 100ps
    degree = 1
    start = int(round(len(MSD)/5.0))
    p = np.polyfit(time[start:],MSD[start:],degree)   # last 4/5 of the dataset
    f = np.polyval(p,time)    
    d_mean_lastpart = p[0]#*10**(-8)        # to get [m^2/s]. [A^2/ps] = 10**(-8)[m^2/s]

    p1 = np.polyfit(time[0:start],MSD[0:start],degree)
    f1 = np.polyval(p1,time)    
    d_mean_firstpart = p1[0]#*10**(8)       # to get [m^2/s]. [A^2/ps] = 10**(-8)[m^2/s]

    p2 = np.polyfit(time[0:int(round(start/2.0))],MSD[0:int(round(start/2.0))],degree)
    f2 = np.polyval(p2,time)    
    d_mean_initpart = p2[0]#*10**(8)        # to get [m^2/s]. [A^2/ps] = 10**(-8)[m^2/s]
    
    p3 = np.polyfit(time[0:100],MSD[0:100],degree)
    f3 = np.polyval(p3,time)
    d_mean_actual = p3[0]#*10**(8)          # to get [m^2/s]. [A^2/ps] = 10**(-8)[m^2/s]

    print "#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#"
    print " Nvalues in estimate 1 = 100"
    print " Nvalues in estimate 2 = %g" % (int(round(start/2.0)))
    print " Nvalues in estimate 3 = %g" % (int(round(start)))
    print " Nvalues in estimate 4 = %g" % (len(time[start:]))
    
    ###############################################################################
    #                             plotting
    #plt.close('all')
    plt.figure()                                                                # Figure 1
    Title = 'Mean square displacement'
    legends = ['$msd_{x}$','$msd_{y}$','$msd_{z}$','$msd$']#['$msd$']#
    Ylabel = 'mean square displacement $ [A^2] $'
    Xlabel = 'time $ [ps] $'
    pltname = 'MSD_bulkwater'
    linestyles = ['--r','--y','--k','o-b']#['--r'] #
    plt.hold(True)
    plt.plot(time,MSDX,linestyles[0])
    plt.plot(time,MSDY,linestyles[1])
    plt.plot(time,MSDZ,linestyles[2])
    plt.plot(time,MSD,linestyles[3],markevery=10)
    plt.title(Title)
    plt.xlabel(Xlabel)
    plt.ylabel(Ylabel)
    plt.legend(legends,loc='lower right')
    plt.hold(False)
    write_to_file(ofile,Title,Xlabel,Ylabel,pltname,legends,linestyles,[time],[MSDX,MSDY,MSDZ,MSD])
    
    
    from mpl_toolkits.axes_grid.inset_locator import zoomed_inset_axes, inset_axes, mark_inset

    plt.figure()
    ax = plt.axes()                                                             # Figure 2
    Title = 'Mean square displacement'
    legends = ['MSD','$f_1 %.3f $' % (d_mean_actual),'$f_2 %.3f $' % (d_mean_initpart),'$f_3 %.3f $' % (d_mean_firstpart),'$f_4 %.3f $' % (d_mean_lastpart)]
    pltname = 'MSD_bulkwater_mean'
    Xlabel = 'time $ [ps] $'
    Ylabel = 'mean square displacement $ [A^2] $'
    linestyles = ['b-','--y','--r','--m','--y']
    
    plt.hold(True)
    plt.plot(time,MSD,linestyles[0])
    plt.plot(time[0:int(round(len(f)/8.0))],f3[0:int(round(len(f)/8.0))],linestyles[1])  # very first part
    plt.plot(time[0:int(round(len(f)/6.0))],f2[0:int(round(len(f)/6.0))],linestyles[2])  # initpart
    plt.plot(time[0:int(round(len(f)/5.0))],f1[0:int(round(len(f)/5.0))],linestyles[3])  # firstpart
    plt.plot(time,f,linestyles[4])                                                       # lastpart

    plt.hold(False)
    plt.xlabel(Xlabel)
    plt.ylabel(Ylabel)
    plt.legend(legends,loc='lower right')
    plt.title(Title)
    
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()
    
    axin = inset_axes(ax, width="30%", height="35%", loc=8)
    plt.hold(True)
    axin.plot(time,MSD,linestyles[0])
    axin.plot(time,f3,linestyles[1],linewidth=3.0)
    axin.plot(time,f2,linestyles[2])
    #axin.plot(time,f1,linestyles[3])
    plt.hold(False)
    axin.set_xlim(500, 520)
    ya = 0.0
    yb = 4.0
    Npoints = 4
    points = np.linspace(ya,yb,Npoints)
    axin.set_ylim(ya, yb)
    axin.set_xticks([])
    axin.set_yticks(points)
    mark_inset(ax, axin, loc1=3, loc2=4, fc="none", ec="0.5")    
    
    write_to_file(ofile,Title,Xlabel,Ylabel,pltname,legends,linestyles,[time],[MSD,f])
    
    print "#################################################################"
    print "# Diffusion estimate f1   : D = %.10f " % (d_mean_actual)    
    print "# Diffusion estimate f2   : D = %.10f " % (d_mean_initpart)
    print "# Diffusion estimate f3   : D = %.10f " % (d_mean_firstpart) 
    print "# Diffusion estimate f4   : D = %.10f " % (d_mean_lastpart)

    plt.figure()                                                                # Figure 3
    Title = 'Time evolution of the local diffusion constant. $ D=(msd/6dt) $'
    legends = ['$D(t)$']; linestyles = ['-y']
    pltname = 'Diffusion_bulkwater'
    Ylabel = 'displacement $ [ %g m^2/s]$' % factor
    plt.hold(True)
    plt.plot(time,D_local,linestyles[0])    
    plt.hold(False)
    plt.legend(legends,loc='upper right')
    plt.title(Title)
    plt.xlabel(Xlabel)
    plt.ylabel(Ylabel)
    write_to_file(ofile,Title,Xlabel,Ylabel,pltname,legends,linestyles,[time],[D_local])

    plt.show()
    ######################## END OF DISPLACEMENT ########################

def storefileinfo(Time,Natoms,frontname,pathtodatabase):
    print "##############################################################"
    print "# Processing input files. Storing information in grand matrix"
    print "#\n"
    grandmatrix = []
    for time in Time:
        path = frontname + '.%r.txt' % (time)
        print "# Processing file: %s " % path
        filename = os.path.abspath(os.path.join(pathtodatabase,path))
        t, Natoms, system_size, matrix, readstructure, entries, types = readfile(filename) # this takes time processing !!
        info = []
        for i in range(Natoms):
            atom = [matrix['id'][i],matrix['type'][i], matrix['x'][i],matrix['y'][i],matrix['z'][i]]
            info.append(atom)
            
        grandmatrix.append(info)
    
    print "##############################################################"
    print "# Done processing..."
    return grandmatrix
    ######################## END OF STOREFILEINFO #######################

def TEST():
    Nestimates = 5
    Nstatefiles = 10
    Natoms = 2
    tstep = 100
    t0 = 100
    for estimate in range(Nestimates):
        stop = Nstatefiles + estimate
        print "#######################"
        print " estimate %g" % estimate
        print "# %g" % stop
        for step in range(Nstatefiles):
            Time = t0*(1 + estimate) + step*tstep
            print " file = %g, step = %g" % (Time, step)
            for i in range(Natoms):
                print "atom"
    ######################## END OF TEST ###############################
    
    
def diffusion(ofile,timestep,Natoms,Nstatefiles,GM,Type,time,systemsize):
    '''
    Diffusion should calculate the diffusion constant for the system by calculating
    the msd for Nstatefiles, starting at statefile 0, and then do the same calculation
    for Nstatefiles from statefile 1, and so on untill we are out of statefiles. Then
    diffusion calculates the mean diffusion coeff and its standard deviation.
    
    '''
    divisor = 3.0                         # Divide msd into divisor number of areas. Use last parts in estimate. e.g divisor = 3, use last 2/3 of data in estimate.
    t0 = time[0]
    for i in range(len(time)):
        time[i] = time[i] - t0
    
    Lx = systemsize[1] -systemsize[0]    
    Ly = systemsize[3] -systemsize[2]
    Lz = systemsize[5] -systemsize[4]
    halfsystemsize = 0.5*(Lx+Ly+Lz)/3.0   # half of average system size
    
    print "# Running estimates on the diffusion coeff..." 
    factor = 10**(-8)                                               # conversionfactor from [A^2/ps] -> [m^2/s]
    degree = 1                                                      # Linear fitting. (polyfit degree 1) 
    Ntotalstatefiles = len(time)                                    # Number of statefiles to use in simulation
    Nestimates = int(np.floor(Ntotalstatefiles/float(Nstatefiles))) # Number of estimates on the diffusion coeff
    
    dt = time[1] - time[0]
    print "# dt=%g [ps]" % dt
    D = []; Dx = []; Dy = []; Dz = []; D_line = [];diffusion = []   # estimates on the diffusion constant
    MeanSquareDisplacement = []
    for estimate in range(Nestimates):
        stop = Nstatefiles + estimate                               # What statefileindex to stop on
        print "######################## estimate nr %g" % (estimate +1)
        initialmx = GM[estimate]

        MSDX = [];MSDY = [];MSDZ = [];MSD = [];dvalues = [];        # Mean square displacement
        for step in range(Nstatefiles):
            nextmx = step + estimate
            counter=0;msdx=0; msdy=0;msdz=0;msd=0
            for i in range(Natoms):
                if (GM[nextmx][i][1] == Type):                  # only track correct atom types
                    ID = GM[nextmx][i][0]                       # atom ID 
                    ix = None
                    for j in range(Natoms):
                        k = initialmx[j][0]
                        if (k == ID):                           # find the same atom at the initial timestep
                            ix = j
                            break
                    
                    dx = GM[nextmx][i][2] - initialmx[ix][2]    # Aangstrom squared
                    dy = GM[nextmx][i][3] - initialmx[ix][3]
                    dz = GM[nextmx][i][4] - initialmx[ix][4]
                    msdx += dx**2
                    msdy += dy**2
                    msdz += dz**2
                    ########################################
                    #        MINIMUM IMAGE CONVENTION      #
                    if (dx < -halfsystemsize): dx = dx + Lx
                    if (dx > halfsystemsize): dx = dx - Lx
                    if (dy < -halfsystemsize): dy = dy + Ly
                    if (dy > halfsystemsize): dy = dy - Ly
                    if (dz < -halfsystemsize): dz = dz + Lz
                    if (dz > halfsystemsize): dz = dz - Lz
                    ########################################

                    msd += dx*dx + dy*dy + dz*dz
                    counter += 1
            dvalues.append(msd/(6*counter*dt))
            MSDX.append(msdx/(6*counter));MSDY.append(msdy/(6*counter));MSDZ.append(msdz/(6*counter));MSD.append(msd/(6*counter));
        
        # Make an estimate on the diffusion constant:
        extra = int(round((stop - estimate)/divisor))
        
        p = np.polyfit(time[(estimate+extra):stop],MSD[extra:],degree)
        px = np.polyfit(time[(estimate+extra):stop],MSDX[extra:],degree)
        py = np.polyfit(time[(estimate+extra):stop],MSDY[extra:],degree)
        pz = np.polyfit(time[(estimate+extra):stop],MSDZ[extra:],degree)
        D_line.append(np.polyval(p,time[estimate:stop]))  # should fit the MSD
        D.append(p[0]); Dx.append(px[0]);Dy.append(py[0]);Dz.append(pz[0]); diffusion.append(dvalues)
        MeanSquareDisplacement.append(MSD)

    #print "# counter=%g" % counter
    meanD = np.mean(D); stdD = np.std(D)
    meanDx = np.mean(Dx); stdDx = np.std(Dz)
    meanDy = np.mean(Dy); stdDy = np.std(Dy)
    meanDz = np.mean(Dz); stdDz = np.std(Dz)
    
    print "################################################################"
    print "# Average of the estimates on the diffusion coeff.\n# Estimate the diffusion coeff for %g time intervalls" % (len(D))
    print "# meanD  = %.4f pm %.5f x %g [m^2/s]" % (meanD,stdD,factor)
    print "# meanDx = %.4f pm %.5f x %g [m^2/s]" % (meanDx,stdDx,factor)
    print "# meanDy = %.4f pm %.5f x %g [m^2/s]" % (meanDy,stdDy,factor)
    print "# meanDz = %.4f pm %.5f x %g [m^2/s]" % (meanDz,stdDz,factor)
    print "################################################################"
    
    lD = len(D)
    Dmean = np.linspace(meanD,meanD,lD)
    Dxmean = np.linspace(meanDx,meanDx,lD)
    Dymean = np.linspace(meanDy,meanDy,lD)
    Dzmean = np.linspace(meanDz,meanDz,lD)
    
    ##########################################################################################
    #                                Plotting
    
    line = np.linspace(1,len(D),len(D))
    Title = ' Estimates on the diffusion coefficient\n$D=%.3f \\pm %.4f [ %g m^2/s]$' % (meanD,stdD,factor)
    legends = ['$D$','$mean(D)$']#['$D$','$D_x$','$D_y$','$D_z$','$mean(D)$','$mean(D_x)$','$mean(D_y)$','$mean(D_z)$']
    linestyles = ['ob','--b']#['ob','vr','dm','xk','--b','--r','--m','--k']
    Xlabel = 'estimate nr'; Ylabel = 'D $[%g m^2/s] $' % factor
    pltname = 'Diffusion_coeff_bulkwater'
    plt.figure()
    plt.hold(True)
    plt.plot(line,D,linestyles[0])
    #plt.plot(line,Dx,linestyles[1])
    #plt.plot(line,Dy,linestyles[2])
    #plt.plot(line,Dz,linestyles[3])
    plt.plot(line,Dmean,linestyles[1])
    #plt.plot(line,Dxmean,linestyles[5])
    #plt.plot(line,Dymean,linestyles[6])
    #plt.plot(line,Dzmean,linestyles[7])
    plt.hold(False)
    plt.title(Title)
    plt.legend(legends,loc='upper left')
    plt.xlabel(Xlabel)
    plt.ylabel(Ylabel)
    Title = ' Estimates on the diffusion coefficient\\n$D=%.3g \\pm %.4g [m^2/s]$' % (meanD*factor,stdD*factor)
    #write_to_file(ofile,Title,Xlabel,Ylabel,pltname,legends,linestyles,[line],[D,Dx,Dy,Dz,Dmean,Dxmean,Dymean,Dzmean])
    write_to_file(ofile,Title,Xlabel,Ylabel,pltname,legends,linestyles,[line],[D,Dmean])
    
    if (Nestimates < 15):
        plt.figure()
        Title = 'Mean square displacement for bulk water'
        Xlabel = 'time $ [ps] $'
        Ylabel = 'msd $  [%g A^2] $' % factor
        legends = []; Legends = []
        pltname = 'MSD_time_bulkwater'
        linestyles = []
        colour1 = ['-rx','-bx','-gx','-yx','-kx','-mx','-rv','-bv','-gv','-yv','-kv','-mv']
        colour2 = ['r','b','g','y','k','m','r','b','g','y','k','m']
        xarray = []; yarray = []
        plt.hold(True)
        for est in range(Nestimates):
            kk = est
            if (est >= len(colour1)): kk = 0
            start = est*Nstatefiles
            stop = Nstatefiles + start
            leg = '$D_{%g}=%.3g$' % (est+1,D[est])
            legend1 = '$estimate_{%g}$' % (est+1)
            legend2 = '$fit_{%g}$' % (est+1)
            legends.append(legend1)
            legends.append(legend2)
            Legends.append(leg)
            Legends.append(legend2)
            style1 = str(colour1[kk])
            style2 = '--' + str(colour2[kk])
            plt.plot(time[start:stop],MeanSquareDisplacement[est],style1)
            xarray.append(time[start:stop])
            yarray.append(MeanSquareDisplacement[est])
            plt.plot(time[start:stop],D_line[est],style2)
            linestyles.append(style1);linestyles.append(style2);
            xarray.append(time[start:stop])
            yarray.append(D_line[est])
        plt.hold(False)
        #annotate('local max', xy=(2, 1), xytext=(3, 1.5),arrowprops=dict(facecolor='black', shrink=0.05),)
        plt.title(Title); plt.xlabel(Xlabel);plt.ylabel(Ylabel); plt.legend(Legends,loc='lower right')
        
        write_to_file(ofile,Title,Xlabel,Ylabel,pltname,legends,linestyles,xarray,yarray)

    sumMSD = 0
    for msd in MeanSquareDisplacement:
        sumMSD += np.array(msd)
    
    len_sumMSD = len(sumMSD)
    start = int(round(len_sumMSD/divisor))
    sumMSD /= len(MeanSquareDisplacement)
    p = np.polyfit(time[start:len_sumMSD],sumMSD[start:],1)
    f = np.polyval(p,time[0:len_sumMSD])
    
    diff = []
    stdev = 0
    for i in range(len_sumMSD):
        value = abs(f[i] - sumMSD[i])
        diff.append(value)
        stdev += value
    
    stdev = (stdev/len_sumMSD)
        
    print "# 2.nd estimate on the diffusion coeff.\n# Taking the average of the msd and then estimate the diffusion coeff:"
    print "# Diffusion constant D=%.4f pm %.5f x %g [m^2/s]" % (p[0],stdev,factor)
    
    plt.figure()
    Title = 'Average of the sum of estimates on D\n$D=%.3f \pm %.4f [ %g m^2/s]$' % (p[0],stdev,factor)
    legends = ['msd','linear fit']; linestyles = ['b-','--r']
    pltname = 'MSD_averaged_sum_fitted_lines'
    plt.hold(True)
    plt.errorbar(time[0:len_sumMSD],sumMSD,yerr=diff)
    plt.plot(time[0:len_sumMSD],f,'--r')
    plt.hold(False)
    plt.title(Title); plt.xlabel(Xlabel); plt.ylabel(Ylabel)
    plt.legend(legends,loc='upper left')
    Title = 'Average of the sum of estimates on D\\newline$D=%.3g \\pm %.4g [m^2/s]$' % (p[0]*factor,stdev*factor)
    write_to_file(ofile,Title,Xlabel,Ylabel,pltname,legends,linestyles,[time[0:len(sumMSD)]],[sumMSD,f])
    
    plt.show()
    
    ######################## END OF DIFFUSION #######################
   
def bulkdensity(ofile,timestep,mass_of_system,Time,frontname,pathtodatabase):
    '''
    masses[i] contains the mass of atom in types[i].
    systemSize = [xlo xhi ylo yhi zlo zhi]
    GM is a grand matrix containing all atom IDs, types and positions.
    The bulkdensity function calculates the density of a bulk system
    at every timestep, and plots the density as a function of time
    mass_of_system
    '''
    avogadro_aangstrom_relation = 10.0/6.022
    density = []; volume = []
    i = 0
    for time in Time:
        path = frontname + '.%r.txt' % (time)
        print "# Processing file: %s " % path
        filename = os.path.abspath(os.path.join(pathtodatabase,path))
        t, Natoms, system_size, matrix, readstructure, entries, types = readfile(filename) # this takes time processing !!
        dx = system_size[1] - system_size[0]  # Aangstrom
        dy = system_size[3] - system_size[2]  # Aangstrom
        dz = system_size[5] - system_size[4]  # Aangstrom
        vol = dx*dy*dz                        # volume of system [Aangstrom^3]
        volume.append(vol)
        density.append((mass_of_system/vol)*avogadro_aangstrom_relation)
        i += 1
    
    t = []
    for time in Time:
        t.append(time*timestep/1000.0) # converts to pico-seconds
    
    if "carb" in frontname:
        dens = 1.7966*10**(-3)                                    # g/cm^3
        benchmarkline = np.linspace(dens,dens,5)        # carbon dioxide benchmarkline for density
        pltname = "Density_cardondioxide"
    elif "water" in frontname:
        dens = 0.998                                              # g/cm^3
        benchmarkline = np.linspace(dens,dens,5)          # water benchmarkline for density
        pltname = "Density_water"
    elif "quartz" in frontname:
        dens = 2.648                                              # g/cm^3
        benchmarkline = np.linspace(dens,dens,5)
        pltname = "Density_quartz"
    elif "portlandite" in frontname:
        dens = 2.24      # g/cm^3
        benchmarkline = np.linspace(dens,dens,5)
        pltname = "Density_portlandite" 
    else:
        print "ERRORERRORERRORERRORERRORERRORERRORERRORERRORERRORERRORERRORERROR"
        print "Error!\n Could not recognize either 'carb', 'water', 'quartz or 'portlandite' in frontname"

    
    start = int(round(0.1*len(density)))    
    meanDens = np.mean(density[start:]); stdDens = np.std(density[start:])
    line = np.linspace(meanDens,meanDens,6)
    Title = 'Density profile\\newline$\\rho = %.4g \\pm %.5g [g/cm^3]$' % (meanDens,stdDens)
    Xlabel = 'time $ [ps] $'; Ylabel = 'density $ [g/cm^3] $'
    legends = ['$\\rho$','$\\rho_{avg}$','$\\rho_{exp}$']; linestyles = ['g--','b-d','r-o']
    tstop = t[-1]
    tstart = t[0]
    dt = tstop - tstart
    line_points = [tstart,(tstart + dt/5), (tstart + 2*dt/5), (tstart + 3*dt/5), (tstart + 4*dt/5), tstop]
    benchmarkline_points = [tstart,(tstart + dt/4), (tstart + dt/2), (tstart + 3*dt/4), tstop]
    plt.figure()
    plt.hold(True)
    plt.plot(t,density,'g--')
    plt.plot(line_points,line,'b-d', linewidth=3.0)
    plt.plot(benchmarkline_points,benchmarkline,'r-o',linewidth=2.0)
    plt.hold(False)
    plt.title(Title)
    plt.legend(legends)
    plt.xlabel(Xlabel); plt.ylabel(Ylabel)

    meanVol = np.mean(volume); stdVol = np.std(volume)
    Title = '$V_{system} = %g \\pm %g A^3$' % (meanVol,stdVol)
    Xlabel = 'time [ps]'; Ylabel = 'volume'
    plt.figure()
    plt.plot(t,volume,'-b')
    plt.title(Title); plt.xlabel(Xlabel);plt.ylabel(Ylabel)
    plt.legend(['volume'],loc='upper left')
        
    plt.show()
    write_to_file(ofile,Title,Xlabel,Ylabel,pltname,legends,linestyles,[t,line_points,benchmarkline_points],[density,line,benchmarkline])

    # END ------------------------#
    
def RDF(ofile,Time,Natoms,systemsize,GM,info_atom1,info_atom2,path):
    '''
    Calculates the average Radial Distribution Function for a set of inputfiles.
    The RDF is calculated as number of Type2 as function of Type1. 
    '''
    Type1 = info_atom1[0]; Type2 = info_atom2[0]
    type1 = info_atom1[1]; type2 = info_atom2[1]
    Nmolec1 = info_atom1[2]; Nmolec2 = info_atom2[2]
    
    Lx = systemsize[1] -systemsize[0]    
    Ly = systemsize[3] -systemsize[2]
    Lz = systemsize[5] -systemsize[4]
    halfsystemsize = 0.5*(Lx+Ly+Lz)/3.0   # half of average system size
    Ncontainers = 100                     # number of atomcontainers/bins
    #dr = 0.1    
    dr = halfsystemsize/Ncontainers       # size of bin [Angstrom]
    
    atomcontainer = np.zeros((Ncontainers,1))
    g = np.zeros((Ncontainers,1))
    r = np.linspace(0,Ncontainers*dr,Ncontainers)
    
    #lTime = len(Time)
    lTime = 5
    print "# Averaging over timesteps..."
    for t in range(lTime):
        print " # Timestep ( %g / %g ) " % (t+1,lTime)
        for i in range(Natoms):
            if (GM[t][i][1] == Type1):
                for j in range(Natoms):
                    if (GM[t][j][1] == Type2):
                        #if (j != i ):
                            dx = GM[t][i][2] - GM[t][j][2]
                            dy = GM[t][i][3] - GM[t][j][3]
                            dz = GM[t][i][4] - GM[t][j][4]
                            
                            dx = minimumImageConvention(dx,halfsystemsize,Lx)
                            dy = minimumImageConvention(dy,halfsystemsize,Ly)
                            dz = minimumImageConvention(dz,halfsystemsize,Lz)
                            
                            dr2 = np.sqrt(dx*dx + dy*dy + dz*dz)
                            container = int(round(dr2/dr))
                            if (container < Ncontainers):
                                atomcontainer[container] += 1

    for i in range(Ncontainers):
        atomcontainer[i] = atomcontainer[i]/lTime
    
    factor = (4.0/3.0)*3.14159265359
    sqNumdens = Nmolec1*(Nmolec2/(Lx*Ly*Lz))
    index_of_max_value = 0
    max_value = 0
    for i in np.arange(1,Ncontainers,1):
        V1 = factor*r[i-1]**3
        V2 = factor*r[i]**3
        dV = (V2 - V1)
        g[i] = (atomcontainer[i]/(dV*sqNumdens))
        if g[i] > max_value:
            max_value = g[i]
            index_of_max_value = i
    
    vmdfilename = 'vmd_rdf_' + type1 + type2 + '.dat'
    fullpath = os.path.join(path,vmdfilename)
    try:
        vmdfile = open(fullpath,'r')
    except:
        vmdfile = None

    plt.figure()
    Title = 'Radial Distribution %s - %s' % (type1,type2)
    Ylabel = 'g(r) ' + type1 + '-' + type2
    Xlabel = 'r [A]'
    plt.hold(True)
    plt.plot(r,g,'-b')                      # plotting the main line
    r_line = np.linspace(r[index_of_max_value],r[index_of_max_value],5)
    vertical_line = np.linspace(0,max_value,5)
    if (vmdfile != None):
    #if (vmdfile == None):        
        xx = []; gg = []; kk = []; gk = []
        legends = ['$g(r)$','$g(r)_{vmd}$','$max$']
        linestyles = ['-b','-r','k--']
        pltname = 'rdf_bulkwater_%s-%s' % (type1,type2)
        for line in vmdfile:
            #aline = vmdfile.readline()
            x,gr,k = line.split()
            xx.append(float(x)); gg.append(float(gr)); kk.append(float(k))
        for i in range(len(g)):
            gk.append(g[i][0])
        write_to_file(ofile,Title,Xlabel,Ylabel,pltname,legends,linestyles,[r,xx,r_line],[gk,gg,vertical_line])
        plt.plot(xx,gg,linestyles[1])   # plotting the vmd-file line

    #if (vmdfile != None):
    if (vmdfile == None):
        legends = ['$g(r)$','$max$']
        linestyles = ['-b','k--']

    position = (max_value/2.0, r[index_of_max_value])
    print position
    print type1, type2
    plt.annotate('r=%.3f' % (r[index_of_max_value]), xy=position,xytext=position, arrowprops=None,style='italic')
    plt.plot(r_line,vertical_line,'k--')
    plt.hold(False)
    plt.axis([0,r[Ncontainers-1],0,4])
    plt.title(Title); plt.xlabel(Xlabel);plt.ylabel(Ylabel);plt.legend(legends,loc='upper right')
    plt.show()
    
    
    ##################### END OF Radial Distribution Calculation ##################################    
    
def simplecalc(frontname,database,initialmx,Natoms,Time,systemsize):

    Lx = systemsize[1] -systemsize[0]    
    Ly = systemsize[3] -systemsize[2]
    Lz = systemsize[5] -systemsize[4]
    
    halfsystemsize = 0.5*(Lx+Ly+Lz)/3.0   # half of average system size
    
    rsquared = 0
    msd = []
    for t in Time:
        name = frontname + '.%r.txt' % t
        filename = os.path.abspath(os.path.join(database,name))
        time,Natoms,system_size,mx,readstructure,entries,types = readfile(filename)
        Noxygen = 0
        for i in range(Natoms):
            if (mx['type'][i] == 1):
                for j in range(Natoms):
                    if (initialmx['id'][j] == mx['id'][i]):
                        dx = mx['x'][i] - initialmx['x'][j]
                        dy = mx['y'][i] - initialmx['y'][j]
                        dz = mx['z'][i] - initialmx['z'][j]
                        ########################################
                        #        MINIMUM IMAGE CONVENTION      #
                        if (dx < -halfsystemsize): dx = dx + Lx
                        if (dx > halfsystemsize): dx = dx - Lx
                        if (dy < -halfsystemsize): dy = dy + Ly
                        if (dy > halfsystemsize): dy = dy - Ly
                        if (dz < -halfsystemsize): dz = dz + Lz
                        if (dz > halfsystemsize): dz = dz - Lz
                        ########################################
                        dr2 = dx*dx + dy*dy + dz*dz
                        Noxygen += 1
                        rsquared += dr2
    
        msd.append(rsquared/Noxygen)
        
    plt.figure()
    plt.plot(Time,msd)
    plt.show()
    
def read_inputfile(inputfilename):
    infile = open(inputfilename,'r') # open for reading
    Nestimates = None; Type = 1; Timestep = None
    InitialFile = None; PathToDatabase = None; arg = None
    for line in infile:
        someline = infile.readline()
        someline = someline.split()
        if (someline[0] == 'Nestimates'): 
            try:
                Nestimates = int(someline[1])
            except:
                print "Nestimates value cannot be converted to int"
        if (someline[0] == 'Type'):
            try:
                Type = int(someline[1])
            except:
                print "Type cannot be converted to int"
            
        if (someline[0] == 'Timestep'):
            try:
                Timestep = float(someline[1])
            except:
                print "Timestep value not recognized!"
        if (someline[0] == 'InitialFile'): InitialFile = someline[1]
        if (someline[0] == 'PathToDatabase'): PathToDatabase = someline[1]
        if (someline[0] == 'arg'): arg = someline[1]
    values = [Nestimates,Type,Timestep,arg,InitialFile,PathToDatabase]
    if (None in values):
        print "Failed to provide input arguments!"
        print values
    return values
    
def minimumImageConvention(delta,halfsystemsize,systemlength):
    ########################################
    #        MINIMUM IMAGE CONVENTION      #
    if (delta < -halfsystemsize): delta = delta + systemlength
    if (delta > halfsystemsize): delta = delta - systemlength
    ########################################
    return delta

def main():
    
    #calculation = "TraPPE"       # diffusion calc for the TraPPE co2 system
    #calculation = "EMP2"         # diffusion calc for the EMP2 co2 system
    #calculation = "FPF"          # diffusion calc for the FPF co2 system
    calculation = "Bulk water"   # diffusion for bulk water
    #calculation = "Portlandite"  # diffusion for Portlandite cube
    
    Nstatefiles = 100           # Number of statefiles to use in every estimate
    timestep = 2.0              # timestep used in MD simulation. [fsec]
    
    if (calculation == "Bulk water"):
        print " Running for bulk water!!!"
        '''
        startingtime = 200000
        pathtodatabase = '/home/goran/lammps-28Jun14/examples/Abel_runs/water/pure_H2O/Small_bulk_system' # bulk water
        initialfile = 'dump.water_nvt_fromRestart.250000.txt'
        frontname = "dump.water_nvt_fromRestart"
        arg = "water_nvt_fromRestart"
        outputfilename = 'diffusion_bulkwater_plotshit.dat'
        Type = 1                    # What atom type to follow
        '''
        startingtime = 250100
        pathtodatabase = "/home/goran/lammps-28Jun14/theproject/water/run3_300K_100ps_nvt" # Bulk water
        initialfile = "dump.water_nvt.250100.txt" # Filename of initial state
        frontname = "dump.water_nvt"
        arg = None
        outputfilename = "diffusion_bulkwater_plotshit.data"
        Type = 1

    if (calculation == "TraPPE"):
        last_timestep = 250000
        startingtime = 200000
        Nsteps = last_timestep - startingtime
        pathtodatabase = '/home/goran/lammps-28Jun14/theproject/carbondioxide/TraPPE/run_2015_02_27'         # bulk carbon dioxide, TraPPE model
        initialfile = 'dump.bulk_carbondioxide_TraPPE.0.txt'
        frontname = "dump.bulk_carbondioxide_TraPPE"
        arg = "bulk_carbondioxide_TraPPE"
        outputfilename = "diffusion_TraPPE_co2_plotshit.dat"
        Type = 1                    # What atom type to follow

    if (calculation == "EMP2"):
        last_timestep = 250000
        startingtime = 200000
        Nsteps = last_timestep - startingtime
        pathtodatabase = '/home/goran/lammps-28Jun14/theproject/carbondioxide/emp2/run_2015_02_27'         # bulk carbon dioxide, EMP2 model
        initialfile = "dump.bulk_carbondioxide_EMPtwo.0.txt"
        frontname = "dump.bulk_carbondioxide_EMPtwo"
        arg = "bulk_carbondioxide_EMPtwo"
        outputfilename = "diffusion_EMP2_co2_plotshit.dat"
        Type = 1                    # What atom type to follow
    
    if (calculation == "FPF"):
        last_timestep = 250000
        startingtime = 200000
        Nsteps = last_timestep - startingtime
        pathtodatabase = '/home/goran/lammps-28Jun14/theproject/carbondioxide/FPF/run_2015_02_26'         # bulk carbon dioxide, FPF model
        initialfile = 'dump.bulk_carbondioxide_FPF.0.txt'
        frontname = "dump.bulk_carbondioxide_FPF"
        arg = "bulk_carbondioxide_FPF"
        outputfilename = 'FPF_co2_diffusion_plotshit.dat'
        outputfilename = "diffusion_FPF_co2_plotshit.dat"
        Type = 1                    # What atom type to follow
        
    if (calculation == "Portlandite"):
        pathtodatabase = '/home/goran/lammps-28Jun14/examples/Abel_runs/portlandite/oldrun'
        initialfile = "dump.portlandite.250000.txt"
        frontname = "dump.portlandite"
        arg = "portlandite"
        outputfilename = "diffusion_Portlandite_co2_plotshit.dat"
        Type = 1                    # What atom type to follow
        
    initialfile = os.path.join(pathtodatabase,initialfile)
    
    
    notsortedfiles, Time = gothroughfiles(pathtodatabase,arg)  # the files are not well enough sorted!
    t0, Natoms, system_size, matrix0, readstructure, entries, types = readfile(initialfile)

    ofile = open(outputfilename,'w')

    #print len(Time), Time[-1]
    #print Time
    #Time = Time[0:1000]
    #print "system size:\n", system_size
    time = []
    start = 0; k = 0;
    for t in Time:
        if t == startingtime:
            start = k
        time.append(t*timestep/1000.0) # convert to pico-seconds
        k += 1
    #######################################################
    #               DENSITY
    
    oxygen = 15.9994                         # g/mol
    hydrogen = 1.00794                       # g/mol
    carbon = 12.0107                         # g/mol
    calcium = 40.078                         # g/mol
    silicon = 28.0855                        # g/mol
    
    if ('water' in initialfile):
        mass = (oxygen + 2*hydrogen)*Natoms/3.0   # g/mol (system)
        Type = 1 # follow oxygen
        system = 'water'
    elif ('carbondioxide' in initialfile):
        mass = (carbon + 2*oxygen)*Natoms/3.0     # g/mol (system)
        Type = 1 # follow carbon
        system = 'water'
    elif ('portlandite' in initialfile):
        mass = (calcium + (oxygen + hydrogen)*2)*Natoms/5.0  # g/mol (system)
        system = 'portlandite'
    elif ('quartz' in initialfile):
        mass = (silicon + 2*oxygen)*Natoms/3.0   # g/mol (system)
        system = 'alpha-quartz'
    else:
        print "failed to characterize system"
    
    #==========================================================================
    #=============== MSD for system, one estimate
    #displacement(ofile,timestep,Natoms,matrix0,Type,Time[start:],frontname,pathtodatabase,system_size)
    
    #==========================================================================
    #=============== MSD for system, len(time)/Nstatefiles estimates
    lTime = len(Time)
    stop = lTime - 1
    stop = start + 500
    print "len(Time) = %g " % (len(Time))
    if (stop > lTime):
        print "Error! stop > len(Time)\n setting stop=len(Time)"
        stop = lTime-1    
    print "start time index = %g" % start
    print "stop time index  = %g" % stop
    GM = storefileinfo(Time,Natoms,frontname,pathtodatabase)                # store fileinfo to a grand matrix
    diffusion(ofile,timestep,Natoms,Nstatefiles,GM,Type,time,system_size)
    
    #==========================================================================
    #====================== Radial Distribution ===============================
    #Carbon Dioxide:
    #Type1 = [1,'C',216] # molID,type,Nmolec
    #Type2 = [2,'O',432] # molID,type,Nmolec

    # Water:
    #Type1 = [1,'O',216] # molID,type,Nmolec
    #Type2 = [2,'H',432] # molID,type,Nmolec

    # Portlandite:
    #Type1 = [1,'Ca',1344]
    #Type2 = [2,'O',2688]

    # Alpha-quartz:
    #Type1 = [1,'O',72] # molID,type,Natoms
    #Type2 = [2,'Si',36] # molID,type,Natoms
    #RDF(ofile,Time[start:],Natoms,system_size,GM,Type1,Type1,pathtodatabase)    # RDF type1 - type1
    #RDF(ofile,Time[start:],Natoms,system_size,GM,Type2,Type2,pathtodatabase)    # RDF type2 - type2
    #RDF(ofile,Time[start:],Natoms,system_size,GM,Type1,Type2,pathtodatabase)    # RDF type1 - type2

    #==========================================================================  
    #======================= Bulk Density =====================================
    #notsort, Time = gothroughfiles(pathtodatabase,arg)
    #bulkdensity(ofile,timestep,mass,Time[start:],frontname,pathtodatabase)
    print "##############################################################"
    print "# Starting timestep:             %g" % Time[start]
    print "# Last timestep:                 %g" % Time[stop]
    print "# Number of files:               %g" % len(time)
    print "# Number of files in estimates:  %g" % Nstatefiles
    print "# Number of estimates:           %g" % (int(np.floor(len(time)/float(Nstatefiles))))
    print "##############################################################"

    ofile.close()
    ######################## END OF MAIN ###########################

if __name__ == "__main__":
    import os
    import numpy as np
    from displacement_ydir import readfile, gothroughfiles
    import matplotlib.pyplot as plt
    from statistics_from_dumpfiles import write_to_file

    """
    font = {'family' : 'sans-serif',
        'weight' : 'bold',
        'size'   : 14}
    plt.rc('font', **font)
    """
    main()