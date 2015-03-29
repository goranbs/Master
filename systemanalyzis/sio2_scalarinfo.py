# -*- coding: utf-8 -*-
"""
Created on Wed Mar  4 16:09:59 2015

@author: goran
"""

def plot_warmup(to_pico_sec): 
    path1 = "/home/goran/lammps-28Jun14/theproject/alpha-quartz/reaxff/Alexandru_sio2_files/slab/run3_warmup_0K_300K_2015_03_05" # slab
    filename1 = "Adump.SiO2_reax_warmup.txt" # slab
    
    fullpath1 = os.path.join(path1,filename1)

    ofile = open(fullpath1,'r')
    ofile.readline(); ofile.readline()

    
    time1 = []; time2 = []; time3 = []   # timestep
    temp1 = []; temp2 = []; temp3 = []   # temperature of system
    Ep1 = []; Ep2 = []; Ep3 = []         # potential energy
    As1 = []; As2 = []; As3 = []         # surface area
    Lz1 = []; Lz2 = []; Lz3 = []         # hight of system
    N1 = []; N2 = []; N3 = []            # number of particles
    Rc1 = []; Rc2 = []; Rc3 = []         # cutoff radius
    
    for line in ofile:
        t,K,Ep,As,Lz,N,Rc = line.split()
        time1.append(int(t)); temp1.append(float(K)); Ep1.append(float(Ep)); As1.append(float(As)); Lz1.append(float(Lz)); N1.append(float(N)); Rc1.append(float(Rc))
    
    ############################################################################
    
    Time1 = []
    for i in range(len(time1)):
        Time1.append(time1[i]*to_pico_sec)
    
    plt.figure()
    Title = "Temperature"
    Xlabel = "time [ps]"
    Ylabel = "temperature [K]"
    Legends = ['$T_{bulk}$','$T_{bulk-onecell}$','$T_{slab}$']
    rc = ['r-']
    plt.plot(Time1,temp1,rc[0],markevery=30)
    plt.title(Title); plt.xlabel(Xlabel); plt.ylabel(Ylabel); plt.legend(Legends)
    
    plt.figure()
    Title = "Potential energy"
    Xlabel = "time [ps]"
    Ylabel = "energy [kcal/mol]"
    Legends = ['$Ep_{bulk}$','$Ep_{bulk-onecell}$','$Ep_{slab}$']
    rc = ['r-']
    plt.plot(Time1,Ep1,rc[0],markevery=30)
    plt.title(Title); plt.xlabel(Xlabel); plt.ylabel(Ylabel); plt.legend(Legends)
    
    plt.show()
    
def main():

    #########################################################################################
    #                        Parameters to set before execution    
    path1 = "/home/goran/lammps-28Jun14/theproject/alpha-quartz/reaxff/Alexandru_sio2_files/bulk/run2_0K_2015_03_04"                            # bulk
    path2 = "/home/goran/lammps-28Jun14/theproject/alpha-quartz/reaxff/Alexandru_sio2_files/slab/run2_0K_2015_03_04"                            # slab
    path3 = "/home/goran/lammps-28Jun14/theproject/alpha-quartz/reaxff/Alexandru_sio2_files/bulk_single_cell/run1_2015_03-05"                   #bulk one-cell
    path4 = "/home/goran/lammps-28Jun14/theproject/alpha-quartz/reaxff/Alexandru_sio2_files/slab/run4_nonperiodic_z"                            # slab nonperiodic Z
    path5 = "/home/goran/lammps-28Jun14/theproject/alpha-quartz/reaxff/Alexandru_sio2_files/bulk/run3_0K_displace_atoms_2015_03_10"             # bulk displaced
    path6 = "/home/goran/lammps-28Jun14/theproject/alpha-quartz/reaxff/Alexandru_sio2_files/bulk_single_cell/run2_0K_displace_atoms_2015_03_10" # bulk one-cell displaced

    path_single_cell = "/home/goran/lammps-28Jun14/theproject/alpha-quartz/reaxff/Alexandru_sio2_files/bulk_single_cell/"
    path_slab = "/home/goran/lammps-28Jun14/theproject/alpha-quartz/reaxff/Alexandru_sio2_files/slab"
    path_slab_nonp = "/home/goran/lammps-28Jun14/theproject/alpha-quartz/reaxff/Alexandru_sio2_files/slab_nonperiodic"

    
    path7 = os.path.join(path_single_cell,"run3_warmup_0-300K")         # Bulk One-cell
    path8 = os.path.join(path_single_cell,"run4_warmup_300-500K")
    path9 = os.path.join(path_single_cell,"run5_warmup_500-1000K")
    path10 = os.path.join(path_single_cell,"run6_warmup_1000-1500K")
    
    path11 = os.path.join(path_slab,"run3_warmup_0-300K")               # 27-layer slab
    path12 = os.path.join(path_slab,"run5_warmup_300-500K")
    path13 = os.path.join(path_slab,"run6_warmup_500-1000K")
    path14 = os.path.join(path_slab,"run7_warmup_1000-1500K")
    
    path15 = os.path.join(path_slab_nonp, "run1_0-300K")                # 27-layer slab non-periodic in z-dir
    path16 = os.path.join(path_slab_nonp, "run2_300-500K")
    path17 = os.path.join(path_slab_nonp, "run3_500-1000K")
    path18 = os.path.join(path_slab_nonp, "run4_1000-1500K")

    path19 = "/home/goran/lammps-28Jun14/theproject/alpha-quartz/reaxff/Alexandru_sio2_files/cleaved/run1_100-0K" # cleaved surface   
    
    path20 = "/home/goran/lammps-28Jun14/theproject/alpha-quartz/reaxff/Alexandru_sio2_files/carbondioxide_SiO2/pure_CO2/run1_0K_nvt_nonperiodic"                        # pure CO2             non-periodic
    path21 = "/home/goran/lammps-28Jun14/theproject/alpha-quartz/reaxff/Alexandru_sio2_files/carbondioxide_SiO2/CO3_1_f/run1_cooling_0K_nvt"     # CO2 on SiO2 surface  non-periodic
    path22 = "/home/goran/lammps-28Jun14/theproject/alpha-quartz/reaxff/Alexandru_sio2_files/carbondioxide_SiO2/CO3_2_f/run1_cooling_0K_nvt"     # CO2 on SiO2 surface  non-periodic
    
    path23 = "/home/goran/lammps-28Jun14/theproject/alpha-quartz/reaxff/Alexandru_sio2_files/carbondioxide_SiO2/CO3_1/run1_0K_nvt_periodic"    
    path24 = "/home/goran/lammps-28Jun14/theproject/alpha-quartz/reaxff/Alexandru_sio2_files/carbondioxide_SiO2/CO3_1_f/run2_0K_nvt_periodic"
    path25 = "/home/goran/lammps-28Jun14/theproject/alpha-quartz/reaxff/Alexandru_sio2_files/carbondioxide_SiO2/CO3_2_f/run2_0K_nvt_periodic"
    path26 = "/home/goran/lammps-28Jun14/theproject/alpha-quartz/reaxff/Alexandru_sio2_files/carbondioxide_SiO2/SiO4_C_O/run1_0K_nvt_periodic"
    path27 = "/home/goran/lammps-28Jun14/theproject/alpha-quartz/reaxff/Alexandru_sio2_files/carbondioxide_SiO2/pure_CO2/run2_0K_nvt_periodic"
    #path28 = "/home/goran/lammps-28Jun14/theproject/alpha-quartz/reaxff/Alexandru_sio2_files/carbondioxide_SiO2/pure_CO2/run2_0K_nvt_nonperiodic"
    
    path28 = "/home/goran/lammps-28Jun14/theproject/alpha-quartz/reaxff/Alexandru_sio2_files/carbondioxide_SiO2/CO3_1/minimize"
    path29 = "/home/goran/lammps-28Jun14/theproject/alpha-quartz/reaxff/Alexandru_sio2_files/carbondioxide_SiO2/CO3_1_f/minimize"
    path30 = "/home/goran/lammps-28Jun14/theproject/alpha-quartz/reaxff/Alexandru_sio2_files/carbondioxide_SiO2/CO3_2_f/minimize"
    path31 = "/home/goran/lammps-28Jun14/theproject/alpha-quartz/reaxff/Alexandru_sio2_files/carbondioxide_SiO2/SiO4_C_O/minimize"

    path32 = "/home/goran/lammps-28Jun14/theproject/alpha-quartz/reaxff/Alexandru_sio2_files/slab/minimize"
    
    
    filename_bulk_warmup = "Adump.SiO2_reax_warmup.txt"
    filename_slab_warmup = "Adump.SiO2_reax_warmup.txt"
        
        
    filename1 = "Adump.SiO2_reax_cooling.txt"               # cooling of the bulk system
    filename2 = "Adump.SiO2_reax_zero_temp.txt"             # bulk system at 0K
    filename3 = "Adump.SiO2_reax_cooling.txt"               # slab system at 0K    
    filename4 = "Adump.SiO2_reax_at_zero_temp.txt"          # cooling of the slab system
    filename5 = "Adump.SiO2_reax_0K_highprecision.txt"      # bulk system at 0K. High precision of potential energy! remember to divide by 100.000
    filename6 = "Adump.SiO2_reax_0K_highprecision.txt"      # slab system at 0K. High precision of potential energy! multiplied by 100.000
    filename7 = "Adump.SiO2_reax_cooldown.txt"              # bulk one-cell system. cooling of the system
    filename8 = "Adump.SiO2_reax_restart_multiplied.txt"    # bulk one-cell system. at 0K. potential energy multiplied by factor 100.000
    filename9 = "Adump.SiO2_reax_nonperiodicZ_cooling.txt"  # Slab nonperiodicZ. cooling of the system
    filename10 = "Adump.SiO2_reax_0K.txt"                   # Slab nonperiodicZ. system at 0K
    filename11 = "Adump.SiO2_reax_displaced.txt"            # Slab 40000 timesteps extra-run
    filename12 = "Adump.SiO2_reax_displaced.txt"            # Slab 40000 timesteps extra-run
    filename13 = "Adump.SiO2_reax_displaced.txt"            # Bulk displaced
    filename14 = "Adump.SiO2_reax_displaced.txt"            # Bulk one-cell displaced
    filename15 = "Adump.carbondioxide_reax.txt"             # CO2      non-periodic
    filename16 = "Adump.SiO2_CO2_reax.txt"                  # CO3_1_f  non-periodic
    filename17 = "Adump.SiO2_CO2_reax.txt"                  # CO3_2_f  non-periodic

    
    somefactor = 100000   # number to divide the potential energy in filename5 and fileneme6 by!
    #Nsteps1 = 300000      # number of timesteps in simulation 1
    #Nsteps2 = 30000       # number of timesteps in simulation 2
    #Nsteps3 = 100000      # number of timesteps in simulation 3
    dt = 0.1              # fsec
    
    #########################################################################################
        
    to_pico_sec = dt/1000   # conversionfactor from timestep to ps
    #plot_warmup(to_pico_sec)
    
    fullpath1 = os.path.join(path1,filename1)    # bulk
    fullpath2 = os.path.join(path1,filename2)    # bulk at 0K
    fullpath3 = os.path.join(path2,filename3)    # slab
    fullpath4 = os.path.join(path2,filename4)    # slab at 0K
    fullpath5 = os.path.join(path1,filename5)    # bulk at 0K highprecision
    fullpath6 = os.path.join(path2,filename6)    # slab at 0K highprecision
    fullpath7 = os.path.join(path3,filename7)    # bulk one-cell at 0K highprecision
    fullpath8 = os.path.join(path3,filename8)    # bulk one-cell at 0K highprecision
    fullpath9 = os.path.join(path4,filename9)    # bulk one-cell at 0K highprecision
    fullpath10 = os.path.join(path4,filename10)  # slab non-periodic in z dir. 
    fullpath11 = os.path.join(path2,filename11)  # slab 
    fullpath12 = os.path.join(path4,filename12)  # slab non-periodic
    fullpath13 = os.path.join(path5,filename13)  # bulk 
    fullpath14 = os.path.join(path6,filename14)  # bulk one-cell
    
    fullpath15 = os.path.join(path7,"Adump.SiO2_reax_300K.txt")   # fullpath of bulk onecell 300K run
    fullpath16 = os.path.join(path8,"Adump.SiO2_reax_500K.txt")   # fullpath of bulk onecell 500K run
    fullpath17 = os.path.join(path9,"Adump.SiO2_reax_1000K.txt")  # fullpath of bulk onecell 1000K run
    fullpath18 = os.path.join(path10,"Adump.SiO2_reax_1500K.txt") # fullpath of bulk onecell 1500K run
    
    fullpath19 = os.path.join(path11,"Adump.SiO2_300K.txt")       # fullpath of slab 300K run
    fullpath20 = os.path.join(path12,"Adump.SiO2_500K.txt")       # fullpath of slab 500K run
    fullpath21 = os.path.join(path13,"Adump.SiO2_1000K.txt")      # fullpath of slab 1000K run
    fullpath22 = os.path.join(path14,"Adump.SiO2_1500K.txt")      # fullpath of slab 1500K run
    
    fullpath23 = os.path.join(path15,"Adump.SiO2_reax_300K.txt")  # fullpath of slab nonperioodic Z 300K run
    fullpath24 = os.path.join(path16,"Adump.SiO2_reax_500K.txt")  # fullpath of slab nonperioodic Z 500K run
    fullpath25 = os.path.join(path17,"Adump.SiO2_reax_1000K.txt") # fullpath of slab nonperioodic Z 1000K run
    fullpath26 = os.path.join(path18,"Adump.SiO2_reax_1500K.txt") # fullpath of slab nonperioodic Z 1500K run
        
    fullpath27 = os.path.join(path19,"Adump.SiO2_cleaved.txt")    # cleaved surface
    fullpath28 = os.path.join(path19,"Adump.SiO2_0K.txt")    # cleaved surface
    
    fullpath29 = os.path.join(path20,filename15) # CO2 system
    fullpath30 = os.path.join(path21,filename16) # CO3_1_f
    fullpath31 = os.path.join(path22,filename17) # CO3_2_f
    
    fullpath32 = os.path.join(path23,"Adump.SiO2_periodic.txt")             # CO3_1
    fullpath33 = os.path.join(path24,"Adump.SiO2_CO2_periodic.txt")         # CO3_1_f
    fullpath34 = os.path.join(path25,"Adump.SiO2_CO2_periodic.txt")         # CO3_2_f
    fullpath35 = os.path.join(path26,"Adump.SiO2_periodic.txt")             # SiO4_C_O   
    fullpath36 = os.path.join(path27,"Adump.carbondioxide_periodic.txt")    # CO2 molecule. periodic system
    #fullpath37 = os.path.join(path28,"Adump.carbondioxide_single_nonperiodic")  # CO2 molecule. nonperiodic system
    
    fullpath37 = os.path.join(path28,"Adump.quartz_n_carbondioxide.txt")   # CO3_1    minimization
    fullpath38 = os.path.join(path29,"Adump.quartz_n_carbondioxide.txt")   # CO3_1_f  minimization
    fullpath39 = os.path.join(path30,"Adump.quartz_n_carbondioxide.txt")   # CO3_2_f  minimization
    fullpath40 = os.path.join(path31,"Adump.quartz_n_carbondioxide.txt")   # SiO4_C_O minimization
    
    fullpath41 = os.path.join(path32,"Adump.quartz_n_carbs.txt")  # minimized energy of slab. periodic
    
    ofile1 = open(fullpath1,'r')
    ofile2 = open(fullpath2,'r')
    ofile3 = open(fullpath3,'r')
    ofile4 = open(fullpath4,'r')
    ofile5 = open(fullpath5,'r')
    ofile6 = open(fullpath6,'r')
    ofile7 = open(fullpath7,'r')
    ofile8 = open(fullpath8,'r')
    ofile9 = open(fullpath9,'r')
    ofile10 = open(fullpath10,'r')
    ofile11 = open(fullpath11,'r')  # slab
    ofile12 = open(fullpath12,'r')  # slab non-periodic
    ofile13 = open(fullpath13,'r')  # bulk
    ofile14 = open(fullpath14,'r')  # bulk one-cell
    
    ofile15 = open(fullpath15,'r')  # bulk one-cell warmup from 0-300K
    ofile16 = open(fullpath16,'r') # bulk one-cell warmup from 300-500K
    ofile17 = open(fullpath17,'r') # bulk one-cell warmup from 500-1000K
    ofile18 = open(fullpath18,'r') # bulk one-cell warmup from 1000-1500K
    
    ofile19 = open(fullpath19,'r')  # slab warmup from 0-300K
    ofile20 = open(fullpath20,'r') # slab warmup from 300-500K
    ofile21 = open(fullpath21,'r') # slab warmup from 500-1000K
    ofile22 = open(fullpath22,'r') # slab warmup from 1000-1500K
    
    ofile23 = open(fullpath23,'r')  # slab nonperiodic Z warmup from 0-300K
    ofile24 = open(fullpath24,'r')  # slab nonperiodic Z warmup from 300-500K
    ofile25 = open(fullpath25,'r')  # slab nonperiodic Z warmup from 500-1000K
    ofile26 = open(fullpath26,'r')  # slab nonperiodic Z warmup from 1000-1500K
    
    ofile27 = open(fullpath27,'r')  # cleaved surface
    ofile28 = open(fullpath28,'r')  # cleaved surface part 2
    
    ofile29 = open(fullpath29,'r')  # CO2
    ofile30 = open(fullpath30,'r')  # CO3_1_f
    ofile31 = open(fullpath31,'r')  # CO3_2_f
    
    ofile32 = open(fullpath32,'r')  # CO3_1
    ofile33 = open(fullpath33,'r')  # CO3_1_f
    ofile34 = open(fullpath34,'r')  # CO3_2_f
    ofile35 = open(fullpath35,'r')  # SiO4_C_O
    ofile36 = open(fullpath36,'r')  # CO2 molecule periodic system
    #ofile37 = open(fullpath37,'r')  # CO2 molecule nonperiodic system
    
    ofile37 = open(fullpath37,'r')  # CO3_1     minimized
    ofile38 = open(fullpath38,'r')  # CO3_1_f   minimized
    ofile39 = open(fullpath39,'r')  # CO3_2_f   minimized
    ofile40 = open(fullpath40,'r')  # SiO4_C_O  minimized
    
    ofile41  = open(fullpath41,'r')  # minimized system. dense27 layer slab. periodic
    
    
    bulk_warmups = [ofile15,ofile16,ofile17,ofile18]
    slab_warmups = [ofile19,ofile20,ofile21]#,ofile22]
    slab_nonp_warmups = [ofile23,ofile24,ofile25,ofile26]
    
    # jump two fires lines:
    ofile1.readline(); ofile2.readline(); ofile3.readline(); ofile4.readline(); ofile5.readline(); ofile6.readline(); ofile7.readline(); ofile8.readline(); ofile9.readline() # skip first line
    ofile1.readline(); ofile2.readline(); ofile3.readline(); ofile4.readline(); ofile5.readline(); ofile6.readline(); ofile7.readline(); ofile8.readline(); ofile9.readline() # skip second line!
    ofile10.readline(); ofile11.readline(); ofile12.readline(); ofile13.readline(); ofile14.readline()
    ofile10.readline(); ofile11.readline();  ofile12.readline(); ofile13.readline(); ofile14.readline()
    ofile15.readline(); ofile16.readline(); ofile17.readline(); ofile18.readline() 
    ofile15.readline(); ofile16.readline(); ofile17.readline(); ofile18.readline()
    ofile19.readline(); ofile20.readline(); ofile21.readline(); ofile22.readline()
    ofile19.readline(); ofile20.readline(); ofile21.readline(); ofile22.readline()
    ofile23.readline(); ofile24.readline(); ofile25.readline(); ofile26.readline()
    ofile23.readline(); ofile24.readline(); ofile25.readline(); ofile26.readline()
    
    ofile27.readline(); ofile27.readline(); ofile28.readline(); ofile28.readline()
    ofile29.readline(); ofile30.readline(); ofile31.readline()
    ofile29.readline(); ofile30.readline(); ofile31.readline()
    ofile32.readline(); ofile33.readline(); ofile34.readline(); ofile35.readline(); ofile36.readline()#; ofile37.readline()
    ofile32.readline(); ofile33.readline(); ofile34.readline(); ofile35.readline(); ofile36.readline()#; ofile37.readline()

    ofile37.readline(); ofile38.readline(); ofile39.readline(); ofile40.readline()
    ofile37.readline(); ofile38.readline(); ofile39.readline(); ofile40.readline()
    
    ofile41.readline(); ofile41.readline()
    
    ###########################################################################
    # 1) bulk. e.g  all lists with end character 1 is for bulk system
    # 2) slab system
    # 3) bulk one-cell system
    # 4) slab system that is nonperiodic in z direction
    
    time1 = []; time2 = []; time3 = []; time4 = []   # timestep
    temp1 = []; temp2 = []; temp3 = []; temp4 = []   # temperature of system
    Ep1 = []; Ep2 = []; Ep3 = []; Ep4 = []           # potential energy
    As1 = []; As2 = []; As3 = []; As4 = []           # surface area
    Lz1 = []; Lz2 = []; Lz3 = []; Lz4 = []           # hight of system
    N1 = []; N2 = []; N3 = []; N4 = []               # number of particles
    Rc1 = []; Rc2 = []; Rc3 = []; Rc4 = []           # cutoff radius
    
    time5 = []; temp5 = []; Ep5 = []; As5 = []; Lz5 = []; N5 = []; Rc5 = []  # cleaved surface    
    time6 = []; temp6 = []; Ep6 = []; As6 = []; Lz6 = []; N6 = []; Rc6 = []  # CO2    
    time7 = []; temp7 = []; Ep7 = []; As7 = []; Lz7 = []; N7 = []; Rc7 = []  # CO3_1_f
    time8 = []; temp8 = []; Ep8 = []; As8 = []; Lz8 = []; N8 = []; Rc8 = []  # CO3_2_f
    time9 = []; temp9 = []; Ep9 = []; As9 = []; Lz9 = []; N9 = []; Rc9 = []         # CO3_1    periodic
    time10 = []; temp10 = []; Ep10 = []; As10 = []; Lz10 = []; N10 = []; Rc10 = []  # CO3_1_f  periodic
    time11 = []; temp11 = []; Ep11 = []; As11 = []; Lz11 = []; N11 = []; Rc11 = []  # CO3_2_f  periodic
    time12 = []; temp12 = []; Ep12 = []; As12 = []; Lz12 = []; N12 = []; Rc12 = []  # SiO4_C_O periodic
    time13 = []; temp13 = []; Ep13 = []; As13 = []; Lz13 = []; N13 = []; Rc13 = []  # CO2
    #time14 = []; temp14 = []; Ep14 = []; As14 = []; Lz14 = []; N14 = []; Rc14 = []  # SiO4_C_O periodic
    
    for line in ofile27:
        t,K,Ep,As,Lz,N,Rc = line.split()
        time5.append(int(t)); temp5.append(float(K)); Ep5.append(float(Ep)/100.0); As5.append(float(As)); Lz5.append(float(Lz)); N5.append(float(N)); Rc5.append(float(Rc))
    for line in ofile28:
        t,K,Ep,As,Lz,N,Rc = line.split()
        time5.append(int(t)); temp5.append(float(K)); Ep5.append(float(Ep)); As5.append(float(As)); Lz5.append(float(Lz)); N5.append(float(N)); Rc5.append(float(Rc))
                        
    for line in ofile1:
        t,K,Ep,As,Lz,N,Rc = line.split()
        time1.append(int(t)); temp1.append(float(K)); Ep1.append(float(Ep)); As1.append(float(As)); Lz1.append(float(Lz)); N1.append(float(N)); Rc1.append(float(Rc))
    for line in ofile2:
        t,K,Ep,As,Lz,N,Rc = line.split()
        time1.append(int(t)); temp1.append(float(K)); Ep1.append(float(Ep)); As1.append(float(As)); Lz1.append(float(Lz)); N1.append(float(N)); Rc1.append(float(Rc))
    for line in ofile3:
        t,K,Ep,As,Lz,N,Rc = line.split()
        time2.append(int(t)); temp2.append(float(K)); Ep2.append(float(Ep)); As2.append(float(As)); Lz2.append(float(Lz)); N2.append(float(N)); Rc2.append(float(Rc))
    for line in ofile4:
        t,K,Ep,As,Lz,N,Rc = line.split()
        time2.append(int(t)); temp2.append(float(K)); Ep2.append(float(Ep)); As2.append(float(As)); Lz2.append(float(Lz)); N2.append(float(N)); Rc2.append(float(Rc))
    for line in ofile5:
        t,K,Ep,As,Lz,N,Rc = line.split()
        time1.append(int(t)); temp1.append(float(K)); Ep1.append(float(Ep)/somefactor); As1.append(float(As)); Lz1.append(float(Lz)); N1.append(float(N)); Rc1.append(float(Rc))
    for line in ofile6:
        t,K,Ep,As,Lz,N,Rc = line.split()
        time2.append(int(t)); temp2.append(float(K)); Ep2.append(float(Ep)/somefactor); As2.append(float(As)); Lz2.append(float(Lz)); N2.append(float(N)); Rc2.append(float(Rc))
    for line in ofile7:
        t,K,Ep,As,Lz,N,Rc = line.split()
        time3.append(int(t)); temp3.append(float(K)); Ep3.append(float(Ep)); As3.append(float(As)); Lz3.append(float(Lz)); N3.append(float(N)); Rc3.append(float(Rc))
    for line in ofile8:
        t,K,Ep,As,Lz,N,Rc = line.split()
        time3.append(int(t)); temp3.append(float(K)); Ep3.append(float(Ep)/somefactor); As3.append(float(As)); Lz3.append(float(Lz)); N3.append(float(N)); Rc3.append(float(Rc))
    for line in ofile9:
        t,K,Ep,As,Lz,N,Rc = line.split()
        time4.append(int(t)); temp4.append(float(K)); Ep4.append(float(Ep)); As4.append(float(As)); Lz4.append(float(Lz)); N4.append(float(N)); Rc4.append(float(Rc))
    for line in ofile10:
        t,K,Ep,As,Lz,N,Rc = line.split()
        time4.append(int(t)); temp4.append(float(K)); Ep4.append(float(Ep)/100); As4.append(float(As)); Lz4.append(float(Lz)); N4.append(float(N)); Rc4.append(float(Rc))
    for line in ofile11: # slab displace
        t,K,Ep,As,Lz,N,Rc = line.split()
        time2.append(int(t)); temp2.append(float(K)); Ep2.append(float(Ep)/100); As2.append(float(As)); Lz2.append(float(Lz)); N2.append(float(N)); Rc2.append(float(Rc))
    for line in ofile29:  # CO2
        t,K,Ep,As,Lz,N,Rc = line.split()
        time6.append(int(t)); temp6.append(float(K)); Ep6.append(float(Ep)/100); As6.append(float(As)); Lz6.append(float(Lz)); N6.append(float(N)); Rc6.append(float(Rc))
    for line in ofile30:  # CO3_1_f
        t,K,Ep,As,Lz,N,Rc = line.split()
        time7.append(int(t)); temp7.append(float(K)); Ep7.append(float(Ep)/100); As7.append(float(As)); Lz7.append(float(Lz)); N7.append(float(N)); Rc7.append(float(Rc))
    for line in ofile31: # CO3_2_f
        t,K,Ep,As,Lz,N,Rc = line.split()
        time8.append(int(t)); temp8.append(float(K)); Ep8.append(float(Ep)/100); As8.append(float(As)); Lz8.append(float(Lz)); N8.append(float(N)); Rc8.append(float(Rc))
    
    for line in ofile32: # CO3_1     periodic
        t,K,Ep,As,Lz,N,Rc = line.split()
        time9.append(int(t)); temp9.append(float(K)); Ep9.append(float(Ep)/100); As9.append(float(As)); Lz9.append(float(Lz)); N9.append(float(N)); Rc9.append(float(Rc))

    for line in ofile33: # CO3_1_f   periodic
        t,K,Ep,As,Lz,N,Rc = line.split()
        time10.append(int(t)); temp10.append(float(K)); Ep10.append(float(Ep)/100); As10.append(float(As)); Lz10.append(float(Lz)); N10.append(float(N)); Rc10.append(float(Rc))
        
    for line in ofile34: # CO3_2_f   periodic
        t,K,Ep,As,Lz,N,Rc = line.split()
        time11.append(int(t)); temp11.append(float(K)); Ep11.append(float(Ep)/100); As11.append(float(As)); Lz11.append(float(Lz)); N11.append(float(N)); Rc11.append(float(Rc))
    
    for line in ofile35: # SiO4_C_O   periodic
        t,K,Ep,As,Lz,N,Rc = line.split()
        time12.append(int(t)); temp12.append(float(K)); Ep12.append(float(Ep)/100); As12.append(float(As)); Lz12.append(float(Lz)); N12.append(float(N)); Rc12.append(float(Rc))
    
    for line in ofile36: # CO2   periodic
        t,K,Ep,As,Lz,N,Rc = line.split()
        time13.append(int(t)); temp13.append(float(K)); Ep13.append(float(Ep)/100); As13.append(float(As)); Lz13.append(float(Lz)); N13.append(float(N)); Rc13.append(float(Rc))
    
    #for line in ofile37: # CO2   non-periodic
    #    t,K,Ep,As,Lz,N,Rc = line.split()
    #    time14.append(int(t)); temp14.append(float(K)); Ep14.append(float(Ep)/100); As14.append(float(As)); Lz14.append(float(Lz)); N14.append(float(N)); Rc14.append(float(Rc))
    
    for line in ofile37:
        t,K,Ep,As,Lz,N,Rc = line.split()
        tt = int(t)
        if (tt <= time9[-1]):
            dt = time9[1] - time9[0]
            tt = time9[-1] + dt
        time9.append(tt); temp9.append(float(K)); Ep9.append(float(Ep)/100.0); As9.append(float(As)); Lz9.append(float(Lz)); N9.append(float(N)); Rc9.append(float(Rc))
    
    for line in ofile38:
        t,K,Ep,As,Lz,N,Rc = line.split()
        tt = int(t)
        if (tt <= time10[-1]):
            dt = time10[1] - time10[0]
            tt = time10[-1] + dt
        time10.append(tt); temp10.append(float(K)); Ep10.append(float(Ep)/100.0); As10.append(float(As)); Lz10.append(float(Lz)); N10.append(float(N)); Rc10.append(float(Rc))
        
    for line in ofile39:
        t,K,Ep,As,Lz,N,Rc = line.split()
        tt = int(t)
        if (tt <= time11[-1]):
            dt = time11[1] - time11[0]
            tt = time11[-1] + dt
        time11.append(tt); temp11.append(float(K)); Ep11.append(float(Ep)/100.0); As11.append(float(As)); Lz11.append(float(Lz)); N11.append(float(N)); Rc11.append(float(Rc))
        
    for line in ofile40:
        t,K,Ep,As,Lz,N,Rc = line.split()
        tt = int(t)
        if (tt <= time12[-1]):
            dt = time12[1] - time12[0]
            tt = time12[-1] + dt
        time12.append(tt); temp12.append(float(K)); Ep12.append(float(Ep)/100.0); As12.append(float(As)); Lz12.append(float(Lz)); N12.append(float(N)); Rc12.append(float(Rc))
    
    
    for line in ofile41:
        t,K,Ep,As,Lz,N,Rc = line.split()
        tt = int(t)
        if (tt <= time2[-1]):
            dt = time2[1] - time2[0]
            tt = time2[-1] + dt
        time2.append(tt); temp2.append(float(K)); Ep2.append(float(Ep)/100.0); As2.append(float(As)); Lz2.append(float(Lz)); N2.append(float(N)); Rc2.append(float(Rc))
    
    '''
    for line in ofile12: # slab non-periodic displace
        t,K,Ep,As,Lz,N,Rc = line.split()
        time4.append(int(t)); temp4.append(float(K)); Ep4.append(float(Ep)/100); As4.append(float(As)); Lz4.append(float(Lz)); N4.append(float(N)); Rc4.append(float(Rc))
    for line in ofile13: # bulk displace
        t,K,Ep,As,Lz,N,Rc = line.split()
        time1.append(int(t)); temp1.append(float(K)); Ep1.append(float(Ep)/100); As1.append(float(As)); Lz1.append(float(Lz)); N1.append(float(N)); Rc1.append(float(Rc))
    for line in ofile14: # bulk one-cell displace
        t,K,Ep,As,Lz,N,Rc = line.split()
        time3.append(int(t)); temp3.append(float(K)); Ep3.append(float(Ep)/100); As3.append(float(As)); Lz3.append(float(Lz)); N3.append(float(N)); Rc3.append(float(Rc))
     '''                              

    time_wb = []; time_ws = []; time_wsnp = []  # time arrays
    temp_wb = []; temp_ws = []; temp_wsnp = []  # temperature arrays
    Ep_wb = []; Ep_ws = []; Ep_wsnp = []        # potential energy
    As_wb = []; As_ws = []; As_wsnp = []        # surface area
    N_wb = []; N_ws = []; N_wsnp = []           # Number of atoms
    Rc_wb = []; Rc_ws = []; Rc_wsnp = []          # cutoff radius in simulation.
    Lz_wb = []; Lz_ws = []; Lz_wsnp = []        # hight of system. z-dir
    
    for afile in bulk_warmups:
        for line in afile:
            try:
                t,K,Ep,As,Lz,N,Rc = line.split()
                time_wb.append(int(t)); temp_wb.append(float(K)); Ep_wb.append(float(Ep)/100); As_wb.append(float(As)); Lz_wb.append(float(Lz)); N_wb.append(float(N)); Rc_wb.append(float(Rc))
            except:
                print line
                print "Error, not able to unpack datas into: t,K,Ep,As,Lz,N,Rc"

    for afile in slab_warmups:
        # Read only every 10th line! 
        i = 0
        for line in afile:
            i += 1
            if (i == 10):
                t,K,Ep,As,Lz,N,Rc = line.split()
                #print t
                time_ws.append(int(t)); temp_ws.append(float(K)); Ep_ws.append(float(Ep)/100); As_ws.append(float(As)); Lz_ws.append(float(Lz)); N_ws.append(float(N)); Rc_ws.append(float(Rc))
                i = 0   
    for line in ofile22:
        t,K,Ep,As,Lz,N,Rc = line.split()
        time_ws.append(int(t)); temp_ws.append(float(K)); Ep_ws.append(float(Ep)/100); As_ws.append(float(As)); Lz_ws.append(float(Lz)); N_ws.append(float(N)); Rc_ws.append(float(Rc))
        
    
    ltime = len(time_wb)
    i = 0
    for afile in slab_nonp_warmups:
        for line in afile:
            t,K,Ep,As,Lz,N,Rc = line.split()
            #print "(%g/%g) times: t_np= %d  t_p= %d  t_b= %d" % (i,ltime,int(t),time_ws[i],time_wb[i])
            time_wsnp.append(int(t)); temp_wsnp.append(float(K)); Ep_wsnp.append(float(Ep)/100); As_wsnp.append(float(As)); Lz_wsnp.append(float(Lz)); N_wsnp.append(float(N)); Rc_wsnp.append(float(Rc))
            i+=1
    ###########################################################################################
    #                    Time for som serious PLOTTING !!!!!!!!!!!!!!!!!!!!!
    
    #print len(Ep_wb), len(N_wb), len(N_ws), len(N_wsnp)  
    '''    
    print len(time9), len(time10), len(time11), len(time12)
    print time9
    print "---------------------------------\n"
    print time10
    print "---------------------------------\n"
    print time11
    print "---------------------------------\n"
    print time12
    print "---------------------------------\n"    
    '''
    Time1 = []; Time2 = []; Time3 = []; Time4 = []; Epp = []
    Time_wb = []; Time_ws = []; Epp_wb = []; Time_wsnp = []
    Time5 = []; Time6 = [];  Time7 = []; Time8 = []; Time9 = []; Time10 = []; Time11 = []; Time12 = []; Time13 = []; Time14 = []
    for i in range(len(time1)):
        Time1.append(time1[i]*to_pico_sec)  # bulk                           1
    for i in range(len(time2)):
        Time2.append(time2[i]*to_pico_sec)  # slab                           2
    for i in range(len(time3)): 
        Time3.append(time3[i]*to_pico_sec)  # bulk one-cell                  3
        Epp.append(Ep3[i]*(N1[i]/N3[i]))
    for i in range(len(time4)):
        Time4.append(time4[i]*to_pico_sec)  # slab non-periodic              4
    for i in range(len(time_wb)):
        Time_wb.append(time_wb[i]*to_pico_sec)      # bulk warmup system     wb
        Epp_wb.append(Ep_wb[i]*(N_ws[i]/N_wb[i]))
    for i in range(len(time_ws)):
        Time_ws.append(time_ws[i]*to_pico_sec)      # slab warmup system     ws
    for i in range(len(time_wsnp)):
        Time_wsnp.append(time_wsnp[i]*to_pico_sec)  # slab nonp warmup system   wsnp
    for i in range(len(time5)):
        Time5.append(time5[i]*to_pico_sec)      # cleaved surface            5
    for i in range(len(time6)):
        Time6.append(time6[i]*to_pico_sec)      # CO2                        6
    for i in range(len(time7)):
        Time7.append(time7[i]*to_pico_sec)      # CO3_1_f                    7
        Time8.append(time7[i]*to_pico_sec)
    
    for i in range(len(time9)):
        Time9.append(time9[i]*to_pico_sec)        # CO3_1                    9
    for i in range(len(time10)):
        Time10.append(time10[i]*to_pico_sec)      # CO3_1_f                  10
    for i in range(len(time11)):
        Time11.append(time11[i]*to_pico_sec)      # CO3_2_f                  11
    for i in range(len(time12)):
        Time12.append(time12[i]*to_pico_sec)      # SiO4_C_O                 12
    for i in range(len(time13)):
        Time13.append(time13[i]*to_pico_sec)      # CO2 periodic             13
    #for i in range(len(time14)):
    #    Time14.append(time14[i]*to_pico_sec)      # CO2 non-periodic         14
        
            
    
    ############################################################################################
    #=========================== Carbon dioxide on SiO2 surface ===============================#
    
    #print len(time7), len(time8)
    #print len(Ep7), len(Ep8)
    #print len(Time9), len(Time10),len(Time11),len(Time12)
    #rint len(Ep9), len(Ep10),len(Ep11),len(Ep12)

    plt.figure()
    Title = "Potential energy"
    Xlabel = "time [ps]"
    Ylabel = "energy [kcal/mol]"
    Legends = ['$Ep_{CO3_{1f,np}}$','$Ep_{CO3_{2f,np}}$','$Ep_{CO3_{1,p}}$','$Ep_{CO3_{1f,p}}$','$Ep_{CO3_{2f_p}}$','$Ep_{SiO4-C-O_{p}}$']
    rc = ['r--','y--','b-','r-','y-','m-']
    lstart = int(round(len(Time6)*0.86))
    mean_ep_co3_1_f = np.mean(Ep7[lstart:])
    std_ep_co3_1_f = np.std(Ep7[lstart:])
    mean_ep_co3_2_f = np.mean(Ep8[lstart:])
    std_ep_co3_2_f = np.std(Ep8[lstart:])
    mline_co3_1f = np.linspace(mean_ep_co3_1_f,mean_ep_co3_1_f,2)
    mline_co3_2f = np.linspace(mean_ep_co3_2_f,mean_ep_co3_2_f,2)
    
    m_co31 = np.mean(Ep9[lstart:])
    std_co31 = np.std(Ep9[lstart:])
    mlco31 = np.linspace(m_co31,m_co31,2)
    
    m_co31f = np.mean(Ep10[lstart:])
    std_co31f = np.std(Ep10[lstart:])
    mlco31f = np.linspace(m_co31f,m_co31f,2)
    
    m_co32f = np.mean(Ep11[lstart:])
    std_co32f = np.std(Ep11[lstart:])
    mlco32f = np.linspace(m_co32f,m_co32f,2)
    
    m_SiO4_C_O = np.mean(Ep12[lstart:])
    std_SiO4_C_O = np.std(Ep12[lstart:])
    mlSiO4_C_O = np.linspace(m_SiO4_C_O,m_SiO4_C_O,2)   
    
    plt.hold(True)
    plt.plot(Time7,Ep7,rc[0])
    plt.plot(Time8,Ep8,rc[1])
    plt.plot(Time9,Ep9,rc[2])
    plt.plot(Time10,Ep10,rc[3])
    plt.plot(Time11,Ep11,rc[4])
    plt.plot(Time12,Ep12,rc[5])
    #plt.plot(Time13,Ep13,'r--')
    plt.plot([Time7[lstart],Time7[-1]], mline_co3_1f,'k--')
    plt.plot([Time7[lstart],Time8[-1]], mline_co3_2f,'k--')
    plt.plot([Time9[lstart],Time9[-1]], mlco31,'k--')
    plt.plot([Time9[lstart],Time9[-1]], mlco31f,'k--')
    plt.plot([Time9[lstart],Time9[-1]], mlco32f,'k--')
    plt.plot([Time9[lstart],Time9[-1]], mlSiO4_C_O,'k--')
    plt.hold(False)
    plt.title(Title); plt.xlabel(Xlabel); plt.ylabel(Ylabel); plt.legend(Legends,loc='center right')
    
    print "###############################################"
    print "# Ep_{CO3_{1f,np}}  = %.2f +- %.3f " % (mean_ep_co3_1_f,std_ep_co3_1_f)
    print "# Ep_{CO3_{2f,np}}  = %.2f +- %.3f " % (mean_ep_co3_2_f,std_ep_co3_2_f)
    print "# Ep_{CO3_{1,p}}}   = %.2f +- %.3f " % (m_co31,std_co31)
    print "# Ep_{CO3_{1f,p}}   = %.2f +- %.3f " % (m_co31f,std_co31f)
    print "# Ep_{CO3_{2f,p}}   = %.2f +- %.3f " % (m_co32f,std_co32f)
    print "# Ep_{SiO4-C-O_{p}} = %.2f +- %.3f " % (m_SiO4_C_O,std_SiO4_C_O)
    
    
    ############################################################################################
    #====================== Surface energy at higher temperatures ==============================
    ############################################################################################

    timesteps = [300000,350000,550000,600000,1100000,1150000,1650000,1700000] # timesteps for different stages during run
    ts = [] # indexes for timesteps for different stages during run
    for timestep in timesteps:
        ts.append(int(round(timestep/1000)))
    
    #print N_wb[0], N_ws[0], len(Time_wb), len(Time_ws), len(As_wb), len(As_ws)
    plt.figure()
    Title = "Potential energy"
    Xlabel = "time [ps]"
    Ylabel = "energy [kcal/mol]"
    Legends = ['%g x $Ep_{bulk,onecell}$' % (N_ws[0]/N_wb[0]),'$Ep_{slab}$','$Ep_{slab,np}$']
    rc = ['b-','r-','y-']
    plt.hold(True)
    plt.plot(Time_wb,Epp_wb,rc[0])
    plt.plot(Time_ws,Ep_ws,rc[1])
    plt.plot(Time_wsnp,Ep_wsnp,rc[2])
    plt.hold(False)
    plt.title(Title); plt.xlabel(Xlabel); plt.ylabel(Ylabel); plt.legend(Legends,loc='center right')
    

    plt.figure()
    start = 50
    cs = 500
    Title = "Temperature"
    Ylabel = "temperature [K]"
    Legends = ['$T_{bulk,onecell}$','$T_{slab}$','$T_{slab,np}$']
    rc = ['b-*','r-d','y-o']
    plt.hold(True)
    plt.plot(Time_wb,temp_wb,rc[0], markevery=40)  # bulk one-cell
    plt.plot(Time_ws,temp_ws,rc[1], markevery=50)  # slab
    plt.plot(Time_wsnp,temp_wsnp,rc[2], markevery=60)  # slab
    plt.plot(Time_ws,np.linspace(0,0,len(time_ws)),'k--')
    plt.plot(np.linspace(0,0,len(temp_ws)),np.linspace(0,temp_ws[0],len(temp_ws)),'k--')
    
    plt.annotate('T=300K', xy=(1.0, 1.0), xytext=(start,cs))
    plt.plot(np.linspace(Time_ws[ts[0]],Time_ws[ts[0]], 10), np.linspace(0,500,10),'k--')
    plt.plot(np.linspace(Time_ws[ts[1]],Time_ws[ts[1]], 10), np.linspace(0,500,10),'k--')
    
    plt.annotate('T=500K', xy=(1.0, 1.0), xytext=(start+30,cs+500))
    plt.plot(np.linspace(Time_ws[ts[2]],Time_ws[ts[2]], 10), np.linspace(0,1000,10),'k--')
    plt.plot(np.linspace(Time_ws[ts[3]],Time_ws[ts[3]], 10), np.linspace(0,1000,10),'k--')

    plt.annotate('T=1000K', xy=(1.0, 1.0), xytext=(start+80,cs+1000))
    plt.plot(np.linspace(Time_ws[ts[4]],Time_ws[ts[4]], 10), np.linspace(0,1500,10),'k--')
    plt.plot(np.linspace(Time_ws[ts[5]],Time_ws[ts[5]], 10), np.linspace(0,1500,10),'k--')

    plt.annotate('T=1500K', xy=(1.0, 1.0), xytext=(start+130,cs+1200))
    plt.plot(np.linspace(Time_ws[ts[6]],Time_ws[ts[6]], 10), np.linspace(0,1600,10),'k--')
    plt.plot(np.linspace(Time_ws[ts[7]-1],Time_ws[ts[7]-1], 10), np.linspace(0,1600,10),'k--')
    plt.hold(False)
    plt.title(Title); plt.xlabel(Xlabel); plt.ylabel(Ylabel); plt.legend(Legends,loc='upper left')
     
    Es_p = []; Es_np = []
    kcalmol2joule = 6.9477E-21
    eqA2msq = 1.0E-20
    for i in range(len(time_ws)):

        As = eqA2msq*As_ws[i]        # surface area of slab. fully periodic system
        As_np = eqA2msq*As_wsnp[i]   # surface area of nonperiodic slab.
        
        Eslab = kcalmol2joule*Ep_ws[i]            # pot E [J] for periodic slab
        Eslab_np = kcalmol2joule*Ep_wsnp[i]       # pot E [J] for nonperiodic slab
        
        Ebulk = kcalmol2joule*Ep_wb[i]*(float(N_ws[i])/N_wb[i])  # pot E [J] for bulk one-cell system 
        
        Es_p.append((Eslab - Ebulk)/(2*As))                      # surface energy for periodic slab and bulk system
        Es_np.append((Eslab_np - Ebulk)/(2*As_np))               # surf e for nonp slab and bulk systm
        
        
    #################################
    # Calculate surface energy at 300K, 500K, 1000K and 1500K:

    mean_Esp_300K = np.mean(Es_p[ts[0]:ts[1]])
    std_Esp_300K = np.std(Es_p[ts[0]:ts[1]])
    mean_Esp_500K = np.mean(Es_p[ts[2]:ts[3]])
    std_Esp_500K = np.std(Es_p[ts[2]:ts[3]])
    mean_Esp_1000K = np.mean(Es_p[ts[4]:ts[5]])
    std_Esp_1000K = np.std(Es_p[ts[4]:ts[5]])
    mean_Esp_1500K = np.mean(Es_p[ts[6]:ts[7]-1])
    std_Esp_1500K = np.std(Es_p[ts[6]:ts[7]-1])
    
    mean_Esnp_300K = np.mean(Es_np[ts[0]:ts[1]])
    std_Esnp_300K = np.std(Es_np[ts[0]:ts[1]])
    mean_Esnp_500K = np.mean(Es_np[ts[2]:ts[3]])
    std_Esnp_500K = np.std(Es_np[ts[2]:ts[3]])
    mean_Esnp_1000K = np.mean(Es_np[ts[4]:ts[5]])
    std_Esnp_1000K = np.std(Es_np[ts[4]:ts[5]])
    mean_Esnp_1500K = np.mean(Es_np[ts[6]:ts[7]-1])
    std_Esnp_1500K = np.std(Es_np[ts[6]:ts[7]-1])
        
    plt.figure()
    a = 0.45    # vertical line start x-hight
    b = 0.55    # vertical line stop x-hight
    c = 0.6    # yloc of annotate text (upper)
    d = 0.4     # yloc of annotate text (lower)
    start = 60
    Title = "Surface energy"
    Ylabel = "energy $[J/m^2]$"
    Legends = ['$Es_{p}$','$Es_{np}$']
    rc = ['r-','y-']
    plt.hold(True)
    plt.plot(Time_ws,Es_p,rc[0])
    plt.plot(Time_wsnp,Es_np,rc[1])

    plt.annotate('T=300K', xy=(1.0, 1.0), xytext=(start,c))
    plt.plot(np.linspace(Time_ws[ts[0]],Time_ws[ts[0]], 10), np.linspace(a,b,10),'k--')
    plt.plot(np.linspace(Time_ws[ts[1]],Time_ws[ts[1]], 10), np.linspace(a,b,10),'k--')
    
    plt.annotate('T=500K', xy=(1.0, 1.0), xytext=((start+20),d))
    plt.plot(np.linspace(Time_ws[ts[2]],Time_ws[ts[2]], 10), np.linspace(a,b,10),'k--')
    plt.plot(np.linspace(Time_ws[ts[3]],Time_ws[ts[3]], 10), np.linspace(a,b,10),'k--')

    plt.annotate('T=1000K', xy=(1.0, 1.0), xytext=((start+70),c))
    plt.plot(np.linspace(Time_ws[ts[4]],Time_ws[ts[4]], 10), np.linspace(a,b,10),'k--')
    plt.plot(np.linspace(Time_ws[ts[5]],Time_ws[ts[5]], 10), np.linspace(a,b,10),'k--')

    plt.annotate('T=1500K', xy=(1.0, 1.0), xytext=((start + 120) ,d))
    plt.plot(np.linspace(Time_ws[ts[6]],Time_ws[ts[6]], 10), np.linspace(a,b,10),'k--')
    plt.plot(np.linspace(Time_ws[ts[7]-1],Time_ws[ts[7]-1], 10), np.linspace(a,b,10),'k--')
    plt.hold(False)
    plt.title(Title); plt.xlabel(Xlabel); plt.ylabel(Ylabel); plt.legend(Legends)
    
    print "##################################################################"
    print "#                   SURFACE ENERGY"
    print "#-----------------------------------------------------------------"
    print "#  300K Es_o,p  =  %.6f +- %.7f [J/m^2]" % (mean_Esp_300K,std_Esp_300K)
    print "#  300K Es_o,np =  %.6f +- %.7f [J/m^2]" % (mean_Esnp_300K,std_Esnp_300K)
    print "#-----------------------------------------------------------------"
    print "#  500K Es_o,p  =  %.6f +- %.7f [J/m^2]" % (mean_Esp_500K,std_Esp_500K)
    print "#  500K Es_o,np =  %.6f +- %.7f [J/m^2]" % (mean_Esnp_500K,std_Esnp_500K)
    print "#-----------------------------------------------------------------"   
    print "# 1000K Es_o,p  =  %.6f +- %.7f [J/m^2]" % (mean_Esp_1000K,std_Esp_1000K)
    print "# 1000K Es_o,np =  %.6f +- %.7f [J/m^2]" % (mean_Esnp_1000K,std_Esnp_1000K)
    print "#-----------------------------------------------------------------"
    print "# 1500K Es_o,p  =  %.6f +- %.7f [J/m^2]" % (mean_Esp_1500K,std_Esp_1500K)
    print "# 1500K Es_o,np =  %.6f +- %.7f [J/m^2]" % (mean_Esnp_1500K,std_Esnp_1500K)
    print "##################################################################"

    ############################################################################################
    # =====================  Surface energy at ~0deg C. ========================================
    ############################################################################################            
    plt.figure()
    Title = "Potential energy"
    Xlabel = "time [ps]"
    Ylabel = "energy [kcal/mol]"
    Legends = ['$Ep_{bulk}$','%g x $Ep_{bulk,onecell}$' % (N1[0]/N3[0]),'$Ep_{slab}$', '$Ep_{slab,non-p}$','$Ep_{cleaved}$']
    rc = ['r-x','b-v','y-o','m-d','k-d']
    plt.hold(True)
    plt.plot(Time1,Ep1,rc[0],markevery=20)  # bulk
    plt.plot(Time3,Epp,rc[2],markevery=20)  # bulk one-cell
    plt.plot(Time2,Ep2,rc[1],markevery=20)  # slab
    plt.plot(Time4,Ep4,rc[3],markevery=25)  # slab non-periodic
    plt.plot(Time5,Ep5,rc[4],markevery=20)  # cleaved
    plt.hold(False)
    plt.title(Title); plt.xlabel(Xlabel); plt.ylabel(Ylabel); plt.legend(Legends)

    print "#HALLOHALLOHALLOHALLOHALLOHALLOHALLOHALLOHALLOHALLOHALLOHALLO#"
    print Epp[-1]
  
    #print len(Time1), len(Time3), len(temp1), len(temp3)
    plt.figure()
    Title = "Temperature"
    Ylabel = "temperature [K]"
    Legends = ['$T_{bulk}$','$T_{bulk,onecell}$','$T_{slab}$','$T_{slab,non-p}$','$T_{cleaved}$']

    rc = ['b-*','r-d','y-o','m--x','k--v']
    plt.hold(True)
    plt.plot(Time1,temp1,rc[0], markevery=20)  # bulk
    plt.plot(Time3,temp3,rc[2], markevery=22)  # bulk one-cell
    plt.plot(Time2,temp2,rc[1], markevery=35)  # slab
    plt.plot(Time4,temp4,rc[3], markevery=30)  # slab non-periodic
    plt.plot(Time5,temp5,rc[4], markevery=50)  # slab non-periodic
    plt.plot(Time2,np.linspace(0,0,len(time2)),'k--')
    plt.plot(np.linspace(0,0,len(temp2)),np.linspace(0,temp2[0],len(temp2)),'k--')
    plt.hold(False)
    plt.title(Title); plt.xlabel(Xlabel); plt.ylabel(Ylabel); plt.legend(Legends)

    #print len(As5), len(As4), len(As3), len(As2)
    Esurf_p = []; Esurf_onecell_p = []; Esurf_nonp = []; Esurf_onecell_nonp = []
    Esurf_oc = []; Esurf_bc = []

    for i in range(len(time4)):

        As = eqA2msq*As2[i]        # surface area of slab. fully periodic system
        As_nonp = eqA2msq*As4[i]   # surface area of slab. non-periodic in z dir
        
        Eslab = kcalmol2joule*Ep2[i]            # pot E [J] for periodic slab
        Eslab_nonp = kcalmol2joule*Ep4[i]       # pot E [J] for non-periodic slab
        
        Ebulk = kcalmol2joule*Ep1[i]*(float(N1[i])/N2[i])          # pot E [J] for bulk
        Ebulk_onecell = kcalmol2joule*Ep3[i]*(float(N1[i])/N3[i])  # pot E [J] for bulk one-cell system 
        
        Esurf_p.append((Eslab - Ebulk)/(2*As))                              # surface energy for periodic slab and bulk system
        
        Esurf_onecell_p.append((Eslab - Ebulk_onecell)/(2*As))              # surface energy for periodic slab and one-cell bulk system
         
        Esurf_nonp.append((Eslab_nonp - Ebulk)/(2*As_nonp))                 # surface energy for non-periodic and bulk system
    
        Esurf_onecell_nonp.append((Eslab_nonp - Ebulk_onecell)/(2*As_nonp)) # surface enerfy for non-periodic and one-cell bulk system
        
        Ecleaved = kcalmol2joule*Ep5[i]
        As_c = eqA2msq*As5[i]
        
        Ebulk = kcalmol2joule*Ep1[i]*(float(N5[i])/N2[i])          # pot E [J] for bulk
        Ebulk_onecell = kcalmol2joule*Ep3[i]*(float(N5[i])/N3[i])  # pot E [J] for bulk one-cell system 
        
        Esurf_bc.append((Ecleaved - Ebulk)/(2*As_c))
        Esurf_oc.append((Ecleaved - Ebulk_onecell)/(2*As_c))
        
        
        #------ End for-loop ----------------------------------------------#

    from mpl_toolkits.axes_grid.inset_locator import zoomed_inset_axes, inset_axes, mark_inset

    plt.figure()
    ax = plt.axes()
    Title = "Surface energy"
    Ylabel = "energy $[J/m^2]$"
    Legends = ['$Es_{b,p}$','$Es_{o,p}$','$Es_{b,np}$','$Es_{o,np}$','$Es_{b,c}$']
    rc = ['r-','r--','y-','y--','m-','m--']
    plt.hold(True)
    plt.plot(Time4,Esurf_p,rc[0])             # Es periodic slab system and bulk system
    plt.plot(Time4,Esurf_onecell_p, rc[1])    # Es periodic slab system and one-cell bulk system 
    plt.plot(Time4,Esurf_nonp,rc[2])          # Es non-periodic slab system and bulk system 
    plt.plot(Time4,Esurf_onecell_nonp,rc[3])  # Es non-periodic slab system and one-cell bulk system 
    plt.plot(Time4,Esurf_bc,rc[4])            # Es cleaved surface. bulk system
    plt.plot(Time4,Esurf_oc,rc[5])            # Es cleaved surface. one-cell system
    plt.hold(False)
    plt.title(Title); plt.xlabel(Xlabel); plt.ylabel(Ylabel); plt.legend(Legends)    

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()
    
    axin = inset_axes(ax, width="30%", height="30%", loc=5)
    plt.hold(True)
    axin.plot(Time4,Esurf_p,rc[0])
    axin.plot(Time4,Esurf_onecell_p,rc[1])
    axin.plot(Time4,Esurf_nonp,rc[2])
    axin.plot(Time4,Esurf_onecell_nonp,rc[3])
    axin.plot(Time4,Esurf_bc,rc[4])
    axin.plot(Time4,Esurf_oc,rc[5])
    plt.hold(False)
    axin.set_xlim(28.0, 43.0)
    ya = 0.432
    yb = 0.439
    Npoints = 6
    points = np.linspace(ya,yb,Npoints)
    axin.set_ylim(ya, yb)
    axin.set_xticks([])
    axin.set_yticks(points)
    mark_inset(ax, axin, loc1=1, loc2=4, fc="none", ec="0.5")
#    mark_inset(ax, axin, loc1=1, loc2=1)
    
    # calculate mean and standard deviations!
    ltime = len(Time1)                  # number of timesteps in total. There are 43ps in total 43/ltime = ps/timestep
    start = int(round(ltime*(30/43.0))) # timestep after 30ps
    mean_op = np.mean(Esurf_onecell_p[start:])
    std_op = np.std(Esurf_onecell_p[start:]) 
    mean_bp = np.mean(Esurf_p[start:])
    std_bp = np.std(Esurf_p[start:])
    mean_onp = np.mean(Esurf_onecell_nonp[start:])
    std_onp = np.std(Esurf_onecell_nonp[start:])
    mean_bnp = np.mean(Esurf_nonp[start:])
    std_bnp = np.std(Esurf_nonp[start:])
    mean_oc = np.mean(Esurf_oc[start:])
    std_oc = np.std(Esurf_oc[start:])  
    mean_bc = np.mean(Esurf_bc[start:])
    std_bc = np.std(Esurf_bc[start:])    
    
    ##########################################################################
    # ---- Print som results:    
    print "#                   SURFACE ENERGY"
    print "# Es_b,p  = %.6f +- %.g [J/m^2]" % (mean_bp,std_bp)
    print "# Es_o,p  = %.6f +- %.g [J/m^2]" % (mean_op,std_op)
    print "# Es_b,np = %.6f +- %.g [J/m^2]" % (mean_bnp,std_bnp)
    print "# Es_o,np = %.6f +- %.g [J/m^2]" % (mean_onp,std_onp)    
    print "# Es_b,c  = %.6f +- %.g [J/m^2]" % (mean_bc,std_bc)
    print "# Es_o,c  = %.6f +- %.g [J/m^2]" % (mean_oc,std_oc)
    print "####################################################################################"
    print "# Number of atoms        Surface area          System hight       cut-off distance  "
    print "# --------------------------------------------------------------------------------- "
    print "# N_o =    %g atoms   As_o =  %g [A^2]   Lz_o =  %g [A]   Rc_o =  %g [A]" % (N3[-1],As3[-1],Lz3[-1],Rc3[-1])
    print "# N_b =  %g atoms   As_b =  %g [A^2]   Lz_b =  %g [A]   Rc_b =  %g [A]" % (N1[-1],As1[-1],Lz1[-1],Rc1[-1])
    print "# N_p =  %g atoms   As_p =  %g [A^2]   Lz_p =  %g [A]   Rc_p =  %g [A]" % (N2[-1],As2[-1],Lz2[-1],Rc2[-1])
    print "# N_np = %g atoms   As_np = %g [A^2]   Lz_np = %g [A]   Rc_np = %g [A]" % (N4[-1],As4[-1],Lz4[-1],Rc4[-1])
    print "# N_c =  %g atoms   As_c =  %g [A^2]   Lz_c =  %.3f [A]   Rc_c =  %g [A]" % (N5[-1],As5[-1],Lz5[-1],Rc5[-1])
    print "####################################################################################"


    ###########################################################################
    # Calculate bindingenergy of carbondioxide to SiO2 surfaces
    #-------------------------------------------------------------------------#
    # Eb = (Ep_{SiO2 + N*(CO2)} - Ep_{cleaved} - Ep_{N*(CO2)})/N_{CO2}

    Ep_cleaved_ = np.mean(Ep5[lstart:])      # kcal/mol
    Ep_cleaved = kcalmol2joule*Ep_cleaved_   # Joule
    std_Ecleaved = kcalmol2joule*np.std(Ep5[lstart:])
    Ep_slab_ = np.mean(Ep3[lstart:])      # kcal/mol
    Ep_slab = kcalmol2joule*Ep_slab_      # Joule
    std_Eslab = kcalmol2joule*np.std(Ep5[lstart:])
    Ep_CO2_8 = []; Ep_CO2_1 = []
    for i in range(len(Time13[lstart:])):
        ii = lstart + i
        Ep_CO2_8.append(8*Ep13[ii])
        Ep_CO2_1.append(1*Ep13[ii])
    
    plt.figure()
    plt.hold(True)
    plt.plot(Time13[lstart:],Ep_CO2_1,'b-x')
    plt.plot(Time13[lstart:],Ep_CO2_8,'b-o')
    plt.plot(Time5[lstart:len(Ep9)],Ep5[lstart:len(Ep9)],'k--d', markevery=3)
    plt.plot(Time9[lstart:],Ep9[lstart:],'r--o', markevery=5)
    plt.plot(Time10[lstart:],Ep10[lstart:],'m--v', markevery=4)
    plt.plot(Time11[lstart:],Ep11[lstart:],'m--d', markevery=3)
    plt.plot(Time12[lstart:],Ep12[lstart:],'r--x', markevery=2)
    plt.hold(False)
    plt.legend(['Ep_CO2_1','Ep_CO2_8','Ep_{c}','Ep_{co3-1}','Ep_{c03-1f}','Ep_{co3-2f}','Ep_{sio4-c-o}'],loc='center left')
    plt.xlabel('time [ps]'); plt.ylabel('potential energy [kcal/mol]')    
    
    
    
    N_CO3_1 = 1            # N CO2 molecules in CO3_1 system
    N_CO3_1f = 8           # N CO2 molecules in CO3_1f system
    N_CO3_2f = 8           # N CO2 molecules in CO3_2f system
    N_SiO4_CO = 1          # N CO2 molecules in SiO4_C_O system
    
    eb_CO3_1 = 0
    eb_CO3_1_f = 0 
    eb_CO3_2_f = 0
    eb_SiO4_C_O = 0
    
    nvalues = len(Time13[lstart:])-1  # number of time values to average over
    for i in range(nvalues):
        j = i + lstart
        eb_CO3_1 += (Ep9[j] - Ep5[j] - N_CO3_1*Ep13[j])/N_CO3_1
        eb_CO3_1_f += (Ep10[j] - Ep5[j] - N_CO3_1*Ep13[j])/N_CO3_1f
        eb_CO3_2_f += (Ep11[j] - Ep5[j] - N_CO3_1*Ep13[j])/N_CO3_2f
        eb_SiO4_C_O += (Ep12[j] - Ep5[j] - N_CO3_1*Ep13[j])/N_SiO4_CO

    eb_CO3_1 = eb_CO3_1/nvalues
    eb_CO3_1_f = eb_CO3_1_f/nvalues
    eb_CO3_2_f = eb_CO3_2_f/nvalues
    eb_SiO4_C_O = eb_SiO4_C_O/nvalues

    

    print "##############################################################"
    print "#  ----------------  Binding energies  ----------------------#"
    print "#  eb_CO3_1    = %.3f [kcal/mol]" % eb_CO3_1
    print "#  eb_CO3_1_f  = %.3f [kcal/mol]" % eb_CO3_1_f
    print "#  eb_CO3_2_f  = %.3f [kcal/mol]" % eb_CO3_2_f
    print "#  eb_SiO4_C_O = %.3f [kcal/mol]" % eb_SiO4_C_O
    print "#  ----------------  Binding energies  ----------------------#"
    print "#  eb_CO3_1    = %g [Joule]" % (eb_CO3_1*kcalmol2joule)
    print "#  eb_CO3_1_f  = %g [Joule]" % (eb_CO3_1_f*kcalmol2joule)
    print "#  eb_CO3_2_f  = %g [Joule]" % (eb_CO3_2_f*kcalmol2joule)
    print "#  eb_SiO4_C_O = %g [Joule]" % (eb_SiO4_C_O*kcalmol2joule)
    '''
    Ep_CO2_p_ = (np.mean(Ep13[lstart:]))  # kcal/mol
    Ep_CO2_p = kcalmol2joule*Ep_CO2_p_    # Joule
    std_EpCO2_p = kcalmol2joule*(np.std(Ep13[lstart:]) )
    
    Ep_CO3_1 = kcalmol2joule*m_co31         # Epot of CO3_1 at 0K
    Ep_CO3_1f = kcalmol2joule*m_co31f       #
    Ep_CO3_2f = kcalmol2joule*m_co32f       #
    Ep_SiO4_CO = kcalmol2joule*m_SiO4_C_O   #
    
    std_Ep_CO3_1 = kcalmol2joule*std_co31        # standard deviations
    std_Ep_CO3_1f = kcalmol2joule*std_co31f
    std_Ep_CO3_2f = kcalmol2joule*std_co32f
    std_Ep_SiO4_CO = kcalmol2joule*std_SiO4_C_O

    std_Eb_CO3_1 = np.sqrt(std_Ep_CO3_1**2 + std_Ecleaved**2 + std_EpCO2_p)      # standard deviations
    std_Eb_CO3_1f = np.sqrt(std_Ep_CO3_1f**2 + std_Ecleaved**2 + std_EpCO2_p)
    std_Eb_CO3_2f = np.sqrt(std_Ep_CO3_2f**2 + std_Ecleaved**2 + std_EpCO2_p)
    std_Eb_SiO4_CO = np.sqrt(std_Ep_SiO4_CO**2 + std_Ecleaved**2 + std_EpCO2_p)
    
    Eb_CO3_1 = (Ep_CO3_1 - Ep_cleaved - N_CO3_1*Ep_CO2_p)/N_CO3_1                # binding energy of carbon dioxide   [Joule]
    Eb_CO3_1f = (Ep_CO3_1f - Ep_cleaved - N_CO3_1f*Ep_CO2_p)/N_CO3_1f
    Eb_CO3_2f = (Ep_CO3_2f - Ep_cleaved - N_CO3_2f*Ep_CO2_p)/N_CO3_2f
    Eb_SiO4_CO = (Ep_SiO4_CO - Ep_cleaved - N_SiO4_CO*Ep_CO2_p)/N_SiO4_CO

    Eb_CO3_1_ = (m_co31 - Ep_cleaved_ - N_CO3_1*Ep_CO2_p_)/N_CO3_1               # binding energy of carbon dioxide  [kcal/mol]
    Eb_CO3_1f_ = (m_co31f - Ep_cleaved_ - N_CO3_1f*Ep_CO2_p_)/N_CO3_1f
    Eb_CO3_2f_ = (m_co32f - Ep_cleaved_ - N_CO3_2f*Ep_CO2_p_)/N_CO3_2f
    Eb_SiO4_CO_ = (m_SiO4_C_O - Ep_cleaved_ - N_SiO4_CO*Ep_CO2_p_)/N_SiO4_CO
    
    Eb_CO3_1_s = (Ep_CO3_1 - Ep_slab - N_CO3_1*Ep_CO2_p)/N_CO3_1                 # binding energy of carbon dioxide   [Joule]
    Eb_CO3_1f_s = (Ep_CO3_1f - Ep_slab - N_CO3_1f*Ep_CO2_p)/N_CO3_1f
    Eb_CO3_2f_s = (Ep_CO3_2f - Ep_slab - N_CO3_2f*Ep_CO2_p)/N_CO3_2f
    Eb_SiO4_CO_s = (Ep_SiO4_CO - Ep_slab - N_SiO4_CO*Ep_CO2_p)/N_SiO4_CO
    
    std_Eb_CO3_1_s = np.sqrt(std_Ep_CO3_1**2 + std_Eslab**2 + std_EpCO2_p)       # standard deviations
    std_Eb_CO3_1f_s = np.sqrt(std_Ep_CO3_1f**2 + std_Eslab**2 + std_EpCO2_p)
    std_Eb_CO3_2f_s = np.sqrt(std_Ep_CO3_2f**2 + std_Eslab**2 + std_EpCO2_p)
    std_Eb_SiO4_CO_s = np.sqrt(std_Ep_SiO4_CO**2 + std_Eslab**2 + std_EpCO2_p)
    
    
#    print Ep_CO3_1, Ep_CO3_1f, Ep_CO3_2f, Ep_SiO4_CO
#    print m_co31, m_co31f, m_co32f, m_SiO4_C_O
#    print Ep_cleaved, Ep_CO2_p
    print "####################################################################################"
    print "#                     Surface binding energy of CO2                                #"
    print "# ---------------------------------------------------------------------------------#"
    print "# --------------------          Cleaved surface          --------------------------#"
    print "# Eb_{CO3-1}    = %g +- %g [J] " % (Eb_CO3_1,std_Eb_CO3_1)
    print "# Eb_{CO3-1-f}  = %g +- %g [J] " % (Eb_CO3_1f,std_Eb_CO3_1f)
    print "# Eb_{CO3-2-f}  = %g +- %g [J] " % (Eb_CO3_2f,std_Eb_CO3_2f)
    print "# Eb_{SiO4-C-O} = %g +- %g [J] " % (Eb_SiO4_CO,std_Eb_SiO4_CO)
    print "# ---------------------------------------------------------------------------------#"
    print "# -------------------         27 layer slab surface        ------------------------#"
    print "# Eb_{CO3-1}    = %g +- %g [J] " % (Eb_CO3_1_s,std_Eb_CO3_1_s)
    print "# Eb_{CO3-1-f}  = %g +- %g [J] " % (Eb_CO3_1f_s,std_Eb_CO3_1f_s)
    print "# Eb_{CO3-2-f}  = %g +- %g [J] " % (Eb_CO3_2f_s,std_Eb_CO3_2f_s)
    print "# Eb_{SiO4-C-O} = %g +- %g [J] " % (Eb_SiO4_CO_s,std_Eb_SiO4_CO_s)
    print "# ---------------------------------------------------------------------------------#"
    #print "# Eb_{CO3-1}    = %g [kcal/mol] " % (Eb_CO3_1_)
    #print "# Eb_{CO3-1-f}  = %g [kcal/mol] " % (Eb_CO3_1f_)
    #print "# Eb_{CO3-2-f}  = %g [kcal/mol] " % (Eb_CO3_2f_)
    #print "# Eb_{SiO4-C-O} = %g [kcal/mol] " % (Eb_SiO4_CO_)
    #print "# ---------------------------------------------------------------------------------#"
    #print "# Eb_{CO3-1}    = %g [J] " % (kcalmol2joule*Eb_CO3_1_)
    #print "# Eb_{CO3-1-f}  = %g [J] " % (kcalmol2joule*Eb_CO3_1f_)
    #print "# Eb_{CO3-2-f}  = %g [J] " % (kcalmol2joule*Eb_CO3_2f_)
    #print "# Eb_{SiO4-C-O} = %g [J] " % (kcalmol2joule*Eb_SiO4_CO_)
    print "#####################################################################################"
    print "# Potential energies:"
    print "# ----------------------------------------------------------------------------------#"
    print "# Ep_CO2_1   = %.3f [kcal/mol] " % (Ep_CO2_p_)
    print "# Ep_Cleaved = %.3f [kcal/mol] " % (Ep_cleaved_)
    print "# Ep_27slab  = %.3f [kcal/mol] " % (Ep_slab_)
    print "# Ep_CO3_1   = %.3f [kcal/mol] " % (m_co31)
    print "# Ep_CO3_1_f = %.3f [kcal/mol] " % (m_co31f)
    print "# Ep_CO3_2_f = %.3f [kcal/mol] " % (m_co32f)
    print "# Ep_SiO4_CO = %.3f [kcal/mol] " % (m_SiO4_C_O)
    print "#####################################################################################"
    '''
    
    print "-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-."
    E_slab = kcalmol2joule*Ep2[-1]
    E_bulk = kcalmol2joule*Ep1[-1]*(float(N1[-1])/N2[-1])
    A_surf = eqA2msq*As2[-1]
    Esurf = (Eslab - Ebulk)/(2*A_surf)
    print " Eslab = %g" % (E_slab)
    print " Ebulk = %g" % (E_bulk)
    print " As    = %g" % (A_surf)
    print " Esurf = %.5f" % (Esurf)
    plt.show()
    
if (__name__ == "__main__"):
    import matplotlib.pyplot as plt
    import numpy as np    
    import os
    main()