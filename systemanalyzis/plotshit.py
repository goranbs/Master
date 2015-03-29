# -*- coding: utf-8 -*-
"""
Created on Thu Dec  4 14:25:13 2014

@author: goran

Plotting program. 

Takes a datafile as argument. The datafile contains information about the plots that
should be made.
title
xlabel
ylabel
legends
linestyles
data-arrays
"""
def printhelp(arg):
    """
    Print information about plotshit
    """
    print "###################################################################"
    print "Name of program:    %s   " % arg
    print "Useage: python %s [path-to-file] [showplots] [saveplots] \n " % arg
    print "%s takes at least one cmd argument [path-to-file]\n " % arg
    print "showplot and saveplot is optional values (yes/no)\n "
    print "The inputfile is a textfile, formatted the following way:\n "
    print " title <figuretitle>\n xlabel <xlabelname>\n ylabel <ylabelname>\n figurename <output-figure-name>\n linestyle graph1 graph2 ... graphn\n legends graph1 graph2 ..... graphn\n x-array1 x11 x12 x13 ......... x1k\n y-array1 y11 y12 y13 ......... y1k\n y-array2 y21 y22 y23 ......... y2k\n .      .   .   .     .........  .\n .      .   .   .     .........  .\n y-arrayn yn1 yn2 yn3 ......... ynk\n "
    string = '\\newline'
    print "You can also specify one x-array for every y-array.\n Then, x-array1 belongs to y-array1, and so on.\n"
    print "Use:\n       %s for newline\n       $ signs for math expressions" % (string)
    print "Whitespaces separeates the input variables,\n this is most often the error encountered in the legends option"
    print "###################################################################"
    
def get_savepath(currentpath,pathgiven,truevalues):
    #print "path is from the submit-folder. That is: %s \n" % currentpath
    saveindatabase = raw_input("Do you want to save the plots in the directory you gave on the cmd (yes/no)? ")
    if (saveindatabase in truevalues):
        fullpathtosavedir = os.path.abspath(os.path.join(currentpath,pathgiven))
        fullpathtosavedir = os.path.dirname(fullpathtosavedir)
        print "save in %s" % fullpathtosavedir
    else:    
        print "Where do you want to save the plots?"
        savepath = raw_input("Give path: ")
        if (savepath[0] == "~"):
            home = os.path.expanduser("~")
            savepath = home + savepath[1:]
            fullpathtosavedir = os.path.abspath(savepath)
        else: 
            fullpathtosavedir = os.path.abspath(os.path.join(currentpath,savepath))
            print "path chosen: %s " % fullpathtosavedir
    
        if not os.path.isdir(fullpathtosavedir): # if the path is unique, do you wish to create it?
            createfolder = raw_input("not a valid path. Do you wish to create a folder here? (yes/no) ")
            if (createfolder in ["Y","y","Yes","yes","true","True"]): 
                print "creating folder...."
                os.makedirs(fullpathtosavedir) # use default mode
            else:
                print "Ok then..."
                fullpathtosavedir = get_savepath(currentpath)
        else: 
            print "The path '%s'  is an existing path," % fullpathtosavedir
            savehere = raw_input("Do you wish to save plots here? (yes/no) ")
            if not (savehere in ["Y","y","Yes","yes","true","True","jepp","Jepp"]):
                print "Ok then..."
                fullpathtosavedir = get_savepath(currentpath)
    
    #print "save in: %s " % fullpathtosavedir
    return fullpathtosavedir


def main(argv):
    """
    The program should be called from the commandline with the full path to the
    input data file should be provided as a commandline argument.
    
    masterplot.py ~/Goran/Master/plotinfo_1.dat
    """
    pathtofile = argv[1]  # path from the directory where this python program is run from
    pathtocurrent = os.getcwd()
    fullpathtofile = os.path.abspath(os.path.join(pathtocurrent,pathtofile))
    #print "pathtofile=     %s" % pathtofile
    #print "pathtocurrent = %s" % pathtocurrent
    #print "fullpathtofile= %s" % fullpathtofile
    truevalues = ["Y","y","Yes","yes","true","True"]
    
    if (argv[1] in ["-h", "h", "Help", "help", "HELP", "--help"]):
        printhelp(argv[0])
        sys.exit(1)
    
    if ((len(argv) < 2) or (len(argv) > 4)):
        sys.stderr.write("Useage:python %s <path-to-plotfile> <showplot=yes/no/None> <saveplot=yes/no/None>" % (argv[0]))
        return 1
    if not os.path.exists(fullpathtofile):
        sys.stderr.write("Error: Database %s was not found!" % (argv[1]))
        return 1
    if ((len(argv) == 3) or (len(argv) == 4)):
        showplot = argv[2]
        if (showplot in truevalues):
            showplot = True
            print "Show plots = True"
        if (showplot in ["N","n","No","no","false","False"]):
            showplot = False
            print "Show plots = False"
        if (len(argv) == 3):
            print "Commandline argument for saveplot not given. saveplot = False"
    if (len(argv) == 4):
        saveplot = argv[3]
        if (saveplot in truevalues):
            saveplot = True
            print "Save plots = True"
            pathtosavedir = get_savepath(pathtocurrent,pathtofile,truevalues)
        if (saveplot in ["N","n","No","no","false","False"]):
            saveplot = False
            print "Save plots = False"
        
                
    
    ########################################################################
    #                     The plotting program                             #
    
    import matplotlib.pyplot as plt
    plt.rc('font',**{'family':'sans-serif','sans-serif':['Helvetica'], 'size' : 13})
    plt.rc('text', usetex=True) 
    
    ofile = open(fullpathtofile,'r')
    

    nlines = 0
    nfigures = 0
    fig = []
    no = 1
    while 1:
        line = ofile.readline()
        if not line: break
        line = line.split()
        if ("title" in line[0]):
            nfigures += 1
            if no > 1:
                fig.append(no)
            no = 1
        else: 
            no += 1
        
        nlines += 1
    fig.append(no)
    
    #print nlines, nfigures, fig
    
    ofile.seek(0) # go back to beginning of file!
    
    plots = []
    
    for fignr in range(nfigures):
        aplot = {"Title": None,
                  "Xlabel": None,
                  "Ylabel" : None,
                  "Plotname" : None,
                  "Fs" : [],
                  "Linestyles" : [],
                  "xarray" : [],
                  "x-arrange" : [],
                  "yarray" : [],
                  "y-arrange" : []}
                  
        for lines in range(fig[fignr]):
            wholeline = ofile.readline()
            line = wholeline.split()
      
            if ("title" in line[0]):                              # Title
                tit = wholeline[(len(line[0])+1):-1]
                aplot["Title"] = tit
            if ("xlabel" in line[0]):                             # Xlabel
                ylab = wholeline[(len(line[0])+1):-1]
                aplot["Xlabel"] = ylab
            if ("ylabel" in line[0]):                             # Ylabel
                ylab = wholeline[(len(line[0])+1):-1]
                aplot["Ylabel"] = ylab
            if ("plotname" in line[0]):                           # Plotname
                aplot["Plotname"] = line[1]
            if ("legends" in line[0]):                            # Legends
                aplot["Legends"] = line[1:]
            if ("linestyles" in line[0]):                         # Linestyles
                aplot["Linestyles"] = line[1:]
            if ("x-array" in line[0]):                            # xarray
                array_nr = int(line[0][7:])
                xarray = [float(x) for x in line[1:]]
                aplot["xarray"].append(xarray)
                aplot["x-arrange"].append(array_nr)
            if ("y-array" in line[0]):                            # yarray
                array_nr = int(line[0][7:])
                yarray = [float(y) for y in line[1:]]
                aplot["yarray"].append(yarray)
                aplot["y-arrange"].append(array_nr)
    
        plots.append(aplot)
    
    ofile.close()
    
    ##########################################################################
    #                              PLOTTING                                  #
    ##########################################################################
    
    
    counter = 0    
    for aplot in plots:
        
        if (len(aplot["Linestyles"]) == 0 ):
            colours = ['b','r','g','y','m','k']
            lines = ['-','--','-x','-d','-o','.','-v']
            aplot["Linestyles"] = []
            for line in lines:
                for colour in colours:
                    aplot["Linestyles"].append(colour + line)
        counter += 1
        print "plotting plot number: %g " % counter
        plt.figure()
        plt.hold(True)
        Nxarrays= len(aplot["xarray"])
        Nyarrays= len(aplot["yarray"])
        if (Nxarrays == Nyarrays):  # new part not tested !!!! could be that x1 is not plotted with y1...
            for i in aplot["x-arrange"]:
                for j in aplot["y-arrange"]:
                    if (j == i):
                        linestyle = aplot["Linestyles"][j-1]
                        plt.plot(aplot["xarray"][i-1], aplot["yarray"][j-1],linestyle)
        elif (Nxarrays == 1):
            for i in range(Nyarrays):
                linestyle = aplot["Linestyles"][i]
                #print "plotting %s yarray number %g" % (aplot["Title"],i)
                plt.plot(aplot["xarray"][0], aplot["yarray"][i],linestyle)
        plt.hold(False)
        plt.title(aplot["Title"])
        plt.xlabel(aplot["Xlabel"]); plt.ylabel(aplot["Ylabel"])
        if (len(aplot["Legends"]) == Nyarrays):
            plt.legend(aplot["Legends"],loc="upper right")
        else:
            print "################################################\n Warning!!!\n Error occured when adding legends!"
            print "number of legends = %g" % (len(aplot["Legends"]))
            print "number of graphs  = %g" % (len(aplot["yarray"]))
            print "\n No legends added to the plot!"
        if (saveplot):
            name = aplot["Plotname"] + '.png'
            print " (%g/%g) Saving plot... plotname: %s" % (counter,len(plots),name)
            fname = os.path.join(pathtosavedir,name)
            plt.savefig(fname,format='png')
    
    plt.show(showplot)
    plt.close('all')
    print "Done saving all the plots, exiting successfully...\n                          .... At least as far as I know :-)"
    
    
if (__name__ == "__main__"):
    import sys, os
    sys.exit(main(sys.argv))

