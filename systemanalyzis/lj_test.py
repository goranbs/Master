# -*- coding: utf-8 -*-
"""
Created on Wed Jan 14 11:15:25 2015

@author: goran
"""

def LAMMPS_lj(epsilon,sigma,r):
    
    E = []
    for i in range(len(r)):
        Enrj = 4*epsilon*((sigma)/r[i])**12 - ((sigma/r[i])**6)
        E.append(Enrj)
    
    return E

    
def Cygan_lj(D0,R0,r):
    E = []
    for i in range(len(r)):
        Enrj = D0*((R0)/r[i])**12 - 2*((R0/r[i])**6)
        E.append(Enrj)
        
    return E
    
def compare(f1,f2,r):
    Title = 'compareing functions f1 and f2'
    Legends = ['f1','f2']
    Xlabel = 'x-axis'
    Ylabel = 'Y-axis'
    plotting(Title,Legends,Xlabel,Ylabel,[r,r],[f1,f2])
    

def plotting(title,legends,xlabel,ylabel,xarrays,yarrays):
    plt.close('all')
    plt.figure()
    plt.hold(True)
    for i in range(len(xarrays)):
        plt.plot(xarrays[i],yarrays[i])

    plt.hold(False)
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.legend(legends)
    plt.show(True)
    plt.close('all')
        
def main():
    
    r0 = 2.0
    r1 = 4.0
    Npoints = 1000
    D = 0.1553 
    R0a = 3.166
    R0b = 3.55
    Title = 'Cygan vs LAMMPS lj potential. D=%g' % D
    Xlabel = 'r'
    Ylabel = 'Energy'
    Legends = ['LAMMPS Ra=%g' % (R0a),'Cygan Ra=%g' % (R0a),'LAMMPS Rb=%g' % (R0b),'Cygan Rb=%g' % (R0b)] 
    
    
    r = np.linspace(r0,r1,Npoints)
    
    Elammps_a = LAMMPS_lj(D,R0a,r)
    Ecygan_a = Cygan_lj(D,R0a,r)
    Elammps_b = LAMMPS_lj(D,R0b,r)
    Ecygan_b = Cygan_lj(D,R0b,r)
    
    plotting(Title,Legends,Xlabel,Ylabel,[r,r,r,r],[Elammps_a,Ecygan_a,Elammps_b,Ecygan_b])    
    
    r = np.linspace(0,np.pi/5.0,100)
    R = 2.5
    S = 2.8
    D = S-R
    f1 = np.cos(np.pi*(r-R)/(S-R))
    f2 = -np.sin(np.pi*(r-R)/(2*D))
    compare(f1,f2,r)
    
if (__name__ == "__main__"):
    import numpy as np
    import matplotlib.pyplot as plt
    main()