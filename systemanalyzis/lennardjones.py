# -*- coding: utf-8 -*-
"""
Created on Wed Jan 28 11:38:24 2015

@author: goran
"""

def LennardJones(r,epsilon,sigma=None,R=None):
    
    
    E1 = np.zeros(np.shape(r))
    E2 = np.zeros(np.shape(r))
    LJ12term = []
    LJ6term = []
    if (R != None):
        for i in range(len(r)):
            lj12term = (R/r[i])**12
            lj6term = 2*(R/r[i])**6
            LJ12term.append(lj12term)
            LJ6term.append(lj6term)
            E1[i] = epsilon*(lj12term - lj6term)
    if ((sigma != None) and (R == None)):
            E2[i] = 4*epsilon*((sigma/r[i])**12 - (sigma/r[i])**6)
        
    if (R == sigma):
        print "Error! R=sigma, R = sigma^(1/6)\n you have to give either R or sigma!"
        return 1
    
    if ((R != None) and (sigma != None)):
        eps_line = np.linspace(-epsilon,0,4)
        sigma_line = np.linspace(0,sigma,4)
        zero_line = np.linspace(0,r[-1],4)
        plt.figure()
        plt.hold(True)
        plt.plot(r,E1,'b-o')
        plt.plot(zero_line,[0,0,0,0],'r-d')
        plt.plot(sigma_line,[0,0,0,0],'--y',linewidth=3.0)
        plt.plot([R,R,R,R],eps_line,'--k',linewidth=3.0)
        plt.annotate('R', xy=(1.0, 1.0), xytext=(1.125, 0.1))
        plt.annotate('$\sigma$', xy=(1.0,1.0), xytext=(0.4, 0.1))
        plt.annotate('$\epsilon$', xy=(1.0, 1.0), xytext=(1.2, -0.25))
        plt.hold(False)        
        plt.title('Lennard-Jones')
        plt.xlabel('radius'); plt.ylabel('Energy')


    plt.figure()
    plt.hold(True)
    plt.plot(r,LJ12term,'r-')
    plt.plot(r,LJ6term,'b-')
    plt.legend(['lj12-term','lj6-term'])
    plt.hold(False)
    plt.show(True)
    plt.close('all')
    
def main():
    eps = 1.0
    sig = 1.0
    R = sig*2**(1/6.0)
    r = np.linspace(0.97,3,100)
    LennardJones(r,eps,sig,R)
    

if (__name__ == "__main__"):
    import matplotlib.pyplot as plt
    #plt.rc('text', usetex=True) 
    import numpy as np
    main()