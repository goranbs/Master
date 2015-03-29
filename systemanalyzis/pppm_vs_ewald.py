# -*- coding: utf-8 -*-
"""
Created on Wed Feb 18 12:02:56 2015

@author: goran
"""



# plot scaling of pppm vs ewald

import numpy as np
import matplotlib.pyplot as plt

N = np.linspace(0,100,100)

pppm = N*np.log(N)

ewald = N**(3.0/2)


plt.figure()
plt.hold(True)
plt.plot(N,pppm)
plt.plot(N,ewald)
plt.hold(False)
plt.legend(['pppm','ewald'],loc='upper left')
