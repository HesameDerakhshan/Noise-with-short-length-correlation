# -*- coding: utf-8 -*-
''' 
This code Genarate two noise:  1.white noise ƞ  2. a noise with shot-length correlation ξ  , v 0.9                                                *
Date: 11.09.1399  
                                                  
Copyleft 2020 H. Derakhshan All left reserved!                                                
                                                                            
Licensed under The GNU General Public License v3.0                
--------------------------------------------------------------------------------
'''
# General Modules

import numpy as np
import matplotlib.pyplot as plt

#   This function gets x argument (descret pionts of the lattice)، Then return  ƞ and ξ noise  
def noise(x):                     
    import math                   # For math functions
    N    = len(x)                 # N is total descret pionts of the lattice
    beta = 0 * 0.5                # Beta 0 White nose, 1 Purple noise,2 Brownian motion
    xi  = np.empty(N)             # Array of the ξ noises 
    eta = np.empty(N)             # Array of the ƞ noises
    Phi  = np.random.rand(int(N/2))*2*math.pi # φ_n uniformly disriuted random number array
    for n in range(N):
        k    = 1   
        xin  = 0                  # ξ_n noise
        etan = 0                  # ƞ_n noise derivative of ξ_n
        while k <= (N/2):
            w     = (2 * math.pi * x[n] * k / N) + Phi[k-1]
            xin  += k**(-beta) * math.cos(w)
            etan += -k**(1-beta) * math.sin(w) 
            k    += 1
        xi[n]  = xin                      # Array of ξ_n noise
        eta[n] = (4 * math.pi / N) * etan # Array of ƞ_n noise
    return eta,xi
    
#   This function gets x and length argument
def autocorr(x,length):
    return np.array([1]+[np.corrcoef(x[:-i], x[i:])[0,1]  \
        for i in range(1, length)])

#main part of code
x1 = np.linspace(0, 1.9, 20)     # Array of discrete points
x2 = np.linspace(2, 9, 8)        # Array of discrete points
x3 = np.linspace(10, 11.9, 20)   
# Array of discrete points
x4 = np.linspace(12, 58, 47)     # Array of discrete points
x5 = np.linspace(59, 60, 11)     # Array of discrete points
x = np.append(x1, x2)            #sum of two arry
x = np.append(x, x3)             #sum of two arry
x = np.append(x, x4)             #sum of two arry
x = np.append(x, x5)             #sum of two arry
xi_list, eta_list, ac_xi_list, ac_eta_list = ([] for i in range(4)) #4 empty array
for i in range(1):
    print(i)
    eta,xi = noise(x)    
    xi_list.append(xi)   # list of xi  
    eta_list.append(eta) # list of eta
#plot
fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2)
fig.suptitle('ξ and ƞ noise')

ax1.acorr(eta_list[0], usevlines=True, normed=True, maxlags=50, lw=1)
ax2.acorr(xi_list[0], usevlines=True, normed=True, maxlags=50, lw=1)   
ax3.scatter(x, eta_list[0], c='blue', s=2, zorder=2.5)
ax4.scatter(x ,xi_list[0] , c='blue', s=2, zorder=2.5)
# plot label/title
ax1.set_xlabel('lag')
ax1.set_ylabel('atc(eta)')
ax2.set_xlabel('lag')
ax2.set_ylabel('atc(xi)')
ax3.set_xlabel('x')
ax3.set_ylabel('eta')
ax4.set_xlabel('x')
ax4.set_ylabel('xi')

ymin = eta.min()
ymax = eta.max()

ax3.set_ylim(ymin-2, ymax+2)
ax4.set_ylim(ymin-2, ymax+2)
ax3.grid(True)
ax4.grid(True) 
plt.show()

