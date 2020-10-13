# -------------------------------------------------
# SIMULATION RESULTS PLOTTING FOR 1D1V ESV SIMU OF 
# ELECTROSTATIC WAKEFIELD EMMISSION BY AN ELECTRON 
# -------------------------------------------------
# Written by Michael J Touati - CLPU - 2020/07/15

import os
import numpy as np
import matplotlib
import math
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from matplotlib import mlab, cm

dir= os.path.dirname("figures/")
if not os.path.exists(dir):
    os.mkdir(dir)

subdir= "results"
os.path.dirname(subdir)
if not os.path.exists(subdir):
      os.mkdir(subdir)


# find the number of spatial grid points (Nvx,Nx) by reading 'results.dat'
t0  = []
v0  = []
file = open('results/fe.dat', 'r')
line = file.readline()
line = line.strip()
array = line.split()
t0.append(float(array[0]))
v0.append(float(array[1]))
counter = 0
for line in file:
    line      = line.strip()
    array     = line.split()
    t0.append(float(array[0]))
    v0.append(float(array[1]))
    counter = counter + 1
    if v0[counter]!=v0[counter-1]:
        Nx = counter
        break
for line in file:
    line      = line.strip()
    array     = line.split()
    t0.append(float(array[0]))
    counter = counter + 1
    if t0[counter]!=t0[counter-1]:
        NxvNx = counter
        break
file.close()
Nxv = NxvNx / Nx

font = {'family': 'non-serif',
        'style':  'normal',
        'color':  'black',
        'weight': 'normal',
        'size': 16,
        }

print('Distribution function plot')
filename = 'results/fe.dat'
N1 = Nxv
N2 = Nx
# cmap = 'gist_ncar'
cmap = 'nipy_spectral'
#cmap = 'gist_stern'
#cmap = 'gnuplot'
#cmap = 'gnuplot2'
#cmap = 'Blues'
#cmap = 'jet'
title = r'$\log_{10}{\displaystyle \left ( \mathbf{f} \left (x,\,v_x,\,t\right)\,(n_0.{v_T}^{-3}) \right )}$'
subdir= "figures/f_log"
os.path.dirname(subdir)
if not os.path.exists(subdir):
      os.mkdir(subdir)
name='figures/f_log/f_log_'

N3 = int(N1*N2)
vx = []
x = []
p = []
N1 = int(N1)
N2 = int(N2)
N3 = int(N3)
X = np.zeros((N1,N2))
VX = np.zeros((N1,N2))
P = np.zeros((N1,N2))
file = open(filename, 'r')
counter = 0
for line in file:
    line      = line.strip()
    array     = line.split()
    vx.append(float(array[1]))
    x.append(float(array[2]))
    p.append(float(array[3]))
    counter = counter + 1
    if counter % N3 == 0:
        time = math.floor(100.*float(array[0]))/100.
        N = int(counter / N3)
        # val=max(abs(np.amax(p[(N-1)*N3:N*N3])),abs(np.amin(p[(N-1)*N3:N*N3])))
        # if (val != 0):
        #     Power=math.floor(np.log(val)/np.log(10))
        # else:
        #     Power = 0
        # p[(N-1)*N3:N*N3] = [pp / (10.**Power) for pp in p[(N-1)*N3:N*N3]]
        for i in range(0,N1):
            for k in range(0,N2):
                X[i][k]            = x[(N-1)*N3+i*N2+k]
                VX[i][k]           = vx[(N-1)*N3+i*N2+k]
                #p[(N-1)*N3+i*N2+k] = np.log(np.max([1.e-15,p[(N-1)*N3+i*N2+k]/np.log(10.)])) # double precision
                P[i][k]            = np.log(p[(N-1)*N3+i*N2+k])/np.log(10.)
        cmap = plt.get_cmap(cmap)
        Maxval = 1+int(np.amax(np.log(p[(int(N)-1)*int(N3):int(N)*int(N3)]))/np.log(10.))
        Minval = -15. #int(np.amin(p[(int(N)-1)*int(N3):int(N)*int(N3)]))
        norm = cm.colors.Normalize(vmax=Maxval, vmin=Minval)
        fig=plt.figure()
        plt.rc('text', usetex=True)
        plt.pcolormesh(X,VX,P,cmap=cmap,norm=norm,vmax=Maxval,vmin=Minval)
        cbar=plt.colorbar()
        cbar.ax.tick_params(labelsize=16)
        plt.title(title+' at '+str(time)+r'$\,\omega_p^{-1}$', fontdict=font)
        plt.xticks(fontsize=16)
        plt.xlabel(r'$x\,(\lambda_\mathrm{Debye})$', fontdict=font)
        plt.xlim([np.amin(x[(N-1)*N3:N*N3]),np.amax(x[(N-1)*N3:N*N3])])
        plt.ylabel(r'$v_x\,(v_\mathrm{T})$', fontdict=font)
        plt.yticks(fontsize=16)
        plt.ylim([np.amin(vx[(N-1)*N3:N*N3]),np.amax(vx[(N-1)*N3:N*N3])])
        fig.savefig(name+str(N)+'.png',bbox_inches='tight')
        plt.close(fig)
file.close()