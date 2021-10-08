# -----------------------------------------------
# SIMULATION RESULTS PLOTTING FOR 1D1V ESVM SIMU  
# -----------------------------------------------
# Initial commit written by Michael J Touati

import os
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from matplotlib import mlab, cm

####################
# Input parameters #
####################

font = {'family': 'non-serif',
        'style':  'normal',
        'color':  'black',
        'weight': 'normal',
        'size': 16,
        }

# cmap_fe = 'gist_ncar'
# cmap_fe = 'nipy_spectral'
#cmap_fe = 'gist_stern'
#cmap_fe = 'gnuplot'
#cmap_fe = 'gnuplot2'
cmap_fe = 'Blues'
#cmap_fe = 'jet'

#############################
# Create figures/ directory #
#############################

dir= os.path.dirname("figures/")
if not os.path.exists(dir):
    os.mkdir(dir)

##########################################
# Find the number of spatial grid points # 
# (Nvx,Nx) by reading 'results/fe.dat'   #
##########################################

print('Search for the number of phase-space cells:')
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
        NvxNx = counter
        break
file.close()
Nvx = int(NvxNx / Nx)

print('* found Nx  = '+str(Nx) +' space cells')
print('* found Nvx = '+str(Nvx)+' velocity cells')
print(' ')

#######################################################
# 1D1V plasma electron distribution function colormap #
#######################################################

print('1D1V plasma electron distribution function phase-space colormap plot')
filename = 'results/fe.dat'
N1 = Nvx
N2 = Nx
title = r'$f_e \left ( x,\,v_x,\,t\right )\,(n_0/v_{T_{e_0}})$'
subdir= "figures/fe"
os.path.dirname(subdir)
if not os.path.exists(subdir):
      os.mkdir(subdir)
name='figures/fe/fe_'

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
        time = int(100.*float(array[0]))/100.
        print(' t (/omega_p) = '+str(time))
        N = int(counter / N3)
        val=max(abs(np.amax(p[(N-1)*N3:N*N3])),abs(np.amin(p[(N-1)*N3:N*N3])))
        if (val != 0):
            Power=int(np.log(val)/np.log(10))
        else:
            Power = 0
        p[(N-1)*N3:N*N3] = [pp / (10.**Power) for pp in p[(N-1)*N3:N*N3]]
        for i in range(0,N1):
            for k in range(0,N2):
                X[i][k]=x[(N-1)*N3+i*N2+k]
                VX[i][k]=vx[(N-1)*N3+i*N2+k]
                P[i][k]=p[(N-1)*N3+i*N2+k]
        cmap = plt.get_cmap(cmap_fe)
        Maxval = np.amax(p[(int(N)-1)*int(N3):int(N)*int(N3)])
        Minval = np.amin(p[(int(N)-1)*int(N3):int(N)*int(N3)])
        if ( Minval == Maxval ) :
            if ( Minval == 0. ) :
                Maxval =  1.
                Minval = -1.
            else :
                Maxval = 1.1*Maxval
                Minval = 0.9*Maxval
        norm = cm.colors.Normalize(vmax=Maxval, vmin=Minval)
        fig=plt.figure()
        plt.rc('text', usetex=True)
        plt.pcolormesh(X,VX,P,cmap=cmap,norm=norm,shading='gouraud')
        cbar=plt.colorbar()
        cbar.ax.tick_params(labelsize=16)
        plt.title(title+' at '+r'$t =$'+str(time)+r'$\,\omega_p^{-1}$', fontdict=font)
        plt.xticks(fontsize=16)
        plt.xlabel(r'$x\,(\lambda_\mathrm{Debye})$', fontdict=font)
        plt.xlim([np.amin(x[(N-1)*N3:N*N3]),np.amax(x[(N-1)*N3:N*N3])])
        plt.ylabel(r'$v_x\,(v_{T_{e_0}})$', fontdict=font)
        plt.yticks(fontsize=16)
        plt.ylim([np.amin(vx[(N-1)*N3:N*N3]),np.amax(vx[(N-1)*N3:N*N3])])
        #plt.axes().set_aspect('equal')
        fig.savefig(name+str(N)+'.png',bbox_inches='tight')
        plt.close(fig)
file.close()
