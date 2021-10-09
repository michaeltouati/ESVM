#######################################################################
##                                                                   ##
##             ElectroStatic Vlasov-Maxwell (ESVM) code              ##
##                                                                   ##
## Copyright © 2015 Michaël J TOUATI                                 ##
##                                                                   ##
## This file is part of ESVM.                                        ##
##                                                                   ##
## ESVM is free software: you can redistribute it and/or modify      ##
## it under the terms of the GNU General Public License as published ##
## by the Free Software Foundation, either version 3 of the License, ##
## or (at your option) any later version.                            ##
##                                                                   ##
## ESVM is distributed in the hope that it will be useful,           ##
## but WITHOUT ANY WARRANTY; without even the implied warranty of    ##
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the     ##
## GNU General Public License for more details.                      ##
##                                                                   ##
## You should have received a copy of the GNU General Public License ##
## along with ESVM. If not, see <https://www.gnu.org/licenses/>.     ##
##                                                                   ##
#######################################################################
## Initial commit written by Michaël J TOUATI - Dec. 2015
import os
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from matplotlib import mlab, cm
import library as lib

####################
# input parameters #
####################

font = {'family': 'non-serif',
        'style':  'normal',
        'color':  'black',
        'weight': 'normal',
        'size': 16,
        }

# cmap_fe_log = 'gist_ncar'
cmap_fe_log = 'nipy_spectral'
#cmap_fe_log = 'gist_stern'
#cmap_fe_log = 'gnuplot'
#cmap_fe_log = 'gnuplot2'
#cmap_fe_log = 'Blues'
#cmap_fe_log = 'jet'

##########
# Script #
##########

dir= os.path.dirname("figures/")
if not os.path.exists(dir):
    os.mkdir(dir)

print(' --------------------------------------------------------------------')
print(' 1D1V plasma electron distribution function phase-space in log. scale')
print(' --------------------------------------------------------------------')
print('  ')

print(' Search for the number of phase-space cells:')
[Nx,Nvx] = lib.search_Nx_Nvx('results/fe.dat')
print(' * found Nx  = '+str(Nx) +' space cells')
print(' * found Nvx = '+str(Nvx)+' velocity cells')
print('  ')

print(' Density plot at :')
filename = 'results/fe.dat'
N1 = Nvx
N2 = Nx

title = r'$\log_{10}{\displaystyle \left ( f_e \left (x,\,v_x,\,t\right)\,(n_0/v_{T_{e_0}}) \right )}$'

subdir= "figures/fe_log"
os.path.dirname(subdir)
if not os.path.exists(subdir):
      os.mkdir(subdir)
name='figures/fe_log/fe_log_'

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
        print(' * t (/omega_p) = '+str(time))
        N = int(counter / N3)
        for i in range(0,N1):
            for k in range(0,N2):
                X[i][k]            = x[(N-1)*N3+i*N2+k]
                VX[i][k]           = vx[(N-1)*N3+i*N2+k]
                P[i][k]            = np.log(p[(N-1)*N3+i*N2+k])/np.log(10.)
        cmap = plt.get_cmap(cmap_fe_log)
        Maxval = 1+int(np.amax(np.log(p[(int(N)-1)*int(N3):int(N)*int(N3)]))/np.log(10.))
        Minval = -15
        if Maxval <= Minval :
            Maxval = -14
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
        fig.savefig(name+str(N)+'.png',bbox_inches='tight')
        plt.close(fig)
file.close()
