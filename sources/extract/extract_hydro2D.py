# -----------------------------------------------
# SIMULATION RESULTS PLOTTING FOR 1D1V ESVM SIMU 
# -----------------------------------------------
# Initial commit written by Michaël J Touati

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

cmap_Ex  = 'jet'
cmap_ne  = 'jet'
cmap_ve  = 'jet'
cmap_je  = 'jet'
cmap_vte = 'jet'

#############################
# Create figures/ directory #
#############################

dir= os.path.dirname("figures/")
if not os.path.exists(dir):
    os.mkdir(dir)

#################################
# functions used in this script #
#################################

def plot_1D_hydro_quantity_2Dmap(N_x, file_name, font_properties, colormap, plot_title,plot_name):
    t = []
    x = []
    p = []
    file = open(file_name, 'r')
    counter = 0
    for line in file:
        line      = line.strip()
        array     = line.split()
        t.append(float(array[0]))
        x.append(float(array[1]))
        p.append(float(array[2]))
        counter += 1
        if counter % N_x == 0:
          N_t = int(counter / N_x)
    file.close()
    T = np.zeros((N_x,N_t))
    X = np.zeros((N_x,N_t))
    P = np.zeros((N_x,N_t))
    for i in range(0,N_x):
        for n in range(0,N_t):
            T[i][n] = t[n*Nx+i]
            X[i][n] = x[n*Nx+i]
            P[i][n] = p[n*Nx+i]
    cmap = plt.get_cmap(colormap)
    Maxval = np.amax(p)
    Minval = np.amin(p)
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
    plt.pcolormesh(X,T,P,cmap=cmap,norm=norm,shading='gouraud')
    cbar=plt.colorbar()
    cbar.ax.tick_params(labelsize=16)
    plt.title(plot_title, fontdict=font_properties)
    plt.xticks(fontsize=16)
    plt.xlabel(r'$x\,(\lambda_\mathrm{Debye})$', fontdict=font_properties)
    plt.xlim([np.amin(x),np.amax(x)])
    plt.ylabel(r'$t\,(/ \omega_p)$', fontdict=font)
    plt.yticks(fontsize=16)
    plt.ylim([np.amin(t),np.amax(t)])
    fig.savefig(name+'.png',bbox_inches='tight')
    plt.close(fig)

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

################################
# 2D space-time colormap plots #
################################

print('2D space-time colormap plots:')
print('* 1D elecrostatic field')
filename = 'results/Ex.dat'
name     ='figures/Ex'
title    = r'$E_x \left ( x,\,t\right )\,(m_e v_{T_{e_0}} \omega_p / e)$'
plot_1D_hydro_quantity_2Dmap(Nx, filename, font, cmap_Ex, title, name)

print('* 1D plasma electron density')
filename = 'results/ne.dat'
title    = r'$n_e \left ( x,\,t\right )\,(n_0)$'
name     ='figures/ne'
plot_1D_hydro_quantity_2Dmap(Nx, filename, font, cmap_ne, title, name)

print('* 1D plasma electron mean velocity')
filename = 'results/ve.dat'
name     ='figures/ve'
cmap     = 'jet'
title    = r'$v_e \left ( x,\,t\right )\,(v_{T_{e_0}})$'
plot_1D_hydro_quantity_2Dmap(Nx, filename, font, cmap_ve, title, name)

filename = 'results/je.dat'
cmap     = 'jet'
title    = r'$j_e \left ( x,\,t\right )\,(n_0 e v_{T_{e_0}})$'
name     ='figures/je'
plot_1D_hydro_quantity_2Dmap(Nx, filename, font, cmap_je, title, name)

print('* 1D plasma electron thermal velocity (standard deviation)')
filename = 'results/vTe.dat'
cmap     = 'jet'
name     ='figures/vte'
title    = r'$v_{T_e} \left ( x,\,t\right )\,(v_{T_{e_0}})$'
plot_1D_hydro_quantity_2Dmap(Nx, filename, font, cmap_vte, title, name)
