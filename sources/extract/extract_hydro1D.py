# ------------------------------------------------------
# SIMULATION RESULTS PLOTTING FOR 1D1V ESVM SIMULATIONS
# ------------------------------------------------------
# Initial commit written by Michaël J Touati - Dec. 2015

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

#############################
# Create figures/ directory #
#############################

dir= os.path.dirname("figures/")
if not os.path.exists(dir):
    os.mkdir(dir)

#################################
# functions used in this script #
#################################

def plot_1D_hydro_quantity_scalar_plot(N_x, file_name, font_properties, y_label, plot_name):
    x = []
    p = []
    X = np.zeros(N_x)
    P = np.zeros(N_x)
    file = open(file_name, 'r')
    counter = 0
    for line in file:
        line      = line.strip()
        array     = line.split()
        x.append(float(array[1]))
        p.append(float(array[2]))
        counter = counter + 1
        if counter % N_x == 0:
            N = int(counter / N_x)
            for i in range(0,N_x):
                X[i] = x[(N-1)*N_x+i]
                P[i] = p[(N-1)*N_x+i]
            time = int(100.*float(array[0]))/100.
            print('  t (/omega_p) = '+str(time))
            fig=plt.figure()
            plt.rc('text', usetex=True)
            plt.plot(X, P, linewidth=2, color = 'black')
            plt.title(r'$t =$'+str(time)+r'$\,\omega_p^{-1}$', fontdict=font_properties)
            plt.xticks(fontsize=16)
            plt.xlabel(r'$x\,(\lambda_\mathrm{Debye})$', fontdict=font)
            plt.xlim([np.amin(X),np.amax(X)])
            plt.ylabel(y_label, fontdict=font)
            plt.yticks(fontsize=16)
            ex_min = np.amin(p)
            ex_max = np.amax(p)
            if ( ex_min == ex_max ):
                if (ex_min == 0.) :
                    ex_min = -1.
                    ex_max =  1.
                else :
                    ex_min = 0.9 * ex_max
                    ex_max = 1.1 * ex_max
            plt.ylim([ex_min,ex_max])
            fig.savefig(name+str(N)+'.png',bbox_inches='tight')
            plt.close(fig)
    file.close()

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

################
# Scalar plots #
################

print('1D scalar plots :')
print('* 1D Electrostatic field')
filename = 'results/Ex.dat'
Ylabel   = r'$E_x \left ( x,\,t\right ) \, \left( m_e \omega_p v_{T_{e_0}} / e \right )$'
subdir   = "figures/Ex"
os.path.dirname(subdir)
if not os.path.exists(subdir):
      os.mkdir(subdir)
name     ='figures/Ex/Ex_'
plot_1D_hydro_quantity_scalar_plot(Nx, filename, font, Ylabel, name)

print('* 1D plasma electron density')
filename = 'results/ne.dat'
Ylabel   = r'$n_e \left ( x,\,t\right ) \, \left( n_0 \right )$'
subdir   = "figures/ne"
os.path.dirname(subdir)
if not os.path.exists(subdir):
      os.mkdir(subdir)
name     ='figures/ne/ne_'
plot_1D_hydro_quantity_scalar_plot(Nx, filename, font, Ylabel, name)

print('* 1D plasma electron mean velocity')
filename = 'results/ve.dat'
Ylabel   = r'$v_e \left ( x,\,t\right ) \, \left( v_{T_{e_0}}\right )$'
subdir   = "figures/ve"
os.path.dirname(subdir)
if not os.path.exists(subdir):
      os.mkdir(subdir)
name    ='figures/ve/ve_'
plot_1D_hydro_quantity_scalar_plot(Nx, filename, font, Ylabel, name)

print('* 1D plasma electron current density')
filename = 'results/je.dat'
Ylabel   = r'$j_e \left ( x,\,t\right ) \, \left( n_0 e v_{T_{e_0}}\right )$'
subdir   = "figures/je"
os.path.dirname(subdir)
if not os.path.exists(subdir):
      os.mkdir(subdir)
name     ='figures/je/je_'
plot_1D_hydro_quantity_scalar_plot(Nx, filename, font, Ylabel, name)

print('* 1D plasma electron thermal velocity (standard deviation)')
filename = 'results/vTe.dat'
Ylabel   = r'$v_{T_e} \left ( x,\,t\right ) \, \left( v_{T_{e_0}}\right )$'
subdir   = "figures/vte"
os.path.dirname(subdir)
if not os.path.exists(subdir):
      os.mkdir(subdir)
name     ='figures/vte/vte_'
plot_1D_hydro_quantity_scalar_plot(Nx, filename, font, Ylabel, name)
