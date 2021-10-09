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
import library as lib

####################
# Input parameters #
####################

font = {'family': 'non-serif',
        'style':  'normal',
        'color':  'black',
        'weight': 'normal',
        'size': 16,
        }

##########
# Script #
##########

dir= os.path.dirname("figures/")
if not os.path.exists(dir):
    os.mkdir(dir)

print(' -------------------------')
print(' 1D hydrodynamic moments :')
print(' -------------------------')
print('  ')


print(' Search for the number of phase-space cells:')
[Nx,Nvx] = lib.search_Nx_Nvx('results/fe.dat')
print(' * found Nx  = '+str(Nx) +' space cells')
print(' * found Nvx = '+str(Nvx)+' velocity cells')
print('  ')

################
# Scalar plots #
################

print(' Scalar plot :')
print(' * 1D Electrostatic field at :')
filename = 'results/Ex.dat'
Ylabel   = r'$E_x \left ( x,\,t\right ) \, \left( m_e \omega_p v_{T_{e_0}} / e \right )$'
subdir   = "figures/Ex"
os.path.dirname(subdir)
if not os.path.exists(subdir):
      os.mkdir(subdir)
name     ='figures/Ex/Ex_'
lib.plot_1D_hydro_quantity_scalar_plot(Nx, filename, font, Ylabel, name)

print(' * 1D Electrostatic potential')
filename = 'results/phi.dat'
Ylabel   = r'$\Phi \left ( x,\,t\right ) \, \left( m_e {v_{T_{e_0}}}^2 / e \right )$'
subdir   = "figures/Phi"
os.path.dirname(subdir)
if not os.path.exists(subdir):
      os.mkdir(subdir)
name     ='figures/Phi/Phi_'
lib.plot_1D_hydro_quantity_scalar_plot(Nx, filename, font, Ylabel, name)

print(' * 1D plasma electron density')
filename = 'results/ne.dat'
Ylabel   = r'$n_e \left ( x,\,t\right ) \, \left( n_0 \right )$'
subdir   = "figures/ne"
os.path.dirname(subdir)
if not os.path.exists(subdir):
      os.mkdir(subdir)
name     ='figures/ne/ne_'
lib.plot_1D_hydro_quantity_scalar_plot(Nx, filename, font, Ylabel, name)

print(' * 1D plasma electron mean velocity')
filename = 'results/ve.dat'
Ylabel   = r'$v_e \left ( x,\,t\right ) \, \left( v_{T_{e_0}}\right )$'
subdir   = "figures/ve"
os.path.dirname(subdir)
if not os.path.exists(subdir):
      os.mkdir(subdir)
name    ='figures/ve/ve_'
lib.plot_1D_hydro_quantity_scalar_plot(Nx, filename, font, Ylabel, name)

print(' * 1D plasma electron current density')
filename = 'results/je.dat'
Ylabel   = r'$j_e \left ( x,\,t\right ) \, \left( n_0 e v_{T_{e_0}}\right )$'
subdir   = "figures/je"
os.path.dirname(subdir)
if not os.path.exists(subdir):
      os.mkdir(subdir)
name     ='figures/je/je_'
lib.plot_1D_hydro_quantity_scalar_plot(Nx, filename, font, Ylabel, name)

print(' * 1D plasma electron thermal velocity (standard deviation)')
filename = 'results/vTe.dat'
Ylabel   = r'$v_{T_e} \left ( x,\,t\right ) \, \left( v_{T_{e_0}}\right )$'
subdir   = "figures/vte"
os.path.dirname(subdir)
if not os.path.exists(subdir):
      os.mkdir(subdir)
name     ='figures/vte/vte_'
lib.plot_1D_hydro_quantity_scalar_plot(Nx, filename, font, Ylabel, name)