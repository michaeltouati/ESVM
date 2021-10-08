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

cmap_Ex  = 'jet'
cmap_ne  = 'jet'
cmap_ve  = 'jet'
cmap_je  = 'jet'
cmap_vte = 'jet'

##########
# Script #
##########

dir= os.path.dirname("figures/")
if not os.path.exists(dir):
    os.mkdir(dir)

print('-------------------------')
print('1D hydrodynamic moments :')
print('-------------------------')
print(' ')

print('Search for the number of phase-space cells:')
[Nx,Nvx] = lib.search_Nx_Nvx('results/fe.dat')
print('* found Nx  = '+str(Nx) +' space cells')
print('* found Nvx = '+str(Nvx)+' velocity cells')
print(' ')

print('2D space-time density plots:')
print('* 1D elecrostatic field')
filename = 'results/Ex.dat'
name     ='figures/Ex'
title    = r'$E_x \left ( x,\,t\right )\,(m_e \omega_p v_{T_{e_0}} / e)$'
lib.plot_1D_hydro_quantity_2Dmap(Nx, filename, font, cmap_Ex, title, name)

print('* 1D elecrostatic potential')
filename = 'results/phi.dat'
name     ='figures/Phi'
title    = r'$\Phi \left ( x,\,t\right )\,(m_e {v_{T_{e_0}}}^2 / e)$'
lib.plot_1D_hydro_quantity_2Dmap(Nx, filename, font, cmap_Ex, title, name)

print('* 1D plasma electron density')
filename = 'results/ne.dat'
title    = r'$n_e \left ( x,\,t\right )\,(n_0)$'
name     ='figures/ne'
lib.plot_1D_hydro_quantity_2Dmap(Nx, filename, font, cmap_ne, title, name)

print('* 1D plasma electron mean velocity')
filename = 'results/ve.dat'
name     ='figures/ve'
cmap     = 'jet'
title    = r'$v_e \left ( x,\,t\right )\,(v_{T_{e_0}})$'
lib.plot_1D_hydro_quantity_2Dmap(Nx, filename, font, cmap_ve, title, name)

filename = 'results/je.dat'
cmap     = 'jet'
title    = r'$j_e \left ( x,\,t\right )\,(n_0 e v_{T_{e_0}})$'
name     ='figures/je'
lib.plot_1D_hydro_quantity_2Dmap(Nx, filename, font, cmap_je, title, name)

print('* 1D plasma electron thermal velocity (standard deviation)')
filename = 'results/vTe.dat'
cmap     = 'jet'
name     ='figures/vte'
title    = r'$v_{T_e} \left ( x,\,t\right )\,(v_{T_{e_0}})$'
lib.plot_1D_hydro_quantity_2Dmap(Nx, filename, font, cmap_vte, title, name)
print(' ')