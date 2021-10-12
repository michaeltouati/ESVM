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

simu_name=lib.get_results_dir()

lib.create_dir('figures/')
lib.create_dir('figures/'+simu_name+'/')

print(' -------------------------')
print(' 1D hydrodynamic moments :')
print(' -------------------------')
print('  ')

print(' Search for the number of phase-space cells:')
[Nx,Nvx] = lib.search_Nx_Nvx('results/'+simu_name+'/fe.dat')
print(' * found Nx  = '+str(Nx) +' space cells')
print(' * found Nvx = '+str(Nvx)+' velocity cells')
print('  ')

print(' 2D space-time density plot:')
print(' * 1D elecrostatic field')
filename = 'results/'+simu_name+'/Ex.dat'
name     ='figures/'+simu_name+'/Ex'
title    = r'$E_x \left ( x,\,t\right )\,(m_e \omega_p v_{T_{e_0}} / e)$'
lib.plot_1D_hydro_quantity_2Dmap(Nx, filename, font, cmap_Ex, title, name)

print(' * 1D elecrostatic potential')
filename = 'results/'+simu_name+'/Phi.dat'
name     ='figures/'+simu_name+'/Phi'
title    = r'$\Phi \left ( x,\,t\right )\,(m_e {v_{T_{e_0}}}^2 / e)$'
lib.plot_1D_hydro_quantity_2Dmap(Nx, filename, font, cmap_Ex, title, name)

print(' * 1D plasma electron density')
filename = 'results/'+simu_name+'/ne.dat'
title    = r'$n_e \left ( x,\,t\right )\,(n_0)$'
name     ='figures/'+simu_name+'/ne'
lib.plot_1D_hydro_quantity_2Dmap(Nx, filename, font, cmap_ne, title, name)

print(' * 1D plasma electron mean velocity')
filename = 'results/'+simu_name+'/ve.dat'
name     ='figures/'+simu_name+'/ve'
cmap     = 'jet'
title    = r'$v_e \left ( x,\,t\right )\,(v_{T_{e_0}})$'
lib.plot_1D_hydro_quantity_2Dmap(Nx, filename, font, cmap_ve, title, name)

filename = 'results/'+simu_name+'/je.dat'
cmap     = 'jet'
title    = r'$j_e \left ( x,\,t\right )\,(n_0 e v_{T_{e_0}})$'
name     ='figures/'+simu_name+'/je'
lib.plot_1D_hydro_quantity_2Dmap(Nx, filename, font, cmap_je, title, name)

print(' * 1D plasma electron thermal velocity (standard deviation)')
filename = 'results/'+simu_name+'/vTe.dat'
cmap     = 'jet'
name     ='figures/'+simu_name+'/vTe'
title    = r'$v_{T_e} \left ( x,\,t\right )\,(v_{T_{e_0}})$'
lib.plot_1D_hydro_quantity_2Dmap(Nx, filename, font, cmap_vte, title, name)
print('  ')
