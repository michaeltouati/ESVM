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
print(' * found Nx  = '+str(Nx) +' space bins')
print(' * found Nvx = '+str(Nvx)+' velocity bins')
print('  ')

################
# Scalar plots #
################

print(' Scalar plot :')
print(' * 1D Electrostatic field at :')
filename = 'results/'+simu_name+'/Ex.dat'
Ylabel   = r'$E_x \left ( x,\,t\right ) \, \left( m_e \omega_p v_{T_{e_0}} / e \right )$'
lib.create_dir('figures/'+simu_name+'/Ex/')
name     ='figures/'+simu_name+'/Ex/Ex_'
lib.plot_1D_hydro_quantity_scalar_plot(Nx, filename, font, Ylabel, name)

print(' * 1D Electrostatic potential')
filename = 'results/'+simu_name+'/Phi.dat'
Ylabel   = r'$\Phi \left ( x,\,t\right ) \, \left( m_e {v_{T_{e_0}}}^2 / e \right )$'
lib.create_dir('figures/'+simu_name+'/Phi/')
name     ='figures/'+simu_name+'/Phi/Phi_'
lib.plot_1D_hydro_quantity_scalar_plot(Nx, filename, font, Ylabel, name)

print(' * 1D plasma electron density')
filename = 'results/'+simu_name+'/ne.dat'
Ylabel   = r'$n_e \left ( x,\,t\right ) \, \left( n_0 \right )$'
lib.create_dir('figures/'+simu_name+'/ne/')
name     ='figures/'+simu_name+'/ne/ne_'
lib.plot_1D_hydro_quantity_scalar_plot(Nx, filename, font, Ylabel, name)

print(' * 1D plasma electron mean velocity')
filename = 'results/'+simu_name+'/ve.dat'
Ylabel   = r'$v_e \left ( x,\,t\right ) \, \left( v_{T_{e_0}}\right )$'
lib.create_dir('figures/'+simu_name+'/ve/')
name    ='figures/'+simu_name+'/ve/ve_'
lib.plot_1D_hydro_quantity_scalar_plot(Nx, filename, font, Ylabel, name)

print(' * 1D plasma electron current density')
filename = 'results/'+simu_name+'/je.dat'
Ylabel   = r'$j_e \left ( x,\,t\right ) \, \left( n_0 e v_{T_{e_0}}\right )$'
lib.create_dir('figures/'+simu_name+'/je/')
name     ='figures/'+simu_name+'/je/je_'
lib.plot_1D_hydro_quantity_scalar_plot(Nx, filename, font, Ylabel, name)

print(' * 1D plasma electron thermal velocity (standard deviation)')
filename = 'results/'+simu_name+'/vTe.dat'
Ylabel   = r'$v_{T_e} \left ( x,\,t\right ) \, \left( v_{T_{e_0}}\right )$'
lib.create_dir('figures/'+simu_name+'/vTe/')
name     ='figures/'+simu_name+'/vTe/vTe_'
lib.plot_1D_hydro_quantity_scalar_plot(Nx, filename, font, Ylabel, name)
