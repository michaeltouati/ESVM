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
"""
Read and plot data t (fs) | x (microns) | F from files:
* Ex.dat
* Phi.dat
* ne.dat
* je.dat
* ve.dat
* vTe.dat
"""
import library as lib

SIMU_NAME=lib.get_results_dir()

print(' -------------------------')
print(' 1D hydrodynamic moments :')
print(' -------------------------')
print('  ')

lib.create_dir('figures/')

SIMU_DIR = 'figures/'+SIMU_NAME+'/'
lib.create_dir(SIMU_DIR)

RES_DIR = 'results/'+SIMU_NAME+'/'

print(' Search for the number of phase-space cells:')
[N_X,N_VX] = lib.search_nx_nvx(RES_DIR+'fe.dat')
print(' * found Nx  = '+str(N_X) +' space bins')
print(' * found Nvx = '+str(N_VX)+' velocity bins')
print('  ')

################
# Scalar plots #
################

print(' Scalar plot :')
print(' * 1D Electrostatic field at :')
EX_RES = RES_DIR+'Ex.dat'
EX_LBL = r'$E_x \left ( x,\,t\right ) \, \left( m_e \omega_p v_{T_{e_0}} / e \right )$'
EX_DIR = SIMU_DIR+'Ex/'
lib.create_dir(EX_DIR)
EX_FIG = EX_DIR+'Ex_'
lib.plot_1d_hydro_quantity_scalar_plot(N_X, EX_RES, EX_LBL, EX_FIG)

print(' * 1D Electrostatic potential')
PHI_RES = RES_DIR+'Phi.dat'
PHI_LBL = r'$\Phi \left ( x,\,t\right ) \, \left( m_e {v_{T_{e_0}}}^2 / e \right )$'
PHI_DIR = SIMU_DIR+'Phi/'
lib.create_dir(PHI_DIR)
PHI_FIG = PHI_DIR+'Phi_'
lib.plot_1d_hydro_quantity_scalar_plot(N_X, PHI_RES, PHI_LBL, PHI_FIG)

print(' * 1D plasma electron density')
NE_RES = RES_DIR+'ne.dat'
NE_LBL = r'$n_e \left ( x,\,t\right ) \, \left( n_0 \right )$'
NE_DIR = SIMU_DIR + 'ne/'
lib.create_dir(NE_DIR)
NE_FIG = NE_DIR + 'ne_'
lib.plot_1d_hydro_quantity_scalar_plot(N_X, NE_RES, NE_LBL, NE_FIG)

print(' * 1D plasma electron mean velocity')
VE_RES = RES_DIR+'ve.dat'
VE_LBL = r'$v_e \left ( x,\,t\right ) \, \left( v_{T_{e_0}}\right )$'
VE_DIR = SIMU_DIR+'ve/'
lib.create_dir(VE_DIR)
VE_FIG = VE_DIR + 've_'
lib.plot_1d_hydro_quantity_scalar_plot(N_X, VE_RES, VE_LBL, VE_FIG)

print(' * 1D plasma electron current density')
JE_RES = RES_DIR+'je.dat'
JE_LBL = r'$j_e \left ( x,\,t\right ) \, \left( n_0 e v_{T_{e_0}}\right )$'
JE_DIR = SIMU_DIR+'je/'
lib.create_dir(JE_DIR)
JE_FIG = JE_DIR + 'je_'
lib.plot_1d_hydro_quantity_scalar_plot(N_X, JE_RES, JE_LBL, JE_FIG)

print(' * 1D plasma electron thermal velocity')
print('   (standard deviation)')
VTE_RES = RES_DIR+'vTe.dat'
VTE_LBL = r'$v_{T_e} \left ( x,\,t\right ) \, \left( v_{T_{e_0}}\right )$'
VTE_DIR = SIMU_DIR + 'vTe/'
lib.create_dir(VTE_DIR)
VTE_FIG = VTE_DIR + 'vTe_'
lib.plot_1d_hydro_quantity_scalar_plot(N_X, VTE_RES, VTE_LBL, VTE_FIG)
