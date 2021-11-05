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


[SIMU_NAME, RES_DIR, SIMU_DIR, N_X, N_VX] = lib.init_plot_hydro()

print(' 2D space-time density plot:')
print(' * 1D elecrostatic field')
EX_RES = RES_DIR  + 'Ex.dat'
EX_FIG = SIMU_DIR + 'Ex'
EX_TTL = r'$E_x \left ( x,\,t\right )\,(m_e \omega_p v_{T_{e_0}} / e)$'
EX_CMP = 'jet'
lib.plot_1d_hydro_quantity_2dmap(N_X, EX_RES, EX_CMP, EX_TTL, EX_FIG)

print(' * 1D elecrostatic potential')
PHI_RES = RES_DIR  + 'Phi.dat'
PHI_FIG = SIMU_DIR + 'Phi'
PHI_TTL = r'$\Phi \left ( x,\,t\right )\,(m_e {v_{T_{e_0}}}^2 / e)$'
PHI_CMP = 'jet'
lib.plot_1d_hydro_quantity_2dmap(N_X, PHI_RES, PHI_CMP, PHI_TTL, PHI_FIG)

print(' * 1D plasma electron density')
NE_RES = RES_DIR  + 'ne.dat'
NE_FIG = SIMU_DIR + 'ne'
NE_TTL = r'$n_e \left ( x,\,t\right )\,(n_0)$'
NE_CMP = 'jet'
lib.plot_1d_hydro_quantity_2dmap(N_X, NE_RES, NE_CMP, NE_TTL, NE_FIG)

print(' * 1D plasma electron mean velocity')
VE_RES = RES_DIR  + 've.dat'
VE_FIG = SIMU_DIR + 've'
VE_TTL = r'$v_e \left ( x,\,t\right )\,(v_{T_{e_0}})$'
VE_CMP = 'jet'
lib.plot_1d_hydro_quantity_2dmap(N_X, VE_RES, VE_CMP, VE_TTL, VE_FIG)

print(' * 1D plasma electron electrical current')
JE_RES = RES_DIR  + 'je.dat'
JE_FIG = SIMU_DIR + 'je'
JE_TTL = r'$j_e \left ( x,\,t\right )\,(n_0 e v_{T_{e_0}})$'
JE_CMP = 'jet'
lib.plot_1d_hydro_quantity_2dmap(N_X, JE_RES, JE_CMP, JE_TTL, JE_FIG)

print(' * 1D plasma electron thermal velocity')
print('   (standard deviation)')
VTE_RES = RES_DIR  + 'vTe.dat'
VTE_FIG = SIMU_DIR + 'vTe'
VTE_TTL = r'$v_{T_e} \left ( x,\,t\right )\,(v_{T_{e_0}})$'
VTE_CMP = 'jet'
lib.plot_1d_hydro_quantity_2dmap(N_X, VTE_RES, VTE_CMP, VTE_TTL, VTE_FIG)
print('  ')
