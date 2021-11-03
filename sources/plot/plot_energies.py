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
Read and plot data t (fs) | Energy (J) from files:
* UK.dat
* UT.dat
* UE.dat
"""
import library as lib

SIMU_NAME=lib.get_results_dir()

print(' ---------------------------')
print(' Total energies scalar plots')
print(' ---------------------------')
print('  ')

lib.create_dir('figures/')

SIMU_DIR = 'figures/'+SIMU_NAME+'/'
lib.create_dir(SIMU_DIR)

RES_DIR = 'results/'+SIMU_NAME+'/'

# total plasma electron kinetic energy
[TK,UK0] = lib.extract_energy_file(RES_DIR+'UK.dat')
T   = TK
N_T = len(TK)
# total plasma electron internal energy
[TT,UT0] = lib.extract_energy_file(RES_DIR+'UT.dat')
if len(TT) < N_T:
    T   = TT
    N_T = len(TT)
# total electrostatic energy
[TE,UE0] = lib.extract_energy_file(RES_DIR+'UE.dat')
if len(TE) < N_T:
    T   = TE
    N_T = len(TE)
#
UK   = []
UT   = []
UE   = []
UTOT = []
for i in range(0,N_T):
    UK.append(UK0[i])
    UT.append(UT0[i])
    UE.append(UE0[i])
    UTOT.append(UE[i]+UT[i]+UK[i])
T_MIN    = 0.
T_MAX    = T[N_T-1]
T_LBL    = r'$t\,(\omega_p^{-1})$'
U_LBL    = r'$\mathrm{Energy}\,$'
U_LBL   += r'$\left (n_0 {\lambda_{\mathrm{Debye}}}^3 m_e {v_{T_{e_0}}}^2 \right )$'
U_FIGLOG = SIMU_DIR+'energy_log.png'
lib.make_scalars_plot_figure(xplot     = T,
                             yplot1    = UK,
                             legend1   = r'$U_{K_e}$',
                             color1    = 'green',
                             yplot2    = UT,
                             legend2   = r'$U_{T_e}$',
                             color2    = 'red',
                             yplot3    = UE,
                             legend3   = r'$U_{E_x}$',
                             color3    = 'blue',
                             yplot4    = UTOT,
                             legend4   = r'$U_\mathrm{tot}$',
                             color4    = 'black',
                             xplot_min = T_MIN,
                             xplot_max = T_MAX,
                             xlabel    = T_LBL,
                             ylabel    = U_LBL,
                             filename  = U_FIGLOG,
                             logx      = False,
                             logy      = True,
                             grid      = False)

U_FIG = SIMU_DIR+'energy.png'
lib.make_scalars_plot_figure(xplot     = T,
                             yplot1    = UK,
                             legend1   = r'$U_{K_e}$',
                             color1    = 'green',
                             yplot3    = UE,
                             legend3   = r'$U_{E_x}$',
                             color3    = 'blue',
                             xplot_min = T_MIN,
                             xplot_max = T_MAX,
                             xlabel    = T_LBL,
                             ylabel    = U_LBL,
                             filename  = U_FIG,
                             logx      = False,
                             logy      = False,
                             grid      = False)
