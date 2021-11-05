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
Read and plot data t (fs) | x (Debye) | vx (vTe0)  | fe (n0/vTe0)
in log scale from the file fe.dat
"""
import numpy as np
import library as lib

SIMU_NAME=lib.get_results_dir()

print(' --------------------------------------------------------------------')
print(' 1D1V plasma electron distribution function phase-space in log. scale')
print(' --------------------------------------------------------------------')
print('  ')

LOGFE_RES = 'results/'+SIMU_NAME+'/fe.dat'

lib.create_dir('figures/')

SIMU_DIR = 'figures/'+SIMU_NAME+'/'
lib.create_dir(SIMU_DIR)

[N_X,N_VX] = lib.search_nx_nvx(LOGFE_RES)

print(' Density plot at :')

LOGFE_DIR = SIMU_DIR+'logfe/'
lib.create_dir(LOGFE_DIR)

LOGFE_FIG  = LOGFE_DIR + 'logfe_'
LOGFE_CMP  = 'nipy_spectral'
LOGFE_TTL  = r'$\log_{10} \displaystyle ( $'
LOGFE_TTL += r'$f_e \left (x,\,v_x,\,t\right)\,(n_0/v_{T_{e_0}}) )$'
X_LBL      = r'$x\,(\lambda_\mathrm{Debye})$'
VX_LBL     = r'$v_x\,(v_{T_{e_0}})$'
LOGFE_MIN  = -15

N_P       = N_VX*N_X
VX_PLT    = []
X_PLT     = []
LOGFE_PLT = []
X_MAP     = np.zeros((N_VX,N_X))
VX_MAP    = np.zeros((N_VX,N_X))
LOGFE_MAP = np.zeros((N_VX,N_X))
with open(LOGFE_RES, 'r', encoding='utf-8') as FILE :
    COUNTER_LOG = 0
    for LINE in FILE:
        ARRAY = LINE.strip().split()
        VX_PLT.append(float(ARRAY[1]))
        X_PLT.append( float(ARRAY[2]))
        LOGFE_PLT.append(float(ARRAY[3]))
        COUNTER_LOG += 1
        if COUNTER_LOG % N_P == 0:
            STR_TIME = f"{float(ARRAY[0]):1.4E}"
            print(' * t (/omega_p) = '+STR_TIME)
            N_T = int(COUNTER_LOG / N_P)
            for i in range(0,N_VX):
                for k in range(0,N_X):
                    INDEX           = (N_T-1)*N_P+i*N_X+k
                    X_MAP[i][k]     =  X_PLT[INDEX]
                    VX_MAP[i][k]    = VX_PLT[INDEX]
                    LOGFE_TMP       = LOGFE_PLT[INDEX]
                    LOGFE_MAP[i][k] = np.log(LOGFE_TMP)/np.log(10.)
            LOGFE_MAX = 1 + int(np.amax(np.log(LOGFE_PLT[(N_T-1)*N_P:N_T*N_P]))/np.log(10.))
            if LOGFE_MAX <= LOGFE_MIN :
                LOGFE_MAX = -14
            LOGFE_TTL2 = LOGFE_TTL+' at '+r'$t =$'+STR_TIME+r'$\,\omega_p^{-1}$'
            LOGFE_FIG2 = LOGFE_FIG+str(N_T)+'.png'
            P_MIN   = (N_T-1)*N_P
            P_MAX   = N_T*N_P
            X_MIN   = np.amin( X_PLT[P_MIN:P_MAX])
            X_MAX   = np.amax( X_PLT[P_MIN:P_MAX])
            VX_MIN  = np.amin(VX_PLT[P_MIN:P_MAX])
            VX_MAX  = np.amax(VX_PLT[P_MIN:P_MAX])
            lib.make_2d_field_pcolormesh_figure(xmap     =  X_MAP,
                                                ymap     = VX_MAP,
                                                zmap     = LOGFE_MAP,
                                                colormap = LOGFE_CMP,
                                                xmap_min =  X_MIN,
                                                xmap_max =  X_MAX,
                                                ymap_min = VX_MIN,
                                                ymap_max = VX_MAX,
                                                zmap_min = LOGFE_MIN,
                                                zmap_max = LOGFE_MAX,
                                                xlabel   =  X_LBL,
                                                ylabel   = VX_LBL,
                                                title    = LOGFE_TTL2,
                                                eq_axis  = False,
                                                fig_file = LOGFE_FIG2)
