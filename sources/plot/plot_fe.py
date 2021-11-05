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
from the file fe.dat
"""
import numpy as np
import library as lib

SIMU_NAME=lib.get_results_dir()

print(' ------------------------------------------------------')
print(' 1D1V plasma electron distribution function phase-space')
print(' ------------------------------------------------------')
print('  ')

FE_RES = 'results/'+SIMU_NAME+'/fe.dat'

lib.create_dir('figures/')

SIMU_DIR = 'figures/'+SIMU_NAME+'/'
lib.create_dir(SIMU_DIR)

[N_X,N_VX] = lib.search_nx_nvx(FE_RES)

print(' Density plot at :')

FE_DIR = SIMU_DIR+'fe/'
lib.create_dir(FE_DIR)

FE_FIG  = 'figures/'+SIMU_NAME+'/fe/fe_'
FE_CMP  = 'Blues'
FE_TTL  = r'$f_e \left ( x,\,v_x,\,t\right )\,$'
FE_TTL += r'$(n_0/v_{T_{e_0}})$'
X_LBL   = r'$x\,(\lambda_\mathrm{Debye})$'
VX_LBL  = r'$v_x\,(v_{T_{e_0}})$'

N_P    = N_VX*N_X
X_PLT  = []
VX_PLT = []
FE_PLT = []
X_MAP  = np.zeros((N_VX,N_X))
VX_MAP = np.zeros((N_VX,N_X))
FE_MAP = np.zeros((N_VX,N_X))
with open(FE_RES, 'r', encoding='utf-8') as FILE :
    COUNTER = 0
    for LINE in FILE:
        ARRAY = LINE.strip().split()
        VX_PLT.append(float(ARRAY[1]))
        X_PLT.append( float(ARRAY[2]))
        FE_PLT.append(float(ARRAY[3]))
        COUNTER += 1
        if COUNTER % N_P == 0:
            STR_TIME = f"{float(ARRAY[0]):1.4E}"
            print(' * t (/omega_p) = '+STR_TIME)
            N_T   = int(COUNTER / N_P)
            for i in range(0,N_VX):
                for k in range(0,N_X):
                    INDEX        = (N_T-1)*N_P+i*N_X+k
                    X_MAP[i][k]  = X_PLT[INDEX]
                    VX_MAP[i][k] = VX_PLT[INDEX]
                    FE_MAP[i][k] = FE_PLT[INDEX]
            P_MIN   = (N_T-1)*N_P
            P_MAX   = N_T*N_P
            FE_TTL2 = FE_TTL+' at '+r'$t =$'+STR_TIME+r'$\,\omega_p^{-1}$'
            FE_FIG2 = FE_FIG+str(N_T)+'.png'
            X_MIN   = np.amin( X_PLT[P_MIN:P_MAX])
            X_MAX   = np.amax( X_PLT[P_MIN:P_MAX])
            VX_MIN  = np.amin(VX_PLT[P_MIN:P_MAX])
            VX_MAX  = np.amax(VX_PLT[P_MIN:P_MAX])
            FE_MAX  = np.amax(FE_PLT[P_MIN:P_MAX])
            FE_MIN  = np.amin(FE_PLT[P_MIN:P_MAX])
            if FE_MIN == FE_MAX :
                if FE_MIN == 0. :
                    FE_MAX =  1.
                    FE_MIN = -1.
                else :
                    FE_MAX = 1.1*FE_MAX
                    FE_MIN = 0.9*FE_MAX
            lib.make_2d_field_pcolormesh_figure(xmap     =  X_MAP,
                                                ymap     = VX_MAP,
                                                zmap     = FE_MAP,
                                                colormap = FE_CMP,
                                                xmap_min =  X_MIN,
                                                xmap_max =  X_MAX,
                                                ymap_min = VX_MIN,
                                                ymap_max = VX_MAX,
                                                zmap_min = FE_MIN,
                                                zmap_max = FE_MAX,
                                                xlabel   =  X_LBL,
                                                ylabel   = VX_LBL,
                                                title    = FE_TTL2,
                                                eq_axis  = False,
                                                fig_file = FE_FIG2)
