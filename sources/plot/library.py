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

def get_results_dir():
    file = open('input-deck','r')
    for line in file:
        line      = line.strip()
        array     = line.split()
        if (array[0] == '#simu') :
            name = array[1]
            break
    file.close()
    return name

def create_dir(name):
    dir_name= os.path.dirname(name)
    if not os.path.exists(dir_name):
        os.mkdir(dir_name)

def extract_energy_file(file_name):
    t0 = []
    u0 = []
    file = open(file_name,'r')
    for line in file:
        line      = line.strip()
        array     = line.split()
        t0.append(float(array[0]))
        u0.append(float(array[1]))
    file.close()
    return [t0,u0]

def search_Nx_Nvx(file_name):
	t0  = []
	v0  = []
	file = open(file_name, 'r')
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
	        Nx0 = counter
	        break
	for line in file:
	    line      = line.strip()
	    array     = line.split()
	    t0.append(float(array[0]))
	    counter = counter + 1
	    if t0[counter]!=t0[counter-1]:
	        NvxNx0 = counter
	        break
	file.close()
	Nvx0 = int(NvxNx0 / Nx0)
	return [Nx0,Nvx0]

def plot_1D_hydro_quantity_2Dmap(N_x, file_name, font_properties, colormap, plot_title,plot_name):
    t = []
    x = []
    p = []
    file = open(file_name, 'r')
    counter = 0
    for line in file:
        line      = line.strip()
        array     = line.split()
        t.append(float(array[0]))
        x.append(float(array[1]))
        p.append(float(array[2]))
        counter += 1
        if counter % N_x == 0:
          N_t = int(counter / N_x)
    file.close()
    T = np.zeros((N_x,N_t))
    X = np.zeros((N_x,N_t))
    P = np.zeros((N_x,N_t))
    for i in range(0,N_x):
        for n in range(0,N_t):
            T[i][n] = t[n*N_x+i]
            X[i][n] = x[n*N_x+i]
            P[i][n] = p[n*N_x+i]
    cmap = plt.get_cmap(colormap)
    Maxval = np.amax(p)
    Minval = np.amin(p)
    if ( Minval == Maxval ) :
        if ( Minval == 0. ) :
            Maxval =  1.
            Minval = -1.
        else :
            Maxval = 1.1*Maxval
            Minval = 0.9*Maxval
    norm = cm.colors.Normalize(vmax=Maxval, vmin=Minval)
    fig=plt.figure()
    plt.rc('text', usetex=True)
    plt.pcolormesh(X,T,P,cmap=cmap,norm=norm,shading='gouraud')
    cbar=plt.colorbar()
    cbar.ax.tick_params(labelsize=16)
    plt.title(plot_title, fontdict=font_properties)
    plt.xticks(fontsize=16)
    plt.xlabel(r'$x\,(\lambda_\mathrm{Debye})$', fontdict=font_properties)
    plt.xlim([np.amin(x),np.amax(x)])
    plt.ylabel(r'$t\,(/ \omega_p)$', fontdict=font_properties)
    plt.yticks(fontsize=16)
    plt.ylim([np.amin(t),np.amax(t)])
    fig.savefig(plot_name+'.png',bbox_inches='tight')
    plt.close(fig)

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
            print('   * t (/omega_p) = '+str(time))
            fig=plt.figure()
            plt.rc('text', usetex=True)
            plt.plot(X, P, linewidth=2, color = 'black')
            plt.title(r'$t =$'+str(time)+r'$\,\omega_p^{-1}$', fontdict=font_properties)
            plt.xticks(fontsize=16)
            plt.xlabel(r'$x\,(\lambda_\mathrm{Debye})$', fontdict=font_properties)
            plt.xlim([np.amin(X),np.amax(X)])
            plt.ylabel(y_label, fontdict=font_properties)
            plt.yticks(fontsize=16)
            ex_min = np.amin(p)
            if (ex_min < 0.) :
                ex_min = 1.1 * ex_min
            else :
                ex_min = 0.9 * ex_min
            ex_max = np.amax(p)
            if (ex_max > 0.) :
                ex_max = 1.1 * ex_max
            else :
                ex_max = 0.9 * ex_max
            if ( ex_min == ex_max ):
                if (ex_min == 0.) :
                    ex_min = -1.
                    ex_max =  1.
                else :
                    ex_min = 0.9 * ex_max
                    ex_max = 1.1 * ex_max
            plt.ylim([ex_min,ex_max])
            fig.savefig(plot_name+str(N)+'.png',bbox_inches='tight')
            plt.close(fig)
    file.close()
