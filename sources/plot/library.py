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
Python functions library for plotting ESVM simulation results
"""
import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm

FONT_SIZE = 16

FONT = {'style':  'normal',
        'color':  'black',
        'weight': 'normal',
        'size': FONT_SIZE,
        }

def get_results_dir():
    """
    Return ESVM simulation name
    """
    with open('input-deck', 'r', encoding='utf-8') as file :
        for line in file:
            line      = line.strip()
            array     = line.split()
            if array[0] == '#simu' :
                name = array[1]
                break
    to_print=' '+name+' ESVM SIMULATION PLOTS :'
    line_to_print= ' '
    for _ in range(0,len(to_print)-1):
        line_to_print+='='
    print(line_to_print)
    print(to_print)
    print(line_to_print)
    return name

def create_dir(name):
    """
    Create directory entitled name
    """
    dir_name= os.path.dirname(name)
    if not os.path.exists(dir_name):
        os.mkdir(dir_name)

def extract_energy_file(file_name):
    """
    Read and return the two columns of an
    ESVM simulation result energy file entitled file_name
    """
    t_0 = []
    u_0 = []
    with open(file_name, 'r', encoding='utf-8') as file :
        for line in file:
            line      = line.strip()
            array     = line.split()
            t_0.append(float(array[0]))
            u_0.append(float(array[1]))
    return [np.array(t_0),np.array(u_0)]

def make_scalars_plot_figure(**kwargs):
    """
    Plot and save a .png image of scalar fields (maximum 10)
    kwargs keys :
    * Scalars  : xplot_min, xplot_max,
                 yplot_min, yplot_max
    * Vectors  : xplot,
                 yplot1, ... yplot10
    * Strings  : xlabel, ylabel,
                 color1, ..., color10,
                 legend1, ..., legend9
                 title, filename
    * Logicals : logx, logy, grid
    """
    fig=plt.figure()
    plt.rc('text', usetex=True)
    for i in range(1,11):
        condition = ('yplot'+str(i) in kwargs)
        condition = condition and ('legend'+str(i) in kwargs)
        condition = condition and ('color'+str(i) in kwargs)
        if condition :
            yplot  = kwargs['yplot' +str(i)]
            legend = kwargs['legend'+str(i)]
            color  = kwargs['color' +str(i)]
            if kwargs['logx'] :
                if kwargs['logy'] :
                    plt.loglog(kwargs['xplot'],yplot,
                               linewidth=2,label=legend,color=color)
                else:
                    plt.semilogx(kwargs['xplot'],yplot,
                                 linewidth=2,label=legend,color=color)
            else :
                if kwargs['logy'] :
                    plt.semilogy(kwargs['xplot'],yplot,
                                 linewidth=2,label=legend,color=color)
                else :
                    plt.plot(kwargs['xplot'],yplot,
                         linewidth=2,label=legend,color=color)
    leg = plt.legend(fontsize=FONT_SIZE,
                     fancybox=True,
                     bbox_to_anchor=[1., 1.],
                     loc='upper left')
    leg.get_frame().set_alpha(0.5)
    if ('xplot_min' in kwargs) and ('xplot_max' in kwargs) :
        plt.xlim([kwargs['xplot_min'],kwargs['xplot_max']])
    if ('yplot_min' in kwargs) and ('yplot_max' in kwargs) :
        plt.ylim([kwargs['yplot_min'],kwargs['yplot_max']])
    if 'title' in kwargs :
        plt.title(kwargs['title'], fontdict=FONT)
    plt.xticks(fontsize=FONT_SIZE)
    plt.xlabel(kwargs['xlabel'], fontdict=FONT)
    plt.ylabel(kwargs['ylabel'], fontdict=FONT)
    plt.yticks(fontsize=FONT_SIZE)
    if kwargs['grid']:
        plt.grid(which='both', axis='both')
    fig.savefig(kwargs['filename'],bbox_inches='tight')
    plt.close(fig)

def search_nx_nvx(file_name):
    """
    Return ESVM simulation phase-space bins numbers
    """
    print(' Search for the number of phase-space cells:')
    t_0  = []
    v_0  = []
    with open(file_name, 'r', encoding='utf-8') as file :
        line = file.readline()
        line = line.strip()
        array = line.split()
        t_0.append(float(array[0]))
        v_0.append(float(array[1]))
        counter = 0
        for line in file:
            line      = line.strip()
            array     = line.split()
            t_0.append(float(array[0]))
            v_0.append(float(array[1]))
            counter += 1
            if v_0[counter]!=v_0[counter-1]:
                nx_0 = counter
                break
        for line in file:
            line      = line.strip()
            array     = line.split()
            t_0.append(float(array[0]))
            counter = counter + 1
            if t_0[counter]!=t_0[counter-1]:
                nvxnx_0 = counter
                break
    nvx_0 = int(nvxnx_0 / nx_0)
    print(' * found Nx  = '+str(nx_0) +' space bins')
    print(' * found Nvx = '+str(nvx_0)+' velocity bins')
    print('  ')
    return [nx_0,nvx_0]

def make_2d_field_pcolormesh_figure(**kwargs):
    """
    Plot and save a .png image of a 2D field map
    kwargs keys :
    * Scalars   : xmap_min, xmap_max,
                  ymap_min, ymap_max,
                  zmap_min, zmap_max
    * 2D arrays : xmap, ymap, zmap
    * Strings   : xlabel, ylabel,
                  colormap, title,
                  fig_file
    * Logical   : eq_axis
    """
    cmap = plt.get_cmap(kwargs['colormap'])
    norm = cm.colors.Normalize(vmin=kwargs['zmap_min'],
                               vmax=kwargs['zmap_max'])
    fig=plt.figure()
    plt.rc('text', usetex=True)
    plt.pcolormesh(kwargs['xmap'],kwargs['ymap'],kwargs['zmap'],
                   cmap=cmap,norm=norm,shading='gouraud')
    cbar=plt.colorbar(format='%.2f')
    cbar.ax.tick_params(labelsize=FONT_SIZE)
    plt.title(kwargs['title'], fontdict=FONT)
    plt.xticks(fontsize=FONT_SIZE)
    plt.xlabel(kwargs['xlabel'], fontdict=FONT)
    plt.xlim([kwargs['xmap_min'],kwargs['xmap_max']])
    plt.ylabel(kwargs['ylabel'], fontdict=FONT)
    plt.yticks(fontsize=FONT_SIZE)
    plt.ylim([kwargs['ymap_min'],kwargs['ymap_max']])
    if kwargs['eq_axis'] :
        plt.gca().set_aspect('equal')
    txt_condition = 'txt_col' in kwargs
    txt_condition = txt_condition and ('txt_posx' in kwargs)
    txt_condition = txt_condition and ('txt_posy' in kwargs)
    txt_condition = txt_condition and ('txt_str'  in kwargs)
    if txt_condition:
        plt.text(kwargs['txt_posx'], kwargs['txt_posy'],
                 kwargs['txt_str'], color=kwargs['txt_col'],fontsize=FONT_SIZE)
    fig.savefig(kwargs['fig_file'],bbox_inches='tight')
    plt.close(fig)

def mesh_grid_1dhydro(n_x0, n_t0, t_in, x_in, p_in):
    """
    Return 2D arrays from 1D arrays : z_in, x_in, p_in
    """
    t_out = np.zeros((n_x0,n_t0))
    x_out = np.zeros((n_x0,n_t0))
    p_out = np.zeros((n_x0,n_t0))
    for i_x in range(0,n_x0):
        for n_it in range(0,n_t0):
            index = (n_it*n_x0) + i_x
            t_out[i_x][n_it] = t_in[index]
            x_out[i_x][n_it] = x_in[index]
            p_out[i_x][n_it] = p_in[index]
    return [t_out,x_out,p_out]

def get_1d_hydro_quantity(n_x, file_name):
    """
    Extract 1D hydro quantity from ESVM simulation results file
    """
    t_out = []
    x_out = []
    p_out = []
    with open(file_name, 'r', encoding='utf-8') as file :
        counter = 0
        for line in file:
            array = line.strip().split()
            t_out.append(float(array[0]))
            x_out.append(float(array[1]))
            p_out.append(float(array[2]))
            counter += 1
            if counter % n_x == 0:
                n_t = int(counter / n_x)
    return [n_t, t_out,x_out,p_out]

def plot_1d_hydro_quantity_2dmap(n_x, file_name, colormap, plot_title, plot_name):
    """
    Read, plot and save a .png image of one scalar fields as a function of x and t
    """
    [n_t, t_plt, x_plt, p_plt] = get_1d_hydro_quantity(n_x, file_name)
    [t_map,x_map,p_map] = mesh_grid_1dhydro(n_x, n_t, t_plt, x_plt, p_plt)
    maxval = np.amax(p_plt)
    minval = np.amin(p_plt)
    if minval == maxval :
        if  minval == 0. :
            maxval =  1.
            minval = -1.
        else :
            maxval = 1.1*maxval
            minval = 0.9*maxval
    make_2d_field_pcolormesh_figure(xmap     = x_map,
                                    ymap     = t_map,
                                    zmap     = p_map,
                                    colormap = colormap,
                                    xmap_min = np.amin(x_plt),
                                    xmap_max = np.amax(x_plt),
                                    ymap_min = np.amin(t_plt),
                                    ymap_max = np.amax(t_plt),
                                    zmap_min = minval,
                                    zmap_max = maxval,
                                    xlabel   = r'$x\,(\lambda_\mathrm{Debye})$',
                                    ylabel   = r'$t\,(/ \omega_p)$',
                                    title    = plot_title,
                                    eq_axis  = False,
                                    fig_file = plot_name+'.png')

def init_plot_hydro():
    """
    Initialize 1D hydro. quantities plot scripts
    """
    simu_name=get_results_dir()
    print(' -------------------------')
    print(' 1D hydrodynamic moments :')
    print(' -------------------------')
    print('  ')
    create_dir('figures/')
    simu_dir = 'figures/'+simu_name+'/'
    create_dir(simu_dir)
    res_dir = 'results/'+simu_name+'/'
    [n_x,n_vx] = search_nx_nvx(res_dir+'fe.dat')
    return [simu_name, res_dir, simu_dir, n_x, n_vx]

def plot_1d_hydro_quantity_scalar_plot(n_x, file_name, y_label, plot_name):
    """
    Read, plot and save a .png image of one scalar fields as a function of x at all damped times
    """
    x_0   = []
    p_0   = []
    with open(file_name, 'r', encoding='utf-8') as file :
        counter = 0
        for line in file:
            array = line.strip().split()
            x_0.append(float(array[1]))
            p_0.append(float(array[2]))
            counter = counter + 1
            if counter % n_x == 0:
                n_t  = int(counter / n_x)
                time = f"{float(array[0]):1.4E}"
                print('   * t (/omega_p) = '+str(time))
                ex_min = np.amin(p_0)
                if ex_min < 0. :
                    ex_min = 1.1 * ex_min
                else :
                    ex_min = 0.9 * ex_min
                ex_max = np.amax(p_0)
                if ex_max > 0. :
                    ex_max = 1.1 * ex_max
                else :
                    ex_max = 0.9 * ex_max
                if ex_min == ex_max:
                    if  ex_min == 0. :
                        ex_min = -1.
                        ex_max =  1.
                    else :
                        ex_min = 0.9 * ex_max
                        ex_max = 1.1 * ex_max
                fig=plt.figure()
                plt.rc('text', usetex=True)
                plt.plot(x_0[(n_t-1)*n_x:n_t*n_x], p_0[(n_t-1)*n_x:n_t*n_x],
                    linewidth=2, color = 'black')
                plt.title(r'$t =$'+str(time)+r'$\,\omega_p^{-1}$', fontdict=FONT)
                plt.xticks(fontsize=FONT_SIZE)
                plt.xlabel(r'$x\,(\lambda_\mathrm{Debye})$', fontdict=FONT)
                plt.xlim([np.amin(x_0),np.amax(x_0)])
                plt.ylabel(y_label, fontdict=FONT)
                plt.yticks(fontsize=FONT_SIZE)
                plt.ylim([ex_min,ex_max])
                fig.savefig(plot_name+str(n_t)+'.png',bbox_inches='tight')
                plt.close(fig)
