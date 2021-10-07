# -----------------------------------------------
# SIMULATION RESULTS PLOTTING FOR 1D1V ESVM SIMU  
# -----------------------------------------------
# Written by Michael J Touati

import os
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from matplotlib import mlab, cm

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
# cmap_fe = 'gist_ncar'
# cmap_fe = 'nipy_spectral'
#cmap_fe = 'gist_stern'
#cmap_fe = 'gnuplot'
#cmap_fe = 'gnuplot2'
cmap_fe = 'Blues'
#cmap_fe = 'jet'

#############################
# Create figures/ directory #
#############################

dir= os.path.dirname("figures/")
if not os.path.exists(dir):
    os.mkdir(dir)

#################################
# functions used in this script #
#################################

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
            T[i][n] = t[n*Nx+i]
            X[i][n] = x[n*Nx+i]
            P[i][n] = p[n*Nx+i]
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
    plt.ylabel(r'$t\,(/ \omega_p)$', fontdict=font)
    plt.yticks(fontsize=16)
    plt.ylim([np.amin(t),np.amax(t)])
    fig.savefig(name+'.png',bbox_inches='tight')
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
            print('  t (/omega_p) = '+str(time))
            fig=plt.figure()
            plt.rc('text', usetex=True)
            plt.plot(X, P, linewidth=2, color = 'black')
            plt.title(r'$t =$'+str(time)+r'$\,\omega_p^{-1}$', fontdict=font_properties)
            plt.xticks(fontsize=16)
            plt.xlabel(r'$x\,(\lambda_\mathrm{Debye})$', fontdict=font)
            plt.xlim([np.amin(X),np.amax(X)])
            plt.ylabel(y_label, fontdict=font)
            plt.yticks(fontsize=16)
            ex_min = np.amin(p)
            ex_max = np.amax(p)
            if ( ex_min == ex_max ):
                if (ex_min == 0.) :
                    ex_min = -1.
                    ex_max =  1.
                else :
                    ex_min = 0.9 * ex_max
                    ex_max = 1.1 * ex_max
            plt.ylim([ex_min,ex_max])
            fig.savefig(name+str(N)+'.png',bbox_inches='tight')
            plt.close(fig)
    file.close()

###############################
# Total energies scalar plots #
###############################

print('Total energies scalar plot')
# total plasma electron kinetic energy 
[tk,uk0] = extract_energy_file('results/UK.dat')
t = tk
# total plasma electron internal energy
[tt,ut0] = extract_energy_file('results/UT.dat')
if len(tt) < len(t):
    t = tt
# total electrostatic energy
[te,ue0] = extract_energy_file('results/UE.dat')
if len(te) < len(t):
    t = te
#
uk   = []
ut   = []
ue   = []
utot = []
for i in range(0,len(t)):
    uk.append(uk0[i])
    ut.append(ut0[i])
    ue.append(ue0[i])
    utot.append(ue[i]+ut[i]+uk[i])
#
fig=plt.figure()
plt.rc('text', usetex=True)
plt.semilogy(t, uk,'green',linewidth=2,label=r'$U_{K_e}$')
plt.semilogy(t, ut,'red',linewidth=2,label=r'$U_{T_e}$')
plt.semilogy(t, ue,'blue',linewidth=2,label=r'$U_{E_x}$')
plt.semilogy(t, utot,'black',linestyle='--',linewidth=2,label='Total')
leg = plt.legend(loc='best',fontsize=16, fancybox=True)
leg.get_frame().set_alpha(0.5)
plt.xticks(fontsize=16)
plt.xlabel(r'$t\,(\omega_p^{-1})$', fontdict=font)
plt.ylabel(r'$\mathrm{Energy}\,\left (n_0 {\lambda_{\mathrm{Debye}}}^3 m_e {v_{T_{e_0}}}^2 \right )$', fontdict=font)
plt.yticks(fontsize=16)
fig.savefig('figures/energy_log.png',bbox_inches='tight')
plt.close(fig)
# in log scale
fig=plt.figure()
plt.rc('text', usetex=True)
plt.plot(t, uk,'green',linewidth=2,label=r'$\mathbf{U_{K}}$')
plt.plot(t, ue,'blue',linewidth=2,label=r'$\mathbf{U_\mathrm{E}}$')
leg = plt.legend(loc='best',fontsize=16, fancybox=True)
leg.get_frame().set_alpha(0.5)
plt.xticks(fontsize=16)
plt.xlabel(r'$t\,(\omega_p^{-1})$', fontdict=font)
plt.ylabel(r'$\mathrm{Energy}\,\left (n_0 {\lambda_{\mathrm{Debye}}}^3 m_e {v_{T_{e_0}}}^2 \right )$', fontdict=font)
plt.yticks(fontsize=16)
fig.savefig('figures/energy.png',bbox_inches='tight')
plt.close(fig)
print(' ')

##########################################
# Find the number of spatial grid points # 
# (Nvx,Nx) by reading 'results/fe.dat'   #
##########################################

print('Search for the number of phase-space cells:')
t0  = []
v0  = []
file = open('results/fe.dat', 'r')
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
        Nx = counter
        break
for line in file:
    line      = line.strip()
    array     = line.split()
    t0.append(float(array[0]))
    counter = counter + 1
    if t0[counter]!=t0[counter-1]:
        NvxNx = counter
        break
file.close()
Nvx = int(NvxNx / Nx)

print('* found Nx  = '+str(Nx) +' space cells')
print('* found Nvx = '+str(Nvx)+' velocity cells')
print(' ')

################################
# 2D space-time colormap plots #
################################

print('2D space-time colormap plots:')
print('* 1D elecrostatic field')
filename = 'results/Ex.dat'
name     ='figures/Ex'
title    = r'$E_x \left ( x,\,t\right )\,(m_e v_{T_{e_0}} \omega_p / e)$'
plot_1D_hydro_quantity_2Dmap(Nx, filename, font, cmap_Ex, title, name)

print('* 1D plasma electron density')
filename = 'results/ne.dat'
title    = r'$n_e \left ( x,\,t\right )\,(n_0)$'
name     ='figures/ne'
plot_1D_hydro_quantity_2Dmap(Nx, filename, font, cmap_ne, title, name)

print('* 1D plasma electron mean velocity')
filename = 'results/ve.dat'
name     ='figures/ve'
cmap     = 'jet'
title    = r'$v_e \left ( x,\,t\right )\,(v_{T_{e_0}})$'
plot_1D_hydro_quantity_2Dmap(Nx, filename, font, cmap_ve, title, name)

filename = 'results/je.dat'
cmap     = 'jet'
title    = r'$j_e \left ( x,\,t\right )\,(n_0 e v_{T_{e_0}})$'
name     ='figures/je'
plot_1D_hydro_quantity_2Dmap(Nx, filename, font, cmap_je, title, name)

print('* 1D plasma electron thermal velocity (standard deviation)')
filename = 'results/vTe.dat'
cmap     = 'jet'
name     ='figures/vte'
title    = r'$v_{T_e} \left ( x,\,t\right )\,(v_{T_{e_0}})$'
plot_1D_hydro_quantity_2Dmap(Nx, filename, font, cmap_vte, title, name)
print(' ')

###############
# Scalar plots #
###############

print('1D scalar plots :')
print('* 1D Electrostatic field')
filename = 'results/Ex.dat'
Ylabel   = r'$E_x \left ( x,\,t\right ) \, \left( m_e \omega_p v_{T_{e_0}} / e \right )$'
subdir   = "figures/Ex"
os.path.dirname(subdir)
if not os.path.exists(subdir):
      os.mkdir(subdir)
name     ='figures/Ex/Ex_'
plot_1D_hydro_quantity_scalar_plot(Nx, filename, font, Ylabel, name)

print('* 1D plasma electron density')
filename = 'results/ne.dat'
Ylabel   = r'$n_e \left ( x,\,t\right ) \, \left( n_0 \right )$'
subdir   = "figures/ne"
os.path.dirname(subdir)
if not os.path.exists(subdir):
      os.mkdir(subdir)
name     ='figures/ne/ne_'
plot_1D_hydro_quantity_scalar_plot(Nx, filename, font, Ylabel, name)

print('* 1D plasma electron mean velocity')
filename = 'results/ve.dat'
Ylabel   = r'$v_e \left ( x,\,t\right ) \, \left( v_{T_{e_0}}\right )$'
subdir   = "figures/ve"
os.path.dirname(subdir)
if not os.path.exists(subdir):
      os.mkdir(subdir)
name    ='figures/ve/ve_'
plot_1D_hydro_quantity_scalar_plot(Nx, filename, font, Ylabel, name)

print('* 1D plasma electron current density')
filename = 'results/je.dat'
Ylabel   = r'$j_e \left ( x,\,t\right ) \, \left( n_0 e v_{T_{e_0}}\right )$'
subdir   = "figures/je"
os.path.dirname(subdir)
if not os.path.exists(subdir):
      os.mkdir(subdir)
name     ='figures/je/je_'
plot_1D_hydro_quantity_scalar_plot(Nx, filename, font, Ylabel, name)

print('* 1D plasma electron thermal velocity (standard deviation)')
filename = 'results/vTe.dat'
Ylabel   = r'$v_{T_e} \left ( x,\,t\right ) \, \left( v_{T_{e_0}}\right )$'
subdir   = "figures/vte"
os.path.dirname(subdir)
if not os.path.exists(subdir):
      os.mkdir(subdir)
name     ='figures/vte/vte_'
plot_1D_hydro_quantity_scalar_plot(Nx, filename, font, Ylabel, name)

#######################################################
# 1D1V plasma electron distribution function colormap #
#######################################################

print('1D1V plasma electron distribution function phase-space colormap plot')
filename = 'results/fe.dat'
N1 = Nvx
N2 = Nx
title = r'$f_e \left ( x,\,v_x,\,t\right )\,(n_0/v_{T_{e_0}})$'
subdir= "figures/fe"
os.path.dirname(subdir)
if not os.path.exists(subdir):
      os.mkdir(subdir)
name='figures/fe/fe_'

N3 = int(N1*N2)
vx = []
x = []
p = []
N1 = int(N1)
N2 = int(N2)
N3 = int(N3)
X = np.zeros((N1,N2))
VX = np.zeros((N1,N2))
P = np.zeros((N1,N2))
file = open(filename, 'r')
counter = 0
for line in file:
    line      = line.strip()
    array     = line.split()
    vx.append(float(array[1]))
    x.append(float(array[2]))
    p.append(float(array[3]))
    counter = counter + 1
    if counter % N3 == 0:
        time = int(100.*float(array[0]))/100.
        print(' t (/omega_p) = '+str(time))
        N = int(counter / N3)
        val=max(abs(np.amax(p[(N-1)*N3:N*N3])),abs(np.amin(p[(N-1)*N3:N*N3])))
        if (val != 0):
            Power=int(np.log(val)/np.log(10))
        else:
            Power = 0
        p[(N-1)*N3:N*N3] = [pp / (10.**Power) for pp in p[(N-1)*N3:N*N3]]
        for i in range(0,N1):
            for k in range(0,N2):
                X[i][k]=x[(N-1)*N3+i*N2+k]
                VX[i][k]=vx[(N-1)*N3+i*N2+k]
                P[i][k]=p[(N-1)*N3+i*N2+k]
        cmap = plt.get_cmap(cmap_fe)
        Maxval = np.amax(p[(int(N)-1)*int(N3):int(N)*int(N3)])
        Minval = np.amin(p[(int(N)-1)*int(N3):int(N)*int(N3)])
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
        plt.pcolormesh(X,VX,P,cmap=cmap,norm=norm,shading='gouraud')
        cbar=plt.colorbar()
        cbar.ax.tick_params(labelsize=16)
        plt.title(title+' at '+r'$t =$'+str(time)+r'$\,\omega_p^{-1}$', fontdict=font)
        plt.xticks(fontsize=16)
        plt.xlabel(r'$x\,(\lambda_\mathrm{Debye})$', fontdict=font)
        plt.xlim([np.amin(x[(N-1)*N3:N*N3]),np.amax(x[(N-1)*N3:N*N3])])
        plt.ylabel(r'$v_x\,(v_{T_{e_0}})$', fontdict=font)
        plt.yticks(fontsize=16)
        plt.ylim([np.amin(vx[(N-1)*N3:N*N3]),np.amax(vx[(N-1)*N3:N*N3])])
        #plt.axes().set_aspect('equal')
        fig.savefig(name+str(N)+'.png',bbox_inches='tight')
        plt.close(fig)
file.close()