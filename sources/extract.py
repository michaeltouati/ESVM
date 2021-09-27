# -------------------------------------------------
# SIMULATION RESULTS PLOTTING FOR 1D1V ESV SIMU OF 
# ELECTROSTATIC WAKEFIELD EMMISSION BY AN ELECTRON 
# -------------------------------------------------
# Written by Michael J Touati - CLPU - 2020/07/15

import os
import numpy as np
import matplotlib
import math
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from matplotlib import mlab, cm

dir= os.path.dirname("figures/")
if not os.path.exists(dir):
    os.mkdir(dir)

subdir= "results"
os.path.dirname(subdir)
if not os.path.exists(subdir):
      os.mkdir(subdir)

print('Energy plot')

filename = 'results/UK.dat'
t  = []
uk = []
file = open(filename,'r')
for line in file:
    line      = line.strip()
    array     = line.split()
    t.append(float(array[0]))
    uk.append(float(array[1]))
file.close()

filename = 'results/UT.dat'
ut = []
file = open(filename,'r')
for line in file:
    line      = line.strip()
    array     = line.split()
    ut.append(float(array[1]))
file.close()

filename = 'results/UE.dat'
ue = []
file = open(filename,'r')
for line in file:
    line      = line.strip()
    array     = line.split()
    ue.append(float(array[1]))
file.close()

utot = []
for i in range(0,len(t)):
	utot.append(ue[i]+ut[i]+uk[i])

fig=plt.figure()
plt.rc('text', usetex=True)
plt.semilogy(t, uk,'green',linewidth=2,label=r'$\mathbf{U_{K}}$')
plt.semilogy(t, ut,'red',linewidth=2,label=r'$\mathbf{U_{T}}$')
plt.semilogy(t, ue,'blue',linewidth=2,label=r'$\mathbf{U_\mathrm{E}}$')
plt.semilogy(t, utot,'black',linestyle='--',linewidth=2,label='total')
leg = plt.legend(loc='best',fontsize=16, fancybox=True)
leg.get_frame().set_alpha(0.5)
font = {'family': 'non-serif',
		'style':  'normal',
        'color':  'black',
        'weight': 'normal',
        'size': 16,
        }
plt.xticks(fontsize=16)
plt.xlabel(r'$t\,(\omega_p^{-1})$', fontdict=font)
plt.ylabel(r'$\mathrm{Energy}\,(n_0 \lambda_D^3 m_e v_T^2 / 2)$', fontdict=font)
plt.yticks(fontsize=16)
fig.savefig('figures/energy_log.png',bbox_inches='tight')
plt.close(fig)

fig=plt.figure()
plt.rc('text', usetex=True)
plt.plot(t, uk,'green',linewidth=2,label=r'$\mathbf{U_{K}}$')
plt.plot(t, ue,'blue',linewidth=2,label=r'$\mathbf{U_\mathrm{E}}$')
leg = plt.legend(loc='best',fontsize=16, fancybox=True)
leg.get_frame().set_alpha(0.5)
font = {'family': 'non-serif',
		'style':  'normal',
        'color':  'black',
        'weight': 'normal',
        'size': 16,
        }
plt.xticks(fontsize=16)
plt.xlabel(r'$t\,(\omega_p^{-1})$', fontdict=font)
plt.ylabel(r'$\mathrm{Energy}\,(n_0 \lambda_D^3 m_e v_T^2 / 2)$', fontdict=font)
plt.yticks(fontsize=16)
fig.savefig('figures/energy.png',bbox_inches='tight')
plt.close(fig)

# find the number of spatial grid points (Nvx,Nx) by reading 'results.dat'
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
        NxvNx = counter
        break
file.close()
Nxv = NxvNx / Nx

font = {'family': 'non-serif',
        'style':  'normal',
        'color':  'black',
        'weight': 'normal',
        'size': 16,
        }

print('Density plot')

filename = 'results/ne.dat'
title = r'$n_e$'
name='figures/ne'
cmap = 'jet'
title = r'$n_e \left ( x,\,t\right )\,(n_0)$'

t = []
x = []
p = []
file = open(filename, 'r')
counter = 0
for line in file:
    line      = line.strip()
    array     = line.split()
    t.append(float(array[0]))
    x.append(float(array[1]))
    p.append(float(array[2]))
    counter += 1
    if counter % Nx == 0:
      Nt = int(counter / Nx)
file.close()
T = np.zeros((Nx,Nt))
X = np.zeros((Nx,Nt))
P = np.zeros((Nx,Nt))
for i in range(0,Nx):
    for n in range(0,Nt):
        T[i][n] = t[n*Nx+i]
        X[i][n] = x[n*Nx+i]
        P[i][n] = p[n*Nx+i]
cmap = plt.get_cmap(cmap)
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
plt.pcolormesh(X,T,P,cmap=cmap,norm=norm,shading='auto')
cbar=plt.colorbar()
cbar.ax.tick_params(labelsize=16)
plt.title(title, fontdict=font)
plt.xticks(fontsize=16)
plt.xlabel(r'$x\,(\lambda_\mathrm{Debye})$', fontdict=font)
plt.xlim([np.amin(x),np.amax(x)])
plt.ylabel(r'$t\,(/ \omega_p)$', fontdict=font)
plt.yticks(fontsize=16)
plt.ylim([np.amin(t),np.amax(t)])
fig.savefig(name+'.png',bbox_inches='tight')
plt.close(fig)

print('Elecric field plot')

filename = 'results/Ex.dat'
title = r'$E_x$'
name='figures/Ex'
cmap = 'jet'
title = r'$E_x \left ( x,\,t\right )\,(m_e v_T \omega_p / e)$'

t = []
x = []
p = []
file = open(filename, 'r')
counter = 0
for line in file:
    line      = line.strip()
    array     = line.split()
    t.append(float(array[0]))
    x.append(float(array[1]))
    p.append(float(array[2]))
    counter += 1
    if counter % Nx == 0:
      Nt = int(counter / Nx)
file.close()
T = np.zeros((Nx,Nt))
X = np.zeros((Nx,Nt))
P = np.zeros((Nx,Nt))
for i in range(0,Nx):
    for n in range(0,Nt):
        T[i][n] = t[n*Nx+i]
        X[i][n] = x[n*Nx+i]
        P[i][n] = p[n*Nx+i]
cmap = plt.get_cmap(cmap)
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
plt.pcolormesh(X,T,P,cmap=cmap,norm=norm,shading='auto')
cbar=plt.colorbar()
cbar.ax.tick_params(labelsize=16)
plt.title(title, fontdict=font)
plt.xticks(fontsize=16)
plt.xlabel(r'$x\,(\lambda_\mathrm{Debye})$', fontdict=font)
plt.xlim([np.amin(x),np.amax(x)])
plt.ylabel(r'$t\,(/ \omega_p)$', fontdict=font)
plt.yticks(fontsize=16)
plt.ylim([np.amin(t),np.amax(t)])
fig.savefig(name+'.png',bbox_inches='tight')
plt.close(fig)

# distribution function

print('Distribution function plot')
filename = 'results/fe.dat'
N1 = Nxv
N2 = Nx
# cmap = 'gist_ncar'
# cmap = 'nipy_spectral'
cmap = 'Blues'
title = r'$\mathbf{f}\,(n_0.{v_T}^{-3})$'
subdir= "figures/f"
os.path.dirname(subdir)
if not os.path.exists(subdir):
      os.mkdir(subdir)
name='figures/f/f_'

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
        time = math.floor(100.*float(array[0]))/100.
        N = int(counter / N3)
        val=max(abs(np.amax(p[(N-1)*N3:N*N3])),abs(np.amin(p[(N-1)*N3:N*N3])))
        if (val != 0):
            Power=math.floor(np.log(val)/np.log(10))
        else:
            Power = 0
        p[(N-1)*N3:N*N3] = [pp / (10.**Power) for pp in p[(N-1)*N3:N*N3]]
        for i in range(0,N1):
            for k in range(0,N2):
                X[i][k]=x[(N-1)*N3+i*N2+k]
                VX[i][k]=vx[(N-1)*N3+i*N2+k]
                P[i][k]=p[(N-1)*N3+i*N2+k]
        cmap = plt.get_cmap(cmap)
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
        plt.pcolormesh(X,VX,P,cmap=cmap,norm=norm,shading='auto')
        cbar=plt.colorbar()
        cbar.ax.tick_params(labelsize=16)
        plt.title(title+' at '+str(time)+r'$\,\omega_p^{-1}$', fontdict=font)
        plt.xticks(fontsize=16)
        plt.xlabel(r'$x\,(\lambda_\mathrm{Debye})$', fontdict=font)
        plt.xlim([np.amin(x[(N-1)*N3:N*N3]),np.amax(x[(N-1)*N3:N*N3])])
        plt.ylabel(r'$v_x\,(v_\mathrm{T})$', fontdict=font)
        plt.yticks(fontsize=16)
        plt.ylim([np.amin(vx[(N-1)*N3:N*N3]),np.amax(vx[(N-1)*N3:N*N3])])
        #plt.axes().set_aspect('equal')
        fig.savefig(name+str(N)+'.png',bbox_inches='tight')
        plt.close(fig)
file.close()


# probes
print('Ex plot probe')

filename = 'results/Ex.dat'
title = r'$E_x$'
subdir= "figures/Ex"
os.path.dirname(subdir)
if not os.path.exists(subdir):
      os.mkdir(subdir)
name='figures/Ex/Ex_'

x = []
p = []
X = np.zeros(Nx)
P = np.zeros(Nx)
file = open(filename, 'r')
counter = 0
for line in file:
    line      = line.strip()
    array     = line.split()
    x.append(float(array[1]))
    p.append(float(array[2]))
    counter = counter + 1
    if counter % Nx == 0:
        N = int(counter / Nx)
        for i in range(0,Nx):
            X[i] = x[(N-1)*Nx+i]
            P[i] = p[(N-1)*Nx+i]
        time = math.floor(100.*float(array[0]))/100.
        fig=plt.figure()
        plt.rc('text', usetex=True)
        plt.plot(X,P)
        plt.title(title+' at '+str(time)+r'$\,\omega_p^{-1}$', fontdict=font)
        plt.xticks(fontsize=16)
        plt.xlabel(r'$x\,(\lambda_\mathrm{Debye})$', fontdict=font)
        plt.xlim([np.amin(X),np.amax(X)])
        plt.ylabel(r'$E_x\,(m_e \omega_p v_T / e)$', fontdict=font)
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

print('ne plot probe')
filename = 'results/ne.dat'
title = r'$n_e$'
subdir= "figures/ne"
os.path.dirname(subdir)
if not os.path.exists(subdir):
      os.mkdir(subdir)
name='figures/ne/ne_'

x = []
p = []
X = np.zeros(Nx)
P = np.zeros(Nx)
file = open(filename, 'r')
counter = 0
for line in file:
    line      = line.strip()
    array     = line.split()
    x.append(float(array[1]))
    p.append(float(array[2]))
    counter = counter + 1
    if counter % Nx == 0:
        N = int(counter / Nx)
        for i in range(0,Nx):
            X[i] = x[(N-1)*Nx+i]
            P[i] = p[(N-1)*Nx+i]
        time = math.floor(100.*float(array[0]))/100.
        fig=plt.figure()
        plt.rc('text', usetex=True)
        plt.plot(X,P)
        plt.title(title+' at '+str(time)+r'$\,\omega_p^{-1}$', fontdict=font)
        plt.xticks(fontsize=16)
        plt.xlabel(r'$x\,(\lambda_\mathrm{Debye})$', fontdict=font)
        plt.xlim([np.amin(X),np.amax(X)])
        plt.ylabel(r'$n_e\,(n_0)$', fontdict=font)
        plt.yticks(fontsize=16)
        ne_min = np.amin(p)
        ne_max = np.amax(p)
        if ( ne_min == ne_max ):
            if (ne_min == 0.) :
                ne_min = -1.
                ne_max =  1.
            else :
                ne_min = 0.9 * ne_max
                ne_max = 1.1 * ne_max
        plt.ylim([ne_min,ne_max])
        fig.savefig(name+str(N)+'.png',bbox_inches='tight')
        plt.close(fig)
file.close()

print('Current density plot probe')
filename = 'results/je.dat'
title = r'$j_e$'
subdir= "figures/je"
os.path.dirname(subdir)
if not os.path.exists(subdir):
      os.mkdir(subdir)
name='figures/je/je_'

x = []
p = []
X = np.zeros(Nx)
P = np.zeros(Nx)
file = open(filename, 'r')
counter = 0
for line in file:
    line      = line.strip()
    array     = line.split()
    x.append(float(array[1]))
    p.append(float(array[2]))
    counter = counter + 1
    if counter % Nx == 0:
        N = int(counter / Nx)
        for i in range(0,Nx):
            X[i] = x[(N-1)*Nx+i]
            P[i] = p[(N-1)*Nx+i]
        time = math.floor(100.*float(array[0]))/100.
        fig=plt.figure()
        plt.rc('text', usetex=True)
        plt.plot(X,P)
        plt.title(title+' at '+str(time)+r'$\,\omega_p^{-1}$', fontdict=font)
        plt.xticks(fontsize=16)
        plt.xlabel(r'$x\,(\lambda_\mathrm{Debye})$', fontdict=font)
        plt.xlim([np.amin(X),np.amax(X)])
        plt.ylabel(r'$j_e\,(n_0 \,e \,v_T)$', fontdict=font)
        plt.yticks(fontsize=16)
        je_min = np.amin(p)
        je_max = np.amax(p)
        if ( je_min == je_max ):
            if (je_min == 0.) :
                je_min = -1.
                je_max =  1.
            else :
                je_min = 0.9 * je_max
                je_max = 1.1 * je_max
        plt.ylim([je_min,je_max])
        fig.savefig(name+str(N)+'.png',bbox_inches='tight')
        plt.close(fig)
file.close()

print('Mean velocity plot probe')
filename = 'results/ve.dat'
title = r'$v_e$'
subdir= "figures/ve"
os.path.dirname(subdir)
if not os.path.exists(subdir):
      os.mkdir(subdir)
name='figures/ve/ve_'

x = []
p = []
X = np.zeros(Nx)
P = np.zeros(Nx)
file = open(filename, 'r')
counter = 0
for line in file:
    line      = line.strip()
    array     = line.split()
    x.append(float(array[1]))
    p.append(float(array[2]))
    counter = counter + 1
    if counter % Nx == 0:
        N = int(counter / Nx)
        for i in range(0,Nx):
            X[i] = x[(N-1)*Nx+i]
            P[i] = p[(N-1)*Nx+i]
        time = math.floor(100.*float(array[0]))/100.
        fig=plt.figure()
        plt.rc('text', usetex=True)
        plt.plot(X,P)
        plt.title(title+' at '+str(time)+r'$\,\omega_p^{-1}$', fontdict=font)
        plt.xticks(fontsize=16)
        plt.xlabel(r'$x\,(\lambda_\mathrm{Debye})$', fontdict=font)
        plt.xlim([np.amin(X),np.amax(X)])
        plt.ylabel(r'$v_e\,(v_T)$', fontdict=font)
        plt.yticks(fontsize=16)
        ve_min = np.amin(p)
        ve_max = np.amax(p)
        if ( ve_min == ve_max ):
            if (ve_min == 0.) :
                ve_min = -1.
                ve_max =  1.
            else :
                ve_min = 0.9 * ve_max
                ve_max = 1.1 * ve_max
        plt.ylim([ve_min,ve_max])
        fig.savefig(name+str(N)+'.png',bbox_inches='tight')
        plt.close(fig)
file.close()

print('Thermal velocity plot probe')
filename = 'results/vTe.dat'
title = r'$v_{T,e}$'
subdir= "figures/vte"
os.path.dirname(subdir)
if not os.path.exists(subdir):
      os.mkdir(subdir)
name='figures/vte/vte_'

x = []
p = []
X = np.zeros(Nx)
P = np.zeros(Nx)
file = open(filename, 'r')
counter = 0
for line in file:
    line      = line.strip()
    array     = line.split()
    x.append(float(array[1]))
    p.append(float(array[2]))
    counter = counter + 1
    if counter % Nx == 0:
        N = int(counter / Nx)
        for i in range(0,Nx):
            X[i] = x[(N-1)*Nx+i]
            P[i] = p[(N-1)*Nx+i]
        time = math.floor(100.*float(array[0]))/100.
        fig=plt.figure()
        plt.rc('text', usetex=True)
        plt.plot(X,P)
        plt.title(title+' at '+str(time)+r'$\,\omega_p^{-1}$', fontdict=font)
        plt.xticks(fontsize=16)
        plt.xlabel(r'$x\,(\lambda_\mathrm{Debye})$', fontdict=font)
        plt.xlim([np.amin(X),np.amax(X)])
        plt.ylabel(r'$v_{T,e}\,(v_T)$', fontdict=font)
        plt.yticks(fontsize=16)
        vTe_min = np.amin(p)
        vTe_max = np.amax(p)
        if ( vTe_min == vTe_max ):
            if (vTe_min == 0.) :
                vTe_min = -1.
                vTe_max =  1.
            else :
                vTe_min = 0.9 * vTe_max
                vTe_max = 1.1 * vTe_max
        plt.ylim([vTe_min,vTe_max])
        fig.savefig(name+str(N)+'.png',bbox_inches='tight')
        plt.close(fig)
file.close()