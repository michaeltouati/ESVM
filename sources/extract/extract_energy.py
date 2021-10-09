# ------------------------------------------------------
# SIMULATION RESULTS PLOTTING FOR 1D1V ESVM SIMULATIONS
# ------------------------------------------------------
# Initial commit written by Michaël J Touati - Dec. 2015

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

dir= os.path.dirname("figures/")
if not os.path.exists(dir):
    os.mkdir(dir)

print(' --------------------------')
print(' Total energies scalar plot')
print(' --------------------------')
print('  ')

# total plasma electron kinetic energy 
[tk,uk0] = lib.extract_energy_file('results/UK.dat')
t = tk
# total plasma electron internal energy
[tt,ut0] = lib.extract_energy_file('results/UT.dat')
if len(tt) < len(t):
    t = tt
# total electrostatic energy
[te,ue0] = lib.extract_energy_file('results/UE.dat')
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
print(' Scalar plot :')
print(' * Total kinetic and electrostatic energies')
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
print(' * All total energies in log. scale')
fig=plt.figure()
plt.rc('text', usetex=True)
plt.plot(t, uk,'green',linewidth=2,label=r'$U_{K_e}$')
plt.plot(t, ue,'blue',linewidth=2,label=r'$U_{E_x}$')
leg = plt.legend(loc='best',fontsize=16, fancybox=True)
leg.get_frame().set_alpha(0.5)
plt.xticks(fontsize=16)
plt.xlabel(r'$t\,(\omega_p^{-1})$', fontdict=font)
plt.ylabel(r'$\mathrm{Energy}\,\left (n_0 {\lambda_{\mathrm{Debye}}}^3 m_e {v_{T_{e_0}}}^2 \right )$', fontdict=font)
plt.yticks(fontsize=16)
fig.savefig('figures/energy.png',bbox_inches='tight')
plt.close(fig)
print('  ')
