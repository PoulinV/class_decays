
# coding: utf-8

# In[ ]:

# import necessary modules
# uncomment to get plots displayed in notebook
# uncomment to get plots displayed in notebook
#%matplotlib inline
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from classy import Class
from scipy.optimize import fsolve
from scipy.interpolate import interp1d
import math


# In[ ]:

# esthetic definitions for the plots
font = {'size'   : 16, 'family':'STIXGeneral'}
axislabelfontsize='large'
matplotlib.rc('font', **font)
matplotlib.mathtext.rcParams['legend.fontsize']='medium'
plt.rcParams["figure.figsize"] = [8.0,6.0]


# In[ ]:

#############################################
#
# User settings controlling the figure aspect
#
z_max_pk = 46000       # highest redshift involved
k_per_decade = 400     # number of k values, controls final resolution
k_min_tau0 = 40.       # this value controls the minimum k value in the figure (it is k_min * tau0)
P_k_max_inv_Mpc =1.0   # this value is directly the maximum k value in the figure in Mpc
tau_num_early = 2000   # number of conformal time values before recombination, controls final resolution
tau_num_late = 200     # number of conformal time values after recombination, controls final resolution
tau_ini = 10.          # first value of conformal time in Mpc
tau_label_Hubble = 20. # value of time at which we want to place the label on Hubble crossing
tau_label_ks = 40.     # value of time at which we want to place the label on sound horizon crossing
tau_label_kd = 230.    # value of time at which we want to place the label on damping scale crossing
#
# Cosmological parameters and other CLASS parameters
#
common_settings = {# which output? transfer functions only
                   'output':'mTk',
                   # LambdaCDM parameters
                   'h':0.67556,
                   'omega_b':0.022032,
                   'omega_cdm':0.12038,
                   'A_s':2.215e-9,
                   'n_s':0.9619,
                   'tau_reio':0.0925,
                   'N_ncdm':1,
                   'N_ur':2.0328,
                   'background_ncdm_distribution':0,
                   'ncdm_fluid_approximation':2,
                   'Number of momentum bins perturbs': 300,
                   'm_ncdm':0.06,
                   # other output and precision parameters
                   'z_max_pk':z_max_pk,
                   'recfast_z_initial':z_max_pk,
                   #'k_step_sub':'0.01',
                   'k_output_values':'0.1,1',
                   'k_per_decade_for_pk':k_per_decade,
                   'k_per_decade_for_bao':k_per_decade,
                   'k_min_tau0':k_min_tau0, # this value controls the minimum k value in the figure
                   'perturb_sampling_stepsize':'0.05',
                   'P_k_max_1/Mpc':P_k_max_inv_Mpc
                   }

###############
#
# call CLASS
#
###############
M = Class()
M.set(common_settings)
M.compute()
#
# define conformal time sampling array
#
times = M.get_current_derived_parameters(['tau_rec','conformal_age'])
tau_rec=times['tau_rec']
tau_0 = times['conformal_age']
tau1 = np.logspace(math.log10(tau_ini),math.log10(tau_rec),tau_num_early)
tau2 = np.logspace(math.log10(tau_rec),math.log10(tau_0),tau_num_late)[1:]
tau2[-1] *= 0.999 # this tiny shift avoids interpolation errors
tau = np.concatenate((tau1,tau2))
tau_num = len(tau)
#
# use table of background and thermodynamics quantitites to define some functions
# returning some characteristic scales
# (of Hubble crossing, sound horizon crossing, etc.) at different time
#
background = M.get_background() # load background table
#print(background.keys())
thermodynamics = M.get_thermodynamics() # load thermodynamics table
#print(thermodynamics.keys())

#
background_tau = background['conf. time [Mpc]'] # read conformal times in background table
background_z = background['z'] # read redshift
background_aH = 2.*math.pi*background['H [1/Mpc]']/(1.+background['z'])/M.h() # read 2pi * aH in [h/Mpc]
background_rho_m_over_r =    (background['(.)rho_b']+background['(.)rho_cdm'])    /(background['(.)rho_g']+background['(.)rho_ur']) # read rho_r / rho_m (to find time of equality)
background_rho_l_over_m =    background['(.)rho_lambda']    /(background['(.)rho_b']+background['(.)rho_cdm']) # read rho_m / rho_lambda (to find time of equality)
thermodynamics_tau = thermodynamics['conf. time [Mpc]'] # read confromal times in thermodynamics table
#
# define a bunch of interpolation functions based on previous quantities
#
background_z_at_tau = interp1d(background_tau,background_z)
background_aH_at_tau = interp1d(background_tau,background_aH)
background_tau_at_mr = interp1d(background_rho_m_over_r,background_tau)
background_tau_at_lm = interp1d(background_rho_l_over_m,background_tau)
#
# infer arrays of characteristic quantitites calculated at values of conformal time in tau array
#
aH = background_aH_at_tau(tau)
#
# infer times of R/M and M/Lambda equalities
#
tau_eq = background_tau_at_mr(1.)
tau_lambda = background_tau_at_lm(1.)
#
# check and inform user whether intiial arbitrary choice of z_max_pk was OK
max_z_needed = background_z_at_tau(tau[0])
if max_z_needed > z_max_pk:
    print(r'must increase the value of z_max_pk = %f to at least '%max_z_needed)
    () + 1  # this strange line is just a trick to stop the script execution there
else:
    print(r'in a next run with the same values of tau, you may decrease z_max_pk from %f to %f '%(z_max_pk,max_z_needed))
#
# get transfer functions at each time and build arrays cs2(tau,k) and phi(tau,k)
#
for i in range(tau_num):
    one_time = M.get_transfer(background_z_at_tau(tau[i])) # transfer functions at each time tau
    if i ==0:   # if this is the first time in the loop: create the arrays (k, cs2, phi)
        k = one_time['k (h/Mpc)']
        k_num = len(k)
        cs2 = np.zeros((tau_num,k_num))
        phi = np.zeros((tau_num,k_num))
    cs2[i,:] = np.log(np.abs(one_time['cs2_ncdm[0]'][:]))
#    cs2[i,:] = one_time['c_n'][:]
    phi[i,:] = one_time['phi'][:]
 
print(one_time.keys()) 



#%%
#all_k = M.get_perturbations()
#print(all_k['scalar'][0].keys()) 
#one_k0 = all_k['scalar'][0]
#one_k1 = all_k['scalar'][1]
#tau0 = one_k0['tau [Mpc]']
#phi0 = one_k0['phi'][:]
#tau1 = one_k1['tau [Mpc]']
#phi1 = one_k1['phi'][:]

#this allows to obtain the evolution of each perturbation at different wavenumbers
#the problem is that we cannot generate an array because the code does an adaptative integration,
# not all k-modes are computed using the same number of time steps, in fact len(tau0) neq len(tau1)
#%%
# find the global extra of cs2(tau,k) and phi(tau,k), used to define color code later
#
#cs2_amp = max(cs2.max(),-cs2.min())
cs2_amp_max = cs2.max()  #it shouldn't be bigger than log(1/3) = -0.4771, otherwise its numerical error
cs2_amp_min = cs2.min()

cs2_amp_min 

phi_amp = max(phi.max(),-phi.min())

#
# reshaping of (k,tau) necessary to call the function 'pcolormesh'
#
K,T = np.meshgrid(k,tau)
#
# inform user of the size of the grids (related to the figure resolution)
#
#print(r'grid size: %f %f %f'%(len(k),len(tau),cs2.shape))
#
#################
#
# start plotting
#
#################
#
#fig = plt.figure(figsize=(18,8))
fig = plt.figure(figsize=(10,8))
#
# plot cs2(k,tau)
#
ax_cs2 = fig.add_subplot(121)
print(r'> Plotting cs2')
fig_cs2 = ax_cs2.pcolormesh(K,T,cs2,cmap='coolwarm',vmin=cs2_amp_min, vmax=cs2_amp_max) #,shading='gouraud')
print(r'> Done')
#
# plot lines (characteristic times and scales)
#
#ax_cs2.axhline(y=tau_rec,color='k',linestyle='-')
ax_cs2.axhline(y=tau_eq,color='k',linestyle='-')
ax_cs2.axhline(y=tau_lambda,color='k',linestyle='-')
ax_cs2.plot(aH,tau,'r-',linewidth=2)
#
# dealing with labels
#
ax_cs2.set_title(r'$\mathrm{log}_{10}(c_s^2)$')
#ax_cs2.text(1.5*k[0],0.9*tau_rec,r'$\mathrm{rec.}$')
ax_cs2.text(1.5*k[0],0.9*tau_eq,r'$\mathrm{R/M} \,\, \mathrm{eq.}$')
ax_cs2.text(1.5*k[0],0.9*tau_lambda,r'$\mathrm{M/L} \,\, \mathrm{eq.}$')
ax_cs2.annotate(r'$\mathrm{Hubble} \,\, \mathrm{cross.}$',
                  xy=(background_aH_at_tau(tau_label_Hubble),tau_label_Hubble),
                  xytext=(0.1*background_aH_at_tau(tau_label_Hubble),0.8*tau_label_Hubble),
                  arrowprops=dict(facecolor='black', shrink=0.05, width=1, headlength=5, headwidth=5))

# dealing with axes
#
ax_cs2.set_xlim(k[0],k[-1])
ax_cs2.set_xscale('log')
ax_cs2.set_yscale('log')
ax_cs2.set_xlabel(r'$k  \,\,\, \mathrm{[h/Mpc]}$')
ax_cs2.set_ylabel(r'$\tau   \,\,\, \mathrm{[Mpc]}$')
ax_cs2.invert_yaxis()
#
# color legend
#
fig.colorbar(fig_cs2)
#
# plot phi(k,tau)
#
#ax_phi = fig.add_subplot(122)
#ax_phi.set_xlim(k[0],k[-1])

#print(r'> Plotting phi')
#fig_phi = ax_phi.pcolormesh(K,T,phi,cmap='coolwarm',vmin=-0., vmax=phi_amp)
#print(r'> Done')
#
# plot lines (characteristic times and scales)
#
#ax_phi.axhline(y=tau_rec,color='k',linestyle='-')
#ax_phi.axhline(y=tau_eq,color='k',linestyle='-')
#ax_phi.axhline(y=tau_lambda,color='k',linestyle='-')
#ax_phi.plot(aH,tau,'r-',linewidth=2)
#
# dealing with labels
#
#ax_phi.set_title(r'$\phi$')
#ax_phi.text(1.5*k[0],0.9*tau_rec,r'$\mathrm{rec.}$')
#ax_phi.text(1.5*k[0],0.9*tau_eq,r'$\mathrm{R/M} \,\, \mathrm{eq.}$')
#ax_phi.text(1.5*k[0],0.9*tau_lambda,r'$\mathrm{M/L} \,\, \mathrm{eq.}$')
#ax_phi.annotate(r'$\mathrm{Hubble} \,\, \mathrm{cross.}$',
#                  xy=(background_aH_at_tau(tau_label_Hubble),tau_label_Hubble),
#                  xytext=(0.1*background_aH_at_tau(tau_label_Hubble),0.8*tau_label_Hubble),
#                  arrowprops=dict(facecolor='black', shrink=0.05, width=1, headlength=5, headwidth=5))

#
# dealing with axes

#ax_phi.set_xscale('log')
#ax_phi.set_yscale('log')
#ax_phi.set_xlabel(r'$k \,\,\, \mathrm{[h/Mpc]}$')
#ax_phi.set_ylabel(r'$\tau \,\,\, \mathrm{[Mpc]}$')
#ax_phi.invert_yaxis()
#
# color legend
#
#fig.colorbar(fig_phi)
#
# produce the PDF
#
plt.savefig('many_times.png',dpi=300)
plt.show()

