# import necessary modules
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
#from matplotlib import pyplot
import numpy as np
from classy import Class
from operator import truediv
from scipy.interpolate import interp1d

import time
start_time = time.time()

#%%

ax_1 = plt.subplot(211)
ax_2 = plt.subplot(212, sharex = ax_1)
plt.subplots_adjust(hspace=0)

ax_1.set_ylim([-0.011,0.011])
ax_2.set_ylim([-0.022,0.022])
ax_1.set_xlim([2,2500])
ax_2.set_xlim([2,2500])

Log10Gamma_dcdm = 1.24

Gamma_dcdm=10**(Log10Gamma_dcdm )
tau =1./(Gamma_dcdm*1.02e-3)
tau

nbins = 300


log10epsilon = -2.16
epsilon = 10**(log10epsilon)
epsilon

kk = np.logspace(-4,0,1000) # k in h/Mpc

Pk1 = [] # P(k) in (Mpc/h)**3
Pk2 = [] # P(k) in (Mpc/h)**3



lTT,DlTT_mean,DlTT_error_minus,DlTT_error_plus,DlTT_bestfit= np.loadtxt("error_Planck/Planck2018_errorTT.txt",unpack=True)
lEE,DlEE_mean,DlEE_error_minus,DlEE_error_plus,DlEE_bestfit= np.loadtxt("error_Planck/Planck2018_errorEE.txt",unpack=True)
lTE,DlTE_mean,DlTE_error_minus,DlTE_error_plus,DlTE_bestfit= np.loadtxt("error_Planck/Planck2018_errorTE.txt",unpack=True)


#%% read data files from reference model

#files1 = ['/home/guillermo/class_MAJORON/dcdm_bestfit_nofluid_cl_lensed.dat']
files1 = ['/Users/gfranco/cloud/output/dcdm_bestfit_nofluid_cl_lensed.dat']

data1 = []
for data_file1 in files1:
    data1.append(np.loadtxt(data_file1))
    
#files2 = ['/home/guillermo/class_MAJORON/dcdm_bestfit_nofluid_z1_pk.dat']
files2 = ['/Users/gfranco/cloud/output/dcdm_bestfit_nofluid_z1_pk.dat']

data2 = []
for data_file2 in files2:
    data2.append(np.loadtxt(data_file2))

#files3 = ['/home/guillermo/class_MAJORON/dcdm_bestfit_nofluid_z2_pk.dat']
files3 = ['/Users/gfranco/cloud/output/dcdm_bestfit_nofluid_z2_pk.dat']

data3 = []
for data_file3 in files3:
    data3.append(np.loadtxt(data_file3))   
    
cl_lens_dcdm = data1[0]
pk_z0_dcdm   = data2[0]
pk_z3_dcdm   = data3[0]

fcl_tt_dcdm_full  = interp1d(cl_lens_dcdm[:,0],cl_lens_dcdm[:,1])
fcl_ee_dcdm_full  = interp1d(cl_lens_dcdm[:,0],cl_lens_dcdm[:,2])
fcl_pp_dcdm_full  = interp1d(cl_lens_dcdm[:,0],cl_lens_dcdm[:,5])

fpk_z0_dcdm_full  = interp1d(pk_z0_dcdm[:,0],pk_z0_dcdm[:,1])
fpk_z3_dcdm_full  = interp1d(pk_z3_dcdm[:,0],pk_z3_dcdm[:,1])

#%% compute our best-fit DCDM

M = Class()

M.set({'output':'tCl,pCl,lCl,mPk',
                   'lensing':'yes',
                   'l_max_scalars':2600,
                   'n_s':0.9673,
                   'ln10^{10}A_s':3.051,
                   'tau_reio':0.0582,
                   'omega_b':0.0224,
                   'H0':67.70,
                   'P_k_max_h/Mpc':1.0,
                   'z_max_pk' : 4.0
                   })
    
M.set({
    'omega_cdm': 0.000001,
    'Omega_ini_dcdm2': 0.26051,
    'Log10_Gamma_dcdm': Log10Gamma_dcdm,
    'M_dcdm': 1,
    'log10_epsilon_dcdm': log10epsilon,
    'background_ncdm_distribution': '0,1',
    'Quadrature strategy': '0,4',
    'N_ncdm': 2,
    'N_ur':2.0328,
    'm_ncdm':'0.06,0',
    'evolver': 0,
    'ncdm_fluid_approximation': 2,
    'ncdm_fluid_trigger_tau_over_tau_k': 25,
    'Number of momentum bins perturbs': '50,300',
    'massive_daughter_perturbations': 'yes',
    'dark_radiation_perturbations': 'yes',
    })
    
M.compute()
h = M.h()
derived = M.get_current_derived_parameters(['sigma8','Omega_m'])


clM_0 = M.lensed_cl(2500)
ll_DCDM_0 = clM_0['ell'][2:]
clTT_DCDM_0 = clM_0['tt'][2:]
clEE_DCDM_0 = clM_0['ee'][2:]

for k in kk:
    Pk1.append(M.pk(k*h,0.)*h**3) # function .pk(k,z)    
    
for k in kk:
    Pk2.append(M.pk(k*h,3.)*h**3) # function .pk(k,z)
    
S8 = derived['sigma8']*np.sqrt(derived['Omega_m']/0.3)

print("sigma8 for DCDM with epsilon=%.3f and tau= %.0f Gyrs is %f" %(epsilon,tau,derived['sigma8']))
print("and Omega_m is %f, so that S8=%f "%(derived['Omega_m'],S8))

M.struct_cleanup()
M.empty()

fTT_DCDM = interp1d(ll_DCDM_0,ll_DCDM_0*(ll_DCDM_0+1.0)*clTT_DCDM_0/(2.0*np.pi))
fEE_DCDM = interp1d(ll_DCDM_0,ll_DCDM_0*(ll_DCDM_0+1.0)*clEE_DCDM_0/(2.0*np.pi))

fpkz0_DCDM = interp1d(kk, Pk1)
fpkz3_DCDM = interp1d(kk, Pk2)


#%% TIME TO PLOT

ax_2.tick_params(axis='both', which='minor', labelsize=12)

ax_1.semilogx(ll_DCDM_0,fTT_DCDM(ll_DCDM_0)/fcl_tt_dcdm_full(ll_DCDM_0)-1,'r',label=r'$\Lambda$DDM  Best-fit')
ax_2.semilogx(ll_DCDM_0,fEE_DCDM(ll_DCDM_0)/fcl_ee_dcdm_full(ll_DCDM_0)-1,'r',label=r'$\Lambda$DDM  Best-fit')

l_cosmic_variance = np.linspace(0,48,1000)
l_cosmic_variance_1 = np.linspace(0,30,1000)
l_cosmic_variance_2 = np.linspace(30,48,2)
slope =np.array([0.15,0.0343])
ax_1.fill_between(l_cosmic_variance_1, -0.15,0.15, color='lightgray' )
ax_1.fill_between(l_cosmic_variance_2, -slope, slope, color='lightgray' )

ax_1.fill_between(lTT, -(DlTT_error_plus)/DlTT_mean, +(DlTT_error_plus)/DlTT_mean, color='lightgray')

ax_2.fill_between(l_cosmic_variance, -0.18,0.18, color='lightgray' )
ax_2.fill_between(lEE, -(DlEE_error_plus)/DlEE_mean, +(DlEE_error_plus)/DlEE_mean, color = 'lightgray')

ax_2.set_xlabel(r'$\mathrm{multipole} \, \ell$',fontsize=15)
ax_1.set_ylabel(r'$\frac{C_\ell^\mathrm{TT}(\mathrm{approx} )}{C_\ell^\mathrm{TT}(\mathrm{full} )} -1$',fontsize=20)
ax_2.set_ylabel(r'$\frac{C_\ell^\mathrm{EE}(\mathrm{approx} )}{C_\ell^\mathrm{EE}(\mathrm{full} )} -1$',fontsize=20)

ax_2.tick_params(axis="x", labelsize=18)
ax_2.tick_params(axis="y", labelsize=18)
ax_1.tick_params(axis="y", labelsize=18)

ax_1.legend(frameon=False,fontsize =15,loc='upper left',borderaxespad=0.)
plt.show()

plt.clf()

#%%
#plt.figure(1)

plt.xscale('log')
plt.xlim(kk[0],kk[-1])
plt.ylim(-0.09,0.09)

plt.xlabel(r'$k \,\,\,\, [h/\mathrm{Mpc}]$', fontsize=15)
plt.ylabel(r'$\frac{P_{\mathrm{approx}}(k)}{P_{\mathrm{full}}(k)}-1$', fontsize=20)

plt.plot(kk,fpkz0_DCDM(kk)/fpk_z0_dcdm_full(kk)-1.0,'r', label=r'$\Lambda$DDM  Best-fit (z = 0)')
plt.plot(kk,fpkz3_DCDM(kk)/fpk_z3_dcdm_full(kk)-1.0,'b', label=r'$\Lambda$DDM  Best-fit (z = 3) ')

k_range_sigma8 = np.linspace(0.1,0.9,1000) #which are the exact wavenumbers probed by DES-Y1?
plt.fill_between(k_range_sigma8, -0.1,0.1, color='lightgray' )

plt.legend(loc='upper left', fontsize=13)


plt.tick_params(axis="x", labelsize=18)
plt.tick_params(axis="y", labelsize=18)


plt.show()
plt.clf()






    