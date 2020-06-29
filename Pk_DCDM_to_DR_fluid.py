
# import necessary modules
import matplotlib.pyplot as plt
import numpy as np
from classy import Class
from operator import truediv

import time
start_time = time.time()

#%%
Gamma_dcdm = 4

tau =1./(Gamma_dcdm*1.02e-3)
tau

m_dcdm = 0.9
epsilon =0.5*(1.-m_dcdm*m_dcdm)
epsilon

kk = np.logspace(-4,np.log10(1),1000) # k in h/Mpc

Pk1 = [] # P(k) in (Mpc/h)**3
Pk2 = [] # P(k) in (Mpc/h)**3

#%%
#set general configuration
common_settings = {'output':'tCl,pCl,lCl,mPk',
                   'lensing':'yes',
                   'l_max_scalars':2600,
                   'n_s':0.9656,
                   'A_s':2.09e-9,
                   'tau_reio':0.0536,
                   'omega_b':0.02229,
                   'h':0.6739,
                   'P_k_max_1/Mpc':1.0
                   }
#Planck 2018 best-fit (TT,TE,EE+lowE+lensing)
#see table 1. in arXiv: 1807.06209v2 
M = Class()
print("~~~~~computing reference~~~~~")
M.set(common_settings)
M.set({
'omega_cdm': 0.000001,
'omega_ini_dcdm2':  0.1197,
'Gamma_dcdm': Gamma_dcdm,
'M_dcdm': 1,
'm_dcdm': m_dcdm,
'background_ncdm_distribution': 1,
'Quadrature strategy': 4,
'N_ncdm': 1,
'evolver': 0,
'ncdm_fluid_approximation': 3,
'Number of momentum bins perturbs': 5000,
'l_max_ncdm':49,
'massive_daughter_perturbations': 'yes',
'dark_radiation_perturbations': 'yes'
})
M.compute()
derived = M.get_current_derived_parameters(['sigma8'])
print("sigma8 for two-body DCDM (full hierarchy with epsilon =%.3f ) is %f"%(epsilon, derived['sigma8']))
# get P(k) at redhsift z=0
h = M.h() # get reduced Hubble for conversions to 1/Mpc

for k in kk:
    Pk1.append(M.pk(k*h,0.)*h**3) # function .pk(k,z)

M.struct_cleanup()
M.empty()

timeafterref=time.time()
#%%

print("~~~~~time =%.f s; computing our code~~~~~"%(timeafterref-start_time))
M.set(common_settings)
M.set({
'omega_cdm': 0.000001,
'omega_ini_dcdm2':  0.1197,
'Gamma_dcdm': Gamma_dcdm,
'M_dcdm': 1,
'm_dcdm': m_dcdm,
'background_ncdm_distribution': 1,
'Quadrature strategy': 4,
'N_ncdm': 1,
'evolver': 0,
'ncdm_fluid_approximation': 2,
'ncdm_fluid_trigger_tau_over_tau_k': 25,
'Number of momentum bins perturbs': 300,
'massive_daughter_perturbations': 'yes',
'dark_radiation_perturbations': 'yes'
 })
    
M.compute()
h = M.h()
derived = M.get_current_derived_parameters(['sigma8'])
print("~~~~~computing our code in %.f s~~~~~"%(time.time()-timeafterref))

for k in kk:
    Pk2.append(M.pk(k*h,0.)*h**3) # function .pk(k,z)

print("sigma8 for two-body DCDM (fluid WDM with epsilon=%.3f) is %f"%(epsilon, derived['sigma8']))

M.struct_cleanup()
M.empty()

#%%
print("~~~~~ready to plot~~~~~")
plt.figure(1)
plt.xscale('log')
plt.xlim(kk[0],kk[-1])
plt.xlabel(r'$k \,\,\,\, [h/\mathrm{Mpc}]$', fontsize=15)
plt.ylabel(r'$P_{\mathrm{approx}}/P_{\mathrm{full}}-1$', fontsize=20)
plt.plot(kk,map(truediv, list(np.array(Pk2) - np.array(Pk1)), Pk1),'b',label=r'$\tau = %.3f \, \mathrm{Gyrs} \, \, \, \varepsilon = %.3f $'%(tau, epsilon))

plt.legend(loc='best', fontsize=13)

plt.show()

plt.clf()


# for very long lifetimes (several hundreds of gyrs, the ones allowed by Planck),
# we get a reasonable agreement in sigma8 between fluid and full hierarchy
# ---> this is for epsilon=0.5 (decay to DR), we should check for smaller values of epsilon

