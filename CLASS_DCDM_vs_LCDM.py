
# import necessary modules
#get_ipython().magic(u'matplotlib inline')
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from classy import Class
from scipy.optimize import fsolve
import math


from scipy.interpolate import interp1d

import time
start_time = time.time()
#%% compute reference

##create plot
ax_1 = plt.subplot(211)
ax_2 = plt.subplot(212, sharex = ax_1)
plt.subplots_adjust(hspace=0)
#ax_1.set_ylim([-0.5,0.5])
#ax_2.set_ylim([-0.3,0.3])

ax_1.set_ylim([-1,1])
ax_2.set_ylim([-0.5,0.5])
ax_1.set_xlim([2,2500])
ax_2.set_xlim([2,2500])

#Gamma_dcdm = 32.679

Gamma_dcdm = np.array([9.8039, 32.679, 98.0392])
tau =1./(Gamma_dcdm*1.02e-3)
tau

nbins = 300

#m_dcdm = np.array([0.1414, 0.6324, 0.8944])

m_dcdm = 0.6324
epsilon =0.5*(1.0-m_dcdm**2)
epsilon
#%%

#set general configuration
common_settings = {'output':'tCl,pCl,lCl,mPk',
                   'lensing':'yes',
                   'l_max_scalars':2600,
                   'n_s':0.9656,
                   'A_s':2.09e-9,
                   'tau_reio':0.0536,
                   'omega_b':0.02229,
                   'h':0.6739
                   }
#Planck 2018 best-fit (TT,TE,EE+lowE+lensing)
#see table 1. in arXiv: 1807.06209v2 
M = Class()

print("~~~~~computing reference~~~~~")
M.set(common_settings)
M.set({
'omega_cdm': 0.1197})
M.compute()
clM = M.lensed_cl(2600)
ll_LCDM = clM['ell'][2:]
clTT_LCDM = clM['tt'][2:]
clEE_LCDM = clM['ee'][2:]

fTT_ref = interp1d(ll_LCDM,clTT_LCDM)
fEE_ref = interp1d(ll_LCDM,clEE_LCDM)

M.struct_cleanup()
M.empty()

timeafterref=time.time()
#%% compute DCDM 

timeafterref=time.time()
t_i=timeafterref
print("~~~~~time =%.f s; computing our code~~~~~"%(timeafterref-start_time))

for i in range(3):
    M.set(common_settings)
    M.set({
    'omega_cdm': 0.00001,
    'Omega_ini_dcdm2': 0.2636,
    'Gamma_dcdm': Gamma_dcdm[i],
    'M_dcdm': 1,
    'm_dcdm': m_dcdm,
    'background_ncdm_distribution': 1,
    'Quadrature strategy': 4,
    'N_ncdm': 1,
    'evolver': 0,
    'ncdm_fluid_approximation': 2,
    'ncdm_fluid_trigger_tau_over_tau_k': 25,
    'Number of momentum bins perturbs': nbins,
    'massive_daughter_perturbations': 'yes',
    'dark_radiation_perturbations': 'yes'
    })

    M.compute()
    print("~~~~~computing our code in %.f s~~~~~"%(time.time()-t_i))
    t_i = time.time()
    
    if i==0:
        clM_0 = M.lensed_cl(2500)
        ll_DCDM_0 = clM_0['ell'][2:]
        clTT_DCDM_0 = clM_0['tt'][2:]
        clEE_DCDM_0 = clM_0['ee'][2:]
    elif i==1:
        clM_1 = M.lensed_cl(2500)
        ll_DCDM_1 = clM_1['ell'][2:]
        clTT_DCDM_1 = clM_1['tt'][2:]
        clEE_DCDM_1 = clM_1['ee'][2:]
    else:
        clM_2 = M.lensed_cl(2500)
        ll_DCDM_2 = clM_2['ell'][2:]
        clTT_DCDM_2 = clM_2['tt'][2:]
        clEE_DCDM_2 = clM_2['ee'][2:]
    
    M.struct_cleanup()
    M.empty()


fTT_ourcode_0 = interp1d(ll_DCDM_0, clTT_DCDM_0)
fEE_ourcode_0 = interp1d(ll_DCDM_0, clEE_DCDM_0)
fTT_ourcode_1 = interp1d(ll_DCDM_1, clTT_DCDM_1)
fEE_ourcode_1 = interp1d(ll_DCDM_1, clEE_DCDM_1)
fTT_ourcode_2 = interp1d(ll_DCDM_2, clTT_DCDM_2)
fEE_ourcode_2 = interp1d(ll_DCDM_2, clEE_DCDM_2)
#%%

print("~~~~~ready to plot~~~~~")
#plot


ax_2.tick_params(axis='both', which='minor', labelsize=12)

#ax_1.semilogx(ll_DCDM_0,fTT_ourcode_0(ll_DCDM_0)/fTT_ref(ll_DCDM_0)-1,'b',label=r'$\varepsilon = 0.49$')
#ax_2.semilogx(ll_DCDM_0,fEE_ourcode_0(ll_DCDM_0)/fEE_ref(ll_DCDM_0)-1,'b',label=r'$\varepsilon = 0.49$')
ax_1.semilogx(ll_DCDM_0,fTT_ourcode_0(ll_DCDM_0)/fTT_ref(ll_DCDM_0)-1,'b',label=r'$\Gamma^{-1} = 100 \, \mathrm{Gyrs}$')
ax_2.semilogx(ll_DCDM_0,fEE_ourcode_0(ll_DCDM_0)/fEE_ref(ll_DCDM_0)-1,'b',label=r'$\Gamma^{-1} = 100 \, \mathrm{Gyrs}$')
#ax_1.semilogx(ll_DCDM_0,fTT_ourcode_1(ll_DCDM_0)/fTT_ref(ll_DCDM_0)-1,'r',label=r'$\varepsilon = 0.3$')
#ax_2.semilogx(ll_DCDM_0,fEE_ourcode_1(ll_DCDM_0)/fEE_ref(ll_DCDM_0)-1,'r',label=r'$\varepsilon = 0.3$')
ax_1.semilogx(ll_DCDM_0,fTT_ourcode_1(ll_DCDM_0)/fTT_ref(ll_DCDM_0)-1,'r',label=r'$\Gamma^{-1} = 30 \, \mathrm{Gyrs}$')
ax_2.semilogx(ll_DCDM_0,fEE_ourcode_1(ll_DCDM_0)/fEE_ref(ll_DCDM_0)-1,'r',label=r'$\Gamma^{-1} = 30 \, \mathrm{Gyrs}$')
#ax_1.semilogx(ll_DCDM_0,fTT_ourcode_2(ll_DCDM_0)/fTT_ref(ll_DCDM_0)-1,'g',label=r'$\varepsilon = 0.1$')
#ax_2.semilogx(ll_DCDM_0,fEE_ourcode_2(ll_DCDM_0)/fEE_ref(ll_DCDM_0)-1,'g',label=r'$\varepsilon = 0.1$')
ax_1.semilogx(ll_DCDM_0,fTT_ourcode_2(ll_DCDM_0)/fTT_ref(ll_DCDM_0)-1,'g',label=r'$\Gamma^{-1} = 10 \, \mathrm{Gyrs}$')
ax_2.semilogx(ll_DCDM_0,fEE_ourcode_2(ll_DCDM_0)/fEE_ref(ll_DCDM_0)-1,'g',label=r'$$\Gamma^{-1} = 10 \, \mathrm{Gyrs}$')
ax_2.set_xlabel(r'$\mathrm{multipole} \, \ell$',fontsize=15)
ax_1.set_ylabel(r'$\frac{C_\ell^\mathrm{TT}(\mathrm{DCDM})}{C_\ell^\mathrm{TT}(\Lambda \mathrm{CDM} )} -1$',fontsize=20)
ax_2.set_ylabel(r'$\frac{C_\ell^\mathrm{EE}(\mathrm{DCDM})}{C_\ell^\mathrm{EE}(\Lambda \mathrm{CDM} )} -1$',fontsize=20)


#ax_1.text(500, -0.4, r'$\Gamma^{-1} = 30 \, \,  \mathrm{Gyrs}$', fontsize =13)
ax_1.text(500, -0.4, r'$\varepsilon = 0.3$', fontsize =13)
ax_1.legend(frameon=False,fontsize =13,loc='best',borderaxespad=0.)
plt.show()

plt.clf()





