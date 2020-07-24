
# import necessary modules
import matplotlib.pyplot as plt
import numpy as np
from classy import Class
from scipy.interpolate import interp1d


Gamma_dcdm =  np.array([98.0392,980.392, 980.392])
Gamma_Gyrs_m1 =Gamma_dcdm*1.02e-3

epsilon = np.array([0.45, 0.45, 0.4])

ax_1 = plt.subplot(211)
ax_2 = plt.subplot(212, sharex = ax_1)
plt.subplots_adjust(hspace=0)



#ax_1.set_ylim([-10,6000])
#ax_2.set_ylim([0.003,100])

ax_1.set_xlim([2,2500])
ax_2.set_xlim([2,2500])


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
#                   'non linear':'halofit',
#                   'P_k_max_1/Mpc':30.0
                   'P_k_max_1/Mpc':1.0
                   }
#Planck 2018 best-fit (TT,TE,EE+lowE+lensing)
#see table 1. in arXiv: 1807.06209v2 
M = Class()


M.set(common_settings)


for i in range(3):
    M.set(common_settings)
    M.set({
    'omega_cdm': 0.00001,
    'Omega_ini_dcdm2': 0.2636,
    'Gamma_dcdm': Gamma_dcdm[i],
    'M_dcdm': 1,
    'epsilon_dcdm': epsilon[i],
    'background_ncdm_distribution': 1,
    'Quadrature strategy': 4,
    'N_ncdm': 1,
    'evolver': 0,
    'ncdm_fluid_approximation': 2,
    'ncdm_fluid_trigger_tau_over_tau_k': 25,
    'Number of momentum bins perturbs': 300,
    'massive_daughter_perturbations': 'no',
    'dark_radiation_perturbations': 'no',
    'mother_dcdm_perturbations': 'no'
    })

    M.compute()
    if i==0:
        clM_3 = M.lensed_cl(2500)
        ll_DCDM_3 = clM_3['ell'][2:]
        clTT_DCDM_3 = clM_3['tt'][2:]
        clEE_DCDM_3 = clM_3['ee'][2:]

    elif i==1:
        clM_4 = M.lensed_cl(2500)
        ll_DCDM_4 = clM_4['ell'][2:]
        clTT_DCDM_4 = clM_4['tt'][2:]
        clEE_DCDM_4 = clM_4['ee'][2:]
    else:
        clM_5 = M.lensed_cl(2500)
        ll_DCDM_5 = clM_5['ell'][2:]
        clTT_DCDM_5 = clM_5['tt'][2:]
        clEE_DCDM_5 = clM_5['ee'][2:]

    M.struct_cleanup()
    M.empty()
    
    
T_cmb = 2.7225 #we change units for Planck
fTT_ourcode_3 = interp1d(ll_DCDM_3, clTT_DCDM_3*(ll_DCDM_3)*(ll_DCDM_3+1)/2/np.pi*(T_cmb*1.e6)**2)
fEE_ourcode_3 = interp1d(ll_DCDM_3, clEE_DCDM_3*(ll_DCDM_3)*(ll_DCDM_3+1)/2/np.pi*(T_cmb*1.e6)**2)
fTT_ourcode_4 = interp1d(ll_DCDM_4, clTT_DCDM_4*(ll_DCDM_4)*(ll_DCDM_4+1)/2/np.pi*(T_cmb*1.e6)**2)
fEE_ourcode_4 = interp1d(ll_DCDM_4, clEE_DCDM_4*(ll_DCDM_4)*(ll_DCDM_4+1)/2/np.pi*(T_cmb*1.e6)**2)
fTT_ourcode_5 = interp1d(ll_DCDM_5, clTT_DCDM_5*(ll_DCDM_5)*(ll_DCDM_5+1)/2/np.pi*(T_cmb*1.e6)**2)
fEE_ourcode_5 = interp1d(ll_DCDM_5, clEE_DCDM_5*(ll_DCDM_5)*(ll_DCDM_5+1)/2/np.pi*(T_cmb*1.e6)**2)

#%%

ax_2.tick_params(axis='both', which='minor', labelsize=12)

ax_2.set_yscale('log')

ax_1.semilogx(ll_DCDM_3,fTT_ourcode_3(ll_DCDM_3),'b',label=r'$\varepsilon = %.2f, \, \, \Gamma = %.1f \, \mathrm{Gyrs}^{-1} $'%(epsilon[0],Gamma_Gyrs_m1[0]))
ax_2.semilogx(ll_DCDM_3,fEE_ourcode_3(ll_DCDM_3),'b',label=r'$\varepsilon = %.2f \, \, \Gamma = %.1f \, \mathrm{Gyrs}^{-1}  $'%(epsilon[0],Gamma_Gyrs_m1[0]))

ax_1.semilogx(ll_DCDM_3,fTT_ourcode_4(ll_DCDM_3),'r',label=r'$\varepsilon = %.2f \, \, \Gamma = %.1f \, \mathrm{Gyrs}^{-1} $'%(epsilon[1],Gamma_Gyrs_m1[1]))
ax_2.semilogx(ll_DCDM_3,fEE_ourcode_4(ll_DCDM_3),'r',label=r'$\varepsilon = %.2f \, \, \Gamma = %.1f \, \mathrm{Gyrs}^{-1} $'%(epsilon[1],Gamma_Gyrs_m1[1]))

ax_1.semilogx(ll_DCDM_3,fTT_ourcode_5(ll_DCDM_3),'r--',label=r'$\varepsilon = %.2f \, \, \Gamma = %.1f \,  \mathrm{Gyrs}^{-1} $'%(epsilon[2],Gamma_Gyrs_m1[2]))
ax_2.semilogx(ll_DCDM_3,fEE_ourcode_5(ll_DCDM_3),'r--',label=r'$\varepsilon = %.2f \, \, \Gamma = %.1f \, \mathrm{Gyrs}^{-1} $'%(epsilon[2],Gamma_Gyrs_m1[2]))

ax_2.set_xlabel(r'$\mathrm{multipole} \, \ell$',fontsize=15)
ax_1.set_ylabel(r'$ [\ell (\ell+1) / 2 \pi] C_\ell^\mathrm{TT}$',fontsize=20)
ax_2.set_ylabel(r'$[\ell (\ell+1) / 2 \pi] C_\ell^\mathrm{EE}$',fontsize=20)

ax_1.legend(fontsize =13,loc='best')

plt.show()

plt.clf()

