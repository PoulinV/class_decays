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

Log10Gamma_dcdm = 1.2476
Gamma_dcdm=10**(Log10Gamma_dcdm )
tau =1./(Gamma_dcdm*1.02e-3)
tau

nbins = 300

log10epsilon = -2.1624
epsilon = 10**(log10epsilon)
epsilon

M_ncdm=0.27
m_ncdm=M_ncdm/3.0


zz = np.linspace(0,1.6,1000) # redshift z

fsigma8_lcdm = [] # f*sigma8
fsigma8_dcdm = [] # f*sigma8
fsigma8_mNu = [] # f*sigma8

#set general configuration
common_settings = {'output':'tCl,pCl,lCl,mPk',
                   'lensing':'yes',
                   'l_max_scalars':2600,
                   'n_s':0.9663,
                   'ln10^{10}A_s':3.045,
                   'tau_reio':0.055,
                   'omega_b':0.02242,
                   '100*theta_s':1.042059,
                   'P_k_max_h/Mpc':1.0,
                   'z_max_pk' : 5.0
                   }
  


#%% compute reference LCDM

M = Class()
#remember that base lcdm model features one massive neutrino of 0.06 eV
print("~~~~~computing reference LCDM~~~~~")
M.set(common_settings)
#we consider three massive degenerate neutrinos 
M.set({
'omega_cdm': 0.1194, #Omega_cdm 0.2626
'N_ncdm':1,
'background_ncdm_distribution': 0,
'N_ur':0.00641,
'deg_ncdm':3,
'm_ncdm':0.02})
M.compute()


h = M.h() # get reduced Hubble for conversions to 1/Mpc

#derived = M.get_current_derived_parameters(['sigma8','Omega_m'])
#print("Omega_m for LCDM is %f" %derived['Omega_m'])


for z in zz:
    fsigma8_lcdm.append(M.scale_independent_growth_factor_f(z)*M.sigma(8.0/h,z)) 
    
M.struct_cleanup()
M.empty()

timeafterref=time.time()

#%% compute best-fit DCDM
print("~~~~~time =%.f s; computing our code~~~~~"%(timeafterref-start_time))

M.set({'output':'tCl,pCl,lCl,mPk',
                   'lensing':'yes',
                   'l_max_scalars':2600,
                   'n_s':0.9673,
                   'ln10^{10}A_s':3.052,
                   'tau_reio':0.0582,
                   'omega_b':0.02240,
                   '100*theta_s':1.042174, #obtained by running class with all these parameters and H0=67.70
                   'P_k_max_h/Mpc':1.0,
                   'z_max_pk' : 2.0})
M.set({
    'omega_cdm': 0.00001,
    'omega_ini_dcdm2': 0.1194,
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
    'dark_radiation_perturbations': 'yes'
    })

M.compute()

h = M.h()
t_i = time.time()
#derived = M.get_current_derived_parameters(['sigma8','Omega_m'])
#print("~~~~~ DCDM computed in in %.f s~~~~~"%(time.time()-timeafterref))

#print("Omega_m for DCDM with epsilon=%.3f and tau= %.0f Gyrs is %f" %(epsilon,tau,derived['Omega_m']))

for z in zz:
    fsigma8_dcdm.append(M.scale_independent_growth_factor_f(z)*M.sigma(8.0/h,z)) 
    

M.struct_cleanup()
M.empty()

#%% compute LCDM+mNU

M.set(common_settings)
M.set({
'omega_cdm': 0.1154,
'N_ncdm':1,
'background_ncdm_distribution': 0,
'N_ur':0.00641,
'deg_ncdm':3,
'm_ncdm':m_ncdm
})
M.compute()

h = M.h()
#derived = M.get_current_derived_parameters(['sigma8','Omega_m'])
#print("~~~~~LCDM+mNU computed in %.f s~~~~~"%(time.time()-t_i))
#print("Omega_m for LCDM with total neutrino mass M_nu=%.2f eV is %f" %(M_ncdm,derived['Omega_m']) )

for z in zz:
    fsigma8_mNu.append(M.scale_independent_growth_factor_f(z)*M.sigma(8.0/h,z)) 
    

M.struct_cleanup()
M.empty()



#%%

print("~~~~~ready to plot~~~~~")

plt.figure(1)
plt.xlim(zz[0],zz[-1])
plt.ylim(0.2,0.7)

plt.xlabel(r'$z$', fontsize=15)
plt.ylabel(r'$f \, \sigma_8$', fontsize=20)

plt.plot(zz, fsigma8_lcdm, 'k', label=r'$\nu\Lambda$CDM Baseline')
plt.plot(zz, fsigma8_dcdm, 'k--', label=r'$\Lambda$DDM  Best-fit')
plt.plot(zz, fsigma8_mNu, 'k:', label=r'$\nu\Lambda\mathrm{CDM} \, \, (M_{\nu} =%.2f \, \mathrm{eV}) $'%M_ncdm)

#plot data points

BOSS_DR12_z = [0.38, 0.51, 0.61]
BOSS_DR12_fs8= [0.497, 0.458, 0.436]
BOSS_DR12_fs8_err = [0.045,0.038, 0.034]

DR14_quasars_z = 1.52
DR14_quasars_fs8 = 0.426
DR14_quasars_fs8_err = 0.077

WiggleZ_z = [0.44, 0.60, 0.73]
WiggleZ_fs8 = [0.413, 0.390, 0.437]
WiggleZ_fs8_err = [0.08, 0.063, 0.072]

SIXdFGS_z=0.067
SIXdFGS_fs8 = 0.423
SIXdFGS_fs8_err = 0.055

FastSound_z = 1.4
FastSound_fs8 = 0.482
FastSound_fs8_err = 0.116

GAMA_z = [0.18, 0.377]
GAMA_fs8 = [0.36, 0.44]
GAMA_fs8_err = [0.09, 0.06]

SDSS_LRG_z=0.3
SDSS_LRG_fs8 = 0.49
SDSS_LRG_fs8_err = 0.08

SDSS_MGS_z = 0.15
SDSS_MGS_fs8 = 0.49
SDSS_MGS_fs8_err = 0.14

VIPERS_z = [0.598, 0.86]
VIPERS_fs8 = [0.55, 0.4]
VIPERS_fs8_err = [0.12, 0.11]

plt.errorbar(BOSS_DR12_z, BOSS_DR12_fs8, yerr=BOSS_DR12_fs8_err, fmt='o', color='red')
plt.errorbar(DR14_quasars_z, DR14_quasars_fs8, yerr = DR14_quasars_fs8_err, fmt='o', color='orange')
plt.errorbar(WiggleZ_z, WiggleZ_fs8, yerr = WiggleZ_fs8_err, fmt='o', color='blue')
plt.errorbar(SIXdFGS_z, SIXdFGS_fs8, yerr = SIXdFGS_fs8_err, fmt='o', color='green')
plt.errorbar(FastSound_z, FastSound_fs8, yerr = FastSound_fs8_err, fmt='_', color='darkblue')
plt.errorbar(GAMA_z, GAMA_fs8, yerr = GAMA_fs8_err, fmt='v', color='darkred')
plt.errorbar(SDSS_LRG_z, SDSS_LRG_fs8, yerr = SDSS_LRG_fs8_err, fmt='^', color='darkcyan')
plt.errorbar(SDSS_MGS_z, SDSS_MGS_fs8, yerr = SDSS_MGS_fs8_err, fmt='s', color='purple')
plt.errorbar(VIPERS_z, VIPERS_fs8, yerr = VIPERS_fs8_err, fmt='D', color='olive')


plt.text(0.36, 0.55, r'BOSS DR12',color='red', fontsize =13)
plt.text(1.3, 0.33, r'DR14 quasars',color='orange', fontsize =13)
plt.text(0.5, 0.3, r'WiggleZ',color='blue', fontsize =13)
plt.text(0.01, 0.34, r'6dFGS',color='green', fontsize =13)
plt.text(1.3, 0.61, r'FastSound',color='darkblue', fontsize =13)
plt.text(0.21, 0.28, r'GAMA',color='darkred', fontsize =13)
plt.text(0.21, 0.59, r'SDSS LRG',color='darkcyan', fontsize =13)
plt.text(0.1, 0.64, r'SDSS MGS',color='purple', fontsize =13)
plt.text(0.62, 0.55, r'VIPERS',color='olive', fontsize =13)

plt.tick_params(axis="x", labelsize=18)
plt.tick_params(axis="y", labelsize=18)


plt.legend(frameon=False, loc='best', fontsize=15, borderaxespad=0.)

plt.show()
plt.clf()
