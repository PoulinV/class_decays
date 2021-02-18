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

Log10Gamma_dcdm = 0.89
Gamma_dcdm=10**(Log10Gamma_dcdm )
tau =1./(Gamma_dcdm*1.02e-3)
tau

nbins = 300

log10epsilon = -2.09
epsilon = 10**(log10epsilon)
epsilon

m_ncdm=0.29

zz = np.linspace(0.01,2.8,1000) # redshift z

bao_scale_lcdm=[]

bao_scale_dcdm=[]
bao_scale_mNu=[]

#set general configuration
common_settings = {'output':'tCl,pCl,lCl,mPk',
                   'lensing':'yes',
                   'l_max_scalars':2600,
                   'n_s':0.9652,
                   'ln10^{10}A_s':3.043,
                   'tau_reio':0.0540,
                   'omega_b':0.02233,
                   'h':0.6737
                   }    


#%% compute reference LCDM

M = Class()
#remember that base lcdm model features one massive neutrino of 0.06 eV
print("~~~~~computing reference LCDM~~~~~")
M.set(common_settings)
M.set({
'omega_cdm': 0.1198,
'N_ncdm':1,
'background_ncdm_distribution': 0,
'N_ur':2.0328,
'm_ncdm':0.06})
M.compute()

for z in zz:
    da = M.angular_distance(z)
    dr = z/M.Hubble(z)
    dv = (da*da*(1.0+z)*(1.0+z)*dr)**(1./3.)
    bao_scale_lcdm.append(dv/M.rs_drag()) 



M.struct_cleanup()
M.empty()



fbao_scale_lcdm = interp1d(zz, bao_scale_lcdm)

timeafterref=time.time()

#%% compute best-fit DCDM
print("~~~~~time =%.f s; computing our code~~~~~"%(timeafterref-start_time))

M.set({'output':'tCl,pCl,lCl,mPk',
                   'lensing':'yes',
                   'l_max_scalars':2600,
                   'n_s':0.9671,
                   'ln10^{10}A_s':3.048,
                   'tau_reio':0.0562,
                   'omega_b':0.02238,
                   'h':0.6765})
    
M.set({
    'omega_cdm': 0.00001,
    'Omega_ini_dcdm2': 0.2611,
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

print("~~~~~ DCDM computed in in %.f s~~~~~"%(time.time()-timeafterref))
t_i = time.time()

for z in zz:
    da = M.angular_distance(z)
    dr = z/M.Hubble(z)
    dv = (da*da*(1.0+z)*(1.0+z)*dr)**(1./3.)
    bao_scale_dcdm.append(dv/M.rs_drag()) 
    
M.struct_cleanup()
M.empty()


#%% compute LCDM+mNU

M.set(common_settings)
M.set({
'Omega_cdm': 0.258,
'N_ncdm':1,
'background_ncdm_distribution': 0,
'N_ur':2.0328,
'm_ncdm':m_ncdm
})
M.compute()

print("~~~~~LCDM+mNU computed in %.f s~~~~~"%(time.time()-t_i))

for z in zz:
    da = M.angular_distance(z)
    dr = z/M.Hubble(z)
    dv = (da*da*(1.0+z)*(1.0+z)*dr)**(1./3.)
    bao_scale_mNu.append(dv/M.rs_drag()) 


M.struct_cleanup()
M.empty()


#%%

print("~~~~~ready to plot~~~~~")

plt.figure(1)
plt.xlim(zz[0],zz[-1])
plt.ylim(0.98,1.02)

plt.xlabel(r'$z$', fontsize=15)
plt.ylabel(r'$(D_V / r_{\rm drag})/(D_V / r_{\rm drag})_{\rm Planck}$', fontsize=20) 

#plt.plot(zz, bao_scale_lcdm, 'k', label=r'$ 2018 \, \, \, \mathrm{Best-fit} \, \, \, \Lambda \mathrm{CDM} $')

plt.plot(zz, map(truediv, bao_scale_dcdm, bao_scale_lcdm), 'g--', label=r'$\Lambda\mathrm{DDM} \, \, \, \mathrm{Best-fit}$')
plt.plot(zz, map(truediv, bao_scale_mNu, bao_scale_lcdm), 'b:', label=r'$\Lambda\mathrm{CDM}+m_{\nu} =%.2f \, \mathrm{eV} $'%m_ncdm)



rd_fid_bossdr12=147.78
BOSS_DR12_z = [0.38, 0.51, 0.61]
BOSS_DR12_bao= [(1477/rd_fid_bossdr12)/fbao_scale_lcdm(0.38), (1877/rd_fid_bossdr12)/fbao_scale_lcdm(0.51),(2140/rd_fid_bossdr12)/fbao_scale_lcdm(0.61) ]
BOSS_DR12_bao_err=[(16/rd_fid_bossdr12)/fbao_scale_lcdm(0.38),(19/rd_fid_bossdr12)/fbao_scale_lcdm(0.51), (22/rd_fid_bossdr12)/fbao_scale_lcdm(0.61)]

plt.errorbar(BOSS_DR12_z, BOSS_DR12_bao, yerr=BOSS_DR12_bao_err, fmt='^', color='red')

plt.text(0.36, 0.96, r'BOSS DR12',color='red', fontsize =13)

plt.tick_params(axis="x", labelsize=18)
plt.tick_params(axis="y", labelsize=18)


plt.axhline(1.0, color='k')

plt.legend(frameon=False, loc='best', fontsize=15, borderaxespad=0.)

plt.show()
plt.clf()






