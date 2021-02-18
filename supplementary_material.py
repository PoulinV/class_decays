
# import necessary modules
import matplotlib.pyplot as plt
import numpy as np
from classy import Class
from scipy.interpolate import interp1d

#%%

Log10Gamma_dcdm = 2.991
Gamma_dcdm=10**(Log10Gamma_dcdm)
Gamma_dcdm
tau =1./(Gamma_dcdm*1.02e-3)
tau


log10epsilon = -0.30103
epsilon = 10**(log10epsilon)
epsilon

#set general configuration
common_settings = {'output':'tCl,pCl,lCl,mPk',
                   'lensing':'yes',
                   'l_max_scalars':2600,
                   'n_s':0.9663,
                   'ln10^{10}A_s':3.045,
                   'tau_reio':0.05447,
                   'omega_b':0.02242,
                   'h':0.6737,
                   'omega_cdm': 0.000001,
                   'P_k_max_1/Mpc':1.0
                   }



kk = np.logspace(-4,np.log10(1),1000) 

Pk1 = [] # P(k) in (Mpc/h)**3
Pk2 = [] # P(k) in (Mpc/h)**3
Pk3 = [] # P(k) in (Mpc/h)**3


#%% compute reference model (DCDM to DR, standard CLASS implementation)
M = Class()

M.set(common_settings)
M.set({
'Omega_ini_dcdm': 0.2636,
'N_ncdm':1,
'background_ncdm_distribution': 0,
'N_ur':2.0328,
'm_ncdm':0.06,
'Gamma_dcdm': Gamma_dcdm,
'evolver': 0,
'dark_radiation_perturbations': 'yes',
'input_verbose': 1,
'background_verbose' :1,
'thermodynamics_verbose' :1,
'perturbations_verbose': 1,
'transfer_verbose': 1,
'primordial_verbose': 1,
'spectra_verbose': 1,
'nonlinear_verbose': 1,
'lensing_verbose' :1,
'output_verbose': 1
})

M.compute()

#derived = M.get_current_derived_parameters(['sigma8','Omega_m','z_eq'])
#S8 = derived['sigma8']*np.sqrt(derived['Omega_m']/0.3)

#print("sigma8 for LCDM is %f"%derived['sigma8'])
#print("and Omega_m is %f, so that S8=%f "%(derived['Omega_m'],S8))

# get P(k) at redhsift z=0
h = M.h() # get reduced Hubble for conversions to 1/Mpc

for k in kk:
    Pk1.append(M.pk(k*h,0.)*h**3) # function .pk(k,z)
  
M.struct_cleanup()
M.empty()

fPk_ref = interp1d(kk,Pk1)



#%% compute two-body decay in limit epsilon->0.5, using fluid approx.
# with and without setting shear equation to zero. 

#WITH NON-ZERO SHEAR

M.set(common_settings)
M.set({
    'Omega_ini_dcdm2': 0.2636,
    'Gamma_dcdm': Gamma_dcdm,
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
    'switch_off_shear_wdm': 'no',
    'input_verbose': 1,
    'background_verbose' :1,
    'thermodynamics_verbose' :1,
    'perturbations_verbose': 1,
    'transfer_verbose': 1,
    'primordial_verbose': 1,
    'spectra_verbose': 1,
    'nonlinear_verbose': 1,
    'lensing_verbose' :1,
    'output_verbose': 1
    })
    
    
M.compute()
h = M.h()


for k in kk:
    Pk2.append(M.pk(k*h,0.)*h**3) # function .pk(k,z)
  

M.struct_cleanup()
M.empty()

fPk_sigYES = interp1d(kk,Pk2)

#WITH ZERO SHEAR

M.set(common_settings)
M.set({
    'Omega_ini_dcdm2': 0.2636,
    'Gamma_dcdm':Gamma_dcdm,
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
    'switch_off_shear_wdm': 'yes',
    'input_verbose': 1,
    'background_verbose' :1,
    'thermodynamics_verbose' :1,
    'perturbations_verbose': 1,
    'transfer_verbose': 1,
    'primordial_verbose': 1,
    'spectra_verbose': 1,
    'nonlinear_verbose': 1,
    'lensing_verbose' :1,
    'output_verbose': 1
    })
    


M.compute()
h = M.h()


for k in kk:
    Pk3.append(M.pk(k*h,0.)*h**3) # function .pk(k,z)
  

M.struct_cleanup()
M.empty()

fPk_sigNOT = interp1d(kk,Pk3)


#%% START TO PLOT

plt.xscale('log')
#plt.yscale('log')
plt.xlim(kk[0],kk[-1])


plt.tick_params(axis="x", labelsize=18)
plt.tick_params(axis="y", labelsize=18)

plt.xlabel(r'$k \,\,\,\, [h/\mathrm{Mpc}]$', fontsize=20)
plt.ylabel(r'$P(k)_{\mathrm{approx}}/P(k)_{\mathrm{ref}}-1$', fontsize=20)

plt.title(r'$ \varepsilon=0.5 \, \, \, \, \, \, \, \Gamma^{-1}= %0.f \ \mathrm{Gyr}$'%tau,fontsize=20)

plt.plot(kk,fPk_sigYES(kk)/fPk_ref(kk)-1.0,'b',label=r'$\dot{\sigma}_{\mathrm{wdm}} \neq 0$')
plt.plot(kk,fPk_sigNOT(kk)/fPk_ref(kk)-1.0,'r',label=r'$\dot{\sigma}_{\mathrm{wdm}} =  0$')


plt.legend(loc='best', fontsize=13)

plt.show()
plt.clf()












