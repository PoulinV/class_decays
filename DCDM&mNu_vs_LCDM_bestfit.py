
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

ax_1.set_ylim([-0.04,0.04])
ax_2.set_ylim([-0.17,0.17])
ax_1.set_xlim([2,2500])
ax_2.set_xlim([2,2500])


#Log10Gamma_dcdm = 0.89
#Log10Gamma_dcdm = 1.658
Log10Gamma_dcdm = 1.4

Gamma_dcdm=10**(Log10Gamma_dcdm )
tau =1./(Gamma_dcdm*1.02e-3)
tau

nbins = 300

#log10epsilon = -2.09
#log10epsilon = -2.529
log10epsilon = -2.3
epsilon = 10**(log10epsilon)
epsilon

#M_ncdm=0.27
M_ncdm=0.9
m_ncdm=M_ncdm/3.0

kk = np.logspace(-4,0,1000) # k in h/Mpc

Pk1 = [] # P(k) in (Mpc/h)**3
Pk2 = [] # P(k) in (Mpc/h)**3
Pk3 = [] # P(k) in (Mpc/h)**3

Pk4 = [] # P(k) in (Mpc/h)**3
Pk5 = [] # P(k) in (Mpc/h)**3
Pk6 = [] # P(k) in (Mpc/h)**3

#set general configuration
common_settings = {'output':'tCl,pCl,lCl,mPk',
                   'lensing':'yes',
                   'l_max_scalars':2600,
                   'n_s':0.9666,
                   'ln10^{10}A_s':3.047,
                   'tau_reio':0.0559,
                   'omega_b':0.02238,
                   'P_k_max_1/Mpc':1.0,
                   'z_max_pk' : 5.0
                   }

lTT,DlTT_mean,DlTT_error_minus,DlTT_error_plus,DlTT_bestfit= np.loadtxt("error_Planck/Planck2018_errorTT.txt",unpack=True)
lEE,DlEE_mean,DlEE_error_minus,DlEE_error_plus,DlEE_bestfit= np.loadtxt("error_Planck/Planck2018_errorEE.txt",unpack=True)
lTE,DlTE_mean,DlTE_error_minus,DlTE_error_plus,DlTE_bestfit= np.loadtxt("error_Planck/Planck2018_errorTE.txt",unpack=True)

#%% compute reference LCDM
#(actually, it's LCDM+free neutrino mass, in runs with Planck, BAO, SN but no prior on S8)
M = Class()


#ONE IN WHICH WE FIX THETA_S
print("~~~~~computing reference LCDM~~~~~")
M.set(common_settings)
#we consider three massive degenerate neutrinos 
M.set({
'omega_cdm': 0.1196, #Omega_cdm 0.2626
'100*theta_s':1.04187,
'N_ncdm':1,
'background_ncdm_distribution': 0,
'N_ur':0.00641,
'deg_ncdm':3,
'm_ncdm':0.02})
M.compute()

derived = M.get_current_derived_parameters(['sigma8','Omega_m','z_eq'])
print("sigma8 for LCDM is %f"%derived['sigma8'])
print("and Omega_m is %f"%derived['Omega_m'])
print("and z_eq is %f"%derived['z_eq'])

clM = M.lensed_cl(2600)
ll_LCDM = clM['ell'][2:]
clTT_LCDM = clM['tt'][2:]
clEE_LCDM = clM['ee'][2:]

fTT_ref = interp1d(ll_LCDM,clTT_LCDM)
fEE_ref = interp1d(ll_LCDM,clEE_LCDM)

M.struct_cleanup()
M.empty()

#AND OTHER IN WHICH WE FIX H0 ############################################

M.set(common_settings)
#we consider three massive degenerate neutrinos 
M.set({
'omega_cdm': 0.1196, #Omega_cdm 0.2626
'H0':67.58, #value corresponding to the reference values above (from 100*theta_s using shooting method)
'N_ncdm':1,
'background_ncdm_distribution': 0,
'N_ur':0.00641,
'deg_ncdm':3,
'm_ncdm':0.02})
M.compute()

# get P(k) at redhsift z=0
h = M.h() # get reduced Hubble for conversions to 1/Mpc

for k in kk:
    Pk1.append(M.pk(k*h,0.)*h**3) # function .pk(k,z)
    
for k in kk:
    Pk4.append(M.pk(k*h,3.)*h**3) # function .pk(k,z)

M.struct_cleanup()
M.empty()

timeafterref=time.time()


#%% compute our best-fit DCDM

t_i=timeafterref
print("~~~~~time =%.f s; computing our code~~~~~"%(timeafterref-start_time))
M.set({'output':'tCl,pCl,lCl,mPk',
                   'lensing':'yes',
                   'l_max_scalars':2600,
                   'n_s':0.9660,
                   'ln10^{10}A_s':3.048,
                   'tau_reio':0.0543,
                   'omega_b':0.02241,
                   'H0':67.56,
                   'P_k_max_1/Mpc':1.0,
                   'z_max_pk' : 5.0})
    
M.set({
    'omega_cdm': 0.00001,
    'omega_ini_dcdm2': 0.1198,
#    'Omega_ini_dcdm2': 0.2625,
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
derived = M.get_current_derived_parameters(['sigma8','Omega_m','z_eq'])
print("~~~~~computing DCDM in %.f s~~~~~"%(time.time()-t_i))
t_i = time.time()

clM_0 = M.lensed_cl(2500)
ll_DCDM_0 = clM_0['ell'][2:]
clTT_DCDM_0 = clM_0['tt'][2:]
clEE_DCDM_0 = clM_0['ee'][2:]

#For the matter power spectrum of the DCDM, should I also fix h? or it doesn't matter
# since H0 is very close to the reference value?
for k in kk:
    Pk2.append(M.pk(k*h,0.)*h**3) # function .pk(k,z)    
    
for k in kk:
    Pk5.append(M.pk(k*h,3.)*h**3) # function .pk(k,z)

print("sigma8 for DCDM with epsilon=%.3f and tau= %.0f Gyrs is %f" %(epsilon,tau,derived['sigma8']))
print("and Omega_m is %f"%derived['Omega_m'])
print("and z_eq is %f"%derived['z_eq'])

M.struct_cleanup()
M.empty()


fTT_DCDM = interp1d(ll_DCDM_0, clTT_DCDM_0)
fEE_DCDM = interp1d(ll_DCDM_0, clEE_DCDM_0)

#%% compute LCDM+neutrino_mass


# FIXING THETA_S ###################################################
print("~~~~~computing LCDM+nu~~~~~")
M.set(common_settings)
M.set({
#'omega_cdm': 0.1155,
'omega_cdm':  0.1196,
'100*theta_s':1.04187,
'N_ncdm':1,
'background_ncdm_distribution': 0,
'N_ur':0.00641,
'deg_ncdm':3,
'm_ncdm':m_ncdm
})
M.compute()

derived = M.get_current_derived_parameters(['sigma8','Omega_m','z_eq'])
clM_1 = M.lensed_cl(2500)
ll_DCDM_1 = clM_1['ell'][2:]
clTT_DCDM_1 = clM_1['tt'][2:]
clEE_DCDM_1 = clM_1['ee'][2:]

print("sigma8 for LCDM (fixed theta_s) with total neutrino mass M_nu=%.1f eV is %f" %(M_ncdm,derived['sigma8']) )
print("and Omega_m is %f"%derived['Omega_m'])
print("and z_eq is %f"%derived['z_eq'])

S8 = derived['sigma8']*np.sqrt(derived['Omega_m']/0.3)

M.struct_cleanup()
M.empty()


fTT_mNU = interp1d(ll_DCDM_1, clTT_DCDM_1)
fEE_mNU = interp1d(ll_DCDM_1, clEE_DCDM_1)


# FIXING H0 ##################################################################

M.set(common_settings)
M.set({
#'omega_cdm': 0.1155,
'omega_cdm':  0.1196,
'H0':67.58,
'N_ncdm':1,
'background_ncdm_distribution': 0,
'N_ur':0.00641,
'deg_ncdm':3,
'm_ncdm':m_ncdm
})
M.compute()

h = M.h()

print("~~~~~computing LCDM+mNU in %.f s~~~~~"%(time.time()-t_i))

for k in kk:
    Pk3.append(M.pk(k*h,0.)*h**3) # function .pk(k,z)    
    
for k in kk:
    Pk6.append(M.pk(k*h,3.)*h**3) # function .pk(k,z)

derived = M.get_current_derived_parameters(['sigma8','Omega_m','z_eq'])

print("sigma8 for LCDM (fixed H0) with total neutrino mass M_nu=%.1f eV is %f" %(M_ncdm,derived['sigma8']) )
print("and Omega_m is %f"%derived['Omega_m'])
print("and z_eq is %f"%derived['z_eq'])
#S8 = derived['sigma8']*np.sqrt(derived['Omega_m']/0.3)

M.struct_cleanup()
M.empty()

#%% TIME TO PLOT
print("~~~~~ready to plot~~~~~")
ax_2.tick_params(axis='both', which='minor', labelsize=12)

ax_1.semilogx(ll_DCDM_0,fTT_DCDM(ll_DCDM_0)/fTT_ref(ll_DCDM_0)-1,'r',label=r'$\Lambda$DDM  Best-fit')
ax_2.semilogx(ll_DCDM_0,fEE_DCDM(ll_DCDM_0)/fEE_ref(ll_DCDM_0)-1,'r',label=r'$\Lambda$DDM  Best-fit')

ax_1.semilogx(ll_DCDM_0,fTT_mNU(ll_DCDM_0)/fTT_ref(ll_DCDM_0)-1,'b',label= r'$\Lambda\mathrm{CDM}+M_{\nu} =%.1f \, \mathrm{eV} $'%M_ncdm)
ax_2.semilogx(ll_DCDM_0,fEE_mNU(ll_DCDM_0)/fEE_ref(ll_DCDM_0)-1,'b',label= r'$\Lambda\mathrm{CDM}+M_{\nu} =%.1f \, \mathrm{eV} $'%M_ncdm)

#Planck error bars
#ax_1.errorbar(lTT, DlTT_mean/DlTT_mean-1, yerr=(DlTT_error_plus)/DlTT_mean, fmt='.',color='r')

l_cosmic_variance = np.linspace(0,48,1000)
l_cosmic_variance_1 = np.linspace(0,30,1000)
l_cosmic_variance_2 = np.linspace(30,48,2)
slope =np.array([0.15,0.0343])
ax_1.fill_between(l_cosmic_variance_1, -0.15,0.15, color='lightgray' )
ax_1.fill_between(l_cosmic_variance_2, -slope, slope, color='lightgray' )

ax_1.fill_between(lTT, -(DlTT_error_plus)/DlTT_mean, +(DlTT_error_plus)/DlTT_mean, color='lightgray')

ax_2.fill_between(l_cosmic_variance, -0.18,0.18, color='lightgray' )
#ax_2.errorbar(lEE, DlEE_mean/DlEE_mean-1, yerr=DlEE_error_plus/DlEE_mean, fmt='.',color='r')
ax_2.fill_between(lEE, -(DlEE_error_plus)/DlEE_mean, +(DlEE_error_plus)/DlEE_mean, color = 'lightgray')

ax_2.set_xlabel(r'$\mathrm{multipole} \, \ell$',fontsize=15)
ax_1.set_ylabel(r'$\frac{C_\ell^\mathrm{TT}}{C_\ell^\mathrm{TT}(\mathrm{ref} )} -1$',fontsize=20)
ax_2.set_ylabel(r'$\frac{C_\ell^\mathrm{EE}}{C_\ell^\mathrm{EE}(\mathrm{ref} )} -1$',fontsize=20)


ax_2.tick_params(axis="x", labelsize=18)
ax_2.tick_params(axis="y", labelsize=18)
ax_1.tick_params(axis="y", labelsize=18)


ax_1.legend(frameon=False,fontsize =15,loc='upper left',borderaxespad=0.)
plt.show()

plt.clf()

#%%
#plt.figure(1)

fig, axes = plt.subplots()

axes.set_xscale('log')
axes.set_xlim(kk[0],kk[-1])
axes.set_ylim(-0.5,0.05)


axes.set_xlabel(r'$k \,\,\,\, [h/\mathrm{Mpc}]$', fontsize=15)
axes.set_ylabel(r'$\frac{P(k)}{P_{ref}(k)}-1$', fontsize=20)

axes.plot(kk,map(truediv, list(np.array(Pk2) - np.array(Pk1)), Pk1),'r')
axes.plot(kk,map(truediv, list(np.array(Pk3) - np.array(Pk1)), Pk1),'b')

axes.plot(kk,map(truediv, list(np.array(Pk5) - np.array(Pk4)), Pk4),'r--')
axes.plot(kk,map(truediv, list(np.array(Pk6) - np.array(Pk4)), Pk4),'b--')


lines = axes.get_lines()

black_line1 = mlines.Line2D([], [], color='black', linestyle='-', label='z = 0')
black_line2 = mlines.Line2D([], [], color='black', linestyle='--', label='z = 3')

legend1 = plt.legend([lines[i] for i in [0,1]], [r'$\Lambda$DDM  Best-fit', r'$\Lambda\mathrm{CDM}+M_{\nu} =%.1f \, \mathrm{eV} $'%M_ncdm], loc='lower left', fontsize=15, frameon=False)
legend2 = plt.legend(handles= [black_line1,black_line2], loc=6, fontsize=15, frameon=False)

axes.add_artist(legend1)
axes.add_artist(legend2)


k_range_desy1 = np.linspace(0.09,0.62,1000) #which are the exact wavenumbers probed by DES-Y1?
plt.fill_between(k_range_desy1, -0.5,0.5, color='lightgray' )

#plt.text(2.0e-4, -0.26, r'$\Omega_{m} =0.31 \, \, \, \, \, \, \, \, \, \, \, \, \sigma_{8} = 0.75  $', fontsize =15)

plt.text(5.0e-3, -0.26, r'$ S_{8} = %.3f  $'%S8, fontsize =15)

# Q: SO, EVEN IF WE ARE PLOTTING RESIDUALS IN P(K) FOR FIXED H0, THE S8 VALUE WE SHOW
# IS FROM THE CALCULATION AT FIXED 100*THETA_S, IS THAT RIGHT?
#ACTUALLY THE S8 FROM THE CALCULATION AT FIXED H0 IS SMALLER THAN S8=0.75, IS IT FINE
# IF S8=0.75 IS NOT REPRESENTING THE POWER SUPPRESSION WE ARE SHOWING?

plt.tick_params(axis="x", labelsize=18)
plt.tick_params(axis="y", labelsize=18)

#axes.set_position([left, bottom, width, length]) 
axes.set_position([0.13,0.13,0.75,0.8])  

plt.show()
plt.clf()



