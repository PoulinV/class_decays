
# import necessary modules
#get_ipython().magic(u'matplotlib inline')
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from classy import Class
from operator import truediv
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

ax_1.set_ylim([-0.05,0.05])
ax_2.set_ylim([-0.02,0.02])
ax_1.set_xlim([2,2500])
ax_2.set_xlim([2,2500])

#Gamma_dcdm = 980
Gamma_dcdm = 4
#convert gamma to Mpc^{-1}
gamma_Mpc=Gamma_dcdm*1.0e3/2.99792458e8
#Gamma_dcdm = np.array([9.8039, 32.679, 98.0392])
tau =1./(Gamma_dcdm*1.02e-3)
tau

nbins = 300


m_dcdm = np.array([0.7, 0.8, 0.9])
epsilon =0.5*(1.-m_dcdm*m_dcdm)
epsilon



#compute characteristic scale factor at decay, t(a_d)=tau
h=0.6739
H0=h*100
Omega_b=0.02229/h**2
Omega_dcdm=0.2636
Omega_m = Omega_b+Omega_dcdm
a_d = ((3./2.)*(H0/Gamma_dcdm)*np.sqrt(Omega_m))**(2./3.)
a_d
#only if tau < age universe, otherwise a_d=1
v_kick=epsilon/(1.0-epsilon)
lambda_fss=3.0*v_kick/(gamma_Mpc*a_d) 
#formula 4.2 from Aoyama et al., works well if tau < age universe, 
#on the contrary it doesn't give good results
k_fss=(np.pi/lambda_fss)
k_fss_wdm=np.zeros(3) # for the one computed with CLASS
suppression=np.zeros(3) #for the asymptote of the matter power suppression
kk = np.logspace(-4,np.log10(1),1000) # k in h/Mpc
Pk1 = [] # P(k) in (Mpc/h)**3
Pk2 = [] # P(k) in (Mpc/h)**3
Pk3 = [] # P(k) in (Mpc/h)**3
Pk4 = []# P(k) in (Mpc/h)**3
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
'omega_cdm': 0.1197})
M.compute()
derived = M.get_current_derived_parameters(['sigma8'])
print("sigma8 for LCDM is %f"%derived['sigma8'])

clM = M.lensed_cl(2600)
ll_LCDM = clM['ell'][2:]
clTT_LCDM = clM['tt'][2:]
clEE_LCDM = clM['ee'][2:]

fTT_ref = interp1d(ll_LCDM,clTT_LCDM)
fEE_ref = interp1d(ll_LCDM,clEE_LCDM)

# get P(k) at redhsift z=0
h = M.h() # get reduced Hubble for conversions to 1/Mpc

for k in kk:
    Pk1.append(M.pk(k*h,0.)*h**3) # function .pk(k,z)

M.struct_cleanup()
M.empty()

timeafterref=time.time()
#%% compute DCDM 
t_i=timeafterref
print("~~~~~time =%.f s; computing our code~~~~~"%(timeafterref-start_time))

for i in range(3):
    M.set(common_settings)
    M.set({
    'omega_cdm': 0.00001,
    'Omega_ini_dcdm2': 0.2636,
    'Gamma_dcdm': Gamma_dcdm,
    'M_dcdm': 1,
    'm_dcdm': m_dcdm[i],
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
    h = M.h()
    derived = M.get_current_derived_parameters(['sigma8','k_fss_wdm', 'age'])
    print("~~~~~computing our code in %.f s~~~~~"%(time.time()-t_i))
    t_i = time.time()
    
    if i==0:
        clM_0 = M.lensed_cl(2500)
        ll_DCDM_0 = clM_0['ell'][2:]
        clTT_DCDM_0 = clM_0['tt'][2:]
        clEE_DCDM_0 = clM_0['ee'][2:]
        for k in kk:
            Pk2.append(M.pk(k*h,0.)*h**3) # function .pk(k,z)
        print("sigma8 for DCDM with epsilon=%.4f is %f" %(epsilon[0],derived['sigma8']) )
        k_fss_wdm[0]=derived['k_fss_wdm']
        suppression[0]=((np.exp(-derived['age']/tau))**2)-1.0
        print("k_fss_wdm for DCDM with epsilon=%.4f is %f Mpc^-1" %(epsilon[0],k_fss_wdm[0]) )
    elif i==1:
        clM_1 = M.lensed_cl(2500)
        ll_DCDM_1 = clM_1['ell'][2:]
        clTT_DCDM_1 = clM_1['tt'][2:]
        clEE_DCDM_1 = clM_1['ee'][2:]
        for k in kk:
            Pk3.append(M.pk(k*h,0.)*h**3) # function .pk(k,z)
        print("sigma8 for DCDM with epsilon=%.4f is %f" %(epsilon[1],derived['sigma8']) )
        k_fss_wdm[1]=derived['k_fss_wdm']
        suppression[1]=((np.exp(-derived['age']/tau))**2)-1.0
        print("k_fss_wdm for DCDM with epsilon=%.4f is %f Mpc^-1" %(epsilon[1],k_fss_wdm[1]) )
    else:
        clM_2 = M.lensed_cl(2500)
        ll_DCDM_2 = clM_2['ell'][2:]
        clTT_DCDM_2 = clM_2['tt'][2:]
        clEE_DCDM_2 = clM_2['ee'][2:]
        for k in kk:
            Pk4.append(M.pk(k*h,0.)*h**3) # function .pk(k,z)
        print("sigma8 for DCDM with epsilon=%.4f is %f" %(epsilon[2],derived['sigma8']) )
        k_fss_wdm[2]=derived['k_fss_wdm']
        suppression[2]=((np.exp(-derived['age']/tau))**2)-1.0
        print("k_fss_wdm for DCDM with epsilon=%.4f is %f Mpc^-1" %(epsilon[2],k_fss_wdm[2]) )
    
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

ax_1.semilogx(ll_DCDM_0,fTT_ourcode_0(ll_DCDM_0)/fTT_ref(ll_DCDM_0)-1,'b',label=r'$\varepsilon = %.4f$'%epsilon[0])
ax_2.semilogx(ll_DCDM_0,fEE_ourcode_0(ll_DCDM_0)/fEE_ref(ll_DCDM_0)-1,'b',label=r'$\varepsilon = %.4f$'%epsilon[0])
#ax_1.semilogx(ll_DCDM_0,fTT_ourcode_0(ll_DCDM_0)/fTT_ref(ll_DCDM_0)-1,'b',label=r'$\Gamma^{-1} = 100 \, \mathrm{Gyrs}$')
#ax_2.semilogx(ll_DCDM_0,fEE_ourcode_0(ll_DCDM_0)/fEE_ref(ll_DCDM_0)-1,'b',label=r'$\Gamma^{-1} = 100 \, \mathrm{Gyrs}$')
ax_1.semilogx(ll_DCDM_0,fTT_ourcode_1(ll_DCDM_0)/fTT_ref(ll_DCDM_0)-1,'r',label=r'$\varepsilon = %.4f$'%epsilon[1])
ax_2.semilogx(ll_DCDM_0,fEE_ourcode_1(ll_DCDM_0)/fEE_ref(ll_DCDM_0)-1,'r',label=r'$\varepsilon = %.4f$'%epsilon[1])
#ax_1.semilogx(ll_DCDM_0,fTT_ourcode_1(ll_DCDM_0)/fTT_ref(ll_DCDM_0)-1,'r',label=r'$\Gamma^{-1} = 30 \, \mathrm{Gyrs}$')
#ax_2.semilogx(ll_DCDM_0,fEE_ourcode_1(ll_DCDM_0)/fEE_ref(ll_DCDM_0)-1,'r',label=r'$\Gamma^{-1} = 30 \, \mathrm{Gyrs}$')
ax_1.semilogx(ll_DCDM_0,fTT_ourcode_2(ll_DCDM_0)/fTT_ref(ll_DCDM_0)-1,'g',label=r'$\varepsilon = %.4f$'%epsilon[2])
ax_2.semilogx(ll_DCDM_0,fEE_ourcode_2(ll_DCDM_0)/fEE_ref(ll_DCDM_0)-1,'g',label=r'$\varepsilon = %.4f$'%epsilon[2])
#ax_1.semilogx(ll_DCDM_0,fTT_ourcode_2(ll_DCDM_0)/fTT_ref(ll_DCDM_0)-1,'g',label=r'$\Gamma^{-1} = 10 \, \mathrm{Gyrs}$')
#ax_2.semilogx(ll_DCDM_0,fEE_ourcode_2(ll_DCDM_0)/fEE_ref(ll_DCDM_0)-1,'g',label=r'$$\Gamma^{-1} = 10 \, \mathrm{Gyrs}$')
ax_2.set_xlabel(r'$\mathrm{multipole} \, \ell$',fontsize=15)
ax_1.set_ylabel(r'$\frac{C_\ell^\mathrm{TT}(\mathrm{DCDM})}{C_\ell^\mathrm{TT}(\Lambda \mathrm{CDM} )} -1$',fontsize=20)
ax_2.set_ylabel(r'$\frac{C_\ell^\mathrm{EE}(\mathrm{DCDM})}{C_\ell^\mathrm{EE}(\Lambda \mathrm{CDM} )} -1$',fontsize=20)


ax_1.text(500, -0.04, r'$\Gamma^{-1} = 392 \, \,  \mathrm{Gyrs}$', fontsize =13)
#ax_1.text(500, -0.4, r'$\varepsilon = 0.3$', fontsize =13)
ax_1.legend(frameon=False,fontsize =13,loc='best',borderaxespad=0.)
plt.show()

plt.clf()

#%%
plt.figure(1)
plt.xscale('log')
plt.xlim(kk[0],kk[-1])
#plt.ylim(-1.5,2.0)
plt.xlabel(r'$k \,\,\,\, [h/\mathrm{Mpc}]$', fontsize=15)
plt.ylabel(r'$P_{\mathrm{DCDM}}/P_{\Lambda\mathrm{CDM}}-1$', fontsize=20)
plt.plot(kk,map(truediv, list(np.array(Pk2) - np.array(Pk1)), Pk1),'b',label=r'$\varepsilon = %.4f$'%epsilon[0])
plt.plot(kk,map(truediv, list(np.array(Pk3) - np.array(Pk1)), Pk1),'r', label=r'$\varepsilon = %.4f$'%epsilon[1])
plt.plot(kk,map(truediv, list(np.array(Pk4) - np.array(Pk1)), Pk1),'g', label = r'$\varepsilon = %.4f$'%epsilon[2])
plt.legend(loc='best', fontsize=13)
#plt.axvline(k_fss[0], color='b', linestyle='dotted')
#plt.axvline(k_fss[1], color='r', linestyle='dotted')
#plt.axvline(k_fss[2], color='g', linestyle='dotted')

plt.axvline(k_fss_wdm[0], color='b', linestyle='dotted')
plt.axvline(k_fss_wdm[1], color='r', linestyle='dotted')
plt.axvline(k_fss_wdm[2], color='g', linestyle='dotted')

#plt.axhline(suppression[0], color='b', linestyle='dotted')
#plt.axhline(suppression[1], color='r', linestyle='dotted')
#plt.axhline(suppression[2], color='g', linestyle='dotted')
plt.text(0.2e-3, -0.08, r'$\Gamma^{-1} = %.1f \, \mathrm{Gyrs} $'%tau, fontsize =13)
#plt.text(0.2e-3, -1.0, r'$\Gamma^{-1} = %.1f \, \mathrm{Gyrs} $'%tau, fontsize =13)
plt.show()

plt.clf()


