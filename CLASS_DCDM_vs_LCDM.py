# import necessary modules
import matplotlib.pyplot as plt
import numpy as np
from classy import Class
from scipy.interpolate import interp1d

plt.rcParams["figure.figsize"] = [4.0,8.0]

import time
start_time = time.time()

plot_type = 1

#lTT,DlTT_mean,DlTT_error_minus,DlTT_error_plus,DlTT_bestfit= np.loadtxt("error_Planck/Planck2018_errorTT.txt",unpack=True)
#lEE,DlEE_mean,DlEE_error_minus,DlEE_error_plus,DlEE_bestfit= np.loadtxt("error_Planck/Planck2018_errorEE.txt",unpack=True)
#lTE,DlTE_mean,DlTE_error_minus,DlTE_error_plus,DlTE_bestfit= np.loadtxt("error_Planck/Planck2018_errorTE.txt",unpack=True)


#%% compute reference


if plot_type ==1:
    Gamma_dcdm = np.array([32.679, 32.679, 32.679,32.679])
    epsilon = np.array([0.4999, 0.1, 0.01, 0.001])

else:
    
    Gamma_dcdm = np.array([3.2679,9.8039, 32.679, 98.0392])
    epsilon =  np.array([0.1,0.1,0.1, 0.1])



tau =1./(Gamma_dcdm*1.02e-3)
tau

nbins = 300

k_fss_wdm=np.zeros(4) # for the one computed with CLASS
suppression=np.zeros(4) #for the asymptote of the matter power suppression
kk = np.logspace(-4,0,1000) # k in h/Mpc

Pk1 = [] # P(k) in (Mpc/h)**3
Pk2 = [] # P(k) in (Mpc/h)**3
Pk3 = [] # P(k) in (Mpc/h)**3
Pk4 = []# P(k) in (Mpc/h)**3
Pk5 = []# P(k) in (Mpc/h)**3

#Pk6 = [] # P(k) in (Mpc/h)**3
#Pk7 = []# P(k) in (Mpc/h)**3
#Pk8 = []# P(k) in (Mpc/h)**3
#Pk9 = []# P(k) in (Mpc/h)**3
#Pk10 = []# P(k) in (Mpc/h)**3



#zhigh=2.

#%%

#set general configuration
common_settings = {'output':'tCl,pCl,lCl, mPk',
                   'lensing':'yes',
                   'l_max_scalars':2600,
                   'n_s':0.9673,
                   'ln10^{10}A_s':3.052,
                   'tau_reio':0.0582,
                   'omega_b': 0.0224,
                   '100*theta_s':1.042168,
                   'P_k_max_h/Mpc':1.0
 #                  'z_max_pk' : 4.0
                   }

#Before I had Planck 2018 best-fit (TT,TE,EE+lowE+lensing)
#see table 1. in arXiv: 1807.06209v2 
#BUT NOW I HAVE BEST-FIT VALUES FROM COMBINED ANALYSIS OF OUR LETTER 

M = Class()

print("~~~~~computing reference~~~~~")
M.set(common_settings)
M.set({
'omega_cdm': 0.1194,
'N_ncdm':1,
'background_ncdm_distribution': 0,
'N_ur':2.0328,
'm_ncdm':0.06})
M.compute()
derived = M.get_current_derived_parameters(['sigma8','Omega_m'])
#print("sigma8 for LCDM is %f"%derived['sigma8'])
#print("Omega_m for LCDM is %f"%derived['Omega_m'])

clM = M.lensed_cl(2600)
#clM = M.raw_cl(2600)
ll_LCDM = clM['ell'][2:]
clTT_LCDM = clM['tt'][2:]
clEE_LCDM = clM['ee'][2:]
clphiphi_LCDM = clM['pp'][2:]


fTT_ref = interp1d(ll_LCDM,clTT_LCDM)
fEE_ref = interp1d(ll_LCDM,clEE_LCDM)
fphiphi_ref = interp1d(ll_LCDM,clphiphi_LCDM)


# get P(k) at redhsift z=0
h = M.h() # get reduced Hubble for conversions to 1/Mpc

for k in kk:
    Pk1.append(M.pk(k*h,0.)*h**3) # function .pk(k,z)
#    Pk6.append(M.pk(k*h,zhigh)*h**3) # function .pk(k,z)


M.struct_cleanup()
M.empty()

fpk_ref = interp1d(kk, Pk1)
#fpk_ref2 = interp1d(kk, Pk6)


timeafterref=time.time()
#%% compute DCDM 
t_i=timeafterref
print("~~~~~time =%.f s; computing our code~~~~~"%(timeafterref-start_time))

for i in range(4):
    M.set(common_settings)
    M.set({
    'omega_cdm': 0.00001,
    'omega_ini_dcdm2':  0.1194,
    'Gamma_dcdm': Gamma_dcdm[i],
    'M_dcdm': 1,
    'epsilon_dcdm': epsilon[i],
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
    derived = M.get_current_derived_parameters(['sigma8','k_fss_wdm', 'age', 'rho0_wdm_over_rho0_m', 'Omega_m'])
    print("~~~~~computing our code in %.f s~~~~~"%(time.time()-t_i))
    t_i = time.time()
    
    if i==0:
        clM_0 = M.lensed_cl(2500)
#        clM_0 = M.raw_cl(2500)
        ll_DCDM_0 = clM_0['ell'][2:]
        clTT_DCDM_0 = clM_0['tt'][2:]
        clEE_DCDM_0 = clM_0['ee'][2:]
        clphiphi_DCDM_0 = clM_0['pp'][2:]
        for k in kk:
            Pk2.append(M.pk(k*h,0.)*h**3) # function .pk(k,z)
#            Pk7.append(M.pk(k*h,zhigh)*h**3) # function .pk(k,z)
        k_fss_wdm[0]=derived['k_fss_wdm']
        om_m=derived['Omega_m']
        om_l=1.0-om_m
        g0=(5.0/2.0)*om_m/((om_m**(4.0/7.0))-om_l+(1.0+0.5*om_m)*(1.0+(1.0/70.0)*om_l))
        suppression[0]=-5.0*derived['rho0_wdm_over_rho0_m']**1.75 #this formula works well only for very high values of Gamma, and not very high epsilon
#        suppression[0]=-1.0+(((1.0-derived['rho0_wdm_over_rho0_m'])**6)/((1.0-(6.0/25.0)*derived['rho0_wdm_over_rho0_m'])**2))*(g0)**(-(6./5.)*derived['rho0_wdm_over_rho0_m'])
        print("k_fss_wdm for DCDM with epsilon=%.4f and Gamma^{-1}=%.1f Gyrs is %f Mpc^-1" %(epsilon[0],tau[0], k_fss_wdm[0]) )
        print("rho_wdm/rho_m for DCDM with epsilon=%.4f  and Gamma^{-1}=%.1f Gyrs  is %f " %(epsilon[0],tau[0],derived['rho0_wdm_over_rho0_m']))
        Omega_dcdm_ini =0.1194/(M.h()**2) 
        analytical_ratio = (Omega_dcdm_ini/derived['Omega_m'])*(1.0-np.exp(-derived['age']/tau[0]))*np.sqrt(1.0-2.0*epsilon[0])
        print("and its analytical formula gives %f"%analytical_ratio)
        
    elif i==1:
        clM_1 = M.lensed_cl(2500)
        ll_DCDM_1 = clM_1['ell'][2:]
        clTT_DCDM_1 = clM_1['tt'][2:]
        clEE_DCDM_1 = clM_1['ee'][2:]
        clphiphi_DCDM_1 = clM_1['pp'][2:]
        for k in kk:
            Pk3.append(M.pk(k*h,0.)*h**3) # function .pk(k,z)
#            Pk8.append(M.pk(k*h,zhigh)*h**3) # function .pk(k,z)
        k_fss_wdm[1]=derived['k_fss_wdm']
        om_m=derived['Omega_m']
        om_l=1.0-om_m
        g0=(5.0/2.0)*om_m/((om_m**(4.0/7.0))-om_l+(1.0+0.5*om_m)*(1.0+(1.0/70.0)*om_l))
        suppression[1]=-5.0*derived['rho0_wdm_over_rho0_m']**1.75  #this formula works well only for very high values of Gamma, and not very high epsilon
#        suppression[1]=-1.0+(((1.0-derived['rho0_wdm_over_rho0_m'])**6)/((1.0-(6.0/25.0)*derived['rho0_wdm_over_rho0_m'])**2))*(g0)**(-(6./5.)*derived['rho0_wdm_over_rho0_m'])
        print("k_fss_wdm for DCDM with epsilon=%.4f  and Gamma^{-1}=%.1f Gyrs  is %f Mpc^-1" %(epsilon[1],tau[1],k_fss_wdm[1]) )
        print("rho_wdm/rho_m for DCDM with epsilon=%.4f  and Gamma^{-1}=%.1f Gyrs  is %f " %(epsilon[1],tau[1],derived['rho0_wdm_over_rho0_m']))
        Omega_dcdm_ini =0.1194/(M.h()**2) 
        analytical_ratio = (Omega_dcdm_ini/derived['Omega_m'])*(1.0-np.exp(-derived['age']/tau[1]))*np.sqrt(1.0-2.0*epsilon[1])
        print("and its analytical formula gives %f"%analytical_ratio)
    elif i==2:
        clM_2 = M.lensed_cl(2500)
        ll_DCDM_2 = clM_2['ell'][2:]
        clTT_DCDM_2 = clM_2['tt'][2:]
        clEE_DCDM_2 = clM_2['ee'][2:]
        clphiphi_DCDM_2 = clM_2['pp'][2:]
        for k in kk:
            Pk4.append(M.pk(k*h,0.)*h**3) # function .pk(k,z)
#            Pk9.append(M.pk(k*h,zhigh)*h**3) # function .pk(k,z)
        k_fss_wdm[2]=derived['k_fss_wdm']
        om_m=derived['Omega_m']
        om_l=1.0-om_m
        g0=(5.0/2.0)*om_m/((om_m**(4.0/7.0))-om_l+(1.0+0.5*om_m)*(1.0+(1.0/70.0)*om_l))
        suppression[2]=-5.0*derived['rho0_wdm_over_rho0_m']**1.75 #this formula works well only for very high values of Gamma, and not very high epsilon
#        suppression[2]=-1.0+(((1.0-derived['rho0_wdm_over_rho0_m'])**6)/((1.0-(6.0/25.0)*derived['rho0_wdm_over_rho0_m'])**2))*(g0)**(-(6./5.)*derived['rho0_wdm_over_rho0_m'])
        print("k_fss_wdm for DCDM with epsilon=%.4f  and Gamma^{-1}=%.1f Gyrs  is %f Mpc^-1" %(epsilon[2],tau[2], k_fss_wdm[2]) )
        print("rho_wdm/rho_m for DCDM with epsilon=%.4f  and Gamma^{-1}=%.1f Gyrs  is %f " %(epsilon[2],tau[2],derived['rho0_wdm_over_rho0_m']))
        Omega_dcdm_ini =0.1194/(M.h()**2) 
        analytical_ratio = (Omega_dcdm_ini/derived['Omega_m'])*(1.0-np.exp(-derived['age']/tau[2]))*np.sqrt(1.0-2.0*epsilon[2])
        print("and its analytical formula gives %f"%analytical_ratio)
    else:
        clM_3 = M.lensed_cl(2500)
        ll_DCDM_3 = clM_3['ell'][2:]
        clTT_DCDM_3 = clM_3['tt'][2:]
        clEE_DCDM_3 = clM_3['ee'][2:]
        clphiphi_DCDM_3 = clM_3['pp'][2:]
        for k in kk:
            Pk5.append(M.pk(k*h,0.)*h**3) # function .pk(k,z)
#            Pk10.append(M.pk(k*h,zhigh)*h**3) # function .pk(k,z)
        k_fss_wdm[3]=derived['k_fss_wdm']
        om_m=derived['Omega_m']
        om_l=1.0-om_m
        g0=(5.0/2.0)*om_m/((om_m**(4.0/7.0))-om_l+(1.0+0.5*om_m)*(1.0+(1.0/70.0)*om_l))
        suppression[3]=-5.0*derived['rho0_wdm_over_rho0_m']**1.75 #this formula works well only for very high values of Gamma, and not very high epsilon
#        suppression[3]=-1.0+(((1.0-derived['rho0_wdm_over_rho0_m'])**6)/((1.0-(6.0/25.0)*derived['rho0_wdm_over_rho0_m'])**2))*(g0)**(-(6./5.)*derived['rho0_wdm_over_rho0_m'])
        print("k_fss_wdm for DCDM with epsilon=%.4f  and Gamma^{-1}=%.1f Gyrs  is %f Mpc^-1" %(epsilon[3],tau[3], k_fss_wdm[3]) )
        print("rho_wdm/rho_m for DCDM with epsilon=%.4f  and Gamma^{-1}=%.1f Gyrs  is %f " %(epsilon[3],tau[3],derived['rho0_wdm_over_rho0_m']))
        Omega_dcdm_ini =0.1194/(M.h()**2) 
        analytical_ratio = (Omega_dcdm_ini/derived['Omega_m'])*(1.0-np.exp(-derived['age']/tau[3]))*np.sqrt(1.0-2.0*epsilon[3])
        print("and its analytical formula gives %f"%analytical_ratio)
        
 
    M.struct_cleanup()
    M.empty()


fTT_ourcode_0 = interp1d(ll_DCDM_0, clTT_DCDM_0)
fEE_ourcode_0 = interp1d(ll_DCDM_0, clEE_DCDM_0)
fphiphi_ourcode_0 = interp1d(ll_DCDM_0, clphiphi_DCDM_0)
fTT_ourcode_1 = interp1d(ll_DCDM_1, clTT_DCDM_1)
fEE_ourcode_1 = interp1d(ll_DCDM_1, clEE_DCDM_1)
fphiphi_ourcode_1 = interp1d(ll_DCDM_1, clphiphi_DCDM_1)
fTT_ourcode_2 = interp1d(ll_DCDM_2, clTT_DCDM_2)
fEE_ourcode_2 = interp1d(ll_DCDM_2, clEE_DCDM_2)
fphiphi_ourcode_2 = interp1d(ll_DCDM_2, clphiphi_DCDM_2)
fTT_ourcode_3 = interp1d(ll_DCDM_3, clTT_DCDM_3)
fEE_ourcode_3 = interp1d(ll_DCDM_3, clEE_DCDM_3)
fphiphi_ourcode_3 = interp1d(ll_DCDM_3, clphiphi_DCDM_3)


fpk_ourcode_0 = interp1d(kk, Pk2)
fpk_ourcode_1 = interp1d(kk, Pk3)
fpk_ourcode_2 = interp1d(kk, Pk4)
fpk_ourcode_3 = interp1d(kk, Pk5)

#fpk_ourcode_4 = interp1d(kk, Pk7)
#fpk_ourcode_5 = interp1d(kk, Pk8)
#fpk_ourcode_6 = interp1d(kk, Pk9)
#fpk_ourcode_7 = interp1d(kk, Pk10)




#%% repeat the computations but neglecting perturbations of the warm daughter 
#for i in range(3):
#    M.set(common_settings)
#    M.set({   'omega_cdm': 0.00001,   'Omega_ini_dcdm2': 0.2636, 'Gamma_dcdm': Gamma_dcdm, 'M_dcdm': 1, 'm_dcdm': m_dcdm[i], 'background_ncdm_distribution': 1, 'Quadrature strategy': 4,
#   'N_ncdm': 1, 'evolver': 0, 'Number of momentum bins perturbs': nbins, 'massive_daughter_perturbations': 'no', 'dark_radiation_perturbations': 'no','mother_dcdm_perturbations': 'yes'})

#    M.compute()
#    h = M.h()
#    if i==0:
#        clM_3 = M.lensed_cl(2500)
#        ll_DCDM_3 = clM_3['ell'][2:]
#        clTT_DCDM_3 = clM_3['tt'][2:]
#        clEE_DCDM_3 = clM_3['ee'][2:]
#        for k in kk:
#            Pk5.append(M.pk(k*h,0.)*h**3) # function .pk(k,z)

#    elif i==1:
#        clM_4 = M.lensed_cl(2500)
#        ll_DCDM_4 = clM_4['ell'][2:]
#        clTT_DCDM_4 = clM_4['tt'][2:]
#        clEE_DCDM_4 = clM_4['ee'][2:]
#        for k in kk:
#            Pk6.append(M.pk(k*h,0.)*h**3) # function .pk(k,z)
#    else:
#        clM_5 = M.lensed_cl(2500)
#        ll_DCDM_5 = clM_5['ell'][2:]
#        clTT_DCDM_5 = clM_5['tt'][2:]
#        clEE_DCDM_5 = clM_5['ee'][2:]
#        for k in kk:
#            Pk7.append(M.pk(k*h,0.)*h**3) # function .pk(k,z)

#    M.struct_cleanup()
#    M.empty()


#fTT_ourcode_3 = interp1d(ll_DCDM_3, clTT_DCDM_3)
#fEE_ourcode_3 = interp1d(ll_DCDM_3, clEE_DCDM_3)
#fTT_ourcode_4 = interp1d(ll_DCDM_4, clTT_DCDM_4)
#fEE_ourcode_4 = interp1d(ll_DCDM_4, clEE_DCDM_4)
#fTT_ourcode_5 = interp1d(ll_DCDM_5, clTT_DCDM_5)
#fEE_ourcode_5 = interp1d(ll_DCDM_5, clEE_DCDM_5)


#ax_1.semilogx(ll_DCDM_0,fTT_ourcode_3(ll_DCDM_0)/fTT_ref(ll_DCDM_0)-1,'b--',label=r'$\varepsilon = %.4f$'%epsilon[0])
#ax_2.semilogx(ll_DCDM_0,fEE_ourcode_3(ll_DCDM_0)/fEE_ref(ll_DCDM_0)-1,'b--',label=r'$\varepsilon = %.4f$'%epsilon[0])

#ax_1.semilogx(ll_DCDM_0,fTT_ourcode_4(ll_DCDM_0)/fTT_ref(ll_DCDM_0)-1,'r--',label=r'$\varepsilon = %.4f$'%epsilon[1])
#ax_2.semilogx(ll_DCDM_0,fEE_ourcode_4(ll_DCDM_0)/fEE_ref(ll_DCDM_0)-1,'r--',label=r'$\varepsilon = %.4f$'%epsilon[1])

#ax_1.semilogx(ll_DCDM_0,fTT_ourcode_5(ll_DCDM_0)/fTT_ref(ll_DCDM_0)-1,'g--',label=r'$\varepsilon = %.4f$'%epsilon[2])
#ax_2.semilogx(ll_DCDM_0,fEE_ourcode_5(ll_DCDM_0)/fEE_ref(ll_DCDM_0)-1,'g--',label=r'$\varepsilon = %.4f$'%epsilon[2])



#plt.figure(1)
#plt.xscale('log')
#plt.xlim(kk[0],kk[-1])
#plt.plot(kk,map(truediv, list(np.array(Pk5) - np.array(Pk1)), Pk1),'b--',label=r'$\varepsilon = %.4f$'%epsilon[0])
#plt.plot(kk,map(truediv, list(np.array(Pk6) - np.array(Pk1)), Pk1),'r--', label=r'$\varepsilon = %.4f$'%epsilon[1])
#plt.plot(kk,map(truediv, list(np.array(Pk7) - np.array(Pk1)), Pk1),'g--', label = r'$\varepsilon = %.4f$'%epsilon[2])

#%%

print("~~~~~ready to plot~~~~~")
#plot

##create plot

ax_1 = plt.subplot(311)
ax_2 = plt.subplot(312, sharex = ax_1)
ax_3 = plt.subplot(313, sharex = ax_2)
plt.subplots_adjust(hspace=0)

ax_1.set_ylim([-0.19,0.19])
ax_2.set_ylim([-0.11,0.11])
ax_3.set_ylim([-0.41,0.41])



ax_1.set_xlim([2,2500])
ax_2.set_xlim([2,2500])
ax_3.set_xlim([2,2500])


ax_3.tick_params(axis='both', which='minor', labelsize=12)

#Planck error bars
# FOR TT
#l_cosmic_variance_1 = np.linspace(0,30,1000)
#l_cosmic_variance_2 = np.linspace(30,48,2)
#slope =np.array([0.15,0.0343])
#ax_1.fill_between(l_cosmic_variance_1, -0.18,0.18, color='lightgray' )
#ax_1.fill_between(l_cosmic_variance_2, -slope, slope, color='lightgray' )
#ax_1.fill_between(lTT, -(DlTT_error_plus)/DlTT_mean, +(DlTT_error_plus)/DlTT_mean, color='lightgray')



if plot_type ==1:
    
    ax_1.semilogx(ll_DCDM_0,fTT_ourcode_0(ll_DCDM_0)/fTT_ref(ll_DCDM_0)-1,'r',label=r'$\varepsilon = %.1f$'%epsilon[0])
    ax_2.semilogx(ll_DCDM_0,fEE_ourcode_0(ll_DCDM_0)/fEE_ref(ll_DCDM_0)-1,'r',label=r'$\varepsilon = %.1f$'%epsilon[0])
    ax_3.semilogx(ll_DCDM_0,fphiphi_ourcode_0(ll_DCDM_0)/fphiphi_ref(ll_DCDM_0)-1,'r',label=r'$\varepsilon = %.1f$'%epsilon[0])
    
    ax_1.semilogx(ll_DCDM_0,fTT_ourcode_1(ll_DCDM_0)/fTT_ref(ll_DCDM_0)-1,'g',label=r'$\varepsilon = %.1f$'%epsilon[1])
    ax_2.semilogx(ll_DCDM_0,fEE_ourcode_1(ll_DCDM_0)/fEE_ref(ll_DCDM_0)-1,'g',label=r'$\varepsilon = %.1f$'%epsilon[1])
    ax_3.semilogx(ll_DCDM_0,fphiphi_ourcode_1(ll_DCDM_0)/fphiphi_ref(ll_DCDM_0)-1,'g',label=r'$\varepsilon = %.1f$'%epsilon[1])

    ax_1.semilogx(ll_DCDM_0,fTT_ourcode_2(ll_DCDM_0)/fTT_ref(ll_DCDM_0)-1,'b',label=r'$\varepsilon = %.2f$'%epsilon[2])
    ax_2.semilogx(ll_DCDM_0,fEE_ourcode_2(ll_DCDM_0)/fEE_ref(ll_DCDM_0)-1,'b',label=r'$\varepsilon = %.2f$'%epsilon[2])
    ax_3.semilogx(ll_DCDM_0,fphiphi_ourcode_2(ll_DCDM_0)/fphiphi_ref(ll_DCDM_0)-1,'b',label=r'$\varepsilon = %.2f$'%epsilon[2])
    
    ax_1.semilogx(ll_DCDM_0,fTT_ourcode_3(ll_DCDM_0)/fTT_ref(ll_DCDM_0)-1,'k',label=r'$\varepsilon = %.3f$'%epsilon[3])
    ax_2.semilogx(ll_DCDM_0,fEE_ourcode_3(ll_DCDM_0)/fEE_ref(ll_DCDM_0)-1,'k',label=r'$\varepsilon = %.3f$'%epsilon[3])
    ax_3.semilogx(ll_DCDM_0,fphiphi_ourcode_3(ll_DCDM_0)/fphiphi_ref(ll_DCDM_0)-1,'k',label=r'$\varepsilon = %.3f$'%epsilon[3])

else:
    
    ax_1.semilogx(ll_DCDM_0,fTT_ourcode_0(ll_DCDM_0)/fTT_ref(ll_DCDM_0)-1,'k',label=r'$\Gamma^{-1} = %.0f \, \mathrm{Gyrs}$'%tau[0])
    ax_2.semilogx(ll_DCDM_0,fEE_ourcode_0(ll_DCDM_0)/fEE_ref(ll_DCDM_0)-1,'k',label=r'$\Gamma^{-1} = %.0f \, \mathrm{Gyrs}$'%tau[0])
    ax_3.semilogx(ll_DCDM_0,fphiphi_ourcode_0(ll_DCDM_0)/fphiphi_ref(ll_DCDM_0)-1,'k',label=r'$\Gamma^{-1} = %.0f \, \mathrm{Gyrs}$'%tau[0])

    ax_1.semilogx(ll_DCDM_0,fTT_ourcode_1(ll_DCDM_0)/fTT_ref(ll_DCDM_0)-1,'b',label=r'$\Gamma^{-1} = %.0f \, \mathrm{Gyrs}$'%tau[1])
    ax_2.semilogx(ll_DCDM_0,fEE_ourcode_1(ll_DCDM_0)/fEE_ref(ll_DCDM_0)-1,'b',label=r'$\Gamma^{-1} = %.0f \, \mathrm{Gyrs}$'%tau[1])
    ax_3.semilogx(ll_DCDM_0,fphiphi_ourcode_1(ll_DCDM_0)/fphiphi_ref(ll_DCDM_0)-1,'b',label=r'$\Gamma^{-1} = %.0f \, \mathrm{Gyrs}$'%tau[1])

    ax_1.semilogx(ll_DCDM_0,fTT_ourcode_2(ll_DCDM_0)/fTT_ref(ll_DCDM_0)-1,'g',label=r'$\Gamma^{-1} = %.0f \, \mathrm{Gyrs}$'%tau[2])
    ax_2.semilogx(ll_DCDM_0,fEE_ourcode_2(ll_DCDM_0)/fEE_ref(ll_DCDM_0)-1,'g',label=r'$$\Gamma^{-1} = %.0f \, \mathrm{Gyrs}$'%tau[2])
    ax_3.semilogx(ll_DCDM_0,fphiphi_ourcode_2(ll_DCDM_0)/fphiphi_ref(ll_DCDM_0)-1,'g',label=r'$$\Gamma^{-1} = %.0f \, \mathrm{Gyrs}$'%tau[2])

    ax_1.semilogx(ll_DCDM_0,fTT_ourcode_3(ll_DCDM_0)/fTT_ref(ll_DCDM_0)-1,'r',label=r'$\Gamma^{-1} = %.0f \, \mathrm{Gyrs}$'%tau[3])
    ax_2.semilogx(ll_DCDM_0,fEE_ourcode_3(ll_DCDM_0)/fEE_ref(ll_DCDM_0)-1,'r',label=r'$\Gamma^{-1} = %.0f \, \mathrm{Gyrs}$'%tau[3])
    ax_3.semilogx(ll_DCDM_0,fphiphi_ourcode_3(ll_DCDM_0)/fphiphi_ref(ll_DCDM_0)-1,'r',label=r'$\Gamma^{-1} = %.0f \, \mathrm{Gyrs}$'%tau[3])




ax_3.set_xlabel(r'$\mathrm{multipole} \, \ell$',fontsize=18)
ax_1.set_ylabel(r'$\frac{C_\ell^\mathrm{TT}(\Lambda \mathrm{DDM})}{C_\ell^\mathrm{TT}(\Lambda \mathrm{CDM} )} -1$',fontsize=23)
ax_2.set_ylabel(r'$\frac{C_\ell^\mathrm{EE}(\Lambda \mathrm{DDM})}{C_\ell^\mathrm{EE}(\Lambda \mathrm{CDM} )} -1$',fontsize=23)
ax_3.set_ylabel(r'$\frac{C_\ell^{\phi \phi}(\Lambda \mathrm{DDM})}{C_\ell^{\phi \phi}(\Lambda \mathrm{CDM})} -1$',fontsize=23)


ax_3.tick_params(axis="x", labelsize=16)
ax_1.tick_params(axis="y", labelsize=16)
ax_2.tick_params(axis="y", labelsize=16)
ax_3.tick_params(axis="y", labelsize=16)

#ax_1.locator_params(axis='y', nbins=5)
#ax_2.locator_params(axis='y', nbins=5)
#ax_3.locator_params(axis='y', nbins=5)


ax_1.yaxis.set_major_locator(plt.MaxNLocator(6))
ax_2.yaxis.set_major_locator(plt.MaxNLocator(6))
ax_3.yaxis.set_major_locator(plt.MaxNLocator(6))


if plot_type == 1:
    ax_1.text(400, -0.1, r'$\Gamma^{-1} = %.0f \,  \mathrm{Gyrs}$'%tau[0], fontsize =16)
    
else:
    ax_1.text(400, -0.1, r'$\varepsilon = %.1f$'%epsilon[0], fontsize =16)
    
ax_1.legend(frameon=False,fontsize =16,loc='lower left',borderaxespad=0., ncol=2)




if plot_type ==1:
    plt.savefig('cl_dcdm_full_several_epsilon.png',dpi=300,bbox_inches='tight')
else:
    plt.savefig('cl_dcdm_full_several_gamma.png',dpi=300,bbox_inches='tight')


plt.show()
plt.clf()

#%%
#plt.figure(2)

#axe_1 = plt.subplot(211)
#axe_2 = plt.subplot(212, sharex = ax_1)
#plt.subplots_adjust(hspace=0)

axe_1 = plt.subplot(111)

#axe_1.set_ylim([-1.5,2])
#axe_2.set_ylim([-0.14,0.23])

axe_1.set_xlim([kk[0],kk[-1]])

#axe_2.set_xlim([kk[0],kk[-1]])



#plt.xlim(kk[0],kk[-1])
#plt.ylim(-1.5,2.0)
#plt.ylabel(r'$P_{\Lambda\mathrm{DDM}}/P_{\Lambda\mathrm{CDM}}-1$', fontsize=17)

#plt.text(0.000025, 0.22, r'$P_{\Lambda\mathrm{DDM}}/P_{\Lambda\mathrm{CDM}}-1$', va='center', rotation='vertical',fontsize=23)
#axe_2.set_xlabel(r'$k \,\,\,\, [h/\mathrm{Mpc}]$', fontsize=18)

axe_1.set_ylabel(r'$P(k) \,\,\,\, [(\mathrm{Mpc}/h)^3] $', fontsize=18)
axe_1.set_xlabel(r'$k \,\,\,\, [h/\mathrm{Mpc}]$', fontsize=18)


if plot_type ==1:
    
    axe_1.loglog(kk,fpk_ref(kk),color='black',label=r'$\Lambda \mathrm{CDM}$')
    axe_1.loglog(kk,fpk_ourcode_0(kk),color='orange',label=r'$\varepsilon = %.1f$'%epsilon[0])
    axe_1.loglog(kk,fpk_ourcode_1(kk),color='cyan', label=r'$\varepsilon = %.1f$'%epsilon[1])
    axe_1.loglog(kk,fpk_ourcode_2(kk),color='dodgerblue', label = r'$\varepsilon = %.2f$'%epsilon[2])
    axe_1.loglog(kk,fpk_ourcode_3(kk),color='steelblue', label = r'$\varepsilon = %.3f$'%epsilon[3])

    
#    axe_1.semilogx(kk,fpk_ourcode_0(kk)/fpk_ref(kk)-1.0,'r',label=r'$\varepsilon = %.1f$'%epsilon[0])
#    axe_1.semilogx(kk,fpk_ourcode_1(kk)/fpk_ref(kk)-1.0,'g', label=r'$\varepsilon = %.1f$'%epsilon[1])
#    axe_1.semilogx(kk,fpk_ourcode_2(kk)/fpk_ref(kk)-1.0,'b', label = r'$\varepsilon = %.2f$'%epsilon[2])
#    axe_1.semilogx(kk,fpk_ourcode_3(kk)/fpk_ref(kk)-1.0,'k', label = r'$\varepsilon = %.3f$'%epsilon[3])

#    axe_2.semilogx(kk,fpk_ourcode_4(kk)/fpk_ref2(kk)-1.0,'r',label=r'$\varepsilon = %.1f$'%epsilon[0])
#    axe_2.semilogx(kk,fpk_ourcode_5(kk)/fpk_ref2(kk)-1.0,'g', label=r'$\varepsilon = %.1f$'%epsilon[1])
#    axe_2.semilogx(kk,fpk_ourcode_6(kk)/fpk_ref2(kk)-1.0,'b', label = r'$\varepsilon = %.2f$'%epsilon[2])
#    axe_2.semilogx(kk,fpk_ourcode_7(kk)/fpk_ref2(kk)-1.0,'k', label = r'$\varepsilon = %.3f$'%epsilon[3])
    
else:
    
    axe_1.semilogx(kk,fpk_ref(kk),color='black',label=r'$\Lambda \mathrm{CDM}$')
    axe_1.semilogx(kk,fpk_ourcode_0(kk),color='orange',label=r'$\Gamma^{-1} = %.0f \, \mathrm{Gyrs}$'%tau[0])
    axe_1.semilogx(kk,fpk_ourcode_1(kk),color='cyan', label=r'$\Gamma^{-1} = %.0f \, \mathrm{Gyrs}$'%tau[1])
    axe_1.semilogx(kk,fpk_ourcode_2(kk),color='dodgerblue', label = r'$\Gamma^{-1} = %.0f \, \mathrm{Gyrs}$'%tau[2])
    axe_1.semilogx(kk,fpk_ourcode_3(kk),color='steelblue', label = r'$\Gamma^{-1} = %.0f \, \mathrm{Gyrs}$'%tau[3])

    
#    axe_1.semilogx(kk,fpk_ourcode_0(kk)/fpk_ref(kk)-1.0,'r',label=r'$\Gamma^{-1} = %.0f \, \mathrm{Gyrs}$'%tau[0])
#    axe_1.semilogx(kk,fpk_ourcode_1(kk)/fpk_ref(kk)-1.0,'g', label=r'$\Gamma^{-1} = %.0f \, \mathrm{Gyrs}$'%tau[1])
#   axe_1.semilogx(kk,fpk_ourcode_2(kk)/fpk_ref(kk)-1.0,'b', label = r'$\Gamma^{-1} = %.0f \, \mathrm{Gyrs}$'%tau[2])
#   axe_1.semilogx(kk,fpk_ourcode_3(kk)/fpk_ref(kk)-1.0,'k', label = r'$\Gamma^{-1} = %.0f \, \mathrm{Gyrs}$'%tau[3])

#    axe_2.semilogx(kk,fpk_ourcode_4(kk)/fpk_ref2(kk)-1.0,'r',label=r'$\Gamma^{-1} = %.0f \, \mathrm{Gyrs}$'%tau[0])
#    axe_2.semilogx(kk,fpk_ourcode_5(kk)/fpk_ref2(kk)-1.0,'g', label=r'$\Gamma^{-1} = %.0f \, \mathrm{Gyrs}$'%tau[1])
#    axe_2.semilogx(kk,fpk_ourcode_6(kk)/fpk_ref2(kk)-1.0,'b', label = r'$\Gamma^{-1} = %.0f \, \mathrm{Gyrs}$'%tau[2])
#    axe_2.semilogx(kk,fpk_ourcode_7(kk)/fpk_ref2(kk)-1.0,'k', label = r'$\Gamma^{-1} = %.0f \, \mathrm{Gyrs}$'%tau[3])



axe_1.legend(frameon=False, loc='lower left', fontsize=16, borderaxespad=0.)

#axe_1.axvline(k_fss_wdm[0], color='r', linestyle='dotted')
#axe_1.axvline(k_fss_wdm[1], color='g', linestyle='dotted')
#axe_1.axvline(k_fss_wdm[2], color='b', linestyle='dotted')
#axe_1.axvline(k_fss_wdm[3], color='k', linestyle='dotted')

#axe_1.tick_params(axis='x', which='both', bottom='False',labelbottom='False')
#axe_1.tick_params(axis="y", labelsize=16)
#axe_2.tick_params(axis="y", labelsize=16)
#axe_2.tick_params(axis="x", labelsize=16)

axe_1.tick_params(axis="y", labelsize=16)
axe_1.tick_params(axis="x", labelsize=16)



#axe_1.axhline(suppression[0], color='r', linestyle='dashed')
#axe_1.axhline(suppression[1], color='g', linestyle='dashed')
#axe_1.axhline(suppression[2], color='b', linestyle='dashed')
#axe_1.axhline(suppression[3], color='k', linestyle='dashed')

if plot_type ==1:
    axe_1.text(0.13e-3, -0.2, r'$\Gamma^{-1} = %.0f \, \mathrm{Gyrs} $'%tau[0], fontsize =16)
#    axe_1.text(0.6e-2, -0.5,  r'$\bf{z=0}$', fontsize =16)
#    axe_2.text(0.6e-2, -0.15,  r'$\bf{z=%.0f}$'%zhigh, fontsize =16)

    
else:
    axe_1.text(0.13e-3, -0.2,  r'$\varepsilon = %.1f$'%epsilon[0], fontsize =16)
#    axe_1.text(0.6e-2, -0.8,  r'$\bf{z=0}$', fontsize =16)
#    axe_2.text(0.6e-2, -0.4,  r'$\bf{z=%.0f}$'%zhigh, fontsize =16)





#plt.tight_layout()


if plot_type ==1:
    plt.savefig('pk_dcdm_full_several_epsilon.png',dpi=300,bbox_inches='tight')
else:
    plt.savefig('pk_dcdm_full_several_gamma.png',dpi=300,bbox_inches='tight')

plt.tight_layout()
plt.show()
plt.clf()

#%%
