
# import necessary modules
import matplotlib.pyplot as plt
import numpy as np
from classy import Class
from operator import truediv
from scipy.interpolate import interp1d

import time
start_time = time.time()

#%%


Gamma_dcdm = 19.6078
tau =1./(Gamma_dcdm*1.02e-3)
tau

kk = np.logspace(-4,np.log10(1),575) # k in h/Mpc

len(kk)
Pk_0p15_sigYES = [] # P(k) in (Mpc/h)**3
Pk_0p2_sigYES  = [] # P(k) in (Mpc/h)**3
Pk_0p25_sigYES = [] # P(k) in (Mpc/h)**3
Pk_0p3_sigYES  = [] # P(k) in (Mpc/h)**3
Pk_0p45_sigYES = [] # P(k) in (Mpc/h)**3
Pk_0p499_sigYES = [] # P(k) in (Mpc/h)**3
Pk_0p49_sigYES = [] # P(k) in (Mpc/h)**3


Pk_0p15_sigNOT = [] # P(k) in (Mpc/h)**3
Pk_0p2_sigNOT  = [] # P(k) in (Mpc/h)**3
Pk_0p25_sigNOT = [] # P(k) in (Mpc/h)**3
Pk_0p3_sigNOT  = [] # P(k) in (Mpc/h)**3
Pk_0p45_sigNOT = [] # P(k) in (Mpc/h)**3
Pk_0p499_sigNOT = [] # P(k) in (Mpc/h)**3
Pk_0p49_sigNOT = [] # P(k) in (Mpc/h)**3


#set general configuration
common_settings = {'output':'tCl,pCl,lCl,mPk',
                   'lensing':'yes',
                   'l_max_scalars':2600,
                   'P_k_max_1/Mpc':1.0,
                   'omega_ini_dcdm2':0.1195,
                   'omega_cdm':0.00001,
                   'Gamma_dcdm':Gamma_dcdm,
                   'n_s':0.9671,
                   'ln10^{10}A_s':3.048,
                   'tau_reio':0.0562,
                   'omega_b':0.02238,
                   'H0':67.65,
                   'N_ncdm':2,
                   'M_dcdm':1,
                   'background_ncdm_distribution': '0, 1',
                   'Quadrature strategy': '0,4',
                   'm_ncdm': '0.06, 0',
                   'N_ur':2.0328,
                   'evolver':0,
                   'l_max_ncdm':17,
                   'Number of momentum bins perturbs': '50, 300',
                   'ncdm_fluid_approximation':2,
                   'ncdm_fluid_trigger_tau_over_tau_k': 25,
                   'massive_daughter_perturbations' : 'yes',
                   'dark_radiation_perturbations': 'yes',
#                   'back_integration_stepsize': 7.0e-4
}

#%% read data files

files1 = ['/home/guillermo/class_MAJORON/output/dcdm_eps0p15_tau50gyrs_nofluid_pk.dat']
data1 = []
for data_file1 in files1:
    data1.append(np.loadtxt(data_file1))
    
files2 = ['/home/guillermo/class_MAJORON/output/dcdm_eps0p2_tau50gyrs_nofluid_pk.dat']
data2 = []


for data_file2 in files2:
    data2.append(np.loadtxt(data_file2))

files3 = ['/home/guillermo/class_MAJORON/output/dcdm_eps0p25_tau50gyrs_nofluid_pk.dat']
data3 = []
for data_file3 in files3:
    data3.append(np.loadtxt(data_file3))


files4 = ['/home/guillermo/class_MAJORON/output/dcdm_eps0p3_tau50gyrs_nofluid_pk.dat']
data4 = []
for data_file4 in files4:
    data4.append(np.loadtxt(data_file4))
    
files5 = ['/home/guillermo/class_MAJORON/output/dcdm_eps0p45_tau50gyrs_nofluid_pk.dat']
data5 = []
for data_file5 in files5:
    data5.append(np.loadtxt(data_file5))

files6 = ['/home/guillermo/class_MAJORON/output/dcdm_eps0p499_tau50gyrs_nofluid_pk.dat']
data6 = []
for data_file6 in files6:
    data6.append(np.loadtxt(data_file6))
    

files7 = ['/home/guillermo/class_MAJORON/output/dcdm_eps0p49_tau50gyrs_nofluid_pk.dat']
data7 = []
for data_file7 in files7:
    data7.append(np.loadtxt(data_file7))
    
          
    
pk_ref_0p15  = data1[0]
pk_ref_0p2   = data2[0]
pk_ref_0p25  = data3[0]
pk_ref_0p3   = data4[0]
pk_ref_0p45  = data5[0]
pk_ref_0p499 = data6[0]
pk_ref_0p49  = data7[0]

M = Class()
#%% EPSILON=0.15 CASE ########################################################

#M.set(common_settings)
#M.set({'epsilon_dcdm': 0.15, 'switch_off_shear_wdm': 'yes'}) 

#M.compute()

#h = M.h()
#for k in kk:
#    Pk_0p15_sigNOT.append(M.pk(k*h,0.)*h**3) # function .pk(k,z)

#M.struct_cleanup()
#M.empty()

############################################################################
#M.set(common_settings)
#M.set({'epsilon_dcdm': 0.15, 'switch_off_shear_wdm': 'no'}) 

#M.compute()

#h = M.h()
#for k in kk:
#    Pk_0p15_sigYES.append(M.pk(k*h,0.)*h**3) # function .pk(k,z)

#M.struct_cleanup()
#M.empty()

#print("~~~~~computing fluid wdm (both full shear and zero shear) for epsilon=0.15 in in %.f s~~~~~"%(time.time()-start_time))
#t_1=time.time()

#%% EPSILON=0.2 CASE ########################################################

#M.set(common_settings)
#M.set({'epsilon_dcdm': 0.2,'switch_off_shear_wdm': 'yes'}) 

#M.compute()


#print("~~~~~computing fluid wdm for epsilon=0.2 in in %.f s~~~~~"%(time.time()-t_1))
#t_2=time.time()

#h = M.h()
#for k in kk:
#    Pk_0p2_sigNOT.append(M.pk(k*h,0.)*h**3) # function .pk(k,z)

#M.struct_cleanup()
#M.empty()


###############################################################################

#M.set(common_settings)
#M.set({'epsilon_dcdm': 0.2,'switch_off_shear_wdm': 'no'}) 

#M.compute()


#h = M.h()
#for k in kk:
#    Pk_0p2_sigYES.append(M.pk(k*h,0.)*h**3) # function .pk(k,z)

#M.struct_cleanup()
#M.empty()

#%% EPSILON=0.25 CASE ########################################################

#M.set(common_settings)
#M.set({'epsilon_dcdm': 0.25, 'switch_off_shear_wdm': 'yes'}) 

#M.compute()


#print("~~~~~computing fluid wdm  for epsilon=0.25 in in %.f s~~~~~"%(time.time()-t_2))
#t_3=time.time()

#h = M.h()
#for k in kk:
#    Pk_0p25_sigNOT.append(M.pk(k*h,0.)*h**3) # function .pk(k,z)

#M.struct_cleanup()
#M.empty()


###########################################################################

#M.set(common_settings)
#M.set({'epsilon_dcdm': 0.25, 'switch_off_shear_wdm': 'no'}) 

#M.compute()

#h = M.h()
#for k in kk:
#    Pk_0p25_sigYES.append(M.pk(k*h,0.)*h**3) # function .pk(k,z)

#M.struct_cleanup()
#M.empty()

#%% EPSILON=0.3 CASE ########################################################

#M.set(common_settings)
#M.set({'epsilon_dcdm': 0.3, 'switch_off_shear_wdm': 'yes'}) 

#M.compute()


#print("~~~~~computing fluid wdm for epsilon=0.3 in in %.f s~~~~~"%(time.time()-t_3))

#h = M.h()
#for k in kk:
#    Pk_0p3_sigNOT.append(M.pk(k*h,0.)*h**3) # function .pk(k,z)

#M.struct_cleanup()
#M.empty()

#######################################################################

#M.set(common_settings)
#M.set({'epsilon_dcdm': 0.3, 'switch_off_shear_wdm': 'no'}) 

#M.compute()

#h = M.h()
#for k in kk:
#    Pk_0p3_sigYES.append(M.pk(k*h,0.)*h**3) # function .pk(k,z)

#M.struct_cleanup()
#M.empty()

#%% EPSILON =0.45 CASE ########################################################


#M.set(common_settings)
#M.set({'epsilon_dcdm': 0.45, 'switch_off_shear_wdm': 'yes'}) 

#M.compute()


#print("~~~~~computing fluid wdm  for epsilon=0.25 in in %.f s~~~~~"%(time.time()-t_2))
#t_3=time.time()

#h = M.h()
#for k in kk:
#    Pk_0p45_sigNOT.append(M.pk(k*h,0.)*h**3) # function .pk(k,z)

#M.struct_cleanup()
#M.empty()


###########################################################################

#M.set(common_settings)
#M.set({'epsilon_dcdm': 0.45, 'switch_off_shear_wdm': 'no'}) 

#M.compute()

#h = M.h()
#for k in kk:
#    Pk_0p45_sigYES.append(M.pk(k*h,0.)*h**3) # function .pk(k,z)

#M.struct_cleanup()
#M.empty()


#%% EPSILON =0.499 CASE ######################################################
#M.set(common_settings)
#M.set({'epsilon_dcdm': 0.499, 'switch_off_shear_wdm': 'yes'}) 

#M.compute()


#print("~~~~~computing fluid wdm  for epsilon=0.25 in in %.f s~~~~~"%(time.time()-t_2))
#t_3=time.time()

#h = M.h()
#for k in kk:
#    Pk_0p499_sigNOT.append(M.pk(k*h,0.)*h**3) # function .pk(k,z)

#M.struct_cleanup()
#M.empty()
###########################################################################

#M.set(common_settings)
#M.set({'epsilon_dcdm': 0.499, 'switch_off_shear_wdm': 'no'}) 

#M.compute()

#h = M.h()
#for k in kk:
#    Pk_0p499_sigYES.append(M.pk(k*h,0.)*h**3) # function .pk(k,z)

#M.struct_cleanup()
#M.empty()

#%% EPSILON =0.49 CASE ######################################################
#M.set(common_settings)
#M.set({'epsilon_dcdm': 0.49, 'switch_off_shear_wdm': 'yes'}) 

#M.compute()


#print("~~~~~computing fluid wdm  for epsilon=0.25 in in %.f s~~~~~"%(time.time()-t_2))
#t_3=time.time()

#h = M.h()
#for k in kk:
#    Pk_0p49_sigNOT.append(M.pk(k*h,0.)*h**3) # function .pk(k,z)

#M.struct_cleanup()
#M.empty()
###########################################################################

#M.set(common_settings)
#M.set({'epsilon_dcdm': 0.49, 'switch_off_shear_wdm': 'no'}) 

#M.compute()

#h = M.h()
#for k in kk:
#    Pk_0p49_sigYES.append(M.pk(k*h,0.)*h**3) # function .pk(k,z)

#M.struct_cleanup()
#M.empty()


#%% START TO PLOT
print("~~~~~ready to plot~~~~~")
plt.figure(1)
plt.xscale('log')
plt.yscale('log')
plt.xlim(kk[0],kk[-1])


plt.xlabel(r'$k \,\,\,\, [h/\mathrm{Mpc}]$', fontsize=15)
#plt.ylabel(r'$P_{\mathrm{approx}}/P_{\mathrm{full}}-1$', fontsize=20)

plt.ylabel(r'$P(k) \, [(\mathrm{Mpc}/h)^3]$', fontsize=20)

plt.title(r"$\varepsilon = 0.49, \, \, \, \, \Gamma^{-1} = 50 \, \mathrm{Gyrs}$",fontsize=15)

#plt.plot(kk,map(truediv, list(np.array(Pk_0p15_sigNOT) - np.array(pk_ref_0p15[:,1])), pk_ref_0p15[:,1]),'b',label=r'$ \sigma_{\mathrm{wdm}} = 0$')

plt.plot(pk_ref_0p49[:,0],pk_ref_0p49[:,1],'b',label=r'Full hierarchy')
plt.plot(kk,Pk_0p49_sigYES,'r',label=r'Fluid wdm with full shear')
plt.plot(kk,Pk_0p49_sigNOT,'g',label=r'Fluid wdm with shear=0')



plt.legend(loc='best', fontsize=13)

plt.show()
plt.clf()

