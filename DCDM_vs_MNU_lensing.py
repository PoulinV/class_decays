
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

#%% read data files

files1 = ['/Users/gfranco/cloud/output/dcdm_bestfit_nofluid_cl_lensed.dat']
data1 = []
for data_file1 in files1:
    data1.append(np.loadtxt(data_file1))
    
files2 = ['/Users/gfranco/cloud/output/dcdm_bestfit_nofluid_z1_pk.dat']
data2 = []
for data_file2 in files2:
    data2.append(np.loadtxt(data_file2))

files3 = ['/Users/gfranco/cloud/output/dcdm_bestfit_nofluid_z2_pk.dat']
data3 = []
for data_file3 in files3:
    data3.append(np.loadtxt(data_file3))   
    
cl_lens_dcdm = data1[0]
pk_z0_dcdm   = data2[0]
pk_z3_dcdm   = data3[0]

fcl_tt_dcdm_full  = interp1d(cl_lens_dcdm[:,0],cl_lens_dcdm[:,1])
fcl_ee_dcdm_full  = interp1d(cl_lens_dcdm[:,0],cl_lens_dcdm[:,2])
fcl_pp_dcdm_full  = interp1d(cl_lens_dcdm[:,0],cl_lens_dcdm[:,5])

fpk_z0_dcdm_full  = interp1d(pk_z0_dcdm[:,0],pk_z0_dcdm[:,1])
fpk_z3_dcdm_full  = interp1d(pk_z3_dcdm[:,0],pk_z3_dcdm[:,1])


#%%

f = plt.figure(figsize=(10,15))

ax_1 = f.add_subplot(311)
ax_2 = f.add_subplot(312, sharex = ax_1)
ax_3 = f.add_subplot(313, sharex = ax_2)
plt.subplots_adjust(hspace=0)

ax_1.set_ylim([-0.04,0.04])
ax_2.set_ylim([-0.21,0.21])
ax_3.set_ylim([-0.19,0.19])
ax_1.set_xlim([2,2500])
ax_2.set_xlim([2,2500])
ax_3.set_xlim([2,2500])


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
                   'n_s':0.9663,
                   'ln10^{10}A_s':3.045,
                   'tau_reio':0.055,
#                   'z_reio':7.732,
                   'omega_b':0.02242,
                   'P_k_max_h/Mpc':1.0,
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
'omega_cdm': 0.1194, #Omega_cdm 0.2626
'100*theta_s':1.042059, #obtained by running class with all these parameters and H0=67.76
#'H0':67.76, 
'N_ncdm':1,
'background_ncdm_distribution': 0,
'N_ur':0.00641,
'deg_ncdm':3,
'm_ncdm':0.02,
#'input_verbose': 1,
#'background_verbose' :1,
#'thermodynamics_verbose' :1,
#'perturbations_verbose': 1,
#'transfer_verbose': 1,
#'primordial_verbose': 1,
#'spectra_verbose': 1,
#'nonlinear_verbose': 1,
#'lensing_verbose' :1,
#'output_verbose': 1
})
M.compute()

derived = M.get_current_derived_parameters(['sigma8','Omega_m','z_eq'])
S8 = derived['sigma8']*np.sqrt(derived['Omega_m']/0.3)

#print("sigma8 for LCDM is %f"%derived['sigma8'])
#print("and Omega_m is %f, so that S8=%f "%(derived['Omega_m'],S8))
#print("and z_eq is %f"%derived['z_eq'])

clM = M.lensed_cl(2600)
ll_LCDM = clM['ell'][2:]
clTT_LCDM = clM['tt'][2:]
clEE_LCDM = clM['ee'][2:]
clphiphi_LCDM = clM['pp'][2:]

# get P(k) at redhsift z=0
h = M.h() # get reduced Hubble for conversions to 1/Mpc

for k in kk:
    Pk1.append(M.pk(k*h,0.)*h**3) # function .pk(k,z)
    
for k in kk:
    Pk4.append(M.pk(k*h,3.)*h**3) # function .pk(k,z)


M.struct_cleanup()
M.empty()

#AND OTHER IN WHICH WE FIX H0 ############################################

#M.set(common_settings)
#we consider three massive degenerate neutrinos 
#M.set({
#'omega_cdm': 0.1196, #Omega_cdm 0.2626
#'H0':67.58, #value corresponding to the reference values above (from 100*theta_s using shooting method)
#'N_ncdm':1,
#'background_ncdm_distribution': 0,
#'N_ur':0.00641,
#'deg_ncdm':3,
#'m_ncdm':0.02})
#M.compute()

# get P(k) at redhsift z=0
#h = M.h() # get reduced Hubble for conversions to 1/Mpc

#for k in kk:
#    Pk1.append(M.pk(k*h,0.)*h**3) # function .pk(k,z)
    
#for k in kk:
#    Pk4.append(M.pk(k*h,3.)*h**3) # function .pk(k,z)

#M.struct_cleanup()
#M.empty()

fTT_ref = interp1d(ll_LCDM,ll_LCDM*(ll_LCDM+1.0)*clTT_LCDM/(2.0*np.pi))
fEE_ref = interp1d(ll_LCDM,ll_LCDM*(ll_LCDM+1.0)*clEE_LCDM/(2.0*np.pi))
fphiphi_ref = interp1d(ll_LCDM,ll_LCDM*(ll_LCDM+1.0)*clphiphi_LCDM/(2.0*np.pi))


fpk_z0_ref  = interp1d(kk,Pk1)
fpk_z3_ref  = interp1d(kk,Pk4)


timeafterref=time.time()


#%% compute our best-fit DCDM

t_i=timeafterref
print("~~~~~time =%.f s; computing our code~~~~~"%(timeafterref-start_time))
M.set({'output':'tCl,pCl,lCl,mPk',
                   'lensing':'yes',
                   'l_max_scalars':2600,
                   'n_s':0.9673,
                   'ln10^{10}A_s':3.052,
                   'tau_reio':0.0582,
#                   'z_reio':8.058,
                   'omega_b':0.02240,
#                   'H0':67.70,
                   '100*theta_s':1.042174, #obtained by running class with all these parameters and H0=67.70
                   'P_k_max_h/Mpc':1.0,
                   'z_max_pk' : 4.0
                   })
    
M.set({
    'omega_cdm': 0.000001,
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
    'dark_radiation_perturbations': 'yes',
#    'input_verbose': 1,
#    'background_verbose' :1,
#    'thermodynamics_verbose' :1,
#    'perturbations_verbose': 1,
#    'transfer_verbose': 1,
#    'primordial_verbose': 1,
#    'spectra_verbose': 1,
#    'nonlinear_verbose': 1,
#    'lensing_verbose' :1,
#    'output_verbose': 1
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
clphiphi_DCDM_0 = clM_0['pp'][2:]


#For the matter power spectrum of the DCDM, should I also fix h? or it doesn't matter
# since H0 is very close to the reference value?
for k in kk:
    Pk2.append(M.pk(k*h,0.)*h**3) # function .pk(k,z)    
    
for k in kk:
    Pk5.append(M.pk(k*h,3.)*h**3) # function .pk(k,z)
    
    
S8 = derived['sigma8']*np.sqrt(derived['Omega_m']/0.3)

#print("sigma8 for DCDM with epsilon=%.3f and tau= %.0f Gyrs is %f" %(epsilon,tau,derived['sigma8']))
#print("and Omega_m is %f, so that S8=%f "%(derived['Omega_m'],S8))
#print("and z_eq is %f"%derived['z_eq'])

M.struct_cleanup()
M.empty()


fTT_DCDM = interp1d(ll_DCDM_0, ll_DCDM_0*(ll_DCDM_0+1.0)*clTT_DCDM_0/(2.0*np.pi))
fEE_DCDM = interp1d(ll_DCDM_0, ll_DCDM_0*(ll_DCDM_0+1.0)*clEE_DCDM_0/(2.0*np.pi))
fphiphi_DCDM = interp1d(ll_DCDM_0,ll_DCDM_0*(ll_DCDM_0+1.0)*clphiphi_DCDM_0/(2.0*np.pi))



#%% compute LCDM+neutrino_mass


# FIXING THETA_S ###################################################
print("~~~~~computing LCDM+nu~~~~~")
M.set(common_settings)
M.set({
'omega_cdm': 0.1154,
'100*theta_s':1.042059,
#'H0':67.76,
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
clphiphi_DCDM_1 = clM_1['pp'][2:]


S8 = derived['sigma8']*np.sqrt(derived['Omega_m']/0.3)

#print("sigma8 for LCDM (fixed theta_s) with total neutrino mass M_nu=%.2f eV is %f" %(M_ncdm,derived['sigma8']) )
#print("and Omega_m is %f, so that S8=%f"%(derived['Omega_m'],S8))
#print("and z_eq is %f"%derived['z_eq'])

h = M.h()
for k in kk:
    Pk3.append(M.pk(k*h,0.)*h**3) # function .pk(k,z)    
    
for k in kk:
    Pk6.append(M.pk(k*h,3.)*h**3) # function .pk(k,z)



M.struct_cleanup()
M.empty()


fTT_mNU = interp1d(ll_DCDM_1,ll_DCDM_1*(ll_DCDM_1+1.0)*clTT_DCDM_1/(2.0*np.pi))
fEE_mNU = interp1d(ll_DCDM_1,ll_DCDM_1*(ll_DCDM_1+1.0)*clEE_DCDM_1/(2.0*np.pi))
fphiphi_mNU = interp1d(ll_DCDM_1,ll_DCDM_1*(ll_DCDM_1+1.0)*clphiphi_DCDM_1/(2.0*np.pi))


fpk_z0_Mnu  = interp1d(kk,Pk3)
fpk_z3_Mnu  = interp1d(kk,Pk6)



# FIXING H0 ##################################################################

#M.set(common_settings)
#M.set({
#'omega_cdm': 0.1155,
#'omega_cdm':  0.1196,
#'H0':67.58,
#'N_ncdm':1,
#'background_ncdm_distribution': 0,
#'N_ur':0.00641,
#'deg_ncdm':3,
#'m_ncdm':m_ncdm
#})
#M.compute()

#h = M.h()

#print("~~~~~computing LCDM+mNU in %.f s~~~~~"%(time.time()-t_i))

#for k in kk:
#    Pk3.append(M.pk(k*h,0.)*h**3) # function .pk(k,z)    
    
#for k in kk:
#    Pk6.append(M.pk(k*h,3.)*h**3) # function .pk(k,z)

#derived = M.get_current_derived_parameters(['sigma8','Omega_m','z_eq'])

#print("sigma8 for LCDM (fixed H0) with total neutrino mass M_nu=%.1f eV is %f" %(M_ncdm,derived['sigma8']) )
#print("and Omega_m is %f"%derived['Omega_m'])
#print("and z_eq is %f"%derived['z_eq'])
#S8 = derived['sigma8']*np.sqrt(derived['Omega_m']/0.3)

#M.struct_cleanup()
#M.empty()

#%% TIME TO PLOT
print("~~~~~ready to plot~~~~~")


# WITH FLUID APPROXIMATION
ax_1.semilogx(ll_DCDM_0,fTT_DCDM(ll_DCDM_0)/fTT_ref(ll_DCDM_0)-1,'red',label=r'$\Lambda$DDM  Best-fit')
ax_2.semilogx(ll_DCDM_0,fEE_DCDM(ll_DCDM_0)/fEE_ref(ll_DCDM_0)-1,'red',label=r'$\Lambda$DDM  Best-fit')
ax_3.semilogx(ll_DCDM_0,fphiphi_DCDM(ll_DCDM_0)/fphiphi_ref(ll_DCDM_0)-1,'red',label=r'$\Lambda$DDM  Best-fit')

# WITH FULL CALCULATION
#ax_1.semilogx(ll_DCDM_0,fcl_tt_dcdm_full(ll_DCDM_0)/fTT_ref(ll_DCDM_0)-1,'r',label=r'$\Lambda$DDM  Best-fit')
#ax_2.semilogx(ll_DCDM_0,fcl_ee_dcdm_full(ll_DCDM_0)/fEE_ref(ll_DCDM_0)-1,'r',label=r'$\Lambda$DDM  Best-fit')
#ax_3.semilogx(ll_DCDM_0,fcl_pp_dcdm_full(ll_DCDM_0)/fphiphi_ref(ll_DCDM_0)-1,'r',label=r'$\Lambda$DDM  Best-fit')



ax_1.semilogx(ll_DCDM_0,fTT_mNU(ll_DCDM_0)/fTT_ref(ll_DCDM_0)-1,'blue',label= r'$\nu\Lambda\mathrm{CDM} \, \, (M_{\nu} =%.2f \, \mathrm{eV})$'%M_ncdm)
ax_2.semilogx(ll_DCDM_0,fEE_mNU(ll_DCDM_0)/fEE_ref(ll_DCDM_0)-1,'blue',label= r'$\nu\Lambda\mathrm{CDM} \, \, (M_{\nu} =%.2f \, \mathrm{eV})$'%M_ncdm)
ax_3.semilogx(ll_DCDM_0,fphiphi_mNU(ll_DCDM_0)/fphiphi_ref(ll_DCDM_0)-1,'blue',label= r'$\nu\Lambda\mathrm{CDM} \, \, (M_{\nu} =%.2f \, \mathrm{eV}) $'%M_ncdm)

#Planck error bars

# FOR TT
l_cosmic_variance_1 = np.linspace(0,30,1000)
l_cosmic_variance_2 = np.linspace(30,48,2)
slope =np.array([0.15,0.0343])
ax_1.fill_between(l_cosmic_variance_1, -0.15,0.15, color='lightgray' )
ax_1.fill_between(l_cosmic_variance_2, -slope, slope, color='lightgray' )
ax_1.fill_between(lTT, -(DlTT_error_plus)/DlTT_mean, +(DlTT_error_plus)/DlTT_mean, color='lightgray')

# FOR EE
l_cosmic_variance = np.linspace(0,48,1000)
ax_2.fill_between(l_cosmic_variance, -0.21,0.21, color='lightgray' )
ax_2.fill_between(lEE, -(DlEE_error_plus)/DlEE_mean, +(DlEE_error_plus)/DlEE_mean, color = 'lightgray')

#FOR PHI-PHI
l_pp_1 =np.linspace(8,20,10)
l_pp_2 =np.linspace(20,39,10)
l_pp_3 =np.linspace(39,65,10)
l_pp_4 =np.linspace(65,100,10)
l_pp_5 =np.linspace(101,144,10)
l_pp_6 =np.linspace(145,198,10)
l_pp_7 =np.linspace(199,263,10)
l_pp_8 =np.linspace(264,338,10)
l_pp_9 =np.linspace(339,425,10)
l_pp_10=np.linspace(426,525,10)
l_pp_11=np.linspace(526,637,10)
l_pp_12=np.linspace(638,762,10)
l_pp_13=np.linspace(763,901,10)
l_pp_14=np.linspace(902,2048,10)

ax_3.fill_between(l_pp_1, -0.2/1.24, 0.2/1.24, color='lightgray' )
ax_3.fill_between(l_pp_2, -0.11/1.40,0.11/1.40, color='lightgray' )
ax_3.fill_between(l_pp_3, -0.08/1.34,0.08/1.34, color='lightgray' )
ax_3.fill_between(l_pp_4, -0.05/1.14,0.05/1.14, color='lightgray' )
ax_3.fill_between(l_pp_5, -0.05/0.904,0.05/0.904, color='lightgray' )
ax_3.fill_between(l_pp_6, -0.06/0.686,0.06/0.686, color='lightgray' )
ax_3.fill_between(l_pp_7, -0.08/0.513,0.08/0.513, color='lightgray' )
ax_3.fill_between(l_pp_8, -0.10/0.382,0.10/0.382, color='lightgray' )
ax_3.fill_between(l_pp_9, -0.13/0.285,0.13/0.285, color='lightgray' )
ax_3.fill_between(l_pp_10,-0.14/0.213,0.14/0.213, color='lightgray' )
ax_3.fill_between(l_pp_11,-0.19/0.160,0.19/0.160, color='lightgray' )
ax_3.fill_between(l_pp_12,-0.23/0.121,0.23/0.121, color='lightgray' )
ax_3.fill_between(l_pp_13,-0.28/0.0934,0.28/0.0934, color='lightgray' )
ax_3.fill_between(l_pp_14,-0.30/0.0518,0.30/0.0518, color='lightgray' )

#Q: WHAT ARE THE MEAN VALUES BY WHICH I SHOULD DIVIDE THE ERROR BARS?
# SHOULD I "SMOOTH" THE LOOK OF THE BOXES?

ax_3.set_xlabel(r'$\mathrm{multipole} \, \ell$',fontsize=15)
ax_1.set_ylabel(r'$\frac{C_\ell^\mathrm{TT}}{C_\ell^\mathrm{TT}(\mathrm{ref} )} -1$',fontsize=20)
ax_2.set_ylabel(r'$\frac{C_\ell^\mathrm{EE}}{C_\ell^\mathrm{EE}(\mathrm{ref} )} -1$',fontsize=20)
ax_3.set_ylabel(r'$\frac{C_\ell^{\phi \phi} }{C_\ell^{\phi \phi} (\mathrm{ref} )} -1$',fontsize=20)


ax_3.tick_params(axis="x", labelsize=18)
ax_1.tick_params(axis="y", labelsize=18)
ax_2.tick_params(axis="y", labelsize=18)
ax_3.tick_params(axis="y", labelsize=18)


ax_1.legend(frameon=False,fontsize =15,loc='upper left',borderaxespad=0.)

plt.show()
plt.clf()
#%%
#plt.figure(1)

fig, axes = plt.subplots()

axes.set_xscale('log')
axes.set_xlim(kk[0],kk[-1])
axes.set_ylim(-0.4,0.05)


axes.set_xlabel(r'$k \,\,\,\, [h/\mathrm{Mpc}]$', fontsize=15)
axes.set_ylabel(r'$\frac{P(k)}{P_{ref}(k)}-1$', fontsize=20)


#DCDM
#WITH FLUID APPROXIMATION
#axes.plot(kk,map(truediv, list(np.array(Pk2) - np.array(Pk1)), Pk1),'r')
#axes.plot(kk,map(truediv, list(np.array(Pk5) - np.array(Pk4)), Pk4),'r--')

#WITH FULL CALCULATION
axes.plot(kk,fpk_z0_dcdm_full(kk)/fpk_z0_ref(kk)-1.0,'red')
axes.plot(kk,fpk_z3_dcdm_full(kk)/fpk_z3_ref(kk)-1.0,'red',linestyle='dashed')


#LCDM+MNU
#axes.plot(kk,map(truediv, list(np.array(Pk3) - np.array(Pk1)), Pk1),'b')
#axes.plot(kk,map(truediv, list(np.array(Pk6) - np.array(Pk4)), Pk4),'b--')

axes.plot(kk,fpk_z0_Mnu(kk)/fpk_z0_ref(kk)-1.0,'blue')
axes.plot(kk,fpk_z3_Mnu(kk)/fpk_z3_ref(kk)-1.0,'blue',linestyle='dashed')



lines = axes.get_lines()

black_line1 = mlines.Line2D([], [], color='black', linestyle='-', label='z = 0')
black_line2 = mlines.Line2D([], [], color='black', linestyle='--', label='z = 3')

legend1 = plt.legend([lines[i] for i in [0,2]], [r'$\Lambda$DDM  Best-fit', r'$\nu\Lambda\mathrm{CDM} \, \, (M_{\nu} =%.2f \, \mathrm{eV}) $'%M_ncdm], loc='lower left', fontsize=15, frameon=False)
legend2 = plt.legend(handles= [black_line1,black_line2], loc=6, fontsize=15, frameon=False)

axes.add_artist(legend1)
axes.add_artist(legend2)


k_range_desy1 = np.linspace(0.09,0.62,1000) #which are the exact wavenumbers probed by DES-Y1?
k_range_sigma8 = np.linspace(0.1,0.9,1000) #which are the exact wavenumbers probed by DES-Y1?

plt.fill_between(k_range_sigma8, -1.0,0.5, color='lightgray' )

#plt.text(2.0e-4, -0.26, r'$\Omega_{m} =0.31 \, \, \, \, \, \, \, \, \, \, \, \, \sigma_{8} = 0.75  $', fontsize =15)

plt.text(5.0e-3, -0.26, r'$ S_{8} = 0.76  $', fontsize =15)

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

