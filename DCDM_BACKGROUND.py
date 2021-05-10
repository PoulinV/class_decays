
# import necessary modules
import matplotlib.pyplot as plt
import numpy as np
from classy import Class
from scipy.interpolate import interp1d
import matplotlib.lines as mlines
    
plt.rcParams["figure.figsize"] = [6.0,15.0]

import time
start_time = time.time()


#tau_gyrs=30 # and 1
#Gamma= 1/(tau_gyrs*1.02e-3)
#Gamma


conversion=2.9979e5 #To convert H from Mpc^-1 to km/s/Mpc

#%%
#set general configuration
common_settings = {'write background':'yes',
                   'write thermodynamics':'yes',
                   'omega_b': 0.0224,
                   'tau_reio':0.0582,
                   'H0': 67.7,
#                   '100*theta_s':1.042168,
                   'N_ur':2.0328
                   }



#LCDM
M = Class()
M.set(common_settings)
M.set({'omega_cdm': 0.1194,
       'background_ncdm_distribution':0,
       'Quadrature strategy': 0,
       'Number of momentum bins perturbs':50,
       'N_ncdm': 1,
       'm_ncdm':0.06})

M.compute()
background = M.get_background() # load background table
background_z = background['z'] # read redshift
back_Omega_cdm= background['(.)rho_cdm']/background['(.)rho_crit']
back_H= background['H [1/Mpc]']
r_of_z = background['comov. dist.']
Om_lcdm = M.Omega_m()
H0_lcdm = 100.0*M.h() # in units of km/s/Mpc
M.struct_cleanup()
M.empty()


Hubble_lcdm = interp1d(background_z,conversion*back_H)
Om_cdm = interp1d(background_z,back_Omega_cdm)
r_of_z_lcdm = interp1d(background_z,r_of_z)



# DCDM MODEL 1

eps1 =0.499999
Gamma1 = 326.797
tau= 1/(Gamma1*1.02e-3)

M = Class()
M.set(common_settings)
M.set({'omega_ini_dcdm2':0.1194,
       'omega_cdm': 0.00001,
       'background_ncdm_distribution': '0,1',
       'Quadrature strategy': '0,4',
       'Number of momentum bins perturbs': '50,300',
#       'back_integration_stepsize': 1e-3, #default is 7e-3, need 1e-3 to make curves look smooth
       'back_integration_stepsize': 7e-3, 
       'M_dcdm': 1,
       'N_ncdm': 2,
       'm_ncdm':'0.06,0',
       'Gamma_dcdm': Gamma1, #lifetime 3 Gyrs
       'epsilon_dcdm': eps1,
   
       })
M.compute()

background = M.get_background() # load background table
background_z = background['z'] # read redshift
#back_Omega_dcdm= background['(.)rho_dcdm']/background['(.)rho_crit']
#back_Omega_dr= background['(.)rho_dr']/background['(.)rho_crit']
#back_Omega_wdm= background['(.)rho_ncdm[1]']/background['(.)rho_crit']
back_H= background['H [1/Mpc]']
r_of_z = background['comov. dist.']
Om_dcdm_1 = M.Omega_m()
H0_dcdm_1 = 100.0*M.h() # in units of km/s/Mpc
M.struct_cleanup()
M.empty()



Hubble_dcdm_1 = interp1d(background_z,conversion*back_H)

r_of_z_dcdm_1 = interp1d(background_z,r_of_z)


# DCDM MODEL 2

M = Class()
M.set(common_settings)
M.set({'omega_ini_dcdm2':0.1194,
       'omega_cdm': 0.00001,
       'background_ncdm_distribution': '0,1',
       'Quadrature strategy': '0,4',
       'Number of momentum bins perturbs': '50,300',
       'back_integration_stepsize': 7e-3, #default is 7e-3, need 1e-3 to make curves look smooth
       'M_dcdm': 1,
       'N_ncdm': 2,
       'm_ncdm':'0.06,0',
       'Gamma_dcdm': 3267.973, #lifetime 0.3 Gyrs
       'epsilon_dcdm': 0.499999,  
       })
M.compute()

background = M.get_background() # load background table
background_z = background['z'] # read redshift
#back_Omega_dcdm= background['(.)rho_dcdm']/background['(.)rho_crit']
#back_Omega_dr= background['(.)rho_dr']/background['(.)rho_crit']
#back_Omega_wdm= background['(.)rho_ncdm[1]']/background['(.)rho_crit']
back_H= background['H [1/Mpc]']
r_of_z = background['comov. dist.']
Om_dcdm_2 = M.Omega_m()
H0_dcdm_2 = 100.0*M.h() # in units of km/s/Mpc
M.struct_cleanup()
M.empty()


Hubble_dcdm_2 = interp1d(background_z,conversion*back_H)

r_of_z_dcdm_2 = interp1d(background_z,r_of_z)



# DCDM MODEL 3

eps3 = 0.2

Gamma3 = Gamma1
M = Class()
M.set(common_settings)
M.set({'omega_ini_dcdm2':0.1194,
       'omega_cdm': 0.00001,
       'background_ncdm_distribution': '0,1',
       'Quadrature strategy': '0,4',
       'Number of momentum bins perturbs': '50,300',
#       'back_integration_stepsize': 1e-3, #default is 7e-3, need 1e-3 to make curves look smooth
       'back_integration_stepsize': 7e-3, 
       'M_dcdm': 1,
       'N_ncdm': 2,
       'm_ncdm':'0.06,0',
       'Gamma_dcdm': Gamma3, #lifetime 3 Gyrs
       'epsilon_dcdm': eps3,     
       })
M.compute()

background = M.get_background() # load background table
background_z = background['z'] # read redshift
#back_Omega_dcdm= background['(.)rho_dcdm']/background['(.)rho_crit']
#back_Omega_dr= background['(.)rho_dr']/background['(.)rho_crit']
#back_Omega_wdm= background['(.)rho_ncdm[1]']/background['(.)rho_crit']
back_H= background['H [1/Mpc]']
r_of_z = background['comov. dist.']
Om_dcdm_3 = M.Omega_m()
H0_dcdm_3 = 100.0*M.h() # in units of km/s/Mpc
M.struct_cleanup()
M.empty()


Hubble_dcdm_3 = interp1d(background_z,conversion*back_H)

r_of_z_dcdm_3 = interp1d(background_z,r_of_z)


# DCDM MODEL 4

M = Class()
M.set(common_settings)
M.set({'omega_ini_dcdm2':0.1194,
       'omega_cdm': 0.00001,
       'background_ncdm_distribution': '0,1',
       'Quadrature strategy': '0,4',
       'Number of momentum bins perturbs': '50,300',
#       'back_integration_stepsize': 1e-3, #default is 7e-3, need 1e-3 to make curves look smooth
       'back_integration_stepsize': 7e-3, 
       'M_dcdm': 1,
       'N_ncdm': 2,
       'm_ncdm':'0.06,0',
      'Gamma_dcdm': 3267.973, #lifetime 0.3 Gyrs
      'epsilon_dcdm': 0.2    
       })
M.compute()

background = M.get_background() # load background table
background_z = background['z'] # read redshift
#back_Omega_dcdm= background['(.)rho_dcdm']/background['(.)rho_crit']
#back_Omega_dr= background['(.)rho_dr']/background['(.)rho_crit']
#back_Omega_wdm= background['(.)rho_ncdm[1]']/background['(.)rho_crit']
back_H= background['H [1/Mpc]']
r_of_z = background['comov. dist.']
Om_dcdm_4 = M.Omega_m()
H0_dcdm_4 = 100.0*M.h() # in units of km/s/Mpc
M.struct_cleanup()
M.empty()


Hubble_dcdm_4 = interp1d(background_z,conversion*back_H)

r_of_z_dcdm_4 = interp1d(background_z,r_of_z)


# DCDM MODEL 5

M = Class()
M.set(common_settings)
M.set({'omega_ini_dcdm2':0.1194,
       'omega_cdm': 0.00001,
       'background_ncdm_distribution': '0,1',
       'Quadrature strategy': '0,4',
       'Number of momentum bins perturbs': '50,300',
       'back_integration_stepsize': 7e-3, #default is 7e-3, need 1e-3 to make curves look smooth
       'M_dcdm': 1,
       'N_ncdm': 2,
       'm_ncdm':'0.06,0',
       'Gamma_dcdm': 32.679, #lifetime 30 Gyrs
       'epsilon_dcdm': 0.1})

M.compute()

background = M.get_background() # load background table
background_z = background['z'] # read redshift
back_Omega_dcdm= background['(.)rho_dcdm']/background['(.)rho_crit']
back_Omega_dr= background['(.)rho_dr']/background['(.)rho_crit']
back_Omega_wdm= background['(.)rho_ncdm[1]']/background['(.)rho_crit']

back_Omega_wdm= background['(.)rho_ncdm[1]']/background['(.)rho_crit']

back_w_wdm= background['(.)p_ncdm[1]']/background['(.)rho_ncdm[1]']

#back_H= background['H [1/Mpc]']
M.struct_cleanup()
M.empty()



Om_dcdm = interp1d(background_z,back_Omega_dcdm)
Om_dr = interp1d(background_z,back_Omega_dr)
Om_wdm = interp1d(background_z,back_Omega_wdm)
w_wdm = interp1d(background_z,back_w_wdm)


#%%
# time to plot omegas and Hubble parameter

zz=np.linspace(0,1e4, 1e6)
#zzz=np.linspace(0,2.5, 1e4)
zzz=np.linspace(0,20, 1e4)


ax_1=plt.subplot(211)
ax_2=plt.subplot(212)
plt.subplots_adjust(hspace=0.25)

ax_1.set_xlim([1e-1,1e4])
ax_1.set_ylim([1e-4,10])

ax_1.loglog(zz,Om_cdm(zz), 'black', label=r'CDM')
ax_1.loglog(zz,Om_dcdm(zz), 'blue', label=r'DCDM')
ax_1.loglog(zz,Om_wdm(zz), 'green', label=r'WDM')
ax_1.loglog(zz,Om_dr(zz), 'red', label=r'DR')

ax_1.set_xlabel(r'$z$', fontsize=18)
ax_1.set_ylabel(r'$\Omega (z)$', fontsize=18)
ax_1.legend(loc='lower right', frameon=False,borderaxespad=0., fontsize=13)


ax_1.tick_params(axis='x',labelsize=17)
ax_1.tick_params(axis='y',labelsize=17)

ax_1.text(100,3,r'$\Gamma^{-1}= 30 \ \mathrm{Gyrs} \ \ \ \ \ \varepsilon=0.1$', fontsize=13)



ax_2.set_xlim([0.05,20])
#ax_2.set_ylim([35,160])


ax_2.semilogx(zzz,Hubble_dcdm_3(zzz)/(1.+zzz) , 'green', linestyle= 'dashed')
ax_2.semilogx(zzz,Hubble_dcdm_4(zzz)/(1.+zzz) , 'red', linestyle= 'dashed')
ax_2.semilogx(zzz,Hubble_dcdm_1(zzz)/(1.+zzz) , 'green', linestyle= 'dashdot')
ax_2.semilogx(zzz,Hubble_dcdm_2(zzz)/(1.+zzz) , 'red',linestyle='dashdot')
ax_2.semilogx(zzz,Hubble_lcdm(zzz)/(1.+zzz) , 'black', label=r'$\Lambda \mathrm{CDM}$')





lines = ax_2.get_lines()

black_line1 = mlines.Line2D([], [], color='black', linestyle='dashed', label=r'$\varepsilon=0.2$')
black_line2 = mlines.Line2D([], [], color='black', linestyle='dashdot', label=r'$ \varepsilon=0.5$')
black_line3 = mlines.Line2D([], [], color='black', linestyle='solid', label=r'$\Lambda$CDM')

color_line1 = mlines.Line2D([], [], color='green', linestyle='solid', label=r'$ \Gamma^{-1} = 3 \ \mathrm{Gyrs} $')
color_line2 = mlines.Line2D([], [], color='red', linestyle='solid', label=r'$\Gamma^{-1} = 0.3 \ \mathrm{Gyrs} $')



legend1 = plt.legend(handles= [black_line3, color_line1,color_line2], loc=(0.1,0.65), fontsize=13, frameon=False)
legend2 = plt.legend(handles= [black_line1,black_line2], loc=(0.4,0.65), fontsize=13, frameon=False)


ax_2.add_artist(legend1)
ax_2.add_artist(legend2)

ax_2.set_xlabel(r'$z$', fontsize=18)
ax_2.set_ylabel(r'$H(z)/(1+z) \ \ \ \ [\mathrm{km}/\mathrm{s}/\mathrm{Mpc}]$', fontsize=18)


#ax_2.legend(loc=(0.15,0.45), frameon=False, borderaxespad=0., fontsize=13)

ax_2.tick_params(axis='x',labelsize=17)
ax_2.tick_params(axis='y',labelsize=17)


plt.show()
plt.clf()


#%% plot the kernel of the lensing potential in terms of redshift


zzz=np.linspace(0,10, 1e4)

r_lss = 13870 #conformal distance to last scattering surface, in units of Mpc

plt.yscale('log')


#print(H0_lcdm)
#print(H0_dcdm_3)
#print(H0_dcdm_1)

kernel_lcdm = (H0_lcdm**2*(3./2.)*Om_lcdm*(1.0+zzz)*(1.0 - (r_of_z_lcdm(zzz)/r_lss)))**2/Hubble_lcdm(zzz)
kernel_dcdm_3 = (H0_dcdm_3**2*(3./2.)*Om_dcdm_3*(1.0+zzz)*(1.0 - (r_of_z_dcdm_3(zzz)/r_lss)))**2/Hubble_dcdm_3(zzz)
kernel_dcdm_1 = (H0_dcdm_1**2*(3./2.)*Om_dcdm_1*(1.0+zzz)*(1.0 - (r_of_z_dcdm_1(zzz)/r_lss)))**2/Hubble_dcdm_1(zzz)


#kernel_lcdm = (H0_lcdm**2*(3./2.)*Om_lcdm*(1.0+zzz)*r_of_z_lcdm(zzz)*(1.0 - (r_of_z_lcdm(zzz)/r_lss)))/Hubble_lcdm(zzz)
#kernel_dcdm_3 = (H0_dcdm_3**2*(3./2.)*Om_dcdm_3*(1.0+zzz)*r_of_z_dcdm_3(zzz)*(1.0 - (r_of_z_dcdm_3(zzz)/r_lss)))/Hubble_dcdm_3(zzz)
#kernel_dcdm_1 = (H0_dcdm_1**2*(3./2.)*Om_dcdm_1*(1.0+zzz)*r_of_z_dcdm_1(zzz)*(1.0 - (r_of_z_dcdm_1(zzz)/r_lss)))/Hubble_dcdm_1(zzz)

plt.plot(zzz,kernel_lcdm, 'black',label=r'$\Lambda \mathrm{CDM}$')
plt.plot(zzz,kernel_dcdm_3 , 'blue',label=r'$\varepsilon=%.1f$'%eps3)
plt.plot(zzz,kernel_dcdm_1, 'red', label=r'$\varepsilon=%.1f$'%eps1)



plt.legend(loc='upper right', frameon=False,borderaxespad=0., fontsize=15)


plt.text(8,1e4,r'$\Gamma^{-1}= %.1f \ \mathrm{Gyrs} $'%tau, fontsize=15)

plt.title(r'$W(z) = \left(3 H_0^2 \Omega_m (1+z)/2 \right)^2 \left(1-r(z)/r_{\ast}\right)^2 / H(z)$',fontsize=15, pad=15)

#plt.title(r'$W(z) = \left(3 H_0^2 \Omega_m (1+z)/2 \right) r(z) \left(1-r(z)/r_{\ast}\right)^2 / H(z)$',fontsize=15, pad=15)

plt.xlabel(r'$z$', fontsize=18)

plt.tick_params(axis='x',labelsize=17)
plt.tick_params(axis='y',labelsize=17)

#OK, SO LENSING KERNEL INDEED BECOMES MORE SUPPRESSED FOR BIGGER VALUES OF EPSILON
# THE KEY FACTOR SEEMS TO BE THE REDUCTION IN OMEGA_M (WITHOUT THIS FACTOR, THE KERNEL COULD EVEN INCREASE WITH RESPECT LCDM)
# NOW I SHOULD MAKE SURE THIS KEEPS BEING THE CASE EVEN IF WE FIX THETA_S, INSTEAD OF H0
# (include H0 in the definition of lensing kernel )

plt.show()
plt.clf()

