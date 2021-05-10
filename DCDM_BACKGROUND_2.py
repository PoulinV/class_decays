
# import necessary modules
import matplotlib.pyplot as plt
import numpy as np
from classy import Class
from scipy.interpolate import interp1d

from sklearn.linear_model import LinearRegression
    
plt.rcParams["figure.figsize"] = [6.0,15.0]

import time
start_time = time.time()

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

# DCDM MODEL 1


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

back_w_wdm= background['(.)p_ncdm[1]']/background['(.)rho_ncdm[1]']

M.struct_cleanup()
M.empty()


w_wdm_1 = interp1d(background_z,back_w_wdm)

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
       'Gamma_dcdm': 32.679, #lifetime 30 Gyrs
       'epsilon_dcdm': 0.01})

M.compute()

background = M.get_background() # load background table
background_z = background['z'] # read redshift

back_w_wdm= background['(.)p_ncdm[1]']/background['(.)rho_ncdm[1]']

M.struct_cleanup()
M.empty()


w_wdm_2 = interp1d(background_z,back_w_wdm)

# DCDM MODEL 3

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
       'epsilon_dcdm': 0.001})

M.compute()

background = M.get_background() # load background table
background_z = background['z'] # read redshift

back_w_wdm= background['(.)p_ncdm[1]']/background['(.)rho_ncdm[1]']

M.struct_cleanup()
M.empty()


w_wdm_3 = interp1d(background_z,back_w_wdm)


# DCDM MODEL 4

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
       'epsilon_dcdm': 0.005})

M.compute()

background = M.get_background() # load background table
background_z = background['z'] # read redshift

back_w_wdm= background['(.)p_ncdm[1]']/background['(.)rho_ncdm[1]']

M.struct_cleanup()
M.empty()


w_wdm_4 = interp1d(background_z,back_w_wdm)


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
       'epsilon_dcdm': 0.05})

M.compute()

background = M.get_background() # load background table
background_z = background['z'] # read redshift

back_w_wdm= background['(.)p_ncdm[1]']/background['(.)rho_ncdm[1]']

M.struct_cleanup()
M.empty()


w_wdm_5 = interp1d(background_z,back_w_wdm)





#%% plot equation of state parameter
zzz=np.linspace(0,5, 1e4)

plt.plot(zzz,w_wdm_1(zzz),'green',label= r'$\varepsilon =0.1$')
plt.plot(zzz,w_wdm_5(zzz),'black',label= r'$\varepsilon =0.05$')
plt.plot(zzz,w_wdm_2(zzz),'blue',label= r'$\varepsilon =0.01$')
plt.plot(zzz,w_wdm_4(zzz),'purple',label= r'$\varepsilon =0.005$')
plt.plot(zzz,w_wdm_3(zzz),'red',label= r'$\varepsilon =0.001$')

plt.yscale('log')

plt.ylabel(r'$w_{\mathrm{wdm}}$', fontsize=18)
plt.xlabel(r'$z$', fontsize=18)

plt.legend(loc='best', frameon=False,borderaxespad=0., fontsize=15)

plt.title(r'$\Gamma^{-1} =  30 \ \mathrm{Gyrs}$', fontsize =18)

plt.tick_params(axis='x',labelsize=17)
plt.tick_params(axis='y',labelsize=17)

plt.show()
plt.clf()


#%%

eps = np.array([0.1,0.05, 0.01, 0.005, 0.001 ])
w_inf = np.array([1.74e-3,4.1e-4,1.46e-5,3.66e-6,1.5e-7])



plt.plot(eps,w_inf , 'o')

plt.yscale('log')
plt.xscale('log')


plt.xlabel(r'$\varepsilon$', fontsize=18)
plt.ylabel(r'$w_{\mathrm{wdm}}(z_{\mathrm{ini}})$', fontsize=18)

plt.title(r'$\Gamma^{-1} =  30 \ \mathrm{Gyrs}$', fontsize =18)


plt.tick_params(axis='x',labelsize=17)
plt.tick_params(axis='y',labelsize=17)


eps = eps.reshape(-1,1)
w_inf = w_inf.reshape(-1,1)

model = LinearRegression()
reg = model.fit(np.log(eps),np.log(w_inf))

print('Intercept: \n',reg.intercept_)
print('Coefficients: \n', reg.coef_)

### NICE, WE GET A SLOPE OF 2, AS WE EXPECTED (FOR A GIVEN LAMBDA AND SMALL VALUES OF EPSILON (EQUAL OR SMALLER THAN 0.1) )

a=np.exp(reg.intercept_)
b = reg.coef_
plt.plot(eps,a*eps**b)


plt.show()
plt.clf()

#%%

eps2 = np.array([0.1,0.05,0.01,0.001])
k_fss = np.array([0.006840,0.014442,0.075168,0.758246]) #all for a lifetime of 30 Gyrs

plt.plot(eps2, k_fss, 'o')

plt.yscale('log')
plt.xscale('log')

plt.xlabel(r'$\varepsilon$', fontsize=18)
plt.ylabel(r'$k_{\mathrm{fs}} \ \ \ [\mathrm{Mpc}^{-1}]$', fontsize=18)


plt.title(r'$\Gamma^{-1} =  30 \ \mathrm{Gyrs}$', fontsize =18)


plt.tick_params(axis='x',labelsize=17)
plt.tick_params(axis='y',labelsize=17)


eps2 = eps2.reshape(-1,1)
k_fss = k_fss.reshape(-1,1)

model = LinearRegression()
reg = model.fit(np.log(eps2),np.log(k_fss))

print('Intercept: \n',reg.intercept_)
print('Coefficients: \n', reg.coef_)

a=np.exp(reg.intercept_)
b = reg.coef_

plt.plot(eps2,a*eps2**b)

plt.show()
plt.clf()

