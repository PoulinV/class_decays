import matplotlib.pyplot as plt
from matplotlib import rc
import numpy as np
from scipy.interpolate import interp1d
from classy import Class
from matplotlib.patches import Rectangle
rc('font',**{'family':'serif','serif':['Times']}, **{'weight':'heavy'})
rc('text', usetex=True)
#%%
#model 1, 30% Matter and 70% Lambda
common_settings = {'write background':'yes',
                   'Omega_b': 0.04,
                   'Omega_cdm':0.26,
                   'H0': 65.2
                   }

M = Class()
M.set(common_settings)
M.compute()
background = M.get_background() # load background table
z_table_1 = background['z'] # read redshift
lum_dist_1 = background['lum. dist.'] #in units of Mpc 
M.struct_cleanup()
M.empty()
lum_dist_1 = interp1d(z_table_1,lum_dist_1)

common_settings = {'write background':'yes',
                   'Omega_b': 0.05,
                   'Omega_cdm':0.95,
                   'H0': 65.2
                   }

#model 2, 100% Matter and 0% Lambda
M = Class()
M.set(common_settings)
M.compute()
background = M.get_background() # load background table
z_table_2 = background['z'] # read redshift
lum_dist_2 = background['lum. dist.'] #in units of Mpc 
M.struct_cleanup()
M.empty()

lum_dist_2 = interp1d(z_table_2,lum_dist_2)


#%%
fig, ax = plt.subplots()
logzc=[3.734,3.779,4.481,4.350,3.896,4.178,3.859,3.685,4.030,3.398,3.547,4.166,3.820,3.679,3.924,4.572,3.891,3.625,4.024,4.130,4.111,4.379,4.418,4.283,3.871,4.189,4.177]
z_low = np.power(10.,logzc)/2.99792e5
mu_low = [34.88,34.77,38.33,37.77,35.59,36.67,35.49,34.34,36.50,32.79,33.73,36.60,35.43,34.00,35.82,39.10,35.53,34.13,36.49,36.87,36.53,37.96,38.09,37.63,35.23,37.31,37.11]
mu_err_low = [0.21,0.15,0.23,0.19,0.16,0.25,0.20,0.14,0.20, 0.16,0.17,0.16, 0.18, 0.18,0.20,0.17,0.20, 0.14,0.21,0.17,0.20,0.15,0.36,0.18,0.25, 0.14, 0.19]

z_high  = [0.43, 0.62,0.57,0.30,0.38,0.43,0.44, 0.50, 0.97,0.48]
mu_high = [41.74, 42.98,42.76,41.38,41.63,42.55,41.95,42.40,44.39,42.45]
mu_err_high = [0.28,0.17,0.19,0.24,0.20,0.25,0.17,0.17,0.30,0.17]

z_theo = np.logspace(-3.0, 1.0, num=1000)
mu_theo_1 = 5.0*np.log10(lum_dist_1(z_theo)) + 25.0
mu_theo_2 = 5.0*np.log10(lum_dist_2(z_theo)) + 25.0

ax.set_xscale('log')
ax.errorbar(z_low, mu_low, yerr=[mu_err_low, mu_err_low],fmt='o',color = 'mediumblue',ecolor='mediumblue',capsize=3,markersize=5.)
ax.errorbar(z_high, mu_high, yerr=[mu_err_high, mu_err_high],fmt='o',color = 'darkgreen',ecolor='darkgreen',capsize=3,markersize=5.)

ax.plot(z_theo,mu_theo_1,color = 'black', linestyle='solid')
ax.plot(z_theo,mu_theo_2,color = 'black', linestyle='dashed')

ax.set_xlabel('redshift z',fontsize=18)
ax.set_ylabel('Distance Modulus m-M',fontsize=18)

ax.text(0.91, 0.88, r'$\Omega_m \ \ \ \Omega_\Lambda$', fontsize=15, transform=plt.gcf().transFigure)
ax.text(0.91, 0.84, r'$0.3 \ \ \ \ 0.7$', fontsize=14, transform=plt.gcf().transFigure)
ax.text(0.91, 0.81, r'$ \ 1.0 \ \ \ \ 0.0$', fontsize=14, transform=plt.gcf().transFigure)


ax.text(0.004,40, 'Low-z',color = 'mediumblue',fontsize=20)
ax.text(0.004,42, 'High-z',color = 'darkgreen',fontsize=20)
ax.tick_params(axis="y", labelsize=16)
ax.tick_params(axis="x", labelsize=16)

ax.add_patch(Rectangle((0.001, 25), 0.0023, 1.5,facecolor = 'red',fill=True))


ax.set_ylim([25,46])
ax.set_xlim([0.002,1.3])
ax.annotate("Hubble (1929)", xy=(0.0035, 26.2), xytext=(0.007, 27), arrowprops=dict(arrowstyle="->",linewidth=2),xycoords='data',textcoords='data',fontsize=15)
plt.show()

#%% #model 2, best-fit LCDM from Planck 2018
common_settings = {'write background':'yes',
                   'omega_b': 0.02233,
                   'omega_cdm':0.1198,
                   '100*theta_s':1.04089,
                   'tau_reio':0.0540, 
                   'N_ur':2.0328, 
                   'N_ncdm':1,
                   'm_ncdm':0.06
#                   'background_ncdm_distribution': '0',
#                   'Quadrature strategy': '0'
                   }
M.compute()
background = M.get_background() # load background table
z_table = background['z'] # read redshift
Omega_cdm    = background['(.)rho_cdm']/background['(.)rho_crit']
Omega_b      = background['(.)rho_b']/background['(.)rho_crit']
Omega_g      = background['(.)rho_g']/background['(.)rho_crit']
Omega_nu     = background['(.)rho_ncdm[0]']/background['(.)rho_crit']
Omega_lambda = background['(.)rho_lambda']/background['(.)rho_crit']
M.struct_cleanup()
M.empty()

Om_cdm_atz= interp1d(z_table,Omega_cdm)
Om_b_atz= interp1d(z_table,Omega_b)
Om_g_atz= interp1d(z_table,Omega_g)
Om_nu_atz= interp1d(z_table,Omega_nu)
Om_lam_atz= interp1d(z_table,Omega_lambda)

z_theo = np.logspace(-8.0, 0.0, num=1000)

plt.tick_params(axis="y", labelsize=16)
plt.tick_params(axis="x", labelsize=16)

plt.xlabel('$\rm{redshift} \ z$',fontsize=18)
plt.ylabel('$\Omega_i (z)$',fontsize=18)

plt.plot(z_theo,Om_cdm_atz(z_theo),color = 'blue', linestyle='solid',label =r'cold dark matter')
plt.plot(z_theo,Om_b_atz(z_theo),color = 'green', linestyle='solid', label =r'baryons')
plt.plot(z_theo,Om_g_atz(z_theo),color = 'red', linestyle='solid', label =r'photons')
plt.plot(z_theo,Om_nu_atz(z_theo),color = 'purple', linestyle='solid', label =r'neutrinos')
plt.plot(z_theo,Om_lam_atz(z_theo),color = 'black', linestyle='solid', label =r'Cosmological cte $\Lambda$')

plt.legend(frameon=False,fontsize =16,loc='lower left',borderaxespad=0., ncol=2)


#plt.ylim([1e-5,1])
#plt.xlim([1e-8,0])

plt.xscale('log')
plt.xscale('log')

plt.show()