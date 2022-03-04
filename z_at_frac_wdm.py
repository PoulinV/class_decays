# import necessary modules
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import interp1d
from classy import Class

#set general configuration
common_settings = {'write background':'yes',
                   'omega_b': 0.02233,
                   'tau_reio':0.0540,
                   '100*theta_s':1.04089,
                   'N_ur':2.0328
                   }

# These LCDM params correspond to Planck 2018 best-fit (TT,TE,EE+lowE+lensing)
# see table 1. in arXiv: 1807.06209v2 

# FIRST MODEL
Log10_Gam_kmsMpc = 1.241 
Log10_eps = -2.16

M = Class()
M.set(common_settings)
M.set({'omega_cdm': 0.00001,
       'omega_ini_dcdm2':0.1198,
       'background_ncdm_distribution': '0,1',
       'Quadrature strategy': '0,4',
       'Number of momentum bins perturbs': '50,50', #for computing S8 is enough with 50 bins instead of 300, since it is unaffected by the small (super-horizon) k
       'back_integration_stepsize': 7e-3, #default is 7e-3, need 1e-3 to make curves look smooth
       'M_dcdm': 1,
       'N_ncdm': 2,
       'm_ncdm':'0.06,0',
       'Log10_Gamma_dcdm':Log10_Gam_kmsMpc ,
       'log10_epsilon_dcdm':Log10_eps})

M.compute()
background = M.get_background() # load background table
z_table = background['z'] # read redshift
rho_wdm= background['(.)rho_ncdm[1]']
w_wdm= background['(.)p_ncdm[1]']/background['(.)rho_ncdm[1]']
rho_wdm = rho_wdm*(1.-3.*w_wdm)
rho_dcdm= background['(.)rho_dcdm']
frac_wdm_table= rho_wdm/(rho_wdm+rho_dcdm)
z_wdm= interp1d(100.0*frac_wdm_table,z_table)
z_at_fwdm_10 = z_wdm(10.0)
print('Model with Log10_Gamma_dcdm=%f,log10_epsilon_dcdm=%f has z_{10} = %f \n'%(Log10_Gam_kmsMpc,Log10_eps,z_at_fwdm_10))
M.struct_cleanup()
M.empty()

# SECOND MODEL
Log10_Gam_kmsMpc = 2.1
Log10_eps = -2.8

M = Class()
M.set(common_settings)
M.set({'omega_cdm': 0.00001,
       'omega_ini_dcdm2':0.1198,
       'background_ncdm_distribution': '0,1',
       'Quadrature strategy': '0,4',
       'Number of momentum bins perturbs': '50,50', #for computing S8 is enough with 50 bins instead of 300, since it is unaffected by the small (super-horizon) k
       'back_integration_stepsize': 7e-3, #default is 7e-3, need 1e-3 to make curves look smooth
       'M_dcdm': 1,
       'N_ncdm': 2,
       'm_ncdm':'0.06,0',
       'Log10_Gamma_dcdm':Log10_Gam_kmsMpc ,
       'log10_epsilon_dcdm':Log10_eps})

M.compute()
background = M.get_background() # load background table
z_table = background['z'] # read redshift
rho_wdm= background['(.)rho_ncdm[1]']
w_wdm= background['(.)p_ncdm[1]']/background['(.)rho_ncdm[1]']
rho_wdm = rho_wdm*(1.-3.*w_wdm)
rho_dcdm= background['(.)rho_dcdm']
frac_wdm_table= rho_wdm/(rho_wdm+rho_dcdm)
z_wdm= interp1d(100.0*frac_wdm_table,z_table)
z_at_fwdm_10 = z_wdm(10.0)
print('Model with Log10_Gamma_dcdm=%f,log10_epsilon_dcdm=%f has z_{10} = %f \n'%(Log10_Gam_kmsMpc,Log10_eps,z_at_fwdm_10))
M.struct_cleanup()
M.empty()

