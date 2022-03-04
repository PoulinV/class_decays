
# import necessary modules
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import rc
import numpy as np
from scipy.interpolate import interp1d
from scipy.interpolate import interp2d
import scipy.optimize
from classy import Class


rc('font',**{'family':'serif','serif':['Times']})
rc('text', usetex=True)
font = {'size': 16, 'family': 'STIXGeneral'}
axislabelfontsize='large'
matplotlib.rc('font', **font)
#matplotlib.mathtext.rcParams['legend.fontsize']='medium'
#plt.rcParams["figure.figsize"] = [8.0,6.0]
plt.rcParams["figure.figsize"] = [5.0,10.0]

log_Gam_1,log_eps_1 = np.loadtxt("LogG_vs_Logeps_baseline_S8_BOSS_1.txt",unpack=True)
log_Gam_2,log_eps_2 = np.loadtxt("LogG_vs_Logeps_baseline_S8_BOSS_2.txt",unpack=True)

log_Gam_1 = log_Gam_1 -3.0
log_Gam_2 = log_Gam_2 -3.0
#%%

z_obs = 0.32

Log10_eps        = np.array([-1.3,-1.6,-1.905,-2.16,-2.5,-2.8,-3.1])
Log10_Gam_kmsMpc = np.array([0.0,0.3,0.6,0.9309,1.241,1.5,1.8,2.1])


##NOTE: Log10_eps = -2.27 seems to give k_fs_wdm = 0.2 h/Mpc, for roughly all the Gamma

frac_wdm_at_z = np.zeros([len(Log10_eps),len(Log10_Gam_kmsMpc)])
z_at_frac_wdm = np.zeros([len(Log10_eps),len(Log10_Gam_kmsMpc)])
k_fs_wdm = np.zeros([len(Log10_eps),len(Log10_Gam_kmsMpc)])

Log10eps_for_kfs0p7 = np.zeros(len(Log10_Gam_kmsMpc))
Log10eps_for_kfs0p2 = np.zeros(len(Log10_Gam_kmsMpc))


S_8 = np.zeros([len(Log10_eps),len(Log10_Gam_kmsMpc)])

#set general configuration
common_settings = {'output':'mPk',
                   'P_k_max_h/Mpc':1.0,
                   'k_output_values':1.0,
                   'write background':'yes',
#                   'omega_b': 0.02233,
#                   'tau_reio':0.0540,
#                   'n_s': 0.9652,
#                   'ln10^{10}A_s': 3.043,
#                   '100*theta_s':1.04089,
                   'omega_b': 0.02242,
                   'tau_reio':0.0571,
                   'n_s': 0.9676,
                   'ln10^{10}A_s': 3.049,
                   'H0':67.72,
                   'N_ur':2.0328
                   }
# These LCDM params correspond to Planck 2018 best-fit (TT,TE,EE+lowE+lensing)
# see table 1. in arXiv: 1807.06209v2 
# changing them to the LCDM best-fit params of the DCDM model (in combined analysis), reduces S8 slightly
count = 0
for i in range(len(Log10_eps)):
    for j in range(len(Log10_Gam_kmsMpc)):       
        M = Class()
        M.set(common_settings)
        M.set({'Omega_cdm': 0.000001,
               'Omega_ini_dcdm2':0.26028,
#               'omega_ini_dcdm2':0.1198,
               'background_ncdm_distribution': '0,1',
               'Quadrature strategy': '0,4',
               'Number of momentum bins perturbs': '25,25', #for computing S8 is enough with 50 bins instead of 300, since it is unaffected by the small (super-horizon) k
               'back_integration_stepsize': 7e-3, #default is 7e-3, need 1e-3 to make curves look smooth
               'M_dcdm': 1,
               'N_ncdm': 2,
               'm_ncdm':'0.06,0',
               'evolver': 0,
               'ncdm_fluid_approximation': 2,
               'ncdm_fluid_trigger_tau_over_tau_k': 25,
               'massive_daughter_perturbations': 'yes',
               'dark_radiation_perturbations': 'yes',
               'Log10_Gamma_dcdm': Log10_Gam_kmsMpc[j],
               'log10_epsilon_dcdm': Log10_eps[i]})

        M.compute()
        background = M.get_background() # load background table
        z_table_back = background['z'] # read redshift
        rho_wdm= background['(.)rho_ncdm[1]']
        w_wdm= background['(.)p_ncdm[1]']/background['(.)rho_ncdm[1]']
        rho_wdm = rho_wdm*(1.-3.*w_wdm)
        rho_dcdm= background['(.)rho_dcdm']
        frac_wdm_table= rho_wdm/(rho_wdm+rho_dcdm)
        frac_wdm= interp1d(z_table_back,frac_wdm_table)
        frac_wdm_at_z[i,j] = 100*frac_wdm(z_obs)  
        derived = M.get_current_derived_parameters(['sigma8', 'Omega_m'])
        S_8[i,j] = derived['sigma8']*np.sqrt(derived['Omega_m']/0.3)
        all_k = M.get_perturbations()
        h = M.h()
        one_k = all_k['scalar'][0]
        z_table_perts = -1.0+1.0/one_k['a']
        kfs_wdm_table = one_k['k_fss_wdm[1]']
        kfs_wdm_atz = interp1d(z_table_perts,kfs_wdm_table)
        k_fs_wdm[i,j] = kfs_wdm_atz(0.0)/h 
        z_wdm= interp1d(100.0*frac_wdm_table,z_table_back)
        z_at_frac_wdm[i,j] = z_wdm(1.0)
        count = count +1
        print('Model Number=%d, Log10_Gamma_kmsMpc=%f,log10_epsilon_dcdm=%f,frac_wdm=%f, z_wdm=%f S8 =%f \n'%(count,Log10_Gam_kmsMpc[j],Log10_eps[i],frac_wdm_at_z[i,j],z_at_frac_wdm[i,j],S_8[i,j]))
        M.struct_cleanup()
        M.empty()
        
#compute values of Gamma and epsilon for which k_fs_wdm = k_max 
kfs_interp = interp2d(Log10_Gam_kmsMpc-3.0,Log10_eps, k_fs_wdm)

for i in range(len(Log10_Gam_kmsMpc)):
    def fun1(x):
        return kfs_interp(Log10_Gam_kmsMpc[i]-3.0,x)-0.7
    
    def fun2(x):
        return kfs_interp(Log10_Gam_kmsMpc[i]-3.0,x)-0.2
    
    Log10eps_for_kfs0p7[i] = scipy.optimize.fsolve(fun1,-2.27)
    Log10eps_for_kfs0p2[i] = scipy.optimize.fsolve(fun2,-2.27)
   

delta_x = abs(Log10_Gam_kmsMpc[1]-Log10_Gam_kmsMpc[0])
delta_y = abs(Log10_eps[0]-Log10_eps[1])

x_start = Log10_Gam_kmsMpc[0]-(delta_x/2.)-3.0
x_end = Log10_Gam_kmsMpc[len(Log10_Gam_kmsMpc)-1]+(delta_x/2.)-3.0
y_start = Log10_eps[len(Log10_eps)-1]-(delta_y/2.)
y_end = Log10_eps[0]+(delta_y/2.)

extent = [x_start, x_end, y_start, y_end]
#%% plot fraction of WDM

fig = plt.figure()
ax = fig.add_subplot(111)
ax.set_xlabel(r'$\mathrm{Log}_{10}(\Gamma / [\mathrm{Gyr}^{-1}])$', fontsize =18)
ax.set_ylabel(r'$\mathrm{Log}_{10}(\varepsilon)$', fontsize = 18)
ax.set_title(r'$z \ \ \rm{at} \ \ \rm{which} \ \ \bar{\rho}_{\mathrm{wdm}}/(\bar{\rho}_{\mathrm{wdm}}+\bar{\rho}_{\mathrm{dcdm}}) = 1 \ \% $', fontsize = 22, pad=15)

ax.xaxis.set_major_locator(plt.MaxNLocator(len(Log10_Gam_kmsMpc)))
ax.yaxis.set_major_locator(plt.MaxNLocator(len(Log10_eps)))

ax.set_yticks([-1.3,-1.6,-1.9,-2.2,-2.5,-2.8,-3.1])
ax.set_yticklabels(['$-1.3$','$-1.6$','$-1.9$','$-2.2$','$-2.5$','$-2.8$','$-3.1$'])

ax.set_xticks([-3.0,-2.7,-2.4,-2.1,-1.8,-1.5,-1.2,-0.9])
ax.set_xticklabels(['$-3.0$','$-2.7$','$-2.4$','$-2.1$','$-1.8$','$-1.5$','$-1.2$','$-0.9$'])


ax.plot(log_Gam_1,log_eps_1,color='black',linestyle='solid')
ax.plot(log_Gam_2,log_eps_2,color='black',linestyle='solid')

ax.plot(Log10_Gam_kmsMpc-3.0,Log10eps_for_kfs0p7, color= 'red',linestyle='dashed')
ax.plot([-3.3,-3.0],[Log10eps_for_kfs0p7[0],Log10eps_for_kfs0p7[0]], color= 'red',linestyle='dashed')
ax.plot([-0.9,-0.6],[Log10eps_for_kfs0p7[len(Log10_Gam_kmsMpc)-1],Log10eps_for_kfs0p7[len(Log10_Gam_kmsMpc)-1]+0.07], color= 'red',linestyle='dashed')
ax.annotate("", xy=(0.75-3.0,Log10eps_for_kfs0p7[2]-0.15), xytext=(0.75-3.0,Log10eps_for_kfs0p7[2]), arrowprops=dict(arrowstyle="->",linewidth=2,color= 'red'),xycoords='data',textcoords='data',fontsize=15)
ax.text(0.8-3.0,Log10eps_for_kfs0p7[2]-0.1, r'$k_{\mathrm{fs}} > k_{\mathrm{NL}}$',color= 'red',fontsize=14 )

ax.plot(Log10_Gam_kmsMpc-3.0,Log10eps_for_kfs0p2, color= 'darkorange',linestyle='dashed')
ax.plot([-3.3,-3.0],[Log10eps_for_kfs0p2[0],Log10eps_for_kfs0p2[0]], color= 'darkorange',linestyle='dashed')
ax.plot([-0.9,-0.6],[Log10eps_for_kfs0p2[len(Log10_Gam_kmsMpc)-1],Log10eps_for_kfs0p2[len(Log10_Gam_kmsMpc)-1]+0.07], color= 'darkorange',linestyle='dashed')
ax.annotate("", xy=(0.75-3.0,Log10eps_for_kfs0p2[4]-0.15), xytext=(0.75-3.0,Log10eps_for_kfs0p2[4]), arrowprops=dict(arrowstyle="->",linewidth=2,color= 'darkorange'),xycoords='data',textcoords='data',fontsize=15)
ax.text(0.8-3.0,Log10eps_for_kfs0p2[4]-0.1, r'$k_{\mathrm{fs}} > k_{\mathrm{max}}$',color= 'darkorange',fontsize=14 )


im = ax.imshow(z_at_frac_wdm, extent=extent, origin='upper', interpolation='bicubic', cmap='GnBu')

    
jump_x = (x_end - x_start) / (2.0 * len(Log10_Gam_kmsMpc))
jump_y = (y_end - y_start) / (2.0 * len(Log10_eps))
x_positions = np.linspace(start=x_start, stop=x_end, num=len(Log10_Gam_kmsMpc), endpoint=False)
y_positions = np.linspace(start=y_end, stop=y_start, num=len(Log10_eps), endpoint=False)

for y_index, y in enumerate(y_positions):
    for x_index, x in enumerate(x_positions):
        label = z_at_frac_wdm[y_index, x_index]
        text_x = x + jump_x
        text_y = y + jump_y - delta_y
        if y_index==2 and x_index==3:
            ax.text(text_x, text_y, round(label,2), color='darkviolet', ha='center', va='center')
        else:
            ax.text(text_x, text_y, round(label,2), color='black', ha='center', va='center')

            
fig.colorbar(im)
plt.show()

#%% plot S8
fig = plt.figure()
ax = fig.add_subplot(111)

ax.set_xlabel(r'$\mathrm{Log}_{10}(\Gamma / [\mathrm{Gyr}^{-1}])$', fontsize =18)
ax.set_ylabel(r'$\mathrm{Log}_{10}(\varepsilon)$', fontsize = 18)
ax.set_title(r'$\bar{\rho}_{\mathrm{wdm}}/(\bar{\rho}_{\mathrm{wdm}}+\bar{\rho}_{\mathrm{dcdm}}) [\%]  \ \ \mathrm{at} \ \ z=0.32$', fontsize = 22, pad=15)
#ax.set_title(r'$S_8 = \sigma_8 \sqrt{\Omega_m/0.3} $', fontsize = 18, pad=15)

ax.xaxis.set_major_locator(plt.MaxNLocator(len(Log10_Gam_kmsMpc)))
ax.yaxis.set_major_locator(plt.MaxNLocator(len(Log10_eps)))

ax.set_yticks([-1.3,-1.6,-1.9,-2.2,-2.5,-2.8,-3.1])
ax.set_yticklabels(['$-1.3$','$-1.6$','$-1.9$','$-2.2$','$-2.5$','$-2.8$','$-3.1$'])

ax.set_xticks([-3.0,-2.7,-2.4,-2.1,-1.8,-1.5,-1.2,-0.9])
ax.set_xticklabels(['$-3.0$','$-2.7$','$-2.4$','$-2.1$','$-1.8$','$-1.5$','$-1.2$','$-0.9$'])

ax.plot(log_Gam_1,log_eps_1,color='black',linestyle='solid')
ax.plot(log_Gam_2,log_eps_2,color='black',linestyle='solid')


ax.plot(Log10_Gam_kmsMpc-3.0,Log10eps_for_kfs0p7, color= 'red',linestyle='dashed')
ax.plot([-3.3,-3.0],[Log10eps_for_kfs0p7[0],Log10eps_for_kfs0p7[0]], color= 'red',linestyle='dashed')
ax.plot([-0.9,-0.6],[Log10eps_for_kfs0p7[len(Log10_Gam_kmsMpc)-1],Log10eps_for_kfs0p7[len(Log10_Gam_kmsMpc)-1]+0.07], color= 'red',linestyle='dashed')
ax.annotate("", xy=(0.75-3.0,Log10eps_for_kfs0p7[2]-0.15), xytext=(0.75-3.0,Log10eps_for_kfs0p7[2]), arrowprops=dict(arrowstyle="->",linewidth=2,color= 'red'),xycoords='data',textcoords='data',fontsize=15)
ax.text(0.8-3.0,Log10eps_for_kfs0p7[2]-0.1, r'$k_{\mathrm{fs}} > k_{\mathrm{NL}}$',color= 'red',fontsize=14 )

ax.plot(Log10_Gam_kmsMpc-3.0,Log10eps_for_kfs0p2, color= 'darkorange',linestyle='dashed')
ax.plot([-3.3,-3.0],[Log10eps_for_kfs0p2[0],Log10eps_for_kfs0p2[0]], color= 'darkorange',linestyle='dashed')
ax.plot([-0.9,-0.6],[Log10eps_for_kfs0p2[len(Log10_Gam_kmsMpc)-1],Log10eps_for_kfs0p2[len(Log10_Gam_kmsMpc)-1]+0.07], color= 'darkorange',linestyle='dashed')
ax.annotate("", xy=(0.75-3.0,Log10eps_for_kfs0p2[4]-0.15), xytext=(0.75-3.0,Log10eps_for_kfs0p2[4]), arrowprops=dict(arrowstyle="->",linewidth=2,color= 'darkorange'),xycoords='data',textcoords='data',fontsize=15)
ax.text(0.8-3.0,Log10eps_for_kfs0p2[4]-0.1, r'$k_{\mathrm{fs}} > k_{\mathrm{max}}$',color= 'darkorange',fontsize=14 )


#im = ax.imshow(S_8, extent=extent, origin='upper', interpolation='bicubic', cmap='GnBu_r')
im = ax.imshow(frac_wdm_at_z, extent=extent, origin='upper', interpolation='bicubic', cmap='GnBu')


jump_x = (x_end - x_start)/(2.0*len(Log10_Gam_kmsMpc))
jump_y = (y_end - y_start)/(2.0*len(Log10_eps))
x_positions = np.linspace(start=x_start, stop=x_end, num=len(Log10_Gam_kmsMpc), endpoint=False)
y_positions = np.linspace(start=y_end, stop=y_start, num=len(Log10_eps), endpoint=False)

for y_index, y in enumerate(y_positions):
    for x_index, x in enumerate(x_positions):
#        label = S_8[y_index, x_index]
        label = frac_wdm_at_z[y_index, x_index]
        text_x = x + jump_x
        text_y = y + jump_y - delta_y
        if y_index==2 and x_index==3:
            ax.text(text_x, text_y, round(label,2), color='darkviolet', ha='center', va='center')
        else:
            ax.text(text_x, text_y, round(label,2), color='black', ha='center', va='center')

            
fig.colorbar(im)
plt.show()
