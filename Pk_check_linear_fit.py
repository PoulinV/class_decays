# import necessary modules
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import numpy as np
from classy import Class
from scipy.interpolate import interp1d
plt.rcParams["figure.figsize"] = [8.0,6.0]

#%%

kk = np.logspace(-4,1,1000) # k in h/Mpc
zhigh=2. 
Pk1 = [] # P(k) in (Mpc/h)**3
Pk2 = [] # P(k) in (Mpc/h)**3
Pk3 = [] # P(k) in (Mpc/h)**3
Pk4 = []# P(k) in (Mpc/h)**3
Pk5 = []# P(k) in (Mpc/h)**3


tau_dcdm_1 = 100 #units of Gyr
Gamma_dcdm_1 = 1/tau_dcdm_1
Gamma_dcdm_1 =Gamma_dcdm_1/1.02e-3


omega_b =0.02233
omega_cdm=0.1198
omega_m = omega_b + omega_cdm
h_lcdm=0.6737
#set general configuration
common_settings = {'output':'tCl,pCl,lCl, mPk',
                   'lensing':'yes',
                   'l_max_scalars':2600,
                   'n_s':0.9652,
                   'ln10^{10}A_s':3.043,
                   'tau_reio':0.0540,
                   'omega_cdm':omega_cdm,
                   'omega_b': omega_b,
#                   '100*theta_s':1.04089,
                   'h':h_lcdm,
                   'P_k_max_h/Mpc':10.0,
                   'z_max_pk' : 4.0,
#                   'nonlinear': 'halofit'
                   }

#%%
#Planck 2018 best-fit (TT,TE,EE+lowE+lensing)

M = Class()

print("~~~~~computing LCDM ~~~~~")
M.set(common_settings)
#assume three massless neutrinos for the moment
M.compute()

h = M.h() # get reduced Hubble for conversions to 1/Mpc

for k in kk:
    Pk1.append(M.pk(k*h,0.)*h**3) # function .pk(k,z)
    Pk2.append(M.pk(k*h,zhigh)*h**3) # function .pk(k,z)

M.struct_cleanup()
M.empty()

fpk_ref_z0 = interp1d(kk, Pk1)
fpk_ref_z2 = interp1d(kk, Pk2)


#%%

M.set(common_settings)
M.set({
    'omega_cdm': 0.00001,
    'omega_ini_dcdm': omega_cdm,
    'Gamma_dcdm': Gamma_dcdm_1,
    })
M.compute()
h = M.h()
print("~~~~~computing DCDM  ~~~~~")
#t_i = time.time()

for k in kk:
    Pk3.append(M.pk(k*h,0.)*h**3) # function .pk(k,z)
    Pk4.append(M.pk(k*h,zhigh)*h**3) # function .pk(k,z)
    
M.struct_cleanup()
M.empty()

fpk_dcdm1_z0 = interp1d(kk, Pk3)
fpk_dcdm1_z2 = interp1d(kk, Pk4)


#%% compute fitting formula
u = omega_b/0.02216
v = h_lcdm/0.6776
w = omega_m/0.14116

alpha = (5.323-1.4644*u-1.391*v) + (-2.055+1.329*u+0.8672*v)*w + (0.2682-0.3509*u)*w**2
beta = 0.9260 + (0.05735-0.02690*v)*w + (-0.01373+0.006713*v)*w**2
gamma = (9.553-0.7860*v) + (0.4884+0.1754*v)*w + (-0.2512+0.07558*v)*w**2

def eps_lin(tau,z): #tau should be given in units of Gyrs
    return alpha*(1/tau)**beta*(1/(1+0.105*z))**gamma

def a(tau,z): #tau should be given in units of Gyrs
    return -0.18+0.7208+(2.027/tau)+3.031/(1.+1.1*z)

def b(tau,z): #tau should be given in units of Gyrs
    return -0.09+0.0120+(2.786/tau)+0.6699/(1.+1.1*z)

def p(tau,z): #tau should be given in units of Gyrs
    return -0.099+1.045+(1.225/tau)+0.2207/(1.+1.1*z)

def q(tau,z): #tau should be given in units of Gyrs
    return -0.056+0.9922+(1.735/tau)+0.2154/(1.+1.1*z)

def eps_nonlin(k,tau,z):  #taking f_dcdm = 1
    #tau should be given in units of Gyrs and k in units of Mpc^-1
    return eps_lin(tau,z)*(1.0+a(tau,z)*k**p(tau,z))/(1.0+b(tau,z)*k**q(tau,z))
#%%

kk = np.logspace(-2.3,1,1000) # k in h/Mpc

fig, axe_1 = plt.subplots()

#axe_1 = plt.subplot(211)
#axe_2 = plt.subplot(212, sharex = ax_1)
#plt.subplots_adjust(hspace=0)

axe_1.set_xlim([kk[0],kk[-1]])

axe_1.set_ylim([0.8,1.05])

axe_1.set_xlabel(r'$k \,\,\,\, [h/\mathrm{Mpc}]$', fontsize=20)
axe_1.set_ylabel(r'$P_{\mathrm{DDM}}(k)/P_{\Lambda\mathrm{CDM}}(k)$',fontsize=20)


plt.title(r'$\tau = %.0f \ \mathrm{Gyrs}$'%tau_dcdm_1,fontsize=20)

axe_1.semilogx(kk,fpk_dcdm1_z0(kk)/fpk_ref_z0(kk),color='purple',linestyle='solid',label=r'$z = 0$')
axe_1.semilogx(kk,fpk_dcdm1_z2(kk)/fpk_ref_z2(kk),color='blue',linestyle='solid', label=r'$z = 2$')

axe_1.axhline(1.-eps_lin(tau_dcdm_1,0), color ='purple',linestyle='dotted')
axe_1.axhline(1.-eps_lin(tau_dcdm_1,zhigh), color='blue', linestyle='dotted')

axe_1.semilogx(kk,1.-eps_nonlin(kk,tau_dcdm_1,0),color='purple',linestyle='dashed')
axe_1.semilogx(kk,1.-eps_nonlin(kk,tau_dcdm_1,zhigh),color='blue',linestyle='dashed')

black_line1 = mlines.Line2D([], [], color='black', linestyle='solid', label=r'CLASS linear')
black_line2 = mlines.Line2D([], [], color='black', linestyle='dotted', label=r'$1-\epsilon_{\mathrm{lin}}$')
black_line3 = mlines.Line2D([], [], color='black', linestyle='dashed', label=r'$1-\epsilon_{\mathrm{nonlin}}$')

legend2 = plt.legend(handles=[black_line1,black_line2, black_line3], loc='lower right', fontsize=18, frameon=False)

axe_1.add_artist(legend2)

axe_1.legend(frameon=False, loc='lower left', fontsize=18, borderaxespad=0.)


axe_1.tick_params(axis="y", labelsize=18)
axe_1.tick_params(axis="x", labelsize=18)


plt.show()
plt.clf()


