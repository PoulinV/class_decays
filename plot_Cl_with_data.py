
# coding: utf-8

# In[ ]:

# import necessary modules
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from classy import Class
from scipy.optimize import fsolve
import math
from matplotlib import rc
import matplotlib.patches as patches

from scipy.interpolate import interp1d
from matplotlib.ticker import FixedLocator
from math import floor
from mpl_toolkits.axes_grid1 import make_axes_locatable

import time
start_time = time.time()

# In[ ]:

# esthetic definitions for the plots

#rc('font',**{'family':'serif','serif':['Times']})
#rc('text', usetex=True)
# matplotlib.rc('font', **font)
#matplotlib.mathtext.rcParams['legend.fontsize']='medium'
#plt.rcParams["figure.figsize"] = [8.0,6.0]


###planck 2018: problem, only the large l are binned. lowell unbinned?
lTT,DlTT_mean,DlTT_error_minus,DlTT_error_plus,DlTT_bestfit= np.loadtxt("error_Planck/Planck2018_errorTT.txt",unpack=True)
lEE,DlEE_mean,DlEE_error_minus,DlEE_error_plus,DlEE_bestfit= np.loadtxt("error_Planck/Planck2018_errorEE.txt",unpack=True)
lTE,DlTE_mean,DlTE_error_minus,DlTE_error_plus,DlTE_bestfit= np.loadtxt("error_Planck/Planck2018_errorTE.txt",unpack=True)
#ref_M0p999= np.loadtxt("output/dcdmdr_G10000_M0p999_100000bins_cl_lensed.dat")
# 1:l     2:TT                     3:EE                     4:TE                     5:BB                     6:phiphi                 7:TPhi                   8:Ephi

##create plot
ax_1 = plt.subplot(211)
ax_2 = plt.subplot(212, sharex = ax_1)
plt.subplots_adjust(hspace=0)
#ax_1.set_ylim([-0.07,0.07])
#ax_2.set_ylim([-0.12,0.12])
ax_1.set_ylim([-0.5,0.5])
ax_2.set_ylim([-0.5,0.5])
ax_1.set_xlim([2,2500])
ax_2.set_xlim([2,2500])
# ax_1.set_ylim([30,30000])
# ax_2.set_ylim([-1,9])
# ax_2.set_ylim([-2,2])

# PUT PLANCK 2015 TT, TE, EE+lowP+lensing best-fit values 
#for omega_b, n_s, A_s, tau_reio, Omega_ini_dcdm2 and H0

##LCDM BESTFIT Planck 2018####
common_settings = {'output':'tCl,pCl,lCl,mPk',
                   'lensing':'yes',
                   'format':'camb',
                   'l_max_scalars':2600,                    
                   'n_s':0.9652,
                   'ln10^{10}A_s':3.043,
                   'tau_reio':0.0540,
                   'omega_b':0.02233,
                   'h':0.6737,
                   }

M = Class()




###choose the value of Gamma and the number of bins in perts
Gamma_dcdm = 9803.9215
tau =1./(Gamma_dcdm*1.02e-3)
#tau of 0.1 Gyrs

nbins = 300
m_dcdm = 0.00001
#epsilon of almost 0.5

#
print("~~~~~computing reference~~~~~")
M.set(common_settings)
M.set({
'omega_cdm': 0.00001,
'omega_ini_dcdm': 0.1198,
'N_ncdm':1,
'background_ncdm_distribution': 0,
'N_ur':2.0328,
'm_ncdm':0.06,
'Gamma_dcdm': Gamma_dcdm,
'evolver': 0,
'dark_radiation_perturbations': 'yes'
})
M.compute()
clM = M.lensed_cl(2600)
ll_LCDM = clM['ell'][2:]
clTT_LCDM = clM['tt'][2:]
clEE_LCDM = clM['ee'][2:]
T_cmb = 2.7225 #we change units for Planck
fTT_ref = interp1d(ll_LCDM,clTT_LCDM*(ll_LCDM)*(ll_LCDM+1)/2/np.pi*(T_cmb*1.e6)**2)
fEE_ref = interp1d(ll_LCDM,clEE_LCDM*(ll_LCDM)*(ll_LCDM+1)/2/np.pi*(T_cmb*1.e6)**2)
#fTT_ref = interp1d(ref_M0p999[:,0],ref_M0p999[:,1]*(T_cmb*1.e6)**2)
#fEE_ref = interp1d(ref_M0p999[:,0],ref_M0p999[:,2]*(T_cmb*1.e6)**2)


M.struct_cleanup()
M.empty()


print("~~~~~compute binned cosmic variance based on ref~~~~~")

def binned_cosmic_variance (result,l_ini,width):
    central_l = l_ini+width/2
    weight_total = 0
    result = 0
    Clb = 0
    for i in range(0,int(width)):
        weight_total += (l_ini+i)*(l_ini+i+1)
        result += 2/(2*(l_ini+float(i))+1)*(l_ini+float(i))*(l_ini+float(i)+1)*(l_ini+float(i))*(l_ini+float(i)+1)*fTT_ref(l_ini+i)*fTT_ref(l_ini+i)
        Clb += (l_ini+float(i))*(l_ini+float(i)+1)*fTT_ref(l_ini+i)
    return np.sqrt(result)/Clb

#
l_min = 2.;
l_max = 2000;
n_step = 100.;
j=0.
step = l_min
width= 25
while step < l_max:
        result = 0.0
        if step < 29:
            width = 1
            step = l_min+j*width
            j+=1
            if step == 29:
                j = 0
                l_min = 30
        else:
            width = 30
            step = l_min+j*width
            j+=1

        ax_1.add_patch(patches.Rectangle((int(step), -1*binned_cosmic_variance(result,int(step),width)), width, 2*binned_cosmic_variance(result,int(step),width),color='g', alpha=0.1))
        ax_2.add_patch(patches.Rectangle((int(step), -1*binned_cosmic_variance(result,int(step),width)), width, 2*binned_cosmic_variance(result,int(step),width),color='g',alpha=0.1))


print("~~~~~print planck error bar around mean, i.e, residuals are 0 for simplicity~~~~~")


#ax_1.errorbar(lTT, DlTT_mean/DlTT_mean-1, yerr=(DlTT_error_plus)/DlTT_mean, fmt='.',color='r')

l_cosmic_variance = np.linspace(0,48,1000)

l_cosmic_variance_1 = np.linspace(0,30,1000)

l_cosmic_variance_2 = np.linspace(30,48,2)
slope =np.array([0.13,0.0343])

ax_1.fill_between(l_cosmic_variance_1, -0.13,0.13, color='lightgray' )
ax_1.fill_between(l_cosmic_variance_2, -slope,slope, color='lightgray' )
ax_1.fill_between(lTT, -(DlTT_error_plus)/DlTT_mean, +(DlTT_error_plus)/DlTT_mean, color='lightgray')

#ax_2.errorbar(lEE, DlEE_mean/DlEE_mean-1, yerr=DlEE_error_plus/DlEE_mean, fmt='.',color='r')

ax_2.fill_between(l_cosmic_variance, -0.13,0.13, color='lightgray' )
ax_2.fill_between(lEE, -(DlEE_error_plus)/DlEE_mean, +(DlEE_error_plus)/DlEE_mean, color = 'lightgray')


timeafterref=time.time()

print("~~~~~time =%.f s; computing our code~~~~~"%(timeafterref-start_time))
M.set(common_settings)
M.set({
'omega_cdm': 0.00001,
'Omega_ini_dcdm2': 0.2639,
'Gamma_dcdm': Gamma_dcdm,
'M_dcdm': 1,
'm_dcdm': m_dcdm,
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
print("~~~~~ computing our code in %.f s~~~~~"%(time.time()-timeafterref))
t2 = time.time()


clM = M.lensed_cl(2500)
ll_LCDM = clM['ell'][2:]
clTT_LCDM = clM['tt'][2:]
clEE_LCDM = clM['ee'][2:]
fTT_ourcode = interp1d(ll_LCDM,clTT_LCDM*(ll_LCDM)*(ll_LCDM+1)/2/np.pi*(T_cmb*1.e6)**2)
fEE_ourcode = interp1d(ll_LCDM,clEE_LCDM*(ll_LCDM)*(ll_LCDM+1)/2/np.pi*(T_cmb*1.e6)**2)

M.struct_cleanup()
M.empty()

#compute evolution but setting gamma=0 in the wdm fluid eqs.
M.set(common_settings)

M.set({
'omega_cdm': 0.00001,
'Omega_ini_dcdm2': 0.2636,
'Gamma_dcdm': Gamma_dcdm,
'M_dcdm': 1,
'm_dcdm': m_dcdm,
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
'switch_off_gamma_in_wdm_perts': 'yes'
})

M.compute()

print("~~~~~ computing our code in %.f s~~~~~"%(time.time()-t2))
clM = M.lensed_cl(2500)
ll_LCDM = clM['ell'][2:]
clTT_LCDM = clM['tt'][2:]
clEE_LCDM = clM['ee'][2:]
fTT_ourcode2 = interp1d(ll_LCDM,clTT_LCDM*(ll_LCDM)*(ll_LCDM+1)/2/np.pi*(T_cmb*1.e6)**2)
fEE_ourcode2 = interp1d(ll_LCDM,clEE_LCDM*(ll_LCDM)*(ll_LCDM+1)/2/np.pi*(T_cmb*1.e6)**2)

M.struct_cleanup()
M.empty()

print("~~~~~ready to plot~~~~~")


ax_1.semilogx(ll_LCDM,fTT_ourcode(ll_LCDM)/fTT_ref(ll_LCDM)-1,'g', label=r'Fluid WDM ')
ax_2.semilogx(ll_LCDM,fEE_ourcode(ll_LCDM)/fEE_ref(ll_LCDM)-1,'g', label = r'Fluid WDM')
ax_1.semilogx(ll_LCDM,fTT_ourcode2(ll_LCDM)/fTT_ref(ll_LCDM)-1,'r--', label=r'Fluid WDM with $\Gamma = 0$')
ax_2.semilogx(ll_LCDM,fEE_ourcode2(ll_LCDM)/fEE_ref(ll_LCDM)-1,'r--', label = r'Fluid WDM $\Gamma = 0$')
ax_2.set_xlabel(r'$\mathrm{multipole} \, \, \ell$',fontsize=15)
ax_1.set_ylabel(r'$\frac{C_\ell^\mathrm{TT}(\mathrm{approx})}{C_\ell^\mathrm{TT}(\mathrm{full} )} -1$',fontsize=20)
ax_2.set_ylabel(r'$\frac{C_\ell^\mathrm{EE}(\mathrm{approx})}{C_\ell^\mathrm{EE}(\mathrm{full} )} -1$',fontsize=20)


ax_2.tick_params(axis="x", labelsize=18)
ax_2.tick_params(axis="y", labelsize=18)
ax_1.tick_params(axis="y", labelsize=18)

ax_1.legend(frameon=False,fontsize = 15,loc='upper left',borderaxespad=0.)

plt.show()
plt.clf()



#

# In[ ]:
