
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
# lTT,DlTT_mean,DlTT_error_minus,DlTT_error_plus,DlTT_bestfit= np.loadtxt("error_Planck/Planck2018_errorTT.txt",unpack=True)
# lEE,DlEE_mean,DlEE_error_minus,DlEE_error_plus,DlEE_bestfit= np.loadtxt("error_Planck/Planck2018_errorEE.txt",unpack=True)
# lTE,DlTE_mean,DlTE_error_minus,DlTE_error_plus,DlTE_bestfit= np.loadtxt("error_Planck/Planck2018_errorTE.txt",unpack=True)
lTT,lminTT,lmaxTT,DlTT,DlTT_error= np.loadtxt("error_Planck/Planck2015_errorTT.txt",unpack=True)
lEE,lminEE,lmaxEE,DlEE,DlEE_error= np.loadtxt("error_Planck/Planck2015_errorEE.txt",unpack=True)
lTTlowl,DlTTlowl,DlTT_errorup_lowl,DlTT_errordown_lowl= np.loadtxt("error_Planck/Planck2015_errorTT_lowl.txt",unpack=True)
lEElowl,DlEElowl,DlEE_errorup_lowl,DlEE_errordown_lowl= np.loadtxt("error_Planck/Planck2015_errorEE_lowl.txt",unpack=True)
#ref_M0p999= np.loadtxt("output/dcdmdr_G10000_M0p999_100000bins_cl_lensed.dat")
# 1:l     2:TT                     3:EE                     4:TE                     5:BB                     6:phiphi                 7:TPhi                   8:Ephi
# np.stack(lTTlowl,lTT)
lTT_all = np.append(lTTlowl,lTT)
lEE_all = np.append(lEElowl,lEE)
fract_TT_down = np.append((DlTT_errordown_lowl)/DlTTlowl,(DlTT_error)/DlTT)
fract_TT_up = np.append((DlTT_errorup_lowl)/DlTTlowl,(DlTT_error)/DlTT)
fract_EE_down = np.append((DlEE_errordown_lowl)/DlEElowl,(DlEE_error)/DlEE)
fract_EE_up = np.append((DlEE_errorup_lowl)/DlEElowl,(DlEE_error)/DlEE)
print(fract_TT_down,fract_TT_up)
#ax_2.errorbar(lEE, DlEE_mean/DlEE_mean-1, yerr=DlEE_error_plus/DlEE_mean, fmt='.',color='r')


print(lTT_all,lEE_all)

##create plot
ax_1 = plt.subplot(311)
ax_2 = plt.subplot(312, sharex = ax_1)
ax_3 = plt.subplot(313)
#plt.subplots_adjust(hspace=0)
#ax_1.set_ylim([-0.07,0.07])
#ax_2.set_ylim([-0.12,0.12])
ax_1.set_ylim([-0.13,0.13])
ax_2.set_ylim([-0.13,0.13])
ax_3.set_ylim([-0.005,0.005])
ax_1.set_xlim([2,2500])
ax_2.set_xlim([2,2500])
ax_3.set_xlim([0.0001,1])
# ax_1.set_ylim([30,30000])
# ax_2.set_ylim([-1,9])
# ax_2.set_ylim([-2,2])

# PUT PLANCK 2015 TT, TE, EE+lowP+lensing best-fit values
#for omega_b, n_s, A_s, tau_reio, Omega_ini_dcdm2 and H0

##LCDM BESTFIT Planck 2018####
common_settings = {'output':'tCl,pCl,lCl,mPk',
                   'P_k_max_1/Mpc':1.5,
                   'lensing':'yes',
                   'format':'camb',
                   'l_max_scalars':2600,
                   'n_s':0.9652,
                   'ln10^{10}A_s':3.043,
                   'tau_reio':0.0540,
                   'omega_b':0.02233,
                   'h':0.6737,
                   'deg_ncdm':'3',
                   'N_ur':'0.00641', ###this is required to achieve Nef':'.046 in the standard case
                   'background_ncdm_distribution':'2',# 0 fermi dirac, 1 massive daughters, 2 decaying neutrinos
                   'N_ncdm':'1',
                   'dark_radiation_perturbations':'yes',
                   'lensing':'yes',
                   'loop_over_background':'yes',
                   'decaying_neutrino_model':'yes'
                   }

M = Class()




###choose the value of Gamma and the number of bins in perts
# Gamma_nu = [2e4,1e5,2e5,5e5,8e5]
# Gamma_nu = [10**3.7,10**3.2,10**2.1,10**1.7,10**1.5]
Gamma_nu = [10**2.7,10**2.2,10**1.8,10**1.5,10**1.3]
# Gamma_nu = [1e3,1e4,9e4,4e5]
m_ncdm = [0.1/3,0.15/3,0.2/3,0.25/3,0.3/3]
# m_ncdm = [0.06/3,0.12/3,0.18/3,0.24/3,0.3/3]
# m_ncdm = [0.3/3,0.5/3,0.7/3,0.9/3]
# Gamma_nu = [10**3.7,10**3]
# m_ncdm = [0.16/3,0.25/3]
#m_ncdm = 0.3
#tau of 0.1 Gyrs
i = 0
colors = ['blue','red','green','cyan','orange','black']
# for m_ncdm in [0.1/3,0.3/3,0.5/3,0.7/3,0.9/3]:
for i in [0,1,2,3,4]:
# for m_ncdm in [0.9/3]:

        print("~~~~~computing reference for m = %.2f eV & Gamma = %.2f km/s/sMpc~~~~~"%(m_ncdm[i],Gamma_nu[i]))
        M.set(common_settings)
        M.set({
        # 'Gamma_neutrinos':'%.5f,%.5f,%.5f'%(Gamma_nu,Gamma_nu,Gamma_nu),
        'Gamma_neutrinos':'%.5f'%(Gamma_nu[i]),
        # m_ncdm':'0.08333
        # m_ncdm':'0.00001
        # m_ncdm':'0.01572
        # m_ncdm':' 0.0501811, 0.0509229, 0.00889605
        # m_ncdm':' 0.1501811, 0.1509229, 0.00889605 #if neutrino hierarchy is set to normal or inverted, the last 2 entries will be overwritten
        # m_ncdm':' 0.0785, 0.1509229, 0.00889605 #if neutrino hierarchy is set to normal or inverted, the last 2 entries will be overwritten
        'm_ncdm':'%.5f'%(m_ncdm[i]),#if neutrino hierarchy is set to normal or inverted, the last 2 entries will be overwritten
        'neutrino_hierarchy':'normal',
        'include_new_term_decay_neutrinos':1
        })
        M.compute()
        clM = M.lensed_cl(2600)
        ll_LCDM = clM['ell'][2:]
        clTT_LCDM = clM['tt'][2:]
        clEE_LCDM = clM['ee'][2:]
        kvec = np.logspace(-4,0,100)
        T_cmb = 2.7225 #we change units for Planck
        fTT_ref = interp1d(ll_LCDM,clTT_LCDM*(ll_LCDM)*(ll_LCDM+1)/2/np.pi*(T_cmb*1.e6)**2)
        fEE_ref = interp1d(ll_LCDM,clEE_LCDM*(ll_LCDM)*(ll_LCDM+1)/2/np.pi*(T_cmb*1.e6)**2)
        pk_ref = []
        for k in kvec:
            pk_ref.append(M.pk(k,0.))
        #fTT_ref = interp1d(ref_M0p999[:,0],ref_M0p999[:,1]*(T_cmb*1.e6)**2)
        #fEE_ref = interp1d(ref_M0p999[:,0],ref_M0p999[:,2]*(T_cmb*1.e6)**2)


        M.struct_cleanup()
        M.empty()
        timeafterref=time.time()


        print("~~~~~time =%.f s; computing our code~~~~~"%(timeafterref-start_time))
        M.set(common_settings)
        M.set({
        'Gamma_neutrinos':'%.5f'%(Gamma_nu[i]),
        # m_ncdm':'0.08333
        # m_ncdm':'0.00001
        # m_ncdm':'0.01572
        # m_ncdm':' 0.0501811, 0.0509229, 0.00889605
        # m_ncdm':' 0.1501811, 0.1509229, 0.00889605 #if neutrino hierarchy is set to normal or inverted, the last 2 entries will be overwritten
        # m_ncdm':' 0.0785, 0.1509229, 0.00889605 #if neutrino hierarchy is set to normal or inverted, the last 2 entries will be overwritten
        'm_ncdm':'%.5f'%(m_ncdm[i]),#if neutrino hierarchy is set to normal or inverted, the last 2 entries will be overwritten
        'neutrino_hierarchy':'degenerate',
        'include_new_term_decay_neutrinos':0
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
        pk_ourcode = []
        for k in kvec:
            pk_ourcode.append(M.pk(k,0.))
        M.struct_cleanup()
        M.empty()




        ax_1.semilogx(ll_LCDM,fTT_ourcode(ll_LCDM)/fTT_ref(ll_LCDM)-1,colors[i], label=r'$\Gamma = %.1f, \sum M_\nu = %.2f$'%(Gamma_nu[i],3*m_ncdm[i]))
        ax_2.semilogx(ll_LCDM,fEE_ourcode(ll_LCDM)/fEE_ref(ll_LCDM)-1,colors[i])
        ax_3.semilogx(kvec,np.array(pk_ourcode)/np.array(pk_ref)-1,colors[i])







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
        width= 30
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

        # ax_1.fill_between(l_cosmic_variance_1, -0.13,0.13, color='lightgray' )
        # ax_1.fill_between(l_cosmic_variance_2, -slope,slope, color='lightgray' )
        # ax_1.fill_between(lTT, -(DlTT_error)/DlTT, +(DlTT_error)/DlTT, color='lightgray')
        # ax_1.fill_between(lTTlowl, (DlTT_errordown_lowl-DlTTlowl)/DlTTlowl, +(DlTT_errorup_lowl-DlTTlowl)/DlTTlowl, color='lightgray')
        # ax_1.fill_between(lTT_all, fract_TT_down, fract_TT_up, color='lightgray')
        # ax_2.fill_between(lEE_all, fract_EE_down, fract_EE_up, color='lightgray')
        ax_1.errorbar(lTT_all,0*lTT_all,yerr=[fract_TT_down, fract_TT_up], color='lightgray',fmt='.')
        ax_2.errorbar(lEE_all,0*lEE_all,yerr=[fract_EE_down, fract_EE_up], color='lightgray',fmt='.')

        # print(lTTlowl, (DlTT_errordown_lowl-DlTTlowl)/DlTTlowl, (DlTT_errorup_lowl-DlTTlowl)/DlTTlowl)
        #ax_2.errorbar(lEE, DlEE_mean/DlEE_mean-1, yerr=DlEE_error_plus/DlEE_mean, fmt='.',color='r')

        # ax_2.fill_between(l_cosmic_variance, -0.13,0.13, color='lightgray' )
        # ax_2.fill_between(lEE, -(DlEE_error)/DlEE, +(DlEE_error)/DlEE, color = 'lightgray')
        # ax_1.fill_between(lEElowl, (DlEElowl-DlEE_errordown_lowl)/DlEElowl, +(DlEE_errorup_lowl-DlEElowl)/DlEElowl, color='lightgray')

        kvec2 = np.logspace(-2,0,100)
        ax_3.fill_between(kvec2,-0.002*kvec2/kvec2,0.002*kvec2/kvec2,color='lightgray')





print("~~~~~ready to plot~~~~~")


ax_1.set_xlabel(r'$\mathrm{multipole} \, \, \ell$',fontsize=15)
ax_3.set_xlabel(r'$k$',fontsize=15)
ax_1.set_ylabel(r'$\frac{C_\ell^\mathrm{TT}(\mathrm{old})}{C_\ell^\mathrm{TT}(\mathrm{new} )} -1$',fontsize=15)
ax_2.set_ylabel(r'$\frac{C_\ell^\mathrm{EE}(\mathrm{old})}{C_\ell^\mathrm{EE}(\mathrm{new} )} -1$',fontsize=15)
ax_3.set_ylabel(r'$\Delta P(k)/P(k)$',fontsize=15)


ax_3.tick_params(axis="x", labelsize=15)
ax_2.tick_params(axis="x", labelsize=15)
ax_1.tick_params(axis="x", labelsize=15)
ax_2.tick_params(axis="y", labelsize=15)
ax_1.tick_params(axis="y", labelsize=15)
ax_3.tick_params(axis="y", labelsize=15)

ax_1.legend(frameon=False,fontsize = 12,loc='upper right',borderaxespad=0.)

plt.show()
plt.clf()



#

# In[ ]:
