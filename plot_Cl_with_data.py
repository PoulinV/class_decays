
# coding: utf-8

# In[ ]:

# import necessary modules
#get_ipython().magic(u'matplotlib inline')
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


###planck 2015
# l_TT_low,Dl_TT_low,err_TT_low= np.loadtxt("error_Planck/COM_PowerSpect_CMB-TT-loL-full_R2.02.txt",unpack=True,usecols=(0,1,2))
# l_TE_low,Dl_TE_low,err_TE_low= np.loadtxt("error_Planck/COM_PowerSpect_CMB-TE-loL-full_R2.02.txt",unpack=True,usecols=(0,1,2))
# l_TT_high,Dl_TT_high,err_TT_high= np.loadtxt("error_Planck/COM_PowerSpect_CMB-TT-hiL-binned_R2.02.txt",unpack=True,usecols=(0,3,4))
# l_TE_high,Dl_TE_high,err_TE_high= np.loadtxt("error_Planck/COM_PowerSpect_CMB-TE-hiL-binned_R2.02.txt",unpack=True,usecols=(0,3,4))
# l_EE_low,Dl_EE_low,err_EE_low= np.loadtxt("error_Planck/COM_PowerSpect_CMB-EE-loL-full_R2.02.txt",unpack=True,usecols=(0,1,2))
# l_EE_high,Dl_EE_high,err_EE_high= np.loadtxt("error_Planck/COM_PowerSpect_CMB-EE-hiL-binned_R2.02.txt",unpack=True,usecols=(0,3,4))
# lmin_phiphi,lmax_phiphi,cl_phiphi,err_phiphi= np.loadtxt("error_Planck/agressive_lensing.csv",unpack=True)

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
ax_1.set_ylim([-0.1,0.1])
ax_2.set_ylim([-0.1,0.1])
# ax_1.set_ylim([30,30000])
# ax_2.set_ylim([-1,9])
# ax_2.set_ylim([-2,2])



##LCDM BESTFIT Planck 2018####
common_settings = {'output':'tCl,pCl,lCl,mPk',
                   'lensing':'yes',
                   'format':'camb',
                   'l_max_scalars':2600
                   }
                   # ,'input_verbose': 1,
                   # 'background_verbose': 1,
                   # 'thermodynamics_verbose': 1,
                   # 'perturbations_verbose': 1,
                   # 'transfer_verbose': 1,
                   # 'primordial_verbose': 1,
                   # 'spectra_verbose': 1,
                   # 'nonlinear_verbose': 1,
                   # 'lensing_verbose': 1,
                   # 'output_verbose': 1}
M = Class()




###choose the value of Gamma and the number of bins in perts
Gamma_dcdm = 10
nbins = 300
m_dcdm = 0.00001
#
print("~~~~~computing reference~~~~~")
M.set(common_settings)
M.set({
'omega_cdm': 0.00001,
'Omega_ini_dcdm': 0.24,
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
    # for i in range(0,int(width)):
    #     weight_total += (l_ini+i)*(l_ini+i+1)
    #     result += np.sqrt(2/(2*(l_ini+float(i))+1))*(l_ini+float(i))*(l_ini+float(i)+1)
    # print l_ini,l_ini+width,result, weight_total
    # return result/weight_total/np.sqrt(2)/np.pi
    for i in range(0,int(width)):
        weight_total += (l_ini+i)*(l_ini+i+1)
        result += 2/(2*(l_ini+float(i))+1)*(l_ini+float(i))*(l_ini+float(i)+1)*(l_ini+float(i))*(l_ini+float(i)+1)*fTT_ref(l_ini+i)*fTT_ref(l_ini+i)
        Clb += (l_ini+float(i))*(l_ini+float(i)+1)*fTT_ref(l_ini+i)
    # print l_ini,l_ini+width,np.sqrt(result), Clb, np.sqrt(result)/Clb
    return np.sqrt(result)/Clb

#
#
#
#
l_min = 2.;
l_max = 2000;
n_step = 100.;
j=0.
step = l_min
width= 25
# while j < n_step:
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


        ax_1.add_patch(
            patches.Rectangle(
                (int(step), -1*binned_cosmic_variance(result,int(step),width)),   # (x,y)
                width,          # width
                2*binned_cosmic_variance(result,int(step),width),          # height
                color='g',
                alpha=0.1
            )
        )

        ax_2.add_patch(
            patches.Rectangle(
                (int(step), -1*binned_cosmic_variance(result,int(step),width)),   # (x,y)
                width,          # width
                2*binned_cosmic_variance(result,int(step),width),          # height
                color='g',
                alpha=0.1
            )
        )



print("~~~~~print planck error bar around mean, i.e, residuals are 0 for simplicity~~~~~")


ax_1.errorbar(lTT, DlTT_mean/DlTT_mean-1, yerr=(DlTT_error_plus)/DlTT_mean, fmt='.',color='r')
#ax_1.errorbar(lTT, DlTT_mean/(DlTT_bestfit)-1, yerr=[(DlTT_error_minus)/DlTT_bestfit, (+DlTT_error_plus)/DlTT_bestfit], fmt='.',color='b',label=r'TT')
ax_2.errorbar(lEE, DlEE_mean/DlEE_mean-1, yerr=DlEE_error_plus/DlEE_mean, fmt='.',color='r')
# ax_2.errorbar(lEE, DlEE_mean/(DlEE_bestfit)-1, yerr=[(-DlEE_error_minus)/DlEE_bestfit, (+DlEE_error_plus)/DlEE_bestfit], fmt='.',color='b',label=r'EE')

timeafterref=time.time()

print("~~~~~time =%.f s; computing our code~~~~~"%(timeafterref-start_time))
M.set(common_settings)
M.set({
'omega_cdm': 0.00001,
'Omega_ini_dcdm2': 0.24,
'Gamma_dcdm': Gamma_dcdm,
'M_dcdm': 1,
'm_dcdm': m_dcdm,
'background_ncdm_distribution': 1,
'Quadrature strategy': 4,
'N_ncdm': 1,
'evolver': 0,
'ncdm_fluid_approximation': 2,
'ncdm_fluid_trigger_tau_over_tau_k': 25,
'Number of momentum bins perturbs': nbins,
'massive_daughter_perturbations': 'yes',
'dark_radiation_perturbations': 'yes'
})

#want verbose? uncomment the following

M.compute()
print("~~~~~done computing our code in %.f s~~~~~"%(time.time()-timeafterref))
print("~~~~~ready to plot~~~~~")

clM = M.lensed_cl(2500)
ll_LCDM = clM['ell'][2:]
clTT_LCDM = clM['tt'][2:]
clEE_LCDM = clM['ee'][2:]
fTT_ourcode = interp1d(ll_LCDM,clTT_LCDM*(ll_LCDM)*(ll_LCDM+1)/2/np.pi*(T_cmb*1.e6)**2)
fEE_ourcode = interp1d(ll_LCDM,clEE_LCDM*(ll_LCDM)*(ll_LCDM+1)/2/np.pi*(T_cmb*1.e6)**2)


# plt.loglog(ll_LCDM,clTT_LCDM*(ll_LCDM)*(ll_LCDM+1)/2/np.pi*(T_cmb*1.e6)**2,lw=2,label=r'$H_0 = %.1f $km/s/Mpc'%(H0))
# fcdm_100 = fcdm*100
# ax_2.semilogx(ll_LCDM,clTT_LCDM/fTT(ll_LCDM)-1,lw=2,label=r'$f_{\rm cdm} = %.2f$'%(fcdm))
ax_1.semilogx(ll_LCDM,fTT_ourcode(ll_LCDM)/fTT_ref(ll_LCDM)-1,lw=2)
ax_2.semilogx(ll_LCDM,fEE_ourcode(ll_LCDM)/fEE_ref(ll_LCDM)-1,lw=2)
ax_2.set_xlabel(r'$\ell$',fontsize=20)
ax_1.set_ylabel(r'$\Delta D_\ell^\mathrm{TT}/D_\ell^\mathrm{TT}$',fontsize=20)
ax_2.set_ylabel(r'$\Delta D_\ell^\mathrm{EE}/D_\ell^\mathrm{EE}$',fontsize=20)
# ax_1.tick_params(fontsize=20)
# ax_2.tick_params(fontsize=20)
# ax_1.tick_params(fontsize=20)
ax_2.tick_params(axis="x", labelsize=20)
ax_2.tick_params(axis="y", labelsize=20)
ax_1.tick_params(axis="y", labelsize=20)
# plt.set_logscalex()
plt.legend(frameon=False,prop={'size':15},loc='lower right',borderaxespad=0.)
# plt.savefig('omcdm_%.1f_scan.pdf'%(fcdm_100), bbox_inches='tight')
plt.show()
# plt.savefig('omcdm_%.2f_scan_v2.pdf'%(omega_cdm), bbox_inches='tight')
plt.clf()

M.struct_cleanup()
#

# In[ ]:
