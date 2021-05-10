# import necessary modules
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import numpy as np
from scipy.interpolate import interp1d

# read Planck errors
lTT,DlTT_mean,DlTT_error_minus,DlTT_error_plus,DlTT_bestfit= np.loadtxt("error_Planck/Planck2018_errorTT.txt",unpack=True)
lEE,DlEE_mean,DlEE_error_minus,DlEE_error_plus,DlEE_bestfit= np.loadtxt("error_Planck/Planck2018_errorEE.txt",unpack=True)
lTE,DlTE_mean,DlTE_error_minus,DlTE_error_plus,DlTE_bestfit= np.loadtxt("error_Planck/Planck2018_errorTE.txt",unpack=True)


#%% READ REFERENCE DATA FILES

files1 = ['/Users/gfranco/class_majoron/output/sup4_2/dcdm_eps_hot_gam_0p1H0_nofluid_cl_lensed.dat']
data1 = []
for data_file1 in files1:
    data1.append(np.loadtxt(data_file1))

files2 = ['/Users/gfranco/class_majoron/output/sup4_2/dcdm_eps_hot_gam_H0_nofluid_cl_lensed.dat']
data2 = []
for data_file2 in files2:
    data2.append(np.loadtxt(data_file2))

files3 = ['/Users/gfranco/class_majoron/output/sup4_2/dcdm_eps_hot_gam_10H0_nofluid_cl_lensed.dat']
data3 = []
for data_file3 in files3:
    data3.append(np.loadtxt(data_file3))

files4 = ['/Users/gfranco/class_majoron/output/sup4_2/dcdm_eps_warm_gam_0p1H0_nofluid_cl_lensed.dat']
data4 = []
for data_file4 in files4:
    data4.append(np.loadtxt(data_file4))

files5 = ['/Users/gfranco/class_majoron/output/sup4_2/dcdm_eps_warm_gam_H0_nofluid_cl_lensed.dat']
data5 = []
for data_file5 in files5:
    data5.append(np.loadtxt(data_file5))
    
files6 = ['/Users/gfranco/class_majoron/output/sup4_2/dcdm_eps_warm_gam_10H0_nofluid_cl_lensed.dat']
data6 = []
for data_file6 in files6:
    data6.append(np.loadtxt(data_file6))
    
files7 = ['/Users/gfranco/class_majoron/output/sup4_2/dcdm_eps_cold_gam_0p1H0_nofluid_cl_lensed.dat']
data7 = []
for data_file7 in files7:
    data7.append(np.loadtxt(data_file7))  
    
files8 = ['/Users/gfranco/class_majoron/output/sup4_2/dcdm_eps_cold_gam_H0_nofluid_cl_lensed.dat']
data8 = []
for data_file8 in files8:
    data8.append(np.loadtxt(data_file8))

files9 = ['/Users/gfranco/class_majoron/output/sup4_2/dcdm_eps_cold_gam_10H0_nofluid_cl_lensed.dat']
data9 = []
for data_file9 in files9:
    data9.append(np.loadtxt(data_file9))
    
files10 = ['/Users/gfranco/class_majoron/output/sup4_2/dcdm_eps_warm2_gam_0p1H0_nofluid_cl_lensed.dat']
data10 = []
for data_file10 in files10:
    data10.append(np.loadtxt(data_file10))

files11 = ['/Users/gfranco/class_majoron/output/sup4_2/dcdm_eps_warm2_gam_H0_nofluid_cl_lensed.dat']
data11 = []
for data_file11 in files11:
    data11.append(np.loadtxt(data_file11))

files12 = ['/Users/gfranco/class_majoron/output/sup4_2/dcdm_eps_warm2_gam_10H0_nofluid_cl_lensed.dat']
data12 = []
for data_file12 in files12:
    data12.append(np.loadtxt(data_file12))



cl1 = data1[0]
clTT_full_hot_g0p1H0 = interp1d(cl1[:,0], cl1[:,1])
clEE_full_hot_g0p1H0 = interp1d(cl1[:,0], cl1[:,2])

cl2 = data2[0]
clTT_full_hot_gH0 = interp1d(cl2[:,0], cl2[:,1])
clEE_full_hot_gH0 = interp1d(cl2[:,0], cl2[:,2])

cl3 = data3[0]
clTT_full_hot_g10H0 = interp1d(cl3[:,0], cl3[:,1])
clEE_full_hot_g10H0 = interp1d(cl3[:,0], cl3[:,2])

cl4 = data4[0]
clTT_full_warm_g0p1H0 = interp1d(cl4[:,0], cl4[:,1])
clEE_full_warm_g0p1H0 = interp1d(cl4[:,0], cl4[:,2])

cl5 = data5[0]
clTT_full_warm_gH0 = interp1d(cl5[:,0], cl5[:,1])
clEE_full_warm_gH0 = interp1d(cl5[:,0], cl5[:,2])

cl6 = data6[0]
clTT_full_warm_g10H0 = interp1d(cl6[:,0], cl6[:,1])
clEE_full_warm_g10H0 = interp1d(cl6[:,0], cl6[:,2])

cl7 = data7[0]
clTT_full_cold_g0p1H0 = interp1d(cl7[:,0], cl7[:,1])
clEE_full_cold_g0p1H0 = interp1d(cl7[:,0], cl7[:,2])

cl8 = data8[0]
clTT_full_cold_gH0 = interp1d(cl8[:,0], cl8[:,1])
clEE_full_cold_gH0 = interp1d(cl8[:,0], cl8[:,2])

cl9 = data9[0]
clTT_full_cold_g10H0 = interp1d(cl9[:,0], cl9[:,1])
clEE_full_cold_g10H0 = interp1d(cl9[:,0], cl9[:,2])

cl10 = data10[0]
clTT_full_warm2_g0p1H0 = interp1d(cl10[:,0], cl10[:,1])
clEE_full_warm2_g0p1H0 = interp1d(cl10[:,0], cl10[:,2])

cl11 = data11[0]
clTT_full_warm2_gH0 = interp1d(cl11[:,0], cl11[:,1])
clEE_full_warm2_gH0 = interp1d(cl11[:,0], cl11[:,2])

cl12 = data12[0]
clTT_full_warm2_g10H0 = interp1d(cl12[:,0], cl12[:,1])
clEE_full_warm2_g10H0 = interp1d(cl12[:,0], cl12[:,2])
#%% READ FLUID DATA FILES
    
files1b = ['/Users/gfranco/class_majoron/output/sup4/dcdm_eps_hot_gam_0p1H0_fluid_cl_lensed.dat']
data1b = []
for data_file1b in files1b:
    data1b.append(np.loadtxt(data_file1b))

files2b = ['/Users/gfranco/class_majoron/output/sup4/dcdm_eps_hot_gam_H0_fluid_cl_lensed.dat']
data2b = []
for data_file2b in files2b:
    data2b.append(np.loadtxt(data_file2b))

files3b = ['/Users/gfranco/class_majoron/output/sup4/dcdm_eps_hot_gam_10H0_fluid_cl_lensed.dat']
data3b = []
for data_file3b in files3b:
    data3b.append(np.loadtxt(data_file3b))
    
files4b = ['/Users/gfranco/class_majoron/output/sup4/dcdm_eps_warm_gam_0p1H0_fluid_cl_lensed.dat']
data4b = []
for data_file4b in files4b:
    data4b.append(np.loadtxt(data_file4b))
    
files5b = ['/Users/gfranco/class_majoron/output/sup4/dcdm_eps_warm_gam_H0_fluid_cl_lensed.dat']
data5b = []
for data_file5b in files5b:
    data5b.append(np.loadtxt(data_file5b))

files6b = ['/Users/gfranco/class_majoron/output/sup4/dcdm_eps_warm_gam_10H0_fluid_cl_lensed.dat']
data6b = []
for data_file6b in files6b:
    data6b.append(np.loadtxt(data_file6b))
    
files7b = ['/Users/gfranco/class_majoron/output/sup4/dcdm_eps_cold_gam_0p1H0_fluid_cl_lensed.dat']
data7b = []
for data_file7b in files7b:
    data7b.append(np.loadtxt(data_file7b))
      
files8b = ['/Users/gfranco/class_majoron/output/sup4/dcdm_eps_cold_gam_H0_fluid_cl_lensed.dat']
data8b = []
for data_file8b in files8b:
    data8b.append(np.loadtxt(data_file8b))  
    
files9b = ['/Users/gfranco/class_majoron/output/sup4/dcdm_eps_cold_gam_10H0_fluid_cl_lensed.dat']
data9b = []
for data_file9b in files9b:
    data9b.append(np.loadtxt(data_file9b))
    
files10b = ['/Users/gfranco/class_majoron/output/sup4/dcdm_eps_warm2_gam_0p1H0_fluid_cl_lensed.dat']
data10b = []
for data_file10b in files10b:
    data10b.append(np.loadtxt(data_file10b))

files11b = ['/Users/gfranco/class_majoron/output/sup4/dcdm_eps_warm2_gam_H0_fluid_cl_lensed.dat']
data11b = []
for data_file11b in files11b:
    data11b.append(np.loadtxt(data_file11b))
    
files12b = ['/Users/gfranco/class_majoron/output/sup4/dcdm_eps_warm2_gam_10H0_fluid_cl_lensed.dat']
data12b = []
for data_file12b in files12b:
    data12b.append(np.loadtxt(data_file12b))
    
    
cl1b = data1b[0]
clTT_approx_hot_g0p1H0 = interp1d(cl1b[:,0], cl1b[:,1])
clEE_approx_hot_g0p1H0 = interp1d(cl1b[:,0], cl1b[:,2])

cl2b = data2b[0]
clTT_approx_hot_gH0 = interp1d(cl2b[:,0], cl2b[:,1])
clEE_approx_hot_gH0 = interp1d(cl2b[:,0], cl2b[:,2])

cl3b = data3b[0]
clTT_approx_hot_g10H0 = interp1d(cl3b[:,0], cl3b[:,1])
clEE_approx_hot_g10H0 = interp1d(cl3b[:,0], cl3b[:,2])

cl4b = data4b[0]
clTT_approx_warm_g0p1H0 = interp1d(cl4b[:,0], cl4b[:,1])
clEE_approx_warm_g0p1H0 = interp1d(cl4b[:,0], cl4b[:,2])

cl5b = data5b[0]
clTT_approx_warm_gH0 = interp1d(cl5b[:,0], cl5b[:,1])
clEE_approx_warm_gH0 = interp1d(cl5b[:,0], cl5b[:,2])

cl6b = data6b[0]
clTT_approx_warm_g10H0 = interp1d(cl6b[:,0], cl6b[:,1])
clEE_approx_warm_g10H0 = interp1d(cl6b[:,0], cl6b[:,2])

cl7b = data7b[0]
clTT_approx_cold_g0p1H0 = interp1d(cl7b[:,0], cl7b[:,1])
clEE_approx_cold_g0p1H0 = interp1d(cl7b[:,0], cl7b[:,2])

cl8b = data8b[0]
clTT_approx_cold_gH0 = interp1d(cl8b[:,0], cl8b[:,1])
clEE_approx_cold_gH0 = interp1d(cl8b[:,0], cl8b[:,2])

cl9b = data9b[0]
clTT_approx_cold_g10H0 = interp1d(cl9b[:,0], cl9b[:,1])
clEE_approx_cold_g10H0 = interp1d(cl9b[:,0], cl9b[:,2])

cl10b = data10b[0]
clTT_approx_warm2_g0p1H0 = interp1d(cl10b[:,0], cl10b[:,1])
clEE_approx_warm2_g0p1H0 = interp1d(cl10b[:,0], cl10b[:,2])

cl11b = data11b[0]
clTT_approx_warm2_gH0 = interp1d(cl11b[:,0], cl11b[:,1])
clEE_approx_warm2_gH0 = interp1d(cl11b[:,0], cl11b[:,2])

cl12b = data12b[0]
clTT_approx_warm2_g10H0 = interp1d(cl12b[:,0], cl12b[:,1])
clEE_approx_warm2_g10H0 = interp1d(cl12b[:,0], cl12b[:,2])



#%% TIME TO PLOT 

ll = cl1[:,0]

# set plot configuration
ax_1 = plt.subplot(211)
ax_2 = plt.subplot(212, sharex = ax_1)
plt.subplots_adjust(hspace=0)

ax_1.set_ylim([-0.07,0.07])
ax_2.set_ylim([-0.11,0.11])
ax_1.set_xlim([2,2500])
ax_2.set_xlim([2,2500])

ax_2.tick_params(axis='both', which='minor', labelsize=12)

#plot each TT and EE
ax_1.semilogx(ll,clTT_approx_hot_g0p1H0(ll)/clTT_full_hot_g0p1H0(ll)-1,'red')
ax_2.semilogx(ll,clEE_approx_hot_g0p1H0(ll)/clEE_full_hot_g0p1H0(ll)-1,'red')

ax_1.semilogx(ll,clTT_approx_hot_gH0(ll)/clTT_full_hot_gH0(ll)-1,'red',linestyle='dotted')
ax_2.semilogx(ll,clEE_approx_hot_gH0(ll)/clEE_full_hot_gH0(ll)-1,'red',linestyle='dotted')

ax_1.semilogx(ll,clTT_approx_hot_g10H0(ll)/clTT_full_hot_g10H0(ll)-1,'red',linestyle='dashed')
ax_2.semilogx(ll,clEE_approx_hot_g10H0(ll)/clEE_full_hot_g10H0(ll)-1,'red',linestyle='dashed')

ax_1.semilogx(ll,clTT_approx_warm_g0p1H0(ll)/clTT_full_warm_g0p1H0(ll)-1,'green')
ax_2.semilogx(ll,clEE_approx_warm_g0p1H0(ll)/clEE_full_warm_g0p1H0(ll)-1,'green')

ax_1.semilogx(ll,clTT_approx_warm_gH0(ll)/clTT_full_warm_gH0(ll)-1,'green',linestyle='dotted')
ax_2.semilogx(ll,clEE_approx_warm_gH0(ll)/clEE_full_warm_gH0(ll)-1,'green',linestyle='dotted')

ax_1.semilogx(ll,clTT_approx_warm_g10H0(ll)/clTT_full_warm_g10H0(ll)-1,'green',linestyle='dashed')
ax_2.semilogx(ll,clEE_approx_warm_g10H0(ll)/clEE_full_warm_g10H0(ll)-1,'green',linestyle='dashed')

ax_1.semilogx(ll,clTT_approx_warm2_g0p1H0(ll)/clTT_full_warm2_g0p1H0(ll)-1,'blue')
ax_2.semilogx(ll,clEE_approx_warm2_g0p1H0(ll)/clEE_full_warm2_g0p1H0(ll)-1,'blue')

ax_1.semilogx(ll,clTT_approx_warm2_gH0(ll)/clTT_full_warm2_gH0(ll)-1,'blue',linestyle='dotted')
ax_2.semilogx(ll,clEE_approx_warm2_gH0(ll)/clEE_full_warm2_gH0(ll)-1,'blue',linestyle='dotted')

ax_1.semilogx(ll,clTT_approx_warm2_g10H0(ll)/clTT_full_warm2_g10H0(ll)-1,'blue',linestyle='dashed')
ax_2.semilogx(ll,clEE_approx_warm2_g10H0(ll)/clEE_full_warm2_g10H0(ll)-1,'blue',linestyle='dashed')

ax_1.semilogx(ll,clTT_approx_cold_g0p1H0(ll)/clTT_full_cold_g0p1H0(ll)-1,'black')
ax_2.semilogx(ll,clEE_approx_cold_g0p1H0(ll)/clEE_full_cold_g0p1H0(ll)-1,'black')

ax_1.semilogx(ll,clTT_approx_cold_gH0(ll)/clTT_full_cold_gH0(ll)-1,'black',linestyle='dotted')
ax_2.semilogx(ll,clEE_approx_cold_gH0(ll)/clEE_full_cold_gH0(ll)-1,'black',linestyle='dotted')

ax_1.semilogx(ll,clTT_approx_cold_g10H0(ll)/clTT_full_cold_g10H0(ll)-1,'black',linestyle='dashed')
ax_2.semilogx(ll,clEE_approx_cold_g10H0(ll)/clEE_full_cold_g10H0(ll)-1,'black',linestyle='dashed')


# red, green, blue and black correspond to eps=0.5, 0.1, 0.01, 0.001, respectively
# IMPROVE LEGEND

lines = ax_1.get_lines()

black_line1 = mlines.Line2D([], [], color='black', linestyle='solid', label=r'$\Gamma = 0.1 H_0$')
black_line2 = mlines.Line2D([], [], color='black', linestyle='dotted', label=r'$\Gamma = H_0$')
black_line3 = mlines.Line2D([], [], color='black', linestyle='dashed', label=r'$\Gamma = 10 H_0$')


legend1 = ax_1.legend([lines[i] for i in [0,3,6,9]], [r'$\varepsilon = 0.5$', r'$\varepsilon = 0.1$',r'$\varepsilon=0.01$',r'$\varepsilon = 0.001$'], loc=(0.8,0.6), fontsize=13, frameon=False)
legend2 = ax_1.legend(handles= [black_line1,black_line2,black_line3], loc=(0.55, 0.65), fontsize=13, frameon=False)


ax_1.add_artist(legend1)
ax_1.add_artist(legend2)


#ax_1.legend(frameon=False,fontsize = 13,loc='upper right',borderaxespad=0.)

#plot cosmic variance and Planck error bars
l_cosmic_variance = np.linspace(0,48,1000)
l_cosmic_variance_1 = np.linspace(0,30,1000)
l_cosmic_variance_2 = np.linspace(30,48,2)
slope =np.array([0.15,0.0343])
ax_1.fill_between(l_cosmic_variance_1, -0.15,0.15, color='lightgray' )
ax_1.fill_between(l_cosmic_variance_2, -slope, slope, color='lightgray' )
ax_1.fill_between(lTT, -(DlTT_error_plus)/DlTT_mean, +(DlTT_error_plus)/DlTT_mean, color='lightgray')

ax_2.fill_between(l_cosmic_variance, -0.18,0.18, color='lightgray' )
ax_2.fill_between(lEE, -(DlEE_error_plus)/DlEE_mean, +(DlEE_error_plus)/DlEE_mean, color = 'lightgray')


#set labels
ax_2.set_xlabel(r'$\mathrm{multipole} \, \ell$',fontsize=15)
ax_1.set_ylabel(r'$\frac{C_\ell^\mathrm{TT}(\mathrm{approx} )}{C_\ell^\mathrm{TT}(\mathrm{full} )} -1$',fontsize=20)
ax_2.set_ylabel(r'$\frac{C_\ell^\mathrm{EE}(\mathrm{approx} )}{C_\ell^\mathrm{EE}(\mathrm{full} )} -1$',fontsize=20)


ax_2.tick_params(axis="x", labelsize=18)
ax_2.tick_params(axis="y", labelsize=18)
ax_1.tick_params(axis="y", labelsize=18)


plt.show()
plt.clf()

#CONCLUSION: In the fluid configuration, with 300 momentum bins for the superhorizon modes,
# the residuals in C_l remain always below Planck error bars. However, using instead 100
# momentum bins can lead to residuals bigger than Planck error bars (even if I checked
# that the precision for computing S_8 remains the same). Thus, in order to speed up the code,
# it can be nice to use some adaptative method, i.e. that selects 100 or 300 momentum bins
# depending on the values of gamma and epsilon 











