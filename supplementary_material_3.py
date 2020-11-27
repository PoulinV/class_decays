# import necessary modules
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import numpy as np
import matplotlib.gridspec as gridspec
from scipy.interpolate import interp1d

#%%
############################################# get LCDM reference #############

files = ['/Users/gfranco/class_majoron/output/lcdm_bestfit_pk.dat']
data = []
for data_file in files:
    data.append(np.loadtxt(data_file))

pkk = data[0]  
fpk_LCDM = interp1d(pkk[:,0], pkk[:,1])

#%%



# WITH NEW FITTING FORMULA FOR THE SOUND SPEED
data = np.zeros([4,3])
data[0,0]= (0.78458/0.78459)-1.0
round(data[0,0]*100,3)
data[0,1]= (0.50912/0.50884)-1.0
round(data[0,1]*100,2)
data[0,2]= (0.18742/0.18675)-1.0
round(data[0,2]*100,2)

data[1,0]= (0.75522/0.75641)-1.0
round(data[1,0]*100,2)
data[1,1]= (0.35763/0.36174)-1.0
round(data[1,1]*100,2)
#data[1,2]= (0.06596/0.07128)-1.0

data[1,2]= np.nan


data[2,0]= (0.79042/0.78917)-1.0
round(data[2,0]*100,2)
data[2,1]= (0.57199/0.56569)-1.0
round(data[2,1]*100,2)
data[2,2]= (0.36529/0.36981)-1.0
round(data[2,2]*100,2)

data[3,0]= (0.82584/0.82608)-1.0
round(data[3,0]*100,2)
data[3,1]= (0.81955/0.82040)-1.0
round(data[3,1]*100,2)
data[3,2]= (0.81093/0.81185)-1.0
round(data[3,2]*100,2)




x_start = -1.5
x_end = 1.5
y_start = -2
y_end = 1.5


extent = [x_start, x_end, y_start, y_end]

fig = plt.figure()
ax = fig.add_subplot(111)

ax.set_xlabel(r'$\mathrm{Log}_{10}(\Gamma / H_0)$', fontsize =15)
ax.set_ylabel(r'$\mathrm{Log}_{10}(\varepsilon)$', fontsize = 15)
ax.set_title(r'$S_8^{\mathrm{approx}}/S_8^{\mathrm{full}} -1$', fontsize = 20, pad=15)

#ax.tick_params(axis="x", labelsize=15)
#ax.tick_params(axis="y", labelsize=15)

plt.locator_params(axis='x', nbins=3)


ax.set_yticks([-1.6,-0.7,0.2,1.1])
ax.set_yticklabels(['$-3$','$-2$','$-1$','$-0.30103$'])


plt.locator_params(axis='y', nbins=4)

im = ax.imshow(data, extent=extent, origin='upper', interpolation='none', cmap='viridis')

fig.colorbar(im)
plt.show()


#%% OBTAIN RESIDUALS IN THE MATTER POWER SPECTRUM TODAY

## READ REFERENCE DATA FILES FULL CALCULATION

files1 = ['/Users/gfranco/class_majoron/output/sup4_2/dcdm_eps_hot_gam_0p1H0_nofluid_z1_pk.dat']
data1 = []
for data_file1 in files1:
    data1.append(np.loadtxt(data_file1))

files2 = ['/Users/gfranco/class_majoron/output/sup4_2/dcdm_eps_hot_gam_H0_nofluid_z1_pk.dat']
data2 = []
for data_file2 in files2:
    data2.append(np.loadtxt(data_file2))

files3 = ['/Users/gfranco/class_majoron/output/sup4_2/dcdm_eps_hot_gam_10H0_nofluid_z1_pk.dat']
data3 = []
for data_file3 in files3:
    data3.append(np.loadtxt(data_file3))

files4 = ['/Users/gfranco/class_majoron/output/sup4_2/dcdm_eps_warm_gam_0p1H0_nofluid_z1_pk.dat']
data4 = []
for data_file4 in files4:
    data4.append(np.loadtxt(data_file4))

files5 = ['/Users/gfranco/class_majoron/output/sup4_2/dcdm_eps_warm_gam_H0_nofluid_z1_pk.dat']
data5 = []
for data_file5 in files5:
    data5.append(np.loadtxt(data_file5))
    
files6 = ['/Users/gfranco/class_majoron/output/sup4_2/dcdm_eps_warm_gam_10H0_nofluid_z1_pk.dat']
data6 = []
for data_file6 in files6:
    data6.append(np.loadtxt(data_file6))
    
files7 = ['/Users/gfranco/class_majoron/output/sup4_2/dcdm_eps_cold_gam_0p1H0_nofluid_z1_pk.dat']
data7 = []
for data_file7 in files7:
    data7.append(np.loadtxt(data_file7))  
    
files8 = ['/Users/gfranco/class_majoron/output/sup4_2/dcdm_eps_cold_gam_H0_nofluid_z1_pk.dat']
data8 = []
for data_file8 in files8:
    data8.append(np.loadtxt(data_file8))

files9 = ['/Users/gfranco/class_majoron/output/sup4_2/dcdm_eps_cold_gam_10H0_nofluid_z1_pk.dat']
data9 = []
for data_file9 in files9:
    data9.append(np.loadtxt(data_file9))
    
files10 = ['/Users/gfranco/class_majoron/output/sup4_2/dcdm_eps_warm2_gam_0p1H0_nofluid_z1_pk.dat']
data10 = []
for data_file10 in files10:
    data10.append(np.loadtxt(data_file10))

files11 = ['/Users/gfranco/class_majoron/output/sup4_2/dcdm_eps_warm2_gam_H0_nofluid_z1_pk.dat']
data11 = []
for data_file11 in files11:
    data11.append(np.loadtxt(data_file11))

files12 = ['/Users/gfranco/class_majoron/output/sup4_2/dcdm_eps_warm2_gam_10H0_nofluid_z1_pk.dat']
data12 = []
for data_file12 in files12:
    data12.append(np.loadtxt(data_file12))



pk1 = data1[0]
pkz0_full_hot_g0p1H0 = interp1d(pk1[:,0], pk1[:,1])

pk2 = data2[0]
pkz0_full_hot_gH0 = interp1d(pk2[:,0], pk2[:,1])

pk3 = data3[0]
pkz0_full_hot_g10H0 = interp1d(pk3[:,0], pk3[:,1])

pk4 = data4[0]
pkz0_full_warm_g0p1H0 = interp1d(pk4[:,0], pk4[:,1])

pk5 = data5[0]
pkz0_full_warm_gH0 = interp1d(pk5[:,0], pk5[:,1])

pk6 = data6[0]
pkz0_full_warm_g10H0 = interp1d(pk6[:,0], pk6[:,1])

pk7 = data7[0]
pkz0_full_cold_g0p1H0 = interp1d(pk7[:,0], pk7[:,1])

pk8 = data8[0]
pkz0_full_cold_gH0 = interp1d(pk8[:,0], pk8[:,1])

pk9 = data9[0]
pkz0_full_cold_g10H0 = interp1d(pk9[:,0], pk9[:,1])

pk10 = data10[0]
pkz0_full_warm2_g0p1H0 = interp1d(pk10[:,0], pk10[:,1])

pk11 = data11[0]
pkz0_full_warm2_gH0 = interp1d(pk11[:,0], pk11[:,1])

pk12 = data12[0]
pkz0_full_warm2_g10H0 = interp1d(pk12[:,0], pk12[:,1])


#%% READ FLUID DATA FILES


files1b = ['/Users/gfranco/class_majoron/output/sup4/dcdm_eps_hot_gam_0p1H0_fluid_z1_pk.dat']
data1b = []
for data_file1b in files1b:
    data1b.append(np.loadtxt(data_file1b))

files2b = ['/Users/gfranco/class_majoron/output/sup4/dcdm_eps_hot_gam_H0_fluid_z1_pk.dat']
data2b = []
for data_file2b in files2b:
    data2b.append(np.loadtxt(data_file2b))

files3b = ['/Users/gfranco/class_majoron/output/sup4/dcdm_eps_hot_gam_10H0_fluid_z1_pk.dat']
data3b = []
for data_file3b in files3b:
    data3b.append(np.loadtxt(data_file3b))

files4b = ['/Users/gfranco/class_majoron/output/sup4/dcdm_eps_warm_gam_0p1H0_fluid_z1_pk.dat']
data4b = []
for data_file4b in files4b:
    data4b.append(np.loadtxt(data_file4b))

files5b = ['/Users/gfranco/class_majoron/output/sup4/dcdm_eps_warm_gam_H0_fluid_z1_pk.dat']
data5b = []
for data_file5b in files5b:
    data5b.append(np.loadtxt(data_file5b))
    
files6b = ['/Users/gfranco/class_majoron/output/sup4/dcdm_eps_warm_gam_10H0_fluid_z1_pk.dat']
data6b = []
for data_file6b in files6b:
    data6b.append(np.loadtxt(data_file6b))
    
files7b = ['/Users/gfranco/class_majoron/output/sup4/dcdm_eps_cold_gam_0p1H0_fluid_z1_pk.dat']
data7b = []
for data_file7b in files7b:
    data7b.append(np.loadtxt(data_file7b))  
    
files8b = ['/Users/gfranco/class_majoron/output/sup4/dcdm_eps_cold_gam_H0_fluid_z1_pk.dat']
data8b = []
for data_file8b in files8b:
    data8b.append(np.loadtxt(data_file8b))

files9b = ['/Users/gfranco/class_majoron/output/sup4/dcdm_eps_cold_gam_10H0_fluid_z1_pk.dat']
data9b = []
for data_file9b in files9b:
    data9b.append(np.loadtxt(data_file9b))
    
files10b = ['/Users/gfranco/class_majoron/output/sup4/dcdm_eps_warm2_gam_0p1H0_fluid_z1_pk.dat']
data10b = []
for data_file10b in files10b:
    data10b.append(np.loadtxt(data_file10b))

files11b = ['/Users/gfranco/class_majoron/output/sup4/dcdm_eps_warm2_gam_H0_fluid_z1_pk.dat']
data11b = []
for data_file11b in files11b:
    data11b.append(np.loadtxt(data_file11b))

files12b = ['/Users/gfranco/class_majoron/output/sup4/dcdm_eps_warm2_gam_10H0_fluid_z1_pk.dat']
data12b = []
for data_file12b in files12b:
    data12b.append(np.loadtxt(data_file12b))



pk1b = data1b[0]
pkz0_approx_hot_g0p1H0 = interp1d(pk1b[:,0], pk1b[:,1])

pk2b = data2b[0]
pkz0_approx_hot_gH0 = interp1d(pk2b[:,0], pk2b[:,1])

pk3b = data3b[0]
pkz0_approx_hot_g10H0 = interp1d(pk3b[:,0], pk3b[:,1])

pk4b = data4b[0]
pkz0_approx_warm_g0p1H0 = interp1d(pk4b[:,0], pk4b[:,1])

pk5b = data5b[0]
pkz0_approx_warm_gH0 = interp1d(pk5b[:,0], pk5b[:,1])

pk6b = data6b[0]
pkz0_approx_warm_g10H0 = interp1d(pk6b[:,0], pk6b[:,1])

pk7b = data7b[0]
pkz0_approx_cold_g0p1H0 = interp1d(pk7b[:,0], pk7b[:,1])

pk8b = data8b[0]
pkz0_approx_cold_gH0 = interp1d(pk8b[:,0], pk8b[:,1])

pk9b = data9b[0]
pkz0_approx_cold_g10H0 = interp1d(pk9b[:,0], pk9b[:,1])

pk10b = data10b[0]
pkz0_approx_warm2_g0p1H0 = interp1d(pk10b[:,0], pk10b[:,1])

pk11b = data11b[0]
pkz0_approx_warm2_gH0 = interp1d(pk11b[:,0], pk11b[:,1])

pk12b = data12b[0]
pkz0_approx_warm2_g10H0 = interp1d(pk12b[:,0], pk12b[:,1])


#%%  TIME TO PLOT
#k = pk1[:,0]



fig, (ax_1, ax_2) = plt.subplots(1, 2, gridspec_kw={'width_ratios': [3, 1]})

plt.subplots_adjust(wspace=0)


k = np.logspace(-4,0,1000) # k in h/Mpc


ax_1.set_xscale('log')

ax_1.set_xlim(k[0],k[-1])
ax_1.set_ylim(-0.075,0.075)

ax_1.plot(k,(pkz0_approx_hot_g0p1H0(k)-pkz0_full_hot_g0p1H0(k))/fpk_LCDM(k),'red')

ax_1.plot(k,(pkz0_approx_hot_gH0(k)-pkz0_full_hot_gH0(k))/fpk_LCDM(k),'red',linestyle='dotted')

ax_1.plot(k,(pkz0_approx_hot_g10H0(k)-pkz0_full_hot_g10H0(k))/fpk_LCDM(k),'red',linestyle='dashed')

ax_1.plot(k,(pkz0_approx_warm_g0p1H0(k)-pkz0_full_warm_g0p1H0(k))/fpk_LCDM(k),'green')

ax_1.plot(k,(pkz0_approx_warm_gH0(k)-pkz0_full_warm_gH0(k))/fpk_LCDM(k),'green',linestyle='dotted')

#ax.plot(k,(pkz0_approx_warm_g10H0(k)-pkz0_full_warm_g10H0(k))/fpk_LCDM(k),'green',linestyle='dashed')
## ERRORS (ON S_8) IN THIS CASE ARE BIG , BUT SHOULDN'T BE A PROBLEM, EXCLUDED BY PLANCK

ax_1.plot(k,(pkz0_approx_warm2_g0p1H0(k)-pkz0_full_warm2_g0p1H0(k))/fpk_LCDM(k),'blue')

ax_1.plot(k,(pkz0_approx_warm2_gH0(k)-pkz0_full_warm2_gH0(k))/fpk_LCDM(k),'blue',linestyle='dotted')

ax_1.plot(k,(pkz0_approx_warm2_g10H0(k)-pkz0_full_warm2_g10H0(k))/fpk_LCDM(k),'blue',linestyle='dashed')

ax_1.plot(k,(pkz0_approx_cold_g0p1H0(k)-pkz0_full_cold_g0p1H0(k))/fpk_LCDM(k),'black')

ax_1.plot(k,(pkz0_approx_cold_gH0(k)-pkz0_full_cold_gH0(k))/fpk_LCDM(k),'black',linestyle='dotted')

ax_1.plot(k,(pkz0_approx_cold_g10H0(k)-pkz0_full_cold_g10H0(k))/fpk_LCDM(k),'black',linestyle='dashed')

ax_1.set_xlabel(r'$k \,\,\,\, [h/\mathrm{Mpc}]$', fontsize=15)
ax_1.set_ylabel(r'$(P_{\mathrm{approx}}-P_{\mathrm{full}})/P_{\Lambda\mathrm{CDM}}$', fontsize=20)
ax_1.set_title(r'$z=0$',fontsize=15)


k_range_sigma8 = np.linspace(0.1,0.9,1000) #which are the exact wavenumbers probed by DES-Y1?
ax_1.fill_between(k_range_sigma8, -0.2,0.2, color='lightgray' )

lines = ax_1.get_lines()

black_line1 = mlines.Line2D([], [], color='black', linestyle='solid', label=r'$\Gamma = 0.1 H_0$')
black_line2 = mlines.Line2D([], [], color='black', linestyle='dotted', label=r'$\Gamma = H_0$')
black_line3 = mlines.Line2D([], [], color='black', linestyle='dashed', label=r'$\Gamma = 10 H_0$')

legend1 = ax_1.legend([lines[i] for i in [0,3,5,8]], [r'$\varepsilon = 0.5$', r'$\varepsilon = 0.1$',r'$\varepsilon=0.01$',r'$\varepsilon = 0.001$'], loc='lower left', fontsize=13, frameon=False)

legend2 = ax_1.legend(handles= [black_line1,black_line2,black_line3], loc='upper left', fontsize=13, frameon=False)


ax_1.add_artist(legend1)
ax_1.add_artist(legend2)

ax_1.tick_params(axis="x", labelsize=16)
ax_1.tick_params(axis="y", labelsize=16)


ax_2.tick_params(axis='both', which='both', bottom='False', top='False', labelbottom='False', right='False', left='False', labelleft='False')

ax_2.text(0.15,0.8,r"$|S_8^{\mathrm{approx}}/S_8^{\mathrm{full}} -1|$", fontsize=13)

l1 = mlines.Line2D([], [], color='red', linestyle='solid', label=r'$0.001 \ \%$')
l2 = mlines.Line2D([], [], color='red', linestyle='dotted', label=r'$0.06 \ \% $')
l3 = mlines.Line2D([], [], color='red', linestyle='dashed', label=r'$0.36 \ \% $')

l4 = mlines.Line2D([], [], color='green', linestyle='solid', label=r'$0.16 \ \%$')
l5 = mlines.Line2D([], [], color='green', linestyle='dotted', label=r'$1.14 \ \% $')

l6 = mlines.Line2D([], [], color='blue', linestyle='solid', label=r'$0.16 \ \%$')
l7 = mlines.Line2D([], [], color='blue', linestyle='dotted', label=r'$1.11 \ \% $')
l8 = mlines.Line2D([], [], color='blue', linestyle='dashed', label=r'$1.22 \ \% $')

l9 = mlines.Line2D([], [], color='black', linestyle='solid', label=r'$0.03 \ \%$')
l10 = mlines.Line2D([], [], color='black', linestyle='dotted', label=r'$0.1 \ \% $')
l11 = mlines.Line2D([], [], color='black', linestyle='dashed', label=r'$0.11 \ \% $')


legend3 = ax_2.legend(handles= [l1,l2,l3,l4,l5,l6,l7,l8,l9,l10,l11], loc=(0.2,0.2), fontsize=13, frameon=False)

ax_2.add_artist(legend3)

plt.show()
plt.clf()


#%%  PLOT PERTURBATIONS, delta_wdm at two different wavenumbers, and for 4 different combinations of gamma and epsilon

#filesA = ['/Users/gfranco/class_majoron/output/sup4_2/dcdm_eps_warm2_gam_10H0_nofluid_perturbations_k1_s.dat']
#dataA = []
#for data_fileA in filesA:
#    dataA.append(np.loadtxt(data_fileA))
    
#filesB = ['/Users/gfranco/class_majoron/output/sup4_2/dcdm_eps_warm2_gam_10H0_nofluid_perturbations_k3_s.dat']
#dataB = []
#for data_fileB in filesB:
#    dataB.append(np.loadtxt(data_fileB))

#%% PLOTS for parameters at the 2sigma limit constraints from Planck+SNIa+BAO
kk = np.logspace(-4,0,1000) # k in h/Mpc


###   read reference data files FULL CALCULATION

filesA = ['/Users/gfranco/class_majoron/output/sup4_2/dcdm_eps_0p002_tau_1gyr_nofluid_z1_pk.dat']
dataA = []
for data_fileA in filesA:
    dataA.append(np.loadtxt(data_fileA))
    

filesB = ['/Users/gfranco/class_majoron/output/sup4_2/dcdm_eps_0p008_tau_32gyrs_nofluid_z1_pk.dat']
dataB = []
for data_fileB in filesB:
    dataB.append(np.loadtxt(data_fileB))
    
filesC = ['/Users/gfranco/class_majoron/output/sup4_2/dcdm_eps_0p0032_tau_10gyrs_nofluid_z1_pk.dat']
dataC = []
for data_fileC in filesC:
    dataC.append(np.loadtxt(data_fileC))
    
pkA = dataA[0]
pkz0_full_eps_0p002_tau_1gyr = interp1d(pkA[:,0], pkA[:,1])

pkB = dataB[0]
pkz0_full_eps_0p008_tau_32gyrs = interp1d(pkB[:,0], pkB[:,1])

pkC = dataC[0]
pkz0_full_eps_0p0032_tau_10gyrs = interp1d(pkC[:,0], pkC[:,1])
###################################################################################################

# read data files FLUID CALCULATION 

filesAa = ['/Users/gfranco/class_majoron/output/sup4/dcdm_eps_0p002_tau_1gyr_fluid_z1_pk.dat']
dataAa = []
for data_fileAa in filesAa:
    dataAa.append(np.loadtxt(data_fileAa))
    

filesBb = ['/Users/gfranco/class_majoron/output/sup4/dcdm_eps_0p008_tau_32gyrs_fluid_z1_pk.dat']
dataBb = []
for data_fileBb in filesBb:
    dataBb.append(np.loadtxt(data_fileBb))
    
filesCc = ['/Users/gfranco/class_majoron/output/sup4/dcdm_eps_0p0032_tau_10gyrs_fluid_z1_pk.dat']
dataCc = []
for data_fileCc in filesCc:
    dataCc.append(np.loadtxt(data_fileCc))
    
pkAa = dataAa[0]
pkz0_approx_eps_0p002_tau_1gyr = interp1d(pkAa[:,0], pkAa[:,1])

pkBb = dataBb[0]
pkz0_approx_eps_0p008_tau_32gyrs = interp1d(pkBb[:,0], pkBb[:,1])

pkCc = dataCc[0]
pkz0_approx_eps_0p0032_tau_10gyrs = interp1d(pkCc[:,0], pkCc[:,1])



plt.xscale('log')

plt.xlim(kk[0],kk[-1])
#plt.ylim(-0.1,+0.1)


#plt.plot(kk,(pkz0_approx_eps_0p002_tau_1gyr(kk)-pkz0_full_eps_0p002_tau_1gyr(kk))/fpk_LCDM(kk),'blue', label=r'$\tau = 1 \, \mathrm{Gyrs}, \, \, \, \varepsilon = 0.002, \, \, \, \, S_8 \, \mathrm{residuals} = 0.37 \%$')
#plt.plot(kk,(pkz0_approx_eps_0p0032_tau_10gyrs(kk)-pkz0_full_eps_0p0032_tau_10gyrs(kk))/fpk_LCDM(kk),'green', label=r'$\tau = 10 \, \mathrm{Gyrs}, \, \, \, \varepsilon = 0.0032, \, \, \, \, S_8 \, \mathrm{residuals} = 0.88 \%$')
#plt.plot(kk,(pkz0_approx_eps_0p008_tau_32gyrs(kk)-pkz0_full_eps_0p008_tau_32gyrs(kk))/fpk_LCDM(kk),'red', label=r'$\tau = 32 \, \mathrm{Gyrs}, \, \, \, \varepsilon = 0.008, \, \, \, \, S_8 \, \mathrm{residuals} = 0.89 \%$')


plt.semilogy(kk,pkz0_approx_eps_0p002_tau_1gyr(kk),'red', label=r'$\tau = 10 \, \mathrm{Gyrs}, \, \, \, \varepsilon = 0.0032 \,, \, \, \, S_8 \, \mathrm{residuals} = 0.88 \%$')
plt.semilogy(kk,pkz0_full_eps_0p002_tau_1gyr(kk),'blue', label=r'$\tau = 10 \, \mathrm{Gyrs}, \, \, \, \varepsilon = 0.0032 \,, \, \, \, S_8 \, \mathrm{residuals} = 0.88 \%$')
plt.semilogy(kk,fpk_LCDM(kk),'green', label=r'LCDM')



plt.xlabel(r'$k \,\,\,\, [h/\mathrm{Mpc}]$', fontsize=15)
plt.ylabel(r'$\frac{P_{\mathrm{approx}}-P_{\mathrm{full}}}{P_{\Lambda\mathrm{CDM}}}$', fontsize=20)
plt.title(r'$z=0$',fontsize=15)
k_range_sigma8 = np.linspace(0.1,0.9,1000) #which are the exact wavenumbers probed by DES-Y1?
plt.fill_between(k_range_sigma8, -10000,10000, color='lightgray' )

plt.legend(loc='upper left', fontsize=13)

plt.tick_params(axis="x", labelsize=18)
plt.tick_params(axis="y", labelsize=18)

plt.show()
plt.clf()



  
    

  
    

   


