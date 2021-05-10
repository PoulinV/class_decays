# import necessary modules
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import interp1d
import matplotlib.lines as mlines
#%% READ CL_SHEAR (LCDM) FROM CLASS

lll=[5.0,6.13,7.51,9.2,11.27,13.81,16.92,20.73,25.4,31.12,38.13,46.73,57.25,70.15,85.95,105.32,129.04,158.11,193.73,237.38,290.85,356.38,436.66,535.03,655.57,803.25,984.21,1205.93,1477.6,1810.48,2218.34,2718.09,3330.42,4080.7,5000.]    

lll[26]


files = ['/Users/gfranco/class_majoron/output/trial_lcdm_cl.dat']
data = []
for data_file in files:
    data.append(np.loadtxt(data_file))
cl_class = data[0]
ll=cl_class[:,0]
cl_shear_class = interp1d(ll,(ll*(ll+1.0))*cl_class[:,4])

### READ CL_SHEAR (LCDM) FROM EUCLID LIKELIHOOD
files = ['/Users/gfranco/montepython_public/data/lensing_LCDM_bin9.dat']
data = []
for data_file in files:
    data.append(np.loadtxt(data_file))
cl_lcdm_bin9 = interp1d(lll,data[0])


files = ['/Users/gfranco/montepython_public/data/lensing_LCDM_bin3.dat']
data = []
for data_file in files:
    data.append(np.loadtxt(data_file))
cl_lcdm_bin3 = interp1d(lll,data[0])

#%% READ CL_SHEAR (DCDM) FROM EUCLID LIKELIHOOD

#files = ['/Users/gfranco/montepython_public/data/lensing_DCDM_gam100_f0p03_bin9.dat']
files = ['/Users/gfranco/montepython_public/data/lensing_DCDM_gam50_f0p1_bin9.dat']
data = []
for data_file in files:
    data.append(np.loadtxt(data_file))   
#cl_dcdm_gamma100_f0p03_bin9= interp1d(lll,data[0])
cl_dcdm_gamma50_f0p1_bin9 = interp1d(lll,data[0])


files = ['/Users/gfranco/montepython_public/data/lensing_DCDM_gam50_f0p1_bin3.dat']
data = []
for data_file in files:
    data.append(np.loadtxt(data_file))
cl_dcdm_gamma50_f0p1_bin3 = interp1d(lll,data[0])
 

files = ['/Users/gfranco/montepython_public/data/lensing_DCDM_gam100_f0p1_bin9.dat']
data = []
for data_file in files:
    data.append(np.loadtxt(data_file))
cl_dcdm_gamma100_f0p1_bin9= interp1d(lll, data[0])


files = ['/Users/gfranco/montepython_public/data/lensing_DCDM_gam100_f0p1_bin3.dat']
data = []
for data_file in files:
    data.append(np.loadtxt(data_file))
cl_dcdm_gamma100_f0p1_bin3= interp1d(lll, data[0])

#files = ['/Users/gfranco/montepython_public/data/lensing_DCDM_gam100_f0p3_bin9.dat']
files = ['/Users/gfranco/montepython_public/data/lensing_DCDM_gam150_f0p1_bin9.dat']
data = []
for data_file in files:
    data.append(np.loadtxt(data_file))
cl_dcdm_gamma150_f0p1_bin9 = interp1d(lll, data[0])
#cl_dcdm_gamma100_f0p3_bin9= interp1d(lll,data[0])


files = ['/Users/gfranco/montepython_public/data/lensing_DCDM_gam150_f0p1_bin3.dat']
data = []
for data_file in files:
    data.append(np.loadtxt(data_file))
cl_dcdm_gamma150_f0p1_bin3 = interp1d(lll, data[0])


#%%  READ Window function W(Z) FROM EUCLID LIKELIHOOD


## FOR DCDM
files = ['/Users/gfranco/montepython_public/data/W_of_z_gam50_f0p1_bin9.dat']
data = []
for data_file in files:
    data.append(np.loadtxt(data_file))  
Wz = data[0]    
Wz_gam50_f0p1_bin9 = interp1d(Wz[:,0],Wz[:,1])

files = ['/Users/gfranco/montepython_public/data/W_of_z_gam100_f0p1_bin9.dat']
data = []
for data_file in files:
    data.append(np.loadtxt(data_file))  
Wz = data[0]    
Wz_gam100_f0p1_bin9 = interp1d(Wz[:,0],Wz[:,1])


files = ['/Users/gfranco/montepython_public/data/W_of_z_gam150_f0p1_bin9.dat']
data = []
for data_file in files:
    data.append(np.loadtxt(data_file))  
Wz = data[0]    
Wz_gam150_f0p1_bin9 = interp1d(Wz[:,0],Wz[:,1])


# AND FOR LCDM

files = ['/Users/gfranco/montepython_public/data/W_of_z_lcdm_bin9.dat']
data = []
for data_file in files:
    data.append(np.loadtxt(data_file))  
Wz = data[0]    
Wz_lcdm_bin9 = interp1d(Wz[:,0],Wz[:,1])



#%% READ PK (DCDM) FROM CLASS

#files= ['/Users/gfranco/class_majoron/output/trial_dcdm_gamma100_f0p03_z1_pk.dat']
files= ['/Users/gfranco/class_majoron/output/trial_dcdm_gamma50_f0p1_z1_pk.dat']
data = []
for data_file in files:
    data.append(np.loadtxt(data_file))
pk = data[0]
#pk_z0_dcdm_f0p03_gam100 = interp1d(pk[:,0],pk[:,1])
pk_z0_dcdm_f0p1_gam50 = interp1d(pk[:,0],pk[:,1])


#files= ['/Users/gfranco/class_majoron/output/trial_dcdm_gamma100_f0p03_z2_pk.dat']
files= ['/Users/gfranco/class_majoron/output/trial_dcdm_gamma50_f0p1_z2_pk.dat']
data = []
for data_file in files:
    data.append(np.loadtxt(data_file))   
pk = data[0]
#pk_z1_dcdm_f0p03_gam100 = interp1d(pk_A[:,0],pk_A[:,1])
pk_z1_dcdm_f0p1_gam50 = interp1d(pk[:,0],pk[:,1])


#files = ['/Users/gfranco/class_majoron/output/trial_dcdm_gamma100_f0p03_background.dat']
files = ['/Users/gfranco/class_majoron/output/trial_dcdm_gamma50_f0p1_background.dat']
data = []
for data_file in files:
    data.append(np.loadtxt(data_file))   
back = data[0]
#Hz_dcdm_f0p03_gam100 = interp1d(back[:,0],back[:,3])
Hz_dcdm_f0p1_gam50 = interp1d(back[:,0],back[:,3])

###########################################################################################

#files = ['/Users/gfranco/class_majoron/output/trial_dcdm_gamma100_f0p1_z1_pk.dat']
files = ['/Users/gfranco/class_majoron/output/trial_dcdm_gamma100_f0p1_z1_pk.dat']
data = []
for data_file in files:
    data.append(np.loadtxt(data_file))
pk = data[0]
#pk_z0_dcdm_f0p1_gam100 = interp1d(pk[:,0],pk[:,1])
pk_z0_dcdm_f0p1_gam100 = interp1d(pk[:,0],pk[:,1])


#files = ['/Users/gfranco/class_majoron/output/trial_dcdm_gamma100_f0p1_z2_pk.dat']
files = ['/Users/gfranco/class_majoron/output/trial_dcdm_gamma100_f0p1_z2_pk.dat']
data = []
for data_file in files:
    data.append(np.loadtxt(data_file))   
pk = data[0]
#pk_z1_dcdm_f0p1_gam100 = interp1d(pk[:,0],pk[:,1])
pk_z1_dcdm_f0p1_gam100 = interp1d(pk[:,0],pk[:,1])


#files = ['/Users/gfranco/class_majoron/output/trial_dcdm_gamma100_f0p1_background.dat']
files = ['/Users/gfranco/class_majoron/output/trial_dcdm_gamma100_f0p1_background.dat']
data  = []
for data_file in files:
    data.append(np.loadtxt(data_file))
back = data[0]
#Hz_dcdm_f0p1_gam100 = interp1d(back[:,0],back[:,3])
Hz_dcdm_f0p1_gam100 = interp1d(back[:,0],back[:,3])

###############################################################################

#files = ['/Users/gfranco/class_majoron/output/trial_dcdm_gamma100_f0p3_z1_pk.dat']
files = ['/Users/gfranco/class_majoron/output/trial_dcdm_gamma150_f0p1_z1_pk.dat']

data = []
for data_file in files:
    data.append(np.loadtxt(data_file))
pk = data[0]
#pk_z0_dcdm_f0p3_gam100 = interp1d(pk[:,0],pk[:,1])
pk_z0_dcdm_f0p1_gam150 = interp1d(pk[:,0],pk[:,1])



#files = ['/Users/gfranco/class_majoron/output/trial_dcdm_gamma100_f0p3_z2_pk.dat']
files = ['/Users/gfranco/class_majoron/output/trial_dcdm_gamma150_f0p1_z2_pk.dat']
data = []
for data_file in files:
    data.append(np.loadtxt(data_file))   
pk = data[0]
#pk_z1_dcdm_f0p3_gam100 = interp1d(pk[:,0],pk[:,1])
pk_z1_dcdm_f0p1_gam150 = interp1d(pk[:,0],pk[:,1])



#files = ['/Users/gfranco/class_majoron/output/trial_dcdm_gamma100_f0p3_background.dat']
files  = ['/Users/gfranco/class_majoron/output/trial_dcdm_gamma150_f0p1_background.dat']
data = []
for data_file in files:
    data.append(np.loadtxt(data_file))
back = data[0]
#Hz_dcdm_f0p3_gam100 = interp1d(back[:,0],back[:,3])
Hz_dcdm_f0p1_gam150 = interp1d(back[:,0],back[:,3])


## READ PK (LCDM) FROM CLASS ####################################################

files = ['/Users/gfranco/class_majoron/output/trial_lcdm_z1_pk.dat']
data = []
for data_file in files:
    data.append(np.loadtxt(data_file))
pk = data[0]
pk_z0_lcdm = interp1d(pk[:,0],pk[:,1])


files = ['/Users/gfranco/class_majoron/output/trial_lcdm_z2_pk.dat']
data = []
for data_file in files:
    data.append(np.loadtxt(data_file))
pk = data[0]
pk_z1_lcdm = interp1d(pk[:,0],pk[:,1])


files = ['/Users/gfranco/class_majoron/output/trial_lcdm_background.dat']
data  = []
for data_file in files:
    data.append(np.loadtxt(data_file))
back = data[0]
Hz_lcdm = interp1d(back[:,0],back[:,3])



#%% PLOT CL_SHEAR


fig, axes = plt.subplots()

axes.set_xlabel(r'$\mathrm{multipole} \, \ell$',fontsize=15)
axes.set_ylabel(r'$ C_\ell^{ii}/ C_\ell^{ii} (\Lambda\mathrm{CDM})-1$', fontsize=20)
#axes.set_ylabel(r'$\ell (\ell+1) \ C_\ell^{\phi \phi}$',fontsize=15)
axes.set_xlim([5,5000])

#LCDM: CLASS VS EUCLID
#axes.semilogx(ll,cl_shear_class(ll),'b',label=r'CLASS, $\bar{z}=2.5, \ \ \Delta z =1$')
#axes.semilogx(lll,cl_lcdm_bin9(lll),'r',label=r'Euclid likelihood, $1.58 <z< 3.5$')

#FROM CLASS
#DCDM: SEVERAL F
#axes.semilogx(ll,cl_shear_class(ll),'k',label=r'$\Lambda$CDM')
#axes.semilogx(lll,cl_lcdm_bin9(lll),'k',label=r'$\Lambda$CDM')

#axes.semilogx(lll,cl_dcdm_gamma100_f0p03_bin9(lll)/cl_lcdm_bin9(lll)-1.0,'r',label=r'$\Lambda$DDM, f=0.03, $\Gamma =100$ km/s/Mpc')
#axes.semilogx(lll,cl_dcdm_gamma100_f0p1_bin9(lll)/cl_lcdm_bin9(lll)-1.0,'g',label=r'$\Lambda$DDM, f=0.1, $\Gamma =100$ km/s/Mpc')
#axes.semilogx(lll,cl_dcdm_gamma100_f0p3_bin9(lll)/cl_lcdm_bin9(lll)-1.0,'b',label=r'$\Lambda$DDM, f=0.3, $\Gamma =100$ km/s/Mpc ')

axes.semilogx(lll,cl_dcdm_gamma50_f0p1_bin9(lll)/cl_lcdm_bin9(lll)-1.0,'r')
axes.semilogx(lll,cl_dcdm_gamma100_f0p1_bin9(lll)/cl_lcdm_bin9(lll)-1.0,'g')
axes.semilogx(lll,cl_dcdm_gamma150_f0p1_bin9(lll)/cl_lcdm_bin9(lll)-1.0,'b')

axes.semilogx(lll,cl_dcdm_gamma50_f0p1_bin3(lll)/cl_lcdm_bin3(lll)-1.0,'r--')
axes.semilogx(lll,cl_dcdm_gamma100_f0p1_bin3(lll)/cl_lcdm_bin3(lll)-1.0,'g--')
axes.semilogx(lll,cl_dcdm_gamma150_f0p1_bin3(lll)/cl_lcdm_bin3(lll)-1.0,'b--')

color_line1 = mlines.Line2D([], [], color='red', linestyle='solid', label=r'$\Gamma = 50 \ \mathrm{km} /\mathrm{s} / \mathrm{Mpc}$')
color_line2 = mlines.Line2D([], [], color='green', linestyle='solid', label=r'$\Gamma = 100 \ \mathrm{km} /\mathrm{s} / \mathrm{Mpc}$')
color_line3 = mlines.Line2D([], [], color='blue', linestyle='solid', label=r'$\Gamma = 150 \ \mathrm{km} /\mathrm{s} / \mathrm{Mpc}$')


black_line1 = mlines.Line2D([], [], color='black', linestyle='solid', label=r'$\mathrm{Bin} \ i=9: \ 1.58 <z< 3.5$')
black_line2 = mlines.Line2D([], [], color='black', linestyle='dashed', label=r'$\mathrm{Bin} \ i=3: \ 0.68 < z < 0.79$')

legend1 = plt.legend(handles= [color_line1,color_line2, color_line3], loc=(0.05,0.75), fontsize=13, frameon=False)
legend2 = plt.legend(handles= [black_line2,black_line1], loc=(0.65,0.1), fontsize=13, frameon=False)

axes.add_artist(legend1)
axes.add_artist(legend2)

plt.title(r'$f_{\mathrm{dcdm}} = 0.1$', fontsize=15)

axes.tick_params(axis="x", labelsize=18)
axes.tick_params(axis="y", labelsize=18)

plt.show()
plt.clf()

#%% PLOT PK

kk = np.logspace(-4,0,1000) # k in h/Mpc

#plt.title(r'$z=0$',fontsize=15)

plt.title(r'$z=1$',fontsize=15)

plt.xlim(kk[0],kk[-1])

plt.xlabel(r'$k \,\,\,\, [h/\mathrm{Mpc}]$', fontsize=15)
plt.ylabel(r'$P(k)/P(k)_{\Lambda\mathrm{CDM}}-1$', fontsize=20)


plt.xscale('log')


#plt.plot(kk,pk_z0_dcdm_f0p03_gam100(kk)/pk_z0_lcdm(kk)-1.0,'r', label=r'$\Lambda \mathrm{DDM}, f=0.03, \Gamma$ =100 km/s/Mpc')
#plt.plot(kk,pk_z0_dcdm_f0p1_gam100(kk)/pk_z0_lcdm(kk)-1.0,'g', label=r'$\Lambda \mathrm{DDM}, f=0.1, \Gamma$ =100 km/s/Mpc')
#plt.plot(kk,pk_z0_dcdm_f0p3_gam100(kk)/pk_z0_lcdm(kk)-1.0,'b', label=r'$\Lambda \mathrm{DDM}, f=0.3, \Gamma$ =100 km/s/Mpc')

#plt.plot(kk,pk_z0_dcdm_f0p1_gam50(kk)/pk_z0_lcdm(kk)-1.0,'r', label=r'$\Lambda \mathrm{DDM}, f=0.1, \Gamma$ = 50 km/s/Mpc')
#plt.plot(kk,pk_z0_dcdm_f0p1_gam100(kk)/pk_z0_lcdm(kk)-1.0,'g', label=r'$\Lambda \mathrm{DDM}, f=0.1, \Gamma$ =100 km/s/Mpc')
#plt.plot(kk,pk_z0_dcdm_f0p1_gam150(kk)/pk_z0_lcdm(kk)-1.0,'b', label=r'$\Lambda \mathrm{DDM}, f=0.1, \Gamma$ =150 km/s/Mpc')


#plt.plot(kk,pk_z1_dcdm_f0p03_gam100(kk)/pk_z1_lcdm(kk)-1.0,'r', label=r'$\Lambda \mathrm{DDM}, f=0.03, \Gamma$ =100 km/s/Mpc')
#plt.plot(kk,pk_z1_dcdm_f0p1_gam100(kk)/pk_z1_lcdm(kk)-1.0,'g', label=r'$\Lambda \mathrm{DDM}, f=0.1, \Gamma$ =100 km/s/Mpc')
#plt.plot(kk,pk_z1_dcdm_f0p3_gam100(kk)/pk_z1_lcdm(kk)-1.0,'b', label=r'$\Lambda \mathrm{DDM}, f=0.3, \Gamma$ =100 km/s/Mpc')


plt.plot(kk,pk_z1_dcdm_f0p1_gam50(kk)/pk_z1_lcdm(kk)-1.0,'r', label=r'$\Lambda \mathrm{DDM}, f=0.1, \Gamma$ =50 km/s/Mpc')
plt.plot(kk,pk_z1_dcdm_f0p1_gam100(kk)/pk_z1_lcdm(kk)-1.0,'g', label=r'$\Lambda \mathrm{DDM}, f=0.1, \Gamma$ =100 km/s/Mpc')
plt.plot(kk,pk_z1_dcdm_f0p1_gam150(kk)/pk_z1_lcdm(kk)-1.0,'b', label=r'$\Lambda \mathrm{DDM}, f=0.1, \Gamma$ =150 km/s/Mpc')


plt.legend(loc='best', fontsize=15)

plt.tick_params(axis="x", labelsize=18)
plt.tick_params(axis="y", labelsize=18)


#plt.show()

plt.clf()


#%% plot expansion rate 

zzz=np.linspace(0,2.5, 1e4)

plt.xlim([0,2.5])

conversion=2.9979e5


#plt.plot(zzz,conversion*Hz_lcdm(zzz)/(1.+zzz) , 'black', label=r'$\Lambda \mathrm{CDM}$')
#plt.plot(zzz,conversion*Hz_dcdm_f0p03_gam100(zzz)/(1.+zzz) , 'green', label=r'$\Lambda \mathrm{DDM}, f=0.03, \Gamma$ =100 km/s/Mpc')
#plt.plot(zzz,conversion*Hz_dcdm_f0p1_gam100(zzz)/(1.+zzz) , 'red', label=r'$\Lambda \mathrm{DDM}, f=0.1, \Gamma$ =100 km/s/Mpc')
#plt.plot(zzz,conversion*Hz_dcdm_f0p3_gam100(zzz)/(1.+zzz) , 'blue', label=r'$\Lambda \mathrm{DDM}, f=0.3, \Gamma$ =100 km/s/Mpc')

plt.plot(zzz,conversion*Hz_lcdm(zzz)/(1.+zzz) , 'black', label=r'$\Lambda \mathrm{CDM}$')
plt.plot(zzz,conversion*Hz_dcdm_f0p1_gam50(zzz)/(1.+zzz) , 'green', label=r'$\Lambda \mathrm{DDM}, f=0.1, \Gamma$ =50 km/s/Mpc')
plt.plot(zzz,conversion*Hz_dcdm_f0p1_gam100(zzz)/(1.+zzz) , 'red', label=r'$\Lambda \mathrm{DDM}, f=0.1, \Gamma$ =100 km/s/Mpc')
plt.plot(zzz,conversion*Hz_dcdm_f0p1_gam150(zzz)/(1.+zzz) , 'blue', label=r'$\Lambda \mathrm{DDM}, f=0.1, \Gamma$ =150 km/s/Mpc')


plt.xlabel(r'$z$', fontsize=18)
plt.ylabel(r'$H(z)/(1+z) \ \ \ \ [\mathrm{km}/\mathrm{s}/\mathrm{Mpc}]$', fontsize=18)


plt.legend(loc='best', frameon=False,borderaxespad=0., fontsize=13)



plt.tick_params(axis='x',labelsize=17)
plt.tick_params(axis='y',labelsize=17)


#plt.show()
plt.clf()


#%% plot window function

zzz=np.linspace(0.05,2.9, 1e4)
plt.xlim([0.05,2.9])
plt.plot(zzz, Wz_lcdm_bin9(zzz), 'black', label=r'$\Lambda \mathrm{CDM}$')
plt.plot(zzz, Wz_gam50_f0p1_bin9(zzz), 'green', label=r'$\Gamma = 50 \ \mathrm{km} /\mathrm{s} / \mathrm{Mpc}$')
plt.plot(zzz, Wz_gam100_f0p1_bin9(zzz), 'red', label=r'$\Gamma = 100 \ \mathrm{km} /\mathrm{s} / \mathrm{Mpc}$')
plt.plot(zzz, Wz_gam150_f0p1_bin9(zzz), 'blue', label=r'$\Gamma = 150 \ \mathrm{km} /\mathrm{s} / \mathrm{Mpc}$')


plt.xlabel(r'$z$', fontsize=18)
plt.ylabel(r'$W_9 (z)$', fontsize=18)


plt.legend(loc='best', frameon=False,borderaxespad=0., fontsize=13)

plt.tick_params(axis='x',labelsize=17)
plt.tick_params(axis='y',labelsize=17)


plt.title(r'$f_{\mathrm{dcdm}} = 0.1$', fontsize=15)

plt.show()
plt.clf()

 



