# import necessary modules
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import interp1d

#%%

files1 = ['/Users/gfranco/class_majoron/output/trial_galaxy_lensing_cl.dat']
data1 = []
for data_file1 in files1:
    data1.append(np.loadtxt(data_file1))
    
    
cl_class = data1[0]
ll=cl_class[:,0]
cl_shear_class = interp1d(ll,(ll*(ll+1.0))*cl_class[:,4])

#%%
files2 = ['/Users/gfranco/montepython_public/data/lensing_LCDM_bin9.dat']
data2 = []
for data_file2 in files2:
    data2.append(np.loadtxt(data_file2))
    

cl_shear_euclid = data2[0]
lll=[5.0,6.13,7.51,9.2,11.27,13.81,16.92,20.73,25.4,31.12,38.13,46.73,57.25,70.15,85.95,105.32,129.04,158.11,193.73,237.38,290.85,356.38,436.66,535.03,655.57,803.25,984.21,1205.93,1477.6,1810.48,2218.34,2718.09,3330.42,4080.7,5000.]
    
cl_shear_euclid = interp1d(lll,cl_shear_euclid)

#%%

files3 = ['/Users/gfranco/class_majoron/output/trial_lcdm_cl.dat']
data3 = []
for data_file3 in files3:
    data3.append(np.loadtxt(data_file3))
    
    
cl_lcdm = data1[0]
ll=cl_lcdm[:,0]
cl_shear_lcdm= interp1d(ll,(ll*(ll+1.0))*cl_lcdm[:,4])

files4 = ['/Users/gfranco/class_majoron/output/trial_dcdm_gamma100_f0p03_cl.dat']
data4 = []
for data_file4 in files4:
    data4.append(np.loadtxt(data_file4))
    
    
cl_dcdm_gamma100_f0p03 = data4[0]
ll=cl_dcdm_gamma100_f0p03[:,0]
cl_dcdm_gamma100_f0p03= interp1d(ll,(ll*(ll+1.0))*cl_dcdm_gamma100_f0p03[:,4])


files5 = ['/Users/gfranco/class_majoron/output/trial_dcdm_gamma100_f0p1_cl.dat']
data5 = []
for data_file5 in files5:
    data5.append(np.loadtxt(data_file5))
    
    
cl_dcdm_gamma100_f0p1 = data5[0]
ll=cl_dcdm_gamma100_f0p1[:,0]
cl_dcdm_gamma100_f0p1= interp1d(ll,(ll*(ll+1.0))*cl_dcdm_gamma100_f0p1[:,4])


files6 = ['/Users/gfranco/class_majoron/output/trial_dcdm_gamma100_f0p3_cl.dat']
data6 = []
for data_file6 in files6:
    data6.append(np.loadtxt(data_file6))
    
    
cl_dcdm_gamma100_f0p3 = data6[0]
ll=cl_dcdm_gamma100_f0p3[:,0]
cl_dcdm_gamma100_f0p3= interp1d(ll,(ll*(ll+1.0))*cl_dcdm_gamma100_f0p3[:,4])




#%%
plt.xlabel(r'$\mathrm{multipole} \, \ell$',fontsize=15)
plt.ylabel(r'$\ell (\ell+1) \ C_\ell^{\phi \phi}$',fontsize=15)


plt.tick_params(axis="x", labelsize=18)
plt.tick_params(axis="y", labelsize=18)
plt.title(r'$\bar{z}=2.5, \ \ \Delta z =1$')

#plt.xlim([3,2500])
#plt.semilogx(ll,cl_shear_class(ll),'b',label=r'CLASS')
#plt.semilogx(lll,cl_shear_euclid(lll),'r',label=r'Euclid likelihood')

plt.semilogx(ll,cl_shear_lcdm(ll),'k',label=r'$\Lambda$CDM')
plt.semilogx(ll,cl_dcdm_gamma100_f0p03(ll),'r',label=r'$\Lambda$DDM, f=0.03, $\Gamma =100$ km/s/Mpc')
plt.semilogx(ll,cl_dcdm_gamma100_f0p1(ll),'g',label=r'$\Lambda$DDM, f=0.1, $\Gamma =100$ km/s/Mpc')
plt.semilogx(ll,cl_dcdm_gamma100_f0p3(ll),'b',label=r'$\Lambda$DDM, f=0.3, $\Gamma =100$ km/s/Mpc ')



plt.legend(loc='best', fontsize=15)


plt.show()

plt.clf()