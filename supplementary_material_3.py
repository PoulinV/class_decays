# import necessary modules
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import numpy as np

#%%




data = np.zeros([4,3])
data[0,0]= (0.78459/0.78458)-1.0
data[0,1]= (0.50884/0.50912)-1.0
data[0,2]= (0.18675/0.18742)-1.0
data[1,0]= (0.75641/0.75568)-1.0
data[1,1]= (0.36174/0.35936)-1.0
#data[1,2]= (0.07128/0.06657)-1.0 #NORMAL FLUID
#data[1,2]= (0.0713/0.0659)-1.0  #USING FUDGE FACTOR F=1.2
#data[1,2]= (0.0713/0.0661)-1.0  #USING FUDGE FACTOR F=1.14
data[1,2]= (0.0713/0.0685)-1.0  #USING FUDGE FACTOR F=0.71
data[2,0]= (0.78917/0.79398)-1.0
#data[2,1]= (0.56569/0.59346)-1.0 #NORMAL FLUID
#data[2,1]= (0.56569/0.57528)-1.0 #USING FUDGE FACTOR F=1.2
#data[2,1]= (0.56569/0.58038)-1.0 #USING FUDGE FACTOR F=1.14
data[2,1]= (0.56569/0.62752)-1.0 #USING FUDGE FACTOR F=0.71
#data[2,2]= (0.36981/0.38949)-1.0 #NORMAL FLUID
#data[2,2]= (0.36981/0.36269)-1.0 #USING FUDGE FACTOR F=1.2
#data[2,2]= (0.36981/0.37014)-1.0 #USING FUDGE FACTOR F=1.14
data[2,2]= (0.36981/0.44168)-1.0 #USING FUDGE FACTOR F=0.71
data[3,0]= (0.82608/0.82592)-1.0
data[3,1]= (0.82040/0.82019)-1.0
data[3,2]= (0.81185/0.81221)-1.0



x_start = -1.5
x_end = 1.5
y_start = -2
y_end = 1.5


extent = [x_start, x_end, y_start, y_end]

fig = plt.figure()
ax = fig.add_subplot(111)

ax.set_xlabel(r'$\mathrm{Log}_{10}(\Gamma / H_0)$', fontsize =15)
ax.set_ylabel(r'$\mathrm{Log}_{10}(\varepsilon)$', fontsize = 15)
ax.set_title(r'$S_8^{\mathrm{full}}/S_8^{\mathrm{approx}} -1$', fontsize = 20)

#ax.tick_params(axis="x", labelsize=15)
#ax.tick_params(axis="y", labelsize=15)

plt.locator_params(axis='x', nbins=3)


ax.set_yticks([-1.6,-0.7,0.2,1.1])
ax.set_yticklabels(['$-3$','$-2$','$-1$','$-0.30103$'])


plt.locator_params(axis='y', nbins=4)




im = ax.imshow(data, extent=extent, origin='upper', interpolation='none', cmap='viridis')
#TRY TO PUT COLORBAR IN LOG SCALE 
fig.colorbar(im)
plt.show()


