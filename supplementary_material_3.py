# import necessary modules
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import numpy as np

#%%



size = 3
data = np.zeros([size,size])
data[0,0]= 0.78459-0.78458
data[0,1]= 0.50884-0.50912
data[0,2]= 0.18675-0.18742
data[1,0]= 0.75641-0.75568
data[1,1]= 0.36174-0.35936
data[1,2]= np.nan
data[2,0]= np.nan
data[2,1]= 0.82040-0.82019
data[2,2]= 0.81185-0.81221

data

x_start = -1.5
x_end = 1.5
y_start = -3.7
y_end = 0.8

extraticks=[-0.301, -1, -3]

extent = [x_start, x_end, y_start, y_end]

fig = plt.figure()
ax = fig.add_subplot(111)

ax.set_xlabel(r'$\mathrm{Log}_{10}(\Gamma / H_0)$', fontsize =15)
ax.set_ylabel(r'$\mathrm{Log}_{10}(\varepsilon)$', fontsize = 15)

ax.set_title(r'$S_8^{\mathrm{full}}-S_8^{\mathrm{approx}}$', fontsize = 20)

#ax.tick_params(axis="x", labelsize=15)
#ax.tick_params(axis="y", labelsize=15)

plt.locator_params(axis='x', nbins=3)

plt.yticks(extraticks)
plt.locator_params(axis='y', nbins=3)




im = ax.imshow(data, extent=extent, origin='upper', interpolation='nearest', cmap='viridis')

fig.colorbar(im)
plt.show()


