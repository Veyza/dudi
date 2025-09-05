import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
fnames = ('./results/narrow_jet_shallow_sd', './results/diffuse_source_steep_sd', './results/diffuse_source_shallow_sd', './results/narrow_jet_steep_sd')
ead = (' narrow jet', ' diffuse source')
sd = ('steep size distribution,', 'shallow size distribution,')
d = []
rad = []
for i in range(4):
	d.append(np.loadtxt(fnames[i]+'.dat', usecols=range(2)))
	rad.append(d[i][:,0])
plt.yscale('log')
plt.xlim(0,110)
plt.plot(rad[0], d[0][:,1], 'k--', label = sd[1] + ead[0])
plt.plot(rad[2], d[2][:,1], 'b--', label = sd[1] + ead[1])
plt.plot(rad[3], d[3][:,1], 'k-', label = sd[0] + ead[0])
plt.plot(rad[1], d[1][:,1], 'b-', label = sd[0] + ead[1])
plt.ylabel('deposited mass, $kg/m^2/s$')
plt.xlabel('distance from the source')
plt.legend(loc = 'upper right')
imname = './results/mass_deposition' + '.png'
plt.savefig(imname)