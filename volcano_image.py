import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
plt.gray()
moments = [200, 400, 600, 800, 1000, 1200, 1400, 1600, 1800]
for i in range(9):
	fname = './results/' + str(i+1) + '.dat'
	d = np.loadtxt(fname, usecols=range(128))
	d = np.rot90(d)
	d = np.where(d > 1e-20, d, 1e-19)
	d = np.log10(d)
	plt.figure(i+1)
	imgplot = plt.imshow(d, vmin=-17, vmax=-4)
	plt.colorbar()
	plt.text(75, 117, ('time = ' + str(moments[i]) + ' s'), color = 'w')
	imname = './results/volcano_' + str(i+1) + '.png'
	plt.savefig(imname)