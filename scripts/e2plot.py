import numpy as np
import matplotlib.pyplot as plt

d = np.loadtxt("./results/E2_profile.dat", usecols=range(2))
t = d[:, 0]
dens = d[:, 1]
d = np.loadtxt("./input_data_files/E2_1.6.txt", usecols=range(2))
hrdt = d[:, 0]
hrddens = d[:, 1]
plt.figure(1)
plt.ylabel("number density of grains $> 1.6\\ \\mu m$")
plt.xlabel("seconds from the moment of the closest approach")
plt.ylim(0, 0.12)
plt.suptitle("HRD number density profile of E2 flyby")
model = plt.plot(t, dens, "k-", label="model number density")
hrd = plt.plot(hrdt, hrddens, "bo", label="HRD measurements")
plt.legend(loc="upper left")
plt.show()
