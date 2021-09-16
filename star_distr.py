import math
import numpy as np
from numpy import histogram
import matplotlib.pyplot as plt

with open("sdss-stars.txt") as f:
    lines = f.readlines()
    vmag = [float(line.split('|')[12]) for line in lines]
    j = [float(line.split('|')[9]) for line in lines]
    h = [float(line.split('|')[10]) for line in lines]
    ks = [float(line.split('|')[11]) for line in lines]
    v_j = [float(line.split('|')[13]) for line in lines]
    pmra = [float(line.split('|')[7]) for line in lines]
    pmdec = [float(line.split('|')[8]) for line in lines]

pm = [np.sqrt(pmra[i]**2+pmdec[i]**2) for i in range(len(pmra))]
j_h = [j[i]-h[i] for i in range(len(h))]
h_ks = [h[i]-ks[i] for i in range(len(h))]

#plt.hist(vmag, bins='auto',color='black')  # plt.hist passes it's arguments to np.histogram
#plt.xlabel('mag')
#plt.ylabel('Number of stars')
#plt.title("Apparent magnitude distribution")

#plt.figure(1)
#plt.hist(v_j, bins='auto',color='black')  # plt.hist passes it's arguments to np.histogram
#plt.xlabel('V - J, mag')
#plt.ylabel('Number of stars')
#plt.title("Color index distribution")

plt.figure(2)
plt.plot(h_ks,j_h,'kp')
plt.title('Color-Color Diagram')
plt.xlabel('H - Ks, mag')
plt.ylabel('J - H, mag')

#plt.figure(3)
#plt.hist(pm, bins='auto',color='black')
#plt.xlabel('Proper motion, mas/yr')
#plt.ylabel('Number of stars')
#plt.title("Proper motion distribution")

plt.show()