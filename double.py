import argparse
import math
import numpy as np
import shapelet
from shapelet import *

from mpl_toolkits.mplot3d.axes3d import Axes3D
from matplotlib import cm
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser()
parser.add_argument("--box", help="fits file name")
parser.add_argument("--s_x", help="size of the box to calculate the star centroids")
parser.add_argument("--s_y", help="magnitude limit for bright end")
parser.add_argument("--phi", help="magnitude limit for faint end")
parser.add_argument("--noise", help=(""))
args = parser.parse_args()

dmag = 1.0
rho = 0.0
alpha = 0

box=int(args.box)
s_x = float(args.s_x)
s_y = float(args.s_y)
phi = float(args.phi)
noise = float(args.noise)

xc = box/2
yc = box/2

def starmodel(box,s_x,s_y,phi,noise,xc,yc,dmag,rho,alpha):
    bkgr = 0
    Jmax = 1
    j1 = stimg_gaussian(box, bkgr, Jmax, s_x, s_y, phi, xc, yc, noise)
    Jmax = Jmax + dmag
    xc = xc + rho*np.cos(alpha)
    yc = yc + rho*np.sin(alpha)
    j2 = stimg_gaussian(box, bkgr, Jmax, s_x, s_y, phi, xc, yc, noise)
    j = j1 + j2
    return j

def ellipt(xc,yc,box,j):
    qxx_1 =0; qxy_1 = 0; qyy_1 = 0; q2 = 0
    for x in range(box):
        for y in range(box):
            qxx_1 = qxx_1 + j[y,x]*(x-xc)**2
            qyy_1 = qyy_1 + j[y,x]*(y-yc)**2
            qxy_1 = qxy_1 + j[y,x]*(x-xc)*(y-yc)
            q2 = q2 + j[x,y]
    qxx = qxx_1 / q2; qyy = qyy_1 / q2; qxy = qxy_1 / q2
    e1 = (qxx - qyy) / (qxx + qyy)
    e2 = 2.0 * qxy / (qxx + qyy)
    e = math.sqrt(e1**2+e2**2)
    return math.sqrt(qxx), math.sqrt(qyy), e


def ellipt0(s_x,s_y):
    qxx = s_x**2; qyy = s_y**2; qxy = 0
    e1 = (qxx - qyy) / (qxx + qyy)
    e2 = 2.0 * qxy / (qxx + qyy)
    e = math.sqrt(e1 ** 2 + e2 ** 2)
    return math.sqrt(qxx), math.sqrt(qyy), e

#print ellipt(xc,yc,box)
#print ellipt0(s_x,s_y)

j = starmodel(box,s_x,s_y,phi,noise,xc,yc,dmag,rho,alpha)

fig = plt.figure()
ax = fig.gca(projection='3d')
X = np.arange(0, box, 1)
Y = np.arange(0, box, 1)
X, Y = np.meshgrid(X, Y)
ax.plot_surface(X, Y, j,rstride=1, cstride=1, cmap=cm.Greys,\
                       linewidth=0, antialiased=False)
plt.xlabel('x')
plt.ylabel('y')
plt.savefig('model.png')
plt.show()

def setka(s_x,s_y):
    j = starmodel(box, s_x, s_y, phi, noise, xc, yc, dmag, rho, alpha)
    s = ellipt(xc, yc, box,j)
    s0 = ellipt0(s_x, s_y)
    sx.append(s0[0])
    Sx.append(s[0])
    sy.append(s0[1])
    Sy.append(s[1])
    E.append(s[2]); e.append(s0[2])
    return s,s0

s_x = 0.5
sx = []; Sx = []; sy = []; Sy = []; E = []; e = []
for i in range(10):
    s_x = s_x + 1
    s_y = 1.5 * s_x
    setka(s_x,s_y)

print e

plt.figure(1)
#plt.axis([0,10,0,10])
plt.plot(sx,Sx,'kp')
plt.xlabel('$s_x$', fontsize=30)
plt.ylabel('$S_x$', fontsize=30)
plt.figure(2)
plt.plot(sy,Sy,'kp')
#plt.axis([0,10,0,10])
plt.xlabel('$s_y$', fontsize=30)
plt.ylabel('$S_y$', fontsize=30)
plt.show()