#!/usr/bin/env python

import argparse
import math
import numpy as np
import shapelet
from shapelet import *

from mpl_toolkits.mplot3d.axes3d import Axes3D
from matplotlib import cm
import matplotlib.pyplot as plt


parser = argparse.ArgumentParser()
parser.add_argument("--file", help="fits file name")
parser.add_argument("--box", help="size of the box to calculate the star centroids")
parser.add_argument("--s_x", help="magnitude limit for the bright end")
parser.add_argument("--s_y", help="magnitude limit for the faint end")
parser.add_argument("--phi", help="")
parser.add_argument("--noise", help=(""))
parser.add_argument("--mode", help=("noise, rho or dmag for ellipt plots"))
args = parser.parse_args()

dmag = 1.0
rho = 0.0

box = int(args.box)
s_x = float(args.s_x)
s_y = float(args.s_y)
phi = math.radians(float(args.phi))
noise = float(args.noise)
mode = str(args.mode)

xc = box/2
yc = box/2

def starmodel(box, s_x, s_y, phi, noise, xc, yc, dmag, rho):
    bkgr = 0
    jmax = 1
    j1 = stimg_gaussian(box, bkgr, jmax, s_x, s_y, phi, xc, yc, noise)
    jmax = jmax + dmag
    xc = xc + rho*np.cos(phi)
    yc = yc + rho*np.sin(phi)
    j2 = stimg_gaussian(box, bkgr, jmax, s_x, s_y, phi, xc, yc, noise)
    j = j1 + j2
    return j

#fig = plt.figure()
#ax = fig.gca(projection='3d')
#X = np.arange(0, box, 1)
#Y = np.arange(0, box, 1)
#X, Y = np.meshgrid(X, Y)
#ax.plot_surface(X, Y, J,rstride=1, cstride=1, cmap=cm.Greys,\
#                       linewidth=0, antialiased=False)
#plt.xlabel('x')
#plt.ylabel('y')
#plt.savefig('model.png')
#plt.show()

def ellipt(xc,yc,box):
    qxx_1 =0; qxy_1 = 0; qyy_1 = 0; q2 = 0
    for x in range(box):
        for y in range(box):
            qxx_1 += J[x,y]*(x-xc)**2
            qyy_1 += J[x,y]*(y-yc)**2
            qxy_1 += J[x,y]*(x-xc)*(y-yc)
            q2 += J[x,y]
    qxx = qxx_1 / q2; qyy = qyy_1 / q2; qxy = qxy_1 / q2;
    e1 = (qxx - qyy) / (qxx + qyy)
    e2 = 2.0 * qxy / (qxx + qyy)
    e = math.sqrt(e1 ** 2 + e2 ** 2)
    return e

def ellipt0(s_x, s_y):
    qxx = s_x**2; qyy = s_y**2; qxy = 0
    e1 = (qxx - qyy) / (qxx + qyy)
    e2 = 2.0 * qxy / (qxx + qyy)
    e = math.sqrt(e1 ** 2 + e2 ** 2)
    return e

ee = []
e0 = []

if mode=="rho":
    rhos = []
    for i in range(35):
        J = starmodel(box, s_x, s_y, phi, noise, xc, yc, dmag, rho)
        rhos.append(rho)
        ee.append(ellipt(xc, yc, box))
        e0.append(ellipt0(s_x, s_y))
        rho += 0.1

    print 'e0 =', e0
    print 'e =', ee

    plt.figure(1)
    plt.plot(rhos, ee, 'kp')
    #plt.axis([0,0.02,0,0.1])
    plt.xlabel('$ rho $', fontsize=30)
    plt.ylabel('$e_{calc}$', fontsize=30)
    plt.show()

elif mode=="dmag":
    dmags = []
    for i in range(35):
        J = starmodel(box, s_x, s_y, phi, noise, xc, yc, dmag, rho)
        dmags.append(dmag)
        ee.append(ellipt(xc, yc, box))
        e0.append(ellipt0(s_x, s_y))
        dmag += 0.1

    print 'e0 =', e0
    print 'e =', ee

    plt.figure(1)
    plt.plot(dmags, ee, 'kp')
    # plt.axis([0,0.02,0,0.1])
    plt.xlabel('$ dmag $', fontsize=30)
    plt.ylabel('$e_{calc}$', fontsize=30)
    plt.show()

elif mode=="noise":
    ns = []

    for i in range(20):
        J = starmodel(box, s_x, s_y, phi, noise, xc, yc, dmag, rho)
        ns.append(noise)
        ee.append(ellipt(xc, yc, box))
        e0.append(ellipt0(s_x, s_y))
        noise += 0.001

    print 'e0 =', e0
    print 'e =', ee

    plt.figure(1)
    plt.plot(ns, ee, 'kp')
    plt.axis([0, 0.02, 0, 0.1])
    plt.xlabel('$S/N$', fontsize=30)
    plt.ylabel('$e_{calc}$', fontsize=30)
    plt.show()

else:
    print "choose one of the following modes: dmag, rho or noise"
