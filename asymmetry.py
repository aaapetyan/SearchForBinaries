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
parser.add_argument("--box", help="size of the box to calculate the star centroids")
parser.add_argument("--s_x", help="magnitude limit for bright end")
parser.add_argument("--s_y", help="magnitude limit for faint end")
parser.add_argument("--phi", help="")
parser.add_argument("--noise", help=(""))

args = parser.parse_args()

box = int(args.box)
s_x = float(args.s_x)
s_y = float(args.s_y)
phi = math.radians(float(args.phi))
noise = float(args.noise)

def photocenter(J,bgJ):
    xc=0.0
    yc=0.0
    F=0.0
    for x in range(J.shape[0]):
        for y in range(J.shape[1]):
            xc += (J[y, x]-bgJ)*x
            yc += (J[y, x]-bgJ)*y
            F += (J[y, x]-bgJ)
    if(F!=0):
        xc = xc/F
        yc = yc/F
    else:
        xc=J.shape[0]/2
        yc=J.shape[1]/2
    return xc,yc

def starmodel(box, s_x, s_y, phi, noise, xc, yc, dmag, rho):
    bkgr = 0
    jmax = 1
    j1 = stimg_gaussian(box, bkgr, jmax, s_x, s_y, 0, xc, yc, noise)
    jmax = 2.512**(-dmag)
    xc2 = xc + rho*np.cos(phi)
    yc2 = yc + rho*np.sin(phi)
    j2 = stimg_gaussian(box, bkgr, jmax, s_x, s_y, 0, xc2, yc2, noise)
#    j2 = 0
    j = j1 + j2
    return j

def asymmetry(box,J):
    sumJ = 0; JJ = 0
    for x in range(box):
        for y in range(box):
            sumJ += math.fabs(J[y,x] - J[box-y-1,box-x-1])
            JJ += J[y, x]
    A = sumJ / JJ
    return A

dmag = 0.5
rho = 3
J = starmodel(box, s_x, s_y, phi, noise, box / 2, box / 2, dmag, rho)

print(asymmetry(box,J))