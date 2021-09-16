# add argparse module
import argparse

import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d.axes3d import Axes3D
from matplotlib.colors import SymLogNorm

# analyse command-line parameters
parser = argparse.ArgumentParser()
parser.add_argument("--file", help="fits file name")
parser.add_argument("--lspmfile", help="path to lspm_n.dat")
parser.add_argument("--id", help="LSPM ID of the star")
parser.add_argument("--box", help="size of the box")
parser.add_argument("--filter", help="sdss filter")
args = parser.parse_args()

# load LSPM catalogue file

flspm = open(args.lspmfile, 'r')
rlspm = flspm.readlines()
flspm.close()

# search for star with --id
for line in rlspm:
    if (line.__contains__(args.id)):
        # print(line)
        lspmname = line[5:16]  # star name - we need it
        lhsname = line[17:22]
        nlttname = line[23:28]

        RAlspm = float(line[127:137])
        Declspm = float(line[138:148])
        pmra = float(line[181:188])  # proper motions - we need it
        pmdec = float(line[189:196])

        Jmag = float(line[227:232])
        Hmag = float(line[233:238])
        Kmag = float(line[239:244])
        Vemag = float(line[245:251])  # visual magnitude - we need it
        VJ = float(line[252:257])

# add astropy modules to operate fits-file and wcs
from astropy.io import fits
from astropy.wcs import WCS

hdulist = fits.open(args.file)  # load fits-file
fd = hdulist[0].data  # access to pixels data
h, w = fd.shape  # size of the image in pixels

W = WCS(hdulist[0].header)  # access to WCS parameters

# add astropy modules to operate time representations
from astropy.time import Time
from astropy.time import TimeDelta

# calculation of the central time of exposure
t = Time(hdulist[0].header['DATE-OBS'])
et = TimeDelta(float(hdulist[0].header['EXPTIME']) / 2.0, format='sec')
t += et

hdulist.close()  # we stop use hdu list of fits-file, all data stored in W and fd

# now we know central exposure time and can calculate position of LSPM-star at epoch of observation
# we consider proper motion component as tangential position
# and we should convert tangential coordinates into equatorial ones

# add math
import math


def RADecFromTang(ksi, eta, RA, Dec):
    secRho = math.sqrt(1 + ksi * ksi + eta * eta)
    t11 = -math.sin(RA)
    t12 = math.cos(RA)
    t13 = 0
    t21 = -math.cos(RA) * math.sin(Dec)
    t22 = -math.sin(RA) * math.sin(Dec)
    t23 = math.cos(Dec)
    t31 = math.cos(RA) * math.cos(Dec)
    t32 = math.sin(RA) * math.cos(Dec)
    t33 = math.sin(Dec)
    x = (ksi * t11 + eta * t21 + t31) / secRho
    y = (ksi * t12 + eta * t22 + t32) / secRho
    z = (ksi * t13 + eta * t23 + t33) / secRho
    ra = math.atan2(y, x)
    dec = math.atan2(z, math.sqrt(x * x + y * y))
    if (ra < 0):
        ra += 2 * math.pi
    return ra, dec


t.format = 'decimalyear'
# print(t.value)                            # t.value - epoch - we need it


RA, Dec = RADecFromTang((t.value - 2000.0) * math.pi * pmra / 648000.0,
                        (t.value - 2000.0) * math.pi * pmdec / 648000.0,
                        math.radians(RAlspm), math.radians(Declspm))
# add numpy


import numpy as np

pix_tl = np.array([[0, 0]], np.float_)  # pixel position of the top left corner
pix_br = np.array([[w, h]], np.float_)  # pixel position of the bottom right
pix_tr = np.array([[w, 0]], np.float_)  # pixel position of the top right
cp_tl = W.wcs_pix2world(pix_tl, 1)  # celestial position of the top left corner
cp_br = W.wcs_pix2world(pix_br, 1)  # celestial position of the bottom right
cp_tr = W.wcs_pix2world(pix_tr, 1)  # celestial position of the bottom right


# next we want to load Gaia DR1 data
# we should calculate equatorial coordinate of the center of the image and
# size of the image in arcmins

def angular_distance(point1, point2):  # function to calculate angular distance between points
    ra1 = math.radians(point1[0][0])
    dec1 = math.radians(point1[0][1])
    ra2 = math.radians(point2[0][0])
    dec2 = math.radians(point2[0][1])
    x1 = math.cos(ra1) * math.cos(dec1)
    y1 = math.sin(ra1) * math.cos(dec1)
    z1 = math.sin(dec1)
    x2 = math.cos(ra2) * math.cos(dec2)
    y2 = math.sin(ra2) * math.cos(dec2)
    z2 = math.sin(dec2)
    return math.degrees(math.acos(x1 * x2 + y1 * y2 + z1 * z2))  # result in degrees


# calculate scales in arcsec/pix
scalex = 3600.0 * angular_distance(cp_tl, cp_tr) / w
scaley = 3600.0 * angular_distance(cp_br, cp_tr) / h

# celestial position of the center of the image
pix_center = np.array([[w / 2, h / 2]], np.float_)
cp_center = W.wcs_pix2world(pix_center, 1)

# field of view
fov = 60.0 * angular_distance(cp_tl, cp_br)

# request line for Gaia DR1
Ugaia = ('http://vizier.u-strasbg.fr/viz-bin/asu-tsv'
         '?-source=I/337/gaia&-c=%f %f&-c.geom=b&-c.u=arcmin&-c.r=%d&-sort=_r&-order=I'
         '&-out.all&-out.max=unlimited&-out.form=| -Separated-Values'
         % (cp_center[0, 0], cp_center[0, 1], fov))

# add requests module
import requests

# send get-request and load reply into the list rgaia
rgaia = requests.get(Ugaia).text.split('\n')


# function to convert rgaia to arrays
def gaiadata(gdata):  # converts a list of stars from CDS to numpy arrays
    eqpos = np.empty((0, 2), np.float_)  # array of pairs of coordinates RA,Dec
    gmag = np.empty(0, np.float_)  # array of magnitudes
    for k in range(len(gdata)):
        if (not (gdata[k].__contains__('#')
                 or gdata[k].__contains__('_r')
                 or gdata[k].__contains__('-|')
                 or gdata[k].__contains__('arcmin')
                 or len(gdata[k]) <= 1)
            ):
            # print rgaia[k]
            line = gdata[k].split('|')
            # equatorial coordinates in deg
            eqpos = np.append(eqpos, [[float(line[4]), float(line[6])]], axis=0)
            # G-magnitudes
            gmag = np.append(gmag, [float(line[55])], axis=0)

    return eqpos, gmag


# eqpos - equatorial coordinates, mags - magnitudes
eqpos, mags = gaiadata(rgaia)

# convert equatorial coordinates to pixel positions
pix = W.wcs_world2pix(eqpos, 1)

pix_lspm = W.wcs_world2pix(np.array([[math.degrees(RA), math.degrees(Dec)]], np.float_), 1)

# set the box size
box = int(args.box)
filter = str(args.filter)


# function to calculate ellipticity
def ellipt(xc, yc, box, J):
    qxx_1 = 0; qxy_1 = 0; qyy_1 = 0; q2 = 0
    for x in range(box):
        for y in range(box):
            qxx_1 += J[y, x] * (x - xc) ** 2
            qyy_1 += J[y, x] * (y - yc) ** 2
            qxy_1 += J[y, x] * (x - xc) * (y - yc)
            q2 += J[y, x]
    qxx = qxx_1 / q2; qyy = qyy_1 / q2; qxy = qxy_1 / q2;
    e1 = (qxx - qyy) / (qxx + qyy)
    e2 = 2.0 * qxy / (qxx + qyy)
    e = math.sqrt(e1 ** 2 + e2 ** 2)
    return e


def photocenter(J, bgJ):
    xc = 0.0; yc = 0.0; F = 0.0
    for x in range(J.shape[0]):
        for y in range(J.shape[1]):
            xc += (J[y, x] - bgJ) * x
            yc += (J[y, x] - bgJ) * y
            F += (J[y, x] - bgJ)
    if (F != 0):
        xc = xc / F
        yc = yc / F
    else:
        xc = J.shape[0] / 2
        yc = J.shape[1] / 2
    return xc, yc


Xt = pix_lspm[0, 0]
Yt = pix_lspm[0, 1]
print(Xt, Yt)

Xt = int(Xt)
Yt = int(Yt)

# if lspm-star is in the frame then calculate ellipticity
if (box < Xt < w - box and box < Yt < h - box):
    J = fd[Yt - box / 2:Yt + box / 2, Xt - box / 2:Xt + box / 2]
    Xc = photocenter(J, 0)[0]
    Yc = photocenter(J, 0)[1]
    Et = ellipt(Xc, Yc, box, J)

    fig = plt.figure()
    Jmin = J.min()
    Jmax = J.max()

    plt.imshow(J, origin='lower', cmap=cm.gray_r, interpolation='none',\
    norm=SymLogNorm(linthresh=(Jmax-Jmin)/100,vmin=Jmin,vmax=Jmax))
    plt.xlabel('x')
    plt.ylabel('y')
    plt.savefig('/home/apetyan/PycharmProjects/images/ourstars/' + filter + '/' + lspmname.rstrip() + '_%f.eps' % (Et))
    plt.close()
    fig.clf()