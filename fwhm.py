import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--file", help="fits file name")
args = parser.parse_args()


from astropy.io import fits
from astropy.wcs import WCS
from astropy.time import Time
from astropy.time import TimeDelta

import math
import numpy as np

hdulist = fits.open(args.file)# load file with the name stored in sys.argv[1] (first command line argumrnt)
fd = hdulist[0].data# access to pixels data
t=Time(hdulist[0].header['DATE-OBS'])
et=TimeDelta(float(hdulist[0].header['EXPTIME'])/2.0,format='sec')
t+=et


w = WCS(hdulist[0].header)# access to WCS parameters
hdulist.close()


pix_tl=np.array([[0,0]], np.float_)# pixel position of the top left corner
pix_br=np.array([[fd.shape[1],fd.shape[0]]], np.float_)# pixel position of the bottom right
pix_tr=np.array([[fd.shape[1],0]], np.float_)# pixel position of the top right
cp_tl = w.wcs_pix2world(pix_tl,1)# celestial position of the top left corner
cp_br = w.wcs_pix2world(pix_br,1)# celestial position of the bottom right
cp_tr = w.wcs_pix2world(pix_tr,1)# celestial position of the bottom right

def angular_distance(point1,point2):
    ra=math.pi*point1[0][0]/180
    dec = math.pi * point1[0][1] / 180
    RA = math.pi * point2[0][0] / 180
    DEC = math.pi * point2[0][1] / 180
    x=math.cos(ra)*math.cos(dec)
    y=math.sin(ra)*math.cos(dec)
    z=math.sin(dec)
    X = math.cos(RA) * math.cos(DEC)
    Y = math.sin(RA) * math.cos(DEC)
    Z = math.sin(DEC)
    return math.acos(x*X+y*Y+z*Z)

scalex=648000*angular_distance(cp_tl,cp_tr)/math.pi/fd.shape[1]
scaley=648000*angular_distance(cp_br,cp_tr)/math.pi/fd.shape[0]

pix_center = np.array([[fd.shape[0]/2,fd.shape[1]/2]], np.float_)

cp_center = w.wcs_pix2world(pix_center,1)# celestial position of the center of the image
center_RA = cp_center[0,0]
center_Dec = cp_center[0,1]
fieldsize = 10800*angular_distance(cp_tl,cp_br)/math.pi


#form get-request to CDS to obtain Gaia DR1 data
Ugaia = ('http://vizier.u-strasbg.fr/viz-bin/asu-tsv'
     '?-source=I/337/gaia&-c=%f %f&-c.geom=b&-c.u=arcmin&-c.r=%d&-sort=_r&-order=I'
         '&-out.all&-out.max=unlimited&-out.form=| -Separated-Values'  % (center_RA,center_Dec, fieldsize))
import requests
#send get-request and load reply into the list rgaia
rgaia = requests.get(Ugaia).text.split('\n')

#convert rgaia to arrays
def gaiadata(gdata):#converts a list of stars from CDS to numpy arrays
    eqpos = np.empty((0,2), np.float_)#array of pairs of coordinates RA,Dec
    gmag = np.empty(0, np.float_)#array of magnitudes
    for k in range(len(gdata)):
        if (not (gdata[k].__contains__('#')
                 or gdata[k].__contains__('_r')
                 or gdata[k].__contains__('-|')
                 or gdata[k].__contains__('arcmin')
                 or len(gdata[k]) <= 1)
            ):
            # print rgaia[k]
            line = gdata[k].split('|')
            #equatorial coordinates in deg
            eqpos = np.append(eqpos, [[float(line[4]), float(line[6])]], axis=0)
            #G-magnitudes
            gmag = np.append(gmag, [float(line[55])], axis=0)

    return eqpos,gmag
eqpos,mags = gaiadata(rgaia)

#convert equatorial coordinates to pixel positions
pix = w.wcs_world2pix(eqpos, 1)
Box=40
magb = 10.0
magf = 14.0

def fwhm(J):
    box = J.shape[0]
    levelJ = np.median(J) + (np.max(J)-np.median(J))/2
    S=0.0
    for y in range(box):
        for x in range(box):
            if(J[y,x]>levelJ):
                S+=1
    return 2*math.sqrt(S/math.pi)

b = np.empty(0, np.float_)#array of fwhm estimates
for k in range(pix.shape[0]):
    if((Box <pix[k,0]< fd.shape[1]-Box) and (Box <pix[k,1]< fd.shape[0]-Box)\
       and (mags[k]>magb) and (mags[k]< magf)):

        J = fd[int(pix[k, 1]) - Box / 2:int(pix[k, 1]) + Box / 2,\
            int(pix[k, 0]) - Box / 2:int(pix[k, 0]) + Box / 2]
        b = np.append(b, [scalex*fwhm(J)], axis=0)

print (str(t)+('|%14.7f|%7.3f' % (float(t.mjd), np.median(b))))
