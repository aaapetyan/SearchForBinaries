feb = open('ellipticity-bin.txt','r')
reb = feb.readlines()
feb.close()

import numpy as np
import argparse
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser()
parser.add_argument("--filter", help="filter as -i- or -r-")
parser.add_argument("--out", help="all or most interesting")
args = parser.parse_args()

import math

eter = []
res = open('list_all.txt', 'a')

for line in reb:
    if(line.__contains__('tar') and (line.__contains__(args.filter))):
        #print(line)
        if(line.__contains__('-u-')):
            sdssf = 'u|'
        if (line.__contains__('-g-')):
            sdssf = 'g|'
        if (line.__contains__('-r-')):
            sdssf = 'r|'
        if (line.__contains__('-i-')):
            sdssf = 'i|'
        if (line.__contains__('-z-')):
            sdssf = 'z|'

        p = line.split('|')
        et = float(p[10])
        er = float(p[15])
        dev_er = float(p[16])
        if((not math.isnan(er)) and (not math.isnan(dev_er))):
            if(0<dev_er<0.3):
                lspmid = line[0:12].replace('|t',' |')
                dline = line[16:78]
                if(dline[0]==' '):
                    dline = dline.lstrip(' ')+'|'
                if(args.out=='all'):
                    res.write(lspmid + dline + sdssf \
                          + '%7.3f|%7.3f|%7.3f|%7.3f|%7.3f|' % \
                          (et, er, dev_er, et - er, abs(et - er) / dev_er)+'\n')
                else:
                    if((et<1) and (et-er>0)and(abs(et-er)/dev_er>3)):
                        eter.append(et-er)
                        print(lspmid+ dline + sdssf \
                             +'%7.3f|%7.3f|%7.3f|%7.3f|%7.3f|' % \
                              (et,er,dev_er,et-er,abs(et-er)/dev_er))

print eter