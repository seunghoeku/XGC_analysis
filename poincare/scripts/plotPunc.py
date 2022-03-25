import os
import math
# git@github.com:seunghoeku/XGC_reader.git
#import xgc_reader
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as tri

import random


#fin = open('../../../build/bumm.punc.theta_psi.txt', 'r')
#fin = open('./bumm.punc.theta_psi.txt', 'r')
fin = open('./bumm.TP.txt', 'r')

XT = []
YT = []
ID = 0

lines = fin.readlines()
xt = []
yt = []
for l in lines[1:] :
    x = l.strip(' \n\t').split(',')
    id = int(x[0])
    xc = float(x[1])
    yc = float(x[2])
    if (id != ID) :
        ID = id
        XT.append(xt)
        YT.append(yt)
        xt = []
        yt = []
    xt.append(xc)
    yt.append(yc)

XT.append(xt)
YT.append(yt)

fig = plt.figure(figsize=[8,12])
N = len(XT)
for i in range(N) :
    x = XT[i]
    y = YT[i]
    plt.scatter(x, y, s=1, marker='x', edgecolor='none')

plt.title('VTKm')
plt.show()
