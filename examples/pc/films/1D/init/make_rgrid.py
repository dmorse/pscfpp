#!/usr/bin/env python

import math
import numpy as np

n_interfaces = 2

rgrid = open('preamble','r')
lines = list(rgrid)
rgrid.close()
grid = list(map(int,lines[14].split())) 
params = list(map(float,lines[8].split()))

paramfile = open('../param','r')
paramlines = list(paramfile)
paramfile.close()
normalVecId = int(paramlines[34].split()[1])
t = float(paramlines[35].split()[1])
T = float(paramlines[36].split()[1])
L = params[normalVecId]

linenum = 15
line = '    {:.9f}    {:.9f}\n'

for x in range(0,grid[0]):
    x_coord = x * params[0] / grid[0]
    rhoW = 0.5*(1+np.tanh(4*(((.5*(T-L))+np.abs(x_coord-(L/2)))/t)));
    rho = 1-rhoW
    if x_coord < (T/2):
        lines.append(line.format(rho*0.1,rho*0.9))
    elif x_coord > (params[0] - (T/2)):
        if n_interfaces % 2 == 1: # n_interfaces is odd
            lines.append(line.format(rho*0.9,rho*0.1))
        else:
            lines.append(line.format(rho*0.1,rho*0.9)) 
    else:
        x = x_coord - (T/2)
        rhoB = 0.4 * np.cos(x * n_interfaces * np.pi / (L-T)) + 0.5
        rhoA = 1-rhoB
        lines.append(line.format(rhoA*rho,rhoB*rho))
    linenum += 1

rgrid = open('c.rf','w')
rgrid.writelines(lines)
rgrid.close()
