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

for y in range(0,grid[1]):
    for x in range(0,grid[0]):
        x_coord = x * params[0] / grid_new[0]
        y_coord = y * params[1] / grid_new[1]
        rhoW = 0.5*(1+np.tanh(4*(((.5*(T-L))+np.abs(y_coord-(L/2)))/t)));
        rho = 1-rhoW
        rhoB = 0.4 * np.cos(x_coord * n_interfaces * np.pi / params[0]) + 0.5
        rhoA = 1-rhoB
        lines.append(line.format(rhoA*rho,rhoB*rho))
        linenum += 1

rgrid = open('c.rf','w')
rgrid.writelines(lines)
rgrid.close()
