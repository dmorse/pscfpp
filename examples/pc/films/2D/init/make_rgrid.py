#!/usr/bin/env python

import math
import numpy as np

rgrid = open('preamble','r')
lines = list(rgrid)
rgrid.close()
grid = list(map(int,lines[14].split()))
params = list(map(float,lines[8].split()))

paramfile = open('../param','r')
paramlines = list(paramfile)
paramfile.close()
T = float(paramlines[36].split()[1])

linenum = 15
line = '    {:.9f}    {:.9f}\n'

for y in range(0,grid[1]):
    for x in range(0,grid[0]):
        y_coord = y * params[1] / grid[1]
        x_coord = x * params[0] / grid[0]
        if y_coord < (T/2) or y_coord > (params[1] - (T/2)):
            lines.append(line.format(0,0))
        elif x_coord < 0.4 or (1 < x_coord and x_coord < 1.4): 
            lines.append(line.format(0.8,0.2))
        else:
            lines.append(line.format(0.2,0.8))
        linenum += 1

rgrid = open('c.rf','w')
rgrid.writelines(lines)
rgrid.close()
