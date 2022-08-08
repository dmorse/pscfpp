#!/usr/bin/env python

import math
import numpy as np

rgrid = open('preamble','r')
lines = list(rgrid)
rgrid.close()
grid_new = [64] 
params = [2.2]
T = 0.2
t = 0.1

linenum = 15
line = '    {:.9f}    {:.9f}\n'

for x in range(0,int(grid_new[0])):
    x_coord = x * params[0] / grid_new[0]
    if x_coord < (T/2) or x_coord > (params[0] - (T/2)):
        lines.append(line.format(0,0))
    elif x_coord < 0.5 or (1.1 < x_coord and x_coord < 1.5): 
        lines.append(line.format(0.8,0.2))
    else:
        lines.append(line.format(0.2,0.8))
    linenum += 1

rgrid = open('c.rf','w')
rgrid.writelines(lines)
rgrid.close()
