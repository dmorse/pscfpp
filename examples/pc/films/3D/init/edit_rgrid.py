#!/usr/bin/env python

import math

rgrid = open('bcc.rgrid','r')
lines = list(rgrid)
rgrid.close()

paramfile = open('../param','r')
paramlines = list(paramfile)
paramfile.close()
params = [float(lines[8].strip())]
grid = list(map(int,paramlines[25].split()[1:])) 

normalVecId = int(paramlines[34].split()[1])
t = float(paramlines[35].split()[1])
T = float(paramlines[36].split()[1])
params.append(params[0] * 2 + T)
L = params[1]

lines_new = []
for i in range(0,15):
    lines_new.append(lines[i])
lines_new[4]  = "          tetragonal\n"
lines_new[6]  = "                   2\n"
lines_new[8]  = "    2.0000000000E+00      {:.10f}E+00\n".format(L)
lines_new[10] = "           P_4%m_m_m\n"
lines_new[14] = "                  {}          {}          {}\n".format(grid[0],grid[1],grid[2])

linenum = 15

for z in range(0,grid[2]):
    for y in range(0,grid[1]):
        for x in range(0,grid[0]):
            z_coord = z * L / grid[2]
            if z_coord < (T/2) or z_coord > (L - (T/2)):
                lines_new.append('    {:.9f}    {:.9f}\n'.format(0,0))
            else:
                lines_new.append(lines[linenum])
                linenum += 1

            if linenum >= len(lines):
                linenum = 15

rgrid = open('c.rf','w')
rgrid.writelines(lines_new)
rgrid.close()
