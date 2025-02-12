#!/usr/bin/env python

# Imports
import math
from pscfpp.param import Composite

# Read data for bcc in the bulk
rgrid = open('bcc.rf','r')
lines = list(rgrid)
rgrid.close()
params = [float(lines[8].strip())]

# Get data from param file
paramfile = Composite("../param")
ifg = paramfile.AmIteratorBasis.ImposedFieldsGenerator
normalVecId = ifg.normalVecId
t = ifg.interfaceThickness
T = ifg.excludedThickness
grid = paramfile.Domain.mesh

# add lattice parameter (going from cubic -> tetragonal)
params.append(params[0] * 2 + T)
L = params[-1]

# Write header for rgrid file
lines_new = []
for i in range(0,15):
    lines_new.append(lines[i])
lines_new[4]  = "          tetragonal\n"
lines_new[6]  = "                   2\n"
lines_new[8]  = "    {:.10e}      {:.10e}\n".format(params[0], params[1])
lines_new[10] = "           P_4%m_m_m\n"
lines_new[14] = "                  {}          {}          {}\n".format(grid[0],grid[1],grid[2])

linenum = 15

# Loop through gridpoints
for z in range(0,grid[2]):
    for y in range(0,grid[1]):
        for x in range(0,grid[0]):
            # Calculate z coordinate
            z_coord = z * L / grid[2]

            # Initialize compositions to zero if in the wall,
            # otherwise use data from bcc in the bulk
            if z_coord < (T/2) or z_coord > (L - (T/2)):
                lines_new.append('    {:.9f}    {:.9f}\n'.format(0,0))
            else:
                lines_new.append(lines[linenum])
                linenum += 1

            # If we have written one full unit cell of bulk data,
            # reset linenum back to the start of the bulk unit cell
            if linenum >= len(lines):
                linenum = 15

# Write data to file
rgrid = open('c.rf','w')
rgrid.writelines(lines_new)
rgrid.close()
