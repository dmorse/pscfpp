#!/usr/bin/env python

# Imports
import math
import numpy as np
from pscfpp.param import Composite

# Number of A/B interfaces in unit cell
n_interfaces = 2

# Get data from "preamble" file
rgrid = open('preamble','r')
lines = list(rgrid)
rgrid.close()
grid_new = list(map(int,lines[14].split())) 
params = list(map(float,lines[8].split()))

# Get data from param file
paramfile = Composite("../param")
ifg = paramfile.AmIteratorBasis.ImposedFieldsGenerator
normalVecId = ifg.normalVecId
t = ifg.interfaceThickness
T = ifg.excludedThickness
L = params[normalVecId]

# Setup to loop through gridpoints
linenum = 15
line = '    {:.9f}    {:.9f}\n'

# Loop through gridpoints
for y in range(0,int(grid_new[1])):
    for x in range(0,grid_new[0]):
        
        # Get x and y coords at this gridpoint
        x_coord = x * params[0] / grid_new[0]
        y_coord = y * params[1] / grid_new[1]

        # Calculate wall and polymer volume fractions
        rhoW = 0.5*(1+np.tanh(4*(((.5*(T-L))+np.abs(y_coord-(L/2)))/t)));
        rho = 1-rhoW

        # Estimate A and B compositions using a cosine
        rhoB = 0.4 * np.cos(x_coord * n_interfaces * np.pi / params[0]) + 0.5
        rhoA = 1-rhoB

        # Add data to lines list
        lines.append(line.format(rhoA*rho,rhoB*rho))
        linenum += 1

# Write lines list to file
rgrid = open('c.rf','w')
rgrid.writelines(lines)
rgrid.close()
