#!/usr/bin/env python

# Imports
import math
from pscfpp.param import Composite

# Number of A/B interfaces between the walls
n_interfaces = 3

# Get data from "preamble" file
rgrid = open('preamble','r')
lines = list(rgrid)
rgrid.close()
grid = list(map(int,lines[14].split())) 
params = list(map(float,lines[8].split()))

# Get data from param file
paramfile = Composite("../param")
ifg = paramfile.AmIteratorBasis.ImposedFieldsGenerator
normalVecId = ifg.normalVecId
t = ifg.interfaceThickness
T = ifg.excludedThickness
L = params[normalVecId]

# Setup to loop over gridpoints
linenum = 15
line = '    {:.9f}    {:.9f}\n'

# Loop over gridpoints
for x in range(0,grid[0]):
    # Calculate x coordinate at this gridpoint
    x_coord = x * params[0] / grid[0]

    # Calculate wall and polymer volume fractions
    rhoW = 0.5*(1+math.tanh(4*(((.5*(T-L))+math.fabs(x_coord-(L/2)))/t)));
    rho = 1-rhoW

    # Use rho and n_interfaces to guess A and B compositions
    if x_coord < (T/2):
        # Initialize points inside of the bottom wall
        lines.append(line.format(rho*0.1,rho*0.9))
    elif x_coord > (params[0] - (T/2)):
        # Initialize top wall (either same or opposite of
        # bottom wall, depending on n_interfaces)
        if n_interfaces % 2 == 1: # n_interfaces is odd
            lines.append(line.format(rho*0.9,rho*0.1))
        else:
            lines.append(line.format(rho*0.1,rho*0.9)) 
    else:
        # Initialize the rest of the unit cell, guessing
        # compositions using a cosine function
        x = x_coord - (T/2)
        rhoB = 0.4 * math.cos(x * n_interfaces * math.pi / (L-T)) + 0.5
        rhoA = 1-rhoB
        lines.append(line.format(rhoA*rho,rhoB*rho))
    linenum += 1

# Write data to file
rgrid = open('c.rf','w')
rgrid.writelines(lines)
rgrid.close()
