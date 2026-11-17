#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

@author: kexinchen
"""

with open('out/c_rf') as f:
    content = f.readlines()

# The header ends with the line that follows the label ngrid (or mesh)
# and gives the mesh dimensions. Its length is not fixed: a field
# without a space group, such as the one written after COMPRESS, has
# no group_name entry.
i = 0
while content[i].split()[:1] not in (['ngrid'], ['mesh']):
    i += 1
mesh = [int(x) for x in content[i+1].split()]

#concentration matrix 
c = []
for line in content[i+2:]:
    if line.split():
        c.append([float(x) for x in line.split()])

# Check that a row was read for every grid point
nGrid = 1
for n in mesh:
    nGrid *= n
if len(c) != nGrid:
    raise Exception('Read ' + str(len(c)) + ' rows for a mesh of '
                    + str(nGrid) + ' grid points')

#Cound the summation of concentration at each grid differ from 1
diff = 0
for i in range(len(c)):
    sum_c = 0
    for j in range(len(c[i])):
        #sum the the concentation at grid i
        sum_c += c[i][j]
    diff += abs(sum_c -1)

print("The sum of inaccuracy of concentration on each grid is ")
print(diff)
        