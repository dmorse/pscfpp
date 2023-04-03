#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

@author: kexinchen
"""

import numpy as np
import pandas as pd

def isfloat(num):
    try:
        float(num)
        return True
    except ValueError:
        return False


with open('in/w.rf') as f:
    content = f.readlines()

content = [x.strip() for x in content]

#w matrix
w = []
for line in content:
    mat =[]
    s = line.split(' ')
    for j in range(len(s)):
        if(isfloat(s[j])):
            mat.append(float(s[j]))
    w.append(mat)

#Remove nx and nm value
w = w[15:]

#Cound the summation of concentration at each grid differ from 1
w1 = []
w2 = []
for i in range(10):
    w1.append((w[i][0]-w[i][1])/2)
    w2.append((w[i][0]+w[i][1])/2)

print("w1")
print(w1)
        
print("w2")
print(w2)
