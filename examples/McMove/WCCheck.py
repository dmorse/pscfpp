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


with open('in/w_rf') as f:
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
for i in range(len(w)):
    w1.append((w[i][0]-w[i][1])/2)
    w2.append((w[i][0]+w[i][1])/2)

#Compressor_out
with open('out/w_out') as f:
    contento = f.readlines()

contento = [x.strip() for x in contento]

#w matrix
wo = []
for line in contento:
    mat =[]
    s = line.split(' ')
    for j in range(len(s)):
        if(isfloat(s[j])):
            mat.append(float(s[j]))
    wo.append(mat)

#Remove nx and nm value
wo = wo[15:]

#Cound the summation of concentration at each grid differ from 1
w1O = []
w2O = []
for i in range(len(w)):
    w1O.append((wo[i][0]-wo[i][1])/2)
    w2O.append((wo[i][0]+wo[i][1])/2)

#Compare W_
W_m = np.sum([x - y for x, y in zip(w1, w1O)])
W_p = np.sum([x - y for x, y in zip(w2, w2O)])
print("Difference of W- between input and compressor output ")
print(float(W_m))
print(float(W_p))

        
