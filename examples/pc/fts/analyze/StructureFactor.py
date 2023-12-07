#Print structure vs kR0
#chen7545@umn.edu 

import numpy as np
import matplotlib.pyplot as plt


file_path = 'out/structureFactor'

data_S = [[] for _ in range(2)]

with open(file_path, 'r') as file:
    next(file)
    for line in file:
        data = list(map(float, line.split()))
        for i in range(2):
            data_S[i].append(data[i])
plt.figure()            
plt.plot(data_S[0], data_S[1],'-o')
plt.xlim([0,12])
plt.xlabel(r'$k R_0$')
plt.ylabel(r'$\frac{S(k)}{\rho_0 N}$')
plt.show()