import numpy as np
import matplotlib.pyplot as plt
#fam_1bc = np.loadtxt('fam_open_1bc_morepoints.csv', delimiter=',')
#fam_2bc = np.loadtxt('fam_open_2bc_morepoints.csv', delimiter=',')
fam_1bc = np.loadtxt('fam_Gd220_1bc.csv', delimiter=',')
fam_2bc = np.loadtxt('fam_Gd220_2bc.csv', delimiter=',')

col1 = fam_1bc[:,3]
col2 = fam_2bc[:,3]


#plt.plot(col1)
#plt.plot(col2)
plt.plot(fam_1bc[:,1], col1)
plt.show()

