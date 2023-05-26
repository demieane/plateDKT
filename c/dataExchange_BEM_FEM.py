# This program prints Hello, world!
import numpy as np
import sys
import scipy as sp



print('Hello, from Python! On the following days I will become more useful :) \n')

#Open a file for reading
#the access mode is r


fileID = open("OUTDATA_FEM_DEBUG.txt")
rowsA = int(fileID.readline())
colsA = int(fileID.readline())
AA = np.empty(shape=(rowsA,colsA))
print(rowsA,colsA) 
print(type(rowsA))

for i in range(0,rowsA):
    for j in range(0,colsA):
        AA[i][j]=float(fileID.readline())
#print(i)
#print(j)

rowsb = int(fileID.readline())
colsb = int(fileID.readline())
rhs = np.empty(shape=(rowsb))
print(rowsb,colsb) 
print(type(rowsb))

for i in range(0,rowsb):
    rhs[i]=float(fileID.readline())

#for i in range(0,10):
 #   var = rhs[i]/10**(-6)
print(AA[0][0])
print(rhs[0])

#The generic, symmetric, Hermitian and positive definite solutions 
# are obtained via calling ?GESV, ?SYSV, ?HESV, and ?POSV 
# routines of LAPACK respectively.

#The solutions are computed using LAPACK routine _gesv.
x = np.linalg.solve(AA, rhs)

print(x[0:10])
#sizeM = int(rowsA/2)
#print(sizeM)
#print(x[sizeM-4:sizeM+10])

np.linalg.inv(AA)

print(np.linalg.cond(AA))


