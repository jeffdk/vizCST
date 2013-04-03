import numpy
from fortranFile import FortranFile

filename = "model/hotModelReal0.9_data.dump"
chunksize = 64

f = FortranFile(filename)

# Line 1:
# NSI, NUI, NSURF
print f.readArray(numpy.int32)
# Line 2:
# rpoe, RADEQUAT, COMEGA
print f.readArray(numpy.float64)
# Line 3:
# Etot, BaryM, AngMom, snp, ROTFUNC
print f.readArray(numpy.float64)

# Data
#  rho,gamma,alpha,omega,ed,angv

print f.readArray(numpy.float64)

