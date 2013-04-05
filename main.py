import numpy
from cstModel import cstModel

modelName = 'c0p0'

filename = "model/hotModelReal0.9_data.dump"
filename = "model/a0.8_" + modelName + "/outdata0.3000"


model = cstModel(filename)

mass = model.integralOverStar(["ed"], lambda x: x, 10e10)
volume = model.integralOverStar([], lambda: 1, 10e10)

print "Average density: ", mass / volume

#model.rawPlot('ed', modelName + " a0.8", 30.0, numpy.log10)