from cstModel import cstModel

modelName = 'c30p10'

filename = "model/hotModelReal0.9_data.dump"
filename = "model/tov_" + modelName + "/outdata0.3000"


model = cstModel(filename)


model.rawPlot('ed', modelName)