import re
import numpy

def getList(dict):
    list = []
    for key in dict.keys():
        list.append(key)
          
    return list

def GaussNoise(DataList, sigma):
    NewData=[]
    for i in DataList:
        NewData.append(numpy.random.normal(i, sigma, 1)[0])
    return NewData