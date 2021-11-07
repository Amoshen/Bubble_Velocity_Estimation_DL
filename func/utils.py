# -*- coding: utf-8 -*-
"""

"""
import scipy.io


def turn(GT):
    dim = GT.shape
    for j in range(0,dim[1]):
        for i in range(0,dim[0]//2):
            temp    = GT[i,j]
            GT[i,j] = GT[dim[0]-1-i,j]
            GT[dim[0]-1-i,j] = temp
    return GT 

def SaveTrainResults(loss,loss_test,SavePath):
    data = {}
    data['loss_train'] = loss
    data['loss_test'] = loss_test
    scipy.io.savemat(SavePath+'TrainLoss',data)

def SaveTestResults(Prediction,GT,SavePath):
    data = {}
    data['GT']      = GT
    data['Prediction'] = Prediction
    scipy.io.savemat(SavePath+'TestResults',data) 
    
    

