#!/usr/bin/env python
import os, sys
from scipy import sin, cos, tan, exp, log
from scipy.optimize import minimize
import numpy as np
from Lsidis import *

datalist = np.array([]).reshape(0,8)

def LOADDATA(filename, oldset = np.array([]).reshape(0,8)):
    dataset = np.loadtxt(filename, skiprows = 19)
    newset = np.vstack((oldset, dataset))
    return newset

def chi2(para):
    global datalist
    sum = 0.0
    for item in datalist:
        if item[4] < 1.5 and item[5] < 0.038 and item[6] < 0.2:#Q2
            print item
            sum += (model_0(x = item[5], Q2 = item[4], z = item[6], Pt = item[7], target = "proton", par = para) - item[1])**2 / (item[2]**2 + item[3]**2)
    return sum

def main():
    global datalist
    datalist = LOADDATA("Data/SIDIS/hermes.proton.zxpt-3D.vmsub.mults_piplus.list")
    res = minimize(chi2, [0.5, 0.5], method='Nelder-Mead', tol=1e-6)
    return res



if __name__ == "__main__":
    sys.exit()
