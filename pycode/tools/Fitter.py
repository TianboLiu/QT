#!/usr/bin/env python
import os, sys
from scipy import sin, cos, tan, exp, log
from scipy.optimize import minimize
import numpy as np
from Lsidis import *

def LOADDATA(filename, oldset = np.array([])):
    dataset = np.loadtxt(filename, skiprows = 20)
    newset = np.vstack((oldset, dataset))
    return newset


















if __name__ == "__main__":
    sys.exit()
