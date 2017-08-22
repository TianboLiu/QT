#!/usr/bin/env/ python
import numpy as np
import os, sys
from scipy import pi as PI
from scipy import sin, cos, exp, log

import lhapdf

pdfs = lhapdf.mkPDF("CJ15lo", 0)
ffpion = lhapdf.mkPDF("DSSFFlo", 211)

def Test():
    print "xf(x): x = 0.1, Q = 3.0"
    print pdfs.xfxQ(2, 0.1, 3.0)
    print pdfs.xfxQ(1, 0.1, 3.0)
    print "zD(z): z = 0.2, Q = 3.0"
    print ffpion.xfxQ(2, 0.2, 3.0)
    print ffpion.xfxQ(1, 0.2, 3.0)
    return




if __name__ == "__main__":
    sys.exit()
