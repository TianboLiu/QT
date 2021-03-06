#!/usr/bin/env/ python
import numpy as np
import os, sys
from scipy import pi as PI
from scipy import sin, cos, exp, log

import lhapdf

pdfs = lhapdf.mkPDF("CJ15lo", 0)
ffpion = lhapdf.mkPDF("DSSFFlo", 211)


def xPDF(x, Q2, target = "proton"):
    xq = [pdfs.xfxQ2(2, x, Q2), pdfs.xfxQ2(1, x, Q2), pdfs.xfxQ2(3, x, Q2), pdfs.xfxQ2(4, x, Q2), pdfs.xfxQ2(5, x, Q2)]
    xqb = [pdfs.xfxQ2(-2, x, Q2), pdfs.xfxQ2(-1, x, Q2), pdfs.xfxQ2(-3, x, Q2), pdfs.xfxQ2(-4, x, Q2), pdfs.xfxQ2(-5, x, Q2)]
    xg = [pdfs.xfxQ2(21, x, Q2),]
    if target == "proton":
        return {"u" : xq[0], "d" : xq[1], "s" : xq[2], "c" : xq[3], "b" : xq[4],\
                "ub" : xqb[0], "db" : xqb[1], "sb" : xqb[2], "cb" : xqb[3], "bb" : xqb[4],\
                "g" : xg[0]}
    elif target == "neutron":
        return {"u" : xq[1], "d" : xq[0], "s" : xq[2], "c" : xq[3], "b" : xq[4],\
                "ub" : xqb[1], "db" : xqb[0], "sb" : xqb[2], "cb" : xqb[3], "bb" : xqb[4],\
                "g" : xg[0]}
    elif target == "deuteron":
        return {"u" : xq[0] + xq[1], "d" : xq[1] + xq[0], "s" : xq[2] * 2, "c" : xq[3] * 2, "b" : xq[4] * 2,\
                "ub" : xqb[0] + xqb[1], "db" : xqb[1] + xqb[0], "sb" : xqb[2] * 2, "cb" : xqb[3] * 2, "bb" : xqb[4] * 2,\
                "g" : xg[0] * 2}
    elif target == "helium-3":
        return {"u" : xq[0] * 2 + xq[1], "d" : xq[1] * 2 + xq[0], "s" : xq[2] * 3, "c" : xq[3] * 3, "b" : xq[4] * 3,\
                "ub" : xqb[0] * 2 + xqb[1], "db" : xqb[1] * 2 + xqb[0], "sb" : xqb[2] * 3, "cb" : xqb[3] * 3, "bb" : xqb[4] * 3,\
                "g" : xg[0] * 3}
    else:
        print "xPDF : target not exists!"
        sys.exit()
        return

def zFF(z, Q2, hadron = "pi+"):
    zq = [ffpion.xfxQ2(2, z, Q2), ffpion.xfxQ2(1, z, Q2), ffpion.xfxQ2(3, z, Q2), ffpion.xfxQ2(4, z, Q2), ffpion.xfxQ2(5, z, Q2)]
    zqb = [ffpion.xfxQ2(-2, z, Q2), ffpion.xfxQ2(-1, z, Q2), ffpion.xfxQ2(-3, z, Q2), ffpion.xfxQ2(-4, z, Q2), ffpion.xfxQ2(-5, z, Q2)]
    zg = [ffpion.xfxQ2(21, z, Q2),]
    if hadron == "pi+":
        return {"u" : zq[0], "d" : zq[1], "s" : zq[2], "c" : zq[3], "b" : zq[4],\
                "ub" : zqb[0], "db" : zqb[1], "sb" : zqb[2], "cb" : zqb[3], "bb" : zqb[4],\
                "g" : zg[0]}
    elif hadron == "pi-":
        return {"u" : zqb[0], "d" : zqb[1], "s" : zqb[2], "c" : zqb[3], "b" : zqb[4],\
                "ub" : zq[0], "db" : zq[1], "sb" : zq[2], "cb" : zq[3], "bb" : zq[4],\
                "g" : zg[0]}
    elif hadron == "pi0":
        return {"u" : 0.5 * (zq[0] + zqb[0]), "d" : 0.5 * (zq[1] + zqb[1]), "s" : 0.5 * (zq[2] + zqb[2]), "c" : zq[3] + zqb[3], "b" : zq[4] + zqb[4],\
                "ub" : 0.5 * (zq[0] + zqb[0]), "db" : 0.5 * (zq[1] + zqb[1]), "sb" : 0.5 * (zq[2] + zqb[2]), "cb" : zq[3] + zqb[3], "bb" : zq[4] + zqb[4],\
                "g" : zg[0]}
    else:
        print "zFF : hadron not exists!"
        sys.exit()
        return

def FT(x, Q2, target = "proton"):
    xf = xPDF(x, Q2, target)
    return (2.0 / 3.0)**2 * (xf["u"] + xf["ub"] + xf["c"] + xf["cb"]) + (1.0 / 3.0)**2 * (xf["d"] + xf["db"] + xf["s"] + xf["sb"] + xf["b"] + xf["bb"])
    
def model_0(x, Q2, z, Pt, par, target = "proton", hadron = "pi+"):#gaussian
    xf = xPDF(x, Q2, target)
    zD = zFF(z, Q2, hadron)
    W = z**2 * par[0] + par[1]
    factor = exp(-Pt**2 / W) / (PI * W * z)
    FUUT = factor * ( (2.0 / 3.0)**2 * (xf["u"] + xf["ub"] + xf["c"] + xf["cb"]) + (1.0 / 3.0)**2 * (xf["d"] + xf["db"] + xf["s"] + xf["sb"] + xf["b"] + xf["bb"]) )
    multi = 2.0 * PI * Pt * FUUT / FT(x, Q2, target)#GeV^-1
    return multi
    
def Test():
    return
                
    
def F2(x, Q2, target = "proton", NC = "gamma"):
    xf = xPDF(x, Q2, target)
    if NC == "gamma":
        return (2.0 / 3.0)**2 * (xf["u"] + xf["ub"] + xf["c"] + xf["cb"]) + (1.0 / 3.0)**2 * (xf["d"] + xf["db"] + xf["s"] + xf["sb"] + xf["b"] + xf["bb"])
    elif NC == "gammaZ":
        return 2.0 * (2.0 / 3.0) * (0.5 - 2.0 * 2.0 / 3.0 * 0.23) * (xf["u"] + xf["ub"] + xf["c"] + xf["cb"]) + 2.0 * (-1.0 / 3.0) * (-0.5 + 2.0 * 1.0 / 3.0 * 0.23) * (xf["d"] + xf["db"] + xf["s"] + xf["sb"] + xf["b"] + xf["bb"])
    elif NC == "Z":
        return ((0.5 - 2.0 * 2.0 / 3.0 * 0.23)**2 + 0.5**2) * (xf["u"] + xf["ub"] + xf["c"] + xf["cb"]) + ((-0.5 + 2.0 * 1.0 / 3.0 * 0.23)**2 + 0.5**2) * (xf["d"] + xf["db"] + xf["s"] + xf["sb"] + xf["b"] + xf["bb"])
    else:
        print "F2 : wrong NC option!"
        sys.exit()

def F1(x, Q2, target = "proton", NC = "gamma"):
    return F2(x, Q2, target, NC) / (2.0 * x)


if __name__ == "__main__":
    sys.exit()
