#!/usr/bin/env python
import numpy as np
import scipy as sp
import phasefit as PF
import maskarray
import transitmask
import fitmap  as FM
from dataio import *
def testphase():
    #infile = 'data/HAT-365-0014102.tfalc' #pure outlier
    infile = 'data/HAT-364-0000605.tfalc' #outlier + flatpart of curve
    #infile = 'data/HAT-317-0001481.ltf' #transit
    time = []; readcolumn(time,1,infile);time=np.array(time)
    mag = []; readcolumn(mag,2,infile);mag=np.array(mag)
    #period = 16.3816800  
    period = 6.1594446
    #period = 7.2491767
    #epoch = 55513.5298468
    epoch = 55508.9172669
    #epoch = 55507.2281505
    #q = 0.0109
    q = 0.0127
    #q = 0.0202
    epoch += period*q/2.
    n0tran=round((time[0]-epoch-0.5*period)/period)
    if n0tran<0: 
        epoch+=n0tran*period
    #n0tran=round((time[0]-epoch-0.5*period)/period)
    #print n0tran,epoch
    #print time.shape
    phaselc = PF.PhaseLc(time,mag,period,epoch,10)
    #phaselc.Foldts()
    phaselc.TransitColor(0.04,400)
    return

def testfitsigle():
    data = np.loadtxt("data/color3.txt")
    #data = np.loadtxt("temp1")
    data = np.swapaxes(data,0,1).copy(order='C')
    fits = FM.Fitmap(data)
    #fits.FitSingleMask()
    #e1 = fits.Error()
    #print e1
    #fits.FitSingleFlat(10)
    #e2 = fits.Error()
    fits.FitTTV()
    #e3 = fits.Error()
    #fits.FitTTV(3)
    #e4 = fits.Error()
    #print e1,e2
    #print e3,e4
    #print e1,e2,e3
    fits.StdOutput()
    return

if __name__ == '__main__':
    #testphase()
    testfitsigle()
