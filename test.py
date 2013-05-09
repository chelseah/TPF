#!/usr/bin/env python
import numpy as np
import scipy as sp
import matplotlib 
from matplotlib import pyplot as plt
import blssimple as BS
import phasefit as PF
import maskarray
import transitmask
import fitmap  as FM
from dataio import *
def testphase():
    #infile = 'data/HAT-365-0014102.tfalc' #pure outlier
    #infile = 'data/HAT-364-0000605.tfalc' #outlier + flatpart of curve
    #infile = 'data/HAT-317-0001481.ltf' #transit
    #infile = 'data/kplr011401755.rmtran'
    #infile = 'data/CBPs/Kepler16/kplr012644769.ltf'
    infile = 'data/kplr009347899.ltf'
    time = []; readcolumn(time,1,infile);time=np.array(time)
    mag = []; readcolumn(mag,7,infile);mag=np.array(mag)
    #period = 16.3816800  
    #period = 6.1594446
    #period = 7.2491767
    #period = 13.8550595
    #period = 110.5253269*2
    period = 185.0
    #epoch = 55513.5298468
    #epoch = 55508.9172669
    #epoch = 55507.2281505
    #epoch = 2454974.3977524
    #epoch = 2454993.4068241
    epoch = 2454900.0
    #q = 0.0109
    #q = 0.0127
    #q = 0.0202
    #q = 0.03
    q = 0.0015
    epoch += period*q/2.
    n0tran=round((time[0]-epoch-0.5*period)/period)
    if n0tran<0: 
        epoch+=n0tran*period
    n0tran=round((time[0]-epoch-0.5*period)/period)
   # print n0tran,epoch
    #print time.shape
    phaselc = PF.PhaseLc(time,mag,period,epoch,200)
    #magbin = np.zeros(200)
    #phases = np.arange(200)*1.0/200.0
    #ntran = int((max(time)-epoch-0.5*period)/period)
    #nbin = int(2*q*400)
    #color = np.zeros([nbin,ntran])
    #phaselc.Foldts()
    #phaselc.OutputPhase(magbin)
    #print magbin
    #phaselc.StandOutputPhase()
    phaselc.TransitColor(0.5,700)
    #phaselc.OutputColor(color)
    #mdat = np.ma.masked_array(color,color==0)

    #plt.imshow(mdat)
    phaselc.StandOutputColor()
    #plt.plot(phases,magbin)
    #plt.show()
    return

def testfitsigle():
    data = np.loadtxt("data/Kepler36/Kepler36b")
    #data = np.loadtxt("temp1")
    data = np.swapaxes(data,0,1).copy(order='C')
    fits = FM.Fitmap(data)
    #fits.FitSingleMask()
    #e1 = fits.Error()
    #print e1
    #fits.FitSingleFlat(10)
    #e2 = fits.Error()
    #fits.FitTTV()
    #e3 = fits.Error()
    fits.FitTTV(5)
    #e4 = fits.Error()
    #print e3
    #print e1,e2
    #print e2,e3,e4
    #print e1,e2,e3
    #fits.StdOutput()    
    return

def testBLS():
    #infile = 'data/Kepler2.fklc'
    #infile = 'data/HAT-154-0001553.tfalc' #Kepler2
    infile = 'data/HAT-154-0028698.tfalc' #KOI531
    time = []; readcolumn(time,2,infile);time=np.array(time)
    mag = []; readcolumn(mag,18,infile);mag=np.array(mag)
    #Kepler2
    #P0 = 2.2047355
    #epoch = 54.35780+54900-P0*200
    #dip = (22.29*0.01)**2.
    #qvar = 0.0733
    #KOI531
    P0 = 3.6875
    epoch = 170.88+54833
    dip = 0.0027
    qvar = 0.012
    pmax = 10.0 
    pmin = 1.0
    fn = 200000 
    bls=BS.BlsSpec(time,mag,pmax,pmin,fn,epoch,qvar,dip)
    bls.GenSpec()
    bls.StandOutput()
    return 
if __name__ == '__main__':
    #testphase()
    #testfitsigle()
    testBLS()
