#!/usr/bin/env python
import numpy as np
import matplotlib 
from matplotlib import pyplot as plt
from matplotlib import cm
def main():
    #infile = 'Kepler36b'
    infile = 'temp'
    #infile = 'data/Kepler36/Kepler36b'
    #infile = 'data/Kepler36/Kepler36b'
    data = np.loadtxt(infile)
    mdat = np.ma.masked_array(data,data==0)
    #mm = np.mean(mdat)
    #cmp = data/data.max()
    #ax = plt.figure()
    print data.shape
    asp = data.shape[1]/data.shape[0]
    print asp
    plt.imshow(mdat,aspect=asp,cmap='bone_r')
    plt.xlabel("Nbin")
    plt.ylabel("Ntran")
    plt.title("Kepler36b-TTV")
    plt.colorbar()
    plt.show()
    return
if __name__=='__main__':
    main()
