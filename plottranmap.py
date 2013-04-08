#!/usr/bin/env python
import numpy as np
import matplotlib 
from matplotlib import pyplot as plt
from matplotlib import cm
def main():
    infile = 'temp'
    data = np.loadtxt(infile)
    mdat = np.ma.masked_array(data,data==0)
    #mm = np.mean(mdat)
    #cmp = data/data.max()
    #ax = plt.figure()
    plt.imshow(mdat)
    plt.xlabel("Nbin")
    plt.ylabel("Ntran")
    plt.title("Transit and NonTTV-model")
    plt.colorbar()
    plt.show()
    return
if __name__=='__main__':
    main()
