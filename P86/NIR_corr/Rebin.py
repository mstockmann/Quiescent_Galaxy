

#!/usr/local/bin/python

#-*- coding: utf-8 -*-

from __future__ import division
import glob
import pyfits as pf
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import os
from scipy.optimize import curve_fit
import sys

sys.path.insert(0, '/Volumes/DataDrive/Quiescent_Galaxy')
import Stokky as st

####################################################


def weighted_mean(F,E,M,Np):
    N_zeros = abs(len(M) - np.count_nonzero(M))
    if N_zeros > Np/3: # if 68 percent is zeros
        M_wm = 1
    else:
        M_wm = 0

    w = 1/(E**2)
    F_wm = np.sum(w*F*M)/np.sum(w)
    E_wm = 1/np.sqrt(np.sum(w))
    return F_wm,E_wm, M_wm

def SwapZerosAndOnes(Map):
    Map_swap = np.zeros(shape=(len(Map[:,0]),len(Map[0,:])))
    for j in range(len(Map[:,0])):
        for i in range(len(Map[0,:])):
            if Map[j,i] == 0:
                Map_swap[j,i] = 1
            else:
                Map_swap[j,i] = 0
    return Map_swap








#######################################
##### Rebin a single 2D spectrum  #####
#######################################


# path = glob.glob('../../../X-shooter/P93/Data/Reduction/*/Combined_OBs/NIR_corr/*_avwCombined_OBx*_sig5_V1_NIR_corr.fits')

path = glob.glob('../../../X-shooter/P86_COSMOS_1/Data/*/Combined_OBs/NIR_corr/*_V1_NIR_corr.fits')


for k in range(len(path)):
    print 'Object: %s' % path[k].split('/')[6]

    Flux, Err, BPM, hd_0, hd_1, hd_2 = st.read_2d_fits(path[k])
    Wave = hd_0['CRVAL1'] + (hd_0['CRPIX1'] - 1 + np.arange(hd_0['NAXIS1']))*hd_0['CDELT1'] # Aangstroem

    # Swap 0 <-> 1 to use 
    BPM_alt = SwapZerosAndOnes(BPM)

    # Rebin 2D spectrum
    F = Flux
    E = Err
    M = BPM_alt

    nm2pix = hd_0['CDELT1']
    New_binsize_Ang = 9
    N_pix = New_binsize_Ang/(nm2pix*10)

    #N_pix = 60.0;# print 'Binning = %s pix (%s nm)' % (N_pix,N_pix*nm2pix)
    N_chunk = int(np.round(len(F[0,:])/N_pix,0)); print 'N_chunk = %s' % N_chunk


    Flux_bin = np.zeros(shape=(len(F[:,0]),N_chunk))
    Errs_bin = np.zeros(shape=(len(F[:,0]),N_chunk))
    Qmap_bin = np.zeros(shape=(len(F[:,0]),N_chunk))

    #eps = 1e-40 # Softening in inverse error w_err
    for i in range(int(N_chunk)):
        # The chunks of Flux and Err, multiplied by the M (bad=0, good=1)
        F_chunk = (F[:,i*N_pix:(i+1)*N_pix])
        E_chunk = (E[:,i*N_pix:(i+1)*N_pix])
        M_chunk = (M[:,i*N_pix:(i+1)*N_pix])

        for j in range(len(F[:,0])):
            # Weighted mean
            Flux_bin[j,i], Errs_bin[j,i], Qmap_bin[j,i] = weighted_mean(F_chunk[j,:],E_chunk[j,:],M_chunk[j,:],N_pix)


    ### Read out the fucking binned spectrum ###
    path_out = path[k].replace('_NIR_corr','_NIRcorr_wmrebin%s' % int(N_pix))
    # path_out = path[k].replace('V1','V1_wmrebin%s' % int(N_pix))

    # Changing the wavelength solution to a rebinnined version
    hd_0['CRVAL1'] = hd_0['CRVAL1']+(N_pix-1)*hd_0['CDELT1']/2
    hd_0['CDELT1'] = hd_0['CDELT1']*N_pix # CDELT1 and CD1_1 are the same
    hd_0['CD1_1'] = hd_0['CD1_1']*N_pix
    hd_0['NAXIS1'] = N_chunk

    st.read_out_3arr_2dfits(path_out,Flux_bin,Errs_bin,Qmap_bin,hd_0,hd_1,hd_2)



###############################
###############################















