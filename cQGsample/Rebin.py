

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

sys.path.insert(0, '../')
import Stokky as st

####################################################


def weighted_mean(F,E,M,Number_of_pixels):
    Number_of_bad_pixels = abs(len(M) - np.count_nonzero(M))
    if Number_of_bad_pixels > Number_of_pixels/3: # If more than 1/3 is bad masked (M=0)
        M_wm = 0 # Bad
    else:
        M_wm = 1 # good

    w = 1/(E**2)
    F_wm = np.sum(w*F*M)/np.sum(w)
    E_wm = 1/np.sqrt(np.sum(w))
    return F_wm, E_wm, M_wm

def SwapZerosAndOnes(Map):
    """ In the BPM: Good=0, Bad => 1
        Here we swap them:
        New BPM: Good = 1, Bad = 0 
    """
    Map_swap = np.zeros(shape=(len(Map[:,0]),len(Map[0,:])))
    for j in range(len(Map[:,0])):
        for i in range(len(Map[0,:])):
            if Map[j,i] == 0:
                Map_swap[j,i] = 1
            else:
                Map_swap[j,i] = 0
    return Map_swap




def make_header_rebin(hd,binsize,Nbins):
    hd['CRVAL1'] = hd['CRVAL1']+(binsize-1)*hd['CDELT1']/2
    hd['CDELT1'] = hd['CDELT1']*binsize # CDELT1 and CD1_1 are the same
    hd['CD1_1']  = hd['CD1_1']*binsize
    hd['NAXIS1'] = Nbins
    return hd



#######################################
##### Rebin a single 2D spectrum  #####
#######################################

path = glob.glob('../../X-shooter/cQGsample/Objects/*/NIR_corr/*_V1_NIR_corr.fits')

for k in range(len(path)):
    print 'Object: %s' % path[k].split('/')[-1].split('_')[0]
    
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


    # ### Read out the fucking binned spectrum ###
    path_out = path[k].replace('_NIR_corr','_NIRcorr_wmrebin%s_test' % int(N_pix))


    
    hd_0 = make_header_rebin(hd_0,N_pix,N_chunk)
    hd_1 = make_header_rebin(hd_1,N_pix,N_chunk)
    hd_2 = make_header_rebin(hd_2,N_pix,N_chunk)

    st.read_out_3arr_2dfits(path_out,Flux_bin,Errs_bin,Qmap_bin,hd_0,hd_1,hd_2)
    raise


###############################
###############################















