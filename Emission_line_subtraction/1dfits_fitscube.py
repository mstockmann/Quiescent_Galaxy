#!/usr/local/bin/python

#-*- coding: utf-8 -*-


from __future__ import division
import glob
# import pyfits as pf
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from scipy.interpolate import splrep,splev
from scipy import ndimage
import mpl_toolkits.axisartist as AA
import sys
from scipy import interpolate

sys.path.insert(0, '..')
import Stokky as st

#------------------------------------------------#
#------------------ Functions -------------------#
#------------------------------------------------#

def masking_emission(W,Arr,mask_em,spec_no):
    ID_1 = np.array([0,np.searchsorted(W,[10400])[0]])
    ID_2 = np.searchsorted(W,[13500,14240])
    ID_3 = np.searchsorted(W,[17919,19641])
    ID_4 = np.array([np.searchsorted(W,[20620])[0],-1])

    Arr[ID_1[0]:ID_1[1],spec_no] = 0
    Arr[ID_2[0]:ID_2[1],spec_no] = 0
    Arr[ID_3[0]:ID_3[1],spec_no] = 0
    Arr[ID_4[0]:ID_4[1],spec_no] = 0

    ID_em = np.searchsorted(W,mask_em[spec_no])
    for jj in range(len(mask_em[spec_no])):
        Arr[mask_em[spec_no][jj],spec_no] = 0
    return Arr

def read_in_emission_mask():
    mask_em = []
    lines = [line.rstrip('\n') for line in open('mask/Spurious_pixel_id.txt')]
    lines = [inside.split(',') for inside in lines]
    for ll in lines:
        tmp = [int(k) for k in ll]
        mask_em.append(tmp)
    return mask_em

def update_header(Arr,hdr,hdr_1d,hdr_no):
    if hdr_no == 0:
        hdr['NAXIS3']   = Arr.shape[0]
        hdr['EXTNAME']  = 'FLUX'
        hdr['DATE']     = '2019-11-11'
        hdr['CONTENT']  = 'Flux [erg/cm2/s/AA]'
        hdr['CRPIX3']   = 1
        hdr['CRVAL3']   = hdr_1d['CRVAL1']*10
        hdr['CDELT3']   = hdr_1d['CDELT1']*10

    if hdr_no == 1:
        hdr['NAXIS3']   = Arr.shape[0]
        hdr['EXTNAME']  = 'ERROR'
        hdr['DATE']     = '2019-11-11'
        hdr['CONTENT']  = 'Error [erg/cm2/s/AA]'
        hdr['CRPIX3']   = 1
        hdr['CRVAL3']   = hdr_1d['CRVAL1']*10
        hdr['CDELT3']   = hdr_1d['CDELT1']*10
    
    if hdr_no == 2:
        hdr['NAXIS3']   = Arr.shape[0]
        hdr['EXTNAME']  = 'BADPIX'
        hdr['DATE']     = '2019-11-11'
        hdr['CONTENT']  = 'bpm, =1 bad'
        hdr['CRPIX3']   = 1
        hdr['CRVAL3']   = hdr_1d['CRVAL1']*10
        hdr['CDELT3']   = hdr_1d['CDELT1']*10
    return hdr

def Swap_0_and_1(Map):
    """ In the Map: Good=1, Bad=0
        Here we swap them:
        New BPM: Good=0, Bad=1 
    """
    Map_swap = np.zeros(len(Map))
    for j in range(len(Map_swap)):
        if Map[j] == 0:
            Map_swap[j] = 1
        else:
            Map_swap[j] = 0
    return Map_swap


def save_3d_cube_fits(path,F,E,M,hdf,hde,hdm):
    ''' 1) Load example fits cube

    '''
    # read in example
    f_ex = fits.open('data_cube_ex/nodding_macs2129_full_nir_9.6.main_box.fits')
    hd0 = f_ex[0].header
    hd1 = f_ex[1].header
    hd2 = f_ex[2].header

    # make new 3d cube
    F_3d = np.zeros((len(F), 1, 1))
    E_3d = np.zeros((len(E), 1, 1))
    M_3d = np.zeros((len(M), 1, 1))

    F_3d[:,0,0] = F
    E_3d[:,0,0] = E
    M_g0_b1 = Swap_0_and_1(M)
    M_3d[:,0,0] = M_g0_b1

    hd0 = update_header(F_3d,hd0,hdf,0)
    hd1 = update_header(F_3d,hd1,hde,1)
    hd2 = update_header(F_3d,hd2,hdm,2)

    fits.writeto(path, F_3d, hd0)
    fits.append(path, E_3d, hd1)
    fits.append(path, M_3d, hd2)

#------------------------------------------------#
#------------------------------------------------#









path_1d = glob.glob('../../X-shooter/cQGsample/Objects/105842/NIR_corr/105842_avwCombined_OBx4_sig5_V1_NIRcorr_wmrebin15_opt.fits')

f1 = fits.open(path_1d[0])
hdf1 = f1[0].header
hde1 = f1[0].header
hdm1 = f1[0].header

F1 = f1[0].data
E1 = f1[1].data
M1 = f1[2].data
W1 = (hdf1['CRVAL1'] + (hdf1['CRPIX1'] - 1 + np.arange(hdf1['NAXIS1']))*hdf1['CDELT1'])*10 


# path_out = 'data_cube_ex/f_test.fits'
# save_3d_cube_fits(path_out,F1,E1,M1,hdf1,hde1,hdm1)




















