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

#------------------------------------------------#
#------------------------------------------------#
'''
path = glob.glob('data_cube_ex/nodding_macs2129_full_nir_9.6.main_box.fits')

f = fits.open(path[0])

print f.info()

F1 = f[0].data
header_cube = f[0].header
print header_cube

# NAXIS=3, NAXIS1=1, NAXIS2=1 and NAXIS3=n_wavelengths


# print f[0].header
# print f[1].header


# # print f[0].header
# raise


path_1d = glob.glob('../../X-shooter/cQGsample/Objects/105842/NIR_corr/105842_avwCombined_OBx4_sig5_V1_NIRcorr_wmrebin15_opt.fits')

f1 = fits.open(path_1d[0])
hdf = f1[0].header
F = f1[0].data
E = f1[1].data
M = f1[2].data
W = (hdf['CRVAL1'] + (hdf['CRPIX1'] - 1 + np.arange(hdf['NAXIS1']))*hdf['CDELT1'])*10 


F_3d = np.zeros((len(F), 1, 1),np.float32)
E_3d = np.zeros((len(E), 1, 1),np.float32)
M_3d = np.zeros((len(M), 1, 1),np.float32)

F_3d[:,0,0] = F
E_3d[:,0,0] = E
M_3d[:,0,0] = M

print F_3d.shape



header_cube['NAXIS3'] = F_3d.shape[0]
print header_cube

print hdf['CRPIX1']
print hdf['CRVAL1']
print hdf['CDELT1']
# print hdf

header_cube['CRPIX3'] = 1
header_cube['CRVAL3'] = hdf['CRVAL1']*10
header_cube['CDELT3'] = hdf['CDELT1']*10
print header_cube


path_out = 'data_cube_ex/test.fits'
fits.writeto(path_out, F_3d, header_cube)
'''


# # Read out error array
# pf.append(path, arr2, hd2)

# # Read out bad pixel map
# pf.append(path, arr3, hd3)


f_test = fits.open('data_cube_ex/test.fits')

print f_test.info()

















# # NIR; rebin, 9 A, optimal ext 
# Path_1d_spectra = glob.glob('../../X-shooter/cQGsample/Objects/*/NIR_corr/*_V1_NIRcorr_wmrebin15_opt.fits')


# mask_emission = read_in_emission_mask()

# for i in range(len(Path_1d_spectra)):
#     target = Path_1d_spectra[i].split('/')[-1].split('_')[0]
#     header += '%s ' % target
#     z_spec = st.load_spec_z('../../X-shooter/cQGsample/Collected_info/z_spec.txt',target)

#     # read in 1d spectrum
#     W, F, E, M, hdf, hde, hdm = st.read_in_1d_fits(Path_1d_spectra[i])

#     # Read out to data-cube














