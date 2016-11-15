#!/usr/local/bin/python

#-*- coding: utf-8 -*-


from __future__ import division
import glob
import pyfits as pf
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


def masking_no_info_area(W,Arr,mask_em,spec_no):
    ID_1 = np.array([0,np.searchsorted(W,[10400])[0]])
    ID_2 = np.searchsorted(W,[13500,14240])
    ID_3 = np.searchsorted(W,[17919,19641])
    ID_4 = np.array([np.searchsorted(W,[20620])[0],-1])

    Arr[ID_1[0]:ID_1[1]] = 0
    Arr[ID_2[0]:ID_2[1]] = 0
    Arr[ID_3[0]:ID_3[1]] = 0
    Arr[ID_4[0]:ID_4[1]] = 0

    ID_em = np.searchsorted(W,mask_em[spec_no])
    for jj in range(len(mask_em[spec_no])):
        Arr[mask_em[spec_no][jj]] = 0
    return Arr    

def read_in_spurious_mask():
    mask_em = []
    lines = [line.rstrip('\n') for line in open('mask/Spurious_pixel_id_wmrebin15_NIR.txt')]
    lines = [inside.split(',') for inside in lines]
    for ll in lines:
        tmp = [int(k) for k in ll]
        mask_em.append(tmp)
    return mask_em

#------------------------------------------------#
#------------------------------------------------#



##### OBS #####
# Notes 171060 remove spurios pixel at 15425





# # NIR; rebin, 9 A, optimal ext 
# Path_1d_spectra = glob.glob('../../X-shooter/cQGsample/Objects/*/NIR_corr/*_V1_NIRcorr_wmrebin15_opt.fits')


# NIR; rebin, 9 A, optimal ext & 2d emission removal + stdcorr
Path_1d_spectra = glob.glob('../../X-shooter/cQGsample/Spectra_analysis/3-Stdcorr_slit_aper/NIR/*_V1_NIRcorr_wmrebin15_er2d_opt_stdcorr.fits')
# print len(Path_1d_spectra)
# print Path_1d_spectra
# raise

plot_yes = 1

mask_spurious = read_in_spurious_mask()
W0, F0, E0, M0, hdf0, hde0, hdm0 = st.read_in_1d_fits(Path_1d_spectra[0])
mask_arr = np.zeros(shape=(len(F0),len(Path_1d_spectra)))+1
lim = len(F0)
print lim

for i in range(3,len(Path_1d_spectra),1):
    target = Path_1d_spectra[i].split('/')[-1].split('_')[0]

    z_spec = st.load_spec_z('../../X-shooter/cQGsample/Collected_info/z_spec_notcalib.txt',target)
    # z_spec = 2.02#1.712

    # read in 1d spectrum
    W, F, E, M, hdf, hde, hdm = st.read_in_1d_fits(Path_1d_spectra[i])
    hdf['REDSHIFT'] = (z_spec, 'Spectroscopic redshift')

    # print hdf['CRVAl1']
    # raise
    # M good = 1, bad =0
    # plt.plot(W,M)
    # plt.show()

    # # Do this in separate code
    # if i == 0:
    #     """ cut the K-band (blocking) off for p86
    #     """
        
    M = masking_no_info_area(W,M[:lim],mask_spurious,i)

    if plot_yes == 1:
        # convert to restframe
        W_rf, F_rf, E_rf = st.convert_2_rest_frame(W[:lim], F[:lim], E[:lim], z_spec)

        
        # remove the masked regions
        F_rf_mask = F_rf[M == 1]
        E_rf_mask = E_rf[M == 1]
        W_rf_mask = W_rf[M == 1]

        plt.figure(figsize=(14,10))
        plt.subplot(3,1,1)
        plt.title('%s' % target)
        plt.plot(W_rf, F_rf, color='black')
        plt.plot(W_rf_mask, F_rf_mask, color='green')

        from scipy.ndimage.filters import gaussian_filter
        F_new_rf_conv = gaussian_filter(F_rf_mask, sigma=2)
        plt.plot(W_rf_mask, F_new_rf_conv, color='purple')

        # plt.plot(W_rf, E_rf, color='purple', alpha=0.2)
        plt.fill_between(W_rf, np.zeros(len(E_rf)), E_rf, color='purple', alpha=0.2,edgecolor='purple')

        st.plot_emission_lines()
        st.plot_absorption_lines()
        plt.xlabel('$\lambda [\AA{}]$')
        plt.ylabel('Flux [erg/s/cm2/$\AA{}$]')
        plt.axis([min(W_rf)-500,max(W_rf)+500,-0.4e-17, 1.3e-17])

        plt.subplot(3,1,2)
        plt.plot(np.arange(len(F_rf)), F_rf, color='black')
        plt.plot(np.arange(len(F_rf)), E_rf, color='purple', alpha=0.2)

        plt.ylim([-0.4e-17, 1.7e-17])

        plt.subplot(3,1,3)
        plt.plot(W, F, color='black')
        plt.plot(W, E, color='purple', alpha=0.2)

        plt.show()

    filename = Path_1d_spectra[i].split('/')[-1].replace('.fits','_bpm.fits') 
    path_out = '../../X-shooter/cQGsample/Spectra_analysis/4-Stefano_BPM_update/%s' % filename
    # st.save_1d_cube_fits(path_out,F,E,M,hdf,hde,hdm)
    # raise





# path = '/Users/mstockmann/PhD_DARK/Data/X-shooter/cQGsample/Spectra_analysis/4-Stefano_BPM_update/105842_V1_NIRcorr_wmrebin15_er2d_opt_bpm.fits'


# f = fits.open(path)
# print f.info()

# W, F, E, M, hdf, hde, hdm = st.read_in_1d_fits(path)

# plt.plot(W/10,F)
# plt.show()








