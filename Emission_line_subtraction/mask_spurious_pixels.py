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


# NIR; rebin, 9 A, optimal ext 
Path_1d_spectra = glob.glob('../../X-shooter/cQGsample/Objects/*/NIR_corr/*_V1_NIRcorr_wmrebin15_opt.fits')
# print len(Path_1d_spectra)

# mask_emission = [[193,194,195,196,197,305,306,307,568,569,570,571,572],[141,142,143,229,230,231,232,233,234,235,236,237,238,810,811,812],[],[],[],[],[],[],[],[],[],[],[],[]]

# print mask_emission

# # f = open('mask/emission_id.txt','r')

# # f.read()



mask_emission = read_in_emission_mask()

header = ''
W0, F0, E0, M0, hdf0, hde0, hdm0 = st.read_in_1d_fits(Path_1d_spectra[0])
mask_arr = np.zeros(shape=(len(F0),len(Path_1d_spectra)))+1

i_start = 0
for i in range(i_start,len(Path_1d_spectra),1):
    target = Path_1d_spectra[i].split('/')[-1].split('_')[0]
    header += '%s ' % target
    z_spec = st.load_spec_z('../../X-shooter/cQGsample/Collected_info/z_spec.txt',target)

    # read in 1d spectrum
    W, F, E, M, hdf, hde, hdm = st.read_in_1d_fits(Path_1d_spectra[i])

    if i == i_start:
        ''' cut the K-band (blocking) off for p86
        '''
        lim = len(F)


    mask_arr[:,i] = np.copy(M[:lim])
    mask_arr = masking_emission(W,mask_arr,mask_emission,i)


    # convert to restframe
    W_rf, F_rf, E_rf = st.convert_2_rest_frame(W[:lim], F[:lim], E[:lim], z_spec)

     
    # remove the masked regions
    F_rf_mask = F_rf[mask_arr[:,i] == 1]#*mask_arr[:,i]
    E_rf_mask = E_rf[mask_arr[:,i] == 1]#*mask_arr[:,i]
    W_rf_mask = W_rf[mask_arr[:,i] == 1]

    plt.figure(figsize=(14,10))
    plt.subplot(3,1,1)
    plt.title('%s' % target)
    # plt.plot(W, F, color='black')
    plt.plot(W_rf, F_rf, color='black')
    plt.plot(W_rf_mask, F_rf_mask, color='green')
    plt.plot(W_rf, E_rf, color='purple', alpha=0.2)

    st.plot_emission_lines()
    st.plot_absorption_lines()
    plt.xlabel('$\lambda [\AA{}]$')
    plt.ylabel('Flux [erg/s/cm2/$\AA{}$]')
    plt.axis([min(W_rf)-500,max(W_rf)+500,-0.4e-17, 1.8e-17])

    plt.subplot(3,1,2)
    plt.plot(np.arange(len(F_rf)), F_rf, color='black')
    plt.plot(np.arange(len(F_rf)), E_rf, color='purple', alpha=0.2)

    plt.ylim([-0.4e-17, 1.8e-17])

    plt.subplot(3,1,3)
    plt.plot(W, F, color='black')
    plt.plot(W, E, color='purple', alpha=0.2)


    plt.show()







# make empty array
np.savetxt('mask/Bad_pixel_mask.txt',mask_arr,header=header,fmt='%d',delimiter='\t')












