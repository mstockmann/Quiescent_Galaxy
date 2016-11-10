

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

sys.path.insert(0, '..')
import Stokky as st

####################################################


def read_in_1d_fits(path):
    """ Read in 1d fits file
        F: flux, E: error, W: wavelength (Angstrom), hd: header

        Returns wavelength, flux, error, header
    """
    data_arr = pf.open(path)
    hdf = data_arr[0].header
    hde = data_arr[0].header
    F = data_arr[0].data
    E = data_arr[1].data
    W = (hdf['CRVAL1'] + (hdf['CRPIX1'] - 1 + np.arange(hdf['NAXIS1']))*hdf['CDELT1'])*10 
    return W, F, E, hdf, hde

def load_spec_z(path,obj):
    data = np.genfromtxt(path,dtype=str)
    Name = data[:,0]
    z_spec = data[:,1]
    loc_name = np.where(obj == Name)[0][0]
    return np.float(z_spec[loc_name])



####################
#### 2D spectra ####
####################

Path_2d_spec_new = glob.glob('../../X-shooter/cQGsample/Objects/*/NIR_corr/*_V1_NIRcorr_wmrebin15_opt.fits')


Path_2d_spec_old = glob.glob('../../X-shooter/cQGsample/Objects/*/NIR_corr/*_V1_NIRcorr_wmrebin15_opt_old.fits')


for i in range(len(Path_2d_spec_new)):

    target = Path_2d_spec_new[i].split('/')[-1].split('_')[0]
    z_spec = load_spec_z('../../X-shooter/cQGsample/Collected_info/z_spec.txt',target)
    # z_spec = 2.69
    W_new, F_new, E_new, hdf_new, hde_new = read_in_1d_fits(Path_2d_spec_new[i])   
    W_old, F_old, E_old, hdf_old, hde_old = read_in_1d_fits(Path_2d_spec_old[i])       

    W_new_rf, F_new_rf, E_new_rf = st.convert_2_rest_frame(W_new, F_new, E_new, z_spec)
    W_old_rf, F_old_rf, E_old_rf = st.convert_2_rest_frame(W_old, F_old, E_old, z_spec)



    # Plot
    from scipy.ndimage.filters import gaussian_filter
    F_new_rf_conv = gaussian_filter(F_new_rf, sigma=2)
    F_new_conv = gaussian_filter(F_new, sigma=2)

    plt.figure(figsize=(14,10))
    plt.subplot(2,1,1)
    plt.plot(W_new, F_new,color='black')
    plt.plot(W_new, F_new_conv,color='lime')
    
    plt.plot(W_new, E_new, color='purple', alpha=0.2)
    plt.xlabel('$\lambda [\AA{}]$')
    plt.ylabel('Flux [erg/s/cm2/$\AA{}$]')
    plt.title('%s' % target)


    plt.subplot(2,1,2)
    plt.plot(W_new_rf, F_new_rf,color='black')
    plt.plot(W_new_rf, F_new_rf_conv,color='lime',linewidth=2)
    plt.plot(W_new_rf, E_new_rf, color='purple', alpha=0.2)

    st.plot_emission_lines()
    st.plot_absorption_lines()

    plt.xlabel('$\lambda [\AA{}]$')
    plt.ylabel('Flux [erg/s/cm2/$\AA{}$]')
    # plt.title('%s' % target)
    plt.axis([min(W_new_rf)-500,max(W_new_rf)+500,-0.4e-17, 1.8e-17])
    plt.show()
    






