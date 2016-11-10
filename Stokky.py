
#!/usr/local/bin/python

#-*- coding: utf-8 -*-

from __future__ import division
import glob
import pyfits as pf
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import os
import string


#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@





def read_in_1d_fits(path):
    """ Read in 1d fits file
    	F: flux, E: error, W: wavelength (Angstrom), hd: header

    	Returns wavelength, flux, error, header
    """
    data_arr = pf.open(path)
    hdf = data_arr[0].header
    hde = data_arr[1].header
    hdm = data_arr[2].header
    F = data_arr[0].data
    E = data_arr[1].data
    M = data_arr[2].data
    W = (hdf['CRVAL1'] + (hdf['CRPIX1'] - 1 + np.arange(hdf['NAXIS1']))*hdf['CDELT1'])*10 
    return W, F, E, M, hdf, hde, hdm




def read_2d_fits(path):
    f = pf.open(path)

    hd1 = f[0].header
    hd2 = f[1].header
    hd3 = f[2].header
    Arr1 = f[0].data
    Arr2 = f[1].data
    Arr3 = f[2].data
    return Arr1, Arr2, Arr3, hd1, hd2, hd3


def read_out_3arr_2dfits(path,arr1,arr2,arr3,hd1,hd2,hd3):
    if not os.path.exists(path):
        # Read out flux array
        pf.writeto(path, arr1, hd1)
        
        # Read out error array
        pf.append(path, arr2, hd2)

        # Read out bad pixel map
        pf.append(path, arr3, hd3)

    else:
        os.system('rm %s' % path)
        # Read out flux array
        pf.writeto(path, arr1, hd1)
        
        # Read out error array
        pf.append(path, arr2, hd2)

        # Read out bad pixel map
        pf.append(path, arr3, hd3)
        print 'file already exists'


def convert_2_rest_frame(W,F,E,z):
    W_rest_frame = W / (1+z)
    F_rest_frame = F * (1+z)
    E_rest_frame = E * (1+z)
    return W_rest_frame, F_rest_frame, E_rest_frame


def load_spec_z(path,obj):
    data = np.genfromtxt(path,dtype=str)
    Name = data[:,0]
    z_spec = data[:,1]
    loc_name = np.where(obj == Name)[0][0]
    return np.float(z_spec[loc_name])


def plot_absorption_lines():
    # Absorption lines
    Mg_5169    = 5168.6
    Mg_5175    = 5175.4
    Mg_5183    = 5183.3
    Fe_5332    = 5331.5
    Fe_5401    = 5401.4

    CaI_4227   = 4226.7
    H_beta  = 4861.3
    H_gamma = 4340.5
    H_delta = 4101.7
    B4000   = 4000
    H_eps   = 3970.0
    K       = 3933.7
    H8      = 3889.1
    H9      = 3835.4
    H_alpha = 6562.8
 
    
    ## Plot absorption lines
    lines = np.array([Mg_5169, Mg_5175, Mg_5183, Fe_5332, Fe_5401, CaI_4227, H_beta, H_gamma, H_delta, B4000, H_eps, K, H8, H9, H_alpha])
    name_list = ['Mg_5169', 'Mg_5175', 'Mg_5183', 'Fe_5332', 'Fe_5401', 'CaI_4227', 'H_beta', 'H_gamma', 'H_delta', 'B4000', 'H_eps', 'K', 'H8', 'H9', 'H_alpha']
    for j in range(len(lines)):
        plt.vlines(lines[j],0,5,color='r')
        plt.text(lines[j],0,name_list[j],rotation=-90)

def plot_emission_lines():
    # Emissions linier
    OIII_3299 = 3299
    OIII_3312 = 3312
    NeV_3346  = 3346
    OII_3726  = 3726.2
    FeI_3750  = 3749.5
    HeI_4471  = 4471
    HeII_4686 = 4686
    OIII_4959 = 4959 # A
    OIII_5007 = 5007 # A
    

    
    lines_em = np.array([OIII_3299, OIII_3312, NeV_3346, OII_3726, FeI_3750 , HeI_4471, HeII_4686, OIII_4959, OIII_5007])
    name_list = ['OIII_3299', 'OIII_3312', 'NeV_3346', 'OII_3726', 'FeI_3750', 'HeI_4471', 'HeII_4686', 'OIII_4959', 'OIII_5007']
    for l in range(len(lines_em)):
       plt.vlines(lines_em[l],0,5,color='b')
       plt.text(lines_em[l],0,name_list[l],rotation=-90)




