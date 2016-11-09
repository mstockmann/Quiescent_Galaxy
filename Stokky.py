
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
	hde = data_arr[0].header
	F = data_arr[0].data
	E = data_arr[1].data
	W = (hdf['CRVAL1'] + (hdf['CRPIX1'] - 1 + np.arange(hdf['NAXIS1']))*hdf['CDELT1'])*10 
	return W, F, E, hdf, hde




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


