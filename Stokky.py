
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




