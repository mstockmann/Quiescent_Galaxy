
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









