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


###########################################


path_spec = glob.glob('../../X-shooter/cQGsample/Spectra_analysis/3-Stdcorr_slit_aper/NIR/*_V1_NIRcorr_wmrebin15_er2d_opt_stdcorr.fits')

# print path_spec

for i in range(len(path_spec)):
	target = path_spec[i].split('/')[-1].split('_')[0]
	W, F, E, M, hdf, hde, hdm = st.read_in_1d_fits(path_spec[i])

	Arr = np.zeros(shape=(len(W),4))
	Arr[:,0] = W
	Arr[:,1] = F
	Arr[:,2] = E
	Arr[:,3] = M

	outfile = path_spec[i].replace('.fits','.dat').replace(target,'_%s' % target)
	np.savetxt(outfile,Arr,fmt='%.10g')
	# raise