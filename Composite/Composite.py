#!/usr/local/bin/python
#-*- coding: utf-8 -*-


from __future__ import division

__all__ = ["Composite cQG"]
__version__ = "0.0.0"
__author__ = "Mikkel Stockmann (mstockmann@dark-cosmology.dk)"
__copyright__ = "Copyright 2016 Mikkel Stockmann"


import numpy as np
import glob
# from scipy import interpolate
import matplotlib.pylab as pl
# from methods import latexify
# latexify()


import pyfits as pf

import matplotlib.pyplot as plt
from astropy.io import fits
import os
from scipy.optimize import curve_fit
import sys

sys.path.insert(0, '../')
from Stokky import *

####################################################

def convert_2_rest_frame(W,F,E,z):
	W_rest_frame = W / (1+z)
	F_rest_frame = F * (1+z)
	E_rest_frame = E * (1+z)
	return W_rest_frame, F_rest_frame, E_rest_frame

def wavelength_range(Warr):
	""" Takes 2d array and returns the min and max wavelength
	"""
	wlmin = []
	wlmax = []
	for jj in range(len(Warr)):
		wlmin.append(min(Warr[jj]))
		wlmax.append(max(Warr[jj]))
	return min(wlmin), max(wlmax)

def common_wavelength(wlarr_old, wlarr_new, fluxarr_old, fill_value = 0.):
	""" Jonathan Selsing
	"""
	from scipy import interpolate
	f = interpolate.interp1d(wlarr_old, fluxarr_old, kind='linear', bounds_error = False, fill_value=fill_value)
	fluxarr_new = f(wlarr_new)
	return fluxarr_new

def create_header(targets):
	string = 'Wavelength'
	for i in range(len(targets)):
		if i < 10:
			string += ' UV-%s' % targets[i]
		else:
			string += ' %s' % targets[i]
	return string


def rebin_1d(Npixels,wavelength,arr,arr1=None,type='median'):
	Binsize = int(len(arr)/Npixels)
	
	Arr_binned = np.zeros(Binsize)
	Wave_binned = wavelength[::Npixels][:-1]

	if type == 'median':
		for ii in range(Binsize):
			Arr_binned[ii] = np.median(arr[ii*Npixels:(ii+1)*Npixels]) 
		return Wave_binned,Arr_binned
	
	elif type == 'wmean':
		err_arr = np.zeros(Binsize)
		for ii in range(Binsize):
			val = arr[ii*Npixels:(ii+1)*Npixels]
			sig = arr1[ii*Npixels:(ii+1)*Npixels]
			
			w_err = 1/sig**2
			Arr_binned[ii], sow = np.average(val, weights=w_err,returned=True)
			err_arr[ii] = 1.0/np.sqrt(sow)
		return Wave_binned, Arr_binned, err_arr


####################################################

P93_spec = glob.glob('../../X-shooter/P93/Data/Reduction/*/Combined_OBs/NIR_corr/*V1_NIR_corr_opt.fits')

P86_spec = glob.glob('../../X-shooter/P86_COSMOS_1/Data/CP*/Combined_OBs/NIR_corr/*V1_NIR_corr_opt.fits')

Path = np.concatenate((P93_spec,P86_spec),axis=0)

### Set Parameters ###
interp_step = 0.05 # should be changed so it is CDELT1/4
n_obj = len(Path)




Wave_arr = []
Flux_arr = []
Errs_arr = []
target_list = []
for ii in range(n_obj):
	target_list.append(Path[ii].split('/')[-1].split('_')[0])

	Wave, Flux, Errs, hdf, hde = read_in_1d_fits(Path[ii])

	print hdf['CDELT1']*10



	z_spec=2.0 # Fit redshifts using ppxf
	W_rest, F_rest, E_rest = convert_2_rest_frame(Wave,Flux,Errs,z_spec)


	Wave_arr.append(W_rest)
	Flux_arr.append(F_rest)
	Errs_arr.append(E_rest)


wl_lower, wl_upper = wavelength_range(Wave_arr)
Wave_common = np.arange(wl_lower, wl_upper+interp_step, interp_step)

print Wave_common, wl_lower, wl_upper
raise

Flux_intp = np.zeros((n_obj,len(Wave_common)))
Errs_intp = np.zeros((n_obj,len(Wave_common)))
BPMs_intp = np.zeros((n_obj,len(Wave_common)))


for kk in range(n_obj):
	Flux_intp[kk] = common_wavelength(Wave_arr[kk], Wave_common, Flux_arr[kk])
	Errs_intp[kk] = common_wavelength(Wave_arr[kk], Wave_common, Errs_arr[kk],fill_value=1.0)
	BPMs_intp[kk] = np.ones(len(Wave_common))


# np.savetxt('data/regularised_flux.dat', np.column_stack((Wave_common, Flux_intp.transpose())), header=create_header(target_list))
# np.savetxt('data/regularised_errs.dat', np.column_stack((Wave_common, Errs_intp.transpose())), header=create_header(target_list))
# print('Saving IGM-corrected regularised data to to data/')

# for i in range(len(Path)):
# 	plt.plot(Wave_common,Flux_intp[i]+np.median(Flux_intp[i])*(i+200000))
# 	# plt.scatter(Wave_arr[i],Flux_arr[i],color='r',s=1)
# 	# plt.axis([3600,3700,min(Flux_intp[i]),max(Flux_intp[i])])
# plt.show()




# Make different combinations
wmean__flux = np.zeros(len(Wave_common))
wmean__errs = np.zeros(len(Wave_common))
mean___flux = np.zeros(len(Wave_common))
gmean__flux = np.zeros(len(Wave_common))
median_flux = np.zeros(len(Wave_common))

 
for i, k in enumerate(Flux_intp.transpose()):
	sigma = (Errs_intp.transpose()[i]*BPMs_intp.transpose()[i])
	val = (Flux_intp.transpose()[i]*BPMs_intp.transpose()[i])


	# Weighted average with std deviation errors
	w_err = 1/sigma**2
	wmean__flux[i], sow = np.average(val, axis=0, weights=w_err, returned=True)
	wmean__errs[i] = 1.0/np.sqrt(sow)

	# # Average
	# mean___flux = np.average(val, axis=0)

	# # geometric mean
	# from scipy.stats.mstats import gmean
	# mask = np.where(val == 0)
	# gmean__flux[i] = gmean(val[mask])

	# # median
	# median_flux[i] = np.median(val)






Wave_rebin, wmean_medrebin = rebin_1d(4*45,Wave_common,wmean__flux,type='median')
Wave_rebin, wmean_mw_rebin, err_mw_rebin = rebin_1d(4*45,Wave_common,wmean__flux,wmean__errs,type='wmean')

plt.plot(Wave_rebin,wmean_medrebin,'b')
plt.plot(Wave_rebin,wmean_mw_rebin,'r')
plt.show()























