#!/usr/local/bin/python

#-*- coding: utf-8 -*-

from __future__ import division
import glob
import pyfits as pf
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import os
from functions import *



def read_2d_fits(path):
    f = pf.open(path)

    hd1 = f[0].header
    hd2 = f[1].header
    hd3 = f[2].header
    Arr1 = f[0].data
    Arr2 = f[1].data
    Arr3 = f[2].data
    return Arr1, Arr2, Arr3, hd1, hd2, hd3


################################################################


path_correct_errs = glob.glob('../../../X-shooter/P93/Data/Reduction/*/Combined_OBs/NIR/*_avwCombined_OBx*_sig5_V3.fits')
path_correct_flux = glob.glob('../../../X-shooter/P93/Data/Reduction/*/Combined_OBs/NIR_Banana/*_NIRBanana_avwCombined_tellcorr_*xOB_sig5_V1.fits')

print 'No of Error / Flux spectra: %s / %s' % (len(path_correct_errs), len(path_correct_flux))


for i in range(len(path_correct_errs)):
	obj = path_correct_errs[i].split('/')[7]
	filename = path_correct_errs[i].split('/')[-1]


	F_e, E_e, Q_e, hd1_e, hd2_e, hd3_e = read_2d_fits(path_correct_errs[i])
	F_f, E_f, Q_f, hd1_f, hd2_f, hd3_f = read_2d_fits(path_correct_flux[i])


	###
	print 
	print F_e.shape, E_e.shape, Q_e.shape
	print F_f.shape, E_f.shape, Q_f.shape

	E_e = np.roll(E_e,1,axis=1)
	F_f = F_f[:-1,:-1]
	E_f = E_f[:-1,:-1]
	Q_f = Q_f[:-1,:-1]	

	print '...'
	print F_e.shape, E_e.shape, Q_e.shape
	print F_f.shape, E_f.shape, Q_f.shape
	print 
	print 
	###


	path_out = '../../../X-shooter/P93/Data/Reduction/%s/Combined_OBs/NIR_corr' % obj
	if not os.path.exists(path_out):
		os.system('mkdir %s' % path_out)
	path_out_name = path_out+'/'+filename.replace('V3','V1_NIR_corr')
	print path_out_name
	read_out_to_fits(path_out_name,F_f,E_e,Q_f,hd1_f,hd2_f,hd3_f)




