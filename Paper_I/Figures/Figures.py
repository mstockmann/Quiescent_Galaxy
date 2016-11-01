#!/usr/local/bin/python

#-*- coding: utf-8 -*-

from __future__ import division
import glob
import pyfits as pf
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from astropy.io import fits
import os
from scipy.optimize import curve_fit
import sys

sys.path.insert(0, '../../../Quiescent_Galaxy')
import Stokky as st


def cm2inch(value):
    return value/2.54

#############################################

border_thickness = 0.5
labelsize = 8
axissize = 7
f_size = 7

mpl.rcParams['axes.linewidth'] = border_thickness #set the value globally

mpl.rc('xtick', labelsize=axissize) 
mpl.rc('ytick', labelsize=axissize)
plt.rcParams['xtick.major.pad']='3'
plt.rcParams['ytick.major.pad']='3'

plt.rc('axes', labelsize=10)
plt.rc('font',family='sans serif',style='normal', variant='normal', stretch='normal', size='14')
plt.rc('mathtext',fontset='stixsans')

###############################################


################
### Figure 2 ###


path_2d = glob.glob('../../../X-shooter/P93/Data/Reduction/*/Combined_OBs/NIR_corr/*_V1_NIRcorr_wmrebin15.fits')

path_1d = glob.glob('../../../X-shooter/P93/Data/Reduction/*/Combined_OBs/NIR_corr/*_V1_NIR_corr_opt.fits')


print path_2d



fig_width  = 18
fig_height = 26
plt.figure(figsize=(cm2inch(fig_width), cm2inch(fig_height)))
plt.subplots_adjust(hspace=0.00)

for ii in range(len(path_2d)):


	# 2d spectrum NIR
	Arr1, Arr2, Arr3, hd1, hd2, hd3 = st.read_2d_fits(path_2d[ii])
	zmin, zmax = -5e-20, 3e-19
	r1, r2 = 25, 50


	# 1d spectrum NIR
	W, F, E, hdf, hde = st.read_in_1d_fits(path_1d[ii])


	plt.subplot(len(path_2d)*2,1,ii*2+1)
	# plt.subplot(20,1,1)
	plt.imshow(Arr1[r1:r2,:], cmap='gray', vmin=zmin, vmax=zmax,origin='lower',aspect='auto')

	plt.subplot(len(path_2d)*2,1,ii*2+2)
	# plt.subplot(20,1,2)
	plt.plot(W,F)
	plt.xlim([min(W),max(W)])


plt.show()















# plt.tight_layout(pad=0.0, w_pad=-0.2, h_pad=0)

# outname = 'Fig2.pdf'
# plt.savefig(outname)
# os.system('pdfcrop %s' % outname)
# os.system('rm %s' % outname)
# os.system('mv %s-crop.pdf %s' % (outname.split('.pdf')[0],outname))


################
################