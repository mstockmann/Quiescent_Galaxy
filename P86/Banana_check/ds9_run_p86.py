#!/usr/local/bin/python

#-*- coding: utf-8 -*-

from __future__ import division
import glob
# import pyfits as pf
import numpy as np
import matplotlib.pyplot as plt
# from astropy.io import fits
import os

########################
########################

period = 'P86_COSMOS_1'

Path_nir = glob.glob('../../../X-shooter/%s/Data/CP*/Combined_OBs/NIR/*_V1_wmrebin15.fits' % (period))

Path_nir_banana = glob.glob('/Volumes/DataDrive/X-shooter/%s/Data/CP*/Combined_OBs/NIR_corr/*_V1_NIRcorr_wmrebin15.fits' % (period))


# print len(Path_nir)
# print len(Path_nir_banana)
# raise

for i in range(len(Path_nir_banana)):
    print Path_nir[i].split('/')[6]
    os.system('ds9 -zscale -cmap bb -frame lock wcs %s %s' % (Path_nir[i],Path_nir_banana[i]))






