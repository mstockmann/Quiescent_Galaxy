#!/usr/local/bin/python

#-*- coding: utf-8 -*-


from __future__ import division
import glob
# import pyfits as pf
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import mpl_toolkits.axisartist as AA
import sys
from scipy import interpolate

sys.path.insert(0, '../../..')
import Stokky as st



path = glob.glob('../../../../X-shooter/cQGsample/Objects/108899/NIR_corr/108899_avwCombined_OBx5_sig5_V1_NIRcorr_wmrebin15_opt.fits')

W, F, E, M, hdf, hde, hdm = st.read_in_1d_fits(path[0])

Arr = np.zeros(shape=(len(W),4))

Arr[:,0] = W
Arr[:,1] = F
Arr[:,2] = E
Arr[:,3] = M

np.savetxt('108899_avwCombined_OBx5_sig5_V1_NIRcorr_wmrebin15_opt.dat',Arr,delimiter=' ',fmt='%g')

