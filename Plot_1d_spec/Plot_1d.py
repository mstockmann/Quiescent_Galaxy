

#!/usr/local/bin/python

#-*- coding: utf-8 -*-

from __future__ import division
import glob
import pyfits as pf
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import os
from scipy.optimize import curve_fit
import sys

sys.path.insert(0, '..')
import Stokky as st

####################################################


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

def load_spec_z(path,obj):
    data = np.genfromtxt(path,dtype=str)
    Name = data[:,0]
    z_spec = data[:,1]
    loc_name = np.where(obj == Name)[0][0]
    return np.float(z_spec[loc_name])

####################
#### 2D spectra ####
####################

Path_2d_spec_new = glob.glob('../../X-shooter/cQGsample/Objects/*/NIR_corr/*_V1_NIRcorr_wmrebin15_opt.fits')


Path_2d_spec_old = glob.glob('../../X-shooter/cQGsample/Objects/*/NIR_corr/*_V1_NIRcorr_wmrebin15_opt_old.fits')
print len(Path_2d_spec_new), len(Path_2d_spec_old) 
raise


i = 0
target = Path_2d_spec[i].split('/')[-1].split('_')[0]
z_spec = load_spec_z('../../X-shooter/cQGsample/Collected_info/z_spec.txt','540713')

W_new, F_new, E_new, hdf_new, hde_new = read_in_1d_fits(Path_2d_spec_new[i])   
Wave, Flux, Errs, hdf_nir, hde_nir = read_in_1d_fits(Path_2d_spec_old[i])       

Wave_rf, F_rf, E_rf = st.convert_2_rest_frame(Wave, Flux, Errs, z_spec)





plt.title('%s' % target)

# plt.plot(Wave, Flux,color='black')
# plt.plot(Wave, Errs,color='purple', alpha=0.2)
# plt.plot(Wave,Flux/Errs)
plt.plot(Wave_rf, F_rf,color='black')
plt.plot(Wave_rf, E_rf,color='purple', alpha=0.2)


plt.show()







