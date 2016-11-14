

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

sys.path.insert(0, '../..')
import Stokky as st

####################################################



############################################
#### mask offset emission in 2d spectra ####


Path_2d_spec = glob.glob('../../../X-shooter/cQGsample/Objects/*/NIR_corr/*_V1_NIRcorr_wmrebin15.fits')



# Check spectra
# for ii in range(len(Path_2d_spec)):
#     os.system('gaia %s' % Path_2d_spec[ii])   


### pixel mask [x,y] ###

# # 105842
# target_no = 0
# mask_array = [[[522,524],[554,556],[570,573],[1091,1093]],[[41,46],[41,46],[40,46],[41,49]]] 

# # 108899
# target_no = 1
# mask_array = [[],[]]

# # 155853
# target_no = 2
# mask_array = [[[765,767]],[[39,43]]] 

# # 171060
# target_no = 3
# mask_array = [[],[]]

# # 171687
# target_no = 4
# mask_array = [[[242,248],[723,726],[730,733]],[[28,32],[30,34],[30,34]]]

# # 230929
# target_no = 5
# mask_array = [[],[]]

# # 239220
# target_no = 6
# mask_array = [[],[]]

# # 250513
# target_no = 7
# mask_array = [[],[]]

# 773654
# target_no = 8
# mask_array = [[],[]]

# 90676
# target_no = 9
# mask_array = [[],[]]

# CP-1243752
# target_no = 10
# mask_array = [[],[]]

# CP-1291751
# target_no = 11
# mask_array = [[],[]]

# CP-540713
# target_no = 12
# mask_array = [[],[]]

# CP-561356
# target_no = 13
# mask_array = [[],[]]


Spec_no = target_no
F, E, M, hdf, hde, hdm = st.read_2d_fits(Path_2d_spec[Spec_no])
# print F.shape
# raise

for jj in range(len(mask_array[0])):
    x_tmp = mask_array[0][jj]
    y_tmp = mask_array[1][jj]

    F[y_tmp[0]:y_tmp[1],x_tmp[0]:x_tmp[1]] = np.median(F[y_tmp[0]:y_tmp[1],x_tmp[0]-20:x_tmp[1]+20])



target = Path_2d_spec[Spec_no].split('/')[-1].split('_')[0]
name_out = target+'_'+'_'.join(Path_2d_spec[Spec_no].split('/')[-1].split('_')[4:])
# print name_out.replace('.fits','_2der.fits')

# .replace('.fits','_2der.fits').replace('_avwCombined_OBx*_sig5','')
path_out = '../../../X-shooter/cQGsample/Spectra_analysis/1-2drebin_emission_rm/%s' % (name_out.replace('.fits','_er2d.fits'))
pf.writeto(path_out, F, hdf)

# Read out error array
pf.append(path_out, E, hde)

# Read out bad pixel map
pf.append(path_out, M, hdm)






## Which spectra have off-set emission

# UV-105842 - have offset emission (x=8322;8340, y=41;44) (x=8555;8584,y=39;45) +(x=16367;16401, y=41;44)
# UV-108899 - HST F814 also indicate bright pixels close to center
#             This cannot be seen in 2d spectrum
# UV-155853 - Flag the bad pixel at 16829 A
# UV-171060 - Flag; 15425, 
# UV-171687 - offset emission, mask
# UV-230929
# UV-239220 - noise lines offset from trace. The trace is quite wide though
# UV-250513 - Double trace centered. Looks quiescent in both cases
# UV-773654
# UV-90676 - ind: see nothing

# CP-1243752 - offset noise lines
# CP-1291751 - offset noise lines
# CP-540713 - offset noise lines
# CP-561356











