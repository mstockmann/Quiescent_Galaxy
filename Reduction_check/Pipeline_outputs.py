
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
import sys

sys.path.insert(0, '/Volumes/DataDrive/Quiescent_Galaxy')
import Stokky as st

#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@
#@#@#@# Science Frames #@#@#@#
#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@

target = 105842
# target = 108899
# target = 90676
# target = 155853
# target = 171060
# target = 171687
# target = 230929
# target = 239220
# target = 250513
# target = 773654
No_OB = len(glob.glob('../../X-shooter/P93/Data/Reduction/%s/OB*' % target))


# # Select the path to final reduction products
# path = glob.glob('../../X-shooter/P93/Data/Reduction/%s/OB*/NIR_BananaCorr/Reduction/Output/xsh_scired_slit_nod_SCI_SLIT_FLUX_MERGE2D_NIR_NIR_x_.fits' % target)

# print '%s, Number of exposures / Finished reductions: %s / %s' % (target,len(path),No_OB)
# # Visual inspection of spectra in gaia
# for i in range(len(path)):
#     os.system("gaia %s" % (path[i]))





# predict
path_predict = glob.glob('../../X-shooter/P93/Data/Reduction/%s/OB*/NIR_BananaCorr/Reduction/Output/xsh_predict_FMTCHK_RESID_TAB_LINES_NIR_NIR_x_.fits' % target)


print 'PREDICT: %s. Number of exposures / predict output: %s / %s' % (target,len(path_predict),No_OB)



path_arc_lamp_pinhole = glob.glob('../../X-shooter/P93/Data/Reduction/%s/OB*/data_NIR/XSHOO.2014-04-13T11:14:01.130.fits' % target)

# print path_arc_lamp_pinhole
# raise

i = 0


g = pf.open(path_arc_lamp_pinhole[i])
print g.info()

Counts = g[0].data
# raise




f = pf.open(path_predict[i])

Wavelength = np.array((f[1].data)['Wavelength'])
ResXmodel  = np.array((f[1].data)['ResidXmodel'])
ResYmodel  = np.array((f[1].data)['ResidYmodel'])
Flag       = np.array((f[1].data)['Flag'])


plt.figure(figsize=(14,6))
plt.subplot(2,2,1)
plt.imshow(Counts,vmin=40,vmax=450,cmap='gnuplot2')
plt.axis([0,2000,0,1000])


plt.subplot(2,2,3)
plt.scatter(Wavelength,ResXmodel,color='b',s=5)
plt.plot(Wavelength,Flag,color='r')
plt.title('Line X residuals')
plt.ylabel('Residuals [pix]')
plt.xlabel('Wavelength [nm]')

plt.subplot(2,2,4)
plt.scatter(Wavelength,ResYmodel,color='b',s=1)
plt.plot(Wavelength,Flag,color='r')
plt.title('Line Y residuals')
plt.ylabel('Residuals [pix]')
plt.xlabel('Wavelength [nm]')

plt.show()

# print Wavelength

# Wave, Arr1, Arr2, hd1, hd2 = st.read_in_1d_fits(path_predict[i])

# plt.plot(Wave,Arr1)
# plt.show()


# orderpos
# mflat
# 2dmap
# wavecal
# flexcomp
# response














