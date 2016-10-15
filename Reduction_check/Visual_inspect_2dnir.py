
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
#@#@#@# Science Frames #@#@#@#
#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@

# target = 105842
# target = 108899
# target = 90676
# target = 155853
# target = 171060
# target = 171687
# target = 230929
# target = 239220
# target = 250513
# target = 773654


# Select the path to final reduction products
path = glob.glob('../../X-shooter/P93/Data/Reduction/%s/OB*/NIR_BananaCorr/Reduction/Output/xsh_scired_slit_nod_SCI_SLIT_FLUX_MERGE2D_NIR_NIR_x_.fits' % target)
No_OB = len(glob.glob('../../X-shooter/P93/Data/Reduction/%s/OB*' % target))



print '%s, Number of exposures / Finished reductions: %s / %s' % (target,len(path),No_OB)
# Visual inspection of spectra in gaia
for i in range(len(path)):
    os.system("gaia %s" % (path[i]))











