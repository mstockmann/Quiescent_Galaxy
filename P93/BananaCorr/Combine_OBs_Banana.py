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


## SCIENCE TARGET ##

# target = '105842'; roll_list = [0,-1,-2,0]
# target = '108899'; roll_list = [-1,0,-1,-3,0]
# target = '155853'; roll_list = [-1,0,0,-3,0]
# target = '171060'; roll_list = [1,1,-2,-2,-1]
# target = '171687'; roll_list = [-3,-3,1,0,0]
# target = '230929'; roll_list = [0,-1,2,0,-2,0]
# target = '239220'; roll_list = [0,1,2,-1,-1]
# target = '250513'; roll_list = [-1,1,0,0,0]
# target = '773654'; roll_list = [1,1,-1,1,3]
target = '90676'; roll_list = [0,2,-1,-1,-2]


Nsigma = 5 # Sigma clipping factor

						
OB_paths = glob.glob('../../../X-shooter/P93/Data/Reduction/%s/OB*/NIR_BananaCorr/sci_tellcorr_flux_merged2d_banana_nir_%s_OB*' % (target,target))


N_OB = len(OB_paths)
Normalisation_const = 1e19

f_list = []
OB_name = []
for jj in range(N_OB):
	OB_name.append(OB_paths[jj].split('/')[-1].split('_')[2].split('x')[-1])

	## Read in data
	exec('Flux_%s, Errs_%s, BPM_%s, Wave, hd_0, hd_1, hd_2 = read_fits_file(OB_paths[jj],Normalisation_const)' % (jj,jj,jj))

	## Create list of array widths (x direction)
	exec('f_list.append(len(Flux_%s[:,0]))' % jj)

	## Set all Errs=0 -> Errs=median() (median in wavelength direction)
	exec('Errs_%s = replace_zeros_with_median_lambda(Errs_%s)' % (jj,jj))


## Data array width equalizing (cut all to smallest width)
flen_min = np.min(f_list)
min_len = np.min([len(Flux_0[0,:]), len(Errs_0[0,:]), len(BPM_0[0,:])])


for ii in range(N_OB):
	exec('tmp%s = Flux_%s[:flen_min,:min_len]' % (ii,ii))
	exec('emp%s = Errs_%s[:flen_min,:min_len]' % (ii,ii))
	exec('bpm%s = BPM_%s[:flen_min,:min_len]' % (ii,ii))

	## Plot wavelength collapsed trace
	exec('plot_collapsed_trace_wavedirection(tmp%s,ii)' % ii) 
plt.title(OB_name[jj].split('/')[-1]+' (None corrected)')
plt.show()

## Align trace (using continuous rolls)
for ii in range(N_OB):
	exec('tmp%s = np.roll(tmp%s,roll_list[ii],axis=0)' % (ii,ii))
	exec('emp%s = np.roll(emp%s,roll_list[ii],axis=0)' % (ii,ii))
	exec('bpm%s = np.roll(bpm%s,roll_list[ii],axis=0)' % (ii,ii))

## Plot the trace
for ii in range(N_OB):
	exec('plot_collapsed_trace_wavedirection(tmp%s,ii)' % ii)
plt.title(OB_name[jj].split('/')[-1]+' (Corrected)')
plt.show()  




#######################
## Combine exposures ##

# Make arrays
# Create input (tmp0,tmp1,...,tmpN) and stack the 2d arrays, 
merge_input_tmp = ','.join(['tmp%s' % i for i in range(N_OB)])
merge_input_emp = ','.join(['emp%s' % i for i in range(N_OB)])
merge_input_bpm = ','.join(['bpm%s' % i for i in range(N_OB)])
exec('tmp_merge = np.dstack((%s))' % merge_input_tmp)
exec('emp_merge = np.dstack((%s))' % merge_input_emp)
exec('bpm_merge = np.dstack((%s))' % merge_input_bpm)

## Median sigma clipping
tmp_sigclip, emp_sigclip, bpm_sigclip = median_sigma_clipping(Nsigma,tmp_merge,emp_merge,bpm_merge,N_OB)

## Remove bad pixel
tmp_sig_bpc, emp_sig_bpc, Bpm_sig_bpc = remove_bad_pixels(tmp_sigclip,emp_sigclip,bpm_sigclip,N_OB)

## Average weigthed mean
F_avw_comb, E_avw_comb = calculate_average_weighted_mean(tmp_sig_bpc,emp_sig_bpc,Normalisation_const)

## Read out combined exposures
path_out = '../../../X-shooter/P93/Data/Reduction/%s/Combined_OBs/NIR_Banana/%s_NIRBanana_avwCombined_tellcorr_%sxOB_sig%s_V1.fits' % (target,target,N_OB,Nsigma)


print 'read out: %s' % path_out
read_out_to_fits(path_out,F_avw_comb,E_avw_comb,Bpm_sig_bpc,hd_0,hd_1,hd_2)






