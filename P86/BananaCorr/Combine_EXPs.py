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
from functions import *

######################
###### Programs ######
######################



###################################
###### Combine 2D EXPOSURES  ######
###################################




# Set parameters
## SCIENCE TARGET ##
# target = 'CP-540713'; roll_list = [[-5,0,5],[-5,0,5],[-5,0,5],[-5,0,5],[-5,0,5],[-5,0,5]]; psf_w1, psf_w2 = 30,40
# target = 'CP-561356'; roll_list = [[-5,5,0],[-5,5,0],[-5,5,0],[-5,5,0],[-5,5,0],[-5,5,0],[-5,5,0]]; psf_w1, psf_w2 = 27,43
# target = 'CP-1243752'; roll_list = [[-5,0,5],[-5,0,5],[0,5],[-5,5,0],[-5,5,0],[-5,5,0]]; psf_w1, psf_w2 = 27,43
# target = 'CP-1291751'; roll_list = [[-5,5,0],[-5,5,0],[-5,5,0],[-5,5,0],[-5,5,0],[-5,5,0]]; psf_w1, psf_w2 = 27,43




Nsigma = 5 # Sigma clipping factor



# General path
#root = '../../../X-shooter/P86_COSMOS_1/Data/%s/OB*/NIR_BananaCorr/sci_tellcorr_flux_merged2d_nir_%s_*.fits' % (target,target)
N_OB = len(glob.glob('../../../X-shooter/P86_COSMOS_1/Data/%s/OB*' % target))
Normalisation_const = 1e19


for kk in range(N_OB):
    path_OB = glob.glob('../../../X-shooter/P86_COSMOS_1/Data/%s/OB*' % target)[kk]
    path_EXP = glob.glob(path_OB+'/NIR_BananaCorr/sci_tellcorr*.fits')
    N_EXP = len(path_EXP)
    print 'N_EXP: %s' % N_EXP

    f_list = []
    for jj in range(N_EXP):

        ## Read in data
        exec('Flux_%s, Errs_%s, BPM_%s, Wave, hd_0, hd_1, hd_2 = read_fits_file(path_EXP[jj],Normalisation_const)' % (jj,jj,jj))

        ## Create list of array widths (x direction)
        exec('f_list.append(len(Flux_%s[:,0]))' % jj)

        ## Set all Errs=0 -> Errs=median() (median in wavelength direction)
        exec('Errs_%s = replace_zeros_with_median_lambda(Errs_%s)' % (jj,jj))

    ## Data array width equalizing (cut all to smallest width)
    flen_min = np.min(f_list)#-1
    
    for ii in range(N_EXP):
        exec('tmp%s = Flux_%s[:flen_min,:]' % (ii,ii))
        exec('emp%s = Errs_%s[:flen_min,:]' % (ii,ii))
        exec('bpm%s = BPM_%s[:flen_min,:]' % (ii,ii))



        ## Plot wavelength collapsed trace
        exec('plot_collapsed_trace_wavedirection(tmp%s,ii)' % ii)
    plt.vlines([psf_w1,psf_w2],-6,6,label='average psf: %s %s' % (psf_w1,psf_w2))  
    plt.title(path_OB.split('/')[-1]+' (None corrected)')
    plt.show()

    ## Align trace (using continuous rolls)
    for ii in range(N_EXP):
        exec('tmp%s = np.roll(tmp%s,roll_list[kk][ii],axis=0)' % (ii,ii))
        exec('emp%s = np.roll(emp%s,roll_list[kk][ii],axis=0)' % (ii,ii))
        exec('bpm%s = np.roll(bpm%s,roll_list[kk][ii],axis=0)' % (ii,ii))

    ## Plot the trace
    for ii in range(N_EXP):
        exec('plot_collapsed_trace_wavedirection(tmp%s,ii)' % ii)
    plt.vlines([psf_w1,psf_w2],-6,6,label='average psf: %s %s' % (psf_w1,psf_w2))
    plt.title(path_OB.split('/')[-1]+' (Corrected)')
    plt.show()  

    #######################
    ## Combine exposures ##

    # Make arrays
    # Create input (tmp0,tmp1,...,tmpN) and stack the 2d arrays, 
    merge_input_tmp = ','.join(['tmp%s' % i for i in range(N_EXP)])
    merge_input_emp = ','.join(['emp%s' % i for i in range(N_EXP)])
    merge_input_bpm = ','.join(['bpm%s' % i for i in range(N_EXP)])
    exec('tmp_merge = np.dstack((%s))' % merge_input_tmp)
    exec('emp_merge = np.dstack((%s))' % merge_input_emp)
    exec('bpm_merge = np.dstack((%s))' % merge_input_bpm)
    
    ## Median sigma clipping
    tmp_sigclip, emp_sigclip, bpm_sigclip = median_sigma_clipping(Nsigma,tmp_merge,emp_merge,bpm_merge,N_EXP)

    ## Remove bad pixel
    tmp_sig_bpc, emp_sig_bpc, Bpm_sig_bpc = remove_bad_pixels(tmp_sigclip,emp_sigclip,bpm_sigclip,N_EXP)

    ## Average weigthed mean
    F_avw_comb, E_avw_comb = calculate_average_weighted_mean(tmp_sig_bpc,emp_sig_bpc,Normalisation_const)

    ## Read out combined exposures
    path_out = path_OB+'/NIR_BananaCorr/%s_avwCombined_%sx%s_sig%s_V1.fits' % (target,N_EXP,path_OB.split('/')[-1],Nsigma)
    print 'read out: %s' % path_out
    read_out_to_fits(path_out,F_avw_comb,E_avw_comb,Bpm_sig_bpc,hd_0,hd_1,hd_2)




### Improvements: Check if trace tilts along the wavelength direction











