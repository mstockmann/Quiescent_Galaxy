#!/usr/local/bin/python

#-*- coding: utf-8 -*-

from __future__ import division
import glob
import pyfits as pf
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import os



##################
### Functions ####



def read_fits_file(path,norm_const):
    f = pf.open(path)
    hd_0 = f[0].header
    hd_1 = f[1].header
    hd_2 = f[2].header

    F = f[0].data*norm_const
    E = f[1].data*norm_const
    B = f[2].data
    W = 10*(hd_0['CRVAL1'] + (hd_0['CRPIX1'] - 1 + np.arange(hd_0['NAXIS1']))*hd_0['CDELT1'])

    return F, E, B, W, hd_0, hd_1, hd_2


def replace_zeros_with_median_lambda(Arr):
    Arr_new = np.array(Arr)
    ID_zero_arr = np.where(Arr == 0)
    Arr_med_wave = np.median(Arr,axis=1)
    Arr_med_wave[np.where(Arr_med_wave == 0)] = np.median(Arr_med_wave)

    for jj in range(len(ID_zero_arr[0])):
        index1 = ID_zero_arr[0][jj]
        index2 = ID_zero_arr[1][jj]
        Arr_new[index1,index2] = Arr_med_wave[index1]
    return Arr_new



def plot_collapsed_trace_wavedirection(data,index):
    min_width = len(data[:,0])
    Fplot = np.zeros(min_width)
    for jj in range(min_width):
        Fplot[jj] = np.median(data[jj,:])
    plt.plot(Fplot,label='OB %s' % (index+1))
    plt.legend(loc='upper left')


def median_sigma_clipping(sigma,flux,err,bpm,nexp):

    flux_median = np.median(flux,axis=2)
    nexp_flux_median = ','.join(['flux_median' for i in range(nexp)])    
    exec('flux_median_arr = np.dstack((%s))' % nexp_flux_median)
    Flux_dist_from_median = abs(flux_median_arr-flux)
    sigma_limit = sigma*err

    ID_sigmaclip_pixel = np.where(sigma_limit<Flux_dist_from_median)

    # Define new arrays
    flux_sig = np.array(flux)
    err_sig = np.array(err)
    bpm_sig = np.array(bpm)

    # Clip away pixel larger than Nsigma by setting their constribution to 0
    # The bpm is set to a arbitrary figure > 0, (cannot be 0)                  
    flux_sig[ID_sigmaclip_pixel] = 0
    err_sig[ID_sigmaclip_pixel]  = 0
    bpm_sig[ID_sigmaclip_pixel]  = 5

    return flux_sig, err_sig, bpm_sig



def remove_bad_pixels(flux,errs,bpm,nexp):
    """ Replace all good pixels with 1, and all bad with 0. We then multiply bpm with 
        flux and err and this removes the bad pixels.

        return: flux, error and a updated bad pixel map where combination with only bad pixel
                have gotten the number 22.
    """
    # Swap zeros and ones
    bpm[bpm==0] = 1
    bpm[bpm>1]  = 0

    # Replace zeros with ones if all combinations of one pixel are bad pixel 
    ID_only_zero = (np.sum(bpm,axis=2) == 0)
    bpm[ID_only_zero] = [1 for i in range(nexp)] # similar to [1,1,1] for nexp=3

    # Correct for bad pixel
    flux_bpc = flux*bpm
    errs_bpc = errs*bpm

    # Update bpm to output format:
    # 0: Good, 11 x nexp: N_EXP bad combinations, 100 x nexp: all bad combinations
    bpm[ID_only_zero] = [100 for i in range(nexp)]
    bpm[bpm==0] = 11
    bpm[bpm==1] = 0
    bpm_final = np.sum(bpm,axis=2)

    return flux_bpc, errs_bpc, bpm_final


def calculate_average_weighted_mean(flux,err,norm_const):
    """ Calculate the statistical average weighted mean 
        
        All errors equal to zero will have an infinite contribution 
        and are controlled by assigning a specific value -9999.99 
        that is, when the weight have been calculated, substituted
        with zero.

        return: average weighted mean and error
    """
    # Calculate the weights and make sure that 1/0 is not happening
    Id_value = -9999.99
    err[err==0] = Id_value

    weight = 1/err**2
    weight_Idvalue = 1/(Id_value)**2
    weight[weight == weight_Idvalue] = 0


    # Substitute zero in  flux*weights for median along the wavelength direction
    weight_2d = np.sum(weight,axis=2)
    weight_2d_nozeros = replace_zeros_with_median_lambda(weight_2d)

    # Substitute zero in flux*weights for median along the wavelength direction
    flux_n_weight_2d = np.sum(flux*weight,axis=2)
    flux_n_weight_2d_nozeros = replace_zeros_with_median_lambda(flux_n_weight_2d)

    # Average weighted mean
    F_ac = flux_n_weight_2d_nozeros/weight_2d_nozeros
    E_ac = 1/np.sqrt(weight_2d_nozeros)

    return F_ac/norm_const, E_ac/norm_const


def read_out_to_fits(path,F,E,B,hd0,hd1,hd2):
    if not os.path.exists(path):
	    #Read out flux array
	    pf.writeto(path,F,hd0)

	    # Read out error array
	    pf.append(path,E,hd1)

	    # Read out bad pixel map
	    pf.append(path,B,hd2)

	    # S/N spectrum
	    # hd3 = hd2
	    # hd3['EXTNAME'] = 'S2N'
	    # S2N = F/E
	    # pf.append(path, S2N, hd3)
    else:
        print 'file already exists'
