

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


####################################################


################################
##### Function definitions #####
################################

def read_in_2d_for_OptExt(path,norm_param=1,rebin=0):
    f = pf.open(path)
    F = f[0].data*norm_param; hd1 = f[0].header
    E = f[1].data*norm_param; hd2 = f[1].header
    M = f[2].data; hd3 = f[2].header
    W = (hd1['CRVAL1'] + (hd1['CRPIX1'] - 1 + np.arange(hd1['NAXIS1']))*hd1['CDELT1'])*10 #
    return F, E, M, W, hd1, hd2, hd3


def optext(M,F,P,sig):
    ''' Input
        M: Bad Pixel Map (Good = 1, Bad = 0)
        F: Flux spectrum
        P: Psf profile (made by collapsing the spectrum along the wavelength direction)
        sig: Error spectrum

        Output: Optimal Extracted F, E, BPM (???)
    '''
    

    # Normalized profile
    P_N = P/sum(P)
    P_N[P_N == 0.0] = np.min(P_N[P_N>0])


    f_1d = np.zeros(len(F[0,:]))
    e_1d = np.zeros(len(F[0,:]))
    m_1d = np.zeros(len(F[0,:]))
    for i in range(len(f_1d)):
        F_i = F[:,i]
        M_i = M[:,i]
        w_i = P_N**2/(sig[:,i]**2)

        Number_of_bad_pixels = abs(len(M_i) - np.count_nonzero(M_i))
        # print '%s/%s' % (Number_of_bad_pixels,len(M_i))
        if Number_of_bad_pixels > len(M_i)/3:
            f_1d[i] = (np.sum(F_i*w_i/P_N))/(np.sum(w_i))
            e_1d[i] = 1/np.sqrt(np.sum(w_i))
            m_1d[i] = 0
        else:
            f_1d[i] = (np.sum(M_i*F_i*w_i/P_N))/(np.sum(M_i*w_i))
            e_1d[i] = 1/np.sqrt(np.sum(w_i))
            m_1d[i] = 1

    return f_1d, e_1d, m_1d



def weighted_mean(F,E):
    w = 1/(E**2)
    F_wm = np.sum(w*F)/np.sum(w)
    E_wm = 1/np.sqrt(np.sum(w))
    return F_wm,E_wm


def _plot_trace_slope(F,W,binsize):
    x_bin = int(len(F[0,:])/binsize)
    Y_max = np.zeros(binsize)
    X = np.zeros(binsize)
    for jj in range(len(Y_max)):
        tmp_hor_med = np.zeros(len(F[:,0]))
        for ii in range(len(tmp_hor_med)):
            tmp_hor_med[ii] = np.median(F[ii,jj*x_bin:(jj+1)*x_bin])

        Y_max[jj] = np.where(np.max(tmp_hor_med)==tmp_hor_med)[0][0]
        X[jj] = W[(jj+1)*x_bin/2]

    plt.scatter(X,Y_max,color='r')
    plt.hlines(np.median(Y_max),0,np.max(X)*3,label='%s' % np.median(Y_max))
    plt.axis([min(X)-500,max(X)+500,0,len(F[:,0])])
    plt.legend(loc='upper right')
    plt.show()


def _plot_Wavelength_Collapsed_Spectrum(F,psf):
    # Horizontal collapse of Flux image
    Flux_hor_med = np.zeros(len(F[:,0]))

    for i in range(len(F[:,0])):
        Flux_hor_med[i] = np.median(F[i,:])
    psf_points = psf[0], psf[1]

    plt.figure(1)
    plt.vlines(psf_points,-Norm_const,Norm_const,label='average psf: %s %s' % (np.round(psf_points[0],3),np.round(psf_points[1],3)))
    plt.plot(Flux_hor_med)
    plt.legend(loc='upper right')
    plt.hlines(0,-1,80,linestyle='--')
    plt.axis([-1,80+1,-2e-19*Norm_const,2.5e-19*Norm_const])
    plt.show()
    return Flux_hor_med[psf[0]:psf[1]], psf[0], psf[1]


#############################################
##### Optimal Extraction (Horne et al.) #####
#############################################


####################
#### 2D spectra ####
####################


N_spec = glob.glob('../../X-shooter/cQGsample/Objects/*')
#Path_2d_spec = glob.glob('../../../X-shooter/P93/Data/Reduction/*/Combined_OBs/NIR_corr/*_V1_NIR_corr.fits')

Path_2d_spec = glob.glob('../../X-shooter/cQGsample/Objects/*/NIR_corr/*_V1_NIRcorr_wmrebin15.fits')

klim_list = [[27, 48],[25, 48],[25, 47],[28, 45],[23, 48],[25, 50],[26, 49],[28, 49],[31, 49],[26, 50],[27, 44],[27, 44],[27, 44],[27, 44]]




rebin = 0 # 0: not rebinned, 1: rebinned
Norm_const = 1 # normalise the spectrum F_norm = F*Norm_const
for ii in range(len(N_spec)):
    # ii=1
    print 'Target: %s' % Path_2d_spec[ii].split('/')[6]

    klim = klim_list[ii]

    # read in spectrum
    Flux, Err, BPM, Wave, hd_0, hd_1, hd_2 = read_in_2d_for_OptExt(Path_2d_spec[ii],Norm_const,rebin)

    # Visual confirm a straight trace
    _plot_trace_slope(Flux,Wave,100)


    # Visual confirm the PSF range
    PSF, r1, r2 = _plot_Wavelength_Collapsed_Spectrum(Flux,klim)

    Flux_psf = Flux[r1:r2,:]
    Err_psf  = Err[r1:r2,:]
    BPM_psf  = BPM[r1:r2,:]

    # Optimal Extract (Horne et al.)
    F_opt, E_opt, M_opt = optext(BPM_psf,Flux_psf,PSF,Err_psf)


    plt.plot(Wave,F_opt,'b')
    plt.show()


    F_opt = F_opt/Norm_const
    E_opt = E_opt/Norm_const

    # Read out the optimal extracted spectrum
    path_out = Path_2d_spec[ii].replace('.fits','_opt.fits')
    print path_out

    if not os.path.exists(path_out):
        # Read out
        pf.writeto(path_out, F_opt, hd_0)
        pf.append(path_out, E_opt, hd_1)
        pf.append(path_out, M_opt, hd_2)
    else:
        os.system('mv %s %s' % (path_out,path_out.replace('.fits','_old.fits')))
        # Read out
        pf.writeto(path_out, F_opt, hd_0)
        pf.append(path_out, E_opt, hd_1)
        pf.append(path_out, M_opt, hd_2)

        print 'file already exists'
    # raise







