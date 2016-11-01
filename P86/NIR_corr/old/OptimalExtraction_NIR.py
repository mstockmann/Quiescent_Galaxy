

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

    
    # P is the profile. The profile can be obtained by collapsing the onto the spatial direction (x)
    # Here we create a normalized profile
    P_N = P/sum(P)
    
    # F: 2D flux spectra that needs to be extracted
    # M: Bad pixel mask. Bad = 0, good = 1
    # w = 1/sig**2. Inverse squared error (from weighted average)

    # Softning factor
    eps = 0#1e-40

    f_1d = np.zeros(len(F[0,:]))
    e_1d = np.zeros(len(sig[0,:]))
    for i in range(len(f_1d)):
        F_i = F[:,i]
        M_i = M[:,i]
        w_i = P_N**2/(sig[:,i]**2)
        #w_i = 1/(sig[:,i]**2+eps)
        #w_i = w_i*(1e40)
    
    
        # If 0/0 happens we substitute with F=0, E=1.
        if all(M_i == 0) == True:
            #f_1d[i] = 0
            #e_1d[i] = 1
            f_1d[i] = (np.sum(F_i*w_i/P_N))/(np.sum(w_i))
            e_1d[i] = 1/np.sqrt(np.sum(w_i))
            
            
            #print 'heyehey %s' % i
        else:
            #f_1d[i] = (np.sum(M_i*F_i*P_N*w_i))/(np.sum(M_i*(P_N**2)*w_i))
            f_1d[i] = (np.sum(M_i*F_i*w_i/P_N))/(np.sum(M_i*w_i))
            
            e_1d[i] = 1/np.sqrt(np.sum(w_i))

    return f_1d,e_1d



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


N_spec = glob.glob('../../../X-shooter/P93/Data/Reduction/*')[:10]
Path_2d_spec = glob.glob('../../../X-shooter/P93/Data/Reduction/*/Combined_OBs/NIR_corr/*_V1_NIR_corr.fits')

klim_list = [[27, 48],[25, 47],[25, 47],[28, 45],[23, 48],[25, 50],[26, 49],[28, 49],[31, 49],[26, 50]]


rebin = 0 # 0: not rebinned, 1: rebinned
Norm_const = 1 # normalise the spectrum F_norm = F*Norm_const
for ii in range(len(N_spec)):
    print 'Target: %s' % Path_2d_spec[ii].split('/')[7]

    klim = klim_list[ii]

    # read in spectrum
    Flux, Err, BPM, Wave, hd_0, hd_1, hd_2 = read_in_2d_for_OptExt(Path_2d_spec[ii],Norm_const,rebin)

    # Visual confirm a straight trace
    _plot_trace_slope(Flux,Wave,100)


    # Visual confirm the PSF range
    PSF, r1, r2 = _plot_Wavelength_Collapsed_Spectrum(Flux,klim)


    # Selecting the psf of the data:
    Flux_psf = Flux[r1:r2,:]
    Err_psf  = Err[r1:r2,:]
    BPM_psf  = BPM[r1:r2,:]


    # Optimal Extract (Horne et al.)
    F_opt, E_opt = optext(BPM_psf,Flux_psf,PSF,Err_psf)



    # Remove unphysical S/N values
    S2N_after = F_opt/E_opt
    S2n_limit = 100

    ID = np.where(abs(S2N_after) > S2n_limit)[0]

    F_opt_s2ncorr = np.array(F_opt)
    E_opt_s2ncorr = np.array(E_opt)

    k0 = 0
    k1 = 0
    k2 = 0
    k3 = 0
    for i in range(len(ID)):
        if 5*np.median(F_opt) < abs(F_opt[ID[i]]):
            F_opt_s2ncorr[ID[i]] = np.median(F_opt)
            k0 += 1
            if F_opt_s2ncorr[ID[i]]/E_opt_s2ncorr[ID[i]] > S2n_limit:
                E_opt_s2ncorr[ID[i]] = np.median(E_opt)
                k3 += 1

        elif abs(E_opt[ID[i]]) < 1/5*np.median(E_opt):
            E_opt_s2ncorr[ID[i]] = np.median(E_opt)
            k1 += 1
        
        else:
            F_opt_s2ncorr[ID[i]] = np.median(F_opt)
            E_opt_s2ncorr[ID[i]] = np.median(E_opt)
            k2 += 1
    # print k0,k3,k1,k2, len(ID)


    F_opt = np.array(F_opt_s2ncorr)
    E_opt = np.array(E_opt_s2ncorr)



    F_opt = F_opt/Norm_const
    E_opt = E_opt/Norm_const

    # Read out the optimal extracted spectrum
    path_out = Path_2d_spec[ii].replace('.fits','_opt.fits')


    if not os.path.exists(path_out):
        # Read out flux array
        pf.writeto(path_out, F_opt, hd_0)
        
        # Read out error array
        pf.append(path_out, E_opt, hd_1)
    else:
        print 'file already exists'



    # # S2n_array = Flux_psf/Err_psf
    # # for i in range(len(S2n_array)):
    # #     plt.scatter(i,np.median(S2n_array[i,:]))

    # # S2n_opt_m = np.median(F_opt/E_opt)
    # # plt.scatter(10,S2n_opt_m,color='r')
    # # plt.show()

    # N_bin = 50
    # F_opt_bin = np.zeros(len(F_opt)/N_bin)
    # E_opt_bin = np.zeros(len(F_opt)/N_bin)
    # W_bin = np.zeros(len(Wave)/N_bin)
    # for i in range(int(len(F_opt)/N_bin)):
    #     F_opt_bin[i] = np.median(F_opt[i*N_bin:(i+1)*N_bin])
    #     W_bin[i] = Wave[(i+1)*N_bin/2]
    # plt.plot(W_bin,F_opt_bin)
    # plt.show()


    # plt.figure(2)
    # plt.plot(Wave[::10], F_opt[::10])
    # plt.show()



'''



# Test
S2N_after = F_opt/E_opt
ID = np.where(abs(S2N_after) > S2n_limit)[0]
print 'ID len: %s' % len(ID)
x = np.arange(len(S2N_after))

plt.plot(S2N_after,color='r')
plt.scatter(x[ID],S2N_after[ID],color='b')
plt.show()
'''




# # Plot

# ###
# plt.figure(figsize=(18,3))

# #rect_2D = [10,10,10,10]
# #rect_1D = [left_h, bottom, 0.17, height]

# plt.subplot(2,1,1)
# #plt.axes(rect_2D)

# # Filps the image?
# plt.imshow(Flux, vmin=-1e-19, vmax=5e-19,cmap='gray',origin='lower')
# #, aspect='auto'
# plt.subplot(2,1,2)
# plt.plot(Wave,F_opt)
# plt.axis([min(Wave),max(Wave),-0.01e-18,3e-18])
# plt.show()




# ### Remember to rescale

# F_opt = F_opt/Norm_const
# E_opt = E_opt/Norm_const

# # Read out the optimal extracted spectrum
# path_out = Path_2d_spec[i].replace('.fits','_opt.fits')


# if not os.path.exists(path_out):
#     # Read out flux array
#     pf.writeto(path_out, F_opt, hd_0)
    
#     # Read out error array
#     pf.append(path_out, E_opt, hd_1)
# else:
#     print 'file already exists'








