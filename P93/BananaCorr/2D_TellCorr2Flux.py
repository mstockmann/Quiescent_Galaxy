#!/usr/local/bin/python

#-*- coding: utf-8 -*-

from __future__ import division
import glob
import pyfits as pf
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import os

########################


def read_2d_fits(path):
    f = pf.open(path)

    hd1 = f[0].header
    hd2 = f[1].header
    hd3 = f[2].header
    Arr1 = f[0].data
    Arr2 = f[1].data
    Arr3 = f[2].data
    return Arr1, Arr2, Arr3, hd1, hd2, hd3



def Flux_Telluric_Correction(F,E,t):
    F_t = np.zeros(shape=(len(F[:,0]),len(F[0,:])))
    E_t  = np.zeros(shape=(len(E[:,0]),len(E[0,:])))

    for j in range(len(F[:,0])):
        F_t[j,:] = F[j,:]/t
        E_t[j,:] = E[j,:]/t
    return F_t, E_t


def read_out_tellcorr_arr(path,target,F,E,Q,hd1,hd2,hd3):
    path_out = '/'.join(path.split('/')[:10])+'/sci_tellcorr_flux_merged2d_banana_nir_%s_%s.fits' % (target,path.split('/')[8])

    if not os.path.exists(path_out):
        pf.writeto(path_out, F, hd1)
        pf.append(path_out, E, hd2)
        pf.append(path_out, Q, hd3)
    else:
        os.system("rm %s" % (path_out))
        pf.writeto(path_out, F, hd1)
        pf.append(path_out, E, hd2)
        pf.append(path_out, Q, hd3)



########################


# Set parameters
# target = '105842' # 4 OB'S
# target = '108899' # 5 OB's
# target = '90676' # 5 OB's
# target = '155853' # 5 OB's
# target = '171060' # 5 OB's
# target = '171687' # 5 OB's
# target = '230929' # 5 OB's
# target = '239220' # 5 OB's
# target = '250513' # 5 OB's
target = 773654 # 5 OB's



root1 = '../../../X-shooter/P93/Data/Reduction/%s/OB*/NIR_BananaCorr/Reduction/Output/xsh_scired_slit_nod_SCI_SLIT_FLUX_MERGE2D_NIR_NIR_x_.fits' % target

root2 = '../../../X-shooter/P93/Data/Reduction/%s/OB*/NIR/Reduction_Tell/TellCorr_%s_OB*' % (target,target)


Spath = glob.glob(root1)
Tpath = glob.glob(root2)

print target, len(Spath), len(Tpath)

# Loop over all OB's
for OB_no in range(len(Tpath)):

    ### Create paths ###
    path_sci  = Spath[OB_no]
    path_tell = Tpath[OB_no]

    ### Read-in data ###
    Flux_arr, Err_arr, Qual_arr, hd_0, hd_1, hd_2 = read_2d_fits(path_sci)

    ### Read-in tell ###
    data_tell   = np.genfromtxt(path_tell)
    Wave        = data_tell[:,0]
    Tell_obs    = data_tell[:,1]
    Tell_model  = data_tell[:,2]

    ### Transmission function ###
    # Because of bad edge removal the telluric 1d spectrum (not corrected for bad edges
    # due to the 1d format) is 1 pixel longer than the wavelength direction of the flux
    # spectrum
    trans = (Tell_obs/Tell_model)#[1:]
    
    # We give negative ratio a zero weight (trans[i] = 1)
    for i in range(len(trans)):
        if trans[i] <= 0:
            trans[i] = 1


    # Telluric Corrected 2D Flux Image
    Flux_arr_trans, Err_arr_trans = Flux_Telluric_Correction(Flux_arr,Err_arr,trans)

    read_out_tellcorr_arr(path_sci,target,Flux_arr_trans,Err_arr_trans,Qual_arr,hd_0,hd_1,hd_2)




