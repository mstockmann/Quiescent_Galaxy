

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


def weighted_mean(F,E):
    w = 1/(E**2)
    F_wm = np.sum(w*F)/np.sum(w)
    E_wm = 1/np.sqrt(np.sum(w))
    return F_wm,E_wm

def SwapZerosAndOnes(Map):
    Map_swap = np.zeros(shape=(len(Map[:,0]),len(Map[0,:])))
    for j in range(len(Map[:,0])):
        for i in range(len(Map[0,:])):
            if Map[j,i] == 0:
                Map_swap[j,i] = 1
            else:
                Map_swap[j,i] = 0
    return Map_swap








#######################################
##### Rebin a single 2D spectrum  #####
#######################################

# P93

# target = 90676
# name = '90676_VIS_avwCombined_OBx5_sig5_V1.fits'

# target = 105842
# name = '105842_VIS_avwCombined_OBx4_sig5_V1.fits'

# target = 108899
# name = '108899_VIS_avwCombined_OBx5_sig5_V1.fits'

# target = 155853
# name = '155853_VIS_avwCombined_OBx5_sig5_V1.fits'

# target = 171060
# name = '171060_VIS_avwCombined_OBx5_sig5_V1.fits'

# target = 171687
# name = '171687_VIS_avwCombined_OBx5_sig5_V1.fits'

target = 230929
name = '230929_NIRBanana_avwCombined_tellcorr_6xOB_sig5_V1.fits'

# target = 239220
# name = '239220_VIS_avwCombined_OBx5_sig5_V1.fits'

# target = 250513
# name = '250513_VIS_avwCombined_OBx5_sig5_V1.fits'

# target = 773654
# name = 'xsh_scired_slit_nod_SCI_SLIT_FLUX_MERGE2D_VIS_VIS_1x2_100k.fits'


print target


# UDS targets
#target = 250513
#name = 'sumnir.fits'
#Image = '/Users/mstockmann/X-Shooter/P83_?/'+name
#f = pf.open(Image)
#F = f[0].data; hd_0 = f[0].header

# STOKKYLOKKY
#path = '/Volumes/STOKKYLOKKY/X-shooter/P87/MACS0454/Combined_OBs/'
#name = 'MACS0454_avwCombined_OBx4_tellcorr_sig5_V1.fits'



# Read in Image

path = '../../Data/Reduction/%s/Combined_OBs/NIR_Banana/' % (target)

Image = path+name

# print Image
# raise

f = pf.open(Image)

Flux = f[0].data; hd_0 = f[0].header
Err  = f[1].data; hd_1 = f[1].header
BPM  = f[2].data; hd_2 = f[2].header


#NCom = f[3].data; hd_3 = f[3].header
# S2N  = f[4].data; hd_4 = f[4].header

Wave = hd_0['CRVAL1'] + (hd_0['CRPIX1'] - 1 + np.arange(hd_0['NAXIS1']))*hd_0['CDELT1'] # Aangstroem




# Swap 0 <-> 1 to use 
BPM_alt = SwapZerosAndOnes(BPM)


# Rebin 2D spectrum
F = Flux
E = Err
M = BPM_alt


nm2pix = hd_0['CDELT1']
New_binsize_Ang = 9
N_pix = New_binsize_Ang/(nm2pix*10)

#N_pix = 60.0;# print 'Binning = %s pix (%s nm)' % (N_pix,N_pix*nm2pix)
N_chunk = int(np.round(len(F[0,:])/N_pix,0)); print 'N_chunk = %s' % N_chunk


Flux_bin = np.zeros(shape=(len(F[:,0]),N_chunk))
Errs_bin = np.zeros(shape=(len(F[:,0]),N_chunk))

#eps = 1e-40 # Softening in inverse error w_err
for i in range(int(N_chunk)):
    # The chunks of Flux and Err, multiplied by the M (bad=0, good=1)
    F_chunk = (F[:,i*N_pix:(i+1)*N_pix])
    E_chunk = (E[:,i*N_pix:(i+1)*N_pix])

    for j in range(len(F[:,0])):
        # Weighted mean
        Flux_bin[j,i], Errs_bin[j,i] = weighted_mean(F_chunk[j,:],E_chunk[j,:])




### Read out the fucking binned spectrum ###
path_out = path
file_out = name.split('.')[0]+'_wmrebin%s.fits' % (int(N_pix))

if not os.path.exists(path_out+file_out):
    # Changing the wavelength solution to a rebinnined version
    hd_0['CRVAL1'] = hd_0['CRVAL1']+(N_pix-1)*hd_0['CDELT1']/2
    hd_0['CDELT1'] = hd_0['CDELT1']*N_pix # CDELT1 and CD1_1 are the same
    hd_0['CD1_1'] = hd_0['CD1_1']*N_pix
    hd_0['NAXIS1'] = N_chunk
    #print hd_0['CRVAL1'],hd_0['CRPIX1'], hd_0['NAXIS1'], hd_0['CDELT1']
    pf.writeto(path_out+file_out, Flux_bin, hd_0)
    
    # Read out error array
    pf.append(path_out+file_out, Errs_bin, hd_1)

    # Read out bad pixel map
    pf.append(path_out+file_out, BPM_alt, hd_2)

    # Number of Combinations
    #pf.append(path_out+file_out, Ncombin, hd_3)

    # S/N spectrum
    # S2N_new = Flux_bin/Errs_bin
    # pf.append(path_out+file_out, S2N_new, hd_4)


else:
    print 'file already exists'


###############################
###############################







##################################################
### Rebin a series of spectra P93 (04.05.2016) ###
'''

target_list = ['CP-540713','CP-561356','CP-1243752','CP-1291751']

New_binsize_Ang = 9


for ii in range(len(target_list)):
    target = target_list[ii]

    path_name = glob.glob('/Volumes/DataDrive/X-shooter/P86_COSMOS_1/%s/Combined_OBs/VIS/%s_VIS_avwCombined_*xOB_sig5_V1.fits' % (target,target))[0]

    # Read in Image
    f = pf.open(path_name)

    Flux = f[0].data; hd_0 = f[0].header
    Err  = f[1].data; hd_1 = f[1].header
    BPM  = f[2].data; hd_2 = f[2].header

    # S2N  = f[4].data; hd_4 = f[4].header
    Wave = hd_0['CRVAL1'] + (hd_0['CRPIX1'] - 1 + np.arange(hd_0['NAXIS1']))*hd_0['CDELT1'] # Aangstroem


    # Swap 0 <-> 1 to use 
    BPM_alt = SwapZerosAndOnes(BPM)


    # Rebin 2D spectrum
    F = Flux
    E = Err
    M = BPM_alt


    nm2pix = hd_0['CDELT1']
    N_pix = New_binsize_Ang/(nm2pix*10)


    #N_pix = 60.0;# print 'Binning = %s pix (%s nm)' % (N_pix,N_pix*nm2pix)
    N_chunk = int(np.round(len(F[0,:])/N_pix,0)); print 'N_chunk = %s' % N_chunk


    Flux_bin = np.zeros(shape=(len(F[:,0]),N_chunk))
    Errs_bin = np.zeros(shape=(len(F[:,0]),N_chunk))

    #eps = 1e-40 # Softening in inverse error w_err
    for i in range(int(N_chunk)):
        # The chunks of Flux and Err, multiplied by the M (bad=0, good=1)
        F_chunk = (F[:,i*N_pix:(i+1)*N_pix])
        E_chunk = (E[:,i*N_pix:(i+1)*N_pix])

        for j in range(len(F[:,0])):
            # Weighted mean
            Flux_bin[j,i], Errs_bin[j,i] = weighted_mean(F_chunk[j,:],E_chunk[j,:])




    ### Read out the fucking binned spectrum ###
    path_out = '/'.join(path_name.split('/')[:-1])+'/'
    file_out = (path_name.split('/')[-1]).split('.')[0]+'_wmrebin%s.fits' % (N_pix)

    print path_out+file_out

    if not os.path.exists(path_out+file_out):
        # Changing the wavelength solution to a rebinnined version
        hd_0['CRVAL1'] = hd_0['CRVAL1']+(N_pix-1)*hd_0['CDELT1']/2
        hd_0['CDELT1'] = hd_0['CDELT1']*N_pix # CDELT1 and CD1_1 are the same
        hd_0['CD1_1'] = hd_0['CD1_1']*N_pix
        hd_0['NAXIS1'] = N_chunk
        #print hd_0['CRVAL1'],hd_0['CRPIX1'], hd_0['NAXIS1'], hd_0['CDELT1']
        pf.writeto(path_out+file_out, Flux_bin, hd_0)
        
        # Read out error array
        pf.append(path_out+file_out, Errs_bin, hd_1)

        # Read out bad pixel map
        pf.append(path_out+file_out, BPM_alt, hd_2)

        # Number of Combinations
        #pf.append(path_out+file_out, Ncombin, hd_3)

        # S/N spectrum
        S2N_new = Flux_bin/Errs_bin
        pf.append(path_out+file_out, S2N_new, hd_0)


    else:
        print 'file already exists'


'''






##################################################
##################################################













































##### Old code chunks #####



###############################
##### Rebin of 1D spectra #####
###############################
'''
target = 105842
path = '/Users/mstockmann/X-Shooter/P93/%s/Combined_OBs/' % target

#name = '105842_avwCombined_OBx4_sig5_V1_opt.fits'
#name = '105842_avwCombined_OBx4_sig5_V1_rebin20_opt.fits'



# Read in image
Image = path+name
f = pf.open(Image)

Flux = f[0].data; hd_0 = f[0].header
Err  = f[1].data; hd_1 = f[1].header
Wave = hd_0['CRVAL1'] + (hd_0['CRPIX1'] - 1 + np.arange(hd_0['NAXIS1']))*hd_0['CDELT1'] # Angstrom


print f.info()
#plt.plot(Wave,Flux,'black')
plt.plot(Wave,Err,'r')
plt.show()
'''


'''

#### 1 d binning ####
nm2pix = hd_0['CDELT1']
N_pix = 16; print 'Binning = %s pix (%s nm)' % (N_pix,N_pix*nm2pix)

N_chunk = int(np.round(len(Flux)/N_pix,0)); print 'N_chunk = %s' % N_chunk



Flux_bin = np.zeros(N_chunk)
Errs_bin = np.zeros(N_chunk)

eps = 0#1e-40 # Softening in inverse error w_err
for i in range(int(N_chunk)):
    F_chunk = Flux[i*N_pix:(i+1)*N_pix]
    E_chunk = Err[i*N_pix:(i+1)*N_pix]

    w_err = 1/(E_chunk**2+eps)
    #w_err = w_err*10**40

    Flux_bin[i] = np.sum(w_err*F_chunk)/np.sum(w_err)
    Errs_bin[i] = 1/np.sqrt(np.sum(w_err))


plt.plot(Flux_bin)
plt.axis([0,1000,-2e-18,4e-18])
plt.show()
'''


###############################
###############################





