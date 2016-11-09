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
import math


# import sys
# sys.path.insert(0, '../NIR/')
# from functions.py import *



######################## 
####### Programs ####### 

def read_fits_file(path,norm_const):
    f = pf.open(path)
    hd_0 = f[0].header
    hd_1 = f[1].header

    F = f[0].data*norm_const
    E = f[1].data*norm_const
    W = 10*(hd_0['CRVAL1'] + (hd_0['CRPIX1'] - 1 + np.arange(hd_0['NAXIS1']))*hd_0['CDELT1'])

    return F, E, W, hd_0, hd_1



def Get_binsize(d_nm):
    """ # The binsize needs to match the resolution of the BC03 models in FAST which is 3A in rest-frame. For z~2
        # galaxies this correspond to 3A*(1+2) = 9A. We set it equal to 10 A, to include higher possible z values.
        # The xshooter data have 0.06 nm/A = 0.6 A/pix, and the BC03 models are thus higher quality than the data.
        # From this we can calculate the number of pixels that we will bin; 9 A / (0.6 A/pix) = 15 pix.
    """
    BC03_z2_res = 9
    A2pix = d_nm*10
    bs = int(BC03_z2_res/A2pix)
    print 'Binsize = BC03_res(z~2)/spec_res = %s A/(%s A/pix) = %s pix' % (BC03_z2_res,A2pix,bs)
    return bs


def read_out_catalogue(output_name,Arr,hdr_str):
    fout = open(output_name,'w')
    
    Ntarget = len(Arr[:,0])
    Nfilter = len(Arr[0,:])
    for ii in range(Ntarget):
        Target_filters = Arr[ii,:]
        if ii == 0:
            fout.write('%s' % hdr_str)

        for jj in range(Nfilter):
            fout.write('%4.6g ' % Target_filters[jj])
        fout.write('\n')
    fout.close()


########################
########################



############################################
###### Create spec catalogue for FAST ######
############################################

Spath = glob.glob('../../../X-shooter/cQGsample/Objects/*/VIS_NIR_Banana/*_rebin0.9nm_opt_V1.fits')

# Spath = [Spath[4],Spath[6]]

Normalisation_const = 1e19


# Read in the Flux, Errs, and Wavelength
Namelist = []
for ii in range(len(Spath)):
    Image = Spath[ii]
    Namelist.append(Image.split('/')[5])
    
    ## Read in data
    exec('Flux_%s, Errs_%s, Wave, hd_0, hd_1 = read_fits_file(Image,Normalisation_const)' % (ii,ii))

    # Here I remove the K-band from P86 target due to its bad quality
    exec('Flux_%s = Flux_%s[:1742]' % (ii,ii))
    exec('Errs_%s = Errs_%s[:1742]' % (ii,ii))
    exec('Wave = Wave[:1742]')

    # exec('print Flux_%s.shape, Errs_%s.shape, Wave.shape' % (ii,ii))


# Create binsize
nm2pix = hd_0['CDELT1']
binsize = Get_binsize(nm2pix)



### Selecting only the good pixel regions ###




Range_array = [5850,13500,14240,17919,19641,20300]
mask_goodpix = (Wave > Range_array[0]) & (Wave < Range_array[1]) | (Range_array[2] < Wave) & (Wave < Range_array[3]) | (Range_array[4] < Wave) & (Wave < Range_array[5])
goodPixels = np.where(mask_goodpix == True)[0]

Wave_orig = np.copy(Wave)
Wave = Wave[goodPixels]
for i in range(len(Spath)):
    exec('plt.plot(Wave_orig, Flux_%s,"b")' % i)


    exec('Flux_%s = Flux_%s[goodPixels]' % (i,i))
    exec('Errs_%s = Errs_%s[goodPixels]' % (i,i))

    # Change negative flux values to zero
    #exec('Flux_%s[Flux_%s < 0] = 0' % (i,i))

    # 5sigma clip the outliers
    exec('medFlux = np.median(Flux_%s)' % (i))


    # exec('upper_lim = Flux_%s > (medFlux+10*Errs_%s)' % (i,i))
    # exec('Id_5sigma = np.where(upper_lim)')
    # exec('Flux_%s[Id_5sigma] = medFlux' % (i))

    # exec('lower_lim = Flux_%s < (medFlux-15*Errs_%s)' % (i,i))
    # exec('Id_5sigma = np.where(lower_lim)')
    # exec('Flux_%s[Id_5sigma] = medFlux' % (i))

    exec('plt.plot(Wave, Flux_%s,"r")' % i)
    plt.ylim([-10,80])
    # plt.show()

# plt.plot(Wave, Flux_0,'r')
# plt.show()


#############################################
#############################################

# We have also tried to use mean, weighted mean, but the median represents the spectra better visually.


Arr_med = np.zeros(shape=(len(Wave),2*len(Spath)+3))
Arr_med[:,1] = Wave
Arr_med[:,2] = Wave*0+1

Arr_wm = np.zeros(shape=(len(Wave),2*len(Spath)+3))
Arr_wm[:,1] = Wave
Arr_wm[:,2] = Wave*0+1

n_itt = int(math.floor(len(Wave)/binsize))

eps = 1e-40
for j in range(n_itt):
    Arr_med[j*binsize:(j+1)*binsize,0] = np.zeros(binsize)+j
    Arr_wm[j*binsize:(j+1)*binsize,0] = np.zeros(binsize)+j
    for u in range(len(Spath)):
        # Median binning
        exec('Arr_med[j*binsize:(j+1)*binsize,'+str(u*2+3)+'] = np.zeros(binsize) + np.median(Flux_'+str(u)+'[j*binsize:(j+1)*binsize])')
        exec('Arr_med[j*binsize:(j+1)*binsize,'+str(u*2+4)+'] = Errs_'+str(u)+'[j*binsize:(j+1)*binsize]')

        # Weigthed Mean binning
        exec('w_err = 1/(Errs_%s[j*binsize:(j+1)*binsize]**2+eps)' %  (u))

        exec('Arr_wm[j*binsize:(j+1)*binsize,%s] = np.zeros(binsize) + np.sum(Flux_%s[j*binsize:(j+1)*binsize]*w_err)/np.sum(w_err)' % (str(u*2+3),str(u)))
        exec('Arr_wm[j*binsize:(j+1)*binsize,%s] = 1/np.sqrt(np.sum(w_err))' % (u*2+4))



# Remove the additional zeros at the end of the file, arising from the binsize not matching the len(wave)
H = len(Wave)%binsize
for i in range(H):
    Arr_med = np.delete(Arr_med,-1,0)
    Arr_wm = np.delete(Arr_wm,-1,0)

plt.subplot(2,1,1)
plt.plot(Wave,Flux_0-Arr_wm[:,3],'gray',alpha=0.3)
plt.plot(Arr_wm[:,1],Arr_wm[:,3],'g')
plt.ylim([-10,20])

plt.subplot(2,1,2)
plt.plot(Wave,Flux_0,'gray',alpha=0.3)
plt.plot(Arr_med[:,1],Arr_med[:,3],'g')
plt.ylim([-10,20])
plt.title('Median binning')

# plt.show()




# Read out the file
Arr = Arr_wm

header_string = '#  bin  wl  trans  F1  E1  F2  E2  F3  E3  F4  E4  F5  E5  F6  E6  F7  E7  F8  E8  F9  E9  F10  E10  F11  E11  F12  E12  F13  E13  F14  E14\n'

read_out_catalogue('Output/cQG_all_bin%snm_BC03_V1.spec' % (nm2pix),Arr_wm,header_string)







##########################################################################
### Plotting of the median and weighted average binning of the spectra ###
'''

i=0
#plt.plot(Wave,Flux_0,'green')



path1 = '/Users/mstockmann/X-Shooter/P93/*/Combined_OBs/*V2_medrebin20_opt.fits'
spath1 = glob.glob(path1)

f = pf.open(spath1[i])
hd_0 = f[0].header
nm2pix = hd_0['CDELT1']
#print hd_0
#print f.info()

Wave = np.array((hd_0['CRVAL1'] + (hd_0['CRPIX1'] - 1 + np.arange(hd_0['NAXIS1']))*hd_0['CDELT1'])*10) #Aangstroem

exec('Flux_%s = f[0].data/1e-19' % (i))
exec('Errs_%s = f[1].data/1e-19' % (i))

plt.plot(Wave,Flux_0,'green')



Wave = Arr_med[:,1]
plt.plot(Wave,Arr_med[:,3],'b')

Wave_wm = Arr_wm[:,1]
plt.plot(Wave_wm,Arr_wm[:,3],'r')


plt.plot(Wave_wm,Arr_wm[:,4],'purple')
#plt.plot(Wave,Errs_1,'green')

plt.ylim([-5,20])
plt.show()

'''

##########################################################################
##########################################################################







