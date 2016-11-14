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
from scipy import ndimage
from astropy.coordinates import SkyCoord
import astropy.units as u

########################
####### Programs #######
########################

def read_in_COSMOS_cat(path):
    f = pf.open(path)
    #f0 = f[0].data; hd0 = f[0].header
    f1 = f[1].data; hd1 = f[1].header
    return f1


def extract_filters(fdata,Type='uJy'):
    c = 299792458 #m / s
    if Type == 'uJy':
        F_U     = fdata.field('u_FLUX_APER2')
        E_U     = fdata.field('u_FLUXERR_APER2')
        F_B     = fdata.field('B_FLUX_APER2')
        E_B     = fdata.field('B_FLUXERR_APER2')
        F_V     = fdata.field('V_FLUX_APER2')
        E_V     = fdata.field('V_FLUXERR_APER2')
        F_R     = fdata.field('r_FLUX_APER2')
        E_R     = fdata.field('r_FLUXERR_APER2')
        F_I     = fdata.field('ip_FLUX_APER2')
        E_I     = fdata.field('ip_FLUXERR_APER2')
        F_z     = fdata.field('zp_FLUX_APER2')
        E_z     = fdata.field('zp_FLUXERR_APER2')
        F_Y     = fdata.field('yHSC_FLUX_APER2')
        E_Y     = fdata.field('yHSC_FLUXERR_APER2')
        F_J     = fdata.field('J_FLUX_APER2')
        E_J     = fdata.field('J_FLUXERR_APER2')
        F_H     = fdata.field('H_FLUX_APER2')
        E_H     = fdata.field('H_FLUXERR_APER2')
        F_Ks    = fdata.field('Ks_FLUX_APER2')
        E_Ks    = fdata.field('Ks_FLUXERR_APER2')
        # Ftot_K  = fdata.field('Ks_FLUX_APER3')
        # Etot_K  = fdata.field('Ks_FLUXERR_APER3')
        F_Ch1   = fdata.field('SPLASH_1_FLUX')
        E_Ch1   = fdata.field('SPLASH_1_FLUX_ERR')
        F_Ch2   = fdata.field('SPLASH_2_FLUX')
        E_Ch2   = fdata.field('SPLASH_2_FLUX_ERR')
        F_Ch3   = fdata.field('SPLASH_3_FLUX')
        E_Ch3   = fdata.field('SPLASH_3_FLUX_ERR')
        F_Ch4   = fdata.field('SPLASH_4_FLUX')
        E_Ch4   = fdata.field('SPLASH_4_FLUX_ERR')
    
        return F_U,E_U,F_B,E_B,F_V,E_V,F_R,E_R,F_I,E_I,F_z,E_z,F_Y,E_Y,F_J,E_J,F_H,E_H,F_Ks,E_Ks#, F_Ch1,E_Ch1,F_Ch2,E_Ch2,F_Ch3,E_Ch3,F_Ch4,E_Ch4
    
    elif Type == 'cgsA':
        # Change Units from uJy -> erg/s/cm2/Hz -> erg/s/cm2/A
        mJy2freq = 10**(-29) # uJy = 10^(-29) erg/s/cm2/Hz
        Lam_band_av = np.array([3562,4458,5477,6186,7506,8962,10200,12520,16450,21470]) # FILTER.RES.SWv5.R300.info + http://casu.ast.cam.ac.uk/surveys-projects/vista/technical/filter-set
        Freq2Lambda = c*10**(10)/(Lam_band_av**2)
        
        # Selecting the Bands and correcting to erg/s/cm2/A
        F_U = (fdata.field('u_FLUX_APER2')*mJy2freq)*Freq2Lambda[0]
        E_U = (fdata.field('u_FLUXERR_APER2')*mJy2freq)*Freq2Lambda[0]
        F_B = (fdata.field('B_FLUX_APER2')*mJy2freq)*Freq2Lambda[1]
        E_B = (fdata.field('B_FLUXERR_APER2')*mJy2freq)*Freq2Lambda[1]
        F_V = (fdata.field('V_FLUX_APER2')*mJy2freq)*Freq2Lambda[2]
        E_V = (fdata.field('V_FLUXERR_APER2')*mJy2freq)*Freq2Lambda[2]
        F_R = (fdata.field('r_FLUX_APER2')*mJy2freq)*Freq2Lambda[3]
        E_R = (fdata.field('r_FLUXERR_APER2')*mJy2freq)*Freq2Lambda[3]
        F_I = (fdata.field('ip_FLUX_APER2')*mJy2freq)*Freq2Lambda[4]
        E_I = (fdata.field('ip_FLUXERR_APER2')*mJy2freq)*Freq2Lambda[4]
        F_z = (fdata.field('zp_FLUX_APER2')*mJy2freq)*Freq2Lambda[5]
        E_z = (fdata.field('zp_FLUXERR_APER2')*mJy2freq)*Freq2Lambda[5]
        F_Y = (fdata.field('yHSC_FLUX_APER2')*mJy2freq)*Freq2Lambda[6]
        E_Y = (fdata.field('yHSC_FLUXERR_APER2')*mJy2freq)*Freq2Lambda[6]
        F_J = (fdata.field('J_FLUX_APER2')*mJy2freq)*Freq2Lambda[7]
        E_J = (fdata.field('J_FLUXERR_APER2')*mJy2freq)*Freq2Lambda[7]
        F_H = (fdata.field('H_FLUX_APER2')*mJy2freq)*Freq2Lambda[8]
        E_H = (fdata.field('H_FLUXERR_APER2')*mJy2freq)*Freq2Lambda[8]
        F_Ks = (fdata.field('Ks_FLUX_APER2')*mJy2freq)*Freq2Lambda[9]
        E_Ks = (fdata.field('Ks_FLUXERR_APER2')*mJy2freq)*Freq2Lambda[9]
        Ftot_K = (fdata.field('Ks_FLUX_APER3')*mJy2freq)*Freq2Lambda[9]
        Etot_K = (fdata.field('Ks_FLUXERR_APER3')*mJy2freq)*Freq2Lambda[9]
   
        return F_U,E_U,F_B,E_B,F_V,E_V,F_R,E_R,F_I,E_I,F_z,E_z,F_Y,E_Y,F_J,E_J,F_H,E_H,F_Ks,E_Ks,Ftot_K,Etot_K

    elif Type == 'ABmag':
        print 'Update first to include same bands as Type="uJy"'
        # F_U = f1.field('u_MAG_APER2')
        # E_U = f1.field('u_MAGERR_APER2')
        # F_B = f1.field('B_MAG_APER2')
        # E_B = f1.field('B_MAGERR_APER2')
        # F_V = f1.field('V_MAG_APER2')
        # E_V = f1.field('V_MAGERR_APER2')
        # F_R = f1.field('r_MAG_APER2')
        # E_R = f1.field('r_MAGERR_APER2')
        # F_I = f1.field('ip_MAG_APER2')
        # E_I = f1.field('ip_MAGERR_APER2')
        # F_z = f1.field('zp_MAG_APER2')
        # E_z = f1.field('zp_MAGERR_APER2')
        # F_J = f1.field('J_MAG_APER2')
        # E_J = f1.field('J_MAGERR_APER2')
        # F_H = f1.field('H_MAG_APER2')
        # E_H = f1.field('H_MAGERR_APER2')
        # F_Ks = f1.field('Ks_MAG_APER2')
        # E_Ks = f1.field('Ks_MAGERR_APER2')
        # Ftot_K = f1.field('Ks_MAG_APER3')
        # Etot_K = f1.field('Ks_MAGERR_APER3')




def check_ID(min_dist,id_min):
    ID_check = np.where(min_dist < 0.000005)[0] # 0.00001 rad = 2.06 arcsec
    if id_min != ID_check:
        print 'Shiiiit'
    else:
        print 'ID match'


def make_cat_array(targetnames,coordinates,dataphot,t):
    # Extract filters
    if t == 'uJy':
        # F_U,E_U,F_B,E_B,F_V,E_V,F_R,E_R,F_I,E_I,F_z,E_z,F_Y,E_Y,F_J,E_J,F_H,E_H,F_Ks,E_Ks, F_Ch1,E_Ch1,F_Ch2,E_Ch2,F_Ch3,E_Ch3,F_Ch4,E_Ch4 = extract_filters(dataphot,'uJy')
        size = 21
        F_U,E_U,F_B,E_B,F_V,E_V,F_R,E_R,F_I,E_I,F_z,E_z,F_Y,E_Y,F_J,E_J,F_H,E_H,F_Ks,E_Ks = extract_filters(dataphot,t)

    if t == 'cgsA':
        size = 23
        F_U,E_U,F_B,E_B,F_V,E_V,F_R,E_R,F_I,E_I,F_z,E_z,F_Y,E_Y,F_J,E_J,F_H,E_H,F_Ks,E_Ks, Ftot_K, Etot_K = extract_filters(dataphot,t)


    # RA and DEC of Photometric data
    alpha_J2000 = Phot_data.field('ALPHA_J2000')*2*np.pi/360 # radians
    delta_J2000 = Phot_data.field('DELTA_J2000')*2*np.pi/360 # radians

    # Array = np.zeros(shape=(len(targetnames),29))
    Array = np.zeros(shape=(len(targetnames),size))
    for ii in range(len(targetnames)):
        ra_rad_xsh  = coordinates.ra.deg[ii]*2*np.pi/360
        dec_rad_xsh = coordinates.dec.deg[ii]*2*np.pi/360

        # Calculate minimum distance to nearest target
        r = np.sqrt((alpha_J2000-ra_rad_xsh)**2*np.cos(dec_rad_xsh)**2+(delta_J2000-dec_rad_xsh)**2)
        ID = np.where(r == np.min(r))[0]
        check_ID(r,ID)

        # Fill array
        if t == 'uJy':
            # Array[ii] = [ii+1,F_U[ID],E_U[ID],F_B[ID],E_B[ID],F_V[ID],E_V[ID],F_R[ID],E_R[ID],F_I[ID],E_I[ID],F_z[ID],E_z[ID],F_Y[ID],E_Y[ID],F_J[ID],E_J[ID],F_H[ID],E_H[ID],F_Ks[ID],E_Ks[ID],F_Ch1[ID],E_Ch1[ID],F_Ch2[ID],E_Ch2[ID],F_Ch3[ID],E_Ch3[ID],F_Ch4[ID],E_Ch4[ID]]
            Array[ii] = [ii+1,F_U[ID],E_U[ID],F_B[ID],E_B[ID],F_V[ID],E_V[ID],F_R[ID],E_R[ID],F_I[ID],E_I[ID],F_z[ID],E_z[ID],F_Y[ID],E_Y[ID],F_J[ID],E_J[ID],F_H[ID],E_H[ID],F_Ks[ID],E_Ks[ID]]
        

        if t == 'cgsA':
            Array[ii] = [ii+1,F_U[ID],E_U[ID],F_B[ID],E_B[ID],F_V[ID],E_V[ID],F_R[ID],E_R[ID],F_I[ID],E_I[ID],F_z[ID],E_z[ID],F_Y[ID],E_Y[ID],F_J[ID],E_J[ID],F_H[ID],E_H[ID],F_Ks[ID],E_Ks[ID],Ftot_K[ID],Etot_K[ID]]


    return Array


def read_out_catalogue(output_name,Arr,hdr_str):
    fout = open(output_name,'w')
    
    Ntarget = len(Arr[:,0])
    Nfilter = len(Arr[0,:])
    for ii in range(Ntarget):
        Target_filters = Arr[ii,:]
        if ii == 0:
            fout.write('%s\n' % hdr_str)

        for jj in range(Nfilter):
            fout.write('%4.6g ' % Target_filters[jj])
        fout.write('\n')
    fout.close()


################################################################################################
################################################################################################

#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#
#@#@#@ Make Photometric catalogue for FAST #@#@#@



######################
### Set parameters ###


plotnr = 0
type_unit = 'cgsA'


################################
### Input target coordinates ###


Coord_hms_dms = ['10:01:03.03 02:01:04.00','10:00:42.39 02:03:39.39','09:58:13.51 02:10:02.09','09:59:09.48 02:20:28.78','09:59:32.89 02:21:02.50','09:57:43.86 02:30:38.31','10:00:50.02 02:46:19.57','09:59:17.34 02:39:11.23','10:02:58.98 02:00:37.77','10:01:57.00 02:16:12.14','10:01:18.02 01:49:06.74','10:00:50.13 01:51:00.87','10:00:17.74 02:17:52.96','09:59:28.68 02:19:00.13']
coords = SkyCoord(Coord_hms_dms, unit=(u.hourangle, u.deg))

Namelist = ['105842','108899','155853','171060','171687','239220','230929','250513','773654','90676','CP-540713','CP-561356','CP-1243752','CP-1291751']



##########################
### Read in Photometry ###
Phot_data = read_in_COSMOS_cat('../../../X-shooter/cQGsample/COSMOS2015_Photometry/UVISTA-dr2_COSMOS_chi2_v7_03_02_15.fits')

# # Read in COSMOS z_phot (La Phare galaxy template fits)
# path_z = '/Volumes/DataDrive/X-shooter/P93/COSMOS2015_Photometry/zphot_ID_RA_DEC.txt'
# data = np.genfromtxt(path_z)
# phot_z_gal = data[3,:]
print 'Photometry read'

####################################################
### Create array with the filters of each target ###

Catalogue_arr = make_cat_array(Namelist,coords,Phot_data,type_unit)
print 'Catalogue created'

############
### Plot ###

# ....


################
### Read out ###
# header_string = '#     id    U_colf     U_cole     B_colf     B_cole     V_colf     V_cole     R_colf     R_cole     I_colf     I_cole     z_colf     z_cole     Y_colf     Y_cole     J_colf     J_cole     H_colf     H_cole     K_colf     K_cole     Ch1_colf     Ch1_cole     Ch2_colf     Ch2_cole     Ch3_colf     Ch3_cole     Ch4_colf     Ch4_cole'
# read_out_path = 'Output/PhotCat_xsh_all_V1_uJy.cat'

if type_unit == 'cgsA':
    header_string = '#     id    U_colf     U_cole     B_colf     B_cole     V_colf     V_cole     R_colf     R_cole     I_colf     I_cole     z_colf     z_cole     Y_colf     Y_cole     J_colf     J_cole     H_colf     H_cole     K_colf     K_cole K_totf K_tote'


    read_out_path = 'Output/PhotCat_xsh_all_V1_cgsA_NoSPLASH.cat'
    read_out_catalogue(read_out_path,Catalogue_arr,header_string)

if type_unit == 'uJy':
    header_string = '#     id    U_colf     U_cole     B_colf     B_cole     V_colf     V_cole     R_colf     R_cole     I_colf     I_cole     z_colf     z_cole     Y_colf     Y_cole     J_colf     J_cole     H_colf     H_cole     K_colf     K_cole'

    read_out_path = 'Output/PhotCat_xsh_all_V1_uJy_NoSPLASH.cat'
    read_out_catalogue(read_out_path,Catalogue_arr,header_string)







