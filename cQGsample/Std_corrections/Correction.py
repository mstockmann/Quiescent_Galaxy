#!/usr/local/bin/python

#-*- coding: utf-8 -*-

from __future__ import division
import glob
import pyfits as pf
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy import *
import os
from scipy.optimize import curve_fit
import sys
from cardelli_unred import cardelli_reddening

# sys.path.insert(0, '../..')
# from Stokky import *

sys.path.insert(0, '../..')
import Stokky as st

####################################################



def correct_for_dust(wavelength, ra, dec):
	"""Query IRSA dust map for E(B-V) value and returns reddening array
	----------
	wavelength : numpy array-like
	    Wavelength values for which to return reddening
	ra : float
	    Right Ascencion in degrees
	dec : float
	    Declination in degrees
	Returns
	-------
	reddening : numpy array
	Notes
	-----
	For info on the dust maps, see http://irsa.ipac.caltech.edu/applications/DUST/
	"""

	from astroquery.irsa_dust import IrsaDust
	import astropy.coordinates as coord
	import astropy.units as u
	C = coord.SkyCoord(ra*u.deg, dec*u.deg, frame='fk5')
	dust_image = IrsaDust.get_images(C, radius=2 *u.deg, image_type='ebv', timeout=60)[0]
	ebv = np.mean(dust_image[0].data[40:42, 40:42])
	r_v = 3.1
	av =  r_v * ebv
	from specutils.extinction import reddening
	# from extinction import reddening
	return reddening(wavelength* u.angstrom, av, r_v=r_v, model='ccm89'), ebv


def read_phot_cat(path):
	data = np.genfromtxt(path)

	id = data[:,0]
	Uf = data[:,1]
	Ue = data[:,2]
	Bf = data[:,3]
	Be = data[:,4]
	Vf = data[:,5]
	Ve = data[:,6]
	Rf = data[:,7]
	Re = data[:,8]
	If = data[:,9]
	Ie = data[:,10]
	zf = data[:,11]
	ze = data[:,12]
	Yf = data[:,13]
	Ye = data[:,14]
	Jf = data[:,15]
	Je = data[:,16]
	Hf = data[:,17]
	He = data[:,18]
	Kf = data[:,19]
	Ke = data[:,20]
	Ktf = data[:,21]
	Kte = data[:,22]
	Phot_list = [Uf, Ue, Bf, Be, Vf, Ve, Rf, Re, If, Ie, zf, ze, Yf, Ye, Jf, Je, Hf, He, Kf, Ke, Ktf, Kte]
	return Phot_list


def read_filter_curves(name,filter_names):
	""" Output: [no_filter],[row],[column], column: [0]: No.row, [1]:lambda/A, [2]: Transmission curve
	"""


	f = open(name,'r')
	files = f.read().split('\n')
	QSO_files = []
	for i in range(len(files)-1):# -1 as there is a \n in the last entry
	    k = files[i].split('\t')
	    k_list = k[0].split(' ',8)
	    k_list = filter(lambda x: len(x)>0, k_list)
	    tmp = []
	    if k_list[1] in filter_names:
	        filter_len = int(k_list[0])
	        for j in range(filter_len):
	            k_tmp = files[i+j+1].split('\t')
	            k_list_tmp = k_tmp[0].split(' ',8)
	            k_list_tmp = filter(lambda x: len(x)>0, k_list_tmp)
	            k_list_tmp = map(float,k_list_tmp)
	            tmp.append(k_list_tmp)
	        QSO_files.append(tmp)
	        #print len(QSO_files)
	return np.array(QSO_files)

def pair_transmission_with_spec(Wspec,Wtrans,Trans):
	ID_wave = np.searchsorted(Wtrans,Wspec)
	WH_syn = Wtrans[ID_wave]
	TH_syn = Trans[ID_wave]
	return WH_syn, TH_syn


####################################################


###############################
### de-reddening of spectra ###
'''
path_NIR_VIS = glob.glob('../../../X-shooter/cQGsample/Objects/*/VIS_NIR_Banana/*_rebin0.9nm_opt_V1.fits')




for i in range(len(path_NIR_VIS)):
	W, F, E, hdf, hde = read_in_1d_fits(path_NIR_VIS[i])

	print hdf['ra'], hdf['dec']

	Reddening, ebv = correct_for_dust(W, np.float(hdf['ra']), np.float(hdf['dec']))
	Fcorr = cardelli_reddening(W, F, ebv)
	raise
'''
###############################



################################
### Make slitloss correction ###

# NIR+VIS
# Spec_path = glob.glob('../../../X-shooter/cQGsample/Objects/*/VIS_NIR_Banana/*_rebin0.9nm_opt_V1.fits')

# NIR 2d emission corrected
Spec_path = glob.glob('../../../X-shooter/cQGsample/Spectra_analysis/2-OptimalExt_er/*_V1_NIRcorr_wmrebin15_er2d_opt.fits')



Phot_path = glob.glob('../../../Quiescent_Galaxy/FAST/Photometric_Cat/Output/PhotCat_xsh_all_V1_cgsA_NoSPLASH.cat')

Trans_path = glob.glob('../../../../Programs/FAST/COSMOS_filters/FILTER.RES.SWv5.R300.txt')



# read in photometry
Photometric_list = read_phot_cat(Phot_path[0])
Hcolf = Photometric_list[16]
Kcolf = Photometric_list[18]
Ktotf = Photometric_list[20]


# Read in the filter curves
data_filter = read_filter_curves(Trans_path[0],'UltraVISTA/H_ultravista_qe_atm.bp')
W_Hspec = np.array(data_filter[0,:,1])
T_Hspec = np.array(data_filter[0,:,2])


for ii in range(len(Spec_path)):



	# read in spec - change when bpm is included in 1d nir+vis
	# data_arr = pf.open(Spec_path[0])
	# hdf = data_arr[0].header
	# hde = data_arr[1].header
	# F = data_arr[0].data#/np.median(data_arr[0].data)
	# E = data_arr[1].data#/np.median(data_arr[0].data)
	# W = (hdf['CRVAL1'] + (hdf['CRPIX1'] - 1 + np.arange(hdf['NAXIS1']))*hdf['CDELT1'])*10 

	W, F, E, M, hdf, hde, hdm = st.read_in_1d_fits(Spec_path[0])

	# pair transmission with spec arr
	WH_pair, TH_pair = pair_transmission_with_spec(W,W_Hspec,T_Hspec)


	# plt.plot(W,F,color='black')
	# plt.plot(WH_pair,TH_pair,color='r')
	# plt.show()

	# make synthetic mag
	H_synthetic = np.sum(F*TH_pair)/np.sum(TH_pair)


	# make slit loss and aperture corrections to the flux spec
	corr_aperture = Ktotf[ii]/Kcolf[ii]

	corr_slitloss = Hcolf[ii]*corr_aperture/H_synthetic
	print corr_slitloss, corr_aperture

	F_corr = F*corr_slitloss#*corr_aperture

	
	Lam_band_av = np.array([3562,4458,5477,6186,7506,8962,10200,12520,16450,21470])

	Phot_F_target = np.zeros(10)
	Phot_E_target = np.zeros(10)
	for uu in range(10):
		Phot_F_target[uu] = Photometric_list[uu*2][ii]
		Phot_E_target[uu] = Photometric_list[uu*2+1][ii]


	# print len(Photometric_list)
	# raise

	# print Photometric_list[::2]
	# print Lam_band_av
	# raise

	
	plt.plot(W,F,color='black')
	plt.plot(W,F_corr,color='red')
	plt.errorbar(Lam_band_av,Phot_F_target*corr_aperture,yerr=Phot_E_target,color='green',marker='s',)
	plt.axis([2500,23000,-0.2e-18,4e-18])
	plt.show()

	# read out
	# NIR+VIS
	# name_out = Spec_path[ii].split('/')[-1].replace('_avwCombined_OBx5_sig5','').replace('opt_','opt_stdcorr_')
	# path = '../../../X-shooter/cQGsample/Spectra_analysis/3-Stdcorr_slit_aper/NIR+VIS/%s' % name_out

	# NIR
	name_out = Spec_path[ii].split('/')[-1].replace('opt','opt_stdcorr')
	path = '../../../X-shooter/cQGsample/Spectra_analysis/3-Stdcorr_slit_aper/NIR/%s' % name_out
	# print path_out
	

	st.read_out_3arr_2dfits(path,F_corr,E,M,hdf,hde,hdm)




	# Change to st.read_out_3arr_2dfits(...) when bpm is included
	# if not os.path.exists(path):
	# 	# Read out flux array
	# 	pf.writeto(path, F_corr, hdf)

	# 	# Read out error array
	# 	pf.append(path, E, hde)

	# else:
	# 	os.system('mv %s' % path.replace('.fits','_old.fits'))
	# 	# Read out flux array
	# 	pf.writeto(path, F_corr, hdf)

	# 	# Read out error array
	# 	pf.append(path, E, hde)
	# 	print 'former file moved to *_old.fits'
	# raise














########################################
### heliocentric velocity correction ###

""" The P93 program was taking over the whole year, and the position of the earth was worst in the case of 108899 and 773654 where the difference in radial velocity (dfits keyword: QC.VRAD.BARYCOR) were maximum 10 km/s.

	If estimating the maximal error associated with this we get for lambda = 20.000 a
	correction of 1+vrad/c = 0.8 Angstrom. As most of our thin absorption lines are at 
	12 this effect can be neglected.

	If one were to correct for the heliocentric reference frame we would have to 
	multiply lambda_corr = (1+v/c)*lambda_uncorr, for each OB that have a different 
	radial velocity.
"""

########################################
########################################


















# for i in range(len(Nspec)):
# 	# i = 10
# 	target = Nspec[i].split('/')[-1]
# 	print target

# 	### Rebinned ###
# 	# CDELT1 = 0.9
# 	#path_NIR = '/Volumes/DataDrive/X-shooter/P93/Data/%s/Combined_OBs/NIR/*V3_wmrebin15_opt.fits' % (target)
# 	#path_VIS = '/Volumes/DataDrive/X-shooter/P93/Data/%s/Combined_OBs/VIS/*V1_wmrebin45_opt.fits' % (target)



# 	# print path_NIR
# 	# print path_VIS
# 	# print len(path_NIR), len(path_VIS)
# 	# raise




# 	Wave_nir, Flux_nir, Errs_nir, hdf_nir, hde_nir = read_in_1d_fits(path_NIR[0])	
# 	Wave_vis, Flux_vis, Errs_vis, hdf_vis, hde_vis = read_in_1d_fits(path_VIS[0])

# 	print 'Object: %s, NIR_CD1 / VIS_CD1 = %s / %s' % (target,hdf_nir['CDELT1'], hdf_vis['CDELT1'])

# 	print len(Wave_nir), len(Flux_nir)


# 	# plt.plot(Wave_nir,Flux_nir)
# 	# plt.show()

# 	# Create a new collective wavelength grid
# 	# Note: We tried the easy way using np.union1d() but it left out some of the overlaps (??)
# 	Wave_vis_match_point = Wave_vis[np.searchsorted(Wave_vis,Wave_nir[0])]
# 	Match_diff = Wave_vis_match_point-Wave_nir[0]
# 	Wave_vis_match = np.array(Wave_vis)-Match_diff

# 	# print type(Wave_vis_match), type(Wave_nir)
# 	# print Wave_vis_match[-1]
# 	# print Wave_nir[430:440]
# 	# print np.searchsorted(Wave_nir,Wave_vis_match[-1])

# 	ID_vis_left = np.where(Wave_vis_match == Wave_nir[0])[0]
# 	ID_nir_right = np.searchsorted(Wave_nir,Wave_vis_match[-1])


# 	Wave_vis_overlap = Wave_vis_match[ID_vis_left:]
# 	Wave_nir_overlap = Wave_nir[:ID_nir_right+1]

# 	Wave_vis_nir = np.concatenate((Wave_vis_match[:ID_vis_left],Wave_nir_overlap,Wave_nir[ID_nir_right+1:]),axis=0)


# 	# Create a new collective Flux grid
# 	Flux_overlap_vis = Flux_vis[ID_vis_left:]
# 	Flux_overlap_nir = Flux_nir[:ID_nir_right+1]

# 	Errs_overlap_vis = Errs_vis[ID_vis_left:]
# 	Errs_overlap_nir = Errs_nir[:ID_nir_right+1]


# 	# Combine overlapping fluxes (When plotting them they look similar to highest points, but they are different when zooming)
# 	Flux_avw = np.zeros(len(Flux_overlap_nir))
# 	Errs_avw = np.zeros(len(Flux_overlap_nir))

# 	for i in range(len(Flux_overlap_vis)):
# 		err_tmp = np.array([Errs_overlap_vis[i],Errs_overlap_nir[i]])
# 		flux_tmp = np.array([Flux_overlap_vis[i],Flux_overlap_nir[i]])

# 		w_err = 1/(err_tmp**2)
# 		Flux_avw[i] = np.sum(w_err*flux_tmp)/np.sum(w_err)
# 		Errs_avw[i] = 1/np.sqrt(np.sum(w_err)) 


# 	Flux_vis_nir = np.concatenate((Flux_vis[:ID_vis_left],Flux_avw,Flux_nir[ID_nir_right+1:]),axis=0)
# 	Errs_vis_nir = np.concatenate((Errs_vis[:ID_vis_left],Errs_avw,Errs_nir[ID_nir_right+1:]),axis=0)



# 	# We have change the length (NAXIS1) and the start point (CRVAL), but we keep the binning
# 	hdf_vis['CRVAL1'] = Wave_vis_nir[0]/10 # in nm (not Angstrom)
# 	hdf_vis['NAXIS1'] = len(Flux_vis_nir)

# 	hde_vis['CRVAL1'] = Wave_vis_nir[0]/10 # in nm (not Angstrom)
# 	hde_vis['NAXIS1'] = len(Flux_vis_nir)
	

# 	# Read out the fucking file_out
# 	path_out = '/'.join(path_NIR[0].split('/')[:-2])+'/VIS_NIR_Banana/'
# 	file_out = '%s_VIS_NIR_avwCombined_OBx5_sig5_rebin%snm_opt_V1.fits' % (target,hde_vis['CDELT1'])
# 	# print path_out+file_out
# 	# raise

# 	if not os.path.exists(path_out+file_out):
# 	    # Read out flux array
# 	    pf.writeto(path_out+file_out, Flux_vis_nir, hdf_vis)
	    
# 	    # Read out error array
# 	    pf.append(path_out+file_out, Errs_vis_nir, hde_vis)

# 	else:
# 	    print 'file already exists'



# 	## Plot ###
# 	plt.scatter(Wave_vis-Match_diff,Flux_vis)
# 	plt.scatter(Wave_nir,Flux_nir,color='r')
# 	plt.scatter(Wave_vis_overlap,Flux_overlap_vis,marker='s',s=70,color='magenta')
# 	plt.scatter(Wave_nir_overlap,Flux_overlap_nir,marker='x',s=70,color='black')
# 	plt.scatter(Wave_nir_overlap,Flux_avw,marker='o',s=70,color='green')
# 	plt.plot(Wave_vis_nir,Flux_vis_nir+1.9e-17,color='black',label='Flux+offset')

# 	plt.scatter(Wave_vis-Match_diff,Flux_vis)


# 	plt.ylim([-0.2e-17,0.8e-17])
# 	plt.xlim([9900,10200])
# 	plt.title('%s' % target)
# 	plt.show()





















# ###########################################################
# ### Combine NIR and VIS for 105842 left and right trace ###

# '''
# Nspec = glob.glob('/Volumes/DataDrive/X-shooter/P93/Data/*/Combined_OBs/VIS/*_V1_wmrebin45.fits')


# for i in range(1):#len(Nspec)):
# 	target = 105842 #Nspec[i].split('/')[6]
# 	print target


# 	### Rebinned ###
# 	# CDELT1 = 1.2
# 	path_NIR = '/Volumes/DataDrive/X-shooter/P93/Data/%s/Combined_OBs/NIR/*V3_wmrebin15_opt.fits' % (target)
# 	Wave_nir, Flux_nir, Errs_nir, hdf_nir, hde_nir = read_in_1d_fits(glob.glob(path_NIR)[0])

# 	# CDELT1 = 1.2
# 	path_VIS = '/Volumes/DataDrive/X-shooter/P93/Data/%s/Combined_OBs/VIS/*V1_wmrebin45_opt_right.fits' % (target)
# 	Wave_vis, Flux_vis, Errs_vis, hdf_vis, hde_vis = read_in_1d_fits(glob.glob(path_VIS)[0])


# 	# Create a new collective wavelength grid
# 	# Note: We tried the easy way using np.union1d() but it left out some of the overlaps (??)
# 	Wave_vis_match_point = Wave_vis[np.searchsorted(Wave_vis,Wave_nir[0])]
# 	Match_diff = Wave_vis_match_point-Wave_nir[0]
# 	Wave_vis_match = Wave_vis-Match_diff


# 	ID_vis_left = np.where(Wave_vis_match == Wave_nir[0])[0]
# 	ID_nir_right = np.where(Wave_nir == Wave_vis_match[-1])[0]

# 	Wave_vis_overlap = Wave_vis_match[ID_vis_left:]
# 	Wave_nir_overlap = Wave_nir[:ID_nir_right+1]

# 	Wave_vis_nir = np.concatenate((Wave_vis_match[:ID_vis_left],Wave_nir_overlap,Wave_nir[ID_nir_right+1:]),axis=0)


# 	# Create a new collective Flux grid
# 	Flux_overlap_vis = Flux_vis[ID_vis_left:]
# 	Flux_overlap_nir = Flux_nir[:ID_nir_right+1]

# 	Errs_overlap_vis = Errs_vis[ID_vis_left:]
# 	Errs_overlap_nir = Errs_nir[:ID_nir_right+1]
	 

# 	# Combine overlapping fluxes (When plotting them they look similar to highest points, but they are different when zooming)
# 	Flux_avw = np.zeros(len(Flux_overlap_vis))
# 	Errs_avw = np.zeros(len(Flux_overlap_vis))

# 	for i in range(len(Flux_overlap_vis)):
# 		err_tmp = np.array([Errs_overlap_vis[i],Errs_overlap_nir[i]])
# 		flux_tmp = np.array([Flux_overlap_vis[i],Flux_overlap_nir[i]])

# 		w_err = 1/(err_tmp**2)
# 		Flux_avw[i] = np.sum(w_err*flux_tmp)/np.sum(w_err)
# 		Errs_avw[i] = 1/np.sqrt(np.sum(w_err)) 


# 	Flux_vis_nir = np.concatenate((Flux_vis[:ID_vis_left],Flux_avw,Flux_nir[ID_nir_right+1:]),axis=0)
# 	Errs_vis_nir = np.concatenate((Errs_vis[:ID_vis_left],Errs_avw,Errs_nir[ID_nir_right+1:]),axis=0)



# 	# We have change the length (NAXIS1) and the start point (CRVAL), but we keep the binning
# 	hdf_vis['CRVAL1'] = Wave_vis_nir[0]/10 # in nm (not Angstrom)
# 	hdf_vis['NAXIS1'] = len(Flux_vis_nir)

# 	hde_vis['CRVAL1'] = Wave_vis_nir[0]/10 # in nm (not Angstrom)
# 	hde_vis['NAXIS1'] = len(Flux_vis_nir)
	

# 	# Read out the fucking file_out

# 	path_out = '/'.join(path_NIR.split('/')[:-2])+'/VIS_NIR/'
# 	file_out = '%s_VIS_NIR_avwCombined_OBx5_sig5_rebin%snm_opt_V1_right.fits' % (target,hde_vis['CDELT1'])
# 	print path_out+file_out
# 	# raise

# 	if not os.path.exists(path_out+file_out):
# 	    # Read out flux array
# 	    pf.writeto(path_out+file_out, Flux_vis_nir, hdf_vis)
	    
# 	    # Read out error array
# 	    pf.append(path_out+file_out, Errs_vis_nir, hde_vis)

# 	else:
# 	    print 'file already exists'




# 	## Plot ###
# 	plt.scatter(Wave_vis_overlap,Flux_overlap_vis,marker='s',s=70,color='magenta')
# 	plt.scatter(Wave_nir_overlap,Flux_overlap_nir,marker='x',s=70,color='black')
# 	plt.scatter(Wave_nir_overlap,Flux_avw,marker='o',s=70,color='green')

# 	plt.plot(Wave_vis_nir,Flux_vis_nir+0.2e-17,color='black',label='Flux+offset') # offset is 
# 	# 

# 	plt.scatter(Wave_vis-Match_diff,Flux_vis)
# 	plt.scatter(Wave_nir,Flux_nir,color='r')
# 	plt.ylim([-0.2e-17,0.8e-17])
# 	plt.xlim([5000,21200])
# 	plt.title('%s' % target)
# 	plt.show()
# '''

# ###########################################################
# ###########################################################















# ###########################################################################
# ### Read in new file and compare wavelength grid to the constructed one ###
# '''
# path_vis_nir = '/Volumes/DataDrive/X-shooter/P93/Data/90676/Combined_OBs/VIS_NIR/90676_VIS_NIR_avwCombined_OBx5_sig5_rebin1.2_opt_V1.fits'
# f = pf.open(path_vis_nir)
# hd = f[0].header

# Flux = f[0].data
# Errs = f[1].data
# Wave = (hd['CRVAL1'] + (hd['CRPIX1'] - 1 + np.arange(hd['NAXIS1']))*hd['CDELT1'])*10 # Angstrom

# print Wave
# print Wave_vis_nir

# F_d = Flux-Flux_vis_nir
# print F_d.all() == 0

# plt.plot(Wave,F_d)
# plt.ylim([-1e-18,1e-18])
# plt.show()
# '''

# ###########################################################################
# ###########################################################################


# #############################################
# ### Plot with bin width instead of points ###
# '''
# binsize = np.float(path_NIR.split('medrebin')[-1].split('_')[0])

# Flux_nir_bin = np.zeros(len(Flux)*binsize)
# Wave_nir_bin = np.zeros(len(Flux)*binsize)
# for i in range(len(Flux)):
# 	lower = int(i*binsize)
# 	upper = int((i+1)*binsize)

# 	Wave_nir_bin[lower:upper] = Wave_nir[i]
# 	Flux_nir_bin[lower:upper] = Flux_nir[i]



# plt.plot(Wave_nir_bin,Flux_nir_bin)
# plt.show()
# '''

# #############################################
# #############################################





























