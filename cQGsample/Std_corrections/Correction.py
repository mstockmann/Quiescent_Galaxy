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

	# Phot_list = [Uf, Ue, Bf, Be, Vf, Ve, Rf, Re, If, Ie, zf, ze, Yf, Ye, Jf, Je, Hf, He, Kf, Ke, Ktf, Kte]
	Phot_list = [Uf, Ue, Bf, Be, Vf, Ve, Rf, Re, If, Ie, zf, ze, Yf, Ye, Jf, Je, Hf, He, Kf, Ke]
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
	data_filter = np.array(QSO_files)

	Trans_arr = []
	for jj in range(len(data_filter)):
		Arr_tmp = np.zeros(shape=(len(data_filter[jj]),3))
		for nn in range(len(data_filter[jj])):
			Arr_tmp[nn,0] = data_filter[jj][nn][0]
			Arr_tmp[nn,1] = data_filter[jj][nn][1] 
			Arr_tmp[nn,2] = data_filter[jj][nn][2]
		Trans_arr.append(Arr_tmp)

	return Trans_arr

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



# Phot_path = glob.glob('../../../Quiescent_Galaxy/FAST/Photometric_Cat/Output/PhotCat_xsh_all_V1_cgsA_NoSPLASH.cat')

Phot_path = glob.glob('../../../Quiescent_Galaxy/FAST/Photometric_Cat/Output/PhotCat_xsh_all_V1_cgsA_NoSPLASH_Aper3.cat')


Trans_path = glob.glob('../../../../Programs/FAST/COSMOS_filters/FILTER.RES.SWv5.R300.txt')



# read in photometry
Photometric_list = read_phot_cat(Phot_path[0])
Uf = Photometric_list[0]
Bf = Photometric_list[2]
Vf = Photometric_list[4] 
Rf = Photometric_list[6]
If = Photometric_list[8]
zf = Photometric_list[10]
Yf = Photometric_list[12]
Jf = Photometric_list[14]
Hf = Photometric_list[16]
Kf = Photometric_list[18]
# phot_array = np.array([Uf, Bf, Vf, Rf, If, zf, Yf, Jf, Hf, Kf])
# Lam_centre = np.array([3562,4458,5477,6186,7506,8962,10200,12520,16450,21470])
phot_array = np.array([Jf, Hf])
Lam_centre = np.array([12520,16450])




Lam_centre
# Read in the filter curves
# Trans_curve_list = ['CAPAK_v2/u_megaprime_sagem.res','CAPAK_v2/B_subaru.res','CAPAK_v2/V_subaru.res','CAPAK_v2/r_subaru.res','CAPAK_v2/i_subaru.res','CAPAK_v2/z_subaru.res','sw_COSMOS/WFCAM_Y.res','UltraVISTA/J_ultravista_qe_atm.bp','UltraVISTA/H_ultravista_qe_atm.bp','UltraVISTA/Ks_ultravista_qe_atm.bp']
Trans_curve_list = ['UltraVISTA/J_ultravista_qe_atm.bp','UltraVISTA/H_ultravista_qe_atm.bp']
# data_filter = read_filter_curves(Trans_path[0],'UltraVISTA/H_ultravista_qe_atm.bp')
data_filter = read_filter_curves(Trans_path[0],Trans_curve_list)




for ii in range(len(Spec_path)):

	W, F, E, M, hdf, hde, hdm = st.read_in_1d_fits(Spec_path[ii])

	# plt.plot(W,F*1e18,'black')
	# plt.plot(W,M,'green')
	# plt.show()
	# plt.scatter(Lam_centre,[Uf[0], Bf[0], Vf[0], Rf[0], If[0], zf[0], Yf[0], Jf[0], Hf[0], Kf[0]],color='r')
	# # plt.plot(data_filter[0][:,1],data_filter[0][:,2],'red')
	# plt.show()


	Synthetic_photometry = np.zeros(len(data_filter))
	for i in range(len(data_filter)):
		W_trans = data_filter[i][:,1]
		T_trans = data_filter[i][:,2]

		ID_wave = np.searchsorted(W_trans,W)
		W_synth = W_trans[ID_wave]
		T_synth = T_trans[ID_wave]

		# make synthetic mag
		Synthetic_photometry[i] = np.sum(F*T_synth)/np.sum(T_synth)

	# print phot_array[:,0]
	# raise
	# Phot_col_arr = np.zeros(len(phot_array))
	# for col in range(len(phot_array)):
	# 	Phot_col_arr[col] = phot_array[0][0]

	f_slit_loss = phot_array[:,0]/Synthetic_photometry

	# Slit loss correction using the total (aper3) H-band
	F_corr = F*f_slit_loss[1]
	print f_slit_loss



	plt.plot(W,F,color='black')
	plt.plot(W,F_corr,color='red')
	# plt.errorbar(Lam_band_av,Phot_F_target*corr_aperture,yerr=Phot_E_target,color='green',marker='s',)
	plt.axis([2500,23000,-0.2e-18,4e-18])
	plt.show()

	# # read out
	# # NIR+VIS
	# # name_out = Spec_path[ii].split('/')[-1].replace('_avwCombined_OBx5_sig5','').replace('opt_','opt_stdcorr_')
	# # path = '../../../X-shooter/cQGsample/Spectra_analysis/3-Stdcorr_slit_aper/NIR+VIS/%s' % name_out

	# # NIR
	# name_out = Spec_path[ii].split('/')[-1].replace('opt','opt_stdcorr')
	# path = '../../../X-shooter/cQGsample/Spectra_analysis/3-Stdcorr_slit_aper/NIR/%s' % name_out
	# # print path_out
	

	# # st.read_out_3arr_2dfits(path,F_corr,E,M,hdf,hde,hdm)















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





































