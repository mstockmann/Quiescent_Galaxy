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
import sys

sys.path.insert(0, '/Volumes/DataDrive/X-shooter/P93/Codes')
from P93_functions import *

####################################################





####################################################

Nspec = glob.glob('../../../X-shooter/P93/Data/Reduction/*/Combined_OBs/VIS/*_V1_wmrebin45.fits')



for i in range(len(Nspec)):
	target = Nspec[i].split('/')[7]
	print target

	### Rebinned ###
	# CDELT1 = 0.9
	#path_NIR = '/Volumes/DataDrive/X-shooter/P93/Data/%s/Combined_OBs/NIR/*V3_wmrebin15_opt.fits' % (target)
	#path_VIS = '/Volumes/DataDrive/X-shooter/P93/Data/%s/Combined_OBs/VIS/*V1_wmrebin45_opt.fits' % (target)

	# CDELT1 = 0.06
	path_NIR = glob.glob('../../../X-shooter/P93/Data/Reduction/%s/Combined_OBs/NIR_corr/*_V1_NIR_corr_opt.fits' % (target))

	path_VIS = glob.glob('../../../X-shooter/P93/Data/Reduction/%s/Combined_OBs/VIS/*sig5_V1_wmrebin3_opt.fits' % (target))

	Wave_nir, Flux_nir, Errs_nir, hdf_nir, hde_nir = read_in_1d_fits(path_NIR[0])	
	Wave_vis, Flux_vis, Errs_vis, hdf_vis, hde_vis = read_in_1d_fits(path_VIS[0])

	print 'Object: %s, NIR_CD1 / VIS_CD1 = %s / %s' % (target,hdf_nir['CDELT1'], hdf_vis['CDELT1'])

	print len(Wave_nir), len(Flux_nir)


	plt.plot(Wave_nir,Flux_nir)
	plt.show()

	# Create a new collective wavelength grid
	# Note: We tried the easy way using np.union1d() but it left out some of the overlaps (??)
	Wave_vis_match_point = Wave_vis[np.searchsorted(Wave_vis,Wave_nir[0])]
	Match_diff = Wave_vis_match_point-Wave_nir[0]
	Wave_vis_match = Wave_vis-Match_diff


	ID_vis_left = np.where(Wave_vis_match == Wave_nir[0])[0]
	ID_nir_right = np.where(Wave_nir == Wave_vis_match[-1])[0]

	Wave_vis_overlap = Wave_vis_match[ID_vis_left:]
	Wave_nir_overlap = Wave_nir[:ID_nir_right+1]

	Wave_vis_nir = np.concatenate((Wave_vis_match[:ID_vis_left],Wave_nir_overlap,Wave_nir[ID_nir_right+1:]),axis=0)


	# Create a new collective Flux grid
	Flux_overlap_vis = Flux_vis[ID_vis_left:]
	Flux_overlap_nir = Flux_nir[:ID_nir_right+1]

	Errs_overlap_vis = Errs_vis[ID_vis_left:]
	Errs_overlap_nir = Errs_nir[:ID_nir_right+1]
	



	# Combine overlapping fluxes (When plotting them they look similar to highest points, but they are different when zooming)
	Flux_avw = np.zeros(len(Flux_overlap_vis))
	Errs_avw = np.zeros(len(Flux_overlap_vis))

	for i in range(len(Flux_overlap_vis)):
		err_tmp = np.array([Errs_overlap_vis[i],Errs_overlap_nir[i]])
		flux_tmp = np.array([Flux_overlap_vis[i],Flux_overlap_nir[i]])

		w_err = 1/(err_tmp**2)
		Flux_avw[i] = np.sum(w_err*flux_tmp)/np.sum(w_err)
		Errs_avw[i] = 1/np.sqrt(np.sum(w_err)) 


	Flux_vis_nir = np.concatenate((Flux_vis[:ID_vis_left],Flux_avw,Flux_nir[ID_nir_right+1:]),axis=0)
	Errs_vis_nir = np.concatenate((Errs_vis[:ID_vis_left],Errs_avw,Errs_nir[ID_nir_right+1:]),axis=0)



	# We have change the length (NAXIS1) and the start point (CRVAL), but we keep the binning
	hdf_vis['CRVAL1'] = Wave_vis_nir[0]/10 # in nm (not Angstrom)
	hdf_vis['NAXIS1'] = len(Flux_vis_nir)

	hde_vis['CRVAL1'] = Wave_vis_nir[0]/10 # in nm (not Angstrom)
	hde_vis['NAXIS1'] = len(Flux_vis_nir)
	

	# Read out the fucking file_out
	path_out = '/'.join(path_NIR.split('/')[:-2])+'/VIS_NIR/'
	file_out = '%s_VIS_NIR_avwCombined_OBx5_sig5_rebin%snm_opt_V1.fits' % (target,hde_vis['CDELT1'])
	print path_out+file_out
	raise

	if not os.path.exists(path_out+file_out):
	    # Read out flux array
	    pf.writeto(path_out+file_out, Flux_vis_nir, hdf_vis)
	    
	    # Read out error array
	    pf.append(path_out+file_out, Errs_vis_nir, hde_vis)

	else:
	    print 'file already exists'

	## Plot ###
	plt.scatter(Wave_vis_overlap,Flux_overlap_vis,marker='s',s=70,color='magenta')
	plt.scatter(Wave_nir_overlap,Flux_overlap_nir,marker='x',s=70,color='black')
	plt.scatter(Wave_nir_overlap,Flux_avw,marker='o',s=70,color='green')

	plt.plot(Wave_vis_nir,Flux_vis_nir+1.9e-17,color='black',label='Flux+offset') # offset is 
	# 

	plt.scatter(Wave_vis-Match_diff,Flux_vis)
	plt.scatter(Wave_nir,Flux_nir,color='r')
	plt.ylim([-0.2e-17,0.8e-17])
	plt.xlim([5000,21200])
	plt.title('%s' % target)
	plt.show()





















###########################################################
### Combine NIR and VIS for 105842 left and right trace ###

'''
Nspec = glob.glob('/Volumes/DataDrive/X-shooter/P93/Data/*/Combined_OBs/VIS/*_V1_wmrebin45.fits')


for i in range(1):#len(Nspec)):
	target = 105842 #Nspec[i].split('/')[6]
	print target


	### Rebinned ###
	# CDELT1 = 1.2
	path_NIR = '/Volumes/DataDrive/X-shooter/P93/Data/%s/Combined_OBs/NIR/*V3_wmrebin15_opt.fits' % (target)
	Wave_nir, Flux_nir, Errs_nir, hdf_nir, hde_nir = read_in_1d_fits(glob.glob(path_NIR)[0])

	# CDELT1 = 1.2
	path_VIS = '/Volumes/DataDrive/X-shooter/P93/Data/%s/Combined_OBs/VIS/*V1_wmrebin45_opt_right.fits' % (target)
	Wave_vis, Flux_vis, Errs_vis, hdf_vis, hde_vis = read_in_1d_fits(glob.glob(path_VIS)[0])


	# Create a new collective wavelength grid
	# Note: We tried the easy way using np.union1d() but it left out some of the overlaps (??)
	Wave_vis_match_point = Wave_vis[np.searchsorted(Wave_vis,Wave_nir[0])]
	Match_diff = Wave_vis_match_point-Wave_nir[0]
	Wave_vis_match = Wave_vis-Match_diff


	ID_vis_left = np.where(Wave_vis_match == Wave_nir[0])[0]
	ID_nir_right = np.where(Wave_nir == Wave_vis_match[-1])[0]

	Wave_vis_overlap = Wave_vis_match[ID_vis_left:]
	Wave_nir_overlap = Wave_nir[:ID_nir_right+1]

	Wave_vis_nir = np.concatenate((Wave_vis_match[:ID_vis_left],Wave_nir_overlap,Wave_nir[ID_nir_right+1:]),axis=0)


	# Create a new collective Flux grid
	Flux_overlap_vis = Flux_vis[ID_vis_left:]
	Flux_overlap_nir = Flux_nir[:ID_nir_right+1]

	Errs_overlap_vis = Errs_vis[ID_vis_left:]
	Errs_overlap_nir = Errs_nir[:ID_nir_right+1]
	 

	# Combine overlapping fluxes (When plotting them they look similar to highest points, but they are different when zooming)
	Flux_avw = np.zeros(len(Flux_overlap_vis))
	Errs_avw = np.zeros(len(Flux_overlap_vis))

	for i in range(len(Flux_overlap_vis)):
		err_tmp = np.array([Errs_overlap_vis[i],Errs_overlap_nir[i]])
		flux_tmp = np.array([Flux_overlap_vis[i],Flux_overlap_nir[i]])

		w_err = 1/(err_tmp**2)
		Flux_avw[i] = np.sum(w_err*flux_tmp)/np.sum(w_err)
		Errs_avw[i] = 1/np.sqrt(np.sum(w_err)) 


	Flux_vis_nir = np.concatenate((Flux_vis[:ID_vis_left],Flux_avw,Flux_nir[ID_nir_right+1:]),axis=0)
	Errs_vis_nir = np.concatenate((Errs_vis[:ID_vis_left],Errs_avw,Errs_nir[ID_nir_right+1:]),axis=0)



	# We have change the length (NAXIS1) and the start point (CRVAL), but we keep the binning
	hdf_vis['CRVAL1'] = Wave_vis_nir[0]/10 # in nm (not Angstrom)
	hdf_vis['NAXIS1'] = len(Flux_vis_nir)

	hde_vis['CRVAL1'] = Wave_vis_nir[0]/10 # in nm (not Angstrom)
	hde_vis['NAXIS1'] = len(Flux_vis_nir)
	

	# Read out the fucking file_out

	path_out = '/'.join(path_NIR.split('/')[:-2])+'/VIS_NIR/'
	file_out = '%s_VIS_NIR_avwCombined_OBx5_sig5_rebin%snm_opt_V1_right.fits' % (target,hde_vis['CDELT1'])
	print path_out+file_out
	# raise

	if not os.path.exists(path_out+file_out):
	    # Read out flux array
	    pf.writeto(path_out+file_out, Flux_vis_nir, hdf_vis)
	    
	    # Read out error array
	    pf.append(path_out+file_out, Errs_vis_nir, hde_vis)

	else:
	    print 'file already exists'




	## Plot ###
	plt.scatter(Wave_vis_overlap,Flux_overlap_vis,marker='s',s=70,color='magenta')
	plt.scatter(Wave_nir_overlap,Flux_overlap_nir,marker='x',s=70,color='black')
	plt.scatter(Wave_nir_overlap,Flux_avw,marker='o',s=70,color='green')

	plt.plot(Wave_vis_nir,Flux_vis_nir+0.2e-17,color='black',label='Flux+offset') # offset is 
	# 

	plt.scatter(Wave_vis-Match_diff,Flux_vis)
	plt.scatter(Wave_nir,Flux_nir,color='r')
	plt.ylim([-0.2e-17,0.8e-17])
	plt.xlim([5000,21200])
	plt.title('%s' % target)
	plt.show()
'''

###########################################################
###########################################################















###########################################################################
### Read in new file and compare wavelength grid to the constructed one ###
'''
path_vis_nir = '/Volumes/DataDrive/X-shooter/P93/Data/90676/Combined_OBs/VIS_NIR/90676_VIS_NIR_avwCombined_OBx5_sig5_rebin1.2_opt_V1.fits'
f = pf.open(path_vis_nir)
hd = f[0].header

Flux = f[0].data
Errs = f[1].data
Wave = (hd['CRVAL1'] + (hd['CRPIX1'] - 1 + np.arange(hd['NAXIS1']))*hd['CDELT1'])*10 # Angstrom

print Wave
print Wave_vis_nir

F_d = Flux-Flux_vis_nir
print F_d.all() == 0

plt.plot(Wave,F_d)
plt.ylim([-1e-18,1e-18])
plt.show()
'''

###########################################################################
###########################################################################


#############################################
### Plot with bin width instead of points ###
'''
binsize = np.float(path_NIR.split('medrebin')[-1].split('_')[0])

Flux_nir_bin = np.zeros(len(Flux)*binsize)
Wave_nir_bin = np.zeros(len(Flux)*binsize)
for i in range(len(Flux)):
	lower = int(i*binsize)
	upper = int((i+1)*binsize)

	Wave_nir_bin[lower:upper] = Wave_nir[i]
	Flux_nir_bin[lower:upper] = Flux_nir[i]



plt.plot(Wave_nir_bin,Flux_nir_bin)
plt.show()
'''

#############################################
#############################################





























