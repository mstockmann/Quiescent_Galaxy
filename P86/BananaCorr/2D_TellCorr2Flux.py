#!/usr/local/bin/python

#-*- coding: utf-8 -*-

from __future__ import division
import glob
import pyfits as pf
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import os

#################################################################
### Correcting each sub science frame for telluric absorption ###


# Set parameters
# target = 'CP-540713'
# target = 'CP-561356' 
# target = 'CP-1243752'
# target = 'CP-1291751'


Spath = glob.glob('../../../X-shooter/P86_COSMOS_1/Data/%s/OB*/NIR_BananaCorr/Reduction/Output_*/xsh_scired_slit_nod_SKY_SLIT_FLUX_MERGE2D_NIR_NIR_x_.fits' % (target))



Tpath = glob.glob('../../../X-shooter/P86_COSMOS_1/Data/%s/OB*/NIR/TellCorr_%s_OB*' % (target,target))

# print Spath
# print Tpath
# print len(Spath), len(Tpath)
# raise




nn=0
for ii in range(len(Tpath)):

	no_output_per_OB = len(glob.glob('/'.join(Tpath[ii].split('/')[:-1])+'/Reduction/Output*'))-1
	print 'Number of output/OB: %s' % no_output_per_OB


	path_tell   = Tpath[ii]
	data_tell   = np.genfromtxt(path_tell)
	Wave        = data_tell[:,0]
	Tell_obs    = data_tell[:,1]
	Tell_model  = data_tell[:,2]

	trans = (Tell_obs/Tell_model)

	# We give negative ratio a zero weight (trans[i] = 1)
	for i in range(len(trans)):
		if trans[i] <= 0:
			trans[i] = 1

	
	
	for jj in range(no_output_per_OB):

		OB_no = nn #jj+(ii)*2
		nn += 1

		path_sci  = Spath[OB_no]

		
		### Read-in data ###
		f = pf.open(path_sci)
		#print f.info()
		hd_0 = f[0].header#; print hd
		hd_1 = f[1].header
		hd_2 = f[2].header
		Flux_arr = f[0].data
		Err_arr  = f[1].data
		Qual_arr = f[2].data

		# Telluric Corrected 2D Flux Image
		Flux_arr_trans = np.zeros(shape=(len(Flux_arr[:,0]),len(Flux_arr[0,:])))
		Err_arr_trans  = np.zeros(shape=(len(Err_arr[:,0]),len(Err_arr[0,:])))
		for i in range(len(Flux_arr[:,0])):
			Flux_arr_trans[i,:] = Flux_arr[i,:]/trans
			Err_arr_trans[i,:]  = Err_arr[i,:]/trans

		# Plot S/N values to test if the transmission array are creating large deficits.
		# mid_trace = int(len(Flux_arr_trans[:,0])/2)
		# S2N = (Flux_arr_trans[mid_trace,:]*1e19)/(Err_arr_trans[mid_trace,:]*1e19)
		# plt.plot(S2N)
		# plt.show()

		# Read out the Telluric corrected Flux spectrum
		path = '/'.join(path_sci.split('/')[:9])+'/'
		output_name = 'sci_tellcorr_flux_merged2d_nir_%s_%s_%s.fits' % (target,path_sci.split('/')[7],(jj+1))
		print 'reading out: ', path+output_name


		if not os.path.exists(path+output_name):
			pf.writeto(path+output_name, Flux_arr_trans, hd_0)
			pf.append(path+output_name, Err_arr_trans, hd_1)
			pf.append(path+output_name, Qual_arr, hd_2)

		else:
			os.system("rm %s" % (path+output_name))
			pf.writeto(path+output_name, Flux_arr_trans, hd_0)
			pf.append(path+output_name, Err_arr_trans, hd_1)
			pf.append(path+output_name, Qual_arr, hd_2)
			




#################################################################
#################################################################








#---# Testing #---#
'''
# Flux (uncorrected) median - collapsed in the vertical direction
Flux_1d = np.zeros(len(Flux_arr[0,:]))
for i in range(len(Flux_arr[0,:])):
	Flux_1d[i] = np.median(Flux_arr[:,i])

# Flux (corrected) median - collapsed in the vertical direction
Flux_1d_trans = np.zeros(len(Flux_arr[0,:]))
Err_1d_trans  = np.zeros(len(Err_arr[0,:]))
for i in range(len(Flux_arr_trans[0,:])):
	Flux_1d_trans[i] = np.median(Flux_arr_trans[:,i])
	Err_1d_trans[i] = np.median(Err_arr_trans[:,i])


plt.errorbar(Wave,Flux_1d_trans,yerr=Err_1d_trans)
#plt.plot(Wave,Flux_1d_trans,'black')
#plt.plot(Wave,Flux_1d,'r')
plt.show()


# Testing the first 4 pixels
print '\nFlux (uncorrected):'
print '[',Flux_arr[0,0],']','[', Flux_arr[0,1],']'
print '[',Flux_arr[1,0],']','[', Flux_arr[1,1],']\n'

print 'Flux (corrected):'
print '[',Flux_arr_trans[0,0],']','[', Flux_arr_trans[0,1],']'
print '[',Flux_arr_trans[1,0],']','[', Flux_arr_trans[1,1],']\n'

print 'Flux (uncorrected) / Transmission function:'
print '[',Flux_arr[0,0]/trans[0],']','[', Flux_arr[0,1]/trans[1],']'
print '[',Flux_arr[1,0]/trans[0],']','[', Flux_arr[1,1]/trans[1],']\n'
'''
#---#---#---#---#---#



