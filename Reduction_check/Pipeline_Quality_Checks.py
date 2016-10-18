
#!/usr/local/bin/python

#-*- coding: utf-8 -*-

from __future__ import division
import glob
import pyfits as pf
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import os
import string
import subprocess
import sys
import matplotlib.gridspec as gridspec

sys.path.insert(0, '/Volumes/DataDrive/Quiescent_Galaxy')
import Stokky as st

################################ Programs #################################


def QC_predict(period,target,reduct_type,No_OB):

	path_predict = glob.glob('../../X-shooter/%s/Data/Reduction/%s/OB*/%s/Reduction/Output/xsh_predict_FMTCHK_RESID_TAB_LINES_NIR_NIR_x_.fits' % (period,target,reduct_type))

	print 'PREDICT: %s. Number of exposures / predict output: %s / %s' % (target,len(path_predict),No_OB)



	path_arc_lamp_pinhole = glob.glob('../../X-shooter/P93/Data/Reduction/105842/OB*/data_NIR/XSHOO.2014-04-13T11:14:01.130.fits')


	g = pf.open(path_arc_lamp_pinhole[0])
	Counts = g[0].data

	for i in range(len(path_predict)):
		
		f = pf.open(path_predict[i])

		Wavelength = np.array((f[1].data)['Wavelength'])
		ResXmodel  = np.array((f[1].data)['ResidXmodel'])
		ResYmodel  = np.array((f[1].data)['ResidYmodel'])
		Flag       = np.array((f[1].data)['Flag'])

		plt.figure(figsize=(11,9))
		gs = gridspec.GridSpec(3, 2)

		plt.subplot(gs[:2, :])
		plt.imshow(Counts,vmin=40,vmax=450,cmap='gnuplot2')
		plt.axis([0,2000,0,1000])
		plt.annotate('Example frame',xy=(100,900),color='white')
		plt.title('xsh_predict')


		plt.subplot(gs[-1, 0])
		plt.scatter(Wavelength,ResXmodel,color='b',s=5)
		plt.plot(Wavelength,Flag,color='r')
		plt.title('Line X residuals')
		plt.ylabel('Residuals [pix]')
		plt.xlabel('Wavelength [nm]')

		plt.subplot(gs[-1,-1])
		plt.scatter(Wavelength,ResYmodel,color='b',s=1,label='Residuals')
		plt.plot(Wavelength,Flag,color='r',label='Flag')
		plt.title('Line Y residuals')
		plt.ylabel('Residuals [pix]')
		plt.xlabel('Wavelength [nm]')

		plt.legend(loc='upper center',fontsize=10)

		path_out = '/'.join(path_predict[i].split('/')[:-3])+'/Pipeline_Quality_Checks'
		if not os.path.exists(path_out):
			os.system('mkdir %s' % path_out)
		plt.savefig(path_out+'/1-QC_predict.png')
	plt.close("all")



def QC_orderpos(period,target,reduct_type,No_OB):
	path_orderpos = glob.glob('../../X-shooter/%s/Data/Reduction/%s/OB*/%s/Reduction/Output/xsh_orderpos_ORDERPOS_RESID_TAB_NIR_NIR_x_.fits' % (period,target,reduct_type))

	print 'ORDERPOS: %s. Number of exposures / orderpos output: %s / %s' % (target,len(path_orderpos),No_OB)

	path_arc_lamp_pinhole = glob.glob('../../X-shooter/P93/Data/Reduction/105842/OB*/data_NIR/XSHOO.2014-04-13T11:14:01.130.fits')
	g = pf.open(path_arc_lamp_pinhole[0])
	Counts = g[0].data

	for i in range(len(path_orderpos)):	
		f = pf.open(path_orderpos[i])

		X = np.array((f[1].data)['X'])
		Y = np.array((f[1].data)['Y'])
		RESX = np.array((f[1].data)['RESX'])

		plt.figure(figsize=(11,9))
		gs = gridspec.GridSpec(2, 2)

		plt.subplot(gs[:1, :])
		plt.imshow(Counts,vmin=40,vmax=450,cmap='coolwarm')
		plt.axis([0,2000,0,1000])
		plt.annotate('Example frame',xy=(100,900),color='white')
		plt.title('xsh_orderpos')
		plt.scatter(Y,X,s=0.2)


		plt.subplot(gs[-1, 0])
		ax1 = plt.gca()
		plt.scatter(X,RESX,color='b',s=1)
		# plt.plot(Wavelength,Flag,color='r')
		plt.title('Trace X Residuals vs X')
		plt.ylabel('X Residuals [pix]')
		plt.xlabel('X [pix]')
		ax1.grid('on')
		ax1.ticklabel_format(useOffset=False)
		plt.xlim([0,max(X)])

		plt.subplot(gs[-1,-1])
		ax2 = plt.gca()
		plt.scatter(Y,RESX,color='b',s=1,label='Residuals')
		# plt.plot(Wavelength,Flag,color='r',label='Flag')
		plt.title('Trace X Residuals vs Y')
		plt.ylabel('X Residuals [pix]')
		plt.xlabel('Y [pix]')
		ax2.grid('on')
		ax2.ticklabel_format(useOffset=False)

		plt.xlim([0,max(Y)])

		plt.legend(loc='upper center',fontsize=10)

		# plt.show()
		# raise
		path_out = '/'.join(path_orderpos[i].split('/')[:-3])+'/Pipeline_Quality_Checks'
		# print path_out
		if not os.path.exists(path_out):
			os.system('mkdir %s' % path_out)
		plt.savefig(path_out+'/2-QC_orderpos.png')
	plt.close("all")




def QC_mflat(period,target,reduct_type,No_OB):
	path_mflat = glob.glob('../../X-shooter/%s/Data/Reduction/%s/OB*/%s/Reduction/Output/xsh_mflat_MFLAT_GRID_BACK_SLIT_NIR_NIR_x_.fits' % (period,target,reduct_type))

	print 'MFLAT: %s. Number of exposures / mflat output: %s / %s' % (target,len(path_mflat),No_OB)

	path_masterflat = glob.glob('../../X-shooter/%s/Data/Reduction/%s/OB*/%s/Reduction/Output/xsh_mflat_MASTER_FLAT_SLIT_NIR_NIR_x_.fits' % (period,target,reduct_type))
	g = pf.open(path_masterflat[0])
	Counts = g[0].data


	for i in range(len(path_mflat)):	
		f = pf.open(path_mflat[i])

		X = np.array((f[1].data)['X'])
		Y = np.array((f[1].data)['Y'])
		Residual = np.array((f[1].data)['Residual'])
		INT = np.array((f[1].data)['INT'])


		plt.figure(figsize=(14,10))
		gs = gridspec.GridSpec(2, 2)

		plt.subplot(gs[:1, :])
		plt.imshow(Counts.T,vmin=-1.0,vmax=3.5,cmap='YlOrBr_r')
		plt.axis([0,2000,0,1000])
		plt.annotate('Example frame',xy=(100,900),color='white')
		plt.title('xsh_mflat')
		# plt.scatter(Y,X,s=0.2)

		plt.subplot(gs[-1, 0])
		ax1 = plt.gca()
		plt.scatter(X,Residual,color='b',s=1)
		# plt.plot(Wavelength,Flag,color='r')
		plt.title('Trace X Residuals vs X')
		plt.ylabel('X Residuals [Counts]')
		plt.xlabel('X [pix]')
		ax1.grid('on')
		ax1.ticklabel_format(useOffset=False)
		plt.xlim([0,max(X)])

		plt.subplot(gs[-1,-1])
		ax2 = plt.gca()
		plt.scatter(Y,Residual,color='b',s=1,label='Residuals')
		# plt.plot(Wavelength,Flag,color='r',label='Flag')
		plt.title('Trace X Residuals vs Y')
		plt.ylabel('X Residuals [Counts]')
		plt.xlabel('Y [pix]')
		ax2.grid('on')
		ax2.ticklabel_format(useOffset=False)
		plt.xlim([0,max(Y)])
		plt.legend(loc='upper center',fontsize=10)

		path_out = '/'.join(path_mflat[i].split('/')[:-3])+'/Pipeline_Quality_Checks'
		# print path_out
		if not os.path.exists(path_out):
			os.system('mkdir %s' % path_out)
		plt.savefig(path_out+'/3-QC_mflat.png')
	plt.close("all")




def QC_2dmap(period,target,reduct_type,No_OB):
	path_2dmap = glob.glob('../../X-shooter/%s/Data/Reduction/%s/OB*/%s/Reduction/Output/xsh_2dmap_WAVE_RESID_TAB_LINES_NIR_NIR_x_.fits' % (period,target,reduct_type))

	print '2DMAP: %s. Number of exposures / 2dmap output: %s / %s' % (target,len(path_2dmap),No_OB)

	path_master_2dmap = glob.glob('../../X-shooter/%s/Data/Reduction/%s/OB*/%s/Reduction/Output/xsh_2dmap_WAVE_MAP_NIR_NIR_x_.fits' % (period,target,reduct_type))
	g = pf.open(path_master_2dmap[0])
	Counts = g[0].data

	for i in range(len(path_2dmap)):	
		f = pf.open(path_2dmap[i])

		Wavelength = np.array((f[1].data)['Wavelength'])
		ResXmodel  = np.array((f[1].data)['ResidXmodel'])
		ResYmodel  = np.array((f[1].data)['ResidYmodel'])
		Flag       = np.array((f[1].data)['Flag'])

		plt.figure(figsize=(11,9))
		#gs = gridspec.GridSpec(1, 2)
		gs = gridspec.GridSpec(2, 2)

		plt.subplot(gs[:1, :])
		plt.imshow(Counts.T,vmin=500,vmax=2000,cmap='winter')
		plt.axis([0,2040,0,1050])
		plt.annotate('Example frame',xy=(100,950),color='white')
		plt.title('xsh_orderpos')
		# plt.scatter(Y,X,s=0.2)


		plt.subplot(gs[-1, 0])
		plt.scatter(Wavelength,ResXmodel,color='b',s=5)
		plt.plot(Wavelength,Flag,color='r')
		plt.title('Line X residuals')
		plt.ylabel('Residuals [pix]')
		plt.xlabel('Wavelength [nm]')

		plt.subplot(gs[-1,-1])
		plt.scatter(Wavelength,ResYmodel,color='b',s=1,label='Residuals')
		plt.plot(Wavelength,Flag,color='r',label='Flag')
		plt.title('Line Y residuals')
		plt.ylabel('Residuals [pix]')
		plt.xlabel('Wavelength [nm]')

		plt.legend(loc='upper center',fontsize=10)

		# plt.show()
		# raise
		path_out = '/'.join(path_2dmap[i].split('/')[:-3])+'/Pipeline_Quality_Checks'
		# print path_out
		if not os.path.exists(path_out):
			os.system('mkdir %s' % path_out)
		plt.savefig(path_out+'/4-QC_2dmap.png')
	plt.close("all")



def QC_wavecal(period,target,reduct_type,No_OB):

	path_wavecal = glob.glob('../../X-shooter/%s/Data/Reduction/%s/OB*/%s/Reduction/Output/xsh_wavecal_TILT_TAB_SLIT_NIR_NIR_x_.fits' % (period,target,reduct_type))

	print 'WAVECAL: %s. Number of exposures / wavecal output: %s / %s' % (target,len(path_wavecal),No_OB)

	path_master_wavecal = glob.glob('../../X-shooter/%s/Data/Reduction/%s/OB*/%s/Reduction/Output/xsh_wavecal_ARC_SLIT_ON_NIR_NIR_x_.fits' % (period,target,reduct_type))
	g = pf.open(path_master_wavecal[0])
	Counts = g[0].data

	for i in range(len(path_wavecal)):	
		f = pf.open(path_wavecal[i])

		Wavelength = np.array((f[1].data)['Wavelength'])
		SPECRES    = np.array((f[1].data)['SPECRES'])
		# ResYmodel  = np.array((f[1].data)['ResidYmodel'])
		# Flag       = np.array((f[1].data)['Flag'])

		plt.figure(figsize=(11,9))
		#gs = gridspec.GridSpec(1, 2)
		gs = gridspec.GridSpec(2, 2)

		plt.subplot(gs[:1, :])
		plt.imshow(Counts.T,vmin=-10,vmax=80,cmap='autumn')
		plt.axis([0,2040,0,1050])
		plt.annotate('Example frame',xy=(100,950),color='white')
		plt.title('xsh_orderpos')
		# plt.scatter(Y,X,s=0.2)


		plt.subplot(gs[-1, :])
		plt.scatter(Wavelength,SPECRES,color='b',s=5)
		# plt.plot(Wavelength,Flag,color='r')
		plt.title('Line X residuals')
		plt.ylabel('Spectral Resolution')
		plt.xlabel('Wavelength [nm]')
		plt.axis([min(Wavelength)-100,max(Wavelength)+100,min(SPECRES)-1000,max(SPECRES)+1000])

		plt.hlines(4500,min(Wavelength)-100,max(Wavelength)+100,linestyle='--',label='X-Shooter NIR\n Minimum Spectral Resolution')

		plt.legend(loc='upper center',fontsize=10)


		path_out = '/'.join(path_wavecal[i].split('/')[:-3])+'/Pipeline_Quality_Checks'
		if not os.path.exists(path_out):
			os.system('mkdir %s' % path_out)
		plt.savefig(path_out+'/5-QC_wavecal.png')
	plt.close("all")


def QC_flexcomp(period,target,reduct_type,No_OB):
	path_flexcomp = glob.glob('../../X-shooter/%s/Data/Reduction/%s/OB*/%s/Reduction/Output/esorex_xsh_flexcomp_NIR.log' % (period,target,reduct_type))

	print 'FLEXCOMP: %s. Number of exposures / flexcomp output: %s / %s' % (target,len(path_flexcomp),No_OB)

	# path_master_flexcomp = glob.glob('../../X-shooter/%s/Data/Reduction/%s/OB*/%s/Reduction/Output/xsh_2dmap_WAVE_MAP_NIR_NIR_x_.fits' % (period,target,reduct_type))
	# g = pf.open(path_master_2dmap[0])
	# Counts = g[0].data

	for i in range(len(path_flexcomp)):
		data = np.genfromtxt(path_flexcomp[i],delimiter='\n',dtype=None)

		Str_SetFiles = "Mean x residual"
		Loc_SetFiles = [x for x in data if Str_SetFiles in x]
		mean_x_res = []
		mean_y_res = []
		for string in Loc_SetFiles:
			mean_x_res.append(float(string.split(';')[1].split(' ')[-1]))
			mean_y_res.append(float(string.split(';')[2].split(' ')[-1]))
			
		mean_x_res = np.array(mean_x_res)
		mean_y_res = np.array(mean_y_res)

		median_x = np.median(mean_x_res)
		median_y = np.median(mean_y_res)

		plt.figure(figsize=(12,4))
		#gs = gridspec.GridSpec(1, 2)
		gs = gridspec.GridSpec(1, 2)


		plt.subplot(gs[-1, 0])
		plt.scatter(np.arange(len(mean_x_res)),mean_x_res,color='b',s=1)
		plt.hlines(median_x,-100,500,color='r')
		plt.hlines(0.05,-100,500,color='g',linestyle='--')
		plt.title('Median of the Mean X Residuals')
		plt.ylabel('Mean X Residuals [pix]')
		plt.xlabel('N')
		plt.axis([0,len(mean_x_res),0,0.4])

		plt.subplot(gs[-1,-1])
		plt.scatter(np.arange(len(mean_y_res)),mean_y_res,color='b',s=1,label='Residuals')
		plt.hlines(median_y,-100,500,color='r',label='median')
		plt.hlines(0.1,-100,500,color='g',linestyle='--',label='Expected (minimum median residual level)')
		plt.title('Median of the Mean Y Residuals')
		plt.ylabel('Mean Y Residuals [pix]')
		plt.xlabel('N')
		plt.axis([0,len(mean_y_res),0,0.4])

		plt.legend(loc='upper center',fontsize=10)

		path_out = '/'.join(path_flexcomp[i].split('/')[:-3])+'/Pipeline_Quality_Checks'
		if not os.path.exists(path_out):
			os.system('mkdir %s' % path_out)
		plt.savefig(path_out+'/6-QC_flexcomp.png')
	plt.close("all")



def QC_response(period,target,reduct_type,No_OB):
	path_flux_cat = glob.glob('../../X-shooter/%s/Data/Reduction/%s/OB*/data*/M.XSHOOTER.2013-07-01T15:48:35.541.fits' % (period,target))

	path_response = glob.glob('../../X-shooter/%s/Data/Reduction/%s/OB*/%s/Reduction/Output/xsh_respon_slit_nod_FLUX_SLIT_FLUX_MERGE1D_NIR_NIR_x_.fits' % (period,target,reduct_type))

	print 'RESPONSE: %s. Number of OBs / Response output: %s / %s' % (target,len(path_flux_cat),len(path_response))

	for i in range(len(path_response)):

		# read in response
		f = pf.open(path_response[i])
		hd = f[1].header

		Flux_obs = f[0].data
		Errs_obs = f[1].data
		Wave_obs = 10*(hd['CRVAL1'] + (hd['CRPIX1'] - 1 + np.arange(hd['NAXIS1']))*hd['CDELT1']) # Aangstroem

		Flux_Star_name = ((f[0].header)['HIERARCH ESO OBS TARG NAME'].split(' ')[0]).split('-')[0]


		# print path_flux_cat[i] + '*.fits'
		# output = subprocess.check_output('dfits %s*.fits | fitsort PRO.CATG | grep FLUX_STD_CATALOG_NIR' % (path_flux_cat[i]), shell=True)
		# file_path = output.split(' ')[0]

		# print output
		# # raise
		# print glob.glob('../../X-shooter/%s/Data/Reduction/%s/OB*/data_NIR/M.XSHOOTER.2013-07-01T15:48:35.541.fits' % (period,target))

		f0 = pf.open(path_flux_cat[0])
		stdcatlist = list((f0[1].data)['name'])
		stdcatlist[6] = stdcatlist[6].upper()
		for kk in range(len(stdcatlist)):
			stdcatlist[kk] = stdcatlist[kk].replace(' ','')


		# Std star in catalogue (changes over OB's for some targets)
		ID_Flux_Star = stdcatlist.index(Flux_Star_name)+2 # find the observed std star ID in std star catalogue
		hd_0 = f0[ID_Flux_Star].header
		data = f0[ID_Flux_Star].data
		Wave_model = data['LAMBDA']*10 # data.field(0) can also be used
		Flux_model = data['FLUX']

		plt.figure(figsize=(10,8))
		factor = 1e-14
		plt.plot(Wave_obs,Flux_obs/factor,'b',label='Observed flux std star')
		plt.plot(Wave_model,Flux_model/factor,'r',label='Constructed model')
		plt.axis([9000,25000,-0.1,2.0])
		plt.title('xsh_response')
		plt.xlabel(r'$\lambda [\AA{}]$')
		plt.ylabel('Flux [10^-14 erg/s/cm2/A]')
		plt.legend(loc='upper center',fontsize=20)

		path_out = '/'.join(path_response[i].split('/')[:-3])+'/Pipeline_Quality_Checks'
		if not os.path.exists(path_out):
			os.system('mkdir %s' % path_out)
		plt.savefig(path_out+'/7-QC_response.png')
	plt.close("all")




#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@
#@#@#@# Science Frames #@#@#@#
#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@

# target = 105842
# target = 108899
# target = 90676
# target = 155853
# target = 171060
# target = 171687
# target = 230929
# target = 239220
# target = 250513
target = 773654

No_OB = len(glob.glob('../../X-shooter/P93/Data/Reduction/%s/OB*' % target))


period = 'P93'
reduct_type = 'NIR_BananaCorr'




###################
### xsh_predict ###
QC_predict(period,target,reduct_type,No_OB)

####################
### xsh_orderpos ###
QC_orderpos(period,target,reduct_type,No_OB)

#################
### xsh_mflat ###
QC_mflat(period,target,reduct_type,No_OB)

#################
### xsh_2dmap ###
QC_2dmap(period,target,reduct_type,No_OB)

###################
### xsh_wavecal ###
QC_wavecal(period,target,reduct_type,No_OB)

####################
### xsh_flexcomp ###
QC_flexcomp(period,target,reduct_type,No_OB)

####################
### xsh_response ###
QC_response(period,target,reduct_type,No_OB)


################################
### scired visual inspection ###

# path = glob.glob('../../X-shooter/%s/Data/Reduction/%s/OB*/%s/Reduction/Output/xsh_scired_slit_nod_SCI_SLIT_FLUX_MERGE2D_NIR_NIR_x_.fits' % (period,target,reduct_type))

# print '%s, Number of exposures / Finished reductions: %s / %s' % (target,len(path),No_OB)
# # Visual inspection of spectra in gaia
# for i in range(len(path)):
#     os.system("gaia %s" % (path[i]))



