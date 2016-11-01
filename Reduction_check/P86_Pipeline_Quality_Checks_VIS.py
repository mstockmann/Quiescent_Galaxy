
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

	path_predict = glob.glob('../../X-shooter/%s/Data/%s/OB*/%s/Reduction/Output/xsh_predict_FMTCHK_RESID_TAB_LINES_*.fits' % (period,target,reduct_type))


	print 'PREDICT: %s. Number of exposures / predict output: %s / %s' % (target,len(path_predict),No_OB)


	for i in range(len(path_predict)):
		
		f = pf.open(path_predict[i])

		Wavelength = np.array((f[1].data)['Wavelength'])
		ResXmodel  = np.array((f[1].data)['ResidXmodel'])
		ResYmodel  = np.array((f[1].data)['ResidYmodel'])
		Flag       = np.array((f[1].data)['Flag'])

		plt.figure(figsize=(11,9))
		gs = gridspec.GridSpec(1, 2)

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
		plt.savefig(path_out + '/1-QC-predict.png')
	plt.close("all")



def QC_orderpos(period,target,reduct_type,No_OB):
	path_orderpos = glob.glob('../../X-shooter/%s/Data/%s/OB*/%s/Reduction/Output/xsh_orderpos_ORDERPOS_RESID_TAB_*.fits' % (period,target,reduct_type))

	print 'ORDERPOS: %s. Number of exposures / orderpos output: %s / %s' % (target,len(path_orderpos),No_OB)

	for i in range(len(path_orderpos)):	
		f = pf.open(path_orderpos[i])

		X = np.array((f[1].data)['X'])
		Y = np.array((f[1].data)['Y'])
		RESX = np.array((f[1].data)['RESX'])

		plt.figure(figsize=(11,9))
		gs = gridspec.GridSpec(1, 2)

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

		path_out = '/'.join(path_orderpos[i].split('/')[:-3])+'/Pipeline_Quality_Checks'

		if not os.path.exists(path_out):
			os.system('mkdir %s' % path_out)
		plt.savefig(path_out + '/2-QC_orderpos.png')
	plt.close("all")




def QC_mflat(period,target,reduct_type,No_OB):
	path_mflat = glob.glob('../../X-shooter/%s/Data/%s/OB*/%s/Reduction/Output/xsh_mflat_*_MFLAT_GRID_BACK_SLIT_*.fits' % (period,target,reduct_type))

	print 'MFLAT: %s. Number of exposures / mflat output: %s / %s' % (target,len(path_mflat),No_OB)

	for i in range(len(path_mflat)):	
		f = pf.open(path_mflat[i])

		X = np.array((f[1].data)['X'])
		Y = np.array((f[1].data)['Y'])
		Residual = np.array((f[1].data)['Residual'])
		INT = np.array((f[1].data)['INT'])


		plt.figure(figsize=(14,10))
		gs = gridspec.GridSpec(1, 2)

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
		if not os.path.exists(path_out):
			os.system('mkdir %s' % path_out)
		plt.savefig(path_out + '/3-QC_mflat.png')
	plt.close("all")




def QC_2dmap(period,target,reduct_type,No_OB):
	path_2dmap = glob.glob('../../X-shooter/%s/Data/%s/OB*/%s/Reduction/Output/xsh_2dmap_WAVE_RESID_TAB_LINES_*.fits' % (period,target,reduct_type))

	print '2DMAP: %s. Number of exposures / 2dmap output: %s / %s' % (target,len(path_2dmap),No_OB)


	for i in range(len(path_2dmap)):	
		f = pf.open(path_2dmap[i])

		Wavelength = np.array((f[1].data)['Wavelength'])
		ResXmodel  = np.array((f[1].data)['ResidXmodel'])
		ResYmodel  = np.array((f[1].data)['ResidYmodel'])
		Flag       = np.array((f[1].data)['Flag'])

		plt.figure(figsize=(11,9))
		#gs = gridspec.GridSpec(1, 2)
		gs = gridspec.GridSpec(1, 2)

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

		if not os.path.exists(path_out):
			os.system('mkdir %s' % path_out)
		plt.savefig(path_out + '/4-QC_2dmap.png')
	plt.close("all")






#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@
#@#@#@# Science Frames #@#@#@#
#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@

# target = 'CP-540713'
# target = 'CP-561356'
# target = 'CP-1243752'
# target = 'CP-1291751'


No_OB = len(glob.glob('../../X-shooter/P86_COSMOS_1/Data/%s/OB*' % target))


period = 'P86_COSMOS_1'
reduct_type = 'VIS'




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

####################
### xsh_response ###
# is run with phase3 response function


################################
### scired visual inspection ###

# path = glob.glob('../../X-shooter/%s/Data/Reduction/%s/OB*/%s/Reduction/Output/xsh_scired_slit_nod_SCI_SLIT_FLUX_MERGE2D_NIR_NIR_x_.fits' % (period,target,reduct_type))

# print '%s, Number of exposures / Finished reductions: %s / %s' % (target,len(path),No_OB)
# # Visual inspection of spectra in gaia
# for i in range(len(path)):
#     os.system("gaia %s" % (path[i]))



