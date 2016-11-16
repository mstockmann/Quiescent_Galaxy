#!/usr/local/bin/python

#-*- coding: utf-8 -*-


from __future__ import division
import glob
import pyfits as pf
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from scipy.interpolate import splrep,splev
from scipy import ndimage
import mpl_toolkits.axisartist as AA
import sys
from scipy import interpolate

sys.path.insert(0, '..')
import Stokky as st


###########################################


path_spec = glob.glob('../../X-shooter/cQGsample/Spectra_analysis/5-Lick_Indices/Lick_Indices.*_V1_NIRcorr_wmrebin15_er2d_opt_stdcorr.dat')

# print len(path_spec)
# raise


# def make_Lick_arrays(path):
# 	data0 = np.genfromtxt(path_spec[0],dtype=str)
# 	arr_names = []
# 	for i in range(len(data0)):
# 		arr_names.append(data0[i][2])
# 		# print arr_names[i]
# 		exec('%s = np.zeros(shape=(len(path_spec),2))' % arr_names[i])
# 	return arr_names
# arr_names = make_Lick_arrays(path_spec[0])



data0 = np.genfromtxt(path_spec[0],dtype=str)
arr_names = []
for i in range(len(data0)):
	arr_names.append(data0[i][2])
	# print arr_names[i]
	exec('%s = np.zeros(shape=(len(path_spec),2))' % arr_names[i])


for jj in range(len(path_spec)):
	# print path_spec[jj]
	data = np.genfromtxt(path_spec[jj],dtype=str)

	for n in range(len(data)):
		for k in range(2):
			# print data[n][k]

			if data[n][k] == '-NaN' or data[n][k] == '**********':
				exec('%s[%s,%s] = %s' % (arr_names[n],jj,k, -99))
			# elif data[n][k] == '**********':
			# 	exec('%s[%s,%s] = %s' % (arr_names[n],jj,k, -99))
			else:
				exec('%s[%s,%s] = %s' % (arr_names[n],jj,k, data[n][k]))


plt.figure(figsize=(12,10))
plt.subplots_adjust(hspace=0.0,wspace=0)

xpl = 4
ypl = 4


################
### Column 1 ###

x_range = [[0.9,3.1],[-1.0,1.5],[-1.0,1.5],[0,6]]
y_range = [[-8,20],[-8,13],[-8,30],[-8,40]]

x_axis = D4000_n
y_axis = [LICK_HB,LICK_HD_F,LICK_HG_F,Hd_Hg]

plt.subplot(4,4,1+4*0)
plt.errorbar(x_axis[:,0], y_axis[0][:,0],xerr=x_axis[:,1],yerr=y_axis[0][:,1],linestyle='None',marker='s',ecolor='b',color='r')
plt.ylabel(r'H$\beta\ [\AA{}]$')
plt.axis([x_range[0][0],x_range[0][1],y_range[0][0], y_range[0][1]])
plt.tick_params(length = 2)



# plt.subplot(4,3,4)
plt.subplot(4,4,1+4*1)
plt.errorbar(x_axis[:,0], y_axis[1][:,0],xerr=x_axis[:,1],yerr=y_axis[1][:,1],linestyle='None',marker='s',ecolor='b',color='r')

plt.ylabel(r'H$\delta_{F}\ [\AA{}]$')
plt.axis([x_range[0][0],x_range[0][1],y_range[1][0], y_range[1][1]])
plt.tick_params(length = 2)

# plt.subplot(4,3,7)
plt.subplot(4,4,1+4*2)
plt.errorbar(x_axis[:,0], y_axis[2][:,0],xerr=x_axis[:,1],yerr=y_axis[2][:,1],linestyle='None',marker='s',ecolor='b',color='r')

plt.ylabel(r'H$\gamma_{F}\ [\AA{}]$')
plt.yticks([0,10,20,30])
plt.axis([x_range[0][0],x_range[0][1],y_range[2][0], y_range[2][1]])
plt.tick_params(length = 2)

# plt.subplot(4,3,10)
plt.subplot(4,4,1+4*3)

plt.errorbar(x_axis[:,0], y_axis[3][:,0],xerr=x_axis[:,1],yerr=y_axis[3][:,1],linestyle='None',marker='s',ecolor='b',color='r')
plt.xlabel('$\mathrm{D}_{\mathrm{n}}4000$',labelpad=10)
plt.ylabel(r'$H\delta_{F} + H\gamma_{F}\ [\AA{}]$')
plt.axis([x_range[0][0],x_range[0][1],y_range[3][0], y_range[3][1]])

plt.tick_params(length = 2)


################
### Column 2 ###

x_axis = Mg1Fe
y_axis = [LICK_HB,LICK_HD_F,LICK_HG_F,Hd_Hg]

plt.subplot(4,4,2+4*0)
plt.errorbar(x_axis[:,0], y_axis[0][:,0],xerr=x_axis[:,1],yerr=y_axis[0][:,1],linestyle='None',marker='s',ecolor='b',color='r')

# plt.ylabel(r'$H\delta_{F} + H\gamma_{F}\ [\AA{}]$')
plt.axis([x_range[1][0],x_range[1][1],y_range[0][0], y_range[0][1]])
plt.tick_params(length = 2)
plt.yticks([])
plt.xticks([])

plt.subplot(4,4,2+4*1)
plt.errorbar(x_axis[:,0], y_axis[1][:,0],xerr=x_axis[:,1],yerr=y_axis[1][:,1],linestyle='None',marker='s',ecolor='b',color='r')
# plt.ylabel(r'$H\delta_{F} + H\gamma_{F}\ [\AA{}]$')
plt.axis([x_range[1][0],x_range[1][1],y_range[1][0], y_range[1][1]])
plt.tick_params(length = 2)
plt.yticks([])
plt.xticks([])

plt.subplot(4,4,2+4*2)
plt.errorbar(x_axis[:,0], y_axis[2][:,0],xerr=x_axis[:,1],yerr=y_axis[2][:,1],linestyle='None',marker='s',ecolor='b',color='r')
# plt.ylabel(r'$H\delta_{F} + H\gamma_{F}\ [\AA{}]$')
plt.axis([x_range[1][0],x_range[1][1],y_range[2][0], y_range[2][1]])
plt.tick_params(length = 2)
plt.yticks([])
plt.xticks([])

plt.subplot(4,4,2+4*3)
plt.errorbar(x_axis[:,0], y_axis[3][:,0],xerr=x_axis[:,1],yerr=y_axis[3][:,1],linestyle='None',marker='s',ecolor='b',color='r')
plt.xlabel('$\mathrm{Mg}_{1}\mathrm{Fe}$',labelpad=10)

plt.axis([x_range[1][0],x_range[1][1],y_range[3][0], y_range[3][1]])
plt.tick_params(length = 2)
plt.yticks([])
plt.xticks([-0.5,0,0.5,1.0,1.5])


################
### Column 3 ###

x_axis = Mg2Fe
y_axis = [LICK_HB,LICK_HD_F,LICK_HG_F,Hd_Hg]

plt.subplot(4,4,3+4*0)
plt.errorbar(x_axis[:,0], y_axis[0][:,0],xerr=x_axis[:,1],yerr=y_axis[0][:,1],linestyle='None',marker='s',ecolor='b',color='r')

# plt.ylabel(r'$H\delta_{F} + H\gamma_{F}\ [\AA{}]$')
plt.axis([x_range[1][0],x_range[1][1],y_range[1][0], y_range[1][1]])
plt.tick_params(length = 2)
plt.yticks([])
plt.xticks([])

plt.subplot(4,4,3+4*1)
plt.errorbar(x_axis[:,0], y_axis[1][:,0],xerr=x_axis[:,1],yerr=y_axis[1][:,1],linestyle='None',marker='s',ecolor='b',color='r')
# plt.ylabel(r'$H\delta_{F} + H\gamma_{F}\ [\AA{}]$')
plt.axis([x_range[1][0],x_range[1][1],y_range[1][0], y_range[1][1]])

plt.tick_params(length = 2)
plt.yticks([])
plt.xticks([])

plt.subplot(4,4,3+4*2)
plt.errorbar(x_axis[:,0], y_axis[2][:,0],xerr=x_axis[:,1],yerr=y_axis[2][:,1],linestyle='None',marker='s',ecolor='b',color='r')
# plt.ylabel(r'$H\delta_{F} + H\gamma_{F}\ [\AA{}]$')
plt.axis([x_range[1][0],x_range[1][1],y_range[2][0], y_range[2][1]])
plt.tick_params(length = 2)
plt.yticks([])
plt.xticks([])


plt.subplot(4,4,3+4*3)
plt.errorbar(x_axis[:,0], y_axis[3][:,0],xerr=x_axis[:,1],yerr=y_axis[3][:,1],linestyle='None',marker='s',ecolor='b',color='r')
plt.xlabel('$\mathrm{Mg}_{2}\mathrm{Fe}$',labelpad=10)
# plt.ylabel(r'$H\delta_{F} + H\gamma_{F}\ [\AA{}]$')
plt.axis([x_range[1][0],x_range[1][1],y_range[3][0], y_range[3][1]])
plt.tick_params(length = 2)
plt.yticks([])
plt.xticks([-0.5,0,0.5,1.0,1.5])



################
### Column 3 ###

x_axis = MgFe_prime
y_axis = [LICK_HB,LICK_HD_F,LICK_HG_F,Hd_Hg]

plt.subplot(4,4,4+4*0)
plt.errorbar(x_axis[:,0], y_axis[0][:,0],xerr=x_axis[:,1],yerr=y_axis[0][:,1],linestyle='None',marker='s',ecolor='b',color='r')

# plt.ylabel(r'$H\delta_{F} + H\gamma_{F}\ [\AA{}]$')
plt.axis([x_range[3][0],x_range[3][1],y_range[1][0], y_range[1][1]])
plt.tick_params(length = 2)
plt.yticks([])
plt.xticks([])

plt.subplot(4,4,4+4*1)
plt.errorbar(x_axis[:,0], y_axis[1][:,0],xerr=x_axis[:,1],yerr=y_axis[1][:,1],linestyle='None',marker='s',ecolor='b',color='r')
# plt.ylabel(r'$H\delta_{F} + H\gamma_{F}\ [\AA{}]$')
plt.axis([x_range[3][0],x_range[3][1],y_range[1][0], y_range[1][1]])

plt.tick_params(length = 2)
plt.yticks([])
plt.xticks([])

plt.subplot(4,4,4+4*2)
plt.errorbar(x_axis[:,0], y_axis[2][:,0],xerr=x_axis[:,1],yerr=y_axis[2][:,1],linestyle='None',marker='s',ecolor='b',color='r')
# plt.ylabel(r'$H\delta_{F} + H\gamma_{F}\ [\AA{}]$')
plt.axis([x_range[3][0],x_range[3][1],y_range[2][0], y_range[2][1]])
plt.tick_params(length = 2)
plt.yticks([])
plt.xticks([])


plt.subplot(4,4,4+4*3)
plt.errorbar(x_axis[:,0], y_axis[3][:,0],xerr=x_axis[:,1],yerr=y_axis[3][:,1],linestyle='None',marker='s',ecolor='b',color='r')
plt.xlabel("$[\mathrm{MgFe}]'$",labelpad=10)
# plt.ylabel(r'$H\delta_{F} + H\gamma_{F}\ [\AA{}]$')
plt.axis([x_range[3][0],x_range[3][1],y_range[3][0], y_range[3][1]])
plt.tick_params(length = 2)
plt.yticks([])
plt.xticks([1,2,3,4,5,6])

plt.show()
















































