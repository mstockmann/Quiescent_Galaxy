#!/usr/local/bin/python

#-*- coding: utf-8 -*-

###############################################################################
###############################################################################
##                                                                           ##
##   CODE:       Velocity dispersion estimate with ppxf                      ##
##                                                                           ##
##   AUTHOR:     Mikkel Stockmann                                            ##
##   DATE:       29 July 2015                                                ##
##   Modified:   29 March 2016                                               ##
##   GROUP:      Dark Cosmology Centre                                       ##
##                                                                           ##
##   DESCRIPTION:                                                            ##
##   The code can be used to estimate the velocity dispersion, and the       ##
##   spread associated with the dispersion.                                  ##
##                                                                           ##
##   Input model: FAST best model / a range of acceptable FAST models        ##
##   Input data : Each individual spectrum                                   ##
##   Perturbation technique: Perturbing model with WGN from data             ##
##                                                                           ##
###############################################################################
###############################################################################

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

# Import ppxf modules
sys.path.insert(0, '/Users/mstockmann/PhD_DARK/Programs/ppxf')
from pPXF import ppxf
from pPXF import ppxf_util as util

sys.path.insert(0, '/Users/mstockmann/PhD_DARK/Programs/ppxf/pPXF')
from AddPackages import spec
from AddPackages import tools
from AddPackages import analysis
from AddPackages import nsdata


# Setting constants
c = 299792.458 # km/s

#------------------------------------------------#
#------------------ Functions -------------------#
#------------------------------------------------#


def readfitsfiles(path,nconst):
    f = pf.open(path)
    hd_0 = f[0].header
    F = f[0].data*nconst
    E = f[1].data*nconst
    W = (hd_0['CRVAL1'] + (hd_0['CRPIX1'] - 1 + np.arange(hd_0['NAXIS1']))*hd_0['CDELT1'])*10 # Aangstroem
    Name = path.split('/')[-1].split('_')[0]
    return W, F, E, Name

def Make_fast_model_path_list(path,N):
    PATHS = []
    for ii in range(N):
        PATHS.append(glob.glob(path % (ii+1))[0])
    return PATHS

def read_in_fast_model(path):
    data = np.genfromtxt(path)
    f = data[:,1]
    w =  data[:,0]
    return f,w 


def extract_highres_part_of_BC03_model(ffit,wfit,fobs,wobs,r1,r2):
    """ This model cut out the high resolution part of the BC03 model, and
        perform a similar cut for the observed spectra. This is both crucial so
        that the velocity scale remains constant when logrebinning and so that 
        the ppxf, that matches them in pixel space works accordingly.
    """
    id1 = np.searchsorted(wfit,r1)
    id2 = np.searchsorted(wfit,r2)
    w_model_cut = wfit[id1:id2]
    f_model_cut = ffit[id1:id2]

    idobs1 = np.searchsorted(wobs,r1)
    idobs2 = np.searchsorted(wobs,r2)
    w_obs_cut = wobs[idobs1:idobs2]
    f_obs_cut = fobs[idobs1:idobs2]

    return f_model_cut, w_model_cut, f_obs_cut, w_obs_cut


def run_logrebin(f,w,v=None,em=False): # =None
    """ Changing from constant wavelength space to velocity space: v = log(lambda)
        can be usefull when extracting velocity information from a spectrum. The flux [erg/s/cm2/A]
        and the wavelength [A] are related thus changing the wavelength to a velocity binned spectrum
        results in a correction to the flux as well to conserve the flux density. 
        Doing this accordingly we use the ppxf module util.log_rebin()

        A note is to start rebinning the smallest of all the arrays (the one with least points), and use the velocity scale
        to rebin the rest of the spectra. We do this to avoid interpolation in the array with the fewest points. The larger
        will with the fixed velocity scale implement some interpolation but using this way we get minimal interpolation.
    """
    wrange = np.array([min(w),max(w)])
    logf, logw, vscale = util.log_rebin(wrange,f,velscale=v,error_mode=em)
    return logf, logw, vscale



def Gaussian_smoothing_to_velscale(obs,err,vscale):
    """ We smooth the high resolution spectrum to fit the lower
        resolution fast fit.

        OBS: - How does it affect the smoothing when the input data
               is binned?
             - Should we filter the noise, try both?
    """

    FWHM_model = c/2000 # R = 2000. See intro in www.mpe.mpg.de/events/alvio05/contributions/bruzual@alvio05.pdf
    FWHM_gal = c/5300 # R = 5300 (X-shooter manual p.5)

    FWHM_dif = np.sqrt(FWHM_model**2 - FWHM_gal**2)
    sigma_dif_pix = FWHM_dif/2.355/vscale # (km/s -> pix)

    obs_out = ndimage.gaussian_filter1d(obs,sigma_dif_pix)
    err_out = ndimage.gaussian_filter1d(err,sigma_dif_pix)
    return obs_out, err_out, sigma_dif_pix


def restframe_log_mask(zshift):   
    Range_array = [13500,14240,17919,19641,20590]
    R_rest = Range_array/(1+zshift)
    R_logrest = np.log(R_rest)

    mask_goodpix = (logLam1 < R_logrest[0]) | (R_logrest[1] < logLam1) & (logLam1 < R_logrest[2]) | (R_logrest[3] < logLam1) & (logLam1 < R_logrest[4]) 

    gp = np.where(mask_goodpix == True)[0]
    bp = np.where(mask_goodpix == False)[0]

    return gp, bp


def systematic_velocity(obs_wave,model_wave):
    """ The relation between velocity shift and wavelength is defined: v/c=Dlam/lambda
        In order to match up the spectra along the wavelength direction we calculate
        shift in wavelength that is between them and include this in the fit.

        We have tried 3 versions:
        1) vsys = 0
        2) vsys = Dlam*c 
        3) vsys = Dlam/lam*c
        They all vary max 4 km/s


    """
    dlam = model_wave[0]-obs_wave[0]
    mean_velocity = (model_wave[0]+obs_wave[0])/2
    vsys = c*dlam/mean_velocity
    return vsys



def linear_random_array(err):
    array = np.zeros(len(err))
    for dd in range(len(array)):
        array[dd] = np.random.normal(0, err[dd],1) # 1 is number of randoms
    return array



#------------------------------------------------#
#------------------------------------------------#




###################
### Input files ###

# Read in none-binned data (slitloss and aperture corrected)
# obs_spec = '/Volumes/DataDrive/X-Shooter/P93/Data/Slitloss_aper_corr_spec_1D/V3_opt/*.fits'
obs_spec = '/Volumes/DataDrive/X-shooter/P93/Data/Slitloss_aper_corr_spec_1D/V3_wmrebin5_opt/*.fits'
Spath = glob.glob(obs_spec)


Fast_model_path = Make_fast_model_path_list('/Users/mstockmann/PhD_DARK/Programs/FAST/FAST_v1.0/xshoot_phot_spec_Filt2_fixed-z/BEST_FITS/COSMOS2015_xshP93_V2_uJy.xshoot_wm_bin15_BC03V3_S2Ncut_%s.fit',10)


spec_z = np.genfromtxt('/Users/mstockmann/PhD_DARK/Programs/ppxf/code/P93/zspec.txt')[:,-1]



### Define Arrays ###
Npert = 1000
Norm_const = 1e19


v_disp_guess = np.zeros(len(Spath))
Name_list = []
for uu in range(len(Spath)):

    ### Read in observed spectrum ###
    wl_obs, Flux, Errs, name_qso = readfitsfiles(Spath[uu],Norm_const)

    ### Read in single model ###
    fl_model, wl_model_obs = read_in_fast_model(Fast_model_path[uu])

    ### Changing to rest-frame wavelength ###
    wl = wl_obs/(1+spec_z[uu])
    wl_model = wl_model_obs/(1+spec_z[uu])

    ### Select the wavelength range of the high resolution part of BC03 model (model and obs) ###
    fl_model, wl_model, Flux, wl = extract_highres_part_of_BC03_model(fl_model,wl_model,Flux,wl,3340,7500)


    ### Logrebin ###
    # model
    model1, logLam_model, velscale_model = run_logrebin(fl_model,wl_model)
    # observed spectrum
    galaxy, logLam1, velscale = run_logrebin(Flux,wl,velscale_model)
    # observed error spectrum (with Johannes modifications)
    noise, logLam1_err, velscale_err = run_logrebin(Errs,wl,velscale_model,True)

    # Here you can add a logrebin chi2 correction ()

    # Smooth high resolution observations to BC03 model resolution ---- How is the instrum disp included? 
    galaxy, noise, sig_pix = Gaussian_smoothing_to_velscale(galaxy,noise,velscale_model)

    # Make good pixel masks in restframe log based on input redshift
    goodPixels, baadPixels = restframe_log_mask(spec_z[uu])


    #@# Run pPXF #@#
    Sum_degree = 4
    Multi_degree = 4
    sigma = 3*velscale_model
    start = [0,sigma]
    dv = 0 # This have been tested to alter the dispersion 1-4 km/s 
    
    pp = ppxf.ppxf(model1, galaxy, noise, velscale_model, start, goodpixels=goodPixels, plot=False, degree=Sum_degree, mdegree=Multi_degree,moments=2, clean=True, vsyst=dv)

    # plt.show()
    z = (spec_z[uu] + 1)*(1 + pp.sol[0]/c) - 1
    v_disp_guess[uu] = pp.sol[1]
    Name_list.append(name_qso)

    Bestfit_ppxf = pp.bestfit
    noise_ppxf = pp.noise



    # print v_disp_guess
    # print Name_list
    # raise




    #############################
    ### Perturbing the models ###

    

    # Converting model to linear space
    #bestfit_ln = np.interp(wl, np.e**logLam1, Bestfit_ppxf)

    v_pdf_list = np.zeros(Npert)
    for j in range(Npert):

        # Create random realisations of the Linear Errors
        Lin_err_pertb = linear_random_array(Errs)        

        # Logrebin the new error realisation
        noise_pert, logLam1_err, velscale_err = run_logrebin(Lin_err_pertb,wl,velscale_model,True)

        # raise 
        # Convolve the error realisation to match the ppxf bestfit model
        noise_pert = ndimage.gaussian_filter1d(noise_pert,sig_pix)
        #
        # Create the data realisation using the "true" data (ppxf bestfit model)
        # + the random generated logrebin + convolved errors.  
        Data_realisation = Bestfit_ppxf+noise_pert
        #
        ###########################################################

        ####################
        ### Running pPXF ###
        #
        #
        ## Set Parameters:
        templates_p = model1 #Bestfit_ppxf # Input model template
        galaxy_p = Data_realisation # The logrebinned, and convolved data
        noise_p = noise_pert # The logrebinned, and convolved noise
        #
        # All other parameters are kept similar to the initial fit (model1, galaxy_1, noise_1)
        #
        #
        # Run ppxf module
        #
        pp = ppxf.ppxf(templates_p, galaxy_p, noise_p, velscale, start, goodpixels=goodPixels, plot=False, degree=Sum_degree, mdegree=Multi_degree,moments=2, clean=True, vsyst=dv)
        #
        v_pdf_list[j] = pp.sol[1]
        #
        ####################


    np.savetxt('./out/Mikkel_March2016/with_pert_filter/vdisp_%s_Np%s.txt' % (name_qso,Npert),v_pdf_list,delimiter='\n')


    #
    # Calculate mean and median 
    sig_median = np.round(np.median(v_pdf_list),0)
    sig_mean = np.round(np.mean(v_pdf_list),0)
    #
    # Plot the velocity distribution
    plt.figure(uu+1)
    plt.vlines(v_disp_guess[uu],-1,100,label='Bestfit')
    y, x, _ = plt.hist(v_pdf_list,label=r'$\sigma_{\mathrm{med}}=%s, \sigma_{\mathrm{mean}}=%s $' % (sig_median,sig_mean))
    plt.xlim([np.min(v_pdf_list)-50,np.max(v_pdf_list)+50])
    plt.ylim([0,y.max()*1.1])
    plt.ylabel('PDF')
    plt.xlabel(r'$\sigma$')
    plt.legend(loc='upper right')


    plt.savefig('./out/Mikkel_March2016/vdisp_%s.png' % (name_qso))
    # plt.show()
    # Read out vel_pdf
    #outpath = 'out/PSigma_PDF_sd-1_md3_Sept2015/v_pdf_list_N%s_%s.txt' % (Npert,name_qso)
    #np.savetxt(outpath,v_pdf_list)
    #





































##########################
#### Old coding stuff ####
##########################



    # ### Logrebin chi2 correction ###
    # binsize = 100
    # len_loop = len(galaxy[:-binsize]);# print len_loop

    # w_center = np.zeros(len_loop)
    # xhisquare = np.zeros(len_loop)
    # for i in range(len_loop):
    #     w1, w2 = logLam1[i], logLam1[i+binsize]
    #     F_i = galaxy[i:i+binsize]
    #     E_i = noise[i:i+binsize]

    #     ID1 = np.searchsorted(logLam_model,w1)
    #     ID2 = np.searchsorted(logLam_model,w2)
        
    #     if not ID1:
    #         Model_i = F_i+E_i

    #     else:
    #         Model_i = model1[ID1:ID2]
    #         Model_i = np.interp(logLam1[i:i+binsize],logLam_model[ID1:ID2], Model_i)

    #     if len(F_i) == len(Model_i) == len(E_i):
    #         xhisquare[i] = np.sum((F_i-Model_i)**2/E_i**2)/binsize
    #         w_center[i] = w1+(w2-w1)/2

    # xhi_med = np.median(xhisquare)
    # rule1 = np.where(xhisquare < 3*xhi_med)
    # xhisquare_n = xhisquare[rule1]
    # w_center_n = w_center[rule1]

    # # rule2 = np.where(xhisquare < 1.7)
    # # xhisq_fitpart = xhisquare[rule2]
    # # logwl_fitpart = w_center[rule2] 
    # # coefs = np.polyfit(logwl_fitpart,xhisq_fitpart, 2)
    # # w_new = np.arange(logLam1[0],logLam1[-1],(logLam1[-1]-logLam1[0])/len(logLam1))
    # # ffit = np.polyval(coefs, w_new)


    # #Make fit
    # coefs = np.polyfit(w_center_n,xhisquare_n, 6)
    # #coefs = np.

    # #print 'hurra: %s,%s,%s' % ((logLam1[-1]-logLam1[0]),len(logLam1),np.round((logLam1[-1]-logLam1[0])/len(logLam1),4))

    # w_new = np.arange(logLam1[0],logLam1[-1],(logLam1[-1]-logLam1[0])/len(logLam1))
    # ffit = np.polyval(coefs, w_new)

    # id_fit_0 = np.where(ffit<1)
    # ffit[id_fit_0] = 1

    # #print len(noise), len(ffit)

    # if len(noise) == len(ffit):
    #     noise = noise/np.sqrt(ffit)
    # else:
    #     if len(noise) > len(ffit):
    #         noise = noise[:len(ffit)]
    #         noise = noise/np.sqrt(ffit)
    #         #print 'hu'
    #     elif len(noise) < len(ffit):
    #         #print len(ffit)
    #         ffit = ffit[:len(noise)]
    #         w_center = w_center[:len(noise)]

    #         #print len(ffit)

    #         noise = noise/np.sqrt(ffit)
    #         #print 'rra'
    #     else:
    #         print 'hey man'

    # #Errs = Errs/ffit**2
    # #Errs = Errs/np.sqrt(2)

    # # #
    # # print np.sqrt(2), np.sqrt(ffit[0]), np.sqrt(ffit[-1])
    # # #
    # # #
    # # ################################################################
    # #
    # plt.hlines(np.ones(2),-1000,8000,linestyle='solid')
    # plt.plot(w_center_n,xhisquare_n,'.',color='r')

    # #print len(w_center), len(ffit), len(xhisquare)

    # plt.plot(w_center,xhisquare/ffit[:-100],'o',color='purple')
    # #plt.plot(w_center[rule1],xhisquare[rule1],'x')
    # plt.plot(w_center,ffit[:-100],color='g',linewidth=2.0)



    # plt.axis([8,9,0,15])
    # plt.show()

    # #
    # ################################################################





























#### Before 01.01.2016
'''

# Parameters

#No_of_spectra = 1 # 1-10

Npert = 100

for j in range(10):
    No_of_spectra = j+1

    ##################################################################
    ### Read in none-binned data (slitloss and aperture corrected) ###
    root = '/Users/mstockmann/X-Shooter/P93/Slitloss_aper_corr_spec_1D/V2_opt/*.fits'
    #root = '/Users/mstockmann/X-Shooter/P93/Slitloss_aper_corr_spec_1D/V2_medrebin20_opt/*.fits'
    Spath = glob.glob(root)
    #
    wl, Flux, Errs = readfitsfiles(Spath,No_of_spectra)
    name_qso = Spath[No_of_spectra-1].split('/')[-1].split('_')[0]
    #
    # Setting negative elements to median of positive errors
    m = np.median(Errs[Errs > 0])
    Errs[Errs <= 0] = m
    #
    ##################################################################


    ############################
    ### Read in single model ###  ------------ When the correct binning is selected use the fixed z spec
    FAST_model_path = '/Users/mstockmann/PhD_DARK/Programs/FAST/FAST_v1.0/xshoot_phot_spec_Filt2/BEST_FITS/COSMOS2015_xshP93_V2_uJy.xshoot_med_ErrNorm_bin18_BC03V1_%s.fit' % (No_of_spectra)
    #
    bestfit_data = np.genfromtxt(FAST_model_path)
    fl_model = bestfit_data[:,1]
    wl_model = bestfit_data[:,0]
    #
    ############################




    #########################################
    ### Changing to rest-frame wavelength ###
    path_z = '/Users/mstockmann/PhD_DARK/Programs/ppxf/code/zspec.txt'
    data_z = np.genfromtxt(path_z)
    z_spec = data_z[:,-1][No_of_spectra-1]
    #
    wl = wl/(1+z_spec)
    wl_model = wl_model/(1+z_spec)
    #
    #########################################


    ##############################################
    ### Remove model points outside data range ###
    # ----- The resolution element changes for the initial  pixels, try to exclude and see how it affects the v-disp
    ID_min = np.searchsorted(wl_model,wl[0]) # alter ID so len(template) > len(galaxy), (ppxf requirement)
    ID_max = np.searchsorted(wl_model,wl[-1])
    wl_model = wl_model[ID_min:ID_max]
    fl_model = fl_model[ID_min:ID_max]
    #
    ##############################################



    #################################################################################
    ### Check that the resolution element is the same as a function of wavelength ###
    ### This have been tested for the data that appears to be constant, where as the fit have funny effect at low wl
    #k = wl_model[1:]-wl_model[:-1]
    #print 'min/max difference between spacing (should the same): ',np.min(k),'/',np.max(k)
    #plt.plot(wl_model[:-1],k)
    #plt.show()
    #
    #################################################################################


    v_disp_list = np.zeros(Npert+1)
    # Loop over the perturbed spectra to obtain individual velocity dispersion estimates
    for i in range(Npert):

        ####################################################################################
        ### Logrebin the spectra (Conserve the Flux, while changing into velocity space) ###
        #
        # We logrebin the spectra with the smallest resolution first
        # to avoid interpolating the spectra with least data points.
        #
        # Log Rebin fit
        lamRange_model = np.array([min(wl_model),max(wl_model)])
        model1, logLam_model, velscale_model = util.log_rebin(lamRange_model, fl_model, velscale=None)
        #
        # Log rebin for galaxy
        lamRange_gal = np.array([min(wl),max(wl)])
        galaxy, logLam1, velscale = util.log_rebin(lamRange_model, Flux, velscale=velscale_model)
        #
        ##### correct######################
        # Log Rebin Error spectrum
        noise, logLam1_err, velscale_err = util.log_rebin(lamRange_model, Errs, velscale=velscale_model)
        #
        #print len(model1), len(galaxy), len(noise)
        #################################
        #plt.plot(logLam1,galaxy)        #
        #plt.plot(logLam_model,model1)   #
        #plt.axis([8,9,-20,20])          #
        #plt.show()                      #
        #################################
        ####################################################################################



        ###########################################################
        ### Perturb the model using the noise from the data WGN ###
        #
        modelp = np.zeros(len(noise))
        for n in range(len(noise)):
            mu, sigma = model1[n], noise[n]
            modelp[n] = np.random.normal(mu, sigma,1)
        #
        #########################
        #plt.plot(model1,'b')    #
        #plt.plot(modelp,'r')    #
        #plt.axis([8,9,-20,20])  #
        #plt.show()              #
        #########################
        #
        ###########################################################




        ################################################################################
        ### Convolving the highest resolution spectra to match the lowest resolution ###
        #
        # I here assume that the Dv = c/R = const in vel space, which I'm not sure is actually the case.
        #
        #
        FWHM_model = c/2000 # R = 2000. See intro in www.mpe.mpg.de/events/alvio05/contributions/bruzual@alvio05.pdf
        FWHM_gal = c/5300 # R = 5300 (X-shooter manual p.5)
        #
        #
        # Find the difference (Why squared?)
        FWHM_dif = np.sqrt(FWHM_model**2 - FWHM_gal**2)
        #
        # Converting from FWHM = sigma*sqrt(8*ln(2)). Also dividing by velscale change the units from km/s -> pix
        sigma_filt = FWHM_dif/2.355/velscale_model
        #
        #
        # Smoothing/convolving the (higher resolution) spectra, in this case the galaxy
        # Flux
        galaxy = ndimage.gaussian_filter1d(galaxy,sigma_filt)
        # Errs
        noise = ndimage.gaussian_filter1d(noise,sigma_filt)
        #
        ###
        ################################################################################



        ########################################
        ### Selecting the good pixel regions ###
        #
        Range_array = [13500,14240,17919,19641,20590]
        R_rest = Range_array/(1+z_spec)
        R_logrest = np.log(R_rest)
        #
        mask_goodpix = (logLam1 < R_logrest[0]) | (R_logrest[1] < logLam1) & (logLam1 < R_logrest[2]) | (R_logrest[3] < logLam1) & (logLam1 < R_logrest[4])
        goodPixels = np.where(mask_goodpix == True)[0]
        baadPixels = np.where(mask_goodpix == False)[0]
        #
        #########################################



        ####################
        ### Running pPXF ###
        #
        # Guess parameters
        sigma = 400
        #
        #
        ## Set Parameters:
        templates = modelp # Input model template
        # galaxy:   The logrebinned, and convolved data
        # noise:    The logrebinned, and convolved noise
        velscale = velscale_model
        start = [0,sigma]
        # goodpixels: The pixels that is included in the fit
        # plot: True if the plot of the fit needs to be showed
        Sum_degree = 8 # Belli et al. 2014
        Multi_degree = 3 # Belli et al. 2014
        #
        #
        # Run ppxf module
        
        print len(templates), len(galaxy)
        
        #pp = ppxf.ppxf(templates,galaxy, noise, velscale, start, goodpixels=goodPixels, plot=True, degree=Sum_degree, mdegree=Multi_degree)
        
        pp = ppxf.ppxf(templates, galaxy, noise, velscale, start, goodpixels=goodPixels, plot=True, degree=4, mdegree=4)#, vsyst=sigma, clean=True, regul=0)
        #
        #
        
        z = (z_spec + 1)*(1 + pp.sol[0]/c) - 1
        if i == 0:
            v_disp_list[i] = z
        else:
            v_disp_list[i] = pp.sol[1]
        
        print 'Best-fitting redshift z: %s' % (z)
        #print 'Best-fitting error on redshift z:',((z_spec + 1)*(1 + pp.sol/c) - 1) - ((z_spec + 1)*(1 + pp.error*np.sqrt(pp.chi2)/c) - 1)
        #
        #obj_spec = spec.resamplespec(wl,np.e**logLam1,pp.galaxy, oversamp=1000)
        #err_spec = spec.resamplespec(wl,np.e**logLam1_err,pp.noise, oversamp=1000)
        # logLam1_err
        #
        #template_fit = spec.resamplespec(wl,np.e**logLam1,pp.bestfit, oversamp =1000)
        #




    binned = root.split('/')[-2]
    outpath = 'out/Pmodel_dregrees4_mdegrees4/vdisp_Pmodel_N%s_%s_%s.txt' % (Npert,binned,name_qso)
    np.savetxt(outpath,v_disp_list)


'''







### Velocity dispersion PDF - model perturbation ###
'''

path_nobin = 'out/Pmodel_dregrees4_mdegrees4/vdisp_Pmodel_N100_*.txt'

All_path = glob.glob(path_nobin)
#print All_path




plt.figure(figsize=(14,8))
for i in range(len(All_path)):
    ipath = All_path[i]
    Npert = ipath.split('/')[-1].split('_')[2].split('N')[-1]

    data = np.genfromtxt(ipath)

    plt.subplot(2,5,i+1)
    if i == 0:
        plt.annotate('N perturbations: %s' % (Npert),xy=(80,300), bbox=dict(boxstyle='square', fc='red', alpha=0.5))
    plt.hist(data,bins=20,label='%s' % (ipath.split('/')[-1].split('_')[-1].split('.')[0]))
    plt.axis([0,1000,0,400])
    plt.legend(loc='upper left')


plt.xlabel('$\sigma_{vel}$', position=(1,0),fontsize=25)
plt.show()

'''





### Plot the velocity dispersions - data pert ###
'''
path_bin = './out/Pdata_dregrees4_mdegrees4/vdisp_Pdata_V2_medrebin20*.txt'
path_nobin = './out/Pdata_dregrees4_mdegrees4/vdisp_Pdata_V2_opt*.txt'

All_path = glob.glob(path_nobin)

plt.figure(figsize=(14,8))
for i in range(len(All_path)):

    ipath = All_path[i]

    data = np.genfromtxt(ipath)

    plt.subplot(2,5,i+1)
    if i == 0:
        plt.annotate('N perturbations: 1000',xy=(80,300), bbox=dict(boxstyle='square', fc='red', alpha=0.5))
    plt.hist(data,bins=20,label='%s' % (ipath.split('/')[-1].split('_')[-1].split('.')[0]))
    plt.axis([0,1000,0,400])
    plt.legend(loc='upper left')


plt.xlabel('$\sigma_{vel}$', position=(1,0),fontsize=25)
plt.show()

'''



'''
# ###################################################
# ### Calculate the chi2 correction to the errors ###
# #
# binsize = 100
# len_loop = len(Flux[:-binsize])
# w_center = np.zeros(len_loop)
# xhisquare = np.zeros(len_loop)
# for i in range(len_loop):
#     w1, w2 = wl[i], wl[i+binsize]
#     F_i = Flux[i:i+binsize]
#     E_i = Errs[i:i+binsize]

#     ID1 = np.searchsorted(wl_model,w1)
#     ID2 = np.searchsorted(wl_model,w2)
#     if not ID1:
#         Model_i = F_i+E_i
#     else:
#         Model_i = fl_model[ID1:ID2]
#         Model_i = np.interp(wl[i:i+binsize],wl_model[ID1:ID2], Model_i)

#     if len(F_i) == len(Model_i) == len(E_i):
#         xhisquare[i] = np.sum((F_i-Model_i)**2/E_i**2)/binsize
#         w_center[i] = w1+(w2-w1)/2

# xhi_med = np.median(xhisquare)
# rule1 = np.where(xhisquare > 5*xhi_med)
# xhisquare[rule1] = 1
# #
# # Make fit
# coefs = np.polyfit(w_center,xhisquare, 2)
# w_new = np.arange(wl[0],wl[-1],(wl[-1]-wl[0])/len(wl))
# ffit = np.polyval(coefs, w_new)
# #ffit = ffit+(1-ffit[0]) # Scale so that f[0] is chi2/dof
# #
# # Calibrate errors to a average chi2~1
# #Errs = Errs/np.sqrt(ffit)
# #Errs = Errs/ffit**2
# Errs = Errs/np.sqrt(2)

# #
# print np.sqrt(2), np.sqrt(ffit[0]), np.sqrt(ffit[-1])
# #
# #
# ################################################################

# print np.median(xhisquare)
# plt.hlines(np.median(xhisquare),-1000,8000,linestyle='solid')
# plt.plot(w_center,xhisquare,'.',color='r')
# plt.plot(w_center[rule1],xhisquare[rule1],'x')
# plt.plot(w_new,ffit,color='g',linewidth=2.0)

# print len(ffit), len(Errs)
# plt.axis([3000,7500,0,5])
# plt.show()

# ###############################################################

'''









