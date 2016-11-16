;#############################################################################
;
; Usage example for the procedure PPXF, which
; implements the Penalized Pixel-Fitting (pPXF) method by
; Cappellari M., & Emsellem E., 2004, PASP, 116, 138.
; The example also shows how to include a library of templates
; and how to mask gas emission lines if present.
; 
; MODIFICATION HISTORY:
;   V1.0: Written by Michele Cappellari, Leiden 11 November 2003
;   V1.1: Log rebin the galaxy spectrum. Show how to correct the velocity
;       for the difference in starting wavelength of galaxy and templates.
;       MC, Vicenza, 28 December 2004
;   V1.11: Included explanation of correction for instrumental resolution.
;       After feedback from David Valls-Gabaud. MC, Venezia, 27 June 2005
;   V2.0: Included example routine to determine the goodPixels vector
;       by masking known gas emission lines. MC, Oxford, 30 October 2008
;   V2.01: Included instructions for high-redshift usage. Thanks to Paul Westoby
;       for useful feedback on this issue. MC, Oxford, 27 November 2008
;   V2.02: Included example for obtaining the best-fitting redshift.
;       MC, Oxford, 14 April 2009
;   V2.1: Bug fix: Force PSF_GAUSSIAN to produce a Gaussian with an odd 
;       number of elements centered on the middle one. Many thanks to 
;       Harald Kuntschner, Eric Emsellem, Anne-Marie Weijmans and 
;       Richard McDermid for reporting problems with small offsets 
;       in systemic velocity. MC, Oxford, 15 February 2010
;   V2.11: Added normalization of galaxy spectrum to avoid numerical
;       instabilities. After feedback from Andrea Cardullo.
;       MC, Oxford, 17 March 2010
;   V2.2: Perform templates convolution in linear wavelength. 
;       This is useful for spectra with large wavelength range. 
;       MC, Oxford, 25 March 2010
;   V2.21: Updated for Coyote Graphics. MC, Oxford, 11 October 2011
;   V2.22: Renamed PPXF_KINEMATICS_EXAMPLE_SAURON to avoid conflict with the 
;       new PPXF_KINEMATICS_EXAMPLE_SDSS. Removed DETERMINE_GOOPIXELS which was 
;       made a separate routine. MC, Oxford, 12 January 2012
;       
;#############################################################################
pro ppxf_kinematics_example_sauron

; Read a galaxy spectrum and define the wavelength range
;
fits_read, 'spectra/NGC4550_SAURON.fits', gal_lin, h1
lamRange1 = sxpar(h1,'CRVAL1') + [0d,sxpar(h1,'CDELT1')*(sxpar(h1,'NAXIS1')-1d)]
FWHM_gal = 4.2 ; SAURON has an instrumental resolution FWHM of 4.2A.

; If the galaxy is at a significant redshift (z > 0.03), one would need to apply 
; a large velocity shift in PPXF to match the template to the galaxy spectrum.
; This would require a large initial value for the velocity (V > 1e4 km/s) 
; in the input parameter START = [V,sig]. This can cause PPXF to stop! 
; The solution consists of bringing the galaxy spectrum roughly to the 
; rest-frame wavelength, before calling PPXF. In practice there is no 
; need to modify the spectrum before the usual LOG_REBIN, given that a 
; red shift corresponds to a linear shift of the log-rebinned spectrum. 
; One just needs to compute the wavelength range in the rest-frame
; and adjust the instrumental resolution of the galaxy observations.
; This is done with the following three commented lines:
;
; z = 1.23 ; Initial estimate of the galaxy redshift
; lamRange1 = lamRange1/(1+z) ; Compute approximate restframe wavelength range
; FWHM_gal = FWHM_gal/(1+z)   ; Adjust resolution in Angstrom

log_rebin, lamRange1, gal_lin, galaxy, logLam1, VELSCALE=velScale
galaxy = galaxy/median(galaxy) ; Normalize spectrum to avoid numerical issues
noise = galaxy*0 + 1           ; Assume constant noise per pixel here

; Read the list of filenames from the Single Stellar Population library
; by Vazdekis (1999, ApJ, 513, 224). A subset of the library is included 
; for this example with permission. See http://purl.org/cappellari/idl 
; for suggestions of more up-to-date stellar libraries.
;
vazdekis = file_search('spectra/Rbi1.30z*.fits',COUNT=nfiles)
FWHM_tem = 1.8 ; Vazdekis spectra have a resolution FWHM of 1.8A.

; Extract the wavelength range and logarithmically rebin one spectrum
; to the same velocity scale of the SAURON galaxy spectrum, to determine
; the size needed for the array which will contain the template spectra.
;
fits_read, vazdekis[0], ssp, h2
lamRange2 = sxpar(h2,'CRVAL1') + [0d,sxpar(h2,'CDELT1')*(sxpar(h2,'NAXIS1')-1d)]
log_rebin, lamRange2, ssp, sspNew, logLam2, VELSCALE=velScale
templates = dblarr(n_elements(sspNew),nfiles)

; Convolve the whole Vazdekis library of spectral templates 
; with the quadratic difference between the SAURON and the 
; Vazdekis instrumental resolution. Logarithmically rebin 
; and store each template as a column in the array TEMPLATES.

; Quadratic sigma difference in pixels Vazdekis --> SAURON
; The formula below is rigorously valid if the shapes of the 
; instrumental spectral profiles are well approximated by Gaussians. 
;
FWHM_dif = SQRT(FWHM_gal^2 - FWHM_tem^2)
sigma = FWHM_dif/2.355/sxpar(h2,'CDELT1') ; Sigma difference in pixels 

; IMPORTANT: To avoid spurious velocity offsets of the templates, the
; NPIXEL keyword in PSF_GAUSSIAN must be an odd integer as done below
;
lsf = psf_Gaussian(NPIXEL=2*ceil(4*sigma)+1, ST_DEV=sigma, /NORM, NDIM=1)
for j=0,nfiles-1 do begin
    fits_read, vazdekis[j], ssp
    ssp = convol(ssp,lsf) ; Degrade template to SAURON resolution
    ; From IDL 8.1 one can use the following commented line instead of the 
    ; above one, and the line with PSF_GAUSSIAN is not needed any more.  
    ; ssp = gauss_smooth(ssp,sigma)    
    log_rebin, lamRange2, ssp, sspNew, VELSCALE=velScale
    templates[*,j] = sspNew/median(sspNew) ; Normalizes templates 
endfor

; The galaxy and the template spectra do not have the same starting wavelength.
; For this reason an extra velocity shift DV has to be applied to the template
; to fit the galaxy spectrum. We remove this artificial shift by using the
; keyword VSYST in the call to PPXF below, so that all velocities are
; measured with respect to DV. This assume the redshift is negligible.
; In the case of a high-redshift galaxy one should de-redshift its 
; wavelength to the rest frame before using the line below (see above).
;
c = 299792.458d
dv = (logLam2[0]-logLam1[0])*c ; km/s

vel = 450d ; Initial estimate of the galaxy velocity in km/s
goodPixels = ppxf_determine_goodPixels(logLam1,lamRange2,vel)

; Uncomment the following line to set goodPixels by hand, instead 
; of using the optional determine_goodPixels routine above.
;
;goodPixels = [range(1,35),range(50,125),range(145,165),range(190,414)]

; Here the actual fit starts. The best fit is plotted on the screen.
; Gas emission lines are excluded from the pPXF fit using the GOODPIXELS keyword.
;
start = [vel, 180d] ; (km/s), starting guess for [V,sigma]
ppxf, templates, galaxy, noise, velScale, start, sol, $
    GOODPIXELS=goodPixels, /PLOT, MOMENTS=4, DEGREE=4, $
    VSYST=dv, ERROR=error    

print, 'Formal errors:    dV    dsigma       dh3       dh4'
print, error*sqrt(sol[6]), FORMAT='(10x,2f10.1,2f10.3)'

; If the galaxy is at significant redshift z and the wavelength has been
; de-redshifted with the three lines "z = 1.23..." near the beginning of 
; this procedure, the best-fitting redshift is now given by the following 
; commented line (equation 2 of Cappellari et al. 2009, ApJ, 704, L34):
;    
;print, 'Best-fitting redshift z:', (z + 1)*(1 + sol[0]/c) - 1

;##############################################
; Expected Output from PPXF_KINEMATICS_EXAMPLE: 
; (after 0.5s on a 2.5GHz Core i5) 
;##############################################
;
;     feval         V     sigma        h3        h4        h5        h6  Chi2/DOF
;        32     437.9     115.1     0.042     0.056     0.000     0.000     0.000
;Nonzero Templates: 5 / 19
;Templates weights:
;    0.00455      0.000     0.0919      0.163      0.150      0.000      0.000      0.000      0.000      0.000
;      0.000      0.000      0.000      0.000      0.000      0.000      0.604      0.000      0.000
;Formal errors:    dV    dsigma       dh3       dh4
;                 2.4       5.1     0.018     0.031
         
end
;------------------------------------------------------------------------------
