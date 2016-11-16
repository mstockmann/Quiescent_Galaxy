;#############################################################################
;
; Usage example for the procedure PPXF, which
; implements the Penalized Pixel-Fitting (pPXF) method by
; Cappellari M., & Emsellem E., 2004, PASP, 116, 138.
; The example also shows how to include a library of templates
; and how to mask gas emission lines if present.
; The example is specialized for a fit to a SDSS spectrum.
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
;   V2.3: Specialized for SDSS spectrum following requests from users.
;       Renamed PPXF_KINEMATICS_EXAMPLE_SDSS. MC, Oxford, 12 January 2012
;       
;#############################################################################
pro ppxf_kinematics_example_sdss

; Read SDSS DR8 galaxy spectrum taken from here http://www.sdss3.org/dr8/
; The spectrum is *already* log rebinned by the SDSS DR8
; pipeline and log_rebin should not be used in this case.
;
t = mrdfits('spectra/NGC3522_SDSS.fits',1,/SILENT)

; Only use the wavelength range in common between galaxy and stellar library.
; Of course for SDSS spectra one should use a more up-to-date library than 
; the one used in this pPXF example, with a more extended wavelength coverage.
; In that way the entire SDSS wavelength region can be fitted.
; See http://purl.org/cappellari/idl for suggestions of stellar libraries.
; 
t = t[where(t.wavelength gt 4795 and t.wavelength lt 5465)]
galaxy = t.flux/median(t.flux)  ; Normalize spectrum to avoid numerical issues
wave = t.wavelength
noise = galaxy*0 + 1           ; Assume constant noise per pixel here

; The velocity step was already chosen by the SDSS pipeline
; and we convert it below to km/s
;
c = 299792.458d ; speed of light in km/s
velScale = alog(wave[1]/wave[0])*c
FWHM_gal = 2.76 ; SDSS has an instrumental resolution FWHM of 2.76A.

; If the galaxy is at a significant redshift (z > 0.03), one would need to apply 
; a large velocity shift in PPXF to match the template to the galaxy spectrum.
; This would require a large initial value for the velocity (V > 1e4 km/s) 
; in the input parameter START = [V,sig]. This can cause PPXF to stop! 
; The solution consists of bringing the galaxy spectrum roughly to the 
; rest-frame wavelength, before calling PPXF. In practice there is no 
; need to modify the spectrum in any way, given that a red shift 
; corresponds to a linear shift of the log-rebinned spectrum. 
; One just needs to compute the wavelength range in the rest-frame
; and adjust the instrumental resolution of the galaxy observations.
; This is done with the following three commented lines:
;
; z = 1.23 ; Initial estimate of the galaxy redshift
; wave = wave/(1+z) ; Compute approximate restframe wavelength
; FWHM_gal = FWHM_gal/(1+z)   ; Adjust resolution in Angstrom

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
; with the quadratic difference between the SDSS and the 
; Vazdekis instrumental resolution. Logarithmically rebin 
; and store each template as a column in the array TEMPLATES.

; Quadratic sigma difference in pixels Vazdekis --> SDSS
; The formula below is rigorously valid if the shapes of the 
; instrumental spectral profiles are well approximated by Gaussians. 
;
FWHM_dif = SQRT(FWHM_gal^2 - FWHM_tem^2)
sigma = FWHM_dif/2.355/sxpar(h2,'CDELT1') ; Sigma difference in pixels 

; IMPORTANT: To avoid spurious velocity offsets of the templates, the
; NPIXEL keyword in PSF_GAUSSIAN must be an odd integer as done below.
;
lsf = psf_Gaussian(NPIXEL=2*ceil(4*sigma)+1, ST_DEV=sigma, /NORM, NDIM=1)
for j=0,nfiles-1 do begin
    fits_read, vazdekis[j], ssp
    ssp = convol(ssp,lsf) ; Degrade template to SDSS resolution
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
dv = (logLam2[0] - alog(wave[0]))*c ; km/s

vel = 1305d ; Initial estimate of the galaxy velocity in km/s
goodPixels = ppxf_determine_goodPixels(alog(wave),lamRange2,vel)

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

;###################################################
; Expected Output from PPXF_KINEMATICS_EXAMPLE_SDSS: 
; (after 0.5s on a 2.5GHz Core i5) 
;###################################################
;
;     feval         V     sigma        h3        h4        h5        h6  Chi2/DOF
;        32    1305.1      97.9    -0.023    -0.002     0.000     0.000     0.000
;Nonzero Templates: 4 / 19
;Templates weights:
;      0.000      0.000      0.000      0.000      0.000      0.000      0.000      0.000      0.000      0.282
;      0.113      0.172      0.000      0.000      0.000      0.000      0.512      0.000      0.000
;Formal errors:    dV    dsigma       dh3       dh4
;                 2.6       3.9     0.022     0.028
         
end
;------------------------------------------------------------------------------
