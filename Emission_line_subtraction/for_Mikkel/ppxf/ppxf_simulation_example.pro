;----------------------------------------------------------------------------
pro ppxf_simulation_example
; 
; Monte Carlo simulation of the kinematical extraction with pPXF. It is useful
; to determine the desired value for the BIAS keyword of the pPXF procedure.
; This procedure generates a plot similar (but not identical) to Figure 6 in
; Cappellari & Emsellem, 2004, PASP, 116, 138.
;
; A rough guideline to determine the BIAS value is the following: choose the *largest* 
; value which make sure that in the range sigma>3*velScale and for (S/N)>30 the true values 
; for the Gauss-Hermite parameters are well within the rms scatter of the measured values.
; See the documentation in the file ppxf.pro for a more accurate description.
;
; V1.0: By Michele Cappellari, Leiden, 28 March 2003
; V1.1: Included in the standard PPXF distribution. After feedback
;   from Alejandro Garcia Bedregal. MC, Leiden, 13 April 2005
; V1.11: Adjust GOODPIXELS according to the size of the convolution kernel.
;   MC, Oxford, 13 April 2010
; V1.12: Use Coyote Graphics (http://www.idlcoyote.com/) by David W. Fanning.  
;   The required routines are now included in NASA IDL Astronomy Library. 
;   MC, Oxford, 29 July 2011

fits_read, 'spectra/Rbi1.30z+0.00t12.59.fits', ssp, h ; Solar metallicitly, Age=12.59 Gyr
lamRange = sxpar(h,'CRVAL1') + [0d,sxpar(h,'CDELT1')*(sxpar(h,'NAXIS1')-1d)]
log_rebin, lamRange, ssp, star, VELSCALE=velScale

; The finite sampling of the observed spectrum is modeled in detail:
; the galaxy spectrum is obtained by oversampling the actual observed spectrum
; to a high resolution. This represent the true spectrum, which is later resampled
; to lower resolution to simulate the observations on the CCD. Similarly, the
; convolution with a well-sampled LOSVD are done on the high-resolution spectrum,
; and later resampled to the observed resolution before fitting with PPXF.

n = n_elements(star)
factor = 10                    ; Oversampling integer factor for an accurate convolution
starNew = rebin(star,n*factor) ; This is the underlying spectrum, known at high resolution
star = rebin(starNew,n)        ; Make sure that the observed spectrum is the integral over the pixels

vel = 0.3d      ; velocity in *pixels* [=V(km/s)/velScale]
h3 = 0.1d       ; Adopted G-H parameters of the LOSVD
h4 = -0.1d
sn = 60d        ; Adopted S/N of the Monte Carlo simulation
m = 300         ; Number of realizations of the simulation
sigmaV = range(0.8d,4d,m) ; Range of sigma in *pixels* [=sigma(km/s)/velScale]

result = fltarr(m,4) ; This will store the results
t = systime(1)
loadct, 12

for i=0,m-1 do begin

    sigma = sigmaV[i]
    dx = ceil(abs(vel)+4.0*sigma)   ; Sample the Gaussian and GH at least to vel+4*sigma
    x = range(dx,-dx,2*dx*factor+1) ; Evaluate the Gaussian using steps of 1/factor pixel.
    w = (x - vel)/sigma
    w2 = w^2
    gauss = exp(-0.5d*w2)/(sqrt(2d*!dpi)*sigma*factor) ; Normalized total(gauss)=1
    h3poly = w*(2d*w2 - 3d)/Sqrt(3d)           ; H3(y)
    h4poly = (w2*(4d*w2 - 12d) + 3d)/Sqrt(24d) ; H4(y)
    losvd = gauss *(1d + h3*h3poly + h4*h4poly)

    galaxy = convol(starNew,losvd,/EDGE_TRUNCATE) ; Convolve the oversampled spectrum
    galaxy = rebin(galaxy,n) ; Integrate spectrum back to original resolution
    noise = galaxy/sn        ; 1sigma error spectrum
    galaxy = galaxy + noise*randomn(s,n) ; Add noise to the galaxy spectrum
    start = [vel+randomu(s), sigma*(1.0+(randomu(s)-0.5)/3)]*velScale ; Convert to km/s

    ppxf, star, galaxy, noise, velScale, start, sol, $
        GOODPIXELS=range(dx,n-dx), MOMENTS=4, BIAS=0.5, /PLOT
    result[i,*] = sol[0:3]

endfor
print, 'Calculation time:', systime(1)-t, ' seconds'

!p.multi = [0,2,2]
cgplot, sigmaV*velScale, result[*,0]-vel*velScale, PSYM=1, /YNOZERO, $
    YRANGE=[-40,40], YTITLE='V - V!Bin!N (km/s)', XTITLE='!4r!3!Bin!N (km/s)'
cgplot, [0,400], [0,0], COLOR='red', THICK=2, /OVERPLOT
cgplot, velscale*3*[1,1], [-40,40], LINE=1, /OVERPLOT
cgplot, sigmaV*velScale, result[*,1]-sigmaV*velScale, PSYM=1, $
    YRANGE=[-40,40], YTITLE='!4r!3 - !4r!3!Bin!N (km/s)', XTITLE='!4r!3!Bin!N (km/s)'
cgplot, [0,400], [0,0], COLOR='red', THICK=2, /OVERPLOT
cgplot, velscale*3*[1,1], [-40,40], LINE=1, /OVERPLOT
cgplot, sigmaV*velScale, result[*,2], PSYM=1, /YNOZERO, $
    YRANGE=[-0.2,0.2]+h3, YTITLE='h!B3!N', XTITLE='!4r!3!Bin!N (km/s)'
cgplot, [0,400], [h3,h3], COLOR='red', THICK=2, /OVERPLOT
cgplot, velscale*3*[1,1], [-0.2,0.2]+h3, LINE=1, /OVERPLOT
cgplot, sigmaV*velScale, result[*,3], PSYM=1, /YNOZERO, $
    YRANGE=[-0.2,0.2]+h4, YTITLE='h!B4!N', XTITLE='!4r!3!Bin!N (km/s)'
cgplot, [0,400], [h4,h4], COLOR=200, THICK=2, /OVERPLOT
cgplot, velscale*3*[1,1], [-0.2,0.2]+h4, LINE=1, /OVERPLOT
!p.multi = 0

end
;----------------------------------------------------------------------------
