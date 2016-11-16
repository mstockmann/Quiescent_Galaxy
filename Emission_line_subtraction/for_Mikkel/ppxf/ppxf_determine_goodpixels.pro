;------------------------------------------------------------------------------
function ppxf_determine_goodPixels, logLam, lamRangeTemp, vel
;
; PPXF_DETERMINE_GOODPIXELS: Example routine to generate the vector of goodPixels 
;     to be used as input keyword for the routine PPXF. This is useful to mask 
;     gas emission lines or atmospheric absorptions. 
;     It can be trivially adapted to mask different lines.
; 
; INPUT PARAMETERS:
; - LOGLAM: Natural logarithm ALOG(wave) of the wavelength in Angstrom 
;     of each pixel of the log rebinned *galaxy* spectrum.
; - LAMRANGETEMP: Two elements vectors [lamMin2,lamMax2] with the minimum and
;     maximum wavelength in Angstrom in the stellar *template* used in PPXF.
; - VEL: Estimate of the galaxy velocity in km/s.
; 
; V1.0: Michele Cappellari, Leiden, 9 September 2005
; V1.01: Made a separate routine and included additional common emission lines. 
;   MC, Oxford 12 January 2012

;        [OII] Hgamma   Hbeta        [OIII]           [OI]        [NII]        [SII]  
lines = [3727, 4341, 4861.3d, 4958.9d, 5006.9d, 6300, 6363, 6548, 6583, 6716, 6731] 
dv = lines*0+800d ; width/2 of masked gas emission region in km/s
c = 299792.458d ; speed of light in km/s

flag = bytarr(n_elements(logLam))

for j=0,n_elements(lines)-1 do $
    flag or= logLam gt alog(lines[j]) + (vel - dv[j])/c $
         and logLam lt alog(lines[j]) + (vel + dv[j])/c

flag or= logLam lt alog(lamRangeTemp[0]) + (vel + 900d)/c ; Mask edges of
flag or= logLam gt alog(lamRangeTemp[1]) + (vel - 900d)/c ; stellar library

return, where(flag eq 0)
end
;------------------------------------------------------------------------------
