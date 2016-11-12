;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;	Anna Gallazzi, 26.11.2010:
;;; 	Added argument DLAMBDA, which needs to be passed when lineindex is called
;;; 	This enters in the computation of the error on the bandpass flux
;;; 	CHECK ERROR MEASUREMENT!!!!!!

;;; 	Anna Gallazzi, 24.06.2011:
;;; 	Normally index was discarded if max(wave)<red[1] or min(wave]>blue[0]
;;; 	This is a bit too restrictive (imagine there is a gap that only affects
;;; 	a minor portion of the red or blue band). I've introduced an extrapolation
;;; 	from min(wave) [max(wave)] to blue(0) [red(1)]
;;;
;;;     Stefano Zibetti, 1.8.2013: 
;;;     Minor optimization: min and max wavelength computed only once
;;;     and stored
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


;*******************************************************************************
function lineindex, wave, flux, err, indexinfo, w, f_blue_cont = f_blue_cont, $
                    err_blue_cont = err_blue_cont, f_red_cont = f_red_cont, $
                    err_red_cont = err_red_cont
;*******************************************************************************
; pass wave, flux and err only in the relevant range in order to speed
; up calculation!!!
;
;*******************************************************************************
; Check to make sure the the continuum windows are within wavelength range
;*******************************************************************************

  index = [0.0, -1.0]


  minwave=min(wave,max=maxwave)
  if maxwave lt indexinfo.red_continuum[0] or $
     minwave gt indexinfo.blue_continuum[1] then begin ;; either of the sidebands are completely out of range
     f_blue_cont = -9.9
     err_blue_cont = -9.9
     f_red_cont = -9.9
     err_red_cont = -9.9
     print,indexinfo.index, 'cannot be measured in lineindex'
     return, index
  end

  if (maxwave lt indexinfo.red_continuum[1]) then begin  ;; part of the red continuum is out of range
     okred=where(wave gt indexinfo.red_continuum[0],complement=outred)
     print,indexinfo.index,' Extrapolated red bandpass'
                                ;stop
                                ;print,wave[okred],flux[okred]
     ;; add the reddest wavelength point by extrapolation
     wave_new=[wave[okred],indexinfo.red_continuum[1]]
     flux_new_red=interpol(flux[okred],wave[okred],wave_new)
     flux=[flux[outred],flux_new_red]
     wave=[wave,indexinfo.red_continuum[1]]
     ;; the last point in the err spec is given by the mean in the band
     err=[err,mean(err[okred])]
                                ;print,wave_new,flux_new_red
                                ;stop
  endif
  if (minwave gt indexinfo.blue_continuum[0]) then begin
     okblue=where(wave lt indexinfo.blue_continuum[1],complement=outblue)
     print,indexinfo.index,' Extrapolated blue bandpass'
                                ;print,wave[okblue],flux[okblue]
     wave_new=[indexinfo.blue_continuum[0],wave[okblue]]
     flux_new_blue=interpol(flux[okblue],wave[okblue],wave_new)
     flux=[flux_new_blue,flux[outblue]]
     wave=[indexinfo.blue_continuum[0],wave]
     err=[mean(err[okblue]),err]
                                ;print,wave_new,flux_new_blue
                                ;stop
  endif


;if max(wave) lt indexinfo.red_continuum[1] or $
;   min(wave) gt indexinfo.blue_continuum[0] then begin
;   f_blue_cont = -9.9
;   err_blue_cont = -9.9
;   f_red_cont = -9.9
;   err_red_cont = -9.9
;   print,indexinfo.index, 'cannot be measured in lineindex'
;   return, index
;end
;*******************************************************************************
; Find average red & blue continuum fluxes
;*******************************************************************************

  wlb = indexinfo.blue_continuum
  f_blue_cont = integral(wave, flux, wlb[0], wlb[1]) / (wlb[1] - wlb[0])
  err_blue_cont = sqrt(quad_integral(wave, err, wlb[0], wlb[1]) / (wlb[1] - wlb[0]))
  wl_blue = (wlb[0] + wlb[1]) / 2.0

  wlr = indexinfo.red_continuum
  f_red_cont = integral(wave, flux, wlr[0], wlr[1]) / (wlr[1] - wlr[0])
  err_red_cont = sqrt(quad_integral(wave, err, wlr[0], wlr[1]) / (wlr[1] - wlr[0]))
  wl_red = (wlr[0] + wlr[1]) / 2.0


;*******************************************************************************
; Compute continuum for the index (different for BH than Lick & DTT)
;*******************************************************************************

  lick = strpos(indexinfo.index, 'Lick_') & lick = lick[0] ne -1
  dtt = strpos(indexinfo.index, 'DTT_') & dtt = dtt[0] ne -1
  bh = strpos(indexinfo.index, 'BH_') & bh = bh[0] ne -1

  if bh then begin
     f_cont = (f_red_cont[0] + f_blue_cont[0]) / 2.0
     err_cont = sqrt(err_red_cont[0]^2 + err_blue_cont[0]^2) / 2.0
;  print, 'ERR_CONT=', err_cont, err_red_cont, err_blue_cont, format='(A9,2X,3(E11.5,1X))'
;  stop
  endif else begin
     m = (f_red_cont[0] - f_blue_cont[0]) / (wl_red - wl_blue)
     m_err = sqrt(err_blue_cont[0]^2 + err_red_cont[0]^2) / (wl_red - wl_blue) 
     b = f_blue_cont[0]
     b_err = err_blue_cont[0] 

     f_cont = m * (wave - wl_blue) + b
     err_cont = sqrt((m_err * (wave - wl_blue))^2 + b_err^2)
  endelse

;*******************************************************************************
; Compute equivalent width 
;*******************************************************************************

  wli = indexinfo.bandpass
  if indexinfo.units eq 'A' then begin
     eqwidth = integral(wave, (1 - flux/f_cont), wli[0], wli[1]) 

     nerri =  (flux / f_cont) * sqrt((err / flux)^2 + (err_cont / f_cont)^2)
     eqwidth_err=sqrt(quad_integral(wave, nerri, wli[0], wli[1]))
     index = [eqwidth, eqwidth_err]
  endif

;*******************************************************************************
; Compute magnitudes
;*******************************************************************************

  if indexinfo.units eq 'mag' then begin

     nflux = integral(wave, flux/f_cont, wli[0], wli[1]) 
     nerri =  (flux / f_cont) * sqrt((err/ flux)^2 + (err_cont / f_cont)^2)
     nerr =  sqrt(quad_integral(wave, nerri, wli[0], wli[1]))
     mag = -2.5 * alog10(nflux / (wli[1] - wli[0])) 
     mag_err = 2.5 / alog(10) * nerr / nflux

     ;;------------------------------------------------------------------
     ;; Added 10/4/04 JB
     ;; Check to see that the magnitude is finite - for the emission line
     ;; subtracted spectrum this might not be the case... I don't know
     ;; what the best solution would be - but set to -99 seems ok...
     ;;------------------------------------------------------------------
     if (size(mag, /n_dim) gt 0) then begin
        bad = where(finite(mag) eq 0, n_bad)
        if (n_bad gt 0) then mag[bad] = -99.9
     endif else if (finite(mag) eq 0) then mag = -99.9
     
     index = [mag, mag_err]
;;  stop
  endif


;*******************************************************************************

  return, index

end
