;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Anna Gallazzi, 22.11.2010
;; function fiber_index changed to allow a larger wavelength window
;;  accounting for the wavelength binning of the data (different from SDSS)
;;  26.11.2010: Added keyword DLAMBDA to fiber_index
;;  	This enters in the computation of the error on the bandpass fluxes
;;  DLAMBDA should be passed as argument to function d4000
;;
;; Anna Gallazzi, 23.06.2011
;; fiber_index: added check on number of good pixels in each of the index
;;  	bandpasses. Index is measured (in lineindex) only if there are at least
;;  	80% of pixels in the pseudo-continuum bands and 90% in the central band
;; d4000: added check on number of good pixels in the two bandpasses (require
;;  	at least 80%); need to pass rest-frame bin width (dlambda)
;;
;;  07.08.2012:
;;  Added keyword INTERP_IDX in order to use the bestfit continuum+smooth residuals
;;  when the SIDE band(s) are affected by CCD gaps.
;;  It can be passed from fiber_lfit.pro, i.e. only for the emission line corrected spectrum
;;  !!! If set, mod_interp and lambda_interp and okpix need to be
;;  passed !!!!
;;
;;  Stefano Zibetti, 01.08.2013: fiber_index --> spaxel_index
;;  Minor changes and a bit of clean-up to use with CALIFA
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


pro read_idxinfo
    common idxinfo,indexpars
    ;------------------------------------------------------------------------------
    ; Read info about indices - store in common block
    ;------------------------------------------------------------------------------
    indexlist=getenv('HOME')+'/idl/idl_jarle/IDL/platefit/etc/indexlist.txt'
    readcol, indexlist, index, wli, wlf, bc_wli, bc_wlf, rc_wli, rc_wlf, units, $
        format = 'A,F,F,F,F,F,F,A', /silent
        
    indexinfo = {iinfo, index: ' ', bandpass: [0.0, 0.0], $
        blue_continuum: [0.0, 0.0], red_continuum: [0.0, 0.0], $
        units: ' '}
        
    nindices = n_elements(index)
    indexpars = replicate(indexinfo, nindices)
    
    indexpars.index = index
    indexpars.bandpass = transpose([[wli], [wlf]])
    indexpars.blue_continuum = transpose([[bc_wli], [bc_wlf]])
    indexpars.red_continuum = transpose([[rc_wli], [rc_wlf]])
    indexpars.units = units
end



;*******************************************************************************
function f_d4000, wl, flux, err, narrow = narrow, redside=fl_r, blueside=fl_b,$
        okpix=okpix, mod_interp=mod_interp, minfrac=minfrac
    ;*******************************************************************************
    ;; mod_interp is the interpolating model, must be sampled on wl, exactly like flux. If not given interpolation is made linearly
    ;; minfrac=minimum fraction of good pixels in either bands, defaults
    ;; to 2./3.
        
    ;; return [index,error]
    d4000 = [0., -1.]             ; default
    
    if (~(keyword_set(minfrac))) then minfrac=2./3.
    
    ;; check array consistency
    npix=n_elements(wl)
    if ((n_elements(flux) ne npix) or (n_elements(err) ne npix)) then begin
        message,'d4000: array dimensions do not match. Returning.',/continue
        return,d4000
    endif
    
    if (keyword_set(mod_interp)) then if (n_elements(mod_interp) ne npix) then begin
        message,'d4000: mod_interp array dimension does not match the spectrum dimension. Returning.',/continue
        return,d4000
    endif
    if (not keyword_set(okpix)) then okpix=findgen(n_elements(wl))
    
    if keyword_set(narrow) then begin
        ;; Balogh et al., 1999, ApJ, 527, 54
        wl_b1 = 3850.0  &  wl_b2 = 3950.0
        wl_r1 = 4000.0  &  wl_r2 = 4100.0
    endif else begin
        ;; Bruzual, 1983, ApJ, 273, 105
        wl_b1 = 3750.0  &  wl_b2 = 3950.0
        wl_r1 = 4050.0  &  wl_r2 = 4250.0
    endelse
    
    ;;;; Check that there are enough pixels in each bandpass and the
    ;;;; spectral coverage is fine
    mn = min(wl[okpix], max=mx)
    if (mn gt wl_b1 or mx lt wl_r2) then begin
        message,'Insufficient spectral covarage to compute D4000, returning.',/continue
        return,d4000
    endif
    
    ok_blue=where(wl[okpix] ge wl_b1 and wl[okpix] le wl_b2,nok_blue)
    blue=where(wl ge wl_b1 and wl le wl_b2,n_blue)
    frac_blue=1.*nok_blue/n_blue
    ok_red=where(wl[okpix] ge wl_r1 and wl[okpix] le wl_r2,nok_red)
    red=where(wl ge wl_r1 and wl le wl_r2,n_red)
    frac_red=1.*nok_red/n_red
    ;; print,frac_blue,frac_red
    if (frac_blue lt minfrac or frac_red lt minfrac) then begin
        message,'Less than minfrac of good pixels in at least one band. Cannot measure D4000, returning.',/continue
        return,d4000
    endif
    fl_r = 0.0 & fl_b=0.0
    if (~(keyword_set(mod_interp))) then begin  ;; simply linearly interpolate over bad pixels
        fl_r = integral(wl[okpix], flux[okpix]*wl[okpix]*wl[okpix], wl_r1, wl_r2)
        fl_b = integral(wl[okpix], flux[okpix]*wl[okpix]*wl[okpix], wl_b1, wl_b2)
        
        fl_r_err =  sqrt(quad_integral(wl[okpix], (err[okpix]*wl[okpix]*wl[okpix]),wl_r1, wl_r2))
        fl_b_err =  sqrt(quad_integral(wl[okpix], (err[okpix]*wl[okpix]*wl[okpix]),wl_b1, wl_b2))
        
        d4000[0] = (wl_b2 - wl_b1) / (wl_r2 - wl_r1)  *  fl_r / fl_b
        d4000[1] = (wl_b2 - wl_b1) / (wl_r2 - wl_r1) * $
            sqrt((fl_r_err/fl_b)^2 + (fl_b_err * fl_r/fl_b^2)^2)
    endif else begin ;; replace bad pixels with the model
        flux_new=mod_interp
        flux_new[okpix]=flux[okpix]
        err_new=mod_interp/mean(flux[okpix]/err[okpix])
        err_new[okpix]=err[okpix]
        
        fl_r = integral(wl, flux_new*wl*wl,  wl_r1, wl_r2)
        fl_b = integral(wl, flux_new*wl*wl,  wl_b1, wl_b2)
        
        fl_r_err =  sqrt(quad_integral(wl, (err_new*wl*wl), wl_r1, wl_r2))
        fl_b_err =  sqrt(quad_integral(wl, (err_new*wl*wl), wl_b1, wl_b2))
        
        d4000[0] = (wl_b2 - wl_b1) / (wl_r2 - wl_r1)  *  fl_r / fl_b
        d4000[1] = (wl_b2 - wl_b1) / (wl_r2 - wl_r1) * $
            sqrt((fl_r_err/fl_b)^2 + (fl_b_err * fl_r/fl_b^2)^2)
    endelse
    return,d4000
end


;*******************************************************************************
function spaxel_index, restwl, flux, err, dlambda_samp = dlambda_samp, $
        okpix=okpix, mod_interp=mod_interp, $
        index_idx=index_idx, minfrac_side=minfrac_side, minfrac_cen=minfrac_cen
    ;*******************************************************************************
    ;; dlambda_samp is the spectral sampling in wavelength, if not given take the average pixel size
    ;; index_idx is an optional array of idx corresponding to the indices to use from indexpars
    ;; mod_interp is the interpolating model used to replace bad pixels.
    ;;            must be sampled on restwl, exactly like flux.
    ;;            If not given interpolation is made linearly during the integration
    ;; minfrac_[side|cen] is the minimum fraction of good pixels in any of the bands, defaults
    ;; to 0.8 for side and 0.9 for cen.
    ;; If the fraction of good pixels in any of the bands is less than minfrac the index cannot be computed
        
    common idxinfo
    
    if (~(keyword_set(minfrac_side))) then minfrac_side=0.8
    if (~(keyword_set(minfrac_cen ))) then minfrac_cen =0.9
    
    if (not keyword_set(index_idx)) then index_idx=indgen(n_elements(indexpars))
    
    n_indexes=n_elements(indexpars[index_idx])
    
    nwl=n_elements(restwl)
    
    if (not keyword_set(okpix)) then okpix=findgen(nwl)
    
    minwl=min(restwl,max=maxwl)
    if (not keyword_set(dlambda_samp)) then dlambda_samp=(maxwl-minwl)/nwl
    
    ;; create output structure
    idx_measure=create_struct('IDX_NAME','','IDX_VALUE',0.,'IDX_ERR',0.)
    output=replicate(idx_measure,n_indexes)
    
    ;*******************************************************************************
    ; Measure line indicies
    ;*******************************************************************************
    
    pixflag=intarr(nwl)+1
    pixflag[okpix]=0
    
    for i = 0, n_indexes - 1 do begin
        locali_all = where(restwl gt indexpars[index_idx[i]].blue_continuum[0] - 3.*dlambda_samp and $
            restwl lt indexpars[index_idx[i]].red_continuum[1] + 3.*dlambda_samp)
        restwl_tmp=restwl[locali_all]
        flux_tmp  =  flux[locali_all]
        err_tmp   =   err[locali_all]
        okpix_tmp =where(pixflag[locali_all] eq 0)
        ;if (strupcase(indexpars[index_idx[i]].index) eq 'LICK_HG_A') then stop
        if locali_all[0] eq -1 then begin
            message,indexpars[index_idx[i]].index+' cannot be measured, not covered',/continue
            index = [0,-1L]
            flag_blue=1 & flag_red=1 & flag_cen=1
        endif else begin
            ;;;;;;;; Check if there are enough pixels defined in each bandpass
            ;;;;;;;; Require to have at least minfrac of expected pixels
            ; Blue continuum bandpass
            ok_blue=where(restwl_tmp[okpix_tmp] ge indexpars[index_idx[i]].blue_continuum[0] and $
                restwl_tmp[okpix_tmp] le indexpars[index_idx[i]].blue_continuum[1],nblue)
            all_blue=where(restwl_tmp ge indexpars[index_idx[i]].blue_continuum[0] and $
                restwl_tmp le indexpars[index_idx[i]].blue_continuum[1],nblue_all)
            if (float(nblue/nblue_all) lt minfrac_side) then begin
                flag_blue=1
                message,indexpars[index_idx[i]].index+': too many gaps in the blue sideband!',/continue
            endif else flag_blue=0
            ; Red continuum bandpass
            ok_red=where(restwl_tmp[okpix_tmp] ge indexpars[index_idx[i]].red_continuum[0] and $
                restwl_tmp[okpix_tmp] le indexpars[index_idx[i]].red_continuum[1],nred)
            all_red=where(restwl_tmp ge indexpars[index_idx[i]].red_continuum[0] and $
                restwl_tmp le indexpars[index_idx[i]].red_continuum[1],nred_all)
            if (float(nred/nred_all) lt minfrac_side) then begin
                flag_red=1
                message,indexpars[index_idx[i]].index+': too many gaps in the red sideband!',/continue
            endif else flag_red=0
            ; Central bandpas
            ok_cen=where(restwl_tmp[okpix_tmp] ge indexpars[index_idx[i]].bandpass[0] and $
                restwl_tmp[okpix_tmp] le indexpars[index_idx[i]].bandpass[1],ncen)
            all_cen=where(restwl_tmp ge indexpars[index_idx[i]].bandpass[0] and $
                restwl_tmp le indexpars[index_idx[i]].bandpass[1],ncen_all)
            if (float(ncen/ncen_all) lt minfrac_cen) then begin
                flag_cen=1
                message,indexpars[index_idx[i]].index+': too many gaps in the central passband!',/continue
            endif else flag_cen=0
            ; Now compute index only if there are enough pixels in each bandpass
            if (flag_blue or flag_red or flag_cen) then begin
                message,indexpars[index_idx[i]].index+' cannot be measured',/continue
                index = [0,-1L]
            endif else begin
                if (keyword_set(mod_interp)) then begin ;; interpolate using the model replacement and then use the bad pix as if they were good
                    restwl_new=restwl_tmp
                    flux_new=mod_interp[locali_all]
                    flux_new[okpix_tmp]=flux_tmp[okpix_tmp]
                    err_new=mod_interp[locali_all]/mean(flux_tmp[okpix_tmp]/err_tmp[okpix_tmp]) ;; same average SNR over the index region
                    err_new[okpix_tmp]=err_tmp[okpix_tmp]
                    
                endif else begin ;; otherwise just forget about the bad pixels, the integrals will be done interpolating linearly over them
                    restwl_new=restwl_tmp[okpix_tmp]
                    flux_new=flux_tmp[okpix_tmp]
                    err_new=err_tmp[okpix_tmp]
                endelse
                ;;stop
                index = lineindex(restwl_new, flux_new, err_new, indexpars[index_idx[i]])
                
            endelse
        endelse
        
        ;;*****************************************************************************
        ;; Store indices in appropriate sturctures
        ;;*****************************************************************************
        output[i].idx_name=strupcase(indexpars[index_idx[i]].index)
        output[i].idx_value=index[0]
        output[i].idx_err=index[1]
    endfor
    
    ;      stop
    return,output
;;*******************************************************************************
;; Measure 4000 A break (integrate in F-nu units!): TBD
;;*******************************************************************************
    
end

