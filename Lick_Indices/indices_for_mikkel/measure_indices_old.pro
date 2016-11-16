;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; 5/9/2015
;;
;; 14 cQG sample:
;; Measure indices on final reduction from Mikkel
;; spectrum needs correction for slitloss and match to phot
;; spectrum not corrected for emission lines (OIII)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro measure_indices_new_redux
    common idxinfo
    DIR='test_spectrum/'

    file=DIR+'108899_avwCombined_OBx5_sig5_V1_NIRcorr_wmrebin15_opt.dat'
    str={lambda:fltarr(1),flux:fltarr(1),err:fltarr(1),qual:intarr(1)}
    nlambda=file_lines(file)
    data=replicate(str,nlambda)
    openr,lun,file,/get_lun
    readf,lun,data
    free_lun,lun
;;; Put wavelength in air
    vactoair,data.lambda

    zgal=2.2292
    normalisation = 1.e18

    lambda_res1=data.lambda[0]/(1.+zgal)
    flux1=data.flux[0]*normalisation
    error1=data.err[0]*normalisation

    lambda_res=lambda_res1
    flux=flux1
    error=error1
    
    set_plot,'x'
    device,decompose=0
    loadct,13
    ;tmp=get_position_arr(0,nx=1,ny=2,xmin=0.1,xmax=0.9,ymin=0.1,ymax=0.9,ygap=0)
    pos=[0.1,0.1,0.9,0.9]
    plot,lambda_res,flux,xstyle=1,ystyle=1,position=pos,yrange=[0.,12.],xrange=[3700.,5500.]
    oplot,lambda_res,error,color=100
;    oplot,lambda_res1,flux1,color=200
;    oplot,lambda_res1,error1,color=160
    ;pos=get_position_arr(1)
    ;plot,lambda_res,flux/error,xstyle=1,ystyle=1,position=pos,/noerase,xrange=[3700.,5500.]
    plot,lambda_res,flux/error,xstyle=1,ystyle=1,/noerase,xrange=[3700.,5500.]
    ;stop

    measure_idx2,flux1,error,lambda_res,outfile
    
end


pro read_idxinfo
    common idxinfo,indexpars
;------------------------------------------------------------------------------
; Read info about indices - store in common block
;------------------------------------------------------------------------------
    indexlist='indexlist.txt'
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
;**************************************************
function Lick,iidx,indexpars,wl,flux,err
;**************************************************

    dlambda= (wl[n_elements(wl)-1]-wl[0])/float(n_elements(wl))
    print,'dlambda= ',dlambda
    wlb = indexpars[iidx].blue_continuum
    f_blue_cont = integral(wl, flux, wlb[0], wlb[1]) / (wlb[1] - wlb[0])
    err_blue_cont = sqrt(integral(wl, (err)^2, wlb[0], wlb[1]) / (wlb[1] - wlb[0]))
    wl_blue = (wlb[0] + wlb[1]) / 2.0

    wlr = indexpars[iidx].red_continuum
    f_red_cont = integral(wl, flux, wlr[0], wlr[1]) / (wlr[1] - wlr[0])
    err_red_cont = sqrt(integral(wl, (err)^2, wlr[0], wlr[1]) / (wlr[1] - wlr[0]))
    wl_red = (wlr[0] + wlr[1]) / 2.0

    m = (f_red_cont[0] - f_blue_cont[0]) / (wl_red - wl_blue)
    m_err = sqrt(err_blue_cont[0]^2 + err_red_cont[0]^2) / (wl_red - wl_blue)
    b = f_blue_cont[0]
    b_err = err_blue_cont[0] 

    f_cont = m * (wl - wl_blue) + b
    err_cont = sqrt((m_err * (wl - wl_blue))^2 + b_err^2)

    ;*******************************************************************************
    ; Compute equivalent width 
    ;*******************************************************************************

    wli = indexpars[iidx].bandpass
    if (indexpars[iidx].units eq 'A') then begin
    	eqwidth = integral(wl, (1 - flux/f_cont), wli[0], wli[1]) 

    	nerri =  (flux / f_cont) * sqrt((err / flux)^2 + (err_cont / f_cont)^2)
    	;eqwidth_err=sqrt(integral(wl, nerri^2, wli[0], wli[1]))
    	eqwidth_err=sqrt(integral(wl, dlambda*nerri^2, wli[0], wli[1]))
    	Lick = [eqwidth, eqwidth_err]
    endif    
    ;*******************************************************************************
    ; Compute magnitudes
    ;*******************************************************************************

    if (indexpars[iidx].units eq 'mag') then begin

    	nflux = integral(wl, (flux/f_cont), wli[0], wli[1]) 
    	nerri =  (flux / f_cont) * sqrt((err / flux)^2 + (err_cont / f_cont)^2)
    	;nerr =  sqrt(integral(wl, nerri^2, wli[0], wli[1])) 
    	nerr =  sqrt(integral(wl, dlambda*nerri^2, wli[0], wli[1]))
	mag = -2.5 * alog10(nflux / (wli[1] - wli[0])) 
    	mag_err = 2.5 / alog(10) * nerr / nflux
    	Lick=[mag,mag_err]
    endif   
    
    return,Lick
end


;*******************************************************************************
function d4000, wl, flux, err, narrow = narrow, redside=fl_r, blueside=fl_b
;*******************************************************************************

   if keyword_set(narrow) then begin
      ;; Balogh et al., 1999, ApJ, 527, 54
      wl_b1 = 3850.0  &  wl_b2 = 3950.0  
      wl_r1 = 4000.0  &  wl_r2 = 4100.0
   endif else begin
      ;; Bruzual, 1983, ApJ, 273, 105
      wl_b1 = 3750.0  &  wl_b2 = 3950.0  
      wl_r1 = 4050.0  &  wl_r2 = 4250.0
   endelse

   mn = min(wl, max=mx)
   d4000 = fltarr(2)
   fl_r = 0.0 & fl_b=0.0
   if (mn lt wl_b1 and mx gt wl_r2) then begin
;      if (abs(mn-2540.7239) lt 1e-3) then stop
      if (n_elements(wl) ne n_elements(flux)) then stop
      wi=(wl_b1-3.) < max(wl[where(wl lt wl_b1)])
      wf=(wl_r2+3.) > min(wl[where(wl gt wl_r2)])
      ;wi=(wl_b1-5.) < max(wl[where(wl lt wl_b1)])
      ;wf=(wl_r2+5.) > min(wl[where(wl gt wl_r2)])
      locali=where(wl gt (wi-5.) and wl lt (wf+5.))
      if (locali[0] ne -1) then begin
      	if (wl_b1 ge min(wl[locali]) and wl_r2 le max(wl[locali])) then begin
	    print,'computing D4000'
      	    ;the flux is defined as the integral of Fnu*dlambda=Flambda*lambda^2*dlambda
      	    fl_r = integral(wl[locali], flux[locali]*wl[locali]*wl[locali], wl_r1, wl_r2) 
      	    fl_b = integral(wl[locali], flux[locali]*wl[locali]*wl[locali], wl_b1, wl_b2) 
      
      	    fl_r_err =  sqrt(integral(wl[locali], (err[locali]*wl[locali]*wl[locali])^2, wl_r1, wl_r2)) 
      	    fl_b_err =  sqrt(integral(wl[locali], (err[locali]*wl[locali]*wl[locali])^2, wl_b1, wl_b2)) 
      
      	    d4000[0] = (wl_b2 - wl_b1) / (wl_r2 - wl_r1)  *  fl_r / fl_b
      	    d4000[1] = (wl_b2 - wl_b1) / (wl_r2 - wl_r1) * $
                 sqrt((fl_r_err/fl_b)^2 + (fl_b_err * fl_r/fl_b^2)^2)
	    ;stop
	endif
      endif
   endif
   return, d4000

end

;**************************************************
function Lick_models,iidx,indexpars,wl,flux
;**************************************************
    wlb = indexpars[iidx].blue_continuum
    f_blue_cont = integral(wl, flux, wlb[0], wlb[1]) / (wlb[1] - wlb[0])
    wl_blue = (wlb[0] + wlb[1]) / 2.0

    wlr = indexpars[iidx].red_continuum
    f_red_cont = integral(wl, flux, wlr[0], wlr[1]) / (wlr[1] - wlr[0])
    wl_red = (wlr[0] + wlr[1]) / 2.0

    m = (f_red_cont[0] - f_blue_cont[0]) / (wl_red - wl_blue)
    b = f_blue_cont[0]

    f_cont = m * (wl - wl_blue) + b

    ;*******************************************************************************
    ; Compute equivalent width 
    ;*******************************************************************************

    wli = indexpars[iidx].bandpass
    if (indexpars[iidx].units eq 'A') then begin
	eqwidth = integral(wl, (1 - flux/f_cont), wli[0], wli[1]) 
	Lick_models = [eqwidth]
    endif    
    ;*******************************************************************************
    ; Compute magnitudes
    ;*******************************************************************************

    if (indexpars[iidx].units eq 'mag') then begin

	nflux = integral(wl, (flux/f_cont), wli[0], wli[1]) 
	mag = -2.5 * alog10(nflux / (wli[1] - wli[0])) 
	Lick_models=[mag]
    endif   
    
    return,Lick_models
end


;*******************************************************************************
function d4000_models, wl, flux, narrow = narrow, redside=fl_r, blueside=fl_b
;*******************************************************************************

   if keyword_set(narrow) then begin
      ;; Balogh et al., 1999, ApJ, 527, 54
      wl_b1 = 3850.0  &  wl_b2 = 3950.0  
      wl_r1 = 4000.0  &  wl_r2 = 4100.0
   endif else begin
      ;; Bruzual, 1983, ApJ, 273, 105
      wl_b1 = 3750.0  &  wl_b2 = 3950.0  
      wl_r1 = 4050.0  &  wl_r2 = 4250.0
   endelse

   mn = min(wl, max=mx)
   fl_r = 0.0 & fl_b=0.0
   if (mn lt wl_b1 and mx gt wl_r2) then begin
;      if (abs(mn-2540.7239) lt 1e-3) then stop
      if (n_elements(wl) ne n_elements(flux)) then stop
      wi=(wl_b1-3.) < max(wl[where(wl lt wl_b1)])
      wf=(wl_r2+3.) > min(wl[where(wl gt wl_r2)])
      ;wi=(wl_b1-5.) < max(wl[where(wl lt wl_b1)])
      ;wf=(wl_r2+5.) > min(wl[where(wl gt wl_r2)])
      locali=where(wl gt (wi-5.) and wl lt (wf+5.))
      if (locali[0] ne -1) then begin
	if (wl_b1 ge min(wl[locali]) and wl_r2 le max(wl[locali])) then begin
	    ;the flux is defined as the integral of Fnu*dlambda=Flambda*lambda^2*dlambda
	    fl_r = integral(wl[locali], flux[locali]*wl[locali]*wl[locali], wl_r1, wl_r2) 
	    fl_b = integral(wl[locali], flux[locali]*wl[locali]*wl[locali], wl_b1, wl_b2) 
    
	    d4000_models = (wl_b2 - wl_b1) / (wl_r2 - wl_r1)  *  fl_r / fl_b
	endif
      endif
   endif
   return, d4000_models

end


pro measure_idx2,flux,err,lambda,fileidx_out
    common idxinfo

    ok=where(flux gt -90. and err gt 0. and err lt 100. and finite(flux) and finite(err),nok)
    ;; Measure absorption indices
    lidx=fltarr(2,34)
    nameidx=strarr(34)
    oki=where(indexpars.blue_continuum[0] ge min(lambda[ok]) and $
    	indexpars.red_continuum[1] le max(lambda[ok]) and (strcmp(indexpars.index,'Lick',4) eq 1 $
	or strcmp(indexpars.index,'H8',2) eq 1 or strcmp(indexpars.index,'H9',2) eq 1 $
	or strcmp(indexpars.index,'BH_HK',5) eq 1),noki)
    print,noki
    
    lick_idx=spaxel_index(lambda,flux,err,okpix=ok,index_idx=oki,minfrac_side=0.8,minfrac_cen=0.8)
    lidx[0,oki]=lick_idx.idx_value[0]
    lidx[1,oki]=lick_idx.idx_err[0]
    nameidx[oki]=lick_idx.idx_name[0]
    
;   for j=0,noki-1 do begin
;	print,'processing index ',indexpars[oki[j]].index
;	nameidx[oki[j]]=indexpars[oki[j]].index
;	wi=(indexpars[oki[j]].blue_continuum[0]-3.) < max(lambda[ok[where(lambda[ok] lt indexpars[oki[j]].blue_continuum[0])]])
;	wf=(indexpars[oki[j]].red_continuum[1]+3.) > min(lambda[ok[where(lambda[ok] gt indexpars[oki[j]].red_continuum[1])]])
;	;allow for missing pixels (5 = 10A)
;	locali=where(lambda[ok] gt (wi-10.) and lambda[ok] lt (wf+10.))
;	if (locali[0] ne -1) then begin
;	    if (indexpars[oki[j]].blue_continuum[0] ge min(lambda[ok[locali]]) and $
;		indexpars[oki[j]].red_continuum[1] le max(lambda[ok[locali]])) then begin
;		lickidx=Lick(oki[j],indexpars,lambda[ok[locali]],flux[ok[locali]],err[ok[locali]])
;		lidx[0,oki[j]]=lickidx[0] & lidx[1,oki[j]]=lickidx[1]
;	    endif else $
;		print,'index outside wavelentgh range'
;	endif
;   endfor
		
    if (min(lambda[ok]) le 3850.) then begin
       	;d4000_n = d4000(lambda[ok], flux[ok], err[ok], /narrow,redside=r_d4000, blueside=b_d4000)
	d4000_n = f_d4000(lambda, flux, err, /narrow,redside=r_d4000, blueside=b_d4000,okpix=ok,minfrac=0.8)
	lidx[0,28]=d4000_n[0] & lidx[1,28]=d4000_n[1]
	print,'D4',d4000_n
	nameidx[28]='D4000_n'
    endif
    
    print,'D4000: ',lidx[0,28],' pm ',lidx[1,28]
    print,'Hdelta: ',lidx[0,21],' pm ',lidx[1,21]
    print,'Hgamma: ',lidx[0,22],' pm ',lidx[1,22]
    print,'Hbeta: ',lidx[0,8],' pm ',lidx[1,8]
    
    ; Measure also composite indices
    ; [Mg1Fe] = 0.6*Mg1 + 0.4*alog10(Fe4531+Fe5015)
    ; [Mg2Fe] = 0.6*Mg2 + 0.4*alog10(Fe4531+Fe5015)
    ; [MgFe]' = sqrt(Mgb*(0.72*Fe5270 + 0.28*Fe5335))
    ; Mgb/<Fe> = Mgb / (0.5*(Fe5270 + Fe5335))
    ; Hd_Hg= Hdelta_A+Hgamma_A
    nameidx[29]='Mg1Fe'
    lidx[0,29]= 0.6*lidx[0,10] + 0.4*alog10(lidx[0,6]+lidx[0,9])
    lidx[1,29]= (0.4/(alog(10)*(lidx[0,6]+lidx[0,9])))^2.
    lidx[1,29]=sqrt(lidx[1,29]*((lidx[1,6])^2.+(lidx[1,9])^2.)+(0.6*lidx[1,10])^2.)
    nameidx[30]='Mg2Fe'
    lidx[0,30]= 0.6*lidx[0,11] + 0.4*alog10(lidx[0,6]+lidx[0,9])
    lidx[1,30]= (0.4/(alog(10)*(lidx[0,6]+lidx[0,9])))^2.
    lidx[1,30]=sqrt(lidx[1,30]*((lidx[1,6])^2.+(lidx[1,9])^2.)+(0.6*lidx[1,11])^2.)
    nameidx[31]='MgFe_prime'
    lidx[0,31]=sqrt(lidx[0,12]*(0.72*lidx[0,13]+0.28*lidx[0,14]))
    lidx[1,31]=(lidx[1,12]/lidx[0,12])^2.+((0.72*lidx[1,13])^2.+ $
    	(0.28*lidx[1,14])^2.)/(0.72*lidx[0,13]+0.28*lidx[0,14])^2.
    lidx[1,31]=0.5*lidx[0,31]*sqrt(lidx[1,31])
    nameidx[32]='Mg_Fe'
    lidx[0,32]=lidx[0,12]/(0.5*(lidx[0,13]+lidx[0,14]))
    lidx[1,32]=(lidx[1,12]/lidx[0,12])^2.+(lidx[1,13]^2.+lidx[1,14]^2.)/(lidx[0,13]^2.+lidx[0,14]^2.)
    lidx[1,32]=lidx[0,32]*sqrt(lidx[1,32])
    nameidx[33]='Hd_Hg'
    lidx[0,33] = lidx[0,21]+lidx[0,22]
    lidx[1,33] = sqrt(lidx[1,21]^2.+lidx[1,22]^2.)

    openw,lun,fileidx_out,/get_lun
    for i=0,22 do printf,lun,format='(2(1x,f10.4),2x,A)',lidx[0:1,i],nameidx[i] 
    for i=25,26 do printf,lun,format='(2(1x,f10.4),2x,A)',lidx[0:1,i],nameidx[i] 
    for i=28,33 do printf,lun,format='(2(1x,f10.4),2x,A)',lidx[0:1,i],nameidx[i] 
    free_lun,lun
    
    
    ;stop
    
end
