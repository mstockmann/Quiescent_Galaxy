@/Users/SteMaGE/tools/idl/gandalf/gandalf.pro
common paths, inst_dir

pro init_gandalf_LDRG
  common paths
  inst_dir='/Users/SteMaGE/tools/idl/LDRG/'
end

function glxlinemask,gndlf_linestr,lnlam_rest,mask_width=mask_width,NaI=NaI
  ;; return the mask array for the line in the galaxy's rest
  ;; frame (including NaI interstellar if required)
  ;; gndlf_linestr : structure with the line setup a la gandalf
  ;; lnlam_rest : ln(lambda) [approximately] corrected to restframe
  ;; (the spectrum is rebinned on this array)

  if (~(keyword_set(mask_width))) then mask_width=1500. ;; km/s default
  if (~(keyword_set(NaI))) then NaI=0 ;; NaI instrinsic absorptions are NOT masked out by default
  cspeed = 299792.458d

  if (NaI) then $
    wline2mask=where(~(gndlf_linestr.name eq 'sky') and (gndlf_linestr.action eq 'm')) $
  else wline2mask=where(~(gndlf_linestr.name eq 'sky') and ~(gndlf_linestr.name eq 'NaI') and (gndlf_linestr.action eq 'm'))
  nlambda_out=n_elements(lnlam_rest)

  mask_bluelim=alog(gndlf_linestr.lambda[wline2mask]*(1.-0.5*mask_width/cspeed))
  mask_redlim=alog(gndlf_linestr.lambda[wline2mask]*(1.+0.5*mask_width/cspeed))
  maskout=intarr(nlambda_out)
  for i=0, nlambda_out-1 do begin
    for k = 0, n_elements(wline2mask) - 1 do begin
      if (lnlam_rest[i] ge mask_bluelim[k] and lnlam_rest[i] le mask_redlim[k]) then maskout[i]=maskout[i] > 1
    endfor
  endfor

  return, maskout
end

function obslinemask,gndlf_linestr,lnlam_rest,zgal,mask_width=mask_width,NaI=NaI
  ;; returns the mask array for the lines in the obs frame (sky and
  ;; Galactic NaI)
  ;; gndlf_linestr : structure with the line setup a la gandalf
  ;; lnlam_rest : ln(lambda) [approximately] corrected to restframe
  ;; (the spectrum is rebinned on this array)
  ;; zgal : redshift of the galaxy (the same assumed to bring the spec
  ;; in restframe)
  if (~(keyword_set(mask_width))) then mask_width=1500. ;; km/s default
  if (~(keyword_set(NaI))) then NaI=0 ;; NaI Galactic absorptions are NOT masked out by default

  cspeed = 299792.458d

  if (NaI) then $
    wline2mask=where(((gndlf_linestr.name eq 'sky') or (gndlf_linestr.name eq 'NaI')) and (gndlf_linestr.action eq 'm')) $
  else wline2mask=where((gndlf_linestr.name eq 'sky') and (gndlf_linestr.action eq 'm'))
  nlambda_out=n_elements(lnlam_rest)

  mask_bluelim=alog(gndlf_linestr.lambda[wline2mask]*(1.-0.5*mask_width/cspeed)/(1.+zgal))
  mask_redlim=alog(gndlf_linestr.lambda[wline2mask]*(1.+0.5*mask_width/cspeed)/(1.+zgal))
  maskout=intarr(nlambda_out)
  for i=0, nlambda_out-1 do begin
    for k = 0, n_elements(wline2mask) - 1 do begin
      if (lnlam_rest[i] ge mask_bluelim[k] and lnlam_rest[i] le mask_redlim[k]) then maskout[i]=maskout[i] > 1
    endfor
  endfor

  return, maskout

end

pro go
  common paths
  INIT_GANDALF_LDRG

  ;  cd,'/Users/data/LDRGs/108899'
  ;  hdr=headfits('108899_avwCombined_OBx5_sig5_V1_NIRcorr_wmrebin15_opt_bpm.fits')
  ;  zgal=sxpar(hdr,'REDSHIFT')
  ;  dlambda_inst=sxpar(hdr,'CDELT3')
  ;  run_gandalf_LDRG,'108899_avwCombined_OBx5_sig5_V1_NIRcorr_wmrebin15_opt_bpm.fits',$
  ;    '108899_avwCombined_OBx5_sig5_V1_NIRcorr_wmrebin15_opt_bpm',zgal,$
  ;    /out_sens,AoN_threshold=2.0,mdegree=3,models='BC03std',$ ;models='MILES',$
  ;    reskern=100,dlambda_inst=dlambda_inst,minlambda=11000.,maxlambda=22000.
  ;
  ;stop
  cd,'/Users/data/LDRGs/105842'
  root_name='105842_V1_NIRcorr_wmrebin15_er2d_opt_bpm'
  hdr=headfits(root_name+'.fits')
  zgal=sxpar(hdr,'REDSHIFT')
  dlambda_inst=sxpar(hdr,'CDELT1')
  run_gandalf_LDRG,root_name+'.fits',$
    root_name,zgal,$
    AoN_threshold=2.0,mdegree=3,models='MILES',$ ;models='BC03std',$ ;
    dlambda_inst=dlambda_inst,minlambda=11000.,maxlambda=22000.,emission_setup_file='emission_lines_setup_105842.txt'

  stop
  cd,'/Users/data/LDRGs/test'
  run_gandalf_LDRG,'nodding_macs2129_full_nir_9.6.main_box.fits',$
    'nodding_macs2129_full_nir_9.6.main_box_2',2.1489,$
    /out_sens,AoN_threshold=2.0,mdegree=3,models='BC03std',$ ;models='MILES',$
    reskern=100,dlambda_inst=9.6,minlambda=12000.,maxlambda=22000.
  stop

  cd,'/Users/data/LDRGs/MACS2129/reduced201509'
  run_gandalf_MACS,'nodding_macs2129_full_nir_9.6.main_box.fits',$
    'nodding_macs2129_full_nir_9.6.main_box',zgal=2.1489,$
    /out_sens,AoN_threshold=2.0,mdegree=3,models='BC03std',$
    reskern=100,dlambda_inst=9.6,maxlambda=25000.
  stop
  run_gandalf_MACS,'MACS2129-2_nir_master_stack_opt_3sigma_5iter_allcor.rebin16.opt.centerforgrad_box.1d.fits',$
    'MACS2129-2_nir_master_stack_opt_3sigma_5iter_allcor.rebin16.opt.centerforgrad_box.1d',zgal=2.1489,$
    /out_sens,AoN_threshold=2.0,mdegree=3,models='BC03std',$
    reskern=100,dlambda_inst=9.6,maxlambda=25000.

  run_gandalf_MACS,'MACS2129-2_nir_master_stack_opt_3sigma_5iter_allcor.rebin16.opt.outerforgrad_box.1d.fits',$
    'MACS2129-2_nir_master_stack_opt_3sigma_5iter_allcor.rebin16.opt.outerforgrad_box.1d',zgal=2.1489,$
    /out_sens,AoN_threshold=2.0,mdegree=3,models='BC03std',$
    reskern=100,dlambda_inst=9.6,maxlambda=25000.

  run_gandalf_MACS,'MACS2129-2_nir_master_stack_opt_3sigma_5iter_allcor.rebin16.opt.lowerforgrad_box.1d.fits',$
    'MACS2129-2_nir_master_stack_opt_3sigma_5iter_allcor.rebin16.opt.lowerforgrad_box.1d',zgal=2.1489,$
    /out_sens,AoN_threshold=2.0,mdegree=3,models='BC03std',$
    reskern=100,dlambda_inst=9.6,maxlambda=25000.

  run_gandalf_MACS,'MACS2129-2_nir_master_stack_opt_3sigma_5iter_allcor.rebin16.opt.upperforgrad_box.1d.fits',$
    'MACS2129-2_nir_master_stack_opt_3sigma_5iter_allcor.rebin16.opt.upperforgrad_box.1d',zgal=2.1489,$
    /out_sens,AoN_threshold=2.0,mdegree=3,models='BC03std',$
    reskern=100,dlambda_inst=9.6,maxlambda=25000.
  stop

end



pro run_gandalf_LDRG,infile,outfile_prefix,zgal,AoN_threshold=AoN_threshold,plot=plot,mdegree=mdegree,$
  cheb_deg=cheb_deg,models=models,$
  minlambda=minlambda,maxlambda=maxlambda,sigma_emline_0=sigma_emline_0,$
  dlambda_inst=dlambda_inst,$
  emission_setup_file=emission_setup_file
  ; zgal is the reference systemic redshift for the galaxy. If not provided look up 'AZMOOTH3C ZGAL' HIERARCH keyword
  ;; AoN_threshold: Amplitude over Noise threshold for valid emission
  ;; line detection, defaults to 2 (not 4 as in GANDALF examples)
  ;; cheb_deg indicates the degree of the chebychev polynomial used to fit the response residual, defaults to 30
  ;
  ; maxlambda: maximum OBSERVED lambda to consider in the cubes; since pipeline V1.4 lambda>7050 are not
  ;     flux-calibrated in the COMBO cubes and have strange behaviours in the error cubes. Defaults to 7000.
  ;     Set to <=0 to apply no cut
  ; minlambda: minimum OBSERVED lambda to consider in the cubes

  common paths
  cspeed = 299792.458d
  ;; lambda ranges at REST
  ;; consider COMBO cubes are "complete" between 3700 and ~7000 AA
  ;    ;; fit kinematics from H+K to last Fe5782
  ;    loglam_min_kin=alog(3910.)
  ;    loglam_max_kin=alog(5815.)
  ;    ;; fit kinematics from approx H+K only
  ;    loglam_min_kin=alog(3800.)
  ;    loglam_max_kin=alog(4200.)
  ;loglam_min_gnd=alog(3910.)
  ;loglam_max_gnd=alog(6700.)
  ;; fit kinematics from approx H+K to Mgb
  loglam_min_kin=alog(3500.) ;alog(3800.)
  loglam_max_kin=alog(7000.) ;alog(6000.);alog(5350.)
  ;loglam_min_gnd=alog(3910.)
  ;loglam_max_gnd=alog(6700.)
  loglam_min_gnd=alog(3500.)
  loglam_max_gnd=alog(7000.)   ;loglam_max_gnd=alog(6740.) ;; stop right after second SII
  if (n_elements(AoN_threshold) eq 0) then AoN_threshold=2.0
  if (n_elements(mdegree) eq 0) then mdegree=3
  if (n_elements(cheb_deg) eq 0) then cheb_deg=30

  if n_elements(dlambda_inst) eq 0 then begin
    dlambda_inst=6. ;; for CALIFA
    ; dlambda_inst=2.8 ;; for MUSE!!
  endif

  if (n_elements(reskern) eq 0) then reskern=100
  if (reskern lt 0) then begin
    reskern=0
    plotfirst=1
  endif
  ker_w=reskern
  ker_w2=fix(ker_w/2)

  ;; minimum SNR required in order to perform GANDALF.
  ;Spaxel's SNR is read from SNR_MAP extension of the cube output
  ;by azmooth3C (with option -k enabled)
  ; if the SNR_MAP extension is not present then run GANDALF on everything which has SMASK value smaller than 100

  if (n_elements(sigma_emline_0) eq 0) then sigma_emline_0=20. ;; initial velocity dispersion for emission lines: use 20. for CALIFA, 10. for MUSE

  if (n_elements(models) eq 0) then models='BC03std'

  if (models eq 'BC03std') then begin
    burst_model_file=inst_dir+'bc_models_subset_v4_0c_bursts_new.fits'
    sigma_mod=(3./2.35)  ;; bc03
    dlambda_key='CD1_1'
  endif

  if (models eq 'BC03ext') then begin
    burst_model_file=inst_dir+'bc_models_subset_v5_bursts_new.fits'
    sigma_mod=(3./2.35)  ;; bc03
    dlambda_key='CD1_1'
  endif

  if (models eq 'MILES') then begin
    burst_model_file=inst_dir+'bc_models_subset_v4_0cMILES_bursts_new.fits'
    sigma_mod=2.5/2.35 ;; MILES
    dlambda_key='CDELT1'
  endif

  if (n_elements(maxlambda) eq 0) then maxlambda=7000.0

  if n_elements(emission_setup_file) eq 0 then emission_setup_file=inst_dir+'emission_lines_setup_std_fullspec_freeall.txt'


  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;   Read data
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  file=infile
  if ~(file_test(file,/regular)) then if (file_test(file+'.gz',/regular)) then file=file+'.gz' else message,'Input file not found '+file
  flux=mrdfits(file,0,hdr_cube)
  ;if (n_elements(zgal) eq 0) then zgal=float(sxpar_hierarch(hdr_cube,'AZMOOTH3C ZGAL'))
  lambda0=sxpar(hdr_cube,'CRVAL1')
  dlambda=sxpar(hdr_cube,'CDELT1')
  nlambda=sxpar(hdr_cube,'NAXIS1')
  lambda=lambda0+findgen(nlambda)*dlambda
  error=mrdfits(file,'ERROR')
  flag=mrdfits(file,'BADPIX')
  smask=mrdfits(file,'SMASK',status=mrdfits_status)
  if mrdfits_status ne 0 then smask=1

  if (maxlambda gt 0.) then begin
    ;; redefine cubes in the relevant range
    maxzcube=fix((maxlambda-lambda0)/dlambda) < (nlambda-1)
  endif
  if (minlambda gt 0.) then begin
    ;; redefine cubes in the relevant range
    minzcube=0 > fix((minlambda-lambda0)/dlambda)
  endif

  if maxlambda gt 0. or minlambda gt 0. then begin
    nlambda=maxzcube-minzcube+1
    lambda0=lambda0+minzcube*dlambda
    lambda=lambda0+findgen(nlambda)*dlambda
    flux=flux[minzcube:maxzcube];*1e19
    error=error[minzcube:maxzcube];*1e19
    flag=flag[minzcube:maxzcube]
  endif

  ;;; TMP
  ;;flag=intarr(nra,ndec,nlambda)
  ;flag[*,*,where(lambda ge 13550 and lambda le 14240)]=1
  ;flag[*,*,where(lambda ge 18000 and lambda le 19500)]=1

  ;stop
  ;;;



  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;   Load model arrays
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  n_met = 4

  zmod = fltarr(n_met)

  for imod=0,n_met-1 do begin
    hdr = headfits(burst_model_file, exten=imod)
    zmod[imod] = sxpar(hdr, 'Z')
  endfor


  ;; Get header information for wavelength arrays
  n_pixmodel = sxpar(hdr, 'NAXIS1')
  ;; 2013-09-06 temporary change to exclude model #2 from the library which is all NULL at all Z
  n_burst = sxpar(hdr, 'NAXIS2') ;; original to restore after problem with lib is fixed
  ;n_burst = sxpar(hdr, 'NAXIS2')-1
  wl0 = sxpar(hdr, 'CRVAL1')
  dw = sxpar(hdr, dlambda_key)
  burst_wl = findgen(n_pixmodel) * dw + wl0

  burst_lib = fltarr(n_pixmodel, n_burst, n_met)
  gburst_lib=burst_lib
  ssp_norms = fltarr(n_burst, n_met)
  if strmatch(burst_model_file, '*burst*') then $
    all_norms = (mrdfits(burst_model_file, n_met));;[*,[0,1,3,4,5,6,7,8,9,10,11,12,13,14]] ;; remove index vector after MILES lib is fixed!!

  ;; 2013-09-06 temporary change to exclude model #2 from the library which is all NULL at all Z
  ;; original to restore after problem in lib is fixed
  for imod=0, n_met-1 do begin
    tmp = mrdfits(burst_model_file, imod)
    burst_lib[*, *, imod] = tmp[*,*]

    if strmatch(burst_model_file, '*burst*') then begin
      ssp_norms[*, imod] = reform(all_norms[imod, *])
    endif

  endfor


  ;; Convolve models to match spectral resolution of the data
  ;rememeber to define sigma_mod above where the library file is defined
  sigma_obs=dlambda_inst/2.35/(1.+zgal)
  sigma_add=sqrt((sigma_obs^2.-sigma_mod^2.) > 0.001)
  sigma_add=sigma_add/dw        ; sigma in pixel (model)
  lsf=psf_Gaussian(NPIXEL=2*ceil(4*sigma_add)+1, ST_DEV=sigma_add, /NORM, NDIM=1)
  for k=0,n_burst-1 do begin
    for kk=0,n_met-1 do begin
      gburst_lib[*,k,kk] = convol(burst_lib[*,k,kk],lsf)
    endfor
  endfor

  ;stop

  ;------------------------------------------------------------------
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ; read in emission-line setup files
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  print,'--> reading in the emission-line setup file, and creating the corresponding structure'
  ;emission_setup_file='/Users/zibetti/idl/CALIFA/emission_lines_setup_std.txt'
  ;emission_setup_file='/Users/SteMaGE/tools/idl/CALIFA/emission_lines_setup_std_fullspec.txt'
  ;emission_setup_file=inst_dir+'emission_lines_setup_std_fullspec_OIII_kin.txt' & message,'line kinematics tied to [OIII] 5006.77',/continue
  eml_file = emission_setup_file
  readcol,eml_file,eml_i,eml_name,eml_lambda,eml_action,eml_kind,eml_a,eml_v,eml_s,eml_fit,$
    f='(i,a,f,a,a,f,f,f,a)',skipline=2,comment='#',/silent
  ; creating emission setup structure, to be used for masking the
  ; emission-line contaminated regions and fitting the emission lines in
  ; GANDALF
  emission_setup_0 = create_struct('i',eml_i,'name',eml_name,'lambda',eml_lambda,'action',eml_action,$
    'kind',eml_kind,'a',eml_a,'v',eml_v,'s',eml_s,'fit',eml_fit)

  spec=reform(flux[*])      ;/1.e-17
  ;; Log-rebin data array
  restwl=lambda/(1.+zgal)
  lamRange=[restwl[0],restwl[nlambda-1]]  ;; won't be modified
  if (n_elements(velscale) gt 0) then undefine,velscale
  log_rebin,lamRange, spec, spec_rebin, loglam_rest, velscale=velscale  ;;velscale will not be modified anymore
  nlambda_out=n_elements(loglam_rest)

  ;; log-rebin the templates as the data
  lamRange_models=[burst_wl[0],burst_wl[n_pixmodel-1]]
  log_rebin, lamRange_models, gburst_lib[*,0,0], ssp_rebin, loglam_ssp, velscale=velscale
  templates=dblarr(n_elements(ssp_rebin), n_burst*n_met)
  Ha_templ=where(loglam_ssp lt alog(6563.+150.) and loglam_ssp gt alog(6563.-150.))
  templates_Ha=dblarr(n_elements(Ha_templ), n_burst*n_met)
  OIII_templ=where(loglam_ssp lt alog(5007.+150.) and loglam_ssp gt alog(5007.-150.))
  templates_OIII=dblarr(n_elements(OIII_templ), n_burst*n_met)
  for k=0,n_burst-1 do begin
    for kk=0,n_met-1 do begin
      log_rebin, lamRange_models, gburst_lib[*,k,kk], ssp_rebin, velscale=velscale
      kkk=n_met*k+kk
      templates[*,kkk]=ssp_rebin
      templates_Ha[*,kkk]=ssp_rebin[Ha_templ]
      templates_OIII[*,kkk]=ssp_rebin[OIII_templ]
    endfor
  endfor

  ;; run ppxf once on the central pixel to set velscale for log_rebin and set initial guess
  noise=reform(error[*])    ;/1.e-17
  flags=reform(flag[*])
  ;stop
  if (n_elements(velscale) gt 0) then undefine,velscale_tmp
  log_rebin,lamRange, noise^2., noise_rebin, velscale=velscale_tmp
  if velscale_tmp ne velscale then message,'velscales do not correspond in spec and noise, pls check!'
  noise_rebin=sqrt(noise_rebin)
  ;;; print,velscale,n_elements(loglam_rest),n_elements(spec),lamRange

  ;; log rebin bad pixels mask
  flags_logrebin=intarr(nlambda_out)
  ii=0
  for jj=0,nlambda_out-1 do begin
    while ((alog(restwl[ii]) lt loglam_rest[jj]) and (ii lt n_elements(restwl)-1)) do begin
      flags_logrebin[jj]=flags_logrebin[jj]>flags[ii]
      ii=ii+1
    endwhile
    if (ii lt n_elements(restwl)-1) then if ~(flags[ii] eq 0) then flags_logrebin[jj]=1
  endfor


  ;; mask regions
  linesetup=emission_setup_0
  not_HaNII=where(emission_setup_0.lambda lt 6547.96 or emission_setup_0.lambda gt 6583.34,complement=HaNII)
  Ha=where(emission_setup_0.name eq 'Ha')

  OIII=where(emission_setup_0.lambda eq 5006.7700,complement=not_OIII)


  restmask=glxlinemask(linesetup,loglam_rest,/NaI)
  obsmask=obslinemask(linesetup,loglam_rest,zgal,/NaI)
  okpix=where(flags_logrebin eq 0 and restmask eq 0 and obsmask eq 0 $
    and loglam_rest gt loglam_min_kin and loglam_rest lt loglam_max_kin,complement=badpix) ;; from H+K to last Fe5782
  ;;        and exp(loglam_rest) gt 3750. and exp(loglam_rest) lt 5500.,complement=badpix)
  okpix_all=where(restmask eq 0 and obsmask eq 0 $
    and loglam_rest gt loglam_min_kin and loglam_rest lt loglam_max_kin,n_kin_all)
  Ha_pix=where(loglam_rest lt alog(6563.+100.) and loglam_rest gt alog(6563.-100.),npix_Ha)
  OIII_pix=where(loglam_rest lt alog(5007.+100.) and loglam_rest gt alog(5007.-100.),npix_OIII)

  dv=(loglam_ssp[0] - loglam_rest[0])*cspeed ;; velocity offset between model and data
  dv_Ha=(loglam_ssp[Ha_templ[0]] - loglam_rest[Ha_pix[0]])*cspeed ;; velocity offset between model and data
  dv_OIII=(loglam_ssp[OIII_templ[0]] - loglam_rest[OIII_pix[0]])*cspeed ;; velocity offset between model and data

  time0_ppxf=systime(/seconds)
  start=[0.,70.]
  ppxf, templates, spec_rebin, noise_rebin, velscale, start, ppxf_sol, $
    moments=2, vsyst=dv, $
    goodpixels=okpix, error=ppxf_error, degree=0,mdegree=mdegree,$
    bestfit=ppxf_bestfit                      ;, /plot, /clean;, /oversample ;, degree=3;,mdegree=3
  time1_ppxf=systime(/seconds)
  print,'pPXF exec. time=',time1_ppxf-time0_ppxf


  ;; Create output structure
  ;; pPXF output
  ;;z_bestfit=fltarr(nra,ndec)
  Vsys=!values.f_nan
  Vsys_err=!values.f_nan
  Vdisp=!values.f_nan
  Vdisp_err=!values.f_nan
  chi2dof=!values.f_nan
  ;; GANDALF output
  net_stellar=fltarr(nlambda_out)
  net_stellar_allsub=fltarr(nlambda_out)
  best_stellar=fltarr(nlambda_out)
  best_nebular_all=fltarr(nlambda_out)
  best_nebular_det=fltarr(nlambda_out)
  gndf_badpix=intarr(nlambda_out)
  gndf_resid_rms=!values.f_nan

  i_lines = where(emission_setup_0.action eq 'm' and emission_setup_0.name ne 'sky'$
    and emission_setup_0.name ne 'NaI',nlines2fit_all)
  i_lines_indep = where(emission_setup_0.action eq 'm' and emission_setup_0.name ne 'sky'$
    and emission_setup_0.name ne 'NaI' and emission_setup_0.kind eq 'l',nlines2fit)
  em_flux=fltarr(nlines2fit)
  em_ampl=fltarr(nlines2fit)
  em_vsys=fltarr(nlines2fit)
  em_disp=fltarr(nlines2fit)
  em_flux_err=fltarr(nlines2fit)
  em_ampl_err=fltarr(nlines2fit)
  em_vsys_err=fltarr(nlines2fit)
  em_disp_err=fltarr(nlines2fit)
  em_AoN=fltarr(nlines2fit)
  k_ha=where(emission_setup_0.name[i_lines_indep] eq 'Ha')
  k_OIII=where(emission_setup_0.lambda[i_lines_indep] eq 5006.7700)


  linesetup.fit[i_lines] = emission_setup_0.fit[i_lines]
  ;; if emission lines too noisy, one can decide to keep the stellar kinematics
  if 1 then begin
    gal_lines=where(emission_setup_0.I[i_lines] lt 90)
    linesetup.fit[i_lines[gal_lines]]='h'
  endif


  ;; prepare to run gandalf
  linesetup.v[i_lines] = ppxf_sol[0] ;; vsys same as stars

  ;;; possible options to constrain sigma of emission lines
  ;; linesetup.s[i_lines] = sigma_emline_0 ;20. ;; essentially no intrinsic broadening, only instrumental
  linesetup.s[i_lines] = ppxf_sol[1] ;; same as stellar sigma
  restmask_em=glxlinemask(linesetup,loglam_rest,/NaI)


  ;; check for lines that should not be fitted
  linecen=value_locate(loglam_rest,alog((linesetup.lambda[i_lines])*(1.+ppxf_sol[0]/cspeed)))
  to_ignore=where(flags_logrebin[linecen] ne 0 or flags_logrebin[linecen+1] ne 0 or flags_logrebin[linecen-1] ne 0 $
    or obsmask[linecen] ne 0 or obsmask[linecen+1] ne 0 or obsmask[linecen-1] ne 0$
    or loglam_rest[linecen] lt loglam_min_gnd or loglam_rest[linecen] gt loglam_max_gnd $
    ;or spec_rebin[linecen] lt 3.*noise_rebin[linecen]
    ,ntoignore)

  linesetup.action[i_lines] = 'f'
  if ntoignore gt 0 then linesetup.action[i_lines[to_ignore]]='i'
  goodlines2fit=where(linesetup.action eq 'f' and emission_setup_0.name ne 'sky'$
    and emission_setup_0.name ne 'NaI' and emission_setup_0.kind eq 'l',ngoodlines2fit)

  okpix=where(flags_logrebin eq 0 and obsmask eq 0 $
    and loglam_rest gt loglam_min_gnd and loglam_rest lt loglam_max_gnd $
    and spec_rebin gt 3.*noise_rebin,nok,complement=badpix)
  time0_gandalf=systime(/seconds)

  gandalf_sol=ppxf_sol
  gandalf_sol[0]=gandalf_sol[0]+dv
  GANDALF, templates, spec_rebin, noise_rebin, velscale, gandalf_sol, linesetup, $
    loglam_rest[0], loglam_rest[1]-loglam_rest[0], GOODPIXELS=okpix, INT_RES=dlambda_inst/(1.+zgal)/2.35,$ ; INT_DISP=(6./(1.+zgal)/2.35/4861.*cspeed), $
    BESTFIT=gandalf_bestfit, EMISSION_TEMPLATES=emission_templates,$
    L0_TEMPL=loglam_ssp[0],$; /PLOT, $; =plotfirst,
    mdegree=mdegree,degree=0, $
    FOR_ERRORS=1, ERROR=gandalf_error;,/quiet

  resid=gandalf_bestfit-spec_rebin
  resid_noise=robust_sigma(resid[okpix],/ZERO)
  ;; IMPORTANT gandalf_sol is an array [4*nlines2fit+mdegree]
  ;; the first 4*nlines2fit are the parameters of the emission lines:
  ;; last mdegree elements are the coefficients of the multiplicative polynomial (not used in what follows)

  jj=0
  for ii=0,nlines2fit-1 do begin
    if linesetup.action[i_lines_indep[ii]] eq 'f' then begin
      em_flux[ii] = gandalf_sol[jj*4+0]
      em_ampl[ii] = gandalf_sol[jj*4+1]
      em_vsys[ii] = gandalf_sol[jj*4+2]
      em_disp[ii] = gandalf_sol[jj*4+3]
      em_flux_err[ii] = gandalf_error[jj*4+0]
      em_ampl_err[ii] = gandalf_error[jj*4+1]
      em_vsys_err[ii] = gandalf_error[jj*4+2]
      em_disp_err[ii] = gandalf_error[jj*4+3]

      ;; smooth error cube to use for noise computation
      ok_flag=intarr(nlambda_out)
      ok_flag[okpix]=1
      noise_kerw2=10
      dloglam_rest=loglam_rest[1]-loglam_rest[0]
      cen_pix=fix((alog(linesetup.LAMBDA[i_lines_indep[ii]]*(1.+gandalf_sol[jj*4+2]/cspeed))-loglam_rest[0])/dloglam_rest)
      tbsmoothed=[]
      nsmooth=0
      if (ok_flag[cen_pix] eq 1) then begin
        tbsmoothed=[tbsmoothed,noise_rebin[cen_pix]]
        nsmooth=nsmooth+1
      endif
      for kk=1,noise_kerw2 do begin
        if (ok_flag[cen_pix+kk] eq 1) then begin
          tbsmoothed=[tbsmoothed,noise_rebin[cen_pix]]
          nsmooth=nsmooth+1
        endif
        if (ok_flag[cen_pix-kk] eq 1) then begin
          tbsmoothed=[tbsmoothed,noise_rebin[cen_pix]]
          nsmooth=nsmooth+1
        endif
      endfor
      if (nsmooth gt 0) then local_noise=median(tbsmoothed) else local_noise=1.e10

      em_AoN[ii]=gandalf_sol[jj*4+1]/local_noise

      jj=jj+1
    endif else begin
      em_flux[ii] = !values.f_nan
      em_ampl[ii] = !values.f_nan
      em_vsys[ii] = !values.f_nan
      em_disp[ii] = !values.f_nan
      em_AoN[ii]  = !values.f_nan
    endelse
  endfor
  spec_emission_detect=dblarr(nlambda_out)
  spec_emission_all=dblarr(nlambda_out)

  i_lines_OII3727=where((emission_setup_0.lambda[goodlines2fit] eq 3726.03 or emission_setup_0.lambda[goodlines2fit] eq 3728.73) $
    and emission_setup_0.kind[goodlines2fit] eq 'l',nlinesOII3727)
  if nlinesOII3727 lt 2 then i_lines_OII3727=[-1,-1]
  ;; smooth error cube to use for AoN cut
  ok_flag=intarr(nlambda_out)
  ok_flag[okpix]=1
  noise_kerw2=10
  dloglam_rest=loglam_rest[1]-loglam_rest[0]
  ii=0
  for ii=0,ngoodlines2fit-1 do begin
    cen_pix=fix((alog(linesetup.LAMBDA[goodlines2fit[ii]]*(1.+gandalf_sol[ii*4+2]/cspeed))-loglam_rest[0])/dloglam_rest)
    tbsmoothed=[]
    nsmooth=0
    if (ok_flag[cen_pix] eq 1) then begin
      tbsmoothed=[tbsmoothed,noise_rebin[cen_pix]]
      nsmooth=nsmooth+1
    endif
    for kk=1,noise_kerw2 do begin
      if (ok_flag[cen_pix+kk] eq 1) then begin
        tbsmoothed=[tbsmoothed,noise_rebin[cen_pix]]
        nsmooth=nsmooth+1
      endif
      if (ok_flag[cen_pix-kk] eq 1) then begin
        tbsmoothed=[tbsmoothed,noise_rebin[cen_pix]]
        nsmooth=nsmooth+1
      endif
    endfor
    if (nsmooth gt 0) then local_noise=median(tbsmoothed) else local_noise=1.e10

    if (ii eq i_lines_OII3727[0] or ii eq i_lines_OII3727[1]) then begin ;; for the OII-3727 doublet check the sum of the two
      if (ii eq i_lines_OII3727[0]) then begin ;; make the sum as the first line of the doublet is encountered
        if ((gandalf_sol[ii*4+1]+gandalf_sol[i_lines_OII3727[1]*4+1])/local_noise ge AoN_threshold) then spec_emission_detect = spec_emission_detect + emission_templates[*,ii]
        spec_emission_all = spec_emission_all + emission_templates[*,ii]
      endif
    endif else begin
      if (gandalf_sol[ii*4+1]/local_noise ge AoN_threshold) then begin
        ;print,linesetup.name[goodlines2fit[ii]],gandalf_sol[ii*4+1],local_noise
        spec_emission_detect = spec_emission_detect + emission_templates[*,ii]
      endif
      spec_emission_all = spec_emission_all + emission_templates[*,ii]
    endelse
  endfor
  ;; fill in remaining output structures
  net_stellar[*] = spec_rebin - spec_emission_detect
  net_stellar_allsub[*] = spec_rebin - spec_emission_all
  best_stellar[*] = (gandalf_bestfit - spec_emission_all)
  best_nebular_all[*] = spec_emission_all
  best_nebular_det[*] = spec_emission_detect
  gndf_badpix[badpix] = 1
  gndf_resid_rms = resid_noise
  Vsys = ppxf_sol[0]
  Vdisp = ppxf_sol[1]
  Vsys_err = ppxf_error[0]*(ppxf_sol[6])
  Vdisp_err = ppxf_error[1]*sqrt(ppxf_sol[6])


  chi2dof = ppxf_sol[6]

  ;; output fits files
  ;;; stellar component
  hdr_cHDU=hdr_cube
  out_stellar=outfile_prefix+'_GNDF_stellar.fits'
  get_date,dte
  fxaddpar, hdr_cHDU,'AUTHOR','GallazZibetti'
  fxaddpar, hdr_cHDU,'DATE',dte, 'Date in yyyy-mm-dd'
  fxaddpar, hdr_cHDU, 'REDSHIFT',zgal
  fxaddpar, hdr_cHDU, 'CONTENT','Net stellar emission=data-det_nebular'
  fxaddpar, hdr_cHDU,'EXTNAME', 'NET_STELLAR'
  fxaddpar, hdr_cHDU,'NAXIS1', nlambda_out
  fxaddpar, hdr_cHDU,'CRPIX1', 1
  fxaddpar, hdr_cHDU,'CRVAL1', loglam_rest[0]
  fxaddpar, hdr_cHDU,'CDELT1', loglam_rest[1]-loglam_rest[0]
  fxaddpar, hdr_cHDU,'CUNIT1', 'ln(lambda/AA)'
  fxaddpar, hdr_cHDU,'z_GAL', zgal,'Galaxy redshift from cat'
  fxaddpar, hdr_cHDU,'z_SYS', zgal,'Redshift of the reference frame'
  fxaddpar, hdr_cHDU,'GD_AONTH', AoN_threshold,'Detection threshold for lines Amplitude/Noise'
  fxaddpar, hdr_cHDU,'GD_EMSET', emission_setup_file,'Emission line setup file'

  mwrfits,net_stellar,out_stellar,hdr_cHDU,/create

  sxdelpar,hdr_cHDU,'SIMPLE'
  fxaddpar, hdr_cHDU, 'CONTENT','Net stellar emission=data-all_nebular'
  fxaddpar, hdr_cHDU,'EXTNAME', 'NET_STELLAR_ALL'
  mwrfits,net_stellar_allsub,out_stellar,hdr_cHDU

  sxdelpar,hdr_cHDU,'SIMPLE'
  fxaddpar, hdr_cHDU, 'CONTENT','Bestfit stellar continuum'
  fxaddpar, hdr_cHDU,'EXTNAME', 'BEST_STELLAR'
  mwrfits,best_stellar,out_stellar,hdr_cHDU

  fxaddpar, hdr_cHDU, 'CONTENT','Masked pixel cube in GANDALF'
  fxaddpar, hdr_cHDU,'EXTNAME', 'BPM_GANDALF'
  mwrfits,gndf_badpix,out_stellar,hdr_cHDU

  fxaddpar,hdr_cHDU,'NAXES',2
  sxdelpar,hdr_cHDU,['NAXIS3','CRVAL3','CDELT3','CRPIX3','CUNIT3','DISPAXIS']

  fxaddpar, hdr_cHDU, 'CONTENT','L.o.S. velocity from pPXF'
  fxaddpar, hdr_cHDU,'EXTNAME', 'VSYS_PPXF'
  mwrfits,[vsys],out_stellar,hdr_cHDU

  fxaddpar, hdr_cHDU, 'CONTENT','Error on L.o.S. velocity from pPXF'
  fxaddpar, hdr_cHDU,'EXTNAME', 'VSYSERR_PPXF'
  mwrfits,[vsys_err],out_stellar,hdr_cHDU

  fxaddpar, hdr_cHDU, 'CONTENT','L.o.S. dispersion from pPXF'
  fxaddpar, hdr_cHDU,'EXTNAME', 'VDISP_PPXF'
  mwrfits,[vdisp],out_stellar,hdr_cHDU

  fxaddpar, hdr_cHDU, 'CONTENT','Error on L.o.S. dispersion from pPXF'
  fxaddpar, hdr_cHDU,'EXTNAME', 'VDISPERR_PPXF'
  mwrfits,[vdisp_err],out_stellar,hdr_cHDU

  fxaddpar, hdr_cHDU, 'CONTENT','Chisquared per dof from pPXF'
  fxaddpar, hdr_cHDU,'EXTNAME', 'CHI2DOF_PPXF'
  mwrfits,[chi2dof],out_stellar,hdr_cHDU

  fxaddpar, hdr_cHDU, 'CONTENT','Residual rms from GANDALF'
  fxaddpar, hdr_cHDU,'EXTNAME', 'GNDF_RMS'
  mwrfits,[gndf_resid_rms],out_stellar,hdr_cHDU

;  fxaddpar, hdr_cHDU, 'CONTENT','FLAGS from run GANDALF'
;  fxaddpar, hdr_cHDU,'EXTNAME', 'GNDF_FLAGS'
;  mwrfits,[gndf_flags],out_stellar,hdr_cHDU

  if keyword_set(noise_resc) then begin
    fxaddpar, hdr_cHDU, 'CONTENT','Noise scale factor based on fit residuals'
    fxaddpar, hdr_cHDU,'EXTNAME', 'NOISE_SCALE_FACTOR'
    mwrfits,gndf_rescnoise,out_stellar,hdr_cHDU
  endif

  ;;; nebular component
  hdr_cHDU=hdr_cube
  out_nebular=outfile_prefix+'_GNDF_nebular.fits'
  get_date,dte
  fxaddpar, hdr_cHDU,'AUTHOR','GallazZibetti'
  fxaddpar, hdr_cHDU,'DATE',dte, 'Date in yyyy-mm-dd'
  fxaddpar, hdr_cHDU, 'REDSHIFT',zgal
  fxaddpar, hdr_cHDU, 'CONTENT','DETECTED Nebular emission from GANDALF'
  fxaddpar, hdr_cHDU,'EXTNAME', 'NEBULAR_DET'
  fxaddpar, hdr_cHDU,'NAXIS1', nlambda_out
  fxaddpar, hdr_cHDU,'CRPIX1', 1
  fxaddpar, hdr_cHDU,'CRVAL1', loglam_rest[0]
  fxaddpar, hdr_cHDU,'CDELT1', loglam_rest[1]-loglam_rest[0]
  fxaddpar, hdr_cHDU,'CUNIT1', 'ln(lambda/AA)'
  fxaddpar, hdr_cHDU,'z_GAL', zgal,'Galaxy redshift from cat'
  fxaddpar, hdr_cHDU,'z_SYS', zgal,'Redshift of the reference frame'

  zlayer=0
  for ii=0,nlines2fit-1 do begin
    if (linesetup.KIND[i_lines_indep[ii]] eq 'l') then begin
      zlayer=zlayer+1
      fxaddpar, hdr_cHDU,'GD_Z__'+string(linesetup.I[i_lines_indep[ii]],format='(I02)'), zlayer,'GANDALF line Z in the cubes'
    endif else begin
      fxaddpar, hdr_cHDU,'GD_Z__'+string(linesetup.I[i_lines_indep[ii]],format='(I02)'), -1,'GANDALF line Z in the cubes'
    endelse
    fxaddpar, hdr_cHDU,'GD_NAM'+string(linesetup.I[i_lines_indep[ii]],format='(I02)'), linesetup.NAME[i_lines_indep[ii]],'GANDALF line name'
    fxaddpar, hdr_cHDU,'GD_LAM'+string(linesetup.I[i_lines_indep[ii]],format='(I02)'), linesetup.LAMBDA[i_lines_indep[ii]],'GANDALF line lambda'
    fxaddpar, hdr_cHDU,'GD_KND'+string(linesetup.I[i_lines_indep[ii]],format='(I02)'), linesetup.KIND[i_lines_indep[ii]],'GANDALF line kind'
    fxaddpar, hdr_cHDU,'GD_A__'+string(linesetup.I[i_lines_indep[ii]],format='(I02)'), linesetup.A[i_lines_indep[ii]],'GANDALF line ratio'
    fxaddpar, hdr_cHDU,'GD_FIT'+string(linesetup.I[i_lines_indep[ii]],format='(I02)'), linesetup.FIT[i_lines_indep[ii]],'GANDALF line fit type'
  endfor
  fxaddpar, hdr_cHDU,'GD_EMSET', emission_setup_file,'Emission line setup file'

  mwrfits,best_nebular_det,out_nebular,hdr_cHDU,/create

  sxdelpar,hdr_cHDU,'SIMPLE'
  fxaddpar, hdr_cHDU, 'CONTENT','ALL Nebular emission from GANDALF'
  fxaddpar, hdr_cHDU,'EXTNAME', 'NEBULAR_ALL'
  mwrfits,best_nebular_all,out_nebular,hdr_cHDU

  sxdelpar,hdr_cHDU,['NAXIS3','CRVAL3','CDELT3','CRPIX3','CUNIT3','DISPAXIS']
  fxaddpar, hdr_cHDU, 'CONTENT','Emission line flux'
  fxaddpar, hdr_cHDU,'EXTNAME', 'EM_LINE_FLUX'
  mwrfits,em_flux,out_nebular,hdr_cHDU

  fxaddpar, hdr_cHDU, 'CONTENT','Emission line amplitude'
  fxaddpar, hdr_cHDU,'EXTNAME', 'EM_LINE_AMPL'
  mwrfits,em_ampl,out_nebular,hdr_cHDU

  fxaddpar, hdr_cHDU, 'CONTENT','Emission line Vsys'
  fxaddpar, hdr_cHDU,'EXTNAME', 'EM_LINE_VSYS'
  mwrfits,em_vsys,out_nebular,hdr_cHDU

  fxaddpar, hdr_cHDU, 'CONTENT','Emission line Vdisp'
  fxaddpar, hdr_cHDU,'EXTNAME', 'EM_LINE_VDISP'
  mwrfits,em_disp,out_nebular,hdr_cHDU

  fxaddpar, hdr_cHDU, 'CONTENT','Emission line flux error'
  fxaddpar, hdr_cHDU,'EXTNAME', 'EM_LINE_FLUX_ERR'
  mwrfits,em_flux_err,out_nebular,hdr_cHDU

  fxaddpar, hdr_cHDU, 'CONTENT','Emission line amplitude error'
  fxaddpar, hdr_cHDU,'EXTNAME', 'EM_LINE_AMPL_ERR'
  mwrfits,em_ampl_err,out_nebular,hdr_cHDU

  fxaddpar, hdr_cHDU, 'CONTENT','Emission line Vsys error'
  fxaddpar, hdr_cHDU,'EXTNAME', 'EM_LINE_VSYS_ERR'
  mwrfits,em_vsys_err,out_nebular,hdr_cHDU

  fxaddpar, hdr_cHDU, 'CONTENT','Emission line Vdisp error'
  fxaddpar, hdr_cHDU,'EXTNAME', 'EM_LINE_VDISP_ERR'
  mwrfits,em_disp_err,out_nebular,hdr_cHDU

  fxaddpar, hdr_cHDU, 'CONTENT','Emission line Amplitude/Noise'
  fxaddpar, hdr_cHDU,'EXTNAME', 'EM_LINE_AoN'
  mwrfits,em_AoN,out_nebular,hdr_cHDU

  stop
end


pro check_GANDALF_LDRG
  common paths

  ;nameroot='nodding_macs2129_full_nir_9.6.main_box'
  ;nameroot='MACS2129-2_nir_master_stack_opt_3sigma_5iter_allcor.rebin16.opt.centerforgrad_box.1d'
  ;nameroot='MACS2129-2_nir_master_stack_opt_3sigma_5iter_allcor.rebin16.opt.outerforgrad_box.1d'
  ;nameroot='MACS2129-2_nir_master_stack_opt_3sigma_5iter_allcor.rebin16.opt.lowerforgrad_box.1d'
  ;nameroot='MACS2129-2_nir_master_stack_opt_3sigma_5iter_allcor.rebin16.opt.upperforgrad_box.1d'
  ;nameroot='nodding_macs2129_full_nir_9.6.main_box_2'
  ;nameroot='108899_avwCombined_OBx5_sig5_V1_NIRcorr_wmrebin15_opt_bpm'
  nameroot='105842_V1_NIRcorr_wmrebin15_er2d_opt_bpm'
  in_stellar=nameroot+'_GNDF_stellar.fits'
  if ~file_test(in_stellar) then begin
    in_stellar=in_stellar+'.gz'
    if ~file_test(in_stellar) then message,'_GNDF_stellar.fits[.gz] not found'
  endif
  in_nebular=nameroot+'_GNDF_nebular.fits'
  if ~file_test(in_nebular) then begin
    in_nebular=in_nebular+'.gz'
    if ~file_test(in_nebular) then message,'_GNDF_nebular.fits[.gz] not found'
  endif
  in_original=nameroot+'.fits'
  if ~file_test(in_original) then begin
    in_original=in_original+'.gz'
    if ~file_test(in_original) then message,prefix+'.fits[.gz] not found'
  endif

  ;; read in cubes
  ;; original
  original_cube=mrdfits(in_original,0,original_hdr)
  lambda0=sxpar(original_hdr,'CRVAL1')
  dlambda=sxpar(original_hdr,'CDELT1')
  nlambda=sxpar(original_hdr,'NAXIS1')
  lambda=lambda0+findgen(nlambda)*dlambda

  originalnoise_cube=mrdfits(in_original,'ERROR',originalnoise_hdr)
  lambda0=sxpar(originalnoise_hdr,'CRVAL1')
  dlambda=sxpar(originalnoise_hdr,'CDELT1')
  nlambda=sxpar(originalnoise_hdr,'NAXIS1')
  lambda=lambda0+findgen(nlambda)*dlambda

  originalbp_cube=mrdfits(in_original,'BADPIX',originalbp_hdr)
  lambda0=sxpar(originalbp_hdr,'CRVAL1')
  dlambda=sxpar(originalbp_hdr,'CDELT1')
  nlambda=sxpar(originalbp_hdr,'NAXIS1')
  lambdabp=lambda0+findgen(nlambda)*dlambda

  originalSNR_ima=0 ;mrdfits(in_original,'SNR_MAP',originalSNR_hdr)

  ;; stellar cube
  netstellar_cube=mrdfits(in_stellar,0,netstellar_hdr)
  lnlambda0=sxpar(netstellar_hdr,'CRVAL1')
  dlnlambda=sxpar(netstellar_hdr,'CDELT1')
  nlnlambda=sxpar(netstellar_hdr,'NAXIS1')
  z_sys=sxpar(netstellar_hdr,'Z_SYS')
  lambda_netstellar=exp(lnlambda0+findgen(nlnlambda)*dlnlambda)*(1.+z_sys)

  beststellar_cube=mrdfits(in_stellar,'BEST_STELLAR',beststellar_hdr)
  lnlambda0=sxpar(beststellar_hdr,'CRVAL1')
  dlnlambda=sxpar(beststellar_hdr,'CDELT1')
  nlnlambda=sxpar(beststellar_hdr,'NAXIS1')
  z_sys=sxpar(netstellar_hdr,'Z_SYS')
  lambda_beststellar=exp(lnlambda0+findgen(nlnlambda)*dlnlambda)*(1.+z_sys)

  vsys_ppxf=mrdfits(in_stellar,'VSYS_PPXF')
  vdisp_ppxf=mrdfits(in_stellar,'VDISP_PPXF')
  ;gndf_flags=mrdfits(in_stellar,'GNDF_FLAGS')
  gndf_rms=mrdfits(in_stellar,'GNDF_RMS')
  chi2dof_ppxf=mrdfits(in_stellar,'CHI2DOF_PPXF')

  if (fxposit(in_stellar,'NOISE_SCALE_FACTOR',errmsg=errmsg) gt 0 and keyword_set(noise_resc)) then begin
    noiserescmap=mrdfits(in_stellar,'NOISE_SCALE_FACTOR',hdr_noiseresc)
    for i=0,nx-1 do begin
      for j=0,ny-1 do begin
        if noiserescmap[i,j] gt 0. then originalnoise_cube[i,j,*]=originalnoise_cube[i,j,*]*noiserescmap[i,j]
      endfor
    endfor
    torescale=where(noiserescmap gt 0.)
    originalSNR_ima[torescale]=originalSNR_ima[torescale]/noiserescmap[torescale]
    chi2dof_ppxf[torescale]=chi2dof_ppxf[torescale]/(noiserescmap[torescale]^2)
  endif

  ;; nebular cube
  bestnebular_cube=mrdfits(in_nebular,0,bestnebular_hdr)
  lnlambda0=sxpar(bestnebular_hdr,'CRVAL1')
  dlnlambda=sxpar(bestnebular_hdr,'CDELT1')
  nlnlambda=sxpar(bestnebular_hdr,'NAXIS1')
  lambda_bestnebular=exp(lnlambda0+findgen(nlnlambda)*dlnlambda)*(1.+z_sys)

  allnebular_cube=mrdfits(in_nebular,'NEBULAR_ALL',allnebular_hdr)
  lnlambda0=sxpar(allnebular_hdr,'CRVAL1')
  dlnlambda=sxpar(allnebular_hdr,'CDELT1')
  nlnlambda=sxpar(allnebular_hdr,'NAXIS1')
  lambda_allnebular=exp(lnlambda0+findgen(nlnlambda)*dlnlambda)*(1.+z_sys)


  ;; mask skylines using the GANDALF procedures
  
  emission_setup_file=sxpar(beststellar_hdr,'GD_EMSET')
  ;emission_setup_file=inst_dir+'emission_lines_setup_std_fullspec_freeall.txt' ;;'emission_lines_setup_std_fullspec_Ha_kin.txt' ;;'/Users/zibetti/idl/CALIFA/emission_lines_setup_std.txt'
  eml_file = emission_setup_file
  readcol,eml_file,eml_i,eml_name,eml_lambda,eml_action,eml_kind,eml_a,eml_v,eml_s,eml_fit,$
    f='(i,a,f,a,a,f,f,f,a)',skipline=2,comment='#',/silent
  emission_setup = create_struct('i',eml_i,'name',eml_name,'lambda',eml_lambda,'action',eml_action,$
    'kind',eml_kind,'a',eml_a,'v',eml_v,'s',eml_s,'fit',eml_fit)
  ;message,'check sky mask!'
  linterp,lambda_allnebular,obslinemask(emission_setup,lnlambda0+findgen(nlnlambda)*dlnlambda,z_sys),lambda,sky_mask

  original_spec=reform(original_cube[*],nlambda)
  originalnoise_spec=reform(originalnoise_cube[*],nlambda)
  original_bpm=reform(originalbp_cube[*],nlambda)
  okdata=where(original_bpm eq 0)
  linterp,lambda_netstellar,reform(netstellar_cube[*],nlnlambda),lambda,netstellar_spec
  linterp,lambda_beststellar,reform(beststellar_cube[*],nlnlambda),lambda,beststellar_spec
  linterp,lambda_bestnebular,reform(bestnebular_cube[*],nlnlambda),lambda,bestnebular_spec
  linterp,lambda_allnebular,reform(allnebular_cube[*],nlnlambda),lambda,allnebular_spec

  xrange=[3700.,7000.]*(1.+z_sys)
  yrange=fltarr(2)
  yrange[0]=-4.*median(originalnoise_spec[okdata])
  yrange[1]=-yrange[0]+max(original_spec[okdata])
  win=window(dimension=[800,800])
  spec_plot=plot(lambda,original_spec,xrange=xrange,yrange=yrange,xstyle=1,ystyle=1,color='Black',thick=1,$
    ;xtitle='$\lambda_{obs}$/$\AA$',$;textoidl('\lambda_{obs}/\AA')
    ytitle='Flux/[$erg cm^{-2} s^{-1} \AA^{-1}$]',font_name='Times',position=[0.15,0.45,0.9,0.9],$
    /nodata,/current)
    
  tit_lab=text(0.5,0.965,nameroot+' - z='+string(z_sys,format='(F5.3)'),/norm,ALIGNMENT=0.5)
  yrange=spec_plot.yrange
  xx=0
  masks=[]
  while xx lt nlambda do begin
    x0=xx
    while (xx lt nlambda-1 and (original_bpm[xx] ne 0 or sky_mask[xx] gt 0)) do xx=xx+1
    if xx gt x0 then begin
      mask=polygon([lambda[x0],lambda[x0],lambda[xx],lambda[xx]],[yrange[0],yrange[1],yrange[1],yrange[0]],/FILL_BACKGROUND, $
        FILL_COLOR='light steel blue',COLOR='light steel blue',/overplot,/data)
      masks=[masks,mask]
    endif
    xx=xx+1
  endwhile
  spec_plot2=plot(lambda,original_spec,color='Black',thick=1,/overplot)
  bestfit_plot=plot(lambda,beststellar_spec+bestnebular_spec,color='orchid',thick=3,/overplot)
  beststellar_plot=plot(lambda,beststellar_spec,color='Lime Green',thick=2,/overplot)
  allnebular_plot=plot(lambda,allnebular_spec,color='dark red',thick=2,/overplot)
  bestnebular_plot=plot(lambda,bestnebular_spec,color='dodger blue',thick=2,/overplot)
  xspan=(spec_plot.xrange)[1]-(spec_plot.xrange)[0]
  stellar_kin=text((spec_plot.xrange)[0]+0.25*xspan,(spec_plot.yrange)[1]*0.90,'$V_{sys,*}$='+$
    string(vsys_ppxf[0],format='(F+6.1)')+'  $\sigma_*$='+string(vdisp_ppxf[0],format='(F5.1)'),/data,/overplot,font_name='Times',font_size=10)

;  _em_line_vsys=mrdfits(in_nebular,'EM_LINE_VSYS')
;  em_line_vsys=_em_line_vsys[21-1] ;; Halpha
;  _em_line_vdisp=mrdfits(in_nebular,'EM_LINE_VDISP')
;  em_line_vdisp=_em_line_vdisp[1-1] ;; Halpha
;
;  nebular_kin=text((spec_plot.xrange)[0]+0.25*xspan,(spec_plot.yrange)[1]*0.81,'$V_{sys,em.lines}$='+$
;    string(em_line_vsys[i,j],format='(F+6.1)')+'  $\sigma_{em.lines}$='+string(em_line_vdisp[i,j],format='(F5.1)'),/data,/overplot,font_name='Times',font_size=10)

    sp_axes=spec_plot.axes
  sp_axes[2].hide=1
  CSFR_axis=axis('X',location='top',$
    COORD_TRANSFORM=[0,1./(1.+z_sys)],$
    title='$\lambda_{rf}/\AA$',target=spec_plot)

  yrange_res=[yrange[0],-yrange[0]]
  res_plot=plot(lambda,beststellar_spec+bestnebular_spec-original_spec,$
    xrange=xrange,yrange=yrange_res,xstyle=1,ystyle=1,position=[0.15,0.10,0.9,0.37],color='Black',thick=2,symbol="dot",sym_thick=2,linestyle="none",$
    xtitle='$\lambda_{obs}$/$\AA$',$
    ;ytitle='$\Delta Flux/[10^{-16} erg cm^{-2} s^{-1} \AA^{-1}$]',$
    ytitle='$\Delta$Flux',$
    font_name='Times',/current,$
    title='rms='+string(gndf_rms[0],format='(E7.1)')+'   $\chi^2_{dof}$='+string(chi2dof_ppxf[0],format='(F4.1)'))

    rms_plot=polygon([lambda,reverse(lambda)],[originalnoise_spec,-reverse(originalnoise_spec)],replicate(-0.1,n_elements(lambda)*2),$
    color='coral',fill_color='coral',/data,/overplot,target=res_plot)
  xx=0
  masks_res=[]
  while xx lt nlambda do begin
    x0=xx
    while (xx lt nlambda-1 and (original_bpm[xx] ne 0 or sky_mask[xx] gt 0)) do xx=xx+1
    if xx gt x0 then begin
      mask=polygon([lambda[x0],lambda[x0],lambda[xx],lambda[xx]],[yrange[0],yrange[1],yrange[1],yrange[0]],/FILL_BACKGROUND, $
        FILL_COLOR='light steel blue',COLOR='light steel blue',/overplot,/data,target=res_plot)
      masks_res=[masks_res,mask]
    endif
    xx=xx+1
  endwhile
  stop
end
