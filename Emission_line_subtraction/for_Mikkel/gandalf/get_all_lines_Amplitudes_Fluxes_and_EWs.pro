FUNCTION DUST_CALZETTI_lambda,ebv,lambda
; This procedure uses the dust model of Calzetti et al. (2000, ApJ,
; 533, 682), and for a given E(B-V) value returns the flux attenuation
; array, which can be used to get reddened templates. Here the spectra
; are assumed to be binned on a ln-rebinned wavelentgh grid as defined
; by input l0_gal,lstep_gal,npix parameters. The input receiding
; velocity vstar, is used to derive the dust reddening in the galaxy
; rest-frame.
; 
; Can be used also to de-reddened the object spectra by the Milky-Way
; dust extinction, using as E(B-V) the opposite of the Schlegel et
; al. values found in NED and vstar = 0.
;
; Initial version kindly provided by S. Kaviray, Oxford, 2006.

; array to hold k(lambda) values
k = fltarr(n_elements(lambda))           

for i=0,n_elements(lambda)-1 do begin
     ; convert wavelength units from angstroms to micrometres
     l = lambda(i)/1e4                   
     ; assign k values
     if (l ge 0.63 and l le 2.2) then k(i) = 2.659*(-1.857+1.040/l)+4.05
     if (l lt 0.63)              then k(i) = 2.659*(-2.156+1.509/l-0.198/l^2+0.011/l^3)+4.05
     if (l gt 2.2)               then k(i) = 0.0
endfor

return,(10^(-0.4*ebv*k))
; this should be then multiplied by the spectrum flux array
END

pro GET_ALL_LINES_AMPLITUDES_FLUXES_AND_EWS,in_file,emission_setup_file

; Given the input .fit SDSS spectrum, read in the Gandalf output
; _fits.fits and _pars.fits files to rederive the values for the
; amplitude and fluxes of all lines, both as observed and dereddened
; by dust extinction. Also used the best model for the continuum fit
; and the wavelength array to compute the lines EW values.
; 

; read in the output _fits.fits file corresponding to the input SDSS
; spectrum. 
;
fits_file = (strsplit(in_file,'.',/extract))[0]+'_fits.fits'
readfits_spec,fits_file,spec=obj,l=lgl,hdr=hdr
bestfit       = mrdfits(fits_file,1,/silent)
emission      = mrdfits(fits_file,2,/silent)
goodpixels    = mrdfits(fits_file,4,/silent)
continuum_fit = bestfit-emission

; read in the output _pars.fits file corresponding to the input SDSS
; spectrum. This will give us the outputs of pPXF and GANDALF
;
pars_file = (strsplit(in_file,'.',/extract))[0]+'_pars.fits'
fit_results = mrdfits(pars_file,1,/silent)

; observed and rest-frame wavelength array, to be used later
;
c = 299792.458d
obs_l    = 10^lgl
rf_l  = obs_l/exp(fit_results.Vel_stars/c)

; read in the emission-line setup file and create the initial emission-line
; structure. Then manipulate it exactly as done in the SDSS wrapper
; to get the emission-line structure as passed to GANDALF (except for
; the v_ and s_guess which we do not care about)
;
eml_file = emission_setup_file
readcol,eml_file,eml_i,eml_name,eml_lambda,eml_action,eml_kind,eml_a,eml_v,eml_s,eml_fit,$
  f='(i,a,f,a,a,f,f,f,a)',skipline=2,comment='#',/silent
emission_setup = create_struct('i',eml_i,'name',eml_name,'lambda',eml_lambda,'action',eml_action,$
                               'kind',eml_kind,'a',eml_a,'v',eml_v,'s',eml_s,'fit',eml_fit)
; emission-line structure manipulation
i = where(emission_setup.action eq 'm' and emission_setup.name ne 'sky')
emission_setup.action[i] = 'f'

; create the output emission-line structure - containing all lines except skylines
;
dummy = emission_setup
emission_setup_output = create_struct('i',dummy.i[i],'name',dummy.name[i],$
                                      'lambda',dummy.lambda[i],'action',dummy.action[i],$
                                      'kind',dummy.kind[i],'a',dummy.a[i],$
                                      'v',dummy.v[i],'s',dummy.s[i],$
                                      'fit',dummy.fit[i])

; create the same emission-line structure that was passed to Gandalf
; and which contains only the lines that were fitted for this
; particular objects, given its redshift.
;
l_rf1 = min(rf_l[where(goodpixels gt 0)])
l_rf2 = max(rf_l[where(goodpixels gt 0)])
; for sigma_mask = 200km/s 
i = where(emission_setup.action eq 'f' and $
          emission_setup.lambda*(1.0d + 200.0/c) ge l_rf1 and $
          emission_setup.lambda*(1.0d - 200.0/c)) 
dummy = emission_setup
emission_setup = create_struct('i',dummy.i[i],'name',dummy.name[i],$
                               'lambda',dummy.lambda[i],'action',dummy.action[i],$
                               'kind',dummy.kind[i],'a',dummy.a[i],$
                               'v',dummy.v[i],'s',dummy.s[i],$
                               'fit',dummy.fit[i])


; Note on indeces:

; i_... - indeces refers to the index of a line in the emission-line
;         structure, going from 0 to n-1
; j_.. - indeces refers to the index of a line in the emission-line
;        setup, which is effectively a label and can have any
;        value. These are the indeces that are used to tie together
;        the lines and the elements of line doublets and multiplets
;
; Note on reddening:
; If you have NOT used dust reddening to adjust the continuum fit, and
; insted used multiplicative polynomials, then uncomment the line 152
; To be made more automatic...

; separate main and satellite lines
;
i_ml       = where(emission_setup.kind eq 'l',complement=i_sl)

; index j of all the lines that we have fitted.
;
j = emission_setup.i

; index j of the the main line in each multiplet to which the satellite
; lines refer to. Will be used later to reassing intrinsic amplitudes
; and fluxes
;
j_ml_sl = fix(strmid(emission_setup.kind[i_sl],1))


n             = n_elements(emission_setup.i)
sol_gas_V     = dblarr(n)
sol_gas_S     = dblarr(n)
sol_gas_A_int = dblarr(n)
sol_gas_F_int = dblarr(n)
sol_gas_A_obs = dblarr(n)
sol_gas_F_obs = dblarr(n)
sol_gas_EW    = dblarr(n)

; compute the reddening attenuation, at the location of all the lines
; we have fitted, for both diffuse and nebular components.
;
reddening_attenuation = dust_calzetti_lambda(fit_results.EBmV[0],emission_setup.lambda)
if n_elements(fit_results.EBmV) eq 2 then begin
   int_reddening_attenuation = dust_calzetti_lambda(fit_results.EBmV[1],emission_setup.lambda)
   reddening_attenuation_emission = reddening_attenuation*int_reddening_attenuation 
endif else reddening_attenuation_emission = reddening_attenuation

; If we used multiplicative polynomials instead of reddening,
; uncomment this line
;
; reddening_attenuation_emission = reddening_attenuation_emission*0.0 + 1.0d

; assign the values of Amplitude and Fluxes (observed and intrinsic)
; as returned by Gandalf
;
sol_gas_A_obs[i_ml] = fit_results.Ampl
sol_gas_F_int[i_ml] = fit_results.Flux

; restore intrinsic amplitudes, for the lines we have fitted and only
; for the main members of doublets/multiplets
;
sol_gas_A_int[i_ml] = fit_results.Ampl/reddening_attenuation_emission[i_ml]

; compute the observed fluxes, for the same lines 
;
sol_gas_F_obs[i_ml] = fit_results.Flux*reddening_attenuation_emission[i_ml]


; loop over the satellites lines and fill the intrisic values of the
; amplitudes and fluxes, by simply scaling the corresponding values
; for the main elements in the doublet multiplet
;
for k=0,n_elements(i_sl)-1 do begin
   ; index i of the main line to which this i_sl[k] refers to
   i_ml_sl = where(emission_setup.i eq j_ml_sl[k])
   ; fixed! affects on the amplitudes
   sol_gas_A_int[i_sl[k]] = sol_gas_A_int[i_ml_sl]*emission_setup.a[i_sl[k]] *(emission_setup.lambda[i_ml_sl]/emission_setup.lambda[i_sl[k]])
   sol_gas_F_int[i_sl[k]] = sol_gas_F_int[i_ml_sl]*emission_setup.a[i_sl[k]] 
endfor

; now simply apply the reddening attenuation to them
;
sol_gas_A_obs[i_sl] = sol_gas_A_int[i_sl]*reddening_attenuation_emission[i_sl]
sol_gas_F_obs[i_sl] = sol_gas_F_int[i_sl]*reddening_attenuation_emission[i_sl]

; Finally repeat the same exercise for the line velocities and
; velocity dispersions. This assumes lines in doublet/multiplets share
; the same line profile
;
sol_gas_V[i_ml] = fit_results.Vel
sol_gas_S[i_ml] = fit_results.Sigma
for k=0,n_elements(i_sl)-1 do begin
   ; index i of the main line to which this i_sl[k] refers to
   i_ml_sl = where(emission_setup.i eq j_ml_sl[k])
   sol_gas_V[i_sl[k]] = sol_gas_V[i_ml_sl]
   sol_gas_S[i_sl[k]] = sol_gas_S[i_ml_sl]
endfor


; ok now get the EW of the lines
;
for i=0,n_elements(sol_gas_EW)-1 do begin
   ; get the continuum level from the stellar continuum fit in a +5 to
   ; +10 sigma buffer region around each line
   rf_l_line  = emission_setup.lambda[i]
   S_line     = sqrt(sol_gas_S[i]^2+61.5^2)
   F_obs_line = sol_gas_F_obs[i]
   j_buffer   = where(abs(rf_l-rf_l_line) lt 10*(S_line/c)*rf_l_line and abs(rf_l-rf_l_line) gt  5*(S_line/c)*rf_l_line)
   C_line     = median(continuum_fit[j_buffer])
   ; EW in \AA (Flux is in erg/s/cm^2,
   ; the continuum level is a flux
   ; spectral density in erg/s/cm^2/\AA)
   sol_gas_EW[i] = F_obs_line/C_line
endfor

; and the A/N of the lines
;
resid = obj - bestfit
resid_noise = robust_sigma(resid[where(goodpixels gt 0)],/zero)
sol_gas_AoN = sol_gas_A_obs/resid_noise

; finally, write out on screen all these outputs, including zeros for
; the lines that were in the original, input emission setup, but that
; were not fitted due their redshift
;
n                 = n_elements(emission_setup_output.i)
sol_gas_V_all     = dblarr(n)
sol_gas_S_all     = dblarr(n)
sol_gas_A_int_all = dblarr(n)
sol_gas_F_int_all = dblarr(n)
reddening_attenuation_emission_all = dblarr(n)
sol_gas_A_obs_all = dblarr(n)
sol_gas_F_obs_all = dblarr(n)
sol_gas_EW_all    = dblarr(n)
sol_gas_AoN_all   = dblarr(n)
;
i_fitted = where(emission_setup_output.action eq 'f' and $
          emission_setup_output.lambda*(1.0d + 200.0/c) ge l_rf1 and $
          emission_setup_output.lambda*(1.0d - 200.0/c)) 
;
sol_gas_V_all     [i_fitted] =  sol_gas_V     
sol_gas_S_all     [i_fitted] =  sol_gas_S     
sol_gas_A_int_all [i_fitted] =  sol_gas_A_int 
sol_gas_F_int_all [i_fitted] =  sol_gas_F_int 
reddening_attenuation_emission_all[i_fitted] = reddening_attenuation_emission
sol_gas_A_obs_all [i_fitted] =  sol_gas_A_obs 
sol_gas_F_obs_all [i_fitted] =  sol_gas_F_obs 
sol_gas_EW_all    [i_fitted] =  sol_gas_EW    
sol_gas_AoN_all   [i_fitted] =  sol_gas_AoN    
;
sol_gas_V     = sol_gas_V_all     
sol_gas_S     = sol_gas_S_all     
sol_gas_A_int = sol_gas_A_int_all 
sol_gas_F_int = sol_gas_F_int_all 
reddening_attenuation_emission = reddening_attenuation_emission_all
sol_gas_A_obs = sol_gas_A_obs_all 
sol_gas_F_obs = sol_gas_F_obs_all 
sol_gas_EW    = sol_gas_EW_all    
sol_gas_AoN   = sol_gas_AoN_all    

; finally print outputs on screen
print,['name','lambda','line type','rel. A','V','Sig','intr. A','intr. F','dust attenuation','obs. A','obs. F','EW','A/N'],$
      f='(8a10,a18,4a10)'
forprint,emission_setup_output.name,emission_setup_output.lambda,emission_setup_output.kind,emission_setup_output.a,$
         sol_gas_V,sol_gas_S,$
         sol_gas_A_int,sol_gas_F_int,$
         reddening_attenuation_emission,$
         sol_gas_A_obs,sol_gas_F_obs,$
         sol_gas_EW,sol_gas_AoN,text=2,f='(a10,f10.2,a10,5f10.2,f18.3,4f10.2)'
stop
end
