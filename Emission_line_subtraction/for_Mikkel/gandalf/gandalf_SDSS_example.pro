PRO GANDALF_SDSS_EXAMPLE, FOR_ERRORS=for_errors

; Example of 'standard' measurement of the flux and kinematics of the
; emission lines in the SDSS spectra, following the approach of
; Tremonti et al. (2004) to tie the recombination and forbidden lines
; to two separate kinematics. We make use also of the same set of
; Bruzual & Charlot (2003) models as stellar templates. Finally we
; include reddening by dust, including two components: a first diffuse
; component affecting both the stellar and emission-line fluxes, and a
; second that is embedded in the emission-line regions and affecting
; only the line fluxes.

; In the following examples I did not always include the foreground
; Galactic extiction. If not, this will be factored in the diffuse
; screen component (which is 99% ok for low redshift).

GANDALF_SDSS,INFILE='spSpec-51908-0277-120.fit', $
  EMISSION_SETUP_FILE='emission_lines_setup_with_Balmer_decrement',$
  LIBFILE='BC03_Tremonti04_templates.dat', MDEGREE=3, /PLOT, FOR_ERRORS=for_errors

GANDALF_SDSS,INFILE='spSpec-51671-0299-123.fit', $
  EMISSION_SETUP_FILE='emission_lines_setup_with_Balmer_decrement',$
  LIBFILE='BC03_Tremonti04_templates.dat', MDEGREE=3, /PLOT, FOR_ERRORS=for_errors

GANDALF_SDSS,INFILE='spSpec-51985-0295-530.fit', $
  EMISSION_SETUP_FILE='emission_lines_setup_with_Balmer_decrement',$
  LIBFILE='BC03_Tremonti04_templates.dat', MDEGREE=3, /PLOT, FOR_ERRORS=for_errors

; Here foreground galactic extinction is also accounted for, by
; passing the Schlegel et al. E(B-V) reddening value that will be used
; to de-redden the spectra prior to the kinematic and emission-line
; analysis.

GANDALF_SDSS,INFILE='spSpec-53149-1421-349.fit', $
  EMISSION_SETUP_FILE='emission_lines_setup_with_Balmer_decrement',$
  LIBFILE='BC03_Tremonti04_templates.dat', MDEGREE=3, /PLOT, FOR_ERRORS=for_errors, EBV_GAL=0.044

; Here is an object with a BLR, which is dealt with by simply adding a
; couple of initially braoder Gaussian Hb and Ha components. Fit could
; still be improved but this shows just a good start

GANDALF_SDSS,INFILE='spSpec-53141-1417-463.fit', $
  EMISSION_SETUP_FILE='emission_lines_setup_with_Balmer_decrement_BLR',$
  LIBFILE='BC03_Tremonti04_templates.dat', MDEGREE=3, /PLOT, FOR_ERRORS=for_errors


END
@gandalf_SDSS.pro
@ppxf.pro
@gandalf.pro
@bvls.pro
@range.pro
@sauron_colormap.pro
@show_Gandalf_fit
@get_all_lines_amplitudes_fluxes_and_EWs

PRO SHOW_GANDALF_FIT_ALL

show_Gandalf_fit,'spSpec-51908-0277-120.fit'
show_Gandalf_fit,'spSpec-51671-0299-123.fit'
show_Gandalf_fit,'spSpec-51985-0295-530.fit'
show_Gandalf_fit,'spSpec-53149-1421-349.fit'
show_Gandalf_fit,'spSpec-53141-1417-463.fit'

END

PRO MAKEPS_ALL

show_Gandalf_fit,'spSpec-51908-0277-120.fit',/ps,outfile='spSpec-51908-0277-120_Gandalf.ps'
show_Gandalf_fit,'spSpec-51671-0299-123.fit',/ps,outfile='spSpec-51671-0299-123_Gandalf.ps'
show_Gandalf_fit,'spSpec-51985-0295-530.fit',/ps,outfile='spSpec-51985-0295-530_Gandalf.ps'
show_Gandalf_fit,'spSpec-53149-1421-349.fit',/ps,outfile='spSpec-53149-1421-349_Gandalf.ps'
show_Gandalf_fit,'spSpec-53141-1417-463.fit',/ps,outfile='spSpec-53149-1421-349_Gandalf.ps'

END

PRO GET_AMPLITUDES_FLUXES_AND_EWS_ALL

get_all_lines_amplitudes_fluxes_and_EWs,'spSpec-51908-0277-120.fit','emission_lines_setup_with_Balmer_decrement'
get_all_lines_amplitudes_fluxes_and_EWs,'spSpec-51671-0299-123.fit','emission_lines_setup_with_Balmer_decrement'
get_all_lines_amplitudes_fluxes_and_EWs,'spSpec-51985-0295-530.fit','emission_lines_setup_with_Balmer_decrement'
get_all_lines_amplitudes_fluxes_and_EWs,'spSpec-53149-1421-349.fit','emission_lines_setup_with_Balmer_decrement'
get_all_lines_amplitudes_fluxes_and_EWs,'spSpec-53141-1417-463.fit','emission_lines_setup_with_Balmer_decrement_BLR'

END

