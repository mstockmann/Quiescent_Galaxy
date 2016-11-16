pro readfits_spec,fitsfile,spec=spec,l=l,hdr=hdr,verbose=verbose
spec = readfits(fitsfile,hdr,/SILENT)
s    = size(spec)
CRVAL1 = sxpar(hdr,'CRVAL1',/SILENT)
CRPIX1 = sxpar(hdr,'CRPIX1',/SILENT) & if CRPIX1 eq 0 then CRPIX1 = 1
CDELT1 = sxpar(hdr,'CDELT1',/SILENT) & if (CDELT1 eq 0) then CDELT1 = sxpar(hdr,'CD1_1',/SILENT)
NAXIS1 = sxpar(hdr,'NAXIS1',/SILENT)

if keyword_set(verbose) then begin
    if s[0] eq 1 then print,'1D-spectrum'
    if s[0] eq 2 then print,'2D-spectrum',s[1],'x',s[2],format='(a12,i5,a3,i5)'
endif
l    = CRVAL1 + CDELT1*(dindgen(NAXIS1) - CRPIX1 +1)
if keyword_set(verbose) then print,l[fix(CRPIX1)-1]+(CRPIX1-fix(CRPIX1))*CDELT1,'should be equal to',CRVAL1,format='(f10.4,a19,f10.4)'

end
