pro show_Gandalf_fit,in_file,ps=ps,ylog=ylog,outfile=outfile
; Simply plots the GANDALF fit to the SDSS spectra, 
; allowing to zoom in two specified area of interest
; M. Sarzi, Hatfield, December 2007
; 
; show_fit,'spSpec-53149-1421-349.fit'     --> will plot on screen
; show_fit,'spSpec-53149-1421-349.fit',/ps --> will make default.ps
; show_fit,'spSpec-53149-1421-349.fit',/ps,outfile='myplot.ps' --> will make myplot.ps
; show_fit,'spSpec-53149-1421-349.fit',/ylog --> use logarithmic y-scale

sauron_colormap

fits_file = (strsplit(in_file,'.',/extract))[0]+'_fits.fits'
pars_file = (strsplit(in_file,'.',/extract))[0]+'_pars.fits'

; read in object and model spectra
readfits_spec,fits_file,spec=obj,l=lgl,hdr=hdr;,/verbose
fit = mrdfits(fits_file,1,/silent)
emi = mrdfits(fits_file,2,/silent)
cfit = fit-emi
res = obj-fit
; read in also the goodpixels array
goodpixels = mrdfits(fits_file,4,/silent)
; get rest-frame wavelength array
pars = mrdfits(pars_file,1,/silent)
c = 299792.458d
l = 10^lgl
l_rf = l/exp(pars.vel_stars/c)

; lambda plotting range
mn_l = min(l,max=mx_l)
lrange = mx_l - mn_l & lrange = [mn_l-0.05*lrange,mx_l+0.05*lrange]
lrange_rf_0 = lrange/(1+pars.vel_stars/c)
; lambda plotting range zoom 1
;lrange_rf_1 = [3700.1,4250.0]
;region_1='[OII]+[NeIII]+H!4d!6 region'
; lambda plotting range zoom 2
lrange_rf_1 = [4750.1,5300.0]
region_1='H!4b!6+[OIII] region'
;
lrange_rf_2 = [4750.1,5300.0]+870
region_2='NaD region'
;
; lambda plotting range zoom 3
lrange_rf_3 = [6250.0,6799.9]
region_3='H!4a!6+[NII] region'
; flux plotting range
ron = robust_sigma(res[where(goodpixels gt 0)],/zero) 
mx_f = 2*(max(cfit)+5*ron)
mn_f = -5*ron
if keyword_set(ylog) then begin 
    mx_f = max(fit)
    mn_f = min(cfit)-ron
    mx_f = mx_f * 10^(0.05*alog10(mx_f/mn_f)) 
endif
frange = [mn_f,mx_f]
; to plot goodpixels
if not keyword_set(ylog) then begin 
    goodpixels[where(goodpixels eq 0)] = -1
    goodpixels = goodpixels*1e6
endif else begin
    goodpixels = goodpixels*1e6
    goodpixels[where(goodpixels eq 0)] = 1e-6
endelse
; plot positions
pos0=[0.070,0.59,0.970,0.98] ; main
pos1=[0.070,0.10,0.337,0.49] ; zoom 1
pos2=[0.387,0.10,0.653,0.49] ; zoom 2
pos3=[0.703,0.10,0.970,0.49] ; zoom 3


; open PS file if necessary
if not keyword_set(ps) then begin
   device,decomp=0 
   window,0,xs=600*1.5,ys=400*1.5
endif else begin 
   if keyword_set(outfile) then openps,file=outfile,xs=29*0.9,ys=21*0.9,/landscape,yoff=29*0.95,xoff=0 else openps,xs=29*0.9,ys=21*0.9,/landscape,yoff=29*0.95,xoff=0
   !p.charsize = 0.9
   !x.charsize = 0.9
   !y.charsize = 0.9
   !p.charthick= 2
   !p.thick= 2
   !x.thick= 2
   !y.thick= 2
endelse
; plot entire wavelength range
xtitle='!4k!6!Dr.f.!N'
ytitle='flux density'
 plot,l_rf,obj,psym=10,xr=lrange_rf_0,/xs,yr=frange,/ys,xtitle=xtitle,ytitle=ytitle,pos=pos0,ylog=ylog
oplot,l_rf,cfit,col=125,lines=3
oplot,l_rf,fit,col=210,thick=3
oplot,l_rf,emi,col=80
oplot,l_rf,res,psym=3
oplot,l_rf,l_rf*0
oplot,l_rf,l_rf*0 + ron,lines=2
loadct,0,/silent & oplot,l_rf,goodpixels,lines=1,col=150 & sauron_colormap
legend,[fits_file],/top,/left,box=0,pos=[pos0[0],pos0[3]],/normal
if n_elements(pars.EBmV) eq 2 then begin
   str1='E(B-V)!Dscreen!N='+string( strmid(strcompress(round((pars.EBmV)[0]*1000.)/1000.,/remove_all),8,6,/reverse))
   str2='E(B-V)!Dintern!N='+string( strmid(strcompress(round((pars.EBmV)[1]*1000.)/1000.,/remove_all),8,6,/reverse))
   legend,[str1,str2],/top,/right,box=0,pos=[pos0[2],pos0[3]],/normal
endif 
if n_elements(pars.EBmV) eq 1 then begin
   str1='E(B-V)!Dscreen!N='+string( strmid(strcompress(round((pars.EBmV)[0]*1000.)/1000.,/remove_all),8,6,/reverse))
   legend,[str1],/top,/right,box=0,pos=[pos0[2],pos0[3]],/normal
endif 
; zoom region 1
ytitle=''
 plot,l_rf,obj,psym=10,xr=lrange_rf_1,/xs,yr=frange,/ys,xtitle=xtitle,ytitle=ytitle,pos=pos1,ylog=ylog,/noerase
oplot,l_rf,cfit,col=125,lines=3
oplot,l_rf,fit,col=210,thick=3
oplot,l_rf,emi,col=80
oplot,l_rf,res,psym=3
oplot,l_rf,l_rf*0
oplot,l_rf,l_rf*0 + ron,lines=2
loadct,0,/silent & oplot,l_rf,goodpixels,lines=1,col=150 & sauron_colormap
legend,[region_1],/top,/left,box=0,pos=[pos1[0],pos1[3]],/normal
; zoom region 2
 plot,l_rf,obj,psym=10,xr=lrange_rf_2,/xs,yr=frange,/ys,xtitle=xtitle,ytitle=ytitle,pos=pos2,ylog=ylog,/noerase
oplot,l_rf,cfit,col=125,lines=3
oplot,l_rf,fit,col=210,thick=3
oplot,l_rf,emi,col=80
oplot,l_rf,res,psym=3
oplot,l_rf,l_rf*0
oplot,l_rf,l_rf*0 + ron,lines=2
loadct,0,/silent & oplot,l_rf,goodpixels,lines=1,col=150 & sauron_colormap
legend,[region_2],/top,/left,box=0,pos=[pos2[0],pos2[3]],/normal
; zoom region 2
 plot,l_rf,obj,psym=10,xr=lrange_rf_3,/xs,yr=frange,/ys,xtitle=xtitle,ytitle=ytitle,pos=pos3,ylog=ylog,/noerase
oplot,l_rf,cfit,col=125,lines=3
oplot,l_rf,fit,col=210,thick=3
oplot,l_rf,emi,col=80
oplot,l_rf,res,psym=3
oplot,l_rf,l_rf*0
oplot,l_rf,l_rf*0 + ron,lines=2
loadct,0,/silent & oplot,l_rf,goodpixels,lines=1,col=150 & sauron_colormap
legend,[region_3],/top,/left,box=0,pos=[pos3[0],pos3[3]],/normal
; close PS file if necessary
if keyword_set(ps) then begin
    !p.charsize = 1
    !x.charsize = 1
    !y.charsize = 1
    !p.charthick= 1
    !p.thick= 1
    !x.thick= 1
    !y.thick= 1    
    closeps
 endif else stop

end

pro openps,FILE=FILE,XSIZE=XSIZE,YSIZE=YSIZE,XOFFSET=XOFFSET,YOFFSET=YOFFSET,_extra=extra

if not keyword_set(FILE)    then FILE = 'default.ps'
if not keyword_set(XSIZE)   then xsize = 14 
if not keyword_set(YSIZE)   then ysize = 14
if not keyword_set(XOFFSET) then xoffset = 1
if not keyword_set(YOFFSET) then yoffset = 1
set_plot,'ps'
device,filename=file,xsize=xsize,ysize=ysize,xoffset=xoffset,yoffset=yoffset,/color,bits=8,/cmyk,_extra=extra
end

pro closeps
device,/close
CASE !VERSION.OS_FAMILY OF
    'Windows': SET_PLOT,'WIN'
    'MacOS': SET_PLOT,'MAC'
    ELSE: SET_PLOT,'X'
ENDCASE
end
