
function diffusion_scale,t,electron=electron,hole=hole
  ;[t]=nanoseconds, drift time
  if keyword_set(hole) then begin
    sigma=sqrt(2.*29.0*t) ;micrometers
  endif else begin
    sigma=sqrt(2.*25.0*t) ;micrometers
  endelse
  return,sigma
end

function repulsion_scale_extended,t,energy,R0,electron=electron,hole=hole
  ;[t]=nanoseconds, drift time
  ;[energy]=keV, interaction energy
  ;[R0]=um, radius of initial charge cloud
  ncharge=energy*1000./2.96 ;number of charge carriers
  if keyword_set(hole) then begin
    eta=(R0^3+ncharge*1.13*t)^(1./3.) ;micrometers
  endif else begin
    eta=(R0^3+ncharge*0.97*t)^(1./3.) ;micrometers ;micrometers
  endelse
  return,eta
end

function charge_distribution_extended,x,t,energy,R0,electron=electron,hole=hole
  ; convolution of diffusion and repulsion
  ;[R0]=um, radius of initial charge cloud
  if keyword_set(hole) then begin
    sigma=diffusion_scale(t,/hole) ;diffusion size scale, micrometers
    eta=repulsion_scale_extended(t,energy,R0,/hole) ;repulsion size scale, micrometers
  endif else begin
    sigma=diffusion_scale(t,/electron) ;diffusion size scale, micrometers
    eta=repulsion_scale_extended(t,energy,R0,/electron) ;repulsion size scale, micrometers
  endelse

  a=x-eta; lower limit of integration
  b=x+eta; upper limit of integration

  density = (0.75/(eta^3))*(0.5)*(eta^2-x^2-sigma^2)*(erf(b/(sqrt(2.)*sigma))-erf(a/(sqrt(2.)*sigma)))
  density = density+(0.75/(eta^3))*(sigma/sqrt(2.*!pi))*((2.*x-a)*exp(-0.5*a^2/sigma^2)-(2.*x-b)*exp(-0.5*b^2/sigma^2))

  return,density
end

PRO photopeak,x,a,f
  f=x
  f=a(0)*exp(-0.5*((x-a(1))/a(2))^2)
  return
end

PRO tailing_unconstrained,x,b,f
  f=x
  f=(b(0)*exp(b(3)*(x-b(1)))+b(4))*(1.0-erf((x-b(1))/(b(2)*sqrt(2.))))*(1.0+erf((x-b(1)+18.0)/(b(2)*sqrt(2.))))
  return
end

PRO tailing3,x,b,f
  f=x
  f=b(0)*(exp((0.47)*(x-b(1)))+(0.079))*(1.0-erf((x-b(1))/(b(2)*sqrt(2.))))*(1.0+erf((x-b(1)+18.0)/(b(2)*sqrt(2.))))
  return
end

PRO tailing_lab,x,b,f
  f=x
  f=(b(0)*exp(b(3)*(x-b(1)))+b(4))*(1.0-erf((x-b(1))/(b(2)*sqrt(2.))))*(1.0+erf((x-b(1)+18.0)/(b(2)*sqrt(2.))))
  return
end

PRO tailing_fixedpeak,x,b,f
COMMON PEAK
  f=x
  f=(b(0)*exp(b(1)*(x-photon_energy))+b(2))*(1.0-erf((x-photon_energy)/(sigma_energy*sqrt(2.))))*(1.0+erf((x-photon_energy+strip_thresh)/(sigma_energy*sqrt(2.))))
  return
end

PRO tailing_fixedpeak2,x,b,f
  COMMON PEAK
  f=x
  f=(b(0)*exp(b(1)*(x-photon_energy))+b(2))*(1.0-erf((x-photon_energy)/(b(3)*sqrt(2.))))*(1.0+erf((x-photon_energy+strip_thresh)/(b(3)*sqrt(2.))))
  return
end

PRO tailpeak_unconstrained,x,a,f
  COMMON PEAK
  f=x
  f=a(0)*exp(-0.5*((x-a(1))/a(2))^2)+(a(3)*exp(a(4)*(x-a(1)))+a(5))*(1.0-erf((x-a(1))/(a(6)*sqrt(2.))))*(1.0+erf((x-a(1)+strip_thresh)/(a(6)*sqrt(2.))))
  return
end

PRO tailing_unconstrained,x,a,f
  COMMON PEAK
  f=x
  f=(a(0)*exp(a(2)*(x-a(1)))+a(3))*(1.0-erf((x-a(1))/(a(4)*sqrt(2.))))*(1.0+erf((x-a(1)+strip_thresh)/(a(4)*sqrt(2.))))
  return
end

PRO tailpeak_unconstrained_v4,x,a,f ;short exponential + step 
  COMMON PEAK
  f=x
  f=a(0)*exp(-0.5*((x-a(1))/a(2))^2)+(a(3)*exp(a(4)*(x-a(1)))+a(5))*(1.0-erf((x-a(1))/(a(6)*sqrt(2.))))*(1.0+erf((x-a(1)+strip_thresh)/(a(6)*sqrt(2.))))
  return
end

PRO tailing_unconstrained_v4,x,a,f ;short exponential + step 
  COMMON PEAK
  f=x
  f=(a(0)*exp(a(2)*(x-a(1)))+a(3))*(1.0-erf((x-a(1))/(a(4)*sqrt(2.))))*(1.0+erf((x-a(1)+strip_thresh)/(a(4)*sqrt(2.))))
  return
end

PRO tailpeak_v4,x,a,f ;short exponential + step 
  COMMON PEAK
  f=x
  f=a(0)*exp(-0.5*((x-a(1))/a(2))^2)+a(3)*(exp((0.44)*(x-a(1)))+(0.075))*(1.0-erf((x-a(1))/((0.87)*a(2)*sqrt(2.))))*(1.0+erf((x-a(1)+strip_thresh)/((0.87)*a(2)*sqrt(2.))))
  ; f=a(0)*exp(-0.5*((x-a(1))/a(2))^2)+a(3)*(exp((0.48)*(x-a(1)))+(0.098))*(1.0-erf((x-a(1))/((0.88)*a(2)*sqrt(2.))))*(1.0+erf((x-a(1)+strip_thresh)/((0.88)*a(2)*sqrt(2.)))) ;60 keV
  ; f=a(0)*exp(-0.5*((x-a(1))/a(2))^2)+a(3)*(exp((0.43)*(x-a(1)))+(0.074))*(1.0-erf((x-a(1))/((0.87)*a(2)*sqrt(2.))))*(1.0+erf((x-a(1)+strip_thresh)/((0.87)*a(2)*sqrt(2.)))) ;1173 keV
  return
end

PRO tailing_v4,x,b,f ;short exponential + step 
  COMMON PEAK
  f=x
  f=b(0)*(exp((0.44)*(x-b(1)))+(0.075))*(1.0-erf((x-b(1))/((0.87)*b(2)*sqrt(2.))))*(1.0+erf((x-b(1)+strip_thresh)/((0.87)*b(2)*sqrt(2.))))
  ; f=b(0)*(exp((0.48)*(x-b(1)))+(0.098))*(1.0-erf((x-b(1))/((0.88)*b(2)*sqrt(2.))))*(1.0+erf((x-b(1)+strip_thresh)/((0.88)*b(2)*sqrt(2.)))) ;60 keV
  ; f=b(0)*(exp((0.43)*(x-b(1)))+(0.074))*(1.0-erf((x-b(1))/((0.87)*b(2)*sqrt(2.))))*(1.0+erf((x-b(1)+strip_thresh)/((0.87)*b(2)*sqrt(2.)))) ;1173 keV
  return
end


PRO tailpeak_unconstrained_v5,x,a,f ;double exponential
COMMON PEAK
  f=x
  f=a(0)*exp(-0.5*((x-a(1))/a(2))^2)+(a(3)*exp(a(4)*(x-a(1)))+a(5)*exp(a(6)*(x-a(1))))*(1.0-erf((x-a(1))/(a(7)*sqrt(2.))))*(1.0+erf((x-a(1)+strip_thresh)/(a(7)*sqrt(2.))))
  return
end

PRO tailing_unconstrained_v5,x,a,f ;double exponential
  COMMON PEAK
  f=x
  f=(a(0)*exp(a(2)*(x-a(1)))+a(3)*exp(a(4)*(x-a(1))))*(1.0-erf((x-a(1))/(a(5)*sqrt(2.))))*(1.0+erf((x-a(1)+strip_thresh)/(a(5)*sqrt(2.))))
  return
end

PRO tailpeak_v5,x,a,f ;double exponential
  COMMON PEAK
  f=x
  f=a(0)*exp(-0.5*((x-a(1))/a(2))^2)+(a(3)*exp((0.520)*(x-a(1)))+a(4)*exp((0.045)*(x-a(1))))*(1.0-erf((x-a(1))/((0.84)*a(2)*sqrt(2.))))*(1.0+erf((x-a(1)+strip_thresh)/((0.84)*a(2)*sqrt(2.))))
 ; f=a(0)*exp(-0.5*((x-a(1))/a(2))^2)+(a(3)*exp((0.562)*(x-a(1)))+a(4)*exp((0.036)*(x-a(1))))*(1.0-erf((x-a(1))/((0.84)*a(2)*sqrt(2.))))*(1.0+erf((x-a(1)+strip_thresh)/((0.84)*a(2)*sqrt(2.)))) ;60 keV
 ; f=a(0)*exp(-0.5*((x-a(1))/a(2))^2)+(a(3)*exp((0.495)*(x-a(1)))+a(4)*exp((0.044)*(x-a(1))))*(1.0-erf((x-a(1))/((0.84)*a(2)*sqrt(2.))))*(1.0+erf((x-a(1)+strip_thresh)/((0.84)*a(2)*sqrt(2.)))) ;1173 keV
  return
end

PRO tailing_v5,x,b,f ;double exponential
  COMMON PEAK
  f=x
  f=(b(0)*exp((0.520)*(x-b(1)))+b(3)*exp((0.045)*(x-b(1))))*(1.0-erf((x-b(1))/((0.84)*b(2)*sqrt(2.))))*(1.0+erf((x-b(1)+strip_thresh)/((0.84)*b(2)*sqrt(2.))))
 ; f=(b(0)*exp((0.562)*(x-b(1)))+b(3)*exp((0.036)*(x-b(1))))*(1.0-erf((x-b(1))/((0.84)*b(2)*sqrt(2.))))*(1.0+erf((x-b(1)+strip_thresh)/((0.84)*b(2)*sqrt(2.)))) ; 60 keV
 ; f=(b(0)*exp((0.495)*(x-b(1)))+b(3)*exp((0.044)*(x-b(1))))*(1.0-erf((x-b(1))/((0.84)*b(2)*sqrt(2.))))*(1.0+erf((x-b(1)+strip_thresh)/((0.84)*b(2)*sqrt(2.)))) ; 1173 keV
  return
end

PRO tailpeak_unconstrained_v6,x,a,f ;short exponential + sloped step
  COMMON PEAK
  f=x
  f=a(0)*exp(-0.5*((x-a(1))/a(2))^2)+(a(3)*exp(a(4)*(x-a(1)))+a(5)*(1.+a(6)*(x-a(1))))*(1.0-erf((x-a(1))/(a(7)*sqrt(2.))))*(1.0+erf((x-a(1)+strip_thresh)/(a(7)*sqrt(2.))))
  return
end

PRO tailing_unconstrained_v6,x,a,f ;short exponential + sloped step
  COMMON PEAK
  f=x
  f=(a(0)*exp(a(2)*(x-a(1)))+a(3)*(1.+a(4)*(x-a(1))))*(1.0-erf((x-a(1))/(a(5)*sqrt(2.))))*(1.0+erf((x-a(1)+strip_thresh)/(a(5)*sqrt(2.))))
  return
end

PRO tailpeak_v6,x,a,f ;short exponential + sloped step
  COMMON PEAK
  f=x
;  f=a(0)*exp(-0.5*((x-a(1))/a(2))^2)+a(3)*(exp((0.51)*(x-a(1)))+(0.13)*(1.+0.028*(x-a(1))))*(1.0-erf((x-a(1))/((0.84)*a(2)*sqrt(2.))))*(1.0+erf((x-a(1)+strip_thresh)/((0.84)*a(2)*sqrt(2.))))
;  f=a(0)*exp(-0.5*((x-a(1))/a(2))^2)+a(3)*(exp((0.55)*(x-a(1)))+(0.15)*(1.+0.025*(x-a(1))))*(1.0-erf((x-a(1))/((0.85)*a(2)*sqrt(2.))))*(1.0+erf((x-a(1)+strip_thresh)/((0.85)*a(2)*sqrt(2.)))) ;60 keV
  f=a(0)*exp(-0.5*((x-a(1))/a(2))^2)+a(3)*(exp((0.49)*(x-a(1)))+(0.13)*(1.+0.028*(x-a(1))))*(1.0-erf((x-a(1))/((0.84)*a(2)*sqrt(2.))))*(1.0+erf((x-a(1)+strip_thresh)/((0.84)*a(2)*sqrt(2.)))) ;1173 keV
  return
end

PRO tailing_v6,x,b,f ;short exponential + sloped step
  COMMON PEAK
  f=x
;  f=b(0)*(exp((0.51)*(x-b(1)))+(0.13)*(1.+0.028*(x-b(1))))*(1.0-erf((x-b(1))/((0.84)*b(2)*sqrt(2.))))*(1.0+erf((x-b(1)+strip_thresh)/((0.84)*b(2)*sqrt(2.))))
;  f=b(0)*(exp((0.55)*(x-b(1)))+(0.15)*(1.+0.025*(x-b(1))))*(1.0-erf((x-b(1))/((0.85)*b(2)*sqrt(2.))))*(1.0+erf((x-b(1)+strip_thresh)/((0.85)*b(2)*sqrt(2.)))) ;60 keV
  f=b(0)*(exp((0.49)*(x-b(1)))+(0.13)*(1.+0.028*(x-b(1))))*(1.0-erf((x-b(1))/((0.84)*b(2)*sqrt(2.))))*(1.0+erf((x-b(1)+strip_thresh)/((0.84)*b(2)*sqrt(2.)))) ;1173 keV
  return
end


pro peak_fit_extended, energy=energy, drift_time=drift_time, extended=extended, ps=ps, eps=eps

loadct,0 
postscr=keyword_set(ps)
epostscr=keyword_set(eps)
if (postscr or epostscr) then begin                                                              
  set_plot,'ps'
  device,/inches,xsize=8.5,ysize=3.0,xoffset=0.0,yoffset=1.00
  if epostscr then $
    device,filename='f7b.eps',/encapsulated,/preview,/color $
  else device,filename='plt_spectra.ps',encapsulated=0,preview=0,/color
  endif else begin
  window, xsize=1200, ysize=400
endelse

COMMON PEAK

photon_energy = 661.66;keV, photopeak energy
if (keyword_set(energy)) then photon_energy=energy

practical_range=(0.8)*(1.e1)*(0.55)*photon_energy*(1.-0.9841/(1.+0.0030*photon_energy))/(5.323)
R0=0.0
if (keyword_set(extended)) then R0=practical_range/2.
print, "range: ", R0

tau=250. ;drift time, nanoseconds
if (keyword_set(drift_time)) then tau=drift_time

strip_thresh = 18. ;keV, energy threshold for triggering neighbor strip
sigma_energy = (2.17+0.65*(photon_energy/1000.)^0.5)/2.35 ;keV, intrisic 1-sigma spectral resolution at photopeak energy
;print,'FWHM: ', 2.35*sigma_energy
print,'sigma: ', sigma_energy




half_pitch = 2000./2. ;micrometers, strip half pitch
emin=float(fix(photon_energy-22.))
emax=emin+30.
dechn=0.2 ; keV
nchn=fix((emax-emin)/dechn)
c1=8 ;start cannel for fitting
c2=136 ;end channel for fitting
spce=intarr(nchn) ;electron photopeak spectrum
spce(*)=0
losse=intarr(nchn) ;electron charge sharing spectrum
losse(*)=0
spch=intarr(nchn) ;hole photopeak spectrum
spch(*)=0
lossh=intarr(nchn) ;hole charge sharing spectrum
lossh(*)=0
echn=emin+dechn*findgen(nchn)+dechn/2.0

ndx=100000 ;number of spatial points sampled
dx=half_pitch/ndx ;micrometers
x=findgen(ndx)*dx


;ce=diffusion_distribution(x,tau,/electron)
;ce=repulsion_distribution(x,tau,photon_energy,/electron)
ce=charge_distribution_extended(x,tau,photon_energy,R0,/electron)

;ch=diffusion_distribution(x,tau,/hole)
;ch=repulsion_distribution(x,tau,photon_energy,/hole)
ch=charge_distribution_extended(x,tau,photon_energy,R0,/hole)


for i=0,ndx-1 do begin
  eloss=total(ce(i:ndx-1))*dx*photon_energy
  if (eloss lt strip_thresh) then begin
    nrge=photon_energy-eloss+sigma_energy*randomn(seed)
    if (eloss gt 0.1) then begin ; charge loss spectrum
      chne=0
      if (nrge ge emin) then begin
        chne=fix((nrge-emin)/dechn)
        losse(chne)+=1
      endif
    endif else begin ; photopeak spectrum
      chne=0
      if (nrge ge emin) then begin
        chne=fix((nrge-emin)/dechn)
        spce(chne)+=1
      endif
    endelse 
  endif
endfor

;print,'electrons peak, tail: ',total(spce),total(losse)

yheight=1.1*max(spce+losse) ;top of plot
ybase=5 ;bottom of plot

!p.multi=[0,2,1]
!x.range=[emin,emax]
!x.type=0
!x.style=1
!y.range=[ybase,yheight]
!y.type=0
!y.style=1
!x.title='Energy [keV]'
!y.title='Counts'

!p.title='Electron Signal'
plot,echn,spce+losse,psym=10,/ylog
oplot,echn,losse,psym=10
oplot,echn,spce,psym=10
oplot,[photon_energy,photon_energy],[ybase,yheight],psym=0,linestyle=2

tfit=fltarr(nchn)
weights=fltarr(nchn)
weights(*)=1./((spce(*)+losse(*))>1)
c=fltarr(4) ;fitting parameters
sigc=c
quality=0.0
c(0)=max(spce)
c(1)=photon_energy
c(2)=1.0
c(3)=max(losse)/4.
;c(4)=0.15*c(3)
;c(4)=0.55
;c(5)=0.15*c(3)
;c(6)=0.045
;c(7)=0.9*c(2)
ffit=CURVEFIT(echn(c1:c2),spce(c1:c2)+losse(c1:c2),weights(c1:c2),c,sigc,CHISQ=quality,FUNCTION_NAME='tailpeak_v6',ITMAX=100,/NODERIVATIVE,TOL=1.0E-5); Constrained model, 5 parameters
;ffit=CURVEFIT(echn(c1:c2),spce(c1:c2)+losse(c1:c2),weights(c1:c2),c,sigc,CHISQ=quality,FUNCTION_NAME='photopeak',ITMAX=100,/NODERIVATIVE,TOL=1.0E-5); Gaussian only, 3 parameters
;ffit=CURVEFIT(echn(c1:c2),spce(c1:c2)+losse(c1:c2),weights(c1:c2),c,sigc,CHISQ=quality,FUNCTION_NAME='tailpeak_unconstrained_v6',/NODERIVATIVE); Unconstrained model, 8 parameters
;print,'electrons Eo, FWHM, Chi^2: ',c(1),c(2),quality
oplot,echn(c1:c2),ffit

;print,'electrons: ',c(1),c(2),total(spce+losse),total(ffit),quality


a=fltarr(3)
a(0:2)=c(0:2)
photopeak,echn(c1:c2),a,pfit
oplot,echn(c1:c2),pfit
b=fltarr(3); b is used for tailing_v6
b(0)=c(3)
b(1:2)=c(1:2)
tailing_v6,echn(c1:c2),b,tfit
;d=fltarr(6) ;d is used for tailing_unconstrained_v5
;d(0)=c(3)
;d(1)=c(1)
;d(2:5)=c(4:7)
;tailing_unconstrained_v6,echn(c1:c2),d,tfit
oplot,echn(c1:c2),tfit

;print,'electrons peak, tail fits: ',total(pfit),total(tfit)

print,'electrons: ',c(1),c(2),total(spce),total(pfit),total(losse),total(tfit),quality
;print,'electrons: ',c(1),c(2),total(spce+losse),total(ffit),quality ;Gaussian fit




for i=0,ndx-1 do begin
  eloss=total(ch(i:ndx-1))*dx*photon_energy
  if (eloss lt strip_thresh) then begin
    nrge=photon_energy-eloss+sigma_energy*randomn(seed)
    if (eloss gt 0.1) then begin ; charge loss spectrum
     chne=0
      if (nrge ge emin) then begin
        chne=fix((nrge-emin)/dechn)
        lossh(chne)+=1
      endif
    endif else begin ; photopeak spectrum
        chne=0
      if (nrge ge emin) then begin
        chne=fix((nrge-emin)/dechn)
        spch(chne)+=1
      endif
    endelse
  endif
endfor


;print,'holes peak, tail: ',total(spch),total(lossh)


!p.title='Hole Signal'
plot,echn,spch+lossh,psym=10,/ylog
oplot,echn,lossh,psym=10
oplot,echn,spch,psym=10
oplot,[photon_energy,photon_energy],[ybase,yheight],psym=0,linestyle=2

hfit=fltarr(nchn)
weights(*)=1./((spch(*)+lossh(*))>1)
c=fltarr(4) ;fitting parameters
sigc=c
quality=0.0
c(0)=max(spch)
c(1)=photon_energy
c(2)=1.0
c(3)=max(lossh)/4.
;c(4)=0.15*c(3)
;c(4)=0.55
;c(5)=0.15*c(3)
;c(6)=0.045
;c(7)=0.9*c(2)
ffit=CURVEFIT(echn(c1:c2),spch(c1:c2)+lossh(c1:c2),weights(c1:c2),c,sigc,CHISQ=quality,FUNCTION_NAME='tailpeak_v6',ITMAX=100,/NODERIVATIVE,TOL=1.0E-5)
;ffit=CURVEFIT(echn(c1:c2),spch(c1:c2)+lossh(c1:c2),weights(c1:c2),c,sigc,CHISQ=quality,FUNCTION_NAME='photopeak',ITMAX=100,/NODERIVATIVE,TOL=1.0E-5)
;ffit=CURVEFIT(echn(c1:c2),spch(c1:c2)+lossh(c1:c2),weights(c1:c2),c,sigc,CHISQ=quality,FUNCTION_NAME='tailpeak_unconstrained_v6',/NODERIVATIVE); Unconstrained model, 7 parameters
;print,'holes Eo, FWHM, Chi^2: ',c(1),c(2),quality
oplot,echn(c1:c2),ffit

;print,'holes: ',c(1),c(2),total(spch+lossh),total(ffit),quality

a=fltarr(3)
a(0:2)=c(0:2)
photopeak,echn(c1:c2),a,pfit
oplot,echn(c1:c2),pfit
b=fltarr(3); b is used for tailing_v6
b(0)=c(3)
b(1:2)=c(1:2)
tailing_v6,echn(c1:c2),b,tfit
;d=fltarr(6);d is used for tailing_unconstrained
;d(0)=c(3)
;d(1)=c(1)
;d(2:5)=c(4:7)
;tailing_unconstrained_v6,echn(c1:c2),d,tfit
oplot,echn(c1:c2),tfit

;print,'holes peak, tail fits: ',total(pfit),total(tfit)
print,'holes: ',c(1),c(2),total(spch),total(pfit),total(lossh),total(tfit),quality
;print,'holes: ',c(1),c(2),total(spch+lossh),total(ffit),quality ;Gaussian fit



if (postscr or epostscr) then begin
  device,/close
  set_plot,'X'
endif 

return
end


