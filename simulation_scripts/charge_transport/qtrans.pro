;*****************************************************************************
; FILE_NAME: QTRANS
;
; PURPOSE: Simulates the charge transport and signal induction on several 
;          strips of a Ge detector
;          
; INCLUDED PROCEDURES/FUNCTIONS:
;         VECTGRID()
;         INITIAL_POS
;         QMOVE()
;         DRIFT_DIR()
;         QTRANS
;
; FULL DOCUMENTATION: Search for ';+****' below          
; REVISION HISTORY: Based on FIELD3D by Steve Boggs
;                   Update by Mark Bandstra UCB/SSL 2005
;                   removed vector weighting field
;       Revised by SEB UCSD 2023
;*****************************************************************************


;+****************************************************************
; NAME: VECTGRID
;
; PURPOSE: Derive the electric and weighting fields at the grid point from the scalar
; potentials at those points.
;
; CALLING SEQUENCE: egrid=vectgrid(eugrid)
; Note in the update QTRANS we only need the vector efield, not wfield. SEB 2023
;
;
; NOTES: -This program will NOT compile unless you first call the
;        procedure DEFINE_COMMON_BLKS.PRO, which defines the
;        common blocks used in FIELD3D. Do this at least once
;        per IDL session: IDL> define_common_blks.
;        Define_common_blks is also called in field3d, so you
;        don't need to manually run it again even if changes are made.
;
;
; EXAMPLE:
; - make a scalar field using field3d
;    IDL> restore,outfile_dir+'efield.idlsave'   ;restore the saved scalar field
;    IDL> egrid=vectgrid(eugrid) ;return the 3D vector field
;    IDL> save,egrid,FILENAME=outfile_dir+'vefield.sav' ;save if desired, not required
;Note this routine used to require removing the Air regions from grid before calling,
;now this is done internal to the routine below. SEB 2023
;
; PROCEDURES CALLED:
; DEFINE_COMMON_BLKS
;
; REVISION HISTORY: Mod by Susan Amrose    UCB/SSL   2002
;       Revised by SEB UCSD 2023
;-****************************************************************

FUNCTION vectgrid,ugrid


;------------------
;GET COMMON BLOCKS
;------------------
resolve_routine,'field3d',/no_recompile,/compile_full_file ;for access to efgeom geometry function
COMMON FLDGRID ; IMAX,JMAX,KMAX,KAIR,DX,DY,DZA,DZB,HEIGHT,NSTRIP
COMMON PIXEL   ; ICNTCT,ISTRIP
COMMON TYPES   ; GND,HV,Air,Semiconductor,Fixed,Variable

;------------------------------------
;SET UP CONSTANTS, REMOVE AIR REGIONS
;------------------------------------
zmax=KMAX-2*KAIR ;Z dimension when Air regions removed
ugrid=temporary(ugrid[*,*,KAIR:KMAX-KAIR-1]) ;Remove Air regions from the scalar potential field
wucgrid=ugrid  ;copy input array
;This field (wucgrid) provides only the scalar potential field.
;SEB 2023: not clear why wucgrid is needed, kept now for legacy

;---------------------------
;ESTABLISH DETECTOR GEOMETRY
;---------------------------

; Complication to this field calculation business
; -At the boundaries of the grid, we only know the potential field on one side
;  (within the grid), but not the other.  Thus, a different formula is used.
; -On the pixel plane, metal electrodes (contact and steering) create
;  additional boundaries, which also call for a different formula.  Thus,
;  the geometry of the pixel plane becomes relevant in this calculation.


;--------------------------------------
;LOCATE METAL-SEMICONDUCTOR BOUNDARIES
;--------------------------------------
; index small       large
; -x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-
;   metal  |  semiconductor  |  metal
;       l_edge    r_edge
;
; l_edge == left edge of a semiconductor region, r_edge == right edge

; The following grids are (IMAX-5) by JMAX and IMAX by (JMAX-5),
; excluding two and three rows/columns on each side of the pixelgeom[] grid.
;
; 0 1 2 3   IMAX-1
; x-x-x-x- . . . -x-x-x-x-x
;    |<- xl_edge[] ->|
;      |<- xr_edge[] ->|
;
; Note that i=2 and i=IMAX-3 are guaranteed not to be a right edge and a left
; edge respectively, by the definition of the geometry in field3d.pro;
; same in the y-direction.

;SEB 2023: I have simplified how fields are calcuated at the metal-semiconductor
;boundary and within the gaps in the pixel plane. The strip edges are now defined only
;as the single point where the contact stops and semiconductor begins. Now the fields
;in the gap are calculated using only the nearest-neighbor points. Previously the 
;third-order fits for gap points could lead to poor estimates and incorrect signs
;when the points straddled the symmetry axis of the gap. This simplified approach 
;produces much more stable electric field values in the gaps.
;--------------------------------------
;LOWER CONTACT PLANE
;--------------------------------------
;lower contact plane:
; First, reconstruct the pixel-plane geometry.
pixelgeom=(efgeom(imax,jmax,kmax,kair))[*,*,kair]
; efgeom() is defined in field3d.pro.
; Define all metal areas as GND here.
; Note that there is no type `air' on the pixel plane.
pixelgeom = GND * (temporary(pixelgeom) ne Semiconductor)

;IF ((ISTRIP-ICNTCT) eq 3) THEN BEGIN
; yl_edge_l=    (pixelgeom[*,3:JMAX-3] eq Semiconductor)$
;   and (pixelgeom[*,1:JMAX-5] eq GND) ;Originally hard-coded for 5-6 gap points!!!
; yr_edge_l=    (pixelgeom[*,2:JMAX-4] eq Semiconductor)$
;   and (pixelgeom[*,4:JMAX-2] eq GND) ;Originally hard-coded for 5-6 gap points!!!
;     print,'Here we are!'
;ENDIF ELSE IF ((ISTRIP-ICNTCT) eq 5) THEN BEGIN
; yl_edge_l=    (pixelgeom[*,3:JMAX-3] eq Semiconductor)$
;   and (pixelgeom[*,0:JMAX-6] eq GND) ;Now hardcoded for 7 gap points!
; yr_edge_l=    (pixelgeom[*,2:JMAX-4] eq Semiconductor)$
;   and (pixelgeom[*,5:JMAX-1] eq GND) ;Now hardcoded for 7 gap points!
;ENDIF ELSE BEGIN
; print,'Contact edges not defined properly!'  
;ENDELSE

yl_edge_l=    (pixelgeom[*,4:JMAX-5] eq GND)$
  and (pixelgeom[*,5:JMAX-4] eq Semiconductor) ;Only 1 point where contact ends
yr_edge_l=    (pixelgeom[*,4:JMAX-5] eq GND)$
  and (pixelgeom[*,3:JMAX-6] eq Semiconductor) ;Only 1 point where contact ends
y_gap_l = (pixelgeom[*,1:JMAX-2] eq Semiconductor) ;Gap points on the pixel plane

  

; Grids to mask out metal interiors on the pixel plane
;  (the name should actually be x_not_metal_interior)
;SEB 2023: Note that 0=lower, 1=upper; and this defines
;the last points on the contacts as nonmetal so that the vector fields
;are calculated for those points.
x_nonmetal=intarr(IMAX,JMAX,2,/nozero) ;last index, 0=lower, 1=upper
x_nonmetal[1:IMAX-2,*,0]=    (pixelgeom[0:IMAX-3,*] eq Semiconductor)$
  or (pixelgeom[2:IMAX-1,*] eq Semiconductor)
; The following 2 statements are simplified by the restrictions on
; the grid definition in field3d.pro.
x_nonmetal[0,*,0]=     (pixelgeom[1,*]    eq Semiconductor)
x_nonmetal[IMAX-1,*,0]=    (pixelgeom[IMAX-2,*]   eq Semiconductor)

y_nonmetal=intarr(IMAX,JMAX,2,/nozero)
y_nonmetal[*,1:JMAX-2,0]=    (pixelgeom[*,0:JMAX-3] eq Semiconductor)$
  or (pixelgeom[*,2:JMAX-1] eq Semiconductor)
y_nonmetal[*,0,0]=     (pixelgeom[*,1]    eq Semiconductor)
y_nonmetal[*,JMAX-1,0]=    (pixelgeom[*,JMAX-2]   eq Semiconductor)


;--------------------------------------
;UPPER CONTACT PLANE
;--------------------------------------
;upper contact plane:
pixelgeom=(efgeom(imax,jmax,kmax,kair))[*,*,kmax-kair-1]
; efgeom() is defined in field3d.pro.
; Define all metal areas as GND here.
; Note that there is no type `air' on the pixel plane.
pixelgeom = GND * (temporary(pixelgeom) ne Semiconductor)

;IF ((ISTRIP-ICNTCT) eq 3) THEN BEGIN
; xl_edge_u=    (pixelgeom[3:IMAX-3,*] eq Semiconductor)$
;  and (pixelgeom[1:IMAX-5,*] eq GND) ;Originally hard-coded for 5-6 gap points!!!
; xr_edge_u=    (pixelgeom[2:IMAX-4,*] eq Semiconductor)$
;  and (pixelgeom[4:IMAX-2,*] eq GND) ;Originally hard-coded for 5-6 gap points!!!
;  print,'Here we are!'
;ENDIF ELSE IF ((ISTRIP-ICNTCT) eq 5) THEN BEGIN
; xl_edge_u=    (pixelgeom[3:IMAX-3,*] eq Semiconductor)$
;  and (pixelgeom[0:IMAX-6,*] eq GND) ;Now hardcoded for 7 gap points!
; xr_edge_u=    (pixelgeom[2:IMAX-4,*] eq Semiconductor)$
;  and (pixelgeom[5:IMAX-1,*] eq GND) ;Now hardcoded for 7 gap points!
;ENDIF ELSE BEGIN
;  print,'Contact edges not defined properly!'
;ENDELSE

xl_edge_u=    (pixelgeom[4:IMAX-5,*] eq GND)$
  and (pixelgeom[5:IMAX-4,*] eq Semiconductor) ;Only 1 point where contact ends
xr_edge_u=    (pixelgeom[4:IMAX-5,*] eq GND)$
  and (pixelgeom[3:IMAX-6,*] eq Semiconductor) ;Only 1 point where contact ends
x_gap_u = (pixelgeom[1:IMAX-2,*] eq Semiconductor) ;Gap points on the pixel plane


; Grids to mask out metal interiors on the pixel plane
;  (the name should actually be x_not_metal_interior)
;these arrays are defined above
x_nonmetal[1:IMAX-2,*,1]=  (pixelgeom[0:IMAX-3,*] eq Semiconductor)$
  or (pixelgeom[2:IMAX-1,*] eq Semiconductor)
; The following 2 statements are simplified by the restrictions on
; the grid definition in field3d.pro.
x_nonmetal[0,*,1]=     (pixelgeom[1,*]    eq Semiconductor)
x_nonmetal[IMAX-1,*,1]=    (pixelgeom[IMAX-2,*]   eq Semiconductor)

y_nonmetal[*,1:JMAX-2,1]=  (pixelgeom[*,0:JMAX-3] eq Semiconductor)$
  or (pixelgeom[*,2:JMAX-1] eq Semiconductor)
y_nonmetal[*,0,1]=     (pixelgeom[*,1]    eq Semiconductor)
y_nonmetal[*,JMAX-1,1]=    (pixelgeom[*,JMAX-2]   eq Semiconductor)

pixelgeom=0b  ; save memory


;--------------------------------------
;CALCULATE ELECTRIC FIELDS
;--------------------------------------
;Now, we are ready to calculate the electric fields.
dwucdx=fltarr(IMAX,JMAX,zmax)
dwucdy=fltarr(IMAX,JMAX,zmax)
dwucdz=fltarr(IMAX,JMAX,zmax)


;--------------------------------------
;X-DERIVATIVES
;--------------------------------------
;Find the x-derivatives.
;Simple 4th-order derivative at ordinary points
dwucdx[2:IMAX-3,*,*]$
  =(-1./DX)*( (2./ 3.)*(wucgrid[3:IMAX-2,*,*]-wucgrid[1:IMAX-4,*,*])$
  -(1./12.)*(wucgrid[4:IMAX-1,*,*]-wucgrid[0:IMAX-5,*,*]) )

;3rd-order fit to field at left boundaries
dwucdx[0:1,*,*]=(-1./DX)*(  3.    *(wucgrid[1:2,*,*]-wucgrid[0:1,*,*])$
  - 1.5   *(wucgrid[2:3,*,*]-wucgrid[0:1,*,*])$
  +(1./3.)*(wucgrid[3:4,*,*]-wucgrid[0:1,*,*]) )

;3rd-order fit to field at right boundaries
dwucdx[IMAX-2:IMAX-1,*,*]$
  =(-1./DX)*(  3.    *(wucgrid[IMAX-2:IMAX-1,*,*]-wucgrid[IMAX-3:IMAX-2,*,*])$
  - 1.5   *(wucgrid[IMAX-2:IMAX-1,*,*]-wucgrid[IMAX-4:IMAX-3,*,*])$
  +(1./3.)*(wucgrid[IMAX-2:IMAX-1,*,*]-wucgrid[IMAX-5:IMAX-4,*,*]) )

; Metal-semiconductor boundaries on pixel plane
; upper plane Semiconductor left edges

;dwucdx[2:IMAX-4,*,ZMAX-1]$  ; First, reset the fields to null.
;  = (xl_edge_u eq 0) * dwucdx[2:IMAX-4,*,ZMAX-1]$
;  + xl_edge_u       * (-1./DX)$ ; Next, calculate the boundaries properly.
;  *(  3.    *(wucgrid[3:IMAX-3,*,ZMAX-1]-$
;  wucgrid[2:IMAX-4,*,ZMAX-1])$
;  - 1.5   *(wucgrid[4:IMAX-2,*,ZMAX-1]-$
;  wucgrid[2:IMAX-4,*,ZMAX-1])$
;  +(1./3.)*(wucgrid[5:IMAX-1,*,ZMAX-1]-$
;  wucgrid[2:IMAX-4,*,ZMAX-1]) )
;;recall that wucgrid has been resized to be zmax elements at start of program
;xl_edge_u=0b
;; upper plane Semiconductor right edges
;dwucdx[3:IMAX-3,*,ZMAX-1]$  ; First, reset the fields to null.
;  = (xr_edge_u eq 0) * dwucdx[3:IMAX-3,*,ZMAX-1]$
;  + xr_edge_u       * (-1./DX)$ ; Next, calculate the boundaries properly.
;  *(  3.    *(wucgrid[3:IMAX-3,*,ZMAX-1]-$
;  wucgrid[2:IMAX-4,*,ZMAX-1])$
;  - 1.5   *(wucgrid[3:IMAX-3,*,ZMAX-1]-$
;  wucgrid[1:IMAX-5,*,ZMAX-1])$
;  +(1./3.)*(wucgrid[3:IMAX-3,*,ZMAX-1]-$
;  wucgrid[0:IMAX-6,*,ZMAX-1]) )
;xr_edge_u=0b

dwucdx[4:IMAX-5,*,ZMAX-1]$  ; First, reset the fields to null.
  = (xl_edge_u eq 0) * dwucdx[4:IMAX-5,*,ZMAX-1]$
  + xl_edge_u       * (-1./DX)$ ; Next, calculate the boundaries properly.
  *(  3.    *(wucgrid[5:IMAX-4,*,ZMAX-1]-$
  wucgrid[4:IMAX-5,*,ZMAX-1])$
  - 1.5   *(wucgrid[6:IMAX-3,*,ZMAX-1]-$
  wucgrid[4:IMAX-5,*,ZMAX-1])$
  +(1./3.)*(wucgrid[7:IMAX-2,*,ZMAX-1]-$
  wucgrid[4:IMAX-5,*,ZMAX-1]) )
xl_edge_u=0b
; upper plane Semiconductor right edges
dwucdx[4:IMAX-5,*,ZMAX-1]$  ; First, reset the fields to null.
  = (xr_edge_u eq 0) * dwucdx[4:IMAX-5,*,ZMAX-1]$
  + xr_edge_u       * (-1./DX)$ ; Next, calculate the boundaries properly.
  *(  3.    *(wucgrid[4:IMAX-5,*,ZMAX-1]-$
  wucgrid[3:IMAX-6,*,ZMAX-1])$
  - 1.5   *(wucgrid[4:IMAX-5,*,ZMAX-1]-$
  wucgrid[2:IMAX-7,*,ZMAX-1])$
  +(1./3.)*(wucgrid[4:IMAX-5,*,ZMAX-1]-$
  wucgrid[1:IMAX-8,*,ZMAX-1]) )
xr_edge_u=0b
 ; lower plane Semiconductor gaps
dwucdx[1:IMAX-2,*,ZMAX-1]$  ; First, reset the fields to null.
   = (x_gap_u eq 0) * dwucdx[1:IMAX-2,*,ZMAX-1]$
   + x_gap_u       * (-1./DX)$ ; Next, calculate the boundaries properly.
   *(  0.5    *(wucgrid[2:IMAX-1,*,ZMAX-1]-$
   wucgrid[0:IMAX-3,*,ZMAX-1]) )
 x_gap_u=0b

; Finally, mask out all the metal interior on both the upper and lower plane.
dwucdx[*,*,ZMAX-1]= dwucdx[*,*,ZMAX-1] * x_nonmetal[*,*,1]
dwucdx[*,*,0]= dwucdx[*,*,0] * x_nonmetal[*,*,0] ;lower plane
x_nonmetal=0b


;--------------------------------------
;Y-DERIVATIVES
;--------------------------------------
;Find the y-derivatives.
;Simple 4th-order derivative at ordinary points
dwucdy[*,2:JMAX-3,*]$
  =(-1./DY)*( (2./ 3.)*(wucgrid[*,3:JMAX-2,*]-wucgrid[*,1:JMAX-4,*])$
  -(1./12.)*(wucgrid[*,4:JMAX-1,*]-wucgrid[*,0:JMAX-5,*]) )

;3rd-order fit to field at left boundaries
dwucdy[*,0:1,*]=(-1./DY)*(  3.    *(wucgrid[*,1:2,*]-wucgrid[*,0:1,*])$
  - 1.5   *(wucgrid[*,2:3,*]-wucgrid[*,0:1,*])$
  +(1./3.)*(wucgrid[*,3:4,*]-wucgrid[*,0:1,*]) )

;3rd-order fit to field at right boundaries
dwucdy[*,JMAX-2:JMAX-1,*]$
  =(-1./DY)*(  3.    *(wucgrid[*,JMAX-2:JMAX-1,*]-wucgrid[*,JMAX-3:JMAX-2,*])$
  - 1.5   *(wucgrid[*,JMAX-2:JMAX-1,*]-wucgrid[*,JMAX-4:JMAX-3,*])$
  +(1./3.)*(wucgrid[*,JMAX-2:JMAX-1,*]-wucgrid[*,JMAX-5:JMAX-4,*]) )

; Metal-semiconductor boundaries on pixel plane
; lower plane Semiconductor left edges

;dwucdy[*,2:JMAX-4,0]$ ; First, reset the fields to null.
;  = (yl_edge_l eq 0) * dwucdy[*,2:JMAX-4,0]$
;  + yl_edge_l       * (-1./DY)$ ; Next, calculate the boundaries properly.
;  *(  3.    *(wucgrid[*,3:JMAX-3,0]-$
;  wucgrid[*,2:JMAX-4,0])$
;  - 1.5   *(wucgrid[*,4:JMAX-2,0]-$
;  wucgrid[*,2:JMAX-4,0])$
;  +(1./3.)*(wucgrid[*,5:JMAX-1,0]-$
;  wucgrid[*,2:JMAX-4,0]) )
;;recall that wucgrid has been resized to be zmax elements
;yl_edge_l=0b
;; lower plane Semiconductor right edges
;dwucdy[*,3:JMAX-3,0]$ ; First, reset the fields to null.
;  = (yr_edge_l eq 0) * dwucdy[*,3:JMAX-3,0]$
;  + yr_edge_l       * (-1./DY)$ ; Next, calculate the boundaries properly.
;  *(  3.    *(wucgrid[*,3:JMAX-3,0]-$
;  wucgrid[*,2:JMAX-4,0])$
;  - 1.5   *(wucgrid[*,3:JMAX-3,0]-$
;  wucgrid[*,1:JMAX-5,0])$
;  +(1./3.)*(wucgrid[*,3:JMAX-3,0]-$
;  wucgrid[*,0:JMAX-6,0]) )
; yr_edge_l=0b
 
 dwucdy[*,4:JMAX-5,0]$ ; First, reset the fields to null.
   = (yl_edge_l eq 0) * dwucdy[*,4:JMAX-5,0]$
   + yl_edge_l       * (-1./DY)$ ; Next, calculate the boundaries properly.
   *(  3.    *(wucgrid[*,5:JMAX-4,0]-$
   wucgrid[*,4:JMAX-5,0])$
   - 1.5   *(wucgrid[*,6:JMAX-3,0]-$
   wucgrid[*,4:JMAX-5,0])$
   +(1./3.)*(wucgrid[*,7:JMAX-2,0]-$
   wucgrid[*,4:JMAX-5,0]) )
 ;recall that wucgrid has been resized to be zmax elements
 yl_edge_l=0b
 ; lower plane Semiconductor right edges
 dwucdy[*,4:JMAX-5,0]$ ; First, reset the fields to null.
   = (yr_edge_l eq 0) * dwucdy[*,4:JMAX-5,0]$
   + yr_edge_l       * (-1./DY)$ ; Next, calculate the boundaries properly.
   *(  3.    *(wucgrid[*,4:JMAX-5,0]-$
   wucgrid[*,3:JMAX-6,0])$
   - 1.5   *(wucgrid[*,4:JMAX-5,0]-$
   wucgrid[*,2:JMAX-7,0])$
   +(1./3.)*(wucgrid[*,4:JMAX-5,0]-$
   wucgrid[*,1:JMAX-8,0]) )
 yr_edge_l=0b
 ; lower plane Semiconductor gaps
 dwucdy[*,1:JMAX-2,0]$ ; First, reset the fields to null.
   = (y_gap_l eq 0) * dwucdy[*,1:JMAX-2,0]$
   + y_gap_l       * (-1./DY)$ ; Next, calculate the boundaries properly.
   *(  0.5    *(wucgrid[*,2:JMAX-1,0]-$
   wucgrid[*,0:JMAX-2,0]) )
 y_gap_l=0b


; Finally, mask out all the metal interior on both the upper and lower plane.
dwucdy[*,*,ZMAX-1]= dwucdy[*,*,ZMAX-1] * y_nonmetal[*,*,1]
dwucdy[*,*,0]= dwucdy[*,*,0] * y_nonmetal[*,*,0] ;lower plane
y_nonmetal=0b


;--------------------------------------
;Z-DERIVATIVES
;--------------------------------------
;Find the z-derivatives.
; ordinary points
dwucdz[*,*,2:zmax-3]$
  =(-1./DZA)*( (2./ 3.)*(wucgrid[*,*,3:zmax-2]-$
  wucgrid[*,*,1:ZMAX-4])$
  -(1./12.)*(wucgrid[*,*,4:ZMAX-1]-$
  wucgrid[*,*,0:ZMAX-5]) )

; bottom boundaries
dwucdz[*,*,0:1]=(-1./DZA)*(  3.    *(wucgrid[*,*,1:2]-wucgrid[*,*,0:1])$
  - 1.5   *(wucgrid[*,*,2:3]-wucgrid[*,*,0:1])$
  +(1./3.)*(wucgrid[*,*,3:4]-wucgrid[*,*,0:1]) )
;check for sign consistency between three-point fit at edges and simple slope using nearest neighbor
Check=(-1./DZA)*(wucgrid[*,*,1:2]-wucgrid[*,*,0:1])
signchange= (abs(check)/check) ne (abs(dwucdz[*,*,0:1])/dwucdz[*,*,0:1])
;if there is a sign change, use the simple slope.
dwucdz[*,*,0:1] = signchange*Check + (1b-signchange)*dwucdz[*,*,0:1]
print,total(signchange),' sign switches for bottom bndry'
; top boundaries
dwucdz[*,*,zmax-2:zmax-1]$
  =(-1./DZA)*(  3.    *(wucgrid[*,*,zmax-2:zmax-1]-$
  wucgrid[*,*,zmax-3:zmax-2])$
  - 1.5   *(wucgrid[*,*,zmax-2:zmax-1]-$
  wucgrid[*,*,zmax-4:zmax-3])$
  +(1./3.)*(wucgrid[*,*,zmax-2:zmax-1]-$
  wucgrid[*,*,zmax-5:zmax-4]))

;check for sign consistency between three-point fit at edges and simple slope using nearest neighbor
Check=(-1./DZA)*(wucgrid[*,*,zmax-2:zmax-1]-$
  wucgrid[*,*,zmax-3:zmax-2])
signchange= (abs(check)/check) ne $
  (abs(dwucdz[*,*,zmax-2:zmax-1])/dwucdz[*,*,zmax-2:zmax-1])
;if there is a sign change, use the simple slope.
dwucdz[*,*,zmax-2:zmax-1] = signchange*Check + (1b-signchange)*dwucdz[*,*,zmax-2:zmax-1]
print,total(signchange),' sign switches for top bndry'

wucgrid=0b

vgrid=fltarr(IMAX,JMAX,zmax,3,/nozero)
vgrid[*,*,*,0] = dwucdx
vgrid[*,*,*,1] = dwucdy
vgrid[*,*,*,2] = dwucdz

dwucdx=0b
dwucdy=0b
dwucdz=0b

RETURN,vgrid
END ;vectgrid()


;*****************************************************************************
; PROCEDURE NAME: INITIAL_POS
;
; PURPOSE: To sort out the many ways in which simulated event 
;           parameters can be chosen
;
; SYNTAX: initial_pos, zin, xin, yin, ein, qin, num, atten, random=random, 
;         event_atten=event_atten, expdist=expdist
; 
; WAYS TO DEFINE INITIAL PARAMS:
;     Every input event in the simulation must be defined by a 3D initial 
;   interaction position (Xin, Yin, Zin) and an energy (ein or atten). There
;   are several different ways to define these using input params.
;   INPUTS are described in the QTRANS main documentation below.
;    
; OUTPUTS: qin - charge associated with ein
;          num - number of total events
;          atten - chosen attenuation
;*****************************************************************************
PRO initial_pos, zin, xin, yin, ein, qin, num, atten, random=random, event_atten=event_atten, expdist=expdist

COMMON FLDGRID
COMMON PIXEL
COMMON DETECTOR
QE=1.60217733d-19 ; elementary charge

;Determine number of objects
expdist=keyword_set(expdist)? expdist:0
num = ( expdist GT 1 )? expdist: $
   (n_elements(zin) > (n_elements(xin)*(1-keyword_set(random)))) $
  > (n_elements(yin)*(1-keyword_set(random)))

IF (num EQ 0) THEN $
  IF (n_elements(ein) GT 1) THEN num = n_elements(ein) $
  ELSE num=2000
                  ;default 2000 objects, 122 keV exp dist

;---------------------------------
;DETERMINE X,Y,Z INITIAL POSITIONS
;---------------------------------

;choose random grid of x,y from midgap_left to midgap_rt of center strip 
xstr='' & ystr='' & xset=0 & yset=0
IF keyword_set(random) OR (n_elements(xin) EQ 0) THEN BEGIN 
  xin=(randomu(seed,num)*istrip +(imax/2-istrip/2))*DX 
  xstr=' X' & xset=1
ENDIF ELSE IF n_elements(xin) NE num THEN message, $
  'X0 must have same number of elements as Z0'

IF keyword_set(random) OR (n_elements(yin) EQ 0) then BEGIN
  yin=(randomu(seed,num)*istrip +(jmax/2-istrip/2))*DY
  ystr='Y ' & yset=1
ENDIF ELSE IF n_elements(yin) NE num THEN message, $
  'Y0 must have same number of elements as Z0'

both=(xset+yset EQ 2)? ',':' '
IF (xset+yset NE 0) THEN print,'Using random'+xstr+both+ystr+'start positions'

;--------------------------
;DETERMINE INITIAL CHARGE
;--------------------------

ne0=n_elements(ein)
CASE 1 OF 
  
  (ne0 EQ 1): BEGIN
    energy_type = strtrim(ein[0],2)+' keV'
    qin=replicate(QE*ein[0]*1e3/2.98d, num)
    atten=152
    ;atten=get_atten(ein[0])
  END

  (ne0 EQ 0): BEGIN
    energy_type = 'by default 122 keV'
    qin=replicate(QE*122.0e3/2.98d, num)
    ein=[122.0]
    atten=152d
    ;atten=get_atten(122)
  END

  ELSE: BEGIN
    IF (ne0 NE num) THEN message,'E0 must have same number of elements as Z0'
    energy_type = 'Input Energies'
    qin=QE*ein*1e3/2.98d
    atten=0
  ENDELSE

ENDCASE

atten=keyword_set(event_atten)? event_atten: atten
print,'Energies are '+energy_type

;---------------------------
;FINALLY GET Z DISTRIBUTION
;---------------------------
IF keyword_set(expdist) OR (n_elements(zin) EQ 0) THEN BEGIN
  MN=0.0 ;cm
  MX=HEIGHT ;cm
;  zin=randomu_exp(num,-1.*atten,MN,MX,/error)  ; -1.* added MEB, 10/18/05 to fix AC/DC illumination
; SEB 3/20/23, commented out randomu_exp until can get around to verifying it
  IF zin[0] EQ 0 THEN BEGIN
    zin=randomu(seed,num)*(MX-MN) + MN
    print,'Attenuation undefined, random Z Dist being used' 
  ENDIF ELSE print,'Exponential Z0, attenuation: ',strtrim(atten,2)+' cm^-1'
endif
;* note that expdist is used if z0 is NOT defined even if expdist is not set.

return
END  ; initial_pos

; Note the electric and weighting fields assume that the bias voltage is
; applied to the electrode at z=HEIGHT. Signals are taken from the AC electrode
; at z=HEIGHT and the DC (ground) electrode at z=0.

;*****************************************************************************
; FUNCTION NAME: DRIFT_DIR
;
; PURPOSE: Find the drift velocity and direction of holes and
;          electrons given the vector electric field.
;
; INPUTS: e - vector electric field
; e[nums,2,3] has 3 dimensions: internaction #, e/h, and the 3 spatial directions
;
; OUTPUT: Drift Velocity (mag and dir) for each hole and electron
;
; NOTES: Default crystal orientation is <100>.
; Drift velocities come from Ottaviani,Canali & Quaranta, 
; IEEE Trans. Sci., Vol NS 22, pp. 192-204, Feb 1975.
;*****************************************************************************  
FUNCTION drift_dir,e
; Calculate the direction (and magnitude) of the charge displacement

COMMON FLDGRID
COMMON STEP	; DTSTEP,NUMITS, dtstep_tolerance,efield_tolerance
; Only DTSTEP and efield_tolerance are used here.
COMMON DRIFT
COMMON FACTOR_CAL
tmp=size(e,/dimensions) & n0=tmp[0] & n1=tmp[1] & n2=tmp[2]

eabs=dblarr(n0,n1,n2,/nozero) ; magnitude of e[]
eabs[*,*,0]=sqrt(total(e^2,3)) ;[nx*ny,2,3] -- SEB 2023, check this command is correct!
eabs[*,*,1]=eabs[*,*,0]
eabs[*,*,2]=eabs[*,*,0]
vd=dblarr(n0,n1,n2)

;get electron velocities
; spline added MEB 2006
vd[*,0,*]= interpol(vd_e,e_e,eabs[*,0,*], /spline) * $
  (-1)*e[*,0,*]/(eabs[*,0,*] > efield_tolerance) * $
  (eabs[*,0,*] gt efield_tolerance)

;SEB converted to the version below 10/12/23 to match Julia SSD package
;vd[*,0,*]= eabs[*,0,*]*(38609./((1.+(eabs[*,0,*]/511.)^(0.805))^(1./0.805)) + 171.) * $
;  (-1)*e[*,0,*]/(eabs[*,0,*] > efield_tolerance) * $
;  (eabs[*,0,*] gt efield_tolerance)

electron_scale_factor=fac_e; 1.0, historically set to 130./85. to match the experimental range in tdif
vd[*,0,*]=electron_scale_factor * temporary(vd[*,0,*])

;get hole velocities
vd[*,1,*]=  interpol(vd_h,e_h,eabs[*,1,*], /spline) * $ 
  e[*,1,*]/(eabs[*,1,*] > efield_tolerance) * $
  (eabs[*,1,*] gt efield_tolerance)

;SEB converted to the version below 10/12/23 to match Julia SSD package
;vd[*,1,*]= eabs[*,1,*]*(61824./((1.+(eabs[*,1,*]/185.)^(0.942))^(1./0.942)) - 0.)  * $
;  e[*,1,*]/(eabs[*,1,*] > efield_tolerance) * $
;  (eabs[*,1,*] gt efield_tolerance)

hole_scale_factor=fac_h; 1.0, historically set to 135./110. to match the experimental range in tdif
vd[*,1,*]=hole_scale_factor * temporary(vd[*,1,*])

;print,'Drift speeds: ',interpol(vd_e,e_e,eabs[3,0,0], /spline),interpol(vd_h,e_h,eabs[3,1,0], /spline)
return, vd*DTSTEP 
END ;drift_dir()


;*****************************************************************************
; FUNCTION NAME: QMOVE
;
; PURPOSE: This is the main function for charge transport
; 
; Inputs and keywords described in QTRANS documentation below
;*****************************************************************************
FUNCTION qmove,x0,y0,z0,q0,egrid,wgridac,wgriddc,$
               trap=trap,XTRAP,nopulse=nopulse,save_path=save_path,$
               multich=multich,numsites=numsites,$
               sim=sim

COMMON PIXEL
COMMON FLDGRID
COMMON DETECTOR
COMMON STEP
COMMON DIR
COMMON TRAPPING

tmp=size(x0,/dimensions) & nx=tmp[0]>1 ;find number of interactions in x0
en2q = 1.60217733d-19 * 1d3 / 2.98d   ;(C/keV) for a 1MeV event
;-------------
;ARRAY SET UP
;-------------
qac=dblarr(nx,numits+3) ;charge induced at the AC electrode in double precision
qdc=qac	;ditto for the DC electrode
x=dblarr(nx,2,3) ;this array will store the current absolute position at each iteration
dqac=dblarr(nx) ;this array keeps track of the signal lost on the AC electrode due to trapping
dqdc=dqac ;ditto for the DC electrode
qcloud=dblarr(nx,2) ;this is the charge remaining in the drifting clouds after charge trapping
dqcloud=dblarr(nx,2) ;this is the charge trapped for a given iteration
tdrift=dblarr(nx,2) ;drift time before collection on an electrode
 
;------------------------------
;SET INITIAL X, Y, Z POSITIONS
;------------------------------
x[*,0,0]=x0 & x[*,1,0]=x0
x[*,0,1]=y0 & x[*,1,1]=y0
x[*,0,2]=z0 & x[*,1,2]=z0

;------------------------------
;SET INITIAL CHARGES
;------------------------------
qcloud[*,0]=q0 & qcloud[*,1]=q0
dqac[*]=double(0.)
dqdc[*]=double(0.)

;------------------------------
;SET INITIAL DRIFT TIME
;------------------------------
tdrift[*,*]=double(0.)

;-----------------------
;RECORD INFO FOR TESTING
;-----------------------

;SEB 2023: THIS STRUCTURE BELOW IS NOT CURRENTLY IMPLEMENTED
;RETAINED FOR FUTURE IMPLEMENTATION
;if keyword_set(save_path) then begin ;save full path for testing
;  follow=where(x[*,0,2] gt 0,nfollow) ;whose path is needed? 
;  paths=dblarr(nfollow,numits+1,2,3) ;store x,y,z pos at each iter
;  paths[*,0,*,*]=x[follow,*,*] 
;  ef=dblarr(nfollow,numits+1,2,3) ;store electric field Ex, Ey, Ez at each iter
;  wfac=dblarr(nfollow,numits+1,2) ;AC electrode weighting field 
;  wfdc=wfac ;DC electrode weighting field
;  tx1=ef  
;  IF keyword_set(multich) THEN wfac2=wfac & wfdc2=wfdc;second AC & DC weighting field
;  print,'Following ',strtrim(nfollow,2),' charges...'
;endif else begin
;  print,'no charges followed...'
;endelse

;------------------
;DEFINE FIELD GRIDS
;------------------
;These will hold the interpolated field values at each iteration 
efield=dblarr(nx,2,3,/nozero) ;vector electric field
wfieldac=dblarr(nx,2,/nozero) ;scalar weighting field for AC electrode
wfielddc=dblarr(nx,2,/nozero) ;scalar weighting field for DC electrode

if keyword_set(multich) then begin
  wfieldac2=wfieldac ;store a second wfield array if multiple channels are active
  wfielddc2=wfielddc
  qac2=qac ;charge array for second AC electrode channel
  qdc2=qdc ;second DC electrode channel
endif

;-------------------------------------
;BEGIN MAIN CHARGE TRANSPORT ITERATION
;-------------------------------------
k=0L	; Iteration counter, retained to avoid infinite looping.
repeat begin
  k=k+1L	; Increment loop counter.
  dxabs=0b	; Free unused memory from previous iteration.

;-----------------------------------------------------------------------------
;Find the electric vector and scalar weighting field at the current positions.
;-----------------------------------------------------------------------------

  for dir=0,2 do begin
    efield[*,*,dir] = interpolate(egrid[*,*,*,dir],	x[*,*,0]/DX,$
							x[*,*,1]/DY,$
							x[*,*,2]/DZA )
  endfor;dir=0,2						
  wfieldac[*,*] = interpolate(wgridac[*,*,*],x[*,*,0]/DX,$
					 x[*,*,1]/DY,$
					 x[*,*,2]/DZA )
	wfielddc[*,*] = interpolate(wgriddc[*,*,*],x[*,*,0]/DX,$
				  x[*,*,1]/DY,$
				  x[*,*,2]/DZA )

;Consider cross-talk effect
  if keyword_set(multich) then begin ;wfields for neighboring electrodes
      wfieldac2[*,*] = interpolate(wgridac[*,*,*],(x[*,*,0]/DX - double(ISTRIP)),$
						x[*,*,1]/DY,$
						x[*,*,2]/DZA )
			wfielddc2[*,*] = interpolate(wgriddc[*,*,*],x[*,*,0]/DX,$
				  (x[*,*,1]/DY - double(ISTRIP)),$
				  x[*,*,2]/DZA )
				  
    ;Our second channel is to the right in X (on the AC electrode side), 
    ;hence its weighting field is the same as wgridac for a position 
    ;ISTRIP*DX microns to the left in X.
    ;SEB 2023: Not sure this makes sense as presented
  endif

;-----------------------
;CHECK FOR OUT OF BOUNDS
;-----------------------
; Null the fields in the x direction where the position is out of bound.
; Note that the keyword `missing' for interpolate() cannot be used,
; because we want to null each field only in one direction.
  within_bound = (x[*,*,0] ge 0.0d) and (x[*,*,0] le DX*(IMAX-1))
  efield[*,*,0]=efield[*,*,0]*within_bound
  within_bound=0b	; Free memory.
; Ditto, in the y direction.
  within_bound = (x[*,*,1] ge 0.0d) and (x[*,*,1] le DY*(JMAX-1))
  efield[*,*,1]=efield[*,*,1]*within_bound
; In the z direction, null the electric field where the position is either
; out of bound or on a surface; the weighting field is unchanged.
  efield[*,*,2]=efield[*,*,2]*((x[*,*,2] gt 0.0d)$
				   and (x[*,*,2] lt HEIGHT))

;-------------------------------------------------
; Find the movement induced by the electric field
;-------------------------------------------------
	txstep=drift_dir(efield)
  dxabs=sqrt(total(txstep^2,3))

;-------------------------------------------------
; If the charge is still moving update drift time
;-------------------------------------------------

tdrift += DTSTEP*(dxabs ge dxstep_tolerance)

;---------------------------------
;CALCULATE INDUCED CHARGE FOR STEP
;---------------------------------
  if keyword_set(trap) then begin 
    qac(*,k)=qcloud(*,0)*wfieldac(*,0)-qcloud(*,1)*wfieldac(*,1)+dqac(*) ;dqac includes the AC signals from charges trapped in previous iterations
    qdc(*,k)=qcloud(*,0)*wfielddc(*,0)-qcloud(*,1)*wfielddc(*,1)+dqdc(*) ;ditto for dqdc and the DC signals
;    dqcloud(*,0)=qcloud(*,0)*(XTHERM(0)/XTRAP(0)) & dqcloud(*,1)=qcloud(*,1)*(XTHERM(1)/XTRAP(1)) ;the charged trapped in this iteration
;    dqcloud(*,0)=qcloud(*,0)*(dxabs(*,0)/XTRAP(0)) & dqcloud(*,1)=qcloud(*,1)*(dxabs(*,1)/XTRAP(1)) ;the charged trapped in this iteration
    dqcloud(*,0)=qcloud(*,0)*(sqrt(dxabs(*,0)^2+XTHERM(0)^2)/XTRAP(0)) & dqcloud(*,1)=qcloud(*,1)*(sqrt(dxabs(*,1)^2+XTHERM(1)^2)/XTRAP(1)) ;the charged trapped in this iteration
    qcloud(*,0)-=dqcloud(*,0) & qcloud(*,1)-=dqcloud(*,1) ;the charge cloud remaining untrapped
    dqac(*)+=dqcloud(*,0)*wfieldac(*,0)-dqcloud(*,1)*wfieldac(*,1) ;the AC signals from the charges trapped at this + previous iterations
    dqdc(*)+=dqcloud(*,0)*wfielddc(*,0)-dqcloud(*,1)*wfielddc(*,1) ;ditto for the DC signals
  endif else begin
    qac(*,k)=qcloud(*,0)*wfieldac(*,0)-qcloud(*,1)*wfieldac(*,1)
    qdc(*,k)=qcloud(*,0)*wfielddc(*,0)-qcloud(*,1)*wfielddc(*,1)
  endelse
 
  ;Consider cross-talk effect
  if keyword_set(multich) then begin
    qac2(*,k)=q0*(wfieldac2(*,0)-wfieldac2(*,1))
    qdc2(*,k)=q0*(wfielddc2(*,0)-wfielddc2(*,1))
  endif
  
;----------------------------
;CHECK FOR BOUNDARY OVERSHO0T
;----------------------------
; Check that the boundary wasn't overshot, and correct if it is.
; testing change
  tx=x+txstep
  too_high = (tx[*,*,2] ge HEIGHT) ;e- hit AC electrode
  too_low  = (tx[*,*,2] le 0.0d) ;hole hit DC electrode

;-------------------------------
;MOVE CHARGES ALONG EFIELD LINES
;-------------------------------
; Now move the charges and ...
  blas_axpy,x,1.0d,txstep
  x[*,*,2]=x[*,*,2]*(1b-too_low-too_high) + HEIGHT*too_high ; + 0.0*too_low (implied)

;----------------------------------------
;TEST FOR TOTAL COLLECTION OF ALL CHARGES
;----------------------------------------
  txstep=0b                    
endrep until (k ge NUMITS) ;SEB 5/25/05

;------------------------------------
;END MAIN CHARGE TRANSPORT ITERATION
;------------------------------------

if (k gt NUMITS)$
 then message,'Number of iterations exceeds expected value! Forced exit.',$
	/continue
 
;-------------------------------------
;COMBINE WAVEFORMS INTO SINGLE EVENTS
;-------------------------------------
;combine waveforms if an evt occurs at more than 1 site (i.e. is compton scattered)
if keyword_set(numsites) then begin
  combine_waveforms,qac,qdc,numsites,index_keep
  nx0=nx
  nx=n_elements(qac[*,0])
  print,'Combined waveforms ... ',nx0,' reduced to ',nx
endif

;------------------------
;CHANGE Q UNITS TO keV!!!
;------------------------
qac = float (temporary(qac[*,0:k-1])/en2q)   ;to save storage space
qdc = float (temporary(qdc[*,0:k-1])/en2q)
IF keyword_set(multich) THEN BEGIN
  qac2 = float (temporary(qac2[*,0:k-1])/en2q)
  qdc2 = float (temporary(qdc2[*,0:k-1])/en2q)
ENDIF

;------------------------
;SAVE DRIFT TIMES in [ns]
;------------------------
tdrift=tdrift*(1.e9)

;-----------------------------------------
;RETURN STRUCTURE OF TOTAL INDUCED CHARGE
;-----------------------------------------
case 1 of
  keyword_set(nopulse):return,{tac:tdrift[*,0], tdc:tdrift[*,1],$ ;drift times in [ns]
                               eac:qac[*,k-1], edc:(-1.)*qdc[*,k-1],$ ;energy collected in [keV]
                               x0:x0,$ ;starting position
                               y0:y0,$ ;starting energy
                               z0:z0,$
                               e0:fltarr(nx)}
  keyword_set(multich): return,{ac:qac[*,0:k-1], dc:qdc[*,0:k-1],$
                                ac2:qac2[*,0:k-1],$
                                dc2:qdc2[*,0:k-1],$
                                x0:x0,  $
                                y0:y0,$
                                z0:z0, $
                                e0:fltarr(nx)} 
  else: return,{ac:qac[*,0:k-1], dc:qdc[*,0:k-1],$ ;charge signals in [keV]
                tac:tdrift[*,0], tdc:tdrift[*,1],$ ;drift times in [ns]
                eac:qac[*,k-1], edc:(-1.)*qdc[*,k-1],$ ;energy collected in [keV]
                x0:x0,$ ;starting position
                y0:y0,$ ;starting energy
                z0:z0,$
                e0:fltarr(nx)}
endcase
END ;qmove


;+*****************************************************************************
; NAME: QTRANS
;
; PURPOSE: To simulate the charge transport and signal induction on an 
;          active AC electrode (X strip) and DC electrode (Y strip) in a 
;          Ge detector given input event parameters. 
;
; SYNTAX: qtrans,z0,x0,y0,e0,q,outfile=outfile,
;           /trap, /ac_irradiate,
;          (through INITIAL_POS): EXPDIST=NUM OR 1, /RANDOM,  
;               EVENT_ATTEN= # in cm^-1 
;
; INPUTS: 
;   z0 - input array of NUM initial Z (depth) positions (0-HEIGHT cm)
;         *if undefined, /EXPDIST is automatically set 	
;   x0 - input array of NUM initial X positions (0-0.8 cm)
;         *if undefined, /RANDOM is set for x0 ONLY.
;   y0 - input array of NUM initial Y positions (0-0.8 cm)
;         *if undefined, /RANDOM is set for y0 ONLY.
;   e0 - input array of NUM initial energies  (keV) 
;  *set e0 to a single energy if the same energy is 
;  desired for all positions
;  ** if undefined,  Default = 122 keV (NOTE: this does NOT set the
;  attenuation for 122 keV! It only affects the final signal
;  amplitude.)
;
; NOTE: The total number of events NUM is defined by EXPDIST (if set and > 1),
; the array size of z0 (if defined), x0 & y0 (if defined and /RANDOM is not
; set), e0 (if defined and > 1 element), or the default = 2000 in that order. 
;    
; OPTIONAL INPUTS:    *KEYWORDS ALWAYS OVERRIDE INPUT ARRAYS*
;
; /trap - Set to include trapping - no trapping is default!
; /lambda - Set to define the effective trapping lengths
; /nopulse - Set to leave the times dependent charge signals (q.ac, q.dc)
;              out of the saved data file, default is to store the signals
; /ac_irradiate - Set this to transform z0 => HEIGHT - z0. If z0 is an
;                   exponential distribution, this places the majority of
;                   events near z0= HEIGHT cm, as if the AC side of the
;                   detector were irradiated. The default is DC
;                   irradiation which places the majority of events near
;                   z0=0.0 cm.
; 
; (from INITIAL_POS procedure)
; EXPDIST=NUM - Set to choose NUM random z0 positions from an exponential
;               distribution. Default is Cathode (AC) irradiation of the
;               detector which places the majority of events near z0=0.0 cm. 
;
; NOTE:
;            -EXPDIST overwrites values (if present) in z0
;            -Set EXPDIST=1 (same as /EXPDIST) if NUM is already 
;             defined by e0, x0, or y0. (NUM set by EXPDIST will 
;             dominate if [SA1]not equal to 1)
;            -Attenuation of exponential is taken from a single e0
;             value of 60, 122, or 662 keV (hard-coded) OR keyword
;             EVENT_ATTEN (see below).
;            -**if no attenuation value is available from above 
;             sources, then NUM random z0 positions are chosen from a FLAT
;             distribution (Default setting of e0 does NOT set the
;             attenuation!! - to use 122 keV atten, set e0=[122] ) 
;
; /RANDOM  -  Set to choose NUM random values for BOTH x0 and y0 from 
;             flat distributions spanning the mid-gap to mid-gap region
;             around the central strip (0.3 - 0.5 cm).
;
; NOTE: 
;            -Random values overwrite elements (if present) in x0 and y0.
;            -The number of elements initially in x0 & y0 is nulled when
;             /RANDOM is set. Thus, the input number of events (NUM) will
;             be taken from EXPDIST (if set), z0 (if defined), e0 (if
;             defined), or default (2000) in that order. 
;            -To randomly set only ONE of x0 or y0, pass in an undefined
;             array for that parameter rather than using /RANDOM manually.
;
; EVENT_ATTEN = ATTEN - Set to ATTEN in cm^-1 to choose Z-positions from
;             an exponential with this attenuation. Applies ONLY if
;             /EXPDIST is set. 
;
; NOTE:
;             -Overrides value determined by e0 (if present). 
;
; OUTPUT: q - a structure with TAGS described below. 

;      n_evt = total number of input events
;      n_time_steps = total number of time steps before the last charge 
;                     hit a detector edge
;
;      The tags of q are:     
;      {
;       AC:  Array[n_evt, n_time_steps] containing the time resolved signal
;             induced on the AC electrode in units of keV. 
;       DC:   Array[n_evt, n_time_steps] containing the time resolved signal
;              induced on the DC electrode in units of keV.
;       TAC: Array[n_evt] containing the drift time for the AC signal [ns]
;       TDC: Array[n_evt] containing the drift time for the DC signal [ns]
;       EAC: Array[n_evt] energy collected for the AC signal [keV]
;       EDC: Array[n_evt] energy collected for the DC signal [keV]
;       X0: Array[n_evt] containing the initial X position of the event (cm).
;       Y0: Array[n_evt] containing the in:qitial Y position of the event (cm).
;       Z0: Array[n_evt] containing the initial Z position of the event (cm).
;       E0: Array[n_evt] containing the initial energy of the event in keV.
;       }
; 
;    - q can be output as a named variable 'q' in the QTRANS call, 
;      or saved in a file. Use OUTFILE = 'outfilename.sav' to save as an 
;      IDL save file or OUTFILE = 'outfilename.fits' to save as a FITS 
;      file (recommended). Files are saved to the hard-coded directory
;      OUTFILE_DIRECTORY = defined in the COMMON BLOCKS.
;
; HISTORY: Based on QTRANS written by S. Boggs (see QTRANS.PRO for more)
;          Mod by S. Amrose UCB/SSL   2002
;          Update by Mark Bandstra UCB/SSL 2005
;          Revised by SEB UCSD 2023
;***************************************************************************** 
PRO qtrans,z0,x0,y0,e0,q,outfile=outfile,$
               trap=trap,lambda=lambda,nopulse=nopulse,$
               ac_irradiate=ac_irradiate,$
               _ref_extra=ex
               ;through INITIAL_POS: EXPDIST=NUM OR 1, /RANDOM,  
               ;EVENT_ATTEN= # in cm^-1 

IF n_params() EQ 0 THEN BEGIN
  print,'syntax- qtrans,z0,x0,y0,e0 (keV),q,outfile=outfile,trap=trap,lambda=lambda,z0_expdist=z0_expdist,event_energy=event_energy,ac_irradiate=ac_irradiate,_ref_extra=ex'
  return
ENDIF
;SEB 2023: event_energy is no longer defined. Restore it?

!except=1
start_time = systime(1)
printf,-2,'Start time: ',systime(0)

;------------------------
;Define the COMMON BLOCKS
;------------------------
define_common_blks	; This defines all numeric parameters.
;*Note that you must run define_common_blks before successfully compiling qtrans!
COMMON FLDGRID
COMMON PIXEL
COMMON DETECTOR
COMMON DIR ;outfile_dir
field_dir=outfile_dir
;efieldfile='efield_det4_15strip.sav' ;Scalar electric field
;efieldfile='efield_asic.sav' ;Scalar electric field
;efieldfile='efield_Samer500V_15strip.sav' ;Scalar electric field

;The files below were for the spare detector, but with a height of 15.0mm
;efieldfile='efield_asic_15strip.sav' ;Scalar electric field
;wfieldfileac='wfieldac_15strip.sav' ;Use scalar AC wfield to calculate dQ (old: vector)
;wfieldfiledc='wfielddc_15strip.sav' ;Use scalar DC wfield

;The files below were modified for 15.1mm (SEB 2/7/24)
;efieldfile='efield_HP414183.sav' ;Scalar electric field
;wfieldfileac='wfieldac_HP414183.sav' ;Use scalar AC wfield to calculate dQ (old: vector)
;wfieldfiledc='wfielddc_HP414183.sav' ;Use scalar DC wfield

;Modified for any detector ID (SEB 2/8/24)
efieldfile='efield_'+det_name+'.sav' ;Scalar electric field
wfieldfileac='wfieldac_'+det_name+'.sav' ;Use scalar AC wfield to calculate dQ (old: vector)
wfieldfiledc='wfielddc_'+det_name+'.sav' ;Use scalar DC wfield


;The files below were for a pure planar 15.0mm
;efieldfile='eplanar_600V.sav' ;Scalar electric field
;wfieldfileac='wplanar_ac.sav' ;Use scalar AC wfield to calculate dQ (old: vector)
;wfieldfiledc='wplanar_dc.sav' ;Use scalar DC wfield


if not keyword_set(outfile) then outfile = 'outfilename.sav'

;Trapping, default is No Trapping
trap=keyword_set(trap)? 1:0

XTRAP=[1000.0d,$  ; Electron effective trapping length [cm]
  1000.0d ]  ; Hole     effective trapping length [cm]
if keyword_set(lambda) then XTRAP=lambda
print,'XTRAP',XTRAP

;Save pulses, default is save them
nopulse=keyword_set(nopulse)? 1:0

;------------------------------
;RESTORE VECTOR & SCALAR FIELDS
;------------------------------
printf,-2,format='("Restoring vector and scalar fields...",$)'
restore,field_dir+efieldfile ;scalar electric field
print,'restored ',field_dir+efieldfile
egrid=vectgrid(eugrid) ;convert to vector electric field
eugrid=0

restore,field_dir+wfieldfileac
print,'restored ',field_dir+wfieldfileac
wgridac=eugrid[*,*,KAIR:KMAX-KAIR-1] ;AC scalar, remove Air sections
eugrid=0

restore,field_dir+wfieldfiledc
print,'restored ',field_dir+wfieldfiledc
wgriddc=eugrid[*,*,KAIR:KMAX-KAIR-1] ;DC scalar, remove Air sections
eugrid=0

;----------------------------
;SORT OUT INITIAL POSITIONS
;----------------------------

initial_pos,z0,x0,y0,e0,q0,num,atten,_extra=ex ;KEYWORDS: EXPDIST=NUM OR 1, /RANDOM,  
                                               ;EVENT_ATTEN= # in cm^-1 
print,'Total Objects: ',strtrim(num,2)

;-----------------------------
;AC or DC IRRADIATION?
;-----------------------------
if keyword_set(ac_irradiate) then BEGIN
  z0=HEIGHT-z0 ;zmax =1.5(cm)
  print,'AC Irradiate'
ENDIF ELSE IF keyword_set(expdist) THEN print,'DC Irradiate'

printf,-2,format=$
 '("Transporting charges/tracing field lines (will take forever)...",$)'

;------------------
;TRANSPORT CHARGES
;------------------
q=qmove(x0,y0,z0,q0,egrid,wgridac,wgriddc,trap=trap,XTRAP,nopulse=nopulse,_extra=ex)
; This is the main procedure for charge transport.

egrid=0b ;clear memory
wgridac=0b
wgriddc=0b

q.e0=e0 ;SEB 2023: set here because e0 is not sent to qmove

printf,-2,'DONE!!'
;------------------------------------------
;RECORD OUTPUT INTO A FITSFILE OR .SAV FILE
;------------------------------------------
if keyword_set(outfile) then BEGIN
  IF (strpos(outfile,'.fits'))[0] NE -1 THEN BEGIN
    mwrfits,q,outfile_dir+outfile,/create 
    outstyle='fits format'
  ENDIF ELSE BEGIN
    save,q,filename=outfile_dir+outfile
    outstyle='IDL Save format'
  ENDELSE

  print,'Saved outfile in '+outstyle+' to: ',outfile_dir+outfile
endif

end_time = systime(1)
printf,-2,'Time taken=',(end_time-start_time)/60.,' minutes.'

END ;qtrans

