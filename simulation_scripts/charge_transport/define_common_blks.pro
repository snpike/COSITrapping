; Current revision: 2.0 03/15/2023 SEB

pro define_common_blks

; In both field3d.rpo and qtrans.pro, lengths are defined in [cm], as cgs formulae are used.

;SEB: The guidlelines below were for the old pixel simulations. Update these for the strips. 
; Guidelines to defining the grid dimensions
; -To take advantage of the reflective symmetry of the pixel plane, one should
;  always model only one quadrant of it.
;
; -As it is implemented currently, the true boundaries of the grid are i=-0.5,
;  i=IMAX-0.5, j=-0.5 and j=JMAX-0.5.  The modelled geometry is assumed to
;  be reflectively symmetric about these lines.  This implies that features
;  positioned at grid boundaries must each be an odd number of grid cells
;  wide, with the central grid cell right outside the outermost grid point;
;  eg, a 5x5 pixel contact placed at the origin will occupy the cells within
;  i < 2 and j < 2.
;
; -Metal-semiconductor boundaries must be placed along grid points, and their
;  potentials should be set fixed to that of the metal;
;  consequently, metallic features must be at least 1 cell wide.
;
; -Due to the way field gradients are calculated in procedure vectgrid in
;  qtrans.pro, semiconductor features must be bound by at least 5 rows/columns
;  of grid points (ie, at least 4 cells wide if it is not at a grid boundary,
;  and 9 cells wide if it is, with the 5th cell at the boundary).
;  These numbers can be reduced to 4, 3 and 7 respectively, but the
;  computation will be quite a bit more complicated.

;(This is the xy-plane at k=KJUMP)
;
;   |<--      IPIX       -->|		"x" grid points
; -x-x-x-x-x-x-x-x-x-x-x-x-x-x-		Pixel boundaries should be at centres
;  |_|   |_____________|   |_|		 of grid cells, while other boundaries
;    |<->|<--ICNTCT -->| ->| |<-	 should be at grid points.
;    IGAP                 ISTEER

;SEB (2/8/24): This is hardcoded for "15 strips", i.e., 13 full strips and 2 half strips.
;It is also hardcoded that the contacts are 7/8 strip pitch, gap is 1/8 strip pitch. 

;DETECTOR SPECIFICATIONS
det_ID='HP415201' ;Unique detector identifier
det_HEIGHT=double(1.54) ;Crystal thickness
det_IMPURITY=6.0e9 ;Impurity concentration [cm-3] in germanium, p-type [-], n-type [+]
det_HV=1200. ;Operational HV setting
det_PITCH= double(0.20);Strip pitch [cm]
;SEB: Create this folder and document the specs above before running fields.



common PIXEL,ICNTCT,ISTRIP
  ; used in all user-level procedures
  ;  (field3d, capcalc, qtrans, qshare, plotfield, and their subroutines)
 
ICNTCT=29 ;number of grid points in a contact strip (including both endpoints)
ISTRIP=32 ;number of grid points in a strip/gap pair
;SEB  put points along both the gap and contact symmetry axis

common FLDGRID,IMAX,JMAX,KMAX,KAIR,DX,DY,DZA,DZB,HEIGHT,NSTRIP
 ; used in field3d, capcalc, qtrans and plotfield (and their subroutines)

;Define the size of the GRID                 
IMAX=449 ;(14*ISTRIP+1)# of x gridpoints        
JMAX=449 ;# of y gridpoints          
KMAX=161 ;# of z gridpoints          
KAIR=5 ;size of air buffer top and bottom

NSTRIP=fix((IMAX-1)/ISTRIP) ;ISTRIP defined above         
DX=double(NSTRIP*det_PITCH)/float(IMAX-1) ;grid spacing in [cm]
DY=double(NSTRIP*det_PITCH)/float(IMAX-1);double(8.e-1)/float(JMAX-1) ;grid spacing in [cm] 
HEIGHT=det_HEIGHT   ;depth of Ge in [cm]  
DZA=HEIGHT/float(kmax-2*kair-1)  ;grid spacing in [cm] (semiconductor)          
DZB=double(1.0)/float(kair)  ;grid spacing in [cm] (air)  
 

; LAPLACE'S EQUATION AND BOUNDARY CONDITIONS
;  The variable pnt(i,j,k) defines the nature of the point: either contact
;  or type of material.
;      1 = Ground
;      2 = High Voltage (at volt, above)
;      4 = Air
;      0 = Germanium, CdZnTe or other semiconductor
;      5 = Steering electrode voltage
;  The variable pnt2(i,j,k) is different.  If pnt2(i,j,k) is 0, then the pot-
;  ential of the point is fixed.  If pnt2(i,j,k) is 1, then the potential can
;  be varied during the relaxation. */
;  Note: It used to be (fixed,variable)=(1,0).  Hubert changed it
;  on 1- 5-2000 in order to facilitate simple matrix operation.
;  See the main loop in field3d.

common TYPES,GND,HV,Air,Semiconductor,Fixed,Variable
  ; used in field3d and qtrans (and their subroutines)
GND=1 & HV=2 & Air=4 & Semiconductor=0
Fixed=0b & Variable=1b
  
common STEP,DTSTEP,NUMITS,dxstep_tolerance,efield_tolerance
; NUMITS and dxstep_tolerance are used in qmove(), and
; DXSTEP and efield_tolerance in drift_dir(); both defined in qtrans.pro.
; Neither is this common block necessary in terms of functionality,
; but these variables are all user-defined parameters.

DTSTEP=0.50e-9 ;[s], the constant time step for tracing charges = 0.5 ns
MAXTIME=600.e-9 ;[s] 600ns, lenghtened from 300ns to account for shaping electronics
NUMITS=long(MAXTIME/DTSTEP)

; The `tolerance' level for errors accumulated from numerical calculations;
; they are needed because numerical calculations seldom yield an exact zero
; when the true answer is zero.
 ; SEB: Note that neither of the following two parameters are still implemented.
dxstep_tolerance=1.e-5 ; should be less than DXSTEP, obviously.
efield_tolerance=1.0e-20


common DETECTOR,EPS0,EPS,HIGH_VOLTAGE,IMPURITY,QELEC,CNDBULK,CNDSURF
  ; All but the last variable are used in field3d,
  ; only EPS0 and EPS are used in capcalc, while
  ; the last variable is used (not implemented) in qmove() in qtrans.pro.
  ; This common block is not really necessary in terms of functionality,
  ; but these variables are all properties of the detector that may change.
  
EPS0=8.854188e-14 ;permittivity of free space [F/cm]
EPS=16.0*EPS0 ;permittivity of germanium [F/cm]
HIGH_VOLTAGE= det_HV ;Operational HV setting
IMPURITY= det_IMPURITY ;impurity concentration [cm-3] in germanium, p-type [-], n-type [+]
QELEC=1.602177e-19 ;electron charge [Coulomb]
CNDBULK=5.0E-12 ;Bulk conductivity [ohm-1 cm-1] (!!NOT RESISTIVITY!!)
CNDSURF=1.0E-16 ;Surface conductivity [ohm-1] (no change for cgs)


common DRIFT,e_e,vd_e,e_h,vd_h
  ; These are the experimental drift velocity data used in DRIFT_DIR
  ; assuming 77K and deflate <100> crystal orientation
  ; electrons: PhysRevB.24.1014
  ; holes: PhysRevB.16.2781

e_e=[0.,1.584349e+1, 1.920242e+1, 2.249333e+1, 2.654293e+1, 3.086379e+1, 3.615322e+1, 4.234916e+1, $
     4.960696e+1, 5.810859e+1, 6.806724e+1, 7.973260e+1, 9.447783e+1, 1.076033e+2, 1.263235e+2, $
     1.479728e+2, 1.733323e+2, 2.030380e+2, 2.378346e+2, 2.785947e+2, 3.263403e+2, 3.822684e+2, $
     4.477816e+2, 5.245223e+2, 6.144149e+2, 7.197133e+2, 8.430577e+2, 9.875409e+2, 1.156786e+3, $
     1.355035e+3, 1.587261e+3, 1.859286e+3, 2.177930e+3, 2.551184e+3, 2.988405e+3, 3.500558e+3, $
     4.100483e+3, 4.803224e+3, 5.626400e+3, 6.590652e+3, 7.720157e+3, 9.043237e+3, 1.007310e+4] ;V/cm

vd_e=[0.,6.523642e+5, 7.780541e+5, 8.982910e+5, 1.033836e+6, 1.180478e+6, 1.351988e+6, 1.576981e+6, $
      1.815190e+6, 2.054703e+6, 2.324969e+6, 2.634788e+6, 2.987068e+6, 3.271897e+6, 3.666246e+6, $
      4.112519e+6, 4.591129e+6, 5.062380e+6, 5.577113e+6, 6.092759e+6, 6.679841e+6, 7.347673e+6, $
      7.943602e+6, 8.490018e+6, 9.077784e+6, 9.681839e+6, 1.016817e+7, 1.064584e+7, 1.117514e+7, $
      1.155616e+7, 1.192045e+7, 1.226121e+7, 1.239846e+7, 1.265728e+7, 1.268675e+7, 1.264626e+7, $
      1.257227e+7, 1.231646e+7, 1.213048e+7, 1.195528e+7, 1.173763e+7, 1.146813e+7, 1.124086e+7] ;cm/s

e_h=[0.,2.059343e+1, 2.311964e+1, 2.497177e+1, 2.730365e+1, 2.980932e+1, 3.211244e+1, 3.549757e+1, $
     4.039217e+1, 4.597356e+1, 5.203879e+1, 5.978833e+1, 6.632609e+1, 7.813868e+1, 8.947440e+1, $
     1.031570e+2, 1.209779e+2, 1.426057e+2, 1.649539e+2, 1.937726e+2, 2.206406e+2, 2.607205e+2, $
     3.100220e+2, 3.732890e+2, 4.497951e+2, 5.316677e+2, 6.344583e+2, 7.513284e+2, 8.828507e+2, $
     1.048552e+3, 1.243468e+3, 1.465494e+3, 1.728485e+3, 2.043784e+3, 2.515696e+3, 2.964580e+3, $
     3.488190e+3, 4.221915e+3, 4.888399e+3, 6.059044e+3, 7.204958e+3, 8.392104e+3, 9.670185e+3] ;V/cm

vd_h=[0.,9.598331e+5, 1.082017e+6, 1.170087e+6, 1.284689e+6, 1.405773e+6, 1.510779e+6, 1.660925e+6, $
      1.865382e+6, 2.082993e+6, 2.289609e+6, 2.546816e+6, 2.744996e+6, 3.079084e+6, 3.360660e+6, $
      3.663116e+6, 4.009684e+6, 4.393719e+6, 4.728817e+6, 5.095278e+6, 5.400447e+6, 5.810310e+6, $
      6.167911e+6, 6.581091e+6, 6.964178e+6, 7.331810e+6, 7.682748e+6, 8.036594e+6, 8.334384e+6, $
      8.618901e+6, 8.926226e+6, 9.182170e+6, 9.443760e+6, 9.669208e+6, 9.980069e+6, 1.017121e+7, $
      1.036044e+7, 1.054591e+7, 1.067864e+7, 1.081339e+7, 1.085859e+7, 1.089705e+7, 1.089789e+7] ;cm/s

common TRAPPING, XTHERM
;XTRAP[] is the effective trapping length in germanium detectors, including thermal motions
;XTHERM[] is the effective path length due to thermal motion during a time step DTSTEP. 
;Note that the thermal velocities sqrt(3kT/m*) are calcuated at 80K, using the conductivity
;effective mass. m*=0.12me for electrons, m*=0.21me for holes. 
;SEB 2/8/24: XTRAP removed to allow passing it through the function calls.
;XTRAP=[1000.0d,$  ; Electron effective trapping length [cm]
;       1000.0d ]  ; Hole     effective trapping length [cm]
;XTRAP=[1000.0d,$  ; Electron effective trapping length [cm]
;       10000000.0d ]  ; Hole     effective trapping length [cm]
XTHERM=[(1.7e7)*DTSTEP,$  ; Electron thermal path length [cm]
        (1.3e7)*DTSTEP ]  ; Hole     thermal path length [cm]  
 
common FACTOR_CAL, fac_e, fac_h
;scale factors to align simulated pulse profiles with data
;SEB: This was implemented in drift_dir to account for discrepancies between simulated and 
;measured CTDs. Hopefully not needed with implementation of the new drift velocities. 
fac_e=1.0 ;electron stretch factor
fac_h=1.0 ;hole stretch factor

COMMON DIR,  outfile_dir,det_name
;Make sure this det_ID folder is created before running field calculations.
outfile_dir='~/Insync/seboggs@physics.ucsd.edu/Google Drive/germanium/detectors/'+det_ID+'/'
det_name=det_ID
end
