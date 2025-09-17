; Jovian magnetosphere/ionosphere mapping program - 2025 version 1
; Written by Marissa Vogt, mvogt@psi.edu
; Last updated September 2025
;
; PURPOSE OF THIS FUNCTION
; This function allows a user to magnetically map a point between Jupiter's middle magnetosphere and the ionosphere. 
; Points can be mapped using the flux equivalence mapping model of Vogt et al. [2011, 2015] or by tracing field lines from a global field model. 
; 
; HOW TO USE THIS FUNCTION
; (Note that this function maps one point at a time but user can call it in a loop to map several points at once.)
; **Before running this program the user should replace the directory_name variable on line 524 below with the appropriate directory for their local machine
;
;   USER INPUTS - REQUIRED
;   * mapping_type: STRING that specifies the type of mapping to do -- from the magnetosphere to the ionosphere, or 
;       from the ionosphere (north or south) to the magnetosphere
;       Valid inputs are "ion_to_mag", "mag_to_ion_north", and "mag_to_ion_south"
;   * sslong: the subsolar longitude, in degrees SIII left-handed, 0-360
;   * input_position1 and input_position2: Coordinates in the magnetosphere to be mapped to the ionosphere OR 
;       coordinates in the ionosphere to be mapped to the magnetosphere.
;       For mapping from the magnetosphere to the ionosphere, input_position1 should be radial distance and input_position2 should be local time. 
;       Radial distance should be given in Jovian radii and range from 15 to 150 Rj (flux mapping) or 2-85 Rj (field line tracing). Local time should be given in decimal format (e.g. 12.5 for 12:30) and range from 0 to 24 hours. 
;       For mapping from the ionosphere to the magnetosphere, input_position1 should be latitude and input_position2 should be longitude.
;       Latitude should be given in degrees, and should be be negative for southern hemisphere points (range -90 to +90 degrees)
;       Longitude should also be given in degrees, in SIII left-handed, and range is 0-360 degrees. 
;
;   OPTIONAL KEYWORDS
;   "MODEL" - This keyword allows the user to choose which internal field model (VIP 4, the Grodent Anomaly Model, VIPAL, JRM09, or JRM33) to use in the mapping
;     The default value is now the JRM33 model (Connerney et al., 2018, GRL) for both north and south as of March 2024.
;         Note: The previous default was JRM09, and before that the default was the Grodent anomaly model in the north and VIP4 in the south, as in Vogt et al. [2011].
;     If the user wishes to use VIP4 (Connerney et al., 1998), for both north and south, the keyword MODEL should be set to 'vip4'
;     If the user wishes to use the Grodent anomaly model (Grodent et al., 2009) in the north, the keyword MODEL should be set to 'gam'
;     If the user wishes to use VIPAL (Hess et al., 2011), for both north and south, the keyword MODEL should be set to 'vipal'
;     If the user wishes to use JRM09 (Connerney et al., 2018), for both north and south, the keyword MODEL may be set to 'jrm09'
;     If the user wishes to use JRM33 (Connerney et al., 2022), for both north and south, the keyword MODEL may be set to 'jrm33'
;   "FIELDLINE_TRACING" - This keyword (added 2016) allows the user to calculate the M-I mapping by tracing fieldlines from a model rather than the Vogt et al. [2011] flux equivalence calculation
;     For fieldline tracing the user MUST specify a field model. The user can specify that the program should use traced field lines from:
;       * VIP4 (north or south) - set model keyword to "vip4"
;       * Grodent anomaly model (north) - set model keyword to "gam"
;       * VIPAL (north or south) - set model keyword to "vipal"
;       * the unpublished Khurana model (for both north and south) - set model keyword to "khurana" (see http://lasp.colorado.edu/home/mop/graphics/code/ for more info about the model)
;       * the unpublished Khurana model with JRM09 (instead of VIP4) - set model keyword to "khurana_jrm09"
;       * JRM09 (north or south) - set model keyword to "jrm09"
;       * JRM33 (north or south) - set model keyword to "jrm33"
;     The program calculates the mapping using field line tracing results that are pre-calculated (as for the flux equivalence mapping).
;     Note that for VIP4, VIPAL, JRM09, and JRM33 the model results are only valid to ~85 Rj. For the Grodent anomaly model results are only valid to ~65 Rj.  
;     All field line tracing results use the Connerney et al. (1981, 1998) current sheet except the Khurana model and Khurana_jrm09.
;
;
;
;   EXPECTED OUTPUT
;   This program outputs an array of length 2 containing the apppropriate
;     For mapping magnetosphere->northern/southern ionosphere, the result is an array [latitude,longitude] in degrees
;     For mapping ionosphere->magnetosphere, the result is an array [radial distance,local time] in Rj and hours (0-24)
;
;   The program returns array of flag values (-999, -998, -996, etc.) if the point cannot be mapped for some reason.
;   Reasons why a point might not be mapped:
;     * for mapping from the ionosphere to the magnetosphere, the ionospheric point is equatorward of the Ganymede footprint (15 Rj)
;     * for mapping from the magnetosphere to the ionosphere, the magnetospheric point is outside of the Joy et al (2002) magnetopause
;     * for mapping with fieldline tracing, if using VIP4, GAM, or VIPAL, point is outside 95 Rj
;     * the program has an error
;
;
; EXAMPLES - All using 180 degrees subsolar longitude
;   1. For mapping from the northern ionosphere to the magnetosphere using the Vogt et al. (2011, 2015) flux equivalence model, 
;   with input position 60 degrees latitude, 180 degrees longitude:
;     a) using the Vogt et al. flux equivalence calculation with the Grodent anomaly model:
;       result = mapping_function_2025('ion_to_mag',180.,60.,180.,'gam')
;       print,result
;       > result = [87.097375       10.276313]  ;;; means point maps to ~87 Rj and ~10.2 LT
;     b) using the Vogt et al. flux equivalence calculation with VIP4:
;       result = mapping_function_2025('ion_to_mag',180.,60.,180.,'vip4')
;       print,result
;       > result = [58.391824       10.110444]  ;;; means point maps to ~58 Rj and ~10.0 LT 
;     c) using the Vogt et al. flux equivalence calculation with JRM33:
;       result = mapping_function_2025('ion_to_mag',180.,60.,180.,'jrm33')
;       print,result
;       > result = [80.663833578430086       10.219791238872102]  ;;; means point maps to ~81 Rj and ~10.2 LT
;       Note that you will get the same result if you call result = mapping_function_2025('ion_to_mag',180.,60.,180.) since JRM33 is the default.
;
;   2. For mapping from the southern ionosphere to the magnetosphere using the Vogt et al. (2011, 2015) flux equivalence model, 
;     with input position -80 degrees latitude, 90 degrees longitude:
;     a) using the Vogt et al. flux equivalence calculation with VIP4:
;       result = mapping_function_2025('ion_to_mag',180.,-80.,90.,'vip4')
;       print,result
;       > result = [81.154380716253343       13.894088008926612]  ;;; means point maps to ~81 Rj and ~13.9 LT
;     b) using the Vogt et al. flux equivalence calculation with VIPAL:
;       result = mapping_function_2025('ion_to_mag',180.,-80.,90.,'vipal')
;       print,result
;       > result = [84.538550971174530       14.676701118754295]  ;;; means point maps to ~84 Rj and ~14.7 LT
;
;   3. For mapping from the magnetosphere to northern ionosphere using the Vogt et al. (2011, 2015) flux equivalence model,
;   with input position 88 Rj and 9.8 LT:
;     a) using the Vogt et al. flux equivalence calculation with the Grodent anomaly model:
;       result = mapping_function_2025('mag_to_ion_north',180.,88.,9.8,'gam')
;       print,result
;       > result = [60.172613419499804       181.95319557045138]  ;;; means point maps to ~60 degrees latitude and ~182 degrees longitude
;     b) using the Vogt et al. flux equivalence calculation with VIP4:
;       result = mapping_function_2025('mag_to_ion_north',180.,88.,9.8,'vip4')
;       print,result
;       > result = [62.564124582609807       179.15705854962323]  ;;; means point maps to ~63 degrees latitude and ~179 degrees longitude
;
;   4. For mapping from the magnetosphere to southern ionosphere using the Vogt et al. (2011, 2015) flux equivalence model,
;   with input position 88 Rj and 9.8 LT:
;     a) using the Vogt et al. flux equivalence calculation with VIPAL:
;       result = mapping_function_2025('mag_to_ion_south',180.,88.,9.8,'vipal')
;       print,result
;       > result = [-88.102481860795521       192.69359460926808]  ;;; means point maps to roughly -88 degrees latitude and ~193 degrees longitude
;     b) using the Vogt et al. flux equivalence calculation with JRM09:
;       result = mapping_function_2025('mag_to_ion_south',180.,88.,9.8,'jrm09')
;       print,result
;       > result = [-87.800807686083388       189.69085808246234]  ;;; means point maps to roughly -88 degrees latitude and ~191 degrees longitude
;     
;   5. For mapping from the ionosphere to the magnetosphere by tracing field lines from a global field model,
;   with input position 57 degrees latitude, 180 degrees longitude:
;     a) Using field lines traced with the Grodent anomaly model as the internal field model:
;       result = mapping_function_2025('ion_to_mag',180.,57.,180,'gam','fieldline_tracing')
;       print,result
;       > result = [36.673553276772196       11.065309982123598]  ;;; means point maps to ~37 Rj and ~11.1 LT
;     b) Using field lines traced from with JRM09 field model as the internal field model:
;       result = mapping_function_2025('ion_to_mag',180.,57.,180,'jrm33','fieldline_tracing')
;       print,result
;       > result = [23.217494240958338       10.990356069761237]  ;;; means point maps to ~23 Rj and ~11 LT
;       
;   6. For mapping from the magnetosphere to the northern ionosphere by tracing field lines from a global field model,
;   with input position 47 Rj and 18.2 LT:
;     a) Using field lines traced with VIP4 as the internal field model:
;       result = mapping_function_2025('mag_to_ion_north',180.,47.,18.2,'vip4','fieldline_tracing')
;       print,result
;       > result = [66.327301805996257       140.63710324913524]  ;;; means point maps to ~66 degrees latitude and ~141 degrees longitude
;     b) Using field lines traced with the Khurana external field model (with VIP4 as the internal field, default):
;       result = mapping_function_2025('mag_to_ion_north',180.,47.,18.2,'khurana','fieldline_tracing')
;       print,result
;       > result = [66.421892253475733       139.31335528325963]  ;;; means point maps to ~66 degrees latitude and ~139 degrees longitude
;     c) Using field lines traced with the Khurana external field model and JRM09 internal field:
;       result = mapping_function_2025('mag_to_ion_north',180.,47.,18.2,'khurana_jrm09','fieldline_tracing')
;       print,result
;       > result = [63.807225224238451       149.44916663382753]  ;;; means point maps to ~64 degrees latitude and ~149 degrees longitude
;
;   7. For mapping from the magnetosphere to the southern ionosphere by tracing field lines from a global field model,
;   with input position 47 Rj and 18.2 LT:
;     a) Using field lines traced with the VIPAL internal field model:
;       result = mapping_function_2025('mag_to_ion_south',180.,47.,18.2,'vipal','fieldline_tracing')
;       print,result
;       > result = [-69.919105658956227       68.309592929705900]  ;;; means point maps to roughly -70 degrees latitude, 68 degrees longitude
;     b) Using field lines traced with the JRM09 internal field model:
;       result = mapping_function_2025('mag_to_ion_south',180.,47.,18.2,'jrm09','fieldline_tracing')
;       print,result
;       > result = [-69.741396488857163       69.864958182304576]  ;;; means point maps to roughly -70 degrees latitude, 70 degrees longitude
;
; INPUT FILES AND PROGRAMS CALLED
;   Calls routines: atan2, point_inside_polygon
;
;   The program reads in data from files called sslong0_gam.sav, sslong10_jrm09.sav, etc., which should be placed in an appropriate directory on the user's machine.
;   These files can be downloaded from bit.ly/mappingfiles2024 (updated links will be posted on marissavogt.com/mapping and can be obtained by emailing mvogt@psi.edu)
;   ***Before running this program the user will have to replace the directory_name variable in the code (see line 524).
;
; ACKNOWLEDGMENTS
;   Thanks to Bertrand Bonfond, Benjamin Palmaerts, and Andrew Steffl for helping to test early versions of this code and the online mapping tool.
;   Very special thanks to Masafumi Imai for assistance implementing the JRM09 field model.
;   Thanks to Rob Wilson for assistance with code speed ups and other improvements. 
;
; WHEN USING RESULTS FROM THIS CODE IN A PAPER OR PRESENTATION, PLEASE CITE:
;   Vogt, M. F., M. G. Kivelson, K. K. Khurana, R. J. Walker, B. Bonfond, D. Grodent, and A. Radioti (2011), 
;   Improved mapping of Jupiter’s auroral features to magnetospheric sources, 
;   J. Geophys. Res., 116, A03220, doi:10.1029/2010JA016148.
;  and
;   Vogt, M. F., E. J. Bunce, M. G. Kivelson, K. K. Khurana, R. J. Walker, A. Radioti, B. Bonfond, and D. Grodent (2015), 
;   Magnetosphere-ionosphere mapping at Jupiter: Quantifying the effects of using different internal field models, 
;   J. Geophys. Res. Space Physics, 10.1002/2014JA020729
;
; Most recent changes:
;   (September 2025)
;     * Major restructuring of code to make it easier to compare IDL and Python versions
;   (January 2025)
;     * Various speed and accuracy improvements
;   (March 2024)
;     * Added JRM33 internal field option for flux mapping and fieldline tracing (with CAN current sheet)
;     * Changed default from JRM09 to JRM33
;   (January 2021) Early 2021 bug fixes:
;     * Fixed small error in local time mapping (~0.1 LT difference)
;   (September 2020)
;     * The fieldline tracing option is now valid for distances as close as 2 Rj (previously the inner limit was 15 Rj; the flux mapping is valid only 15-150 Rj)
;     * Speed improvements - the program now loads in the pre-calculated contour files as IDL save files rather than text files
;   (May 2020) I've recently made some minor changes to my CAN (Connerney et al. 1981) current sheet code, which affects all mapping options except
;       those using the Khurana model. The expected changes should be very minor (a few RJ/fraction of an hour LT in mapping
;       ionosphere->magnetosphere, or a fraction of a degree in lat/lon mapping magnetosphere->ionosphere).
;       All mapping options except JRM09 now use the CAN current sheet model with the VIP4 dipole values (tilt/longitude). The JRM09 flux mapping and
;       field line tracing options use the CAN sheet with the JRM09 dipole values.
;   (June 2019) Two changes:
;     * Fixed error in local time for fieldline tracing option
;     * Added the option to use the Khurana field model (kk2009) with JRM09 as the internal (planetary) field for field line tracing. The Khurana model normally uses VIP4 as the internal field.
;       For results using the Khurana field model with VIP4 (default) use "khurana".
;       For results using the Khurana field model with JRM09 (new) use "khurana_jrm09"
;   (March 2019) Several big changes in 2019 edition:
;     * Adds the JRM09 internal field model as an option for mapping using flux equivalence OR field line tracing (Connerney et al., GRL 2018)
;     * The mapping (flux equivalence or field line tracing) is now done for Jupiter's oblate, flattened, surface rather than to a constant radial distance
;       of 0.95 RJ. The Vogt et al. (2011, 2015) published models used a spherical surface with radius 0.95 RJ and estimated this introduced a ~1% error in the calculated flux at the surface.
;       This change does slightly shift the mapping results because the starting 15 RJ reference contour (Ganymede footprint) shifts slightly.
;     * Changes the way that the mapping function is called.
;     * Changes the default mapping from GAM (north) and VIP4 (south) to JRM09 for both hemispheres.
;     * Fixes an error that caused sharp discontinuities in mapping at certain subsolar longitudes
;       (e.g., the mapping at subsolar longitude 180 differed significantly from the mapping of the same point at subsolar longitude 181).
;       Thanks to Benjamin Palmaerts for finding this error.
;     * Other miscellaneous minor bug fixes.
;   (January 2017) Added Khurana unpublished global field model as option for mapping using field line tracing, small bug fixes
;   (September 2016) Added ability to map using (pre-calculated results of) field line tracing using a global field model. Default is still to use flux equivalence from the Vogt et al. [2011] model.
;   (August 2014) Error fixed in southern hemisphere mapping - not related to this routine but will affect results of using this routine
;   (10 Feb. 2014) Minor bug fixes to prevent error for SSlong 0-10
;   (4 Feb. 2014) Major fixes in how the interpolation is done for cases where (subsolar longitudes mod 10) != 0



FUNCTION identify_ion_to_mag_match,xpt,ypt,x1,x2,x3,x4,y1,y2,y3,y4

  ;Rotate ionospheric contour lat/lon positions into cartesian coordinates for easy interpolation
  ;define min/max of each of the four (x,y) points in each ionospheric box
  xmin = x1*0.d
  xmax = x1*0.d
  ymin = x1*0.d
  ymax = x1*0.d
  for j=0, n_elements(x1) - 1 do begin
    xmin(j) = min([x1(j), x2(j), x3(j), x4(j)]);
    xmax(j) = max([x1(j), x2(j), x3(j), x4(j)]);
    ymin(j) = min([y1(j), y2(j), y3(j), y4(j)]);
    ymax(j) = max([y1(j), y2(j), y3(j), y4(j)]);
  endfor
  ;quick test to make sure that the input ionospheric positions (xpt and ypt) fix inside any box
  filematch = where(xmin le xpt and xmax ge xpt and ymin le ypt and ymax ge ypt)
  insidept = 0
  ; now loop through the number of boxes that xpt and ypt fit into and do the point inside polygon test
  filematchpts = n_elements(filematch)
  if filematch(0) ne -1 then begin
    if filematchpts eq 1 then begin
      filematch = filematch(0)
      insidept = 1
    endif else begin
      for k=0, filematchpts-1 do begin
        insidept1 = point_inside_polygon(xpt,ypt,[x1(filematch(k)),x3(filematch(k)),x4(filematch(k)),x2(filematch(k))],[y1(filematch(k)),y3(filematch(k)),y4(filematch(k)),y2(filematch(k))])
        insidept2 = point_inside_polygon(xpt,ypt,[x1(filematch(k)),x4(filematch(k)),x3(filematch(k)),x2(filematch(k))],[y1(filematch(k)),y4(filematch(k)),y3(filematch(k)),y2(filematch(k))])
        insidept = max([insidept1,insidept2])
        if insidept eq 1 then begin
          filematch = filematch(k)
          k = k+10000
        endif
      endfor
    endelse
  endif

  if insidept eq 0 then filematch = -1
  return,filematch
  
END
 
 

 
 
 
 
FUNCTION mapping_interpolation_ion_to_mag,xpt,ypt,filematch,rstepsize,r1,ltime1,lat1,long1,r2,ltime2,lat2,long2,r3,ltime3,lat3,long3,r4,ltime4,lat4,long4,x1,x2,x3,x4,y1,y2,y3,y4
    
;    x1_match = x1(filematch);sin(lat1(filematch)*!dtor)*cos(long1(filematch)*!dtor);
;    y1_match = y1(filematch);sin(lat1(filematch)*!dtor)*sin(long1(filematch)*!dtor);
;    x2_match = x2(filematch);sin(lat2(filematch)*!dtor)*cos(long2(filematch)*!dtor);
;    y2_match = y2(filematch);sin(lat2(filematch)*!dtor)*sin(long2(filematch)*!dtor);
;    x3_match = x3(filematch);sin(lat3(filematch)*!dtor)*cos(long3(filematch)*!dtor);
;    y3_match = y3(filematch);sin(lat3(filematch)*!dtor)*sin(long3(filematch)*!dtor);
;    x4_match = x4(filematch);sin(lat4(filematch)*!dtor)*cos(long4(filematch)*!dtor);
;    y4_match = y4(filematch);sin(lat4(filematch)*!dtor)*sin(long4(filematch)*!dtor);
    
    x1prime = x1(filematch) - x1(filematch);
    x2prime = x2(filematch) - x1(filematch)
    x3prime = x3(filematch) - x1(filematch);
    x4prime = x4(filematch) - x1(filematch);
    xpt_prime = xpt - x1(filematch);
    y1prime = y1(filematch) - y1(filematch);
    y2prime = y2(filematch) - y1(filematch);
    y3prime = y3(filematch) - y1(filematch);
    y4prime = y4(filematch) - y1(filematch);
    ypt_prime = ypt - y1(filematch);
    
    theta = atan2(y3prime,x3prime);
    
    cos_theta = cos(theta)
    sin_theta = sin(theta)
    a1 = x1prime*cos_theta + y1prime*sin_theta;
    a2 = x2prime*cos_theta + y2prime*sin_theta;
    a3 = x3prime*cos_theta + y3prime*sin_theta;
    a4 = x4prime*cos_theta + y4prime*sin_theta;
    apt = xpt_prime*cos_theta + ypt_prime*sin_theta;
    b1 = -1.0d*x1prime*sin_theta + y1prime*cos_theta;
    b2 = -1.0d*x2prime*sin_theta + y2prime*cos_theta;
    b3 = -1.0d*x3prime*sin_theta + y3prime*cos_theta;
    b4 = -1.0d*x4prime*sin_theta + y4prime*cos_theta;
    bpt = -1.0d*xpt_prime*sin_theta + ypt_prime*cos_theta;
    
    c = b2 - ((b4-b2)/(a4-a2))*a2;
    b5 = ((b4-b2)/(a4-a2))*apt + c;
    r_interpol = r1(filematch) + bpt*rstepsize/b5;
    
    dist1 = (apt-a1)*(apt-a1) + (bpt-b1)*(bpt-b1);
    dist2 = (apt-a2)*(apt-a2) + (bpt-b2)*(bpt-b2);
    dist3 = (apt-a3)*(apt-a3) + (bpt-b3)*(bpt-b3);
    dist4 = (apt-a4)*(apt-a4) + (bpt-b4)*(bpt-b4);
    
    ltime1 = ltime1(filematch);
    ltime2 = ltime2(filematch);
    ltime3 = ltime3(filematch);
    ltime4 = ltime4(filematch);
      
    if ltime1 lt 0.d then ltime1 = ltime1 + 24.d
    if ltime2 lt 0.d then ltime2 = ltime2 + 24.d
    if ltime3 lt 0.d then ltime3 = ltime3 + 24.d
    if ltime4 lt 0.d then ltime4 = ltime4 + 24.d
  
    if (abs(ltime3 - ltime1) lt 12.0d) then begin
      ltime_interpol = ltime1 + ((ltime3 - ltime1)/(a3-a1))*(apt-a1);
    endif else begin
      ltime_interpol = ltime1 + ((24.0d -abs(ltime3 - ltime1))/(a3-a1))*(apt-a1);
    endelse
  
    ltime_interpol = ltime_interpol mod 24.d
  
  return,[r_interpol,ltime_interpol]

END

  






FUNCTION is_inside_joy_magnetopause,rj,loctime_rad
  pexpanded = 0.039d;
  p_root = pexpanded^(-.25d)
  amage = -0.134d + 0.488d*p_root;
  bmage = -0.581d - 0.225d*p_root;
  cmage = -0.186d - 0.016d*p_root;
  dmage = -0.014d + 0.096d*pexpanded;
  emage = -0.814d - 0.811d*pexpanded;
  fmage = -0.050d + 0.168d*pexpanded;
  ;loctime_rad = loctime*lt_to_rad;
  r = rj;
  y_point = r*sin(loctime_rad);
  if (r gt 200.0) then r = 200d
  xplot = (-1d)*r*cos(loctime_rad)/120d;d
  bplot = dmage + fmage*xplot;
  aplot = emage;
  cplot = amage + bmage*xplot + cmage*(xplot*xplot);
  sqrt_bsq_minus4ac = sqrt(bplot*bplot - 4d*aplot*cplot)
  yplotplus =  (-1d*bplot + sqrt_bsq_minus4ac)/(2d*aplot);
  yplotminus = (-1d*bplot - sqrt_bsq_minus4ac)/(2d*aplot);
  yplotplus = -120d*yplotplus;
  yplotminus = -120d*yplotminus;
  if (y_point le yplotplus and y_point ge yplotminus) then return,1 else return,0
end







FUNCTION mag_to_ion_mapping,north_or_south,xpt,ypt,r1,ltime1,lat1,long1,r2,ltime2,lat2,long2,r3,ltime3,lat3,long3,r4,ltime4,lat4,long4
  ;x1,x2,x3,x4,y1,y2,y3,y4,xpt,ypt,deg_to_rad

  lt_to_rad = !dpi/12d
  ;Rotate ionospheric contour lat/lon positions into cartesian coordinates for easy interpolation
  ltime1_rad = ltime1*lt_to_rad
  ltime2_rad = ltime2*lt_to_rad
  ltime3_rad = ltime3*lt_to_rad
  ltime4_rad = ltime4*lt_to_rad
  x1 = r1*cos(ltime1_rad);
  y1 = r1*sin(ltime1_rad);
  x2 = r2*cos(ltime2_rad);
  y2 = r2*sin(ltime2_rad);
  x3 = r3*cos(ltime3_rad);
  y3 = r3*sin(ltime3_rad);
  x4 = r4*cos(ltime4_rad);
  y4 = r4*sin(ltime4_rad);

  lat1_rad = lat1*!dtor
  lat2_rad = lat2*!dtor
  lat3_rad = lat3*!dtor
  lat4_rad = lat4*!dtor
  long1_rad = long1*!dtor
  long2_rad = long2*!dtor
  long3_rad = long3*!dtor
  long4_rad = long4*!dtor
  x1ion = sin(lat1_rad)*cos(long1_rad);
  y1ion = sin(lat1_rad)*sin(long1_rad);
  x2ion = sin(lat2_rad)*cos(long2_rad);
  y2ion = sin(lat2_rad)*sin(long2_rad);
  x3ion = sin(lat3_rad)*cos(long3_rad);
  y3ion = sin(lat3_rad)*sin(long3_rad);
  x4ion = sin(lat4_rad)*cos(long4_rad);
  y4ion = sin(lat4_rad)*sin(long4_rad);

  xmin = x1*0.d
  xmax = x1*0.d
  ymin = x1*0.d
  ymax = x1*0.d

  for j=0, n_elements(x1) - 1 do begin
    xmin(j) = min([x1(j), x2(j), x3(j), x4(j)]);
    xmax(j) = max([x1(j), x2(j), x3(j), x4(j)]);
    ymin(j) = min([y1(j), y2(j), y3(j), y4(j)]);
    ymax(j) = max([y1(j), y2(j), y3(j), y4(j)]);
  endfor

  filematch = where(xmin le xpt and xmax ge xpt and ymin le ypt and ymax ge ypt)
  filematchpts = n_elements(filematch)
  if filematchpts gt 1 then begin
    insidept = 0
    for k=0, filematchpts-1 do begin
      insidept1 = point_inside_polygon(xpt,ypt,[x1(filematch(k)),x3(filematch(k)),x4(filematch(k)),x2(filematch(k))],[y1(filematch(k)),y3(filematch(k)),y4(filematch(k)),y2(filematch(k))])
      insidept2 = point_inside_polygon(xpt,ypt,[x1(filematch(k)),x4(filematch(k)),x3(filematch(k)),x2(filematch(k))],[y1(filematch(k)),y4(filematch(k)),y3(filematch(k)),y2(filematch(k))])
      insidept = max([insidept1,insidept2])
      if insidept eq 1 then begin
        filematch = filematch(k)
        k = k+10000
      endif
    endfor
    if insidept eq 0 and filematchpts gt 1 then begin
      for k=0, filematchpts-1 do begin
        insidept = point_inside_polygon(xpt+1.d-8,ypt+1.d-8,[x1(filematch(k)),x3(filematch(k)),x4(filematch(k)),x2(filematch(k))],[y1(filematch(k)),y3(filematch(k)),y4(filematch(k)),y2(filematch(k))])
        if insidept eq 1 then begin
          filematch = filematch(k)
          k = k+10000
        endif
      endfor
    endif
    if insidept eq 0 then filematch = -1
  endif

  x1ion = x1ion(filematch);
  x2ion = x2ion(filematch);
  x3ion = x3ion(filematch);
  x4ion = x4ion(filematch);
  y1ion = y1ion(filematch);
  y2ion = y2ion(filematch);
  y3ion = y3ion(filematch);
  y4ion = y4ion(filematch);

  x1prime = x1(filematch) - x1(filematch);
  x2prime = x2(filematch) - x1(filematch);
  x3prime = x3(filematch) - x1(filematch);
  x4prime = x4(filematch) - x1(filematch);
  xpt_prime = xpt - x1(filematch);
  y1prime = y1(filematch) - y1(filematch);
  y2prime = y2(filematch) - y1(filematch);
  y3prime = y3(filematch) - y1(filematch);
  y4prime = y4(filematch) - y1(filematch);
  ypt_prime = ypt - y1(filematch);

  theta = atan2(y3prime,x3prime)

  cos_theta = cos(theta)
  sin_theta = sin(theta)
  a1 = x1prime*cos_theta + y1prime*sin_theta;
  a2 = x2prime*cos_theta + y2prime*sin_theta;
  a3 = x3prime*cos_theta + y3prime*sin_theta;
  a4 = x4prime*cos_theta + y4prime*sin_theta;
  apt = xpt_prime*cos_theta + ypt_prime*sin_theta;
  b1 = -1.0d*x1prime*sin_theta + y1prime*cos_theta;
  b2 = -1.0d*x2prime*sin_theta + y2prime*cos_theta;
  b3 = -1.0d*x3prime*sin_theta + y3prime*cos_theta;
  b4 = -1.0d*x4prime*sin_theta + y4prime*cos_theta;
  bpt = -1.0d*xpt_prime*sin_theta + ypt_prime*cos_theta;

  x1prime = x1ion - x1ion;
  x2prime = x2ion - x1ion;
  x3prime = x3ion - x1ion;
  x4prime = x4ion - x1ion;
  y1prime = y1ion - y1ion;
  y2prime = y2ion - y1ion;
  y3prime = y3ion - y1ion;
  y4prime = y4ion - y1ion;

  thetaion = atan2(y3prime,x3prime);

  cos_thetaion = cos(thetaion)
  sin_thetaion = sin(thetaion)
  a1ion = x1prime*cos_thetaion + y1prime*sin_thetaion;
  a2ion = x2prime*cos_thetaion + y2prime*sin_thetaion;
  a3ion = x3prime*cos_thetaion + y3prime*sin_thetaion;
  a4ion = x4prime*cos_thetaion + y4prime*sin_thetaion;
  b1ion = -1.0d*x1prime*sin_thetaion + y1prime*cos_thetaion;
  b2ion = -1.0d*x2prime*sin_thetaion + y2prime*cos_thetaion;
  b3ion = -1.0d*x3prime*sin_thetaion + y3prime*cos_thetaion;
  b4ion = -1.0d*x4prime*sin_thetaion + y4prime*cos_thetaion;

  aption = a1ion + (a3ion-a1ion)*(apt-a1)/(a3-a1);
  bption = b1ion + (b2ion-b1ion)*(bpt-b1)/(b2-b1);

  xption = aption*cos_thetaion - bption*sin_thetaion + x1ion;
  yption = aption*sin_thetaion + bption*cos_thetaion + y1ion;

  if north_or_south eq 'north' then begin
    lat_interpol = 90. - (asin(sqrt(xption*xption + yption*yption)))*180./!dpi; // will need to be fixed for southern hemisphere
  endif else begin
    lat_interpol = (asin(sqrt(xption*xption + yption*yption)))*180.d/!dpi - 90.; // will need to be changed for southern hemisphere
  endelse
  long_interpol = (atan2(yption,xption))*180.d/!dpi;

  long_interpol = (360.0 - long_interpol + 360.0) mod 360; // convert to LH

  return,[lat_interpol,long_interpol]
END










FUNCTION mapping_function_2025,mapping_type,sslong2,input_position1,input_position2,MODEL,FIELDLINE_TRACING

  directory_name = 'Documents/mapping_function_files/idl_save_files/' ;USER SHOULD REPLACE THIS STRING WITH APPROPRIATE DIRECTORY NAME.
  ion_to_mag = 0
  mag_to_ion_north = 0
  mag_to_ion_south = 0
  rj = 0.d
  loctime = 0.d
  latpt_input = 0.d
  longpt_input = 0.d
  deg_to_rad = !dpi/180d
  lt_to_rad = !dpi/12d

  badpt = 0
  if mapping_type eq 'ion_to_mag' then begin 
    ion_to_mag = 1
    latpt_input = input_position1
    longpt_input = input_position2
  endif else if mapping_type eq 'mag_to_ion_north' then begin
    mag_to_ion_north = 1 
    rj = input_position1
    loctime = input_position2
  endif else if mapping_type eq 'mag_to_ion_south' then begin
    mag_to_ion_south = 1
    rj = input_position1
    loctime = input_position2
  endif else begin
    message,'ERROR. Please check the mapping type entered. Valid inputs are the following strings: ion_to_mag, mag_to_ion_north, mag_to_ion_south'
    badpt = 1
  endelse
  
  if rj eq 15.d then rj = 15.0001d
  
  latpt = latpt_input
  longpt = longpt_input
  rj_input = rj
  loctime_input = loctime
  
  sslong = sslong2
  
  ;test if user-specified point is within the model validity range (15-150 Rj in the magnetosphere), local time given in decimal format 0-24 hours,
  if ion_to_mag eq 1 then begin 
    if (abs(latpt) gt 90d) then badpt = 1
    if (longpt gt 360d) then longpt = longpt mod 360d
    if longpt lt 0d then longpt = 360d + (longpt mod 360d)
  endif else begin
    if ion_to_mag ne 1 and mag_to_ion_north ne 1 and mag_to_ion_south ne 1 then badpt = 1
  endelse
  if (mag_to_ion_north eq 1 or mag_to_ion_south eq 1) and (rj gt 150d or rj lt 15d or loctime lt 0d or loctime gt 24d) then badpt = 1
  if rj lt 15d and rj ge 2d and keyword_set(fieldline_tracing) then badpt = 0
  
  model_filename_north = '_jrm33'
  model_filename_south = '_jrm33'
  if keyword_set(model) and keyword_set(fieldline_tracing) eq 0 then begin
    if strcmp(model,'VIP4',/fold_case) then begin
      model_filename_north = '_vip4'
      model_filename_south = '_vip4'
    endif else if strcmp(model,'gam',/fold_case) then begin
      model_filename_north = '_gam'
      if mag_to_ion_south eq 1 then begin
        message,'Error. Grodent anomaly model can only be used for mapping to northern hemisphere (not south).'
        badpt = 1
      endif
    endif else if strcmp(model,'VIPAL',/fold_case) then begin
      model_filename_north = '_vipal'
      model_filename_south = '_vipal'
    endif else if strcmp(model,'JRM09',/fold_case) then begin
      model_filename_north = '_jrm09'
      model_filename_south = '_jrm09'
    endif else if strcmp(model,'JRM33',/fold_case) then begin
      model_filename_north = '_jrm33'
      model_filename_south = '_jrm33'
    endif else begin
      message,'ERROR. Please check model name. Valid inputs are: gam, vipal, vip4, jrm09, and jrm33. The khurana and khurana_jrm09 model options are only available for the fieldline tracing option.'
      badpt = 1
    endelse
 endif
  
  rmin_limit = 15d
  if keyword_set(fieldline_tracing) then begin
    use_fieldline_tracing = 1
    rmin_limit = 2d
    if model eq 'khurana' or model eq 'KHURANA' or model eq 'Khurana' then begin
      model_filename_north = 'kk2009'
      model_filename_south = 'kk2009'
    endif else if model eq 'khurana_jrm09' or model eq 'KHURANA_JRM09' then begin
      model_filename_north = 'kk2009ext_jrm09int'
      model_filename_south = 'kk2009ext_jrm09int'
    endif else if model eq 'vip4' or model eq 'VIP4' then begin
      model_filename_north = 'vip4'
      model_filename_south = 'vip4'
    endif else if model eq 'vipal' or model eq 'VIPAL' then begin
      model_filename_north = 'vipal'
      model_filename_south = 'vipal'
    endif else if model eq 'jrm09' or model eq 'JRM09' then begin
      model_filename_north = 'jrm09'
      model_filename_south = 'jrm09'
    endif else if model eq 'jrm33' or model eq 'JRM33' then begin
      model_filename_north = 'jrm33'
      model_filename_south = 'jrm33'
    endif else if model eq 'gam' or model eq 'GAM' then begin
      model_filename_north = 'gam'
      if mag_to_ion_south eq 1 then begin
        message,'Error. Grodent anomaly model can only be used for mapping to northern hemisphere (not south).'
        badpt = 1
      endif
    endif else begin
        message,'ERROR. Please check model name. Valid inputs are: gam, vipal, vip4, jrm09, jrm33, khurana, and khurana_jrm09. The khurana and khurana_jrm09 model options are only available for the fieldline tracing option'
        badpt = 1
    endelse
  endif else use_fieldline_tracing = 0
  
  if ion_to_mag eq 0 then begin ; test if magnetosphere point given by user is beyond joy et al magnetopause, if yes then return flag value
    loctime_rad = loctime*lt_to_rad;
    is_inside = is_inside_joy_magnetopause(rj,loctime_rad)
    if is_inside eq 0 then badpt = 1
    rj = double(rj)
    loctime = double(loctime)
  endif else begin
    latpt = double(latpt)
    longpt = double(longpt)
  endelse
  
  sslong_old = sslong + 360.d mod 360.d
  sslong = 360.d - sslong + 360.d
  sslong = sslong mod 360.d
    
  nelements = -1;
  filematch = -1;
  
  r_interpol = -996.d
  ltime_interpol = -996.d
  lat_interpol = -996.d
  long_interpol = -996.d
  
  ;field line tracing for VIP4, GAM, VIPAL not valid outside 100 Rj
  if keyword_set(fieldline_tracing) and model_filename_north ne '_kk2009' then begin
    if rj gt 85.d then badpt = 1
  endif
  
  if badpt eq 0 then begin
    if ion_to_mag eq 1 then begin ; if user wants to map from the ionosphere to the magnetosphere
      longpt = 360.0d - longpt + 360.0d; // covert to right handed
      longpt = longpt mod 360.0d

      latpt_start = latpt;

      xpt = sin((90.d -latpt)*deg_to_rad)*cos(longpt*deg_to_rad);
      ypt = sin((90.d -latpt)*deg_to_rad)*sin(longpt*deg_to_rad);

      ;loop northern hemisphere
      if (latpt ge 0.d) then begin
        latpt = 90.d - latpt; // convert to colatitude

        if sslong mod 10.d eq 0.d then begin ; mapping from northern ionosphere to magnetosphere
          sslong1 = sslong;
          sslong = sslong_old;
          sslongfix = fix(sslong)
          
          if sslong ne 0. and use_fieldline_tracing eq 0 then filetext = directory_name + 'sslong'+strtrim(sslongfix,2)+model_filename_north+'.sav'
          if sslong eq 0. and use_fieldline_tracing eq 0 then filetext = directory_name + 'sslong360'+model_filename_north+'.sav'
          if sslong ne 0. and use_fieldline_tracing eq 1 then filetext = directory_name + 'fieldline_tracing_all_r_'+model_filename_north+'_north_cml'+strtrim(sslongfix,2)+'.sav' 
          if sslong eq 0. and use_fieldline_tracing eq 1 then filetext = directory_name + 'fieldline_tracing_all_r_'+model_filename_north+'_north_cml360.sav'                   
                        
          sslong = sslong1;
          restore,filetext
          
          x1 = sin(lat1*!dtor)*cos(long1*!dtor);
          y1 = sin(lat1*!dtor)*sin(long1*!dtor);
          x2 = sin(lat2*!dtor)*cos(long2*!dtor);
          y2 = sin(lat2*!dtor)*sin(long2*!dtor);
          x3 = sin(lat3*!dtor)*cos(long3*!dtor);
          y3 = sin(lat3*!dtor)*sin(long3*!dtor);
          x4 = sin(lat4*!dtor)*cos(long4*!dtor);
          y4 = sin(lat4*!dtor)*sin(long4*!dtor);
          filematch = identify_ion_to_mag_match(xpt,ypt,x1,x2,x3,x4,y1,y2,y3,y4)
                             
          rstepsize = 5.d
          
          ;===========
          ;For fieldline tracing only - test if point is inside of 15 Rj but outside of 2 Rj
          ;===========
          if filematch(0) eq -1 and use_fieldline_tracing eq 1 then begin
            if sslong ne 0. then filetext = directory_name + 'fieldline_tracing_2to15rj_'+model_filename_north+'_north_cml'+strtrim(sslongfix,2)+'.sav'
            if sslong eq 0. then filetext = directory_name + 'fieldline_tracing_2to15rj_'+model_filename_north+'_north_cml360.sav'
            sslong = sslong1;
            restore,filetext
            x1 = sin(lat1*!dtor)*cos(long1*!dtor);
            y1 = sin(lat1*!dtor)*sin(long1*!dtor);
            x2 = sin(lat2*!dtor)*cos(long2*!dtor);
            y2 = sin(lat2*!dtor)*sin(long2*!dtor);
            x3 = sin(lat3*!dtor)*cos(long3*!dtor);
            y3 = sin(lat3*!dtor)*sin(long3*!dtor);
            x4 = sin(lat4*!dtor)*cos(long4*!dtor);
            y4 = sin(lat4*!dtor)*sin(long4*!dtor);
            filematch = identify_ion_to_mag_match(xpt,ypt,x1,x2,x3,x4,y1,y2,y3,y4)
            rstepsize = 0.5
          endif
          ;=====
          ;End only inside of 15 Rj part     
          ;=====


          ;==========
          ;Do interpolation to calculate mapping of ionospheric point within box
          ;==========
          if filematch(0) ne -1 then begin
            interpolation_output = mapping_interpolation_ion_to_mag(xpt,ypt,filematch,rstepsize,r1,ltime1,lat1,long1,r2,ltime2,lat2,long2,r3,ltime3,lat3,long3,r4,ltime4,lat4,long4,x1,x2,x3,x4,y1,y2,y3,y4)
            r_interpol = interpolation_output(0)
            ltime_interpol = interpolation_output(1)
            ;If filematch works fine, do final test to see if final mapping point is inside the magnetosphere or not ; from joy et al 2002
            is_inside = is_inside_joy_magnetopause(r_interpol,ltime_interpol*lt_to_rad)
            if (is_inside eq 0) then begin
              ; beyond 150 RJ or expanded magnetosphere
              r_interpol = -999.d
              ltime_interpol = -999.d
            endif
          endif else begin ;if unable to find match for point in the files         
            x1_array = sin(lat1*!dtor)*cos(long1*!dtor)
            y1_array = sin(lat1*!dtor)*sin(long1*!dtor)
            ;test to see if point is inside of the 15 Rj contour
            jloop = n_elements(where(r1 eq rmin_limit))
            jloop = jloop - 1;           
            c = 0;
            jloop_start = jloop;
            jloop = jloop;
            nvert = jloop+1;
            y1_match = y1_array(where(r1 eq rmin_limit))
            x1_match = x1_array(where(r1 eq rmin_limit))
            c = point_inside_polygon(xpt,ypt,x1_match,y1_match)
            inside15 = c;
            if (inside15 eq 1) then begin ;maps to outside 150 RJ or Jovian magnetopause
              r_interpol = -999.d
              ltime_interpol = -999.d
            endif else begin; maps to inside of 15 RJ (or 2 RJ for fieldline tracing)
              r_interpol = -998.d
              ltime_interpol = -998.d
            endelse
          endelse;end bad file match
                    
          ;end of loop sslong div by 10
        endif else begin ;if sslong div by 10 != 0, still mapping from northern ionosphere to magnetosphere
;          print,'hi'
          sslong1 = sslong;
          sslong = sslong_old;
          
          sslongmatch_min = dindgen(36)*10.d
          sslongmatch_max = (dindgen(36)+1.d)*10.d
          sslongmatch_file = [360,(indgen(35)+1)*10]
          
          sslongmatch1 = where(sslong ge sslongmatch_min and sslong lt sslongmatch_max)
          
          filetext = directory_name + 'sslong'+strtrim(sslongmatch_file(sslongmatch1(0)),2)+model_filename_north+'.sav'
          sslong = sslong1;
          
          if keyword_set(fieldline_tracing) then begin 
            first_mapping = mapping_function_2025('ion_to_mag',sslongmatch_file(sslongmatch1(0)),latpt_input,longpt_input,model,'fieldline_tracing')
            second_mapping = mapping_function_2025('ion_to_mag',(sslongmatch_file(sslongmatch1(0))+10) mod 360,latpt_input,longpt_input,model,'fieldline_tracing')
          endif else begin
            first_mapping = mapping_function_2025('ion_to_mag',sslongmatch_file(sslongmatch1(0)),latpt_input,longpt_input,model)
            second_mapping = mapping_function_2025('ion_to_mag',(sslongmatch_file(sslongmatch1(0))+10) mod 360,latpt_input,longpt_input,model)
          endelse
                    
          
          if abs(first_mapping(0)) ge 200. then begin
            r_interpol = second_mapping(0)
            ltime_interpol = second_mapping(1)
          endif else if abs(second_mapping(0)) ge 200. then begin
            r_interpol = first_mapping(0)
            ltime_interpol = first_mapping(1)
          endif else begin
            sslong_first = float(sslongmatch_file(sslongmatch1(0)))
            if sslong_first eq 0 or sslong_first eq 360 then sslong_first = 0.
            first_x = first_mapping(0)*cos(first_mapping(1)*lt_to_rad)
            first_y = first_mapping(0)*sin(first_mapping(1)*lt_to_rad)
            second_x = second_mapping(0)*cos(second_mapping(1)*lt_to_rad)
            second_y = second_mapping(0)*sin(second_mapping(1)*lt_to_rad)
            x_interpol = first_x + (second_x-first_x)*((360d -sslong)-sslong_first)/10.d
            y_interpol = first_y + (second_y-first_y)*((360d -sslong)-sslong_first)/10.d
            r_interpol = sqrt(x_interpol*x_interpol + y_interpol*y_interpol)
            ltime_interpol = (atan2(y_interpol,x_interpol)*12.d/!dpi + 48.) mod 24.
          endelse
          
        endelse ;end sslong div 10 != 0 loop for northern ionosphere to magnetosphere
        
        
        
        
        
      endif else begin ;// begin ionosphere->magnetosphere southern loop
      
     
      
        latpt = 90.d + latpt; // convert to colatitude
        
        if sslong mod 10.d eq 0.d then begin ; mapping southern ionosphere to magnetosphere
          sslong1 = sslong;
          sslong = sslong_old;
          sslongfix = fix(sslong)
          
          if sslong ne 0. and use_fieldline_tracing eq 0 then filetext = directory_name + 'sslong'+strtrim(sslongfix,2)+'s'+model_filename_south+'.sav'
          if sslong eq 0. and use_fieldline_tracing eq 0 then filetext = directory_name + 'sslong360s'+model_filename_south+'.sav'
          if sslong ne 0. and use_fieldline_tracing eq 1 then filetext = directory_name + 'fieldline_tracing_all_r_'+model_filename_south+'_south_cml'+strtrim(sslongfix,2)+'.sav'
          if sslong eq 0. and use_fieldline_tracing eq 1 then filetext = directory_name + 'fieldline_tracing_all_r_'+model_filename_south+'_south_cml360.sav'
          
          sslong = sslong1
          
          restore,filetext
                
          x1 = sin(lat1*!dtor)*cos(long1*!dtor);
          y1 = sin(lat1*!dtor)*sin(long1*!dtor);
          x2 = sin(lat2*!dtor)*cos(long2*!dtor);
          y2 = sin(lat2*!dtor)*sin(long2*!dtor);
          x3 = sin(lat3*!dtor)*cos(long3*!dtor);
          y3 = sin(lat3*!dtor)*sin(long3*!dtor);
          x4 = sin(lat4*!dtor)*cos(long4*!dtor);
          y4 = sin(lat4*!dtor)*sin(long4*!dtor);
          filematch = identify_ion_to_mag_match(xpt,ypt,x1,x2,x3,x4,y1,y2,y3,y4)

          rstepsize = 5.d

          ;===========
          ;For fieldline tracing only - test if point is inside of 15 Rj but outside of 2 Rj
          ;===========
          if filematch(0) eq -1 and use_fieldline_tracing eq 1 then begin
            if sslong ne 0. then filetext = directory_name + 'fieldline_tracing_2to15rj_'+model_filename_north+'_south_cml'+strtrim(sslongfix,2)+'.sav'
            if sslong eq 0. then filetext = directory_name + 'fieldline_tracing_2to15rj_'+model_filename_north+'_south_cml360.sav'
            sslong = sslong1;
            restore,filetext
            x1 = sin(lat1*!dtor)*cos(long1*!dtor);
            y1 = sin(lat1*!dtor)*sin(long1*!dtor);
            x2 = sin(lat2*!dtor)*cos(long2*!dtor);
            y2 = sin(lat2*!dtor)*sin(long2*!dtor);
            x3 = sin(lat3*!dtor)*cos(long3*!dtor);
            y3 = sin(lat3*!dtor)*sin(long3*!dtor);
            x4 = sin(lat4*!dtor)*cos(long4*!dtor);
            y4 = sin(lat4*!dtor)*sin(long4*!dtor);
            filematch = identify_ion_to_mag_match(xpt,ypt,x1,x2,x3,x4,y1,y2,y3,y4)
            rstepsize = 0.5
          endif
          ;=====
          ;End only inside of 15 Rj part
          ;=====

          ;==========
          ;Do interpolation to calculate mapping of ionospheric point within box
          ;==========
          if filematch(0) ne -1 then begin
            interpolation_output = mapping_interpolation_ion_to_mag(xpt,ypt,filematch,rstepsize,r1,ltime1,lat1,long1,r2,ltime2,lat2,long2,r3,ltime3,lat3,long3,r4,ltime4,lat4,long4,x1,x2,x3,x4,y1,y2,y3,y4)
            r_interpol = interpolation_output(0)
            ltime_interpol = interpolation_output(1)
            is_inside = is_inside_joy_magnetopause(r_interpol,ltime_interpol*lt_to_rad) ;If filematch works fine, do final test to see if final mapping point is inside the magnetosphere or not ; from joy et al 2002
            if (is_inside eq 0) then begin
              ; beyond 150 RJ or expanded magnetosphere
              r_interpol = -999.d
              ltime_interpol = -999.d
            endif
          endif else begin          ;if unable to find match for point in the files
            x1_array = sin(lat1*!dtor)*cos(long1*!dtor)
            y1_array = sin(lat1*!dtor)*sin(long1*!dtor)
            ;test to see if point is inside of the 15 Rj contour
            jloop = n_elements(where(r1 eq rmin_limit)) - 1
            c = 0;
            jloop_start = jloop;
            nvert = jloop+1;
            y1_match = y1_array(where(r1 eq rmin_limit))
            x1_match = x1_array(where(r1 eq rmin_limit))
            c = point_inside_polygon(xpt,ypt,x1_match,y1_match)
            inside15 = c;
            if (inside15 eq 1) then begin ;maps to outside 150 RJ or Jovian magnetopause
              r_interpol = -999.d
              ltime_interpol = -999.d
            endif else begin; maps to inside of 15 RJ (or 2 RJ for fieldline tracing)
              r_interpol = -998.d
              ltime_interpol = -998.d
            endelse
          endelse; end bad file match

          ;end of loop sslong div by 10
        endif else begin ;if sslong div by 10 != 0 still southern ionosphere to magnetosphere
          sslong1 = sslong;
          sslong = sslong_old;
          
          sslongmatch_min = dindgen(36)*10.d
          sslongmatch_max = (dindgen(36)+1.d)*10.d
          sslongmatch_file = [360,(indgen(35)+1)*10]
          
          sslongmatch1 = where(sslong ge sslongmatch_min and sslong lt sslongmatch_max)
          
          filetext = directory_name + 'sslong'+strtrim(sslongmatch_file(sslongmatch1(0)),2)+'s'+model_filename_south+'.sav'
          sslong = sslong1;
          
          if keyword_set(fieldline_tracing) then begin
            first_mapping = mapping_function_2025('ion_to_mag',sslongmatch_file(sslongmatch1(0)),latpt_input,longpt_input,model,'fieldline_tracing')
            second_mapping = mapping_function_2025('ion_to_mag',(sslongmatch_file(sslongmatch1(0))+10) mod 360,latpt_input,longpt_input,model,'fieldline_tracing')
          endif else begin
            first_mapping = mapping_function_2025('ion_to_mag',sslongmatch_file(sslongmatch1(0)),latpt_input,longpt_input,model)
            second_mapping = mapping_function_2025('ion_to_mag',(sslongmatch_file(sslongmatch1(0))+10) mod 360,latpt_input,longpt_input,model)
          endelse
          
          sslong_first = float(sslongmatch_file(sslongmatch1(0)))
          if sslong_first eq 0 or sslong_first eq 360 then sslong_first = 0.
          
          if abs(first_mapping(0)) ge 200. then begin
            r_interpol = second_mapping(0)
            ltime_interpol = seconmd_mapping(1)
          endif else if abs(second_mapping(0)) ge 200. then begin
            r_interpol = first_mapping(0)
            ltime_interpol = first_mapping(1)
          endif else begin
            first_x = first_mapping(0)*cos(first_mapping(1)*lt_to_rad)
            first_y = first_mapping(0)*sin(first_mapping(1)*lt_to_rad)
            second_x = second_mapping(0)*cos(second_mapping(1)*lt_to_rad)
            second_y = second_mapping(0)*sin(second_mapping(1)*lt_to_rad)
            x_interpol = first_x + (second_x-first_x)*((360d -sslong)-sslong_first)/10.d
            y_interpol = first_y + (second_y-first_y)*((360d -sslong)-sslong_first)/10.d
            r_interpol = sqrt(x_interpol^2d + y_interpol^2d)
            ltime_interpol = (atan2(y_interpol,x_interpol)*(12.d/!dpi) + 48.) mod 24.
          endelse
          
        endelse ;end sslong div 10 != 0 loop
        
        
        
      endelse ;end southern hemisphere ion to mag mapping
      
      
      

      
    endif else if mag_to_ion_north eq 1 then begin
      xpt = rj*cos(loctime_rad);
      ypt = rj*sin(loctime_rad);
      ;latpt = 90. - latpt; // convert to colatitude
      
      if sslong mod 10.d eq 0.d then begin ; magnetosphere to ionosphere north
        sslong1 = sslong;
        sslong = sslong_old;
        sslongfix = fix(sslong)
        
        ;REPLACE DIRECTORY NAME HERE
        if sslong ne 0. then filetext = directory_name + 'sslong'+strtrim(sslongfix,2)+model_filename_north+'.sav'
        if sslong eq 0. then filetext = directory_name + 'sslong360'+model_filename_north+'.sav'
        if sslong ne 0. and use_fieldline_tracing eq 1 then filetext = directory_name + 'fieldline_tracing_all_r_'+model_filename_north+'_north_cml'+strtrim(sslongfix,2)+'.sav'
        if sslong eq 0. and use_fieldline_tracing eq 1 then filetext = directory_name + 'fieldline_tracing_all_r_'+model_filename_north+'_north_cml360.sav'
        if sslong ne 0. and use_fieldline_tracing eq 1 and rj lt 15. then filetext = directory_name + 'fieldline_tracing_2to15rj_'+model_filename_north+'_north_cml'+strtrim(sslongfix,2)+'.sav'
        if sslong eq 0. and use_fieldline_tracing eq 1 and rj lt 15. then filetext = directory_name + 'fieldline_tracing_2to15rj_'+model_filename_north+'_north_cml360.sav'
        
        sslong = sslong1;
        
        restore,filetext
        
        mag_to_ion_mapping_output = mag_to_ion_mapping('north',xpt,ypt,r1,ltime1,lat1,long1,r2,ltime2,lat2,long2,r3,ltime3,lat3,long3,r4,ltime4,lat4,long4)
        lat_interpol = mag_to_ion_mapping_output(0)
        long_interpol = mag_to_ion_mapping_output(1)

      endif else begin; if sslong % 10 ne 0 still magnetosphere to ionosphere north
        sslong1 = sslong;
        sslong = sslong_old;
        
        sslongmatch_min = dindgen(36)*10.d
        sslongmatch_max = (dindgen(36)+1.d)*10.d
        sslongmatch_file = [360,(indgen(35)+1)*10]
        
        sslongmatch1 = where(sslong ge sslongmatch_min and sslong lt sslongmatch_max)
        
        ;REPLACE DIRECTORY NAME HERE
        filetext = directory_name + 'sslong'+strtrim(sslongmatch_file(sslongmatch1(0)),2)+model_filename_north+'.sav'
        sslong = sslong1;
        
        if keyword_set(fieldline_tracing) then begin
          first_mapping = mapping_function_2025('mag_to_ion_north',sslongmatch_file(sslongmatch1(0)),rj_input,loctime_input,model,'fieldline_tracing')
          second_mapping = mapping_function_2025('mag_to_ion_north',(sslongmatch_file(sslongmatch1(0))+10) mod 360,rj_input,loctime_input,model,'fieldline_tracing')
        endif else begin
          first_mapping = mapping_function_2025('mag_to_ion_north',sslongmatch_file(sslongmatch1(0)),rj_input,loctime_input,model)
          second_mapping = mapping_function_2025('mag_to_ion_north',(sslongmatch_file(sslongmatch1(0))+10) mod 360,rj_input,loctime_input,model)
        endelse
        
        sslong_first = float(sslongmatch_file(sslongmatch1(0)))
        if sslong_first eq 0 or sslong_first eq 360 then sslong_first = 0.
        

        if abs(first_mapping(0)) ge 200. then begin
          lat_interpol = second_mapping(0)
          long_interpol = second_mapping(1)
        endif else if abs(second_mapping(0)) ge 200. then begin
          lat_interpol = first_mapping(0)
          long_interpol = first_mapping(1)
        endif else begin
          first_x = abs(sin((90.-first_mapping(0))*deg_to_rad))*cos(first_mapping(1)*deg_to_rad)
          first_y = abs(sin((90.-first_mapping(0))*deg_to_rad))*sin(first_mapping(1)*deg_to_rad)
          second_x = abs(sin((90.-second_mapping(0))*deg_to_rad))*cos(second_mapping(1)*deg_to_rad)
          second_y = abs(sin((90.-second_mapping(0))*deg_to_rad))*sin(second_mapping(1)*deg_to_rad)
          x_interpol = first_x + (second_x-first_x)*((360.-sslong)-sslong_first)/10.d
          y_interpol = first_y + (second_y-first_y)*((360.-sslong)-sslong_first)/10.d
          lat_interpol = 90. - asin(sqrt(x_interpol^2.d + y_interpol^2.d))*180.d/!dpi
          long_interpol = (atan2(y_interpol,x_interpol)*180.d/!dpi + 720.) mod 360.          
        endelse

       
      endelse ; end ionosphere to magnetosphere north sslong != 10

      
      
      
    endif else if mag_to_ion_south eq 1 then begin ;end mag to ion north
    
      xpt = rj*cos(loctime_rad);
      ypt = rj*sin(loctime_rad);
      latpt = 90. - latpt; // convert to colatitude
      ;lat_interpol = 90.0 - lat_interpol;
      
      if sslong mod 10.d eq 0.d then begin
        sslong1 = sslong;
        sslong = sslong_old;
        sslongfix = fix(sslong)
        
        ;REPLACE DIRECTORY NAME HERE
        if sslong ne 0. then filetext = directory_name + 'sslong'+strtrim(sslongfix,2)+'s'+model_filename_south+'.sav'
        if sslong eq 0. then filetext = directory_name + 'sslong360s'+model_filename_south+'.sav'
        if sslong ne 0. and use_fieldline_tracing eq 1 then filetext = directory_name + 'fieldline_tracing_all_r_'+model_filename_south+'_south_cml'+strtrim(sslongfix,2)+'.sav'
        if sslong eq 0. and use_fieldline_tracing eq 1 then filetext = directory_name + 'fieldline_tracing_all_r_'+model_filename_south+'_south_cml360.sav'
        if sslong ne 0. and use_fieldline_tracing eq 1 and rj lt 15. then filetext = directory_name + 'fieldline_tracing_2to15rj_'+model_filename_south+'_south_cml'+strtrim(sslongfix,2)+'.sav'
        if sslong eq 0. and use_fieldline_tracing eq 1 and rj lt 15. then filetext = directory_name + 'fieldline_tracing_2to15rj_'+model_filename_south+'_south_cml360.sav'

        sslong = sslong1;
        
        restore,filetext
        
        mag_to_ion_mapping_output = mag_to_ion_mapping('south',xpt,ypt,r1,ltime1,lat1,long1,r2,ltime2,lat2,long2,r3,ltime3,lat3,long3,r4,ltime4,lat4,long4)
        lat_interpol = mag_to_ion_mapping_output(0)
        long_interpol = mag_to_ion_mapping_output(1)

      endif else begin; if sslong % 10 ne 0, mag to ion south
        sslong1 = sslong;
        sslong = sslong_old;
        
        sslongmatch_min = dindgen(36)*10.d
        sslongmatch_max = (dindgen(36)+1.d)*10.d
        sslongmatch_file = [360,(indgen(35)+1)*10]
        
        sslongmatch1 = where(sslong ge sslongmatch_min and sslong lt sslongmatch_max)
        
        filetext = directory_name + 'sslong'+strtrim(sslongmatch_file(sslongmatch1(0)),2)+'s'+model_filename_south+'.sav'
        sslong = sslong1;
        
        if keyword_set(fieldline_tracing) then begin
          first_mapping = mapping_function_2025('mag_to_ion_south',sslongmatch_file(sslongmatch1(0)),rj_input,loctime_input,model,'fieldline_tracing')
          second_mapping = mapping_function_2025('mag_to_ion_south',(sslongmatch_file(sslongmatch1(0))+10) mod 360,rj_input,loctime_input,model,'fieldline_tracing')
        endif else begin          
          first_mapping = mapping_function_2025('mag_to_ion_south',sslongmatch_file(sslongmatch1(0)),rj_input,loctime_input,model)
          second_mapping = mapping_function_2025('mag_to_ion_south',(sslongmatch_file(sslongmatch1(0))+10) mod 360,rj_input,loctime_input,model)
        endelse
        
        sslong_first = float(sslongmatch_file(sslongmatch1(0)))
        if sslong_first eq 0 or sslong_first eq 360 then sslong_first = 0.
                
        if abs(first_mapping(0)) ge 200. then begin
          lat_interpol = second_mapping(0)
          long_interpol = second_mapping(1)
        endif else if abs(second_mapping(0)) ge 200. then begin
          lat_interpol = first_mapping(0)
          long_interpol = first_mapping(1)
        endif else begin
          first_x = sin((90.d - abs(first_mapping(0)))*deg_to_rad)*cos(first_mapping(1)*deg_to_rad)
          first_y = sin((90.d - abs(first_mapping(0)))*deg_to_rad)*sin(first_mapping(1)*deg_to_rad)
          second_x = sin((90.d - abs(second_mapping(0)))*deg_to_rad)*cos(second_mapping(1)*deg_to_rad)
          second_y = sin((90.d - abs(second_mapping(0)))*deg_to_rad)*sin(second_mapping(1)*deg_to_rad)
          x_interpol = first_x + (second_x-first_x)*((360.-sslong)-sslong_first)/10.d
          y_interpol = first_y + (second_y-first_y)*((360.-sslong)-sslong_first)/10.d
          lat_interpol = asin(sqrt(x_interpol^2.d + y_interpol^2.d))*180.d/!dpi - 90.d
          long_interpol = (atan2(y_interpol,x_interpol)*180.d/!dpi + 720.) mod 360.
        endelse
        
        
      endelse ; end mod 10 ne 0
      
      
    endif ;end mag to ion south mapping
    
    
    
  endif
  sslong = sslong2

  
  if (ion_to_mag eq 1 and abs(r_interpol) lt 15. and use_fieldline_tracing eq 0) or (ion_to_mag eq 1 and abs(r_interpol) lt 2. and use_fieldline_tracing eq 1) then begin
    r_interpol = 999d
    ltime_interpol = 999d
  endif
  if ion_to_mag eq 1 then begin
    return,double([r_interpol,ltime_interpol])
  endif else begin
    return,double([lat_interpol,long_interpol])  
  endelse
  
END