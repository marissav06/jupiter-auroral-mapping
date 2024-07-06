; Jovian magnetosphere/ionosphere mapping program - 2024 version 1
; Written by Marissa Vogt, mvogt@bu.edu - now at the Planetary Science Institute (mvogt@psi.edu)
; Last updated July 2024
;
; PURPOSE OF THIS FUNCTION
; This function allows a user to magnetically map a point between Jupiter's middle magnetosphere and the ionosphere. 
; Points can be mapped using the flux equivalence mapping model of Vogt et al. [2011, 2015] or by tracing field lines from a global field model. 
; 
; HOW TO USE THIS FUNCTION
; (Note that this function maps one point at a time but user can call it in a loop to map several points at once.)
; **Before running this program the user should replace the directory_name variable on line 233 below with the appropriate directory for their local machine
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
;     User MUST specify a field model. The user can specify that the program should use traced field lines from:
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
;   EXAMPLES - All using 180 degrees subsolar longitude
;   1. For mapping from the northern ionosphere to the magnetosphere using the Vogt et al. (2011, 2015) flux equivalence model, 
;   with input position 60 degrees latitude, 180 degrees longitude:
;     a) using the Vogt et al. flux equivalence calculation with the Grodent anomaly model:
;       result = mapping_function_2024('ion_to_mag',180.,60.,180.,'gam')
;       print,result
;       > result = [87.097377,10.260316]  ;;; means point maps to ~87 Rj and ~10.2 LT
;     b) using the Vogt et al. flux equivalence calculation with VIP4:
;       result = mapping_function_2024('ion_to_mag',180.,60.,180.,'vip4')
;       print,result
;       > result = [58.391824,10.059707]  ;;; means point maps to ~58 Rj and ~10.0 LT 
;     c) using the Vogt et al. flux equivalence calculation with VIPAL:
;       result = mapping_function_2024('ion_to_mag',180.,60.,180.,'vipal')
;       print,result
;       > result = [82.094336,9.8960662]  ;;; means point maps to ~82 Rj and ~9.9 LT
;     d) using the Vogt et al. flux equivalence calculation with JRM33:
;       result = mapping_function_2024('ion_to_mag',180.,60.,180.,'jrm09')
;       print,result
;       > result = [80.663835545003778,10.219790900448057]  ;;; means point maps to ~81 Rj and ~10.2 LT
;       You will get the same result if you call result = mapping_function_2024('ion_to_mag',180.,60.,180.) since JRM33 is the default.
;
;   2. For mapping from the southern ionosphere to the magnetosphere using the Vogt et al. (2011, 2015) flux equivalence model, 
;     with input position -80 degrees latitude, 90 degrees longitude:
;     a) using the Vogt et al. flux equivalence calculation with VIP4:
;       result = mapping_function_2024('ion_to_mag',180.,-80.,90.,'vip4')
;       print,result
;       > result = [81.154379,13.894088]  ;;; means point maps to ~81 Rj and ~13.8 LT
;     b) using the Vogt et al. flux equivalence calculation with VIPAL:
;       result = mapping_function_2024('ion_to_mag',180.,-80.,90.,'vipal')
;       print,result
;       > result = [84.538550,14.790648]  ;;; means point maps to ~84 Rj and ~14.8 LT
;     c) using the Vogt et al. flux equivalence calculation with JRM09:
;       result = mapping_function_2024('ion_to_mag',180.,-80.,90.,'jrm09')
;       print,result
;       > result = [82.122500,14.700538]  ;;; means point maps to ~82 Rj and ~14.7 LT
;
;   3. For mapping from the magnetosphere to northern ionosphere using the Vogt et al. (2011, 2015) flux equivalence model,
;   with input position 88 Rj and 9.8 LT:
;     a) using the Vogt et al. flux equivalence calculation with the Grodent anomaly model:
;       result = mapping_function_2024('mag_to_ion_north',180.,88.,9.8,'gam')
;       print,result
;       > result = [60.172613,181.95319]  ;;; means point maps to ~60 degrees latitude and ~182 degrees longitude
;     b) using the Vogt et al. flux equivalence calculation with VIP4:
;       result = mapping_function_2024('mag_to_ion_north',180.,88.,9.8,'vip4')
;       print,result
;       > result = [62.564124,179.15706]  ;;; means point maps to ~63 degrees latitude and ~179 degrees longitude
;     c) using the Vogt et al. flux equivalence calculation with VIPAL:
;       result = mapping_function_2024('mag_to_ion_north',180.,88.,9.8,'vipal')
;       print,result
;       > result = [60.563943,180.32903]  ;;; means point maps to ~61 degrees latitude and ~180 degrees longitude
;     d) using the Vogt et al. flux equivalence calculation with JRM09:
;       result = mapping_function_2024('mag_to_ion_north',180.,88.,9.8,'jrm33')
;       print,result
;       > result = [60.655248788443210,181.46104506560255]  ;;; means point maps to ~61 degrees latitude and ~181 degrees longitude
;       You will get the same result if you call result = mapping_function_2024('mag_to_ion_north',180.,88.,9.8) instead.
;
;   4. For mapping from the magnetosphere to southern ionosphere using the Vogt et al. (2011, 2015) flux equivalence model,
;   with input position 88 Rj and 9.8 LT:
;     a) using the Vogt et al. flux equivalence calculation with VIP4:
;       result = mapping_function_2024('mag_to_ion_south',180.,88.,9.8,'vip4')
;       print,result
;       > result = [-88.276616,138.03374]  ;;; means point maps to roughly -88 degrees latitude and ~138 degrees longitude
;     b) using the Vogt et al. flux equivalence calculation with VIPAL:
;       result = mapping_function_2024('mag_to_ion_south',180.,88.,9.8,'vipal')
;       print,result
;       > result = [-88.102482,192.69359]  ;;; means point maps to roughly -88 degrees latitude and ~193 degrees longitude
;     c) using the Vogt et al. flux equivalence calculation with JRM09:
;       result = mapping_function_2024('mag_to_ion_south',180.,88.,9.8,'jrm09')
;       print,result
;       > result = [-87.800808,189.69086]  ;;; means point maps to roughly -88 degrees latitude and ~191 degrees longitude
;     
;   5. For mapping from the ionosphere to the magnetosphere by tracing field lines from a global field model,
;   with input position 57 degrees latitude, 180 degrees longitude:
;     a) Using field lines traced from the Grodent anomaly model:
;       result = mapping_function_2024('ion_to_mag',180.,57.,180,'gam','fieldline_tracing')
;       print,result
;       > result = [36.673552,11.088121]  ;;; means point maps to ~37 Rj and ~11.1 LT
;     b) Using field lines traced from the JRM09 field model:
;       result = mapping_function_2024('ion_to_mag',180.,57.,180,'jrm09','fieldline_tracing')
;       print,result
;       > result = [29.313409,11.084563]  ;;; means point maps to ~29 Rj and ~11 LT
;       
;   6. For mapping from the magnetosphere to the northern ionosphere by tracing field lines from a global field model,
;   with input position 47 Rj and 18.2 LT:
;     a) Using field lines traced from the Grodent anomaly model:
;       result = mapping_function_2024('mag_to_ion_north',180.,47.,18.2,'gam','fieldline_tracing')
;       print,result
;       > result = [61.475270,150.65616]  ;;; means point maps to ~61 degrees latitude and ~150 degrees longitude
;     b) Using field lines traced from the VIP4 field model:
;       result = mapping_function_2024('mag_to_ion_north',180.,47.,18.2,'vip4','fieldline_tracing')
;       print,result
;       > result = [66.327302,140.63710]  ;;; means point maps to ~66 degrees latitude and ~141 degrees longitude
;     c) Using field lines traced from the Khurana field model (with VIP4, default):
;       result = mapping_function_2024('mag_to_ion_north',180.,47.,18.2,'khurana','fieldline_tracing')
;       print,result
;       > result = [66.421892,139.31335]  ;;; means point maps to ~66 degrees latitude and ~139 degrees longitude
;     d) Using field lines traced from the Khurana field model with JRM09:
;       result = mapping_function_2024('mag_to_ion_north',180.,47.,18.2,'khurana_jrm09','fieldline_tracing')
;       print,result
;       > result = [63.807225,149.44917]  ;;; means point maps to ~64 degrees latitude and ~149 degrees longitude
;
;   7. For mapping from the magnetosphere to the southern ionosphere by tracing field lines from a global field model,
;   with input position 47 Rj and 18.2 LT:
;     a) Using field lines traced from VIPAL:
;       result = mapping_function_2024('mag_to_ion_south',180.,47.,18.2,'vipal','fieldline_tracing')
;       print,result
;       > result = [-69.919107,68.309593]  ;;; means point maps to roughly -70 degrees latitude, 68 degrees longitude
;     b) Using field lines traced from the JRM09 field model:
;       result = mapping_function_2024('mag_to_ion_south',180.,47.,18.2,'jrm09','fieldline_tracing')
;       print,result
;       > result = [-69.741398,69.864959]  ;;; means point maps to roughly -70 degrees latitude, 70 degrees longitude
;
; INPUT FILES AND PROGRAMS CALLED
;   Calls routines: atan2, point_inside_polygon
;
;   The program reads in data from files called sslong0_gam.sav, sslong10_jrm09.sav, etc., which should be placed in an appropriate directory on the user's machine.
;   ***Before running this program the user will have to replace the directory_name variable in the code (see line 233).
;
; ACKNOWLEDGMENTS
;   Thanks to Bertrand Bonfond, Benjamin Palmaerts, and Andrew Steffl for helping to test early versions of this code and the online mapping tool.
;   Very special thanks to Masafumi Imai for assistance implementing the JRM09 field model.
;
; WHEN USING RESULTS FROM THIS CODE IN A PAPER OR PRESENTATION, PLEASE CITE
;   Vogt, M. F., M. G. Kivelson, K. K. Khurana, R. J. Walker, B. Bonfond, D. Grodent, and A. Radioti (2011), 
;   Improved mapping of Jupiter’s auroral features to magnetospheric sources, 
;   J. Geophys. Res., 116, A03220, doi:10.1029/2010JA016148.
;
;   Vogt, M. F., E. J. Bunce, M. G. Kivelson, K. K. Khurana, R. J. Walker, A. Radioti, B. Bonfond, and D. Grodent (2015), 
;   Magnetosphere-ionosphere mapping at Jupiter: Quantifying the effects of using different internal field models, 
;   J. Geophys. Res. Space Physics, 10.1002/2014JA020729
;
;
; Most recent changes:
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

FUNCTION mapping_function_2024,mapping_type,sslong2,input_position1,input_position2,MODEL,FIELDLINE_TRACING

  directory_name = 'Documents/mapping_function_files/idl_save_files/' ;USER SHOULD REPLACE THIS STRING WITH APPROPRIATE DIRECTORY NAME.
  ion_to_mag = 0
  mag_to_ion_north = 0
  mag_to_ion_south = 0
  rj = 0.d
  loctime = 0.d
  latpt_input = 0.d
  longpt_input = 0.d

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
    print,'ERROR. Please check the mapping type entered. Valid inputs are the following strings: ion_to_mag, mag_to_ion_north, mag_to_ion_south'
    badpt = 1
  endelse
  
  if rj eq 15.d then rj = 15.0001d
  
  latpt = latpt_input
  longpt = longpt_input
  rj_input = rj
  loctime_input = loctime
  
  sslong = sslong2
  
  badpt = 0
  ;test if user-specified point is within the model validity range (15-150 Rj in the magnetosphere), local time given in decimal format 0-24 hours,
  if ion_to_mag eq 1 and (abs(latpt) gt 90d) then badpt = 1
  if ion_to_mag eq 1 and (longpt gt 360d) then longpt = longpt mod 360d
  if ion_to_mag eq 1 and longpt lt 0d then longpt = 360d + (longpt mod 360d)
  if (mag_to_ion_north eq 1 or mag_to_ion_south eq 1) and (rj gt 150d or rj lt 15d or loctime lt 0d or loctime gt 24d) then badpt = 1
  if ion_to_mag ne 1 and mag_to_ion_north ne 1 and mag_to_ion_south ne 1 then badpt = 1
  if rj lt 15d and rj ge 2d and keyword_set(fieldline_tracing) then badpt = 0
  
  model_filename_north = '_jrm33'
  model_filename_south = '_jrm33'
  if keyword_set(model) and keyword_set(fieldline_tracing) eq 0 then begin
    if model eq 'vip4' or model eq 'VIP4' then begin
      model_filename_north = '_vip4'
      model_filename_south = '_vip4'
    endif else if model eq 'gam' or model eq 'GAM' then begin
      model_filename_north = '_gam'
      if mag_to_ion_south eq 1 then begin
        print,'Error. Grodent anomaly model can only be used for mapping to northern hemisphere (not south).'
        badpt = 1
      endif
    endif else if model eq 'vipal' or model eq 'VIPAL' then begin
      model_filename_north = '_vipal'
      model_filename_south = '_vipal'
    endif else if model eq 'jrm09' or model eq 'JRM09' then begin
      model_filename_north = '_jrm09'
      model_filename_south = '_jrm09'
    endif else if model eq 'jrm33' or model eq 'JRM33' then begin
      model_filename_north = '_jrm33'
      model_filename_south = '_jrm33'
    endif else if model eq 'khurana' or model eq 'KHURANA' or model eq 'Khurana' then begin
      model_filename_north = '_kk2009'
      model_filename_south = '_kk2009'
    endif else begin
      print,'ERROR. Please check model name. Valid inputs are: gam, vipal, vip4, jrm09, jrm33, khurana*, and khurana_jrm09*'
      print,'*khurana and khurana_jrm09 model options are only available for the fieldline tracing option'
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
        print,'Error. Grodent anomaly model can only be used for mapping to northern hemisphere (not south).'
        badpt = 1
      endif
    endif else begin
        print,'ERROR. Please check model name. Valid inputs are: gam, vipal, vip4, jrm09, jrm33, khurana, and khurana_jrm09'
        print,'*khurana and khurana_jrm09 model options are only available for the fieldline tracing option'
        badpt = 1
    endelse
  endif else use_fieldline_tracing = 0
  
  if ion_to_mag eq 0 then begin ; test if magnetosphere point given by user is beyond joy et al magnetopause, if yes then return flag value
    pexpanded = 0.039;
    amage = -0.134 + 0.488*pexpanded^(-.25d);
    bmage = -0.581 - 0.225*pexpanded^(-.25d);
    cmage = -0.186 - 0.016*pexpanded^(-.25d);
    dmage = -0.014 + 0.096*pexpanded;
    emage = -0.814 - 0.811*pexpanded;
    fmage = -0.050 + 0.168*pexpanded;
    xplot = dindgen(10000)
    iloop = 0;
    while (iloop lt 10000) do begin
      xplot[iloop] = (iloop-200d)/120d;
      iloop = iloop + 1d;
    endwhile
    loctime_rad = loctime*!dpi/12d;
    r = rj;
    xplot = (-1.0)*r*cos(loctime_rad)/120d;
    y_point = r*sin(loctime_rad);
    if (r gt 200.0) then begin
      r = 200.0;
    endif
    xplot = (-1d)*r*cos(loctime_rad)/120d;d
    bplot = dmage + fmage*xplot;
    aplot = emage;
    cplot = amage + bmage*xplot + cmage*(xplot*xplot);
    yplotplus =  (-1.*bplot + sqrt(bplot*bplot - 4.0*aplot*cplot))/(2.0*aplot); // help
    yplotminus = (-1.*bplot - sqrt(bplot*bplot - 4.0*aplot*cplot))/(2.0*aplot);
    yplotplus = -120.*yplotplus;
    yplotminus = -120.*yplotminus;
    is_inside = 0;
    if (y_point lt yplotplus and y_point gt yplotminus) then begin
      is_inside = 1;
    endif
    if is_inside eq 0 then badpt = 1
    rj = float(rj)
    loctime = float(loctime)
  endif else begin
    latpt = float(latpt)
    longpt = float(longpt)
  endelse
  
  sslong_old = sslong + 360.d mod 360.d
  sslong = 360.d - sslong + 360.d
  sslong = sslong mod 360.d
  
  longpt = 360.0d - longpt + 360.0d; // covert to right handed
  longpt = longpt mod 360.0d
  
  latpt_start = latpt;
  
  xpt = sin((90.d -latpt)*!dpi/180.d)*cos(longpt*!dpi/180.d);
  ypt = sin((90.d -latpt)*!dpi/180.d)*sin(longpt*!dpi/180.d);
  
  nelements = -1;
  filematch = -1;
  
  ;flat
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
          
          x1 = sin(lat1*!dpi/180.d)*cos(long1*!dpi/180.d);
          y1 = sin(lat1*!dpi/180.d)*sin(long1*!dpi/180.d);
          x2 = sin(lat2*!dpi/180.d)*cos(long2*!dpi/180.d);
          y2 = sin(lat2*!dpi/180.d)*sin(long2*!dpi/180.d);
          x3 = sin(lat3*!dpi/180.d)*cos(long3*!dpi/180.d);
          y3 = sin(lat3*!dpi/180.d)*sin(long3*!dpi/180.d);
          x4 = sin(lat4*!dpi/180.d)*cos(long4*!dpi/180.d);
          y4 = sin(lat4*!dpi/180.d)*sin(long4*!dpi/180.d);
          
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
          badfilematch = 1
          
          filematchpts = n_elements(filematch)
          insidept = 0
          if total(filematch) ne -1 then begin
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
          endif   
          
          
          rstepsize = 5.d
          ;===========
          ;For fieldline tracing only - test if point is inside of 15 Rj but outside of 2 Rj
          ;===========
          if insidept eq 0 and use_fieldline_tracing eq 1 then begin
            if sslong ne 0. then filetext = directory_name + 'fieldline_tracing_2to15rj_'+model_filename_north+'_north_cml'+strtrim(sslongfix,2)+'.sav'
            if sslong eq 0. then filetext = directory_name + 'fieldline_tracing_2to15rj_'+model_filename_north+'_north_cml360.sav'
            sslong = sslong1;
            restore,filetext
            x1 = sin(lat1*!dpi/180.d)*cos(long1*!dpi/180.d);
            y1 = sin(lat1*!dpi/180.d)*sin(long1*!dpi/180.d);
            x2 = sin(lat2*!dpi/180.d)*cos(long2*!dpi/180.d);
            y2 = sin(lat2*!dpi/180.d)*sin(long2*!dpi/180.d);
            x3 = sin(lat3*!dpi/180.d)*cos(long3*!dpi/180.d);
            y3 = sin(lat3*!dpi/180.d)*sin(long3*!dpi/180.d);
            x4 = sin(lat4*!dpi/180.d)*cos(long4*!dpi/180.d);
            y4 = sin(lat4*!dpi/180.d)*sin(long4*!dpi/180.d);
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
            badfilematch = 1
            filematchpts = n_elements(filematch)
            insidept = 0
            if total(filematch) ne -1 then begin
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
            endif
            rstepsize = 0.5
          endif
          ;=====
          ;End only inside of 15 Rj part     
          ;=====

          if insidept eq 0 then filematch = -1

          r1_array = r1
          x1_array = x1
          y1_array = y1
          if total(filematch) ne -1 then begin
            badfilematch = 0
            x1 = x1(filematch)
            x2 = x2(filematch)
            x3 = x3(filematch)
            x4 = x4(filematch)
            y1 = y1(filematch)
            y2 = y2(filematch)
            y3 = y3(filematch)
            y4 = y4(filematch)
            
            x1prime = x1 - x1;
            x2prime = x2 - x1;
            x3prime = x3 - x1;
            x4prime = x4 - x1;
            xpt_prime = xpt - x1;
            y1prime = y1 - y1;
            y2prime = y2 - y1;
            y3prime = y3 - y1;
            y4prime = y4 - y1;
            ypt_prime = ypt - y1;
            
            theta = atan2(y3prime,x3prime);
            
            a1 = x1prime*cos(theta) + y1prime*sin(theta);
            a2 = x2prime*cos(theta) + y2prime*sin(theta);
            a3 = x3prime*cos(theta) + y3prime*sin(theta);
            a4 = x4prime*cos(theta) + y4prime*sin(theta);
            apt = xpt_prime*cos(theta) + ypt_prime*sin(theta);
            b1 = -1.0d*x1prime*sin(theta) + y1prime*cos(theta);
            b2 = -1.0d*x2prime*sin(theta) + y2prime*cos(theta);
            b3 = -1.0d*x3prime*sin(theta) + y3prime*cos(theta);
            b4 = -1.0d*x4prime*sin(theta) + y4prime*cos(theta);
            bpt = -1.0d*xpt_prime*sin(theta) + ypt_prime*cos(theta);
            
            c = b2 - ((b4-b2)/(a4-a2))*a2;
            b5 = ((b4-b2)/(a4-a2))*apt + c;
            r_interpol = r1_array(filematch) + bpt*rstepsize/b5;
            
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
          endif
          
          ;if unable to find match for point in the files
          if (filematch eq -1) then begin
            ;test to see if point is inside of the 15 Rj contour
            jloop = n_elements(where(r1_array eq rmin_limit))
            
            jloop = jloop - 1;
            
            c = 0;
            iloop = 0;
            jloop_start = jloop;
            jloop = jloop;
            nvert = jloop+1;
            y1_match = y1_array(where(r1_array eq rmin_limit))
            x1_match = x1_array(where(r1_array eq rmin_limit))
            c = point_inside_polygon(xpt,ypt,x1_match,y1_match)
            inside15 = c;
            
            if (inside15 eq 1) then begin ;maps to outside 150 RJ or Jovian magnetopause
              r_interpol = -999.d
              ltime_interpol = -999.d
            endif else begin; maps to inside of 15 RJ (or 2 RJ for fieldline tracing)
              r_interpol = -998.d
              ltime_interpol = -998.d
            endelse
          endif else begin ;filematch works fine
            ;determines whether final mapping point is inside the magnetosphere or not ; from joy et al 2002
            pexpanded = 0.039;
            amage = -0.134d + 0.488d*pexpanded^(-.25d);
            bmage = -0.581d - 0.225d*pexpanded^(-.25d);
            cmage = -0.186d - 0.016d*pexpanded^(-.25d);
            dmage = -0.014d + 0.096d*pexpanded;
            emage = -0.814d - 0.811d*pexpanded;
            fmage = -0.050d + 0.168d*pexpanded;
            
            xplot = dindgen(10000)
            iloop = 0;
            while (iloop lt 10000d) do begin
              xplot(iloop) = (iloop-200.0d)/120.0d;
              iloop = iloop + 1.0d;
            endwhile
            
            loctime = ltime_interpol*!dpi/12.0d;
            r = r_interpol;
            xplot = -1.0d*r*cos(loctime)/120.d;
            y_point = r*sin(loctime);
            if r gt 200.d then r = 200.d
            
            xplot = -1.0d*r*cos(loctime)/120.d;
            
            bplot = dmage + fmage*xplot;
            aplot = emage;
            cplot = amage + bmage*xplot + cmage*(xplot*xplot);
            yplotplus =  (-1.d*bplot + sqrt(bplot*bplot - 4.0d*aplot*cplot))/(2.0d*aplot); // help
            yplotminus = (-1.d*bplot - sqrt(bplot*bplot - 4.0d*aplot*cplot))/(2.0d*aplot);
            yplotplus = -120.d*yplotplus;
            yplotminus = -120.d*yplotminus;
            
            is_inside = 0;
            if (y_point lt yplotplus and y_point gt yplotminus) then is_inside = 1
            
            if (is_inside eq 0) then begin
              ; beyond 150 RJ or expanded magnetosphere
              r_interpol = -999.d
              ltime_interpol = -999.d
            endif else begin ;misc other error confusing!
              ;r_interpol = -997.d
              ;ltime_interpol = -997.d
            endelse
          endelse ; end bad file match
          
          
          
          
          ;end of loop sslong div by 10
        endif else begin ;if sslong div by 10 != 0, still mapping from northern ionosphere to magnetosphere
;          print,'hi'
          sslong1 = sslong;
          sslong = sslong_old;
          
          sslongmatch_min = dindgen(36)*10.d
          sslongmatch_max = (dindgen(36)+1.d)*10.d
          sslongmatch_file = [360,(indgen(35)+1)*10]
          
          sslongmatch1 = where(sslong ge sslongmatch_min and sslong lt sslongmatch_max)
          match1 = sslongmatch_min(sslongmatch1)
          
          filetext = directory_name + 'sslong'+strtrim(sslongmatch_file(sslongmatch1(0)),2)+model_filename_north+'.sav'
          sslong = sslong1;
          
          if keyword_set(fieldline_tracing) then begin 
            first_mapping = mapping_function_2024('ion_to_mag',sslongmatch_file(sslongmatch1(0)),latpt_input,longpt_input,model,fieldline_tracing)
            second_mapping = mapping_function_2024('ion_to_mag',(sslongmatch_file(sslongmatch1(0))+10) mod 360,latpt_input,longpt_input,model,fieldline_tracing)
          endif else begin
            first_mapping = mapping_function_2024('ion_to_mag',sslongmatch_file(sslongmatch1(0)),latpt_input,longpt_input,model)
            second_mapping = mapping_function_2024('ion_to_mag',(sslongmatch_file(sslongmatch1(0))+10) mod 360,latpt_input,longpt_input,model)
          endelse
                    
          sslong_first = float(sslongmatch_file(sslongmatch1(0)))
          if sslong_first eq 0 or sslong_first eq 360 then sslong_first = 0.
          
          first_x = first_mapping(0)*cos(first_mapping(1)*!dpi/12.d)
          first_y = first_mapping(0)*sin(first_mapping(1)*!dpi/12.d)
          second_x = second_mapping(0)*cos(second_mapping(1)*!dpi/12.d)
          second_y = second_mapping(0)*sin(second_mapping(1)*!dpi/12.d)
          x_interpol = first_x + (second_x-first_x)*((360.-sslong)-sslong_first)/10.d
          y_interpol = first_y + (second_y-first_y)*((360.-sslong)-sslong_first)/10.d
          r_interpol = sqrt(x_interpol^2.d + y_interpol^2.d)
          ltime_interpol = (atan2(y_interpol,x_interpol)*12.d/!dpi + 48.) mod 24.
          if abs(first_mapping(0)) ge 200. or abs(second_mapping(0)) ge 200. then begin
            r_interpol = first_mapping(0)
            ltime_interpol = first_mapping(1)
          endif
          
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
                   
          x1 = sin(lat1*!dpi/180.d)*cos(long1*!dpi/180.d);
          y1 = sin(lat1*!dpi/180.d)*sin(long1*!dpi/180.d);
          x2 = sin(lat2*!dpi/180.d)*cos(long2*!dpi/180.d);
          y2 = sin(lat2*!dpi/180.d)*sin(long2*!dpi/180.d);
          x3 = sin(lat3*!dpi/180.d)*cos(long3*!dpi/180.d);
          y3 = sin(lat3*!dpi/180.d)*sin(long3*!dpi/180.d);
          x4 = sin(lat4*!dpi/180.d)*cos(long4*!dpi/180.d);
          y4 = sin(lat4*!dpi/180.d)*sin(long4*!dpi/180.d);
          
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
          badfilematch = 1
          
          filematchpts = n_elements(filematch)
          insidept = 0
          if total(filematch) ne -1 then begin
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
          endif
          
          rstepsize = 5.d
          ;===========
          ;For fieldline tracing only - test if point is inside of 15 Rj but outside of 2 Rj
          ;===========
          if insidept eq 0 and use_fieldline_tracing eq 1 then begin
            if sslong ne 0. then filetext = directory_name + 'fieldline_tracing_2to15rj_'+model_filename_south+'_south_cml'+strtrim(sslongfix,2)+'.sav'
            if sslong eq 0. then filetext = directory_name + 'fieldline_tracing_2to15rj_'+model_filename_south+'_south_cml360.sav'
            restore,filetext
            sslong = sslong1;
            x1 = sin(lat1*!dpi/180.d)*cos(long1*!dpi/180.d);
            y1 = sin(lat1*!dpi/180.d)*sin(long1*!dpi/180.d);
            x2 = sin(lat2*!dpi/180.d)*cos(long2*!dpi/180.d);
            y2 = sin(lat2*!dpi/180.d)*sin(long2*!dpi/180.d);
            x3 = sin(lat3*!dpi/180.d)*cos(long3*!dpi/180.d);
            y3 = sin(lat3*!dpi/180.d)*sin(long3*!dpi/180.d);
            x4 = sin(lat4*!dpi/180.d)*cos(long4*!dpi/180.d);
            y4 = sin(lat4*!dpi/180.d)*sin(long4*!dpi/180.d);
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
            badfilematch = 1            
            filematchpts = n_elements(filematch)
            insidept = 0
            if total(filematch) ne -1 then begin
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
            endif
            rstepsize = 0.5
          endif
          ;=====
          ;End only inside of 15 Rj part
          ;=====

          if insidept eq 0 then filematch = -1
          
          r1_array = r1
          x1_array = x1
          y1_array = y1
          
          if total(filematch) ne -1 then begin
            badfilematch = 0
            x1 = x1(filematch)
            x2 = x2(filematch)
            x3 = x3(filematch)
            x4 = x4(filematch)
            y1 = y1(filematch)
            y2 = y2(filematch)
            y3 = y3(filematch)
            y4 = y4(filematch)
                        
            x1prime = x1 - x1;
            x2prime = x2 - x1;
            x3prime = x3 - x1;
            x4prime = x4 - x1;
            xpt_prime = xpt - x1;
            y1prime = y1 - y1;
            y2prime = y2 - y1;
            y3prime = y3 - y1;
            y4prime = y4 - y1;
            ypt_prime = ypt - y1;

            theta = atan2(y3prime,x3prime);
            
            a1 = x1prime*cos(theta) + y1prime*sin(theta);
            a2 = x2prime*cos(theta) + y2prime*sin(theta);
            a3 = x3prime*cos(theta) + y3prime*sin(theta);
            a4 = x4prime*cos(theta) + y4prime*sin(theta);
            apt = xpt_prime*cos(theta) + ypt_prime*sin(theta);
            b1 = -1.0d*x1prime*sin(theta) + y1prime*cos(theta);
            b2 = -1.0d*x2prime*sin(theta) + y2prime*cos(theta);
            b3 = -1.0d*x3prime*sin(theta) + y3prime*cos(theta);
            b4 = -1.0d*x4prime*sin(theta) + y4prime*cos(theta);
            bpt = -1.0d*xpt_prime*sin(theta) + ypt_prime*cos(theta);
            
            c = b2 - ((b4-b2)/(a4-a2))*a2;
            b5 = ((b4-b2)/(a4-a2))*apt + c;
            r_interpol = r1_array(filematch) + bpt*rstepsize/b5;
            
            ;all new
            x_array = [x1,x2,x3,x4]
            y_array = [y1,y2,y3,y4]
            r_array = [r1(filematch),r2(filematch),r3(filematch),r4(filematch)]
            lt_array = [ltime1(filematch),ltime2(filematch),ltime3(filematch),ltime4(filematch)]
            min_x = where(x_array eq min(x_array))
            min_y = where(y_array eq min(y_array))
            max_x = where(x_array eq max(x_array))
            max_y = where(y_array eq max(y_array))
            deltax = x_array(max_x) - x_array(min_x)
            deltay = y_array(max_x) - y_array(min_x)

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
            
            if ((dist1 le dist2 and dist1 le dist4) or ((dist3 le dist2 and dist3 le dist4))) then begin
              if (abs(ltime3 - ltime1) lt 12.0d) then begin
                ltime_interpol = ltime1 + ((ltime3 - ltime1)/(a3-a1))*(apt-a1);
              endif else begin
                ltime_interpol = ltime1 + ((24.0d -abs(ltime3 - ltime1))/(a3-a1))*(apt-a1);
              endelse
            endif else begin
              if (abs(ltime4 - ltime2) lt 12.0d) then begin
                ltime_interpol = ltime2 + ((ltime4 - ltime2)/(a4-a2))*(apt-a2);
              endif else begin
                ltime_interpol = ltime2 + ((24.0d -abs(ltime4 - ltime2))/(a4-a2))*(apt-a2);
              endelse
            endelse
          endif
          ;end of loop sslong div by 10
        endif else begin ;if sslong div by 10 != 0 still southern ionosphere to magnetosphere
          sslong1 = sslong;
          sslong = sslong_old;
          
          sslongmatch_min = dindgen(36)*10.d
          sslongmatch_max = (dindgen(36)+1.d)*10.d
          sslongmatch_file = [360,(indgen(35)+1)*10]
          
          sslongmatch1 = where(sslong ge sslongmatch_min and sslong lt sslongmatch_max)
          match1 = sslongmatch_min(sslongmatch1)
          
          filetext = directory_name + 'sslong'+strtrim(sslongmatch_file(sslongmatch1(0)),2)+'s'+model_filename_south+'.sav'
          sslong = sslong1;
          
          if keyword_set(fieldline_tracing) then begin
            first_mapping = mapping_function_2024('ion_to_mag',sslongmatch_file(sslongmatch1(0)),latpt_input,longpt_input,model,fieldline_tracing)
            second_mapping = mapping_function_2024('ion_to_mag',(sslongmatch_file(sslongmatch1(0))+10) mod 360,latpt_input,longpt_input,model,fieldline_tracing)
          endif else begin
            first_mapping = mapping_function_2024('ion_to_mag',sslongmatch_file(sslongmatch1(0)),latpt_input,longpt_input,model)
            second_mapping = mapping_function_2024('ion_to_mag',(sslongmatch_file(sslongmatch1(0))+10) mod 360,latpt_input,longpt_input,model)
          endelse
          
          sslong_first = float(sslongmatch_file(sslongmatch1(0)))
          if sslong_first eq 0 or sslong_first eq 360 then sslong_first = 0.
          
          first_x = first_mapping(0)*cos(first_mapping(1)*!dpi/12.d)
          first_y = first_mapping(0)*sin(first_mapping(1)*!dpi/12.d)
          second_x = second_mapping(0)*cos(second_mapping(1)*!dpi/12.d)
          second_y = second_mapping(0)*sin(second_mapping(1)*!dpi/12.d)
          x_interpol = first_x + (second_x-first_x)*((360.-sslong)-sslong_first)/10.d
          y_interpol = first_y + (second_y-first_y)*((360.-sslong)-sslong_first)/10.d
          r_interpol = sqrt(x_interpol^2.d + y_interpol^2.d)
          ltime_interpol = (atan2(y_interpol,x_interpol)*12.d/!dpi + 48.) mod 24.
          if abs(first_mapping(0)) ge 400. or abs(second_mapping(0)) ge 400. then begin
            r_interpol = first_mapping(0)
            ltime_interpol = first_mapping(1)
          endif
          filematch = 1
          
          
          
        endelse ;end sslong div 10 != 0 loop
        
        
        
        ;if unable to find match for point in the files
        if (filematch eq -1) then begin
        
          ;test to see if point is inside of the 15 Rj contour
          jloop = n_elements(where(r1_array eq rmin_limit))
          
          jloop = jloop - 1;
          
          c = 0;
          iloop = 0;
          jloop_start = jloop;
          jloop = jloop;
          nvert = jloop+1;
          y1_match = y1_array(where(r1_array eq rmin_limit))
          x1_match = x1_array(where(r1_array eq rmin_limit))
          c = point_inside_polygon(xpt,ypt,x1_match,y1_match)
          
          inside15 = c;
          
          if (inside15 eq 1) then begin ;maps to outside 150 RJ or Jovian magnetopause
            r_interpol = -999.d
            ltime_interpol = -999.d
          endif else begin; maps to inside of 15 RJ (or 2 Rj for fieldline tracing)
            r_interpol = -998.d
            ltime_interpol = -998.d
          endelse
        endif else begin ;filematch works fine
          ;determines whether final mapping point is inside the magnetosphere or not ; from joy et al 2002
          pexpanded = 0.039;
          amage = -0.134d + 0.488d*pexpanded^(-.25d);
          bmage = -0.581d - 0.225d*pexpanded^(-.25d);
          cmage = -0.186d - 0.016d*pexpanded^(-.25d);
          dmage = -0.014d + 0.096d*pexpanded;
          emage = -0.814d - 0.811d*pexpanded;
          fmage = -0.050d + 0.168d*pexpanded;
          
          xplot = dindgen(10000)
          iloop = 0;
          while (iloop lt 10000d) do begin
            xplot(iloop) = (iloop-200.0d)/120.0d;
            iloop = iloop + 1.0d;
          endwhile
          
          loctime = ltime_interpol*!dpi/12.0d;
          r = r_interpol;
          xplot = -1.0d*r*cos(loctime)/120.d;
          y_point = r*sin(loctime);
          if r gt 200.d then r = 200.d
          
          xplot = -1.0d*r*cos(loctime)/120.d;
          
          bplot = dmage + fmage*xplot;
          aplot = emage;
          cplot = amage + bmage*xplot + cmage*(xplot*xplot);
          yplotplus =  (-1.d*bplot + sqrt(bplot*bplot - 4.0d*aplot*cplot))/(2.0d*aplot); // help
          yplotminus = (-1.d*bplot - sqrt(bplot*bplot - 4.0d*aplot*cplot))/(2.0d*aplot);
          yplotplus = -120.d*yplotplus;
          yplotminus = -120.d*yplotminus;
          
          is_inside = 0;
          if (y_point lt yplotplus and y_point gt yplotminus) then is_inside = 1
          
          if (is_inside eq 0) then begin
            ; beyond 150 RJ or expanded magnetosphere
            r_interpol = -999.d
            ltime_interpol = -999.d
          endif else begin ;misc other error confusing!
          endelse
        endelse ; end bad file match
        
        
        
        
        
        
        
        
        
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
        
        x1 = r1*cos(ltime1*!dpi/12.0);
        y1 = r1*sin(ltime1*!dpi/12.0);
        x2 = r2*cos(ltime2*!dpi/12.0);
        y2 = r2*sin(ltime2*!dpi/12.0);
        x3 = r3*cos(ltime3*!dpi/12.0);
        y3 = r3*sin(ltime3*!dpi/12.0);
        x4 = r4*cos(ltime4*!dpi/12.0);
        y4 = r4*sin(ltime4*!dpi/12.0);
        
        x1ion = sin(lat1*!dpi/180.d)*cos(long1*!dpi/180.d);
        y1ion = sin(lat1*!dpi/180.d)*sin(long1*!dpi/180.d);
        x2ion = sin(lat2*!dpi/180.d)*cos(long2*!dpi/180.d);
        y2ion = sin(lat2*!dpi/180.d)*sin(long2*!dpi/180.d);
        x3ion = sin(lat3*!dpi/180.d)*cos(long3*!dpi/180.d);
        y3ion = sin(lat3*!dpi/180.d)*sin(long3*!dpi/180.d);
        x4ion = sin(lat4*!dpi/180.d)*cos(long4*!dpi/180.d);
        y4ion = sin(lat4*!dpi/180.d)*sin(long4*!dpi/180.d);
        
        nelements = nelements + 1;
        
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
        
        x1 = x1(filematch);
        x2 = x2(filematch);
        x3 = x3(filematch);
        x4 = x4(filematch);
        y1 = y1(filematch);
        y2 = y2(filematch);
        y3 = y3(filematch);
        y4 = y4(filematch);
        
        x1ion = x1ion(filematch);
        x2ion = x2ion(filematch);
        x3ion = x3ion(filematch);
        x4ion = x4ion(filematch);
        y1ion = y1ion(filematch);
        y2ion = y2ion(filematch);
        y3ion = y3ion(filematch);
        y4ion = y4ion(filematch);
        
        x1prime = x1 - x1;
        x2prime = x2 - x1;
        x3prime = x3 - x1;
        x4prime = x4 - x1;
        xpt_prime = xpt - x1;
        y1prime = y1 - y1;
        y2prime = y2 - y1;
        y3prime = y3 - y1;
        y4prime = y4 - y1;
        ypt_prime = ypt - y1;
        
        theta = atan2(y3prime,x3prime)
        
        a1 = x1prime*cos(theta) + y1prime*sin(theta);
        a2 = x2prime*cos(theta) + y2prime*sin(theta);
        a3 = x3prime*cos(theta) + y3prime*sin(theta);
        a4 = x4prime*cos(theta) + y4prime*sin(theta);
        apt = xpt_prime*cos(theta) + ypt_prime*sin(theta);
        b1 = -1.0d*x1prime*sin(theta) + y1prime*cos(theta);
        b2 = -1.0d*x2prime*sin(theta) + y2prime*cos(theta);
        b3 = -1.0d*x3prime*sin(theta) + y3prime*cos(theta);
        b4 = -1.0d*x4prime*sin(theta) + y4prime*cos(theta);
        bpt = -1.0d*xpt_prime*sin(theta) + ypt_prime*cos(theta);
        
        x1prime = x1ion - x1ion;
        x2prime = x2ion - x1ion;
        x3prime = x3ion - x1ion;
        x4prime = x4ion - x1ion;
        y1prime = y1ion - y1ion;
        y2prime = y2ion - y1ion;
        y3prime = y3ion - y1ion;
        y4prime = y4ion - y1ion;
        
        thetaion = atan2(y3prime,x3prime);
        
        a1ion = x1prime*cos(thetaion) + y1prime*sin(thetaion);
        a2ion = x2prime*cos(thetaion) + y2prime*sin(thetaion);
        a3ion = x3prime*cos(thetaion) + y3prime*sin(thetaion);
        a4ion = x4prime*cos(thetaion) + y4prime*sin(thetaion);
        b1ion = -1.0d*x1prime*sin(thetaion) + y1prime*cos(thetaion);
        b2ion = -1.0d*x2prime*sin(thetaion) + y2prime*cos(thetaion);
        b3ion = -1.0d*x3prime*sin(thetaion) + y3prime*cos(thetaion);
        b4ion = -1.0d*x4prime*sin(thetaion) + y4prime*cos(thetaion);
        
        aption = a1ion + (a3ion-a1ion)*(apt-a1)/(a3-a1);
        bption = b1ion + (b2ion-b1ion)*(bpt-b1)/(b2-b1);
        
        xption = aption*cos(thetaion) - bption*sin(thetaion) + x1ion;
        yption = aption*sin(thetaion) + bption*cos(thetaion) + y1ion;
        
        lat_interpol = 90. - (asin(sqrt(xption*xption + yption*yption)))*180./!dpi; // will need to be fixed for southern hemisphere
        long_interpol = (atan2(yption,xption))*180.d/!dpi;
        
        long_interpol = (360.0 - long_interpol + 360.0) mod 360; // convert to LH
      endif else begin; if sslong % 10 ne 0 still magnetosphere to ionosphere north
        sslong1 = sslong;
        sslong = sslong_old;
        
        sslongmatch_min = dindgen(36)*10.d
        sslongmatch_max = (dindgen(36)+1.d)*10.d
        sslongmatch_file = [360,(indgen(35)+1)*10]
        
        sslongmatch1 = where(sslong ge sslongmatch_min and sslong lt sslongmatch_max)
        match1 = sslongmatch_min(sslongmatch1)
        
        ;REPLACE DIRECTORY NAME HERE
        filetext = directory_name + 'sslong'+strtrim(sslongmatch_file(sslongmatch1(0)),2)+model_filename_north+'.sav'
        sslong = sslong1;
        
        if keyword_set(fieldline_tracing) then begin
          first_mapping = mapping_function_2024('mag_to_ion_north',sslongmatch_file(sslongmatch1(0)),rj_input,loctime_input,model,fieldline_tracing)
          second_mapping = mapping_function_2024('mag_to_ion_north',(sslongmatch_file(sslongmatch1(0))+10) mod 360,rj_input,loctime_input,model,fieldline_tracing)
        endif else begin
          first_mapping = mapping_function_2024('mag_to_ion_north',sslongmatch_file(sslongmatch1(0)),rj_input,loctime_input,model)
          second_mapping = mapping_function_2024('mag_to_ion_north',(sslongmatch_file(sslongmatch1(0))+10) mod 360,rj_input,loctime_input,model)
        endelse
        
        sslong_first = float(sslongmatch_file(sslongmatch1(0)))
        if sslong_first eq 0 or sslong_first eq 360 then sslong_first = 0.
        
        first_x = abs(sin((90.-first_mapping(0))*!dpi/180.))*cos(first_mapping(1)*!dpi/180.d)
        first_y = abs(sin((90.-first_mapping(0))*!dpi/180.))*sin(first_mapping(1)*!dpi/180.d)
        second_x = abs(sin((90.-second_mapping(0))*!dpi/180.))*cos(second_mapping(1)*!dpi/180.d)
        second_y = abs(sin((90.-second_mapping(0))*!dpi/180.))*sin(second_mapping(1)*!dpi/180.d)
        x_interpol = first_x + (second_x-first_x)*((360.-sslong)-sslong_first)/10.d
        y_interpol = first_y + (second_y-first_y)*((360.-sslong)-sslong_first)/10.d
        lat_interpol = 90. - asin(sqrt(x_interpol^2.d + y_interpol^2.d))*180.d/!dpi
        long_interpol = (atan2(y_interpol,x_interpol)*180.d/!dpi + 720.) mod 360.
        if abs(first_mapping(0)) ge 400. or abs(second_mapping(0)) ge 400. then begin
          lat_interpol = first_mapping(0)
          long_interpol = first_mapping(1)
        endif
        
        
        
        
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
        
        nelements = -1;
        filematch = -1;
        
        r1_array6 = r1;
        r2_array6 = r2;
        r3_array6 = r3;
        r4_array6 = r4;
        ltime1_array6 = ltime1;
        ltime2_array6 = ltime2;
        ltime3_array6 = ltime3;
        ltime4_array6 = ltime4;
        
        x1 = r1*cos(ltime1*!dpi/12.0);
        y1 = r1*sin(ltime1*!dpi/12.0);
        x2 = r2*cos(ltime2*!dpi/12.0);
        y2 = r2*sin(ltime2*!dpi/12.0);
        x3 = r3*cos(ltime3*!dpi/12.0);
        y3 = r3*sin(ltime3*!dpi/12.0);
        x4 = r4*cos(ltime4*!dpi/12.0);
        y4 = r4*sin(ltime4*!dpi/12.0);
        
        x1_array6 = x1;
        x2_array6 = x2;
        x3_array6 = x3;
        x4_array6 = x4;
        y1_array6 = y1;
        y2_array6 = y2;
        y3_array6 = y3;
        y4_array6 = y4;
        
        x1ion = sin(lat1*!dpi/180.d)*cos(long1*!dpi/180.d);
        y1ion = sin(lat1*!dpi/180.d)*sin(long1*!dpi/180.d);
        x2ion = sin(lat2*!dpi/180.d)*cos(long2*!dpi/180.d);
        y2ion = sin(lat2*!dpi/180.d)*sin(long2*!dpi/180.d);
        x3ion = sin(lat3*!dpi/180.d)*cos(long3*!dpi/180.d);
        y3ion = sin(lat3*!dpi/180.d)*sin(long3*!dpi/180.d);
        x4ion = sin(lat4*!dpi/180.d)*cos(long4*!dpi/180.d);
        y4ion = sin(lat4*!dpi/180.d)*sin(long4*!dpi/180.d);
        
        x1_arrayion6 = x1ion;
        x2_arrayion6 = x2ion;
        x3_arrayion6 = x3ion;
        x4_arrayion6 = x4ion;
        y1_arrayion6 = y1ion;
        y2_arrayion6 = y2ion;
        y3_arrayion6 = y3ion;
        y4_arrayion6 = y4ion;
        
        nelements = nelements + 1;
        
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
        
        x1 = x1_array6(filematch);
        x2 = x2_array6(filematch);
        x3 = x3_array6(filematch);
        x4 = x4_array6(filematch);
        y1 = y1_array6(filematch);
        y2 = y2_array6(filematch);
        y3 = y3_array6(filematch);
        y4 = y4_array6(filematch);
        
        x1ion = x1_arrayion6(filematch);
        x2ion = x2_arrayion6(filematch);
        x3ion = x3_arrayion6(filematch);
        x4ion = x4_arrayion6(filematch);
        y1ion = y1_arrayion6(filematch);
        y2ion = y2_arrayion6(filematch);
        y3ion = y3_arrayion6(filematch);
        y4ion = y4_arrayion6(filematch);
        
        x1prime = x1 - x1;
        x2prime = x2 - x1;
        x3prime = x3 - x1;
        x4prime = x4 - x1;
        xpt_prime = xpt - x1;
        y1prime = y1 - y1;
        y2prime = y2 - y1;
        y3prime = y3 - y1;
        y4prime = y4 - y1;
        ypt_prime = ypt - y1;
        
        theta = atan2(y3prime,x3prime)
        
        a1 = x1prime*cos(theta) + y1prime*sin(theta);
        a2 = x2prime*cos(theta) + y2prime*sin(theta);
        a3 = x3prime*cos(theta) + y3prime*sin(theta);
        a4 = x4prime*cos(theta) + y4prime*sin(theta);
        apt = xpt_prime*cos(theta) + ypt_prime*sin(theta);
        b1 = -1.0*x1prime*sin(theta) + y1prime*cos(theta);
        b2 = -1.0*x2prime*sin(theta) + y2prime*cos(theta);
        b3 = -1.0*x3prime*sin(theta) + y3prime*cos(theta);
        b4 = -1.0*x4prime*sin(theta) + y4prime*cos(theta);
        bpt = -1.0*xpt_prime*sin(theta) + ypt_prime*cos(theta);
        
        x1prime = x1ion - x1ion;
        x2prime = x2ion - x1ion;
        x3prime = x3ion - x1ion;
        x4prime = x4ion - x1ion;
        y1prime = y1ion - y1ion;
        y2prime = y2ion - y1ion;
        y3prime = y3ion - y1ion;
        y4prime = y4ion - y1ion;
        
        thetaion = atan2(y3prime,x3prime);
        
        a1ion = x1prime*cos(thetaion) + y1prime*sin(thetaion);
        a2ion = x2prime*cos(thetaion) + y2prime*sin(thetaion);
        a3ion = x3prime*cos(thetaion) + y3prime*sin(thetaion);
        a4ion = x4prime*cos(thetaion) + y4prime*sin(thetaion);
        b1ion = -1.0*x1prime*sin(thetaion) + y1prime*cos(thetaion);
        b2ion = -1.0*x2prime*sin(thetaion) + y2prime*cos(thetaion);
        b3ion = -1.0*x3prime*sin(thetaion) + y3prime*cos(thetaion);
        b4ion = -1.0*x4prime*sin(thetaion) + y4prime*cos(thetaion);
        
        aption = a1ion + (a3ion-a1ion)*(apt-a1)/(a3-a1);
        bption = b1ion + (b2ion-b1ion)*(bpt-b1)/(b2-b1);
        
        xption = aption*cos(thetaion) - bption*sin(thetaion) + x1ion;
        yption = aption*sin(thetaion) + bption*cos(thetaion) + y1ion;
        
        lat_interpol = (asin(sqrt(xption*xption + yption*yption)))*180.d/!dpi - 90.; // will need to be changed for southern hemisphere
        long_interpol = (atan2(yption,xption))*180.d/!dpi;
        
        long_interpol = (360.0 - long_interpol + 360.0) mod 360.0; // convert to LH
      endif else begin; if sslong % 10 ne 0, mag to ion south
        sslong1 = sslong;
        sslong = sslong_old;
        
        sslongmatch_min = dindgen(36)*10.d
        sslongmatch_max = (dindgen(36)+1.d)*10.d
        sslongmatch_file = [360,(indgen(35)+1)*10]
        
        sslongmatch1 = where(sslong ge sslongmatch_min and sslong lt sslongmatch_max)
        match1 = sslongmatch_min(sslongmatch1)
        
        filetext = directory_name + 'sslong'+strtrim(sslongmatch_file(sslongmatch1(0)),2)+'s'+model_filename_south+'.sav'
        sslong = sslong1;
        
        if keyword_set(fieldline_tracing) then begin
          first_mapping = mapping_function_2024('mag_to_ion_south',sslongmatch_file(sslongmatch1(0)),rj_input,loctime_input,model,fieldline_tracing)
          second_mapping = mapping_function_2024('mag_to_ion_south',(sslongmatch_file(sslongmatch1(0))+10) mod 360,rj_input,loctime_input,model,fieldline_tracing)
        endif else begin          
          first_mapping = mapping_function_2024('mag_to_ion_south',sslongmatch_file(sslongmatch1(0)),rj_input,loctime_input,model)
          second_mapping = mapping_function_2024('mag_to_ion_south',(sslongmatch_file(sslongmatch1(0))+10) mod 360,rj_input,loctime_input,model)
        endelse
        
        sslong_first = float(sslongmatch_file(sslongmatch1(0)))
        if sslong_first eq 0 or sslong_first eq 360 then sslong_first = 0.
        
        first_x = sin((90.d - abs(first_mapping(0)))*!dpi/180.d)*cos(first_mapping(1)*!dpi/180.d)
        first_y = sin((90.d - abs(first_mapping(0)))*!dpi/180.d)*sin(first_mapping(1)*!dpi/180.d)
        second_x = sin((90.d - abs(second_mapping(0)))*!dpi/180.d)*cos(second_mapping(1)*!dpi/180.d)
        second_y = sin((90.d - abs(second_mapping(0)))*!dpi/180.d)*sin(second_mapping(1)*!dpi/180.d)
        x_interpol = first_x + (second_x-first_x)*((360.-sslong)-sslong_first)/10.d
        y_interpol = first_y + (second_y-first_y)*((360.-sslong)-sslong_first)/10.d
        lat_interpol = asin(sqrt(x_interpol^2.d + y_interpol^2.d))*180.d/!dpi - 90.d
        long_interpol = (atan2(y_interpol,x_interpol)*180.d/!dpi + 720.) mod 360.
        
        if abs(first_mapping(0)) ge 200. or abs(second_mapping(0)) ge 200. then begin
          lat_interpol = first_mapping(0)
          long_interpol = first_mapping(1)
        endif
        
        
        
      endelse ; end mod 10 ne 0
      
      
    endif ;end mag to ion south mapping
    
    
    
  endif
  
  if (ion_to_mag eq 1 and abs(r_interpol) lt 15. and use_fieldline_tracing eq 0) or (ion_to_mag eq 1 and abs(r_interpol) lt 2. and use_fieldline_tracing eq 1) then begin
    r_interpol = 999d
    ltime_interpol = 999d
  endif
  if ion_to_mag eq 1 then return,double([r_interpol,ltime_interpol])
  if ion_to_mag eq 0 then return,double([lat_interpol,long_interpol])
  
  
  
  
  sslong = sslong2
  
  
END