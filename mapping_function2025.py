"""
; Jovian magnetosphere/ionosphere mapping program
; Written by Marissa Vogt, mvogt@psi.edu
; Initially translated into Python from IDL in 2025
;
; PURPOSE OF THIS FUNCTION
; This function allows a user to magnetically map a point between Jupiter's middle magnetosphere and the ionosphere.
; Points can be mapped using the flux equivalence mapping model of Vogt et al. [2011, 2015] or by tracing field lines from a global field model.
;
; HOW TO USE THIS FUNCTION
;   (Note that this function maps one point at a time but user can call it in a loop to map several points at once.)
;   **Before running this program the user should replace the directory_name variable on line 484 below with the appropriate directory for their local machine
;
;   USER INPUTS - REQUIRED
;   * mapping_type: STRING that specifies the type of mapping to do -- from the magnetosphere to the ionosphere, or from the ionosphere (north or south) to the magnetosphere; valid inputs are "ion_to_mag", "mag_to_ion_north", and "mag_to_ion_south"
;   * sslong: the subsolar longitude, in degrees SIII left-handed, 0-360
;   * input_position1 and input_position2: Coordinates in the magnetosphere to be mapped to the ionosphere OR coordinates in the ionosphere to be mapped to the magnetosphere.
;
;   For mapping from the magnetosphere to the ionosphere, input_position1 should be radial distance and input_position2 should be local time.
;     Radial distance should be given in Jovian radii and range from 15 to 150 Rj (flux mapping) or 2-85 Rj (field line tracing). Local time should be given in decimal format (e.g. 12.5 for 12:30) and range from 0 to 24 hours.
;   For mapping from the ionosphere to the magnetosphere, input_position1 should be latitude and input_position2 should be longitude.
;     Latitude should be given in degrees, and should be be negative for southern hemisphere points (range -90 to +90 degrees)
;     Longitude should also be given in degrees, in SIII left-handed, and range is 0-360 degrees.
;
;   OPTIONAL KEYWORDS
;     "MODEL" - This keyword allows the user to choose which internal field model (VIP 4, the Grodent Anomaly Model, VIPAL, JRM09, or JRM33) to use in the mapping
;       The default value is now the JRM33 model (Connerney et al., 2018, GRL) for both north and south as of March 2024.
;       Note: The previous default was JRM09, and before that the default was the Grodent anomaly model in the north and VIP4 in the south, as in Vogt et al. [2011].
;       If the user wishes to use VIP4 (Connerney et al., 1998), for both north and south, the keyword MODEL should be set to 'vip4'
;       If the user wishes to use the Grodent anomaly model (Grodent et al., 2009) in the north, the keyword MODEL should be set to 'gam'
;       If the user wishes to use VIPAL (Hess et al., 2011), for both north and south, the keyword MODEL should be set to 'vipal'
;       If the user wishes to use JRM09 (Connerney et al., 2018), for both north and south, the keyword MODEL may be set to 'jrm09'
;       If the user wishes to use JRM33 (Connerney et al., 2022), for both north and south, the keyword MODEL may be set to 'jrm33'
;     "FIELDLINE_TRACING" - This keyword (added 2016) allows the user to calculate the M-I mapping by tracing fieldlines from a model rather than the Vogt et al. [2011] flux equivalence calculation
;       For fieldline tracing, the user MUST specify a field model. The user can specify that the program should use traced field lines from:
;       * VIP4 (north or south) - set model keyword to "vip4"
;       * Grodent anomaly model (north) - set model keyword to "gam"
;       * VIPAL (north or south) - set model keyword to "vipal"
;       * the unpublished Khurana model (for both north and south) - set model keyword to "khurana" (see http://lasp.colorado.edu/home/mop/graphics/code/ for more info about the model)
;       * the unpublished Khurana model with JRM09 (instead of VIP4) - set model keyword to "khurana_jrm09"
;       * JRM09 (north or south) - set model keyword to "jrm09"
;       * JRM33 (north or south) - set model keyword to "jrm33"
;       The program calculates the mapping using field line tracing results that are pre-calculated (as for the flux equivalence mapping).
;       Note that for VIP4, VIPAL, JRM09, and JRM33 the model results are only valid to ~85 Rj. For the Grodent anomaly model results are only valid to ~65 Rj.
;       All field line tracing results use the Connerney et al. (1981, 1998) current sheet except the Khurana model and Khurana_jrm09.
;
;   EXPECTED OUTPUT
;     This program outputs an array of length 2 containing the apppropriate
;     For mapping magnetosphere->northern/southern ionosphere, the result is an array [latitude,longitude] in degrees
;     For mapping ionosphere->magnetosphere, the result is an array [radial distance,local time] in Rj and hours (0-24)   ;
;     The program returns array of flag values (-999, -998, -996, etc.) if the point cannot be mapped for some reason.
;     Reasons why a point might not be mapped (partial list):
;       * (flag value -998) for mapping from the ionosphere to the magnetosphere, the ionospheric point is equatorward of the Ganymede footprint (15 Rj)
;       * (flag value -999) for mapping from the magnetosphere to the ionosphere, the magnetospheric point is outside of the Joy et al (2002) magnetopause
;       * for mapping with fieldline tracing, if using VIP4, GAM, or VIPAL, point is outside 95 Rj
;       * (flag value -996, +999) the program has an error
;
;
;EXAMPLES - All using 180 degrees subsolar longitude
;   1. For mapping from the northern ionosphere to the magnetosphere using the Vogt et al. (2011, 2015) flux equivalence model, 
;   with input position 60 degrees latitude, 180 degrees longitude:
;     a) using the Vogt et al. flux equivalence calculation with the Grodent anomaly model:
;       >>> mapping_function_2025('ion_to_mag',180.,60.,180.,'gam')
;       [87.09737712545181, 10.276312966473846] #;;; means point maps to ~87 Rj and ~10.2 LT
;     b) using the Vogt et al. flux equivalence calculation with VIP4:
;       >>> mapping_function_2025('ion_to_mag',180.,60.,180.,'vip4')
;       [58.39182405613658, 10.11044368113447]  #;;; means point maps to ~58 Rj and ~10.0 LT 
;     c) using the Vogt et al. flux equivalence calculation with JRM33:
;       >>> mapping_function_2025('ion_to_mag',180.,60.,180.,'jrm33')
;       [80.66383554500378, 10.219790900448057]  #;;; means point maps to ~81 Rj and ~10.2 LT
;       Note that you will get the same result if you call mapping_function_2025('ion_to_mag',180.,60.,180.) since JRM33 is the default.
;
;   2. For mapping from the southern ionosphere to the magnetosphere using the Vogt et al. (2011, 2015) flux equivalence model, 
;     with input position -80 degrees latitude, 90 degrees longitude:
;     a) using the Vogt et al. flux equivalence calculation with VIP4:
;       >>> mapping_function_2025('ion_to_mag',180.,-80.,90.,'vip4')
;       [81.15437910056733, 13.894087896699473]  ;;; means point maps to ~81 Rj and ~13.9 LT
;     b) using the Vogt et al. flux equivalence calculation with VIPAL:
;       >>> mapping_function_2025('ion_to_mag',180.,-80.,90.,'vipal')
;       [84.53854971089659, 14.676701021671343]  ;;; means point maps to ~84 Rj and ~14.7 LT
;
;   3. For mapping from the magnetosphere to northern ionosphere using the Vogt et al. (2011, 2015) flux equivalence model,
;   with input position 88 Rj and 9.8 LT:
;     a) using the Vogt et al. flux equivalence calculation with the Grodent anomaly model:
;       >>> mapping_function_2025('mag_to_ion_north',180.,88.,9.8,'gam')
;       [60.17261324723077, 181.9531950372126]  ;;; means point maps to ~60 degrees latitude and ~182 degrees longitude
;     b) using the Vogt et al. flux equivalence calculation with VIP4:
;       >>> mapping_function_2025('mag_to_ion_north',180.,88.,9.8,'vip4')
;       [62.56412437913066, 179.15705792711594]  ;;; means point maps to ~63 degrees latitude and ~179 degrees longitude
;
;   4. For mapping from the magnetosphere to southern ionosphere using the Vogt et al. (2011, 2015) flux equivalence model,
;   with input position 88 Rj and 9.8 LT:
;     a) using the Vogt et al. flux equivalence calculation with VIPAL:
;       >>> mapping_function_2025('mag_to_ion_south',180.,88.,9.8,'vipal')
;       [-88.10248171861885, 192.693606399617]  ;;; means point maps to roughly -88 degrees latitude and ~193 degrees longitude
;     b) using the Vogt et al. flux equivalence calculation with JRM09:
;       >>> mapping_function_2025('mag_to_ion_south',180.,88.,9.8,'jrm09')
;       [-87.80080754181463, 189.69086791528775]  ;;; means point maps to roughly -88 degrees latitude and ~191 degrees longitude
;     
;   5. For mapping from the ionosphere to the magnetosphere by tracing field lines from a global field model,
;   with input position 57 degrees latitude, 180 degrees longitude:
;     a) Using field lines traced with the Grodent anomaly model as the internal field model:
;       >>> mapping_function_2025('ion_to_mag',180.,57.,180,'gam',fieldline_tracing=True)
;       [36.67355187553826, 11.065309757055015]  ;;; means point maps to ~37 Rj and ~11.1 LT
;     b) Using field lines traced from with JRM09 field model as the internal field model:
;       >>> mapping_function_2025('ion_to_mag',180.,57.,180,'jrm33',fieldline_tracing=True)
;       [23.217493716310603, 10.990355840432812]  ;;; means point maps to ~23 Rj and ~11 LT
;       
;   6. For mapping from the magnetosphere to the northern ionosphere by tracing field lines from a global field model,
;   with input position 47 Rj and 18.2 LT:
;     a) Using field lines traced with VIP4 as the internal field model (with Connerney current sheet):
;       >>> mapping_function_2025('mag_to_ion_north',180.,47.,18.2,'vip4',fieldline_tracing=True)
;       [66.32729873645206, 140.63710492193388]  ;;; means point maps to ~66 degrees latitude and ~141 degrees longitude
;     b) Using field lines traced with the Khurana external field model (with VIP4 as the internal field, default):
;       >>> mapping_function_2025('mag_to_ion_north',180.,47.,18.2,'khurana',fieldline_tracing=True)
;       [66.4218894103459, 139.31335650863275]  ;;; means point maps to ~66 degrees latitude and ~139 degrees longitude
;     c) Using field lines traced with the Khurana external field model and JRM09 internal field:
;       >>> mapping_function_2025('mag_to_ion_north',180.,47.,18.2,'khurana_jrm09',fieldline_tracing=True)
;       [63.80722147877524, 149.44916559654052]  ;;; means point maps to ~64 degrees latitude and ~149 degrees longitude
;
;   7. For mapping from the magnetosphere to the southern ionosphere by tracing field lines from a global field model,
;   with input position 47 Rj and 18.2 LT:
;     a) Using field lines traced with the VIPAL internal field model:
;       >>> mapping_function_2025('mag_to_ion_south',180.,47.,18.2,'vipal',fieldline_tracing=True)
;       [-69.91910731406122, 68.30960191789723]  ;;; means point maps to roughly -70 degrees latitude, 68 degrees longitude
;     b) Using field lines traced with the JRM09 internal field model:
;       >>> mapping_function_2025('mag_to_ion_south',180.,47.,18.2,'jrm09',fieldline_tracing=True)
;       [-69.74139804140242, 69.86496638780068]  ;;; means point maps to roughly -70 degrees latitude, 70 degrees longitude
;
; The program reads in data from files called sslong0_gam.txt, sslong10_jrm09.txt, etc., which should be placed in an appropriate directory on the user's machine.
; These files can be downloaded from bit.ly/mappingfiles2024 (updated links will be posted on marissavogt.com/mapping and can be obtained by emailing mvogt@psi.edu)
; ***Before running this program the user will have to replace the directory_name variable in the code (see line 484).
;
; ACKNOWLEDGMENTS
;   Thanks to Bertrand Bonfond, Benjamin Palmaerts, and Andrew Steffl for helping to test early versions of this code and the online mapping tool.
;   Very special thanks to Masafumi Imai for assistance implementing the JRM09 field model.
;   Thanks to Rob Wilson for assistance with code speed ups and other improvements.
;
; WHEN USING RESULTS FROM THIS CODE IN A PAPER OR PRESENTATION, PLEASE CITE
;   Vogt, M. F., M. G. Kivelson, K. K. Khurana, R. J. Walker, B. Bonfond, D. Grodent, and A. Radioti (2011),
;   Improved mapping of Jupiter’s auroral features to magnetospheric sources,
;   J. Geophys. Res., 116, A03220, doi:10.1029/2010JA016148.
;  and
;   Vogt, M. F., E. J. Bunce, M. G. Kivelson, K. K. Khurana, R. J. Walker, A. Radioti, B. Bonfond, and D. Grodent (2015),
;   Magnetosphere-ionosphere mapping at Jupiter: Quantifying the effects of using different internal field models,
;   J. Geophys. Res. Space Physics, 10.1002/2014JA020729
;
; Most recent changes:
; (September 2025)
; * Initial Python version developed (following major rehaul of IDL code format)
"""

import math
import numpy as np
import pandas
from shapely.geometry import Point, Polygon

def is_inside_joy_magnetopause(rj,loctime_rad):
    pexpanded = 0.039
    p_root = pow(pexpanded,-0.25)
    amage = -0.134 + 0.488*p_root
    bmage = -0.581 - 0.225*p_root
    cmage = -0.186 - 0.016*p_root
    dmage = -0.014 + 0.096*pexpanded
    emage = -0.814 - 0.811*pexpanded
    fmage = -0.050 + 0.168*pexpanded
    r = rj
    y_point = r*math.sin(loctime_rad)
    if (r > 200.0):
        r = 200.0
    xplot = (-1.0)*r*math.cos(loctime_rad)/120.0
    bplot = dmage + fmage*xplot
    aplot = emage
    cplot = amage + bmage*xplot + cmage*(xplot*xplot)
    sqrt_bsq_minus4ac = math.sqrt(bplot*bplot - 4.0*aplot*cplot)
    yplotplus =(-1.0*bplot + sqrt_bsq_minus4ac)/(2.0*aplot)
    yplotminus = (-1.0*bplot - sqrt_bsq_minus4ac)/(2.0*aplot)
    yplotplus = -120.0*yplotplus
    yplotminus = -120.0*yplotminus
    is_inside = 0
    if (y_point <= yplotplus and y_point >= yplotminus):
        is_inside = 1
    return is_inside


def identify_ion_to_mag_match(xpt,ypt,x1,x2,x3,x4,y1,y2,y3,y4):
    xmin = x1*0.0
    xmax = x1*0.0
    ymin = x1*0.0
    ymax = x1*0.0

    for j in range(len(x1)):
        xmin[j] = min([x1[j], x2[j], x3[j], x4[j]]);
        xmax[j] = max([x1[j], x2[j], x3[j], x4[j]]);
        ymin[j] = min([y1[j], y2[j], y3[j], y4[j]]);
        ymax[j] = max([y1[j], y2[j], y3[j], y4[j]]);

    filematch = np.where((xmin <= xpt) & (xmax >= xpt) & (ymin <= ypt) & (ymax >= ypt))
    emptymatch = filematch[0].size == 0

    filematchpts = filematch[0].size
    filematch_array = filematch[0]
    insidept = 0
    if emptymatch == False:
        insidept = 0
        k = 0
        while k < filematchpts:
            #insidept1 = point_inside_polygon[xpt,ypt,[x1[filematch[k]],x3[filematch[k]],x4[filematch[k]],x2[filematch[k]]],[y1[filematch[k]],y3[filematch[k]],y4[filematch[k]],y2[filematch[k]]]]
            #insidept2 = point_inside_polygon[xpt,ypt,[x1[filematch[k]],x4[filematch[k]],x3[filematch[k]],x2[filematch[k]]],[y1[filematch[k]],y4[filematch[k]],y3[filematch[k]],y2[filematch[k]]]]
            polygon_coords1 = [(x1[filematch_array[k]], y1[filematch_array[k]]), (x3[filematch_array[k]], y3[filematch_array[k]]), (x4[filematch_array[k]], y4[filematch_array[k]]), (x2[filematch_array[k]], y2[filematch_array[k]])]
            polygon_coords2 = [(x1[filematch_array[k]], y1[filematch_array[k]]), (x4[filematch_array[k]], y4[filematch_array[k]]), (x3[filematch_array[k]], y3[filematch_array[k]]), (x2[filematch_array[k]], y2[filematch_array[k]])]
            polygon1 = Polygon(polygon_coords1)
            polygon2 = Polygon(polygon_coords2)
            # Check if the point is within the polygon
            point = Point(xpt,ypt)
            is_inside1 = polygon1.contains(point)
            is_inside2 = polygon2.contains(point)          
            insidept = (is_inside1 or is_inside2)
            if insidept == True:
                filematch = filematch_array[k]
                k = k+10000
            k+=1

    if insidept == 0:
        filematch = -1

    return filematch


def mapping_interpolation_ion_to_mag(xpt,ypt,filematch,rstepsize,r1,ltime1,lat1,long1,r2,ltime2,lat2,long2,r3,ltime3,lat3,long3,r4,ltime4,lat4,long4,x1,x2,x3,x4,y1,y2,y3,y4):

    x1prime = x1[filematch] - x1[filematch]
    x2prime = x2[filematch] - x1[filematch]
    x3prime = x3[filematch] - x1[filematch]
    x4prime = x4[filematch] - x1[filematch]
    xpt_prime = xpt - x1[filematch]
    y1prime = y1[filematch] - y1[filematch]
    y2prime = y2[filematch] - y1[filematch]
    y3prime = y3[filematch] - y1[filematch]
    y4prime = y4[filematch] - y1[filematch]
    ypt_prime = ypt - y1[filematch]

    theta = math.atan2(y3prime,x3prime)

    cos_theta = math.cos(theta)
    sin_theta = math.sin(theta)
    a1 = x1prime*cos_theta + y1prime*sin_theta
    a2 = x2prime*cos_theta + y2prime*sin_theta
    a3 = x3prime*cos_theta + y3prime*sin_theta
    a4 = x4prime*cos_theta + y4prime*sin_theta
    apt = xpt_prime*cos_theta + ypt_prime*sin_theta
    b1 = -1.0*x1prime*sin_theta + y1prime*cos_theta
    b2 = -1.0*x2prime*sin_theta + y2prime*cos_theta
    b3 = -1.0*x3prime*sin_theta + y3prime*cos_theta
    b4 = -1.0*x4prime*sin_theta + y4prime*cos_theta
    bpt = -1.0*xpt_prime*sin_theta + ypt_prime*cos_theta

    c = b2 - ((b4-b2)/(a4-a2))*a2
    b5 = ((b4-b2)/(a4-a2))*apt + c
    r_interpol = r1[filematch] + bpt*rstepsize/b5

    dist1 = (apt-a1)*(apt-a1) + (bpt-b1)*(bpt-b1)
    dist2 = (apt-a2)*(apt-a2) + (bpt-b2)*(bpt-b2)
    dist3 = (apt-a3)*(apt-a3) + (bpt-b3)*(bpt-b3)
    dist4 = (apt-a4)*(apt-a4) + (bpt-b4)*(bpt-b4)

    ltime1 = ltime1[filematch]
    ltime2 = ltime2[filematch]
    ltime3 = ltime3[filematch]
    ltime4 = ltime4[filematch]

    if ltime1 < 0.0:
        ltime1 = ltime1 + 24.0
    if ltime2 < 0.0:
        ltime2 = ltime2 + 24.0
    if ltime3 < 0.0:
        ltime3 = ltime3 + 24.0
    if ltime4 < 0.0:
        ltime4 = ltime4 + 24.0

    if (abs(ltime3 - ltime1) < 12.0):
        ltime_interpol = ltime1 + ((ltime3 - ltime1)/(a3-a1))*(apt-a1)
    else:
        ltime_interpol = ltime1 + ((24.0 -abs(ltime3 - ltime1))/(a3-a1))*(apt-a1)
    ltime_interpol = ltime_interpol % 24.0

    return [r_interpol,ltime_interpol]



def mag_to_ion_mapping(north_or_south,xpt,ypt,r1,ltime1,lat1,long1,r2,ltime2,lat2,long2,r3,ltime3,lat3,long3,r4,ltime4,lat4,long4):

    lt_to_rad = math.pi/12.0
    #;Rotate ionospheric contour lat/lon positions into cartesian coordinates for easy interpolation
    ltime1_rad = ltime1*lt_to_rad
    ltime2_rad = ltime2*lt_to_rad
    ltime3_rad = ltime3*lt_to_rad
    ltime4_rad = ltime4*lt_to_rad
    x1 = r1*np.cos(ltime1_rad)
    y1 = r1*np.sin(ltime1_rad)
    x2 = r2*np.cos(ltime2_rad)
    y2 = r2*np.sin(ltime2_rad)
    x3 = r3*np.cos(ltime3_rad)
    y3 = r3*np.sin(ltime3_rad)
    x4 = r4*np.cos(ltime4_rad)
    y4 = r4*np.sin(ltime4_rad)

    lat1_rad = np.radians(lat1)
    lat2_rad = np.radians(lat2)
    lat3_rad = np.radians(lat3)
    lat4_rad = np.radians(lat4)
    long1_rad = np.radians(long1)
    long2_rad = np.radians(long2)
    long3_rad = np.radians(long3)
    long4_rad = np.radians(long4)
    x1ion = np.sin(lat1_rad)*np.cos(long1_rad)
    y1ion = np.sin(lat1_rad)*np.sin(long1_rad)
    x2ion = np.sin(lat2_rad)*np.cos(long2_rad)
    y2ion = np.sin(lat2_rad)*np.sin(long2_rad)
    x3ion = np.sin(lat3_rad)*np.cos(long3_rad)
    y3ion = np.sin(lat3_rad)*np.sin(long3_rad)
    x4ion = np.sin(lat4_rad)*np.cos(long4_rad)
    y4ion = np.sin(lat4_rad)*np.sin(long4_rad)

    xmin = x1*0.0
    xmax = x1*0.0
    ymin = x1*0.0
    ymax = x1*0.0

    for j in range(len(x1)):
        xmin[j] = min([x1[j], x2[j], x3[j], x4[j]]);
        xmax[j] = max([x1[j], x2[j], x3[j], x4[j]]);
        ymin[j] = min([y1[j], y2[j], y3[j], y4[j]]);
        ymax[j] = max([y1[j], y2[j], y3[j], y4[j]]);

    filematch = np.where((xmin <= xpt) & (xmax >= xpt) & (ymin <= ypt) & (ymax >= ypt))
    emptymatch = filematch[0].size == 0

    filematchpts = filematch[0].size
    filematch_array = filematch[0]
    insidept = 0
    
    if emptymatch == False:
        insidept = 0
        k = 0
        while k < filematchpts:
            
            #insidept1 = point_inside_polygon[xpt,ypt,[x1[filematch[k]],x3[filematch[k]],x4[filematch[k]],x2[filematch[k]]],[y1[filematch[k]],y3[filematch[k]],y4[filematch[k]],y2[filematch[k]]]]
            #insidept2 = point_inside_polygon[xpt,ypt,[x1[filematch[k]],x4[filematch[k]],x3[filematch[k]],x2[filematch[k]]],[y1[filematch[k]],y4[filematch[k]],y3[filematch[k]],y2[filematch[k]]]]
            polygon_coords1 = [(x1[filematch_array[k]], y1[filematch_array[k]]), (x3[filematch_array[k]], y3[filematch_array[k]]), (x4[filematch_array[k]], y4[filematch_array[k]]), (x2[filematch_array[k]], y2[filematch_array[k]])]
            polygon_coords2 = [(x1[filematch_array[k]], y1[filematch_array[k]]), (x4[filematch_array[k]], y4[filematch_array[k]]), (x3[filematch_array[k]], y3[filematch_array[k]]), (x2[filematch_array[k]], y2[filematch_array[k]])]
            polygon1 = Polygon(polygon_coords1)
            polygon2 = Polygon(polygon_coords2)
            # Check if the point is within the polygon
            point = Point(xpt,ypt)
            is_inside1 = polygon1.contains(point)
            is_inside2 = polygon2.contains(point)          
            insidept = (is_inside1 or is_inside2)
            if insidept == True:
                filematch = filematch_array[k]
                k = k+10000
            k+=1
        if insidept == False & filematchpts > 1:
            k = 0
            while k < filematchpts:
                polygon_coords = [(x1[filematch_array[k]], y1[filematch_array[k]]), (x3[filematch_array[k]], y3[filematch_array[k]]), (x4[filematch_array[k]], y4[filematch_array[k]]), (x2[filematch_array[k]], y2[filematch_array[k]])]
                polygon = Polygon(polygon_coords)
                # Check if the point is within the polygon
                point = Point(xpt+1e-8,ypt+1e-8)
                is_inside = polygon.contains(point)
                #insidept = point_inside_polygon[xpt+1e-8,ypt+1e-8,[x1[filematch[k]],x3[filematch[k]],x4[filematch[k]],x2[filematch[k]]],[y1[filematch[k]],y3[filematch[k]],y4[filematch[k]],y2[filematch[k]]]]
                if is_inside == True:
                    filematch = filematch_array[k]
                    k = k+10000
                k = k+1             

    if insidept == 0:
        filematch = -1

    x1 = x1[filematch]
    x2 = x2[filematch]
    x3 = x3[filematch]
    x4 = x4[filematch]
    y1 = y1[filematch]
    y2 = y2[filematch]
    y3 = y3[filematch]
    y4 = y4[filematch]

    x1ion = x1ion[filematch]
    x2ion = x2ion[filematch]
    x3ion = x3ion[filematch]
    x4ion = x4ion[filematch]
    y1ion = y1ion[filematch]
    y2ion = y2ion[filematch]
    y3ion = y3ion[filematch]
    y4ion = y4ion[filematch]

    x1prime = x1 - x1
    x2prime = x2 - x1
    x3prime = x3 - x1
    x4prime = x4 - x1
    xpt_prime = xpt - x1
    y1prime = y1 - y1
    y2prime = y2 - y1
    y3prime = y3 - y1
    y4prime = y4 - y1
    ypt_prime = ypt - y1

    theta = math.atan2(y3prime,x3prime)

    cos_theta = math.cos(theta)
    sin_theta = math.sin(theta)
    a1 = x1prime*cos_theta + y1prime*sin_theta
    a2 = x2prime*cos_theta + y2prime*sin_theta
    a3 = x3prime*cos_theta + y3prime*sin_theta
    a4 = x4prime*cos_theta + y4prime*sin_theta
    apt = xpt_prime*cos_theta + ypt_prime*sin_theta
    b1 = -1.0*x1prime*sin_theta + y1prime*cos_theta
    b2 = -1.0*x2prime*sin_theta + y2prime*cos_theta
    b3 = -1.0*x3prime*sin_theta + y3prime*cos_theta
    b4 = -1.0*x4prime*sin_theta + y4prime*cos_theta
    bpt = -1.0*xpt_prime*sin_theta + ypt_prime*cos_theta

    x1prime = x1ion - x1ion
    x2prime = x2ion - x1ion
    x3prime = x3ion - x1ion
    x4prime = x4ion - x1ion
    y1prime = y1ion - y1ion
    y2prime = y2ion - y1ion
    y3prime = y3ion - y1ion
    y4prime = y4ion - y1ion

    thetaion = math.atan2(y3prime,x3prime);

    cos_thetaion = math.cos(thetaion)
    sin_thetaion = math.sin(thetaion)
    a1ion = x1prime*cos_thetaion + y1prime*sin_thetaion
    a2ion = x2prime*cos_thetaion + y2prime*sin_thetaion
    a3ion = x3prime*cos_thetaion + y3prime*sin_thetaion
    a4ion = x4prime*cos_thetaion + y4prime*sin_thetaion
    b1ion = -1.0*x1prime*sin_thetaion + y1prime*cos_thetaion
    b2ion = -1.0*x2prime*sin_thetaion + y2prime*cos_thetaion
    b3ion = -1.0*x3prime*sin_thetaion + y3prime*cos_thetaion
    b4ion = -1.0*x4prime*sin_thetaion + y4prime*cos_thetaion

    aption = a1ion + (a3ion-a1ion)*(apt-a1)/(a3-a1)
    bption = b1ion + (b2ion-b1ion)*(bpt-b1)/(b2-b1)

    xption = aption*cos_thetaion - bption*sin_thetaion + x1ion
    yption = aption*sin_thetaion + bption*cos_thetaion + y1ion

    if north_or_south == 'north':
        lat_interpol = 90.0 - (math.asin(math.sqrt(xption*xption + yption*yption)))*180.0/math.pi #// will need to be fixed for southern hemisphere
    else:
        lat_interpol = (math.asin(math.sqrt(xption*xption + yption*yption)))*180.0/math.pi - 90.0 #// will need to be changed for southern hemisphere
    long_interpol = (math.atan2(yption,xption))*180.0/math.pi

    long_interpol = (360.0 - long_interpol + 360.0) % 360#; // convert to LH

    return [lat_interpol,long_interpol]




def mapping_function_2025(mapping_type,sslong2,input_position1,input_position2,model='jrm33',fieldline_tracing=False):

    directory_name = '/Users/marissav/Documents/mapping_function_files/' #;USER SHOULD REPLACE THIS STRING WITH APPROPRIATE DIRECTORY NAME.
    ion_to_mag = 0
    mag_to_ion_north = 0
    mag_to_ion_south = 0
    rj = 0.0
    loctime = 0.0
    latpt_input = 0.0
    longpt_input = 0.0
    deg_to_rad = math.pi/180.0
    lt_to_rad = math.pi/12.0

    badpt = 0
    if mapping_type == 'ion_to_mag':
        ion_to_mag = 1
        latpt_input = input_position1
        longpt_input = input_position2
    elif mapping_type == 'mag_to_ion_north':
        mag_to_ion_north = 1
        rj = input_position1
        loctime = input_position2
    elif mapping_type == 'mag_to_ion_south':
        mag_to_ion_south = 1
        rj = input_position1
        loctime = input_position2
    else:
        print,'ERROR.0 Please check the mapping type entered. Valid inputs are the following strings: ion_to_mag, mag_to_ion_north, mag_to_ion_south'
        badpt = 1

    if rj == 15.0:
        rj = 15.0001

    latpt = latpt_input
    longpt = longpt_input
    rj_input = rj
    loctime_input = loctime

    sslong = sslong2

    #;test if user-specified point is within the model validity range (15-150 Rj in the magnetosphere), local time given in decimal format 0-24 hours,
    if ion_to_mag == 1:
        if (abs(latpt) > 90.0):
            badpt = 1
        if (longpt > 360.0):
            longpt = longpt % 360.0
        if longpt < 0.0: longpt = 360.0 + (longpt % 360.0)
    else:
        if (ion_to_mag != 1) & (mag_to_ion_north != 1) & (mag_to_ion_south != 1):
            badpt = 1
    if ((mag_to_ion_north == 1) | (mag_to_ion_south == 1)) & ((rj > 150.0) | (rj < 15.0) | (loctime < 0.0) | (loctime > 24.0)):
        badpt = 1
    if (rj < 15.0) & (rj >= 2.0) & (fieldline_tracing != False):
        badpt = 0

    model_filename_north = '_jrm33'
    model_filename_south = '_jrm33'
    if (str.lower(model) == 'jrm33') & (fieldline_tracing == False):
        model_filename_north = '_jrm33'
        model_filename_south = '_jrm33'
    elif fieldline_tracing == False:
        if str.lower(model) == 'vip4':
            model_filename_north = '_vip4'
            model_filename_south = '_vip4'
        elif str.lower(model) == 'gam':
            model_filename_north = '_gam'
            if mag_to_ion_south == 1:
                print,'Error. Grodent anomaly model can only be used for mapping to northern hemisphere (not south).'
                badpt = 1
        elif str.lower(model) == 'vipal':
            model_filename_north = '_vipal'
            model_filename_south = '_vipal'
        elif str.lower(model) == 'jrm09':
            model_filename_north = '_jrm09'
            model_filename_south = '_jrm09'
        else:
            print,'ERROR. Please check model name. Valid inputs are: gam, vipal, vip4, jrm09, and jrm33. The khurana and khurana_jrm09 model options are only available for the fieldline tracing option.'
            badpt = 1

    rmin_limit = 15.0

    if fieldline_tracing == True:
        rmin_limit = 2.0
        if str.lower(model) == 'khurana':
            model_filename_north = 'kk2009'
            model_filename_south = 'kk2009'
        elif str.lower(model) == 'khurana_jrm09':
            model_filename_north = 'kk2009ext_jrm09int'
            model_filename_south = 'kk2009ext_jrm09int'
        elif str.lower(model) == 'vip4':
            model_filename_north = 'vip4'
            model_filename_south = 'vip4'
        elif str.lower(model) == 'vipal':
            model_filename_north = 'vipal'
            model_filename_south = 'vipal'
        elif str.lower(model) == 'jrm09':
            model_filename_north = 'jrm09'
            model_filename_south = 'jrm09'
        elif str.lower(model) == 'jrm33':
            model_filename_north = 'jrm33'
            model_filename_south = 'jrm33'
        elif str.lower(model) == 'gam':
            model_filename_north = 'gam'
            if mag_to_ion_south == 1:
                print,'Error. Grodent anomaly model can only be used for mapping to northern hemisphere (not south).'
                badpt = 1
        else:
            print,'ERROR. Please check model name. Valid inputs are: gam, vipal, vip4, jrm09, jrm33, khurana, and khurana_jrm09.'
            badpt = 1


    #; test if magnetosphere point given by user is beyond joy et al magnetopause, if yes then return flag value
    if ion_to_mag == 0:
        loctime_rad = loctime*lt_to_rad
        is_inside = is_inside_joy_magnetopause(rj,loctime_rad)
        if is_inside == 0: badpt = 1

    sslong_old = (sslong + 360.0) % 360.0
    sslong = 360.0 - sslong + 360.0
    sslong = sslong % 360.0

    nelements = -1
    filematch = -1

    r_interpol = -996.0
    ltime_interpol = -996.0
    lat_interpol = -996.0
    long_interpol = -996.0

    #fieldline tracing for VIP4, GAM, VIPAL not valid outside 100 Rj
    if fieldline_tracing & (model_filename_north != '_kk2009') & (rj > 85.0):
        badpt = 1

    #===========
    #begin main mapping loop
    #===========
    if badpt == 0:
        if ion_to_mag == 1:# ; if user wants to map from the ionosphere to the magnetosphere
            longpt = 360.0 - longpt + 360.0# // covert to right handed
            longpt = longpt % 360.0
            latpt_start = latpt
            xpt = math.sin((90.0 -latpt)*deg_to_rad)*math.cos(longpt*deg_to_rad);
            ypt = math.sin((90.0 -latpt)*deg_to_rad)*math.sin(longpt*deg_to_rad);
            #;loop northern hemisphere
            if (latpt >= 0.0):
                latpt = 90.0 - latpt#; // convert to colatitude
                if sslong % 10.0 == 0.0:#; mapping from northern ionosphere to magnetosphere
                    sslong1 = sslong;
                    sslong = sslong_old;
                    sslongfix = int(sslong)
                    sslong_str = str(sslongfix)

                    if (sslong != 0.0) & (fieldline_tracing == False):
                        filetext = directory_name + 'sslong'+sslong_str+model_filename_north+'.txt'
                    if (sslong == 0.0) & (fieldline_tracing == False):
                        filetext = directory_name + 'sslong360'+model_filename_north+'.txt'
                    if (sslong != 0.0) & (fieldline_tracing == True):
                        filetext = directory_name + 'fieldline_tracing_all_r_'+model_filename_north+'_north_cml'+sslong_str+'.txt'
                    if (sslong == 0.0) & (fieldline_tracing == True):
                        filetext = directory_name + 'fieldline_tracing_all_r_'+model_filename_north+'_north_cml360.txt'

                    sslong = sslong1;

                    full_file = pandas.read_table(filetext,delim_whitespace=True, header=None)
                    r1 = np.array(full_file.iloc[:,0].values.tolist())
                    ltime1 = np.array(full_file.iloc[:,1].values.tolist())
                    lat1 = np.array(full_file.iloc[:,2].values.tolist())
                    long1 = np.array(full_file.iloc[:,3].values.tolist())
                    r2 = np.array(full_file.iloc[:,4].values.tolist())
                    ltime2 = np.array(full_file.iloc[:,5].values.tolist())
                    lat2 = np.array(full_file.iloc[:,6].values.tolist())
                    long2 = np.array(full_file.iloc[:,7].values.tolist())
                    if fieldline_tracing == True: #the fieldline tracing files incorrectly flipped the third and fourth position column sets (R, ltime, lat, and lon)
                        r4 = np.array(full_file.iloc[:,8].values.tolist())
                        ltime4 = np.array(full_file.iloc[:,9].values.tolist())
                        lat4 = np.array(full_file.iloc[:,10].values.tolist())
                        long4 = np.array(full_file.iloc[:,11].values.tolist())
                        r3 = np.array(full_file.iloc[:,12].values.tolist())
                        ltime3 = np.array(full_file.iloc[:,13].values.tolist())
                        lat3 = np.array(full_file.iloc[:,14].values.tolist())
                        long3 = np.array(full_file.iloc[:,15].values.tolist())
                    else:
                        r3 = np.array(full_file.iloc[:,8].values.tolist())
                        ltime3 = np.array(full_file.iloc[:,9].values.tolist())
                        lat3 = np.array(full_file.iloc[:,10].values.tolist())
                        long3 = np.array(full_file.iloc[:,11].values.tolist())
                        r4 = np.array(full_file.iloc[:,12].values.tolist())
                        ltime4 = np.array(full_file.iloc[:,13].values.tolist())
                        lat4 = np.array(full_file.iloc[:,14].values.tolist())
                        long4 = np.array(full_file.iloc[:,15].values.tolist())

                    lat1_rad = lat1*math.pi/180.0
                    lat2_rad = lat2*math.pi/180.0
                    lat3_rad = lat3*math.pi/180.0
                    lat4_rad = lat4*math.pi/180.0
                    long1_rad = long1*math.pi/180.0
                    long2_rad = long2*math.pi/180.0
                    long3_rad = long3*math.pi/180.0
                    long4_rad = long4*math.pi/180.0
                    x1 = np.sin(lat1_rad)*np.cos(long1_rad)
                    y1 = np.sin(lat1_rad)*np.sin(long1_rad)
                    x2 = np.sin(lat2_rad)*np.cos(long2_rad)
                    y2 = np.sin(lat2_rad)*np.sin(long2_rad)
                    x3 = np.sin(lat3_rad)*np.cos(long3_rad)
                    y3 = np.sin(lat3_rad)*np.sin(long3_rad)
                    x4 = np.sin(lat4_rad)*np.cos(long4_rad)
                    y4 = np.sin(lat4_rad)*np.sin(long4_rad)
    
                    filematch = identify_ion_to_mag_match(xpt,ypt,x1,x2,x3,x4,y1,y2,y3,y4)                    

                    rstepsize = 5.0
                    #;===========
                    #;For fieldline tracing only - test if point is inside of 15 Rj but outside of 2 Rj
                    #;===========
                    if (filematch == -1) & (fieldline_tracing == True):
                        if (sslong != 0.0):
                            filetext = directory_name + 'fieldline_tracing_2to15rj_'+model_filename_north+'_north_cml'+sslong_str+'.txt'
                        elif sslong == 0.0:
                            filetext = directory_name + 'fieldline_tracing_2to15rj_'+model_filename_north+'_north_cml360.txt'
                        sslong = sslong1

                        full_file = pandas.read_table(filetext,delim_whitespace=True, header=None)
                        r1 = np.array(full_file.iloc[:,0].values.tolist())
                        ltime1 = np.array(full_file.iloc[:,1].values.tolist())
                        lat1 = np.array(full_file.iloc[:,2].values.tolist())
                        long1 = np.array(full_file.iloc[:,3].values.tolist())
                        r2 = np.array(full_file.iloc[:,4].values.tolist())
                        ltime2 = np.array(full_file.iloc[:,5].values.tolist())
                        lat2 = np.array(full_file.iloc[:,6].values.tolist())
                        long2 = np.array(full_file.iloc[:,7].values.tolist())
                        if fieldline_tracing == True:
                            r4 = np.array(full_file.iloc[:,8].values.tolist())
                            ltime4 = np.array(full_file.iloc[:,9].values.tolist())
                            lat4 = np.array(full_file.iloc[:,10].values.tolist())
                            long4 = np.array(full_file.iloc[:,11].values.tolist())
                            r3 = np.array(full_file.iloc[:,12].values.tolist())
                            ltime3 = np.array(full_file.iloc[:,13].values.tolist())
                            lat3 = np.array(full_file.iloc[:,14].values.tolist())
                            long3 = np.array(full_file.iloc[:,15].values.tolist())
                        else:
                            r3 = np.array(full_file.iloc[:,8].values.tolist())
                            ltime3 = np.array(full_file.iloc[:,9].values.tolist())
                            lat3 = np.array(full_file.iloc[:,10].values.tolist())
                            long3 = np.array(full_file.iloc[:,11].values.tolist())
                            r4 = np.array(full_file.iloc[:,12].values.tolist())
                            ltime4 = np.array(full_file.iloc[:,13].values.tolist())
                            lat4 = np.array(full_file.iloc[:,14].values.tolist())
                            long4 = np.array(full_file.iloc[:,15].values.tolist())
                        x1 = np.sin(lat1_rad)*np.cos(long1_rad)
                        y1 = np.sin(lat1_rad)*np.sin(long1_rad)
                        x2 = np.sin(lat2_rad)*np.cos(long2_rad)
                        y2 = np.sin(lat2_rad)*np.sin(long2_rad)
                        x3 = np.sin(lat3_rad)*np.cos(long3_rad)
                        y3 = np.sin(lat3_rad)*np.sin(long3_rad)
                        x4 = np.sin(lat4_rad)*np.cos(long4_rad)
                        y4 = np.sin(lat4_rad)*np.sin(long4_rad)

                        filematch = identify_ion_to_mag_match(xpt,ypt,x1,x2,x3,x4,y1,y2,y3,y4)
                        rstepsize = 0.5

                    #;=====
                    #;End only inside of 15 Rj part
                    #;=====

                    #if insidept == 0:
                    #    filematch[0].size = 0
                    
                    #emptymatch = filematch[0].size == 0
                    if filematch != -1:
                        interpolation_output = mapping_interpolation_ion_to_mag(xpt,ypt,filematch,rstepsize,r1,ltime1,lat1,long1,r2,ltime2,lat2,long2,r3,ltime3,lat3,long3,r4,ltime4,lat4,long4,x1,x2,x3,x4,y1,y2,y3,y4)
                        r_interpol = interpolation_output[0]
                        ltime_interpol = interpolation_output[1]

                         #filematch works fine, determines whether final mapping point is inside the magnetosphere or not ; from joy et al 2002
                        is_inside = is_inside_joy_magnetopause(r_interpol,ltime_interpol*lt_to_rad)
                        if (is_inside == 0): #; beyond 150 RJ or expanded magnetosphere
                            r_interpol = -999.0
                            ltime_interpol = -999.0
                    else: #;if unable to find match for point in the files, test to see if point is inside of the 15 Rj contour
                        jloop = len(np.where(r1 == rmin_limit)) - 1
                        c = 0
                        jloop_start = jloop
                        nvert = jloop+1
                        y1_match = y1[np.where(r1 == rmin_limit)]
                        x1_match = x1[np.where(r1 == rmin_limit)]
                        #c = point_inside_polygon(xpt,ypt,x1_match,y1_match)
                        polygon_coords1 = np.zeros([len(x1_match),2])
                        for m in range(len(x1_match)):
                            polygon_coords1[m,0] = x1_match[m]
                            polygon_coords1[m,1] = y1_match[m]
                        polygon1 = Polygon(polygon_coords1)
                        # Check if the point is within the polygon
                        point = Point(xpt,ypt)
                        inside15 = polygon1.contains(point)
                        if inside15 == True:# ;maps to outside 150 RJ or Jovian magnetopause
                            r_interpol = -999.0
                            ltime_interpol = -999.0
                        else:# maps to inside of 15 RJ (or 2 RJ for fieldline tracing)
                            r_interpol = -998.0
                            ltime_interpol = -998.0

                    #end bad file match

                #;end of loop sslong div by 10
                else: #;if sslong div by 10 != 0, still mapping from northern ionosphere to magnetosphere
                    sslong1 = sslong
                    sslong = sslong_old

                    sslongmatch_min = np.arange(36)*10.0
                    sslongmatch_max = np.arange(1,37)*10.0
                    sslongmatch_file = np.arange(36)*10
                    sslongmatch_file[0] = 360

                    sslongmatch1 = np.array(np.where((sslong >= sslongmatch_min) & (sslong < sslongmatch_max)))

                    sslong = sslong1

                    if fieldline_tracing == True:
                        first_mapping = mapping_function_2025('ion_to_mag',sslongmatch_file[sslongmatch1[0]],latpt_input,longpt_input,model,True)
                        second_mapping = mapping_function_2025('ion_to_mag',(sslongmatch_file[sslongmatch1[0]]+10) % 360,latpt_input,longpt_input,model,True)
                    else:
                        first_mapping = mapping_function_2025('ion_to_mag',sslongmatch_file[sslongmatch1[0]],latpt_input,longpt_input,model,False)
                        second_mapping = mapping_function_2025('ion_to_mag',(sslongmatch_file[sslongmatch1[0]]+10) % 360,latpt_input,longpt_input,model,False)

                    sslong_first = float(sslongmatch_file[sslongmatch1[0]])
                    if (sslong_first == 0) | (sslong_first == 360):
                        sslong_first = 0.

                    if abs(first_mapping[0]) >= 200.0:
                        r_interpol = second_mapping[0]
                        ltime_interpol = second_mapping[1]
                    elif abs(second_mapping[0]) >= 200.0:
                        r_interpol = first_mapping[0]
                        ltime_interpol = first_mapping[1]
                    else:
                        first_x = first_mapping[0]*math.cos(first_mapping[1]*math.pi/12.0)
                        first_y = first_mapping[0]*math.sin(first_mapping[1]*math.pi/12.0)
                        second_x = second_mapping[0]*math.cos(second_mapping[1]*math.pi/12.0)
                        second_y = second_mapping[0]*math.sin(second_mapping[1]*math.pi/12.0)
                        x_interpol = first_x + (second_x-first_x)*((360.0-sslong)-sslong_first)/10.0
                        y_interpol = first_y + (second_y-first_y)*((360.0-sslong)-sslong_first)/10.0
                        r_interpol = math.sqrt(x_interpol*x_interpol + y_interpol*y_interpol)
                        ltime_interpol = (math.atan2(y_interpol,x_interpol)*12.0/math.pi + 48.) % 24.0

                #;end sslong div 10 != 0 loop for northern ionosphere to magnetosphere

            else: #;// begin ionosphere->magnetosphere southern loop

                latpt = 90.0 + latpt#; // convert to colatitude

                if sslong % 10.0 == 0.0:#; mapping southern ionosphere to magnetosphere
                    sslong1 = sslong;
                    sslong = sslong_old;
                    sslongfix = int(sslong)
                    sslong_str = str(sslongfix)

                    if (sslong != 0.0) & (fieldline_tracing == False):
                        filetext = directory_name + 'sslong'+sslong_str+'s'+model_filename_south+'.txt'
                    if (sslong == 0.0) & (fieldline_tracing == False):
                        filetext = directory_name + 'sslong360s'+model_filename_south+'.txt'
                    if (sslong != 0.0) & (fieldline_tracing == True):
                        filetext = directory_name + 'fieldline_tracing_all_r_'+model_filename_south+'_south_cml'+sslong_str+'.txt'
                    if (sslong == 0.0) & (fieldline_tracing == True):
                        filetext = directory_name + 'fieldline_tracing_all_r_'+model_filename_south+'_south_cml360.txt'

                    sslong = sslong1

                    full_file = pandas.read_table(filetext,delim_whitespace=True, header=None)
                    r1 = np.array(full_file.iloc[:,0].values.tolist())
                    ltime1 = np.array(full_file.iloc[:,1].values.tolist())
                    lat1 = np.array(full_file.iloc[:,2].values.tolist())
                    long1 = np.array(full_file.iloc[:,3].values.tolist())
                    r2 = np.array(full_file.iloc[:,4].values.tolist())
                    ltime2 = np.array(full_file.iloc[:,5].values.tolist())
                    lat2 = np.array(full_file.iloc[:,6].values.tolist())
                    long2 = np.array(full_file.iloc[:,7].values.tolist())
                    if fieldline_tracing == True:
                        r4 = np.array(full_file.iloc[:,8].values.tolist())
                        ltime4 = np.array(full_file.iloc[:,9].values.tolist())
                        lat4 = np.array(full_file.iloc[:,10].values.tolist())
                        long4 = np.array(full_file.iloc[:,11].values.tolist())
                        r3 = np.array(full_file.iloc[:,12].values.tolist())
                        ltime3 = np.array(full_file.iloc[:,13].values.tolist())
                        lat3 = np.array(full_file.iloc[:,14].values.tolist())
                        long3 = np.array(full_file.iloc[:,15].values.tolist())
                    else:
                        r3 = np.array(full_file.iloc[:,8].values.tolist())
                        ltime3 = np.array(full_file.iloc[:,9].values.tolist())
                        lat3 = np.array(full_file.iloc[:,10].values.tolist())
                        long3 = np.array(full_file.iloc[:,11].values.tolist())
                        r4 = np.array(full_file.iloc[:,12].values.tolist())
                        ltime4 = np.array(full_file.iloc[:,13].values.tolist())
                        lat4 = np.array(full_file.iloc[:,14].values.tolist())
                        long4 = np.array(full_file.iloc[:,15].values.tolist())
                    
                    lat1_rad = lat1*math.pi/180.0
                    lat2_rad = lat2*math.pi/180.0
                    lat3_rad = lat3*math.pi/180.0
                    lat4_rad = lat4*math.pi/180.0
                    long1_rad = long1*math.pi/180.0
                    long2_rad = long2*math.pi/180.0
                    long3_rad = long3*math.pi/180.0
                    long4_rad = long4*math.pi/180.0
                    x1 = np.sin(lat1_rad)*np.cos(long1_rad)
                    y1 = np.sin(lat1_rad)*np.sin(long1_rad)
                    x2 = np.sin(lat2_rad)*np.cos(long2_rad)
                    y2 = np.sin(lat2_rad)*np.sin(long2_rad)
                    x3 = np.sin(lat3_rad)*np.cos(long3_rad)
                    y3 = np.sin(lat3_rad)*np.sin(long3_rad)
                    x4 = np.sin(lat4_rad)*np.cos(long4_rad)
                    y4 = np.sin(lat4_rad)*np.sin(long4_rad)

                    filematch = identify_ion_to_mag_match(xpt,ypt,x1,x2,x3,x4,y1,y2,y3,y4)

                    rstepsize = 5.0
                    #;===========
                    #;For fieldline tracing only - test if point is inside of 15 Rj but outside of 2 Rj
                    #;===========
                    if (filematch == -1) & (fieldline_tracing == True):
                        if sslong != 0.0:
                            filetext = directory_name + 'fieldline_tracing_2to15rj_'+model_filename_south+'_south_cml'+sslong_str+'.txt'
                        else:
                            filetext = directory_name + 'fieldline_tracing_2to15rj_'+model_filename_south+'_south_cml360.txt'
                        restore,filetext
                        sslong = sslong1;

                        full_file = pandas.read_table(filetext,delim_whitespace=True, header=None)
                        r1 = np.array(full_file.iloc[:,0].values.tolist())
                        ltime1 = np.array(full_file.iloc[:,1].values.tolist())
                        lat1 = np.array(full_file.iloc[:,2].values.tolist())
                        long1 = np.array(full_file.iloc[:,3].values.tolist())
                        r2 = np.array(full_file.iloc[:,4].values.tolist())
                        ltime2 = np.array(full_file.iloc[:,5].values.tolist())
                        lat2 = np.array(full_file.iloc[:,6].values.tolist())
                        long2 = np.array(full_file.iloc[:,7].values.tolist())
                        if fieldline_tracing == True:
                            r4 = np.array(full_file.iloc[:,8].values.tolist())
                            ltime4 = np.array(full_file.iloc[:,9].values.tolist())
                            lat4 = np.array(full_file.iloc[:,10].values.tolist())
                            long4 = np.array(full_file.iloc[:,11].values.tolist())
                            r3 = np.array(full_file.iloc[:,12].values.tolist())
                            ltime3 = np.array(full_file.iloc[:,13].values.tolist())
                            lat3 = np.array(full_file.iloc[:,14].values.tolist())
                            long3 = np.array(full_file.iloc[:,15].values.tolist())
                        else:
                            r3 = np.array(full_file.iloc[:,8].values.tolist())
                            ltime3 = np.array(full_file.iloc[:,9].values.tolist())
                            lat3 = np.array(full_file.iloc[:,10].values.tolist())
                            long3 = np.array(full_file.iloc[:,11].values.tolist())
                            r4 = np.array(full_file.iloc[:,12].values.tolist())
                            ltime4 = np.array(full_file.iloc[:,13].values.tolist())
                            lat4 = np.array(full_file.iloc[:,14].values.tolist())
                            long4 = np.array(full_file.iloc[:,15].values.tolist())
                            
                        x1 = np.sin(lat1_rad)*np.cos(long1_rad)
                        y1 = np.sin(lat1_rad)*np.sin(long1_rad)
                        x2 = np.sin(lat2_rad)*np.cos(long2_rad)
                        y2 = np.sin(lat2_rad)*np.sin(long2_rad)
                        x3 = np.sin(lat3_rad)*np.cos(long3_rad)
                        y3 = np.sin(lat3_rad)*np.sin(long3_rad)
                        x4 = np.sin(lat4_rad)*np.cos(long4_rad)
                        y4 = np.sin(lat4_rad)*np.sin(long4_rad)

                        filematch = identify_ion_to_mag_match(xpt,ypt,x1,x2,x3,x4,y1,y2,y3,y4)
                        rstepsize = 0.5

                    #;=====
                    #;End only inside of 15 Rj part
                    #;=====


                    if filematch != -1:
                        interpolation_output = mapping_interpolation_ion_to_mag(xpt,ypt,filematch,rstepsize,r1,ltime1,lat1,long1,r2,ltime2,lat2,long2,r3,ltime3,lat3,long3,r4,ltime4,lat4,long4,x1,x2,x3,x4,y1,y2,y3,y4)
                        r_interpol = interpolation_output[0]
                        ltime_interpol = interpolation_output[1]
                        is_inside = is_inside_joy_magnetopause(r_interpol,ltime_interpol*lt_to_rad)#filematch works fine, determines whether final mapping point is inside the magnetosphere or not ; from joy et al 2002
                        if (is_inside == 0): #; beyond 150 RJ or expanded magnetosphere
                            r_interpol = -999.0
                            ltime_interpol = -999.0
                    else: #;if unable to find match for point in the files, test to see if point is inside of the 15 Rj contour
                        jloop = len(np.where(r1 == rmin_limit)) - 1
                        c = 0
                        jloop_start = jloop
                        nvert = jloop+1
                        y1_match = y1[np.where(r1 == rmin_limit)]
                        x1_match = x1[np.where(r1 == rmin_limit)]
                        #c = point_inside_polygon(xpt,ypt,x1_match,y1_match)
                        polygon_coords1 = np.zeros([len(x1_match),2])
                        for m in range(len(x1_match)):
                            polygon_coords1[m,0] = x1_match[m]
                            polygon_coords1[m,1] = y1_match[m]
                        polygon1 = Polygon(polygon_coords1)
                        # Check if the point is within the polygon
                        point = Point(xpt,ypt)
                        inside15 = polygon1.contains(point)
                        if inside15 == True:# ;maps to outside 150 RJ or Jovian magnetopause
                            r_interpol = -999.0
                            ltime_interpol = -999.0
                        else:# maps to inside of 15 RJ (or 2 RJ for fieldline tracing)
                            r_interpol = -998.0
                            ltime_interpol = -998.0

                    #;end of loop sslong div by 10
                else: #;if sslong div by 10 != 0 still southern ionosphere to magnetosphere
                    sslong1 = sslong;
                    sslong = sslong_old;

                    sslongmatch_min = np.arange(36)*10.0
                    sslongmatch_max = np.arange(1,37)*10.0
                    sslongmatch_file = np.arange(36)*10
                    sslongmatch_file[0] = 360

                    sslongmatch1 = np.array(np.where((sslong >= sslongmatch_min) & (sslong < sslongmatch_max)))
                    sslong = sslong1

                    if fieldline_tracing == True:
                        first_mapping = mapping_function_2025('ion_to_mag',sslongmatch_file[sslongmatch1[0]],latpt_input,longpt_input,model,True)
                        second_mapping = mapping_function_2025('ion_to_mag',(sslongmatch_file[sslongmatch1[0]]+10) % 360,latpt_input,longpt_input,model,True)
                    else:
                        first_mapping = mapping_function_2025('ion_to_mag',sslongmatch_file[sslongmatch1[0]],latpt_input,longpt_input,model,False)
                        second_mapping = mapping_function_2025('ion_to_mag',(sslongmatch_file[sslongmatch1[0]]+10) % 360,latpt_input,longpt_input,model,False)

                    sslong_first = float(sslongmatch_file[sslongmatch1[0]])
                    if (sslong_first == 0) | (sslong_first == 360):
                        sslong_first = 0.

                    if abs(first_mapping[0]) >= 200.0:
                        r_interpol = second_mapping[0]
                        ltime_interpol = second_mapping[1]
                    elif abs(second_mapping[0]) >= 200.0:
                        r_interpol = first_mapping[0]
                        ltime_interpol = first_mapping[1]
                    else:
                        first_x = first_mapping[0]*np.cos(first_mapping[1]*np.pi/12.0)
                        first_y = first_mapping[0]*np.sin(first_mapping[1]*np.pi/12.0)
                        second_x = second_mapping[0]*np.cos(second_mapping[1]*np.pi/12.0)
                        second_y = second_mapping[0]*np.sin(second_mapping[1]*np.pi/12.0)
                        x_interpol = first_x + (second_x-first_x)*((360.0-sslong)-sslong_first)/10.0
                        y_interpol = first_y + (second_y-first_y)*((360.0-sslong)-sslong_first)/10.0
                        r_interpol = np.sqrt(x_interpol*x_interpol + y_interpol*y_interpol)
                        ltime_interpol = (math.atan2(y_interpol,x_interpol)*12.0/np.pi + 48.) % 24.0            

                #;end sslong div 10 != 0 loop                 
            #endelse ;end southern hemisphere ion to mag mapping
        #end ion to mag mapping
        elif mag_to_ion_north == 1:
            xpt = rj*np.cos(loctime_rad)
            ypt = rj*np.sin(loctime_rad)

            if sslong % 10.0 == 0.0: #; magnetosphere to ionosphere north
                sslong1 = sslong
                sslong = sslong_old
                sslongfix = int(sslong)
                sslong_str = str(sslongfix)

                if sslong != 0.0:
                    filetext = directory_name + 'sslong'+sslong_str+model_filename_north+'.txt'
                if sslong == 0.0:
                    filetext = directory_name + 'sslong360'+model_filename_north+'.txt'
                if (sslong != 0.0) & (fieldline_tracing == True):
                    filetext = directory_name + 'fieldline_tracing_all_r_'+model_filename_north+'_north_cml'+sslong_str+'.txt'
                if (sslong == 0.0) & (fieldline_tracing == True):
                    filetext = directory_name + 'fieldline_tracing_all_r_'+model_filename_north+'_north_cml360.txt'
                if (sslong != 0.0) & (fieldline_tracing == True) & (rj < 15.0):
                    filetext = directory_name + 'fieldline_tracing_2to15rj_'+model_filename_north+'_north_cml'+sslong_str+'.txt'
                if (sslong == 0.0) & (fieldline_tracing == True) & (rj < 15.0):
                    filetext = directory_name + 'fieldline_tracing_2to15rj_'+model_filename_north+'_north_cml360.txt'

                sslong = sslong1
                full_file = pandas.read_table(filetext,delim_whitespace=True, header=None)
                r1 = np.array(full_file.iloc[:,0].values.tolist())
                ltime1 = np.array(full_file.iloc[:,1].values.tolist())
                lat1 = np.array(full_file.iloc[:,2].values.tolist())
                long1 = np.array(full_file.iloc[:,3].values.tolist())
                r2 = np.array(full_file.iloc[:,4].values.tolist())
                ltime2 = np.array(full_file.iloc[:,5].values.tolist())
                lat2 = np.array(full_file.iloc[:,6].values.tolist())
                long2 = np.array(full_file.iloc[:,7].values.tolist())
                if fieldline_tracing == True:
                    r4 = np.array(full_file.iloc[:,8].values.tolist())
                    ltime4 = np.array(full_file.iloc[:,9].values.tolist())
                    lat4 = np.array(full_file.iloc[:,10].values.tolist())
                    long4 = np.array(full_file.iloc[:,11].values.tolist())
                    r3 = np.array(full_file.iloc[:,12].values.tolist())
                    ltime3 = np.array(full_file.iloc[:,13].values.tolist())
                    lat3 = np.array(full_file.iloc[:,14].values.tolist())
                    long3 = np.array(full_file.iloc[:,15].values.tolist())
                else:
                    r3 = np.array(full_file.iloc[:,8].values.tolist())
                    ltime3 = np.array(full_file.iloc[:,9].values.tolist())
                    lat3 = np.array(full_file.iloc[:,10].values.tolist())
                    long3 = np.array(full_file.iloc[:,11].values.tolist())
                    r4 = np.array(full_file.iloc[:,12].values.tolist())
                    ltime4 = np.array(full_file.iloc[:,13].values.tolist())
                    lat4 = np.array(full_file.iloc[:,14].values.tolist())
                    long4 = np.array(full_file.iloc[:,15].values.tolist())


                mag_to_ion_mapping_output = mag_to_ion_mapping('north',xpt,ypt,r1,ltime1,lat1,long1,r2,ltime2,lat2,long2,r3,ltime3,lat3,long3,r4,ltime4,lat4,long4)
                lat_interpol = mag_to_ion_mapping_output[0]
                long_interpol = mag_to_ion_mapping_output[1]

            elif sslong % 10 != 0:  #still magnetosphere to ionosphere north
                sslong1 = sslong
                sslong = sslong_old

                sslongmatch_min = np.arange(36)*10.0
                sslongmatch_max = np.arange(1,37)*10.0
                sslongmatch_file = np.arange(36)*10
                sslongmatch_file[0] = 360

                sslongmatch1 = np.array(np.where((sslong >= sslongmatch_min) & (sslong < sslongmatch_max)))
                sslong = sslong1

                if fieldline_tracing: 
                    first_mapping = mapping_function_2025('mag_to_ion_north',sslongmatch_file[sslongmatch1[0]],rj_input,loctime_input,model,True)
                    second_mapping = mapping_function_2025('mag_to_ion_north',(sslongmatch_file[sslongmatch1[0]]+10) % 360,rj_input,loctime_input,model,True)
                else:
                    first_mapping = mapping_function_2025('mag_to_ion_north',sslongmatch_file[sslongmatch1[0]],rj_input,loctime_input,model,False)
                    second_mapping = mapping_function_2025('mag_to_ion_north',(sslongmatch_file[sslongmatch1[0]]+10) % 360,rj_input,loctime_input,model,False)

                sslong_first = (sslongmatch_file[sslongmatch1[0]])
                if sslong_first == 0 or sslong_first == 360:
                    sslong_first = 0.0

                if abs(first_mapping[0]) >= 200.0:
                    lat_interpol = second_mapping[0]
                    long_interpol = second_mapping[1]
                elif abs(second_mapping[0]) >= 200.0:
                    lat_interpol = first_mapping[0]
                    long_interpol = first_mapping[1]
                else:
                    first_x = abs(math.sin((90.0 -first_mapping[0])*deg_to_rad))*math.cos(first_mapping[1]*deg_to_rad)
                    first_y = abs(math.sin((90.0 -first_mapping[0])*deg_to_rad))*math.sin(first_mapping[1]*deg_to_rad)
                    second_x = abs(math.sin((90.0 -second_mapping[0])*deg_to_rad))*math.cos(second_mapping[1]*deg_to_rad)
                    second_y = abs(math.sin((90.0 -second_mapping[0])*deg_to_rad))*math.sin(second_mapping[1]*deg_to_rad)
                    x_interpol = first_x + (second_x-first_x)*((360.0-sslong)-sslong_first)/10.0
                    y_interpol = first_y + (second_y-first_y)*((360.0-sslong)-sslong_first)/10.0
                    lat_interpol = 90.0- math.asin(math.sqrt(x_interpol*x_interpol + y_interpol*y_interpol))*180.0/math.pi
                    long_interpol = (math.atan2(y_interpol,x_interpol)*180.0/math.pi + 720.0) % 360.0


            #endelse ; end ionosphere to magnetosphere north sslong != 10

        elif mag_to_ion_south == 1: #;end mag to ion north

            xpt = rj*np.cos(loctime_rad)
            ypt = rj*np.sin(loctime_rad)
            latpt = 90.0 - latpt #// convert to colatitude

            if sslong % 10.0 == 0.0:
                sslong1 = sslong
                sslong = sslong_old
                sslongfix = int(sslong)
                sslong_str = str(sslongfix)


                if sslong != 0.0:
                    filetext = directory_name + 'sslong'+sslong_str+'s'+model_filename_south+'.txt'
                if sslong == 0.0:
                    filetext = directory_name + 'sslong360s'+model_filename_south+'.txt'
                if (sslong != 0.0) & (fieldline_tracing == True):
                    filetext = directory_name + 'fieldline_tracing_all_r_'+model_filename_south+'_south_cml'+sslong_str+'.txt'
                if (sslong == 0.0) & (fieldline_tracing == True):
                    filetext = directory_name + 'fieldline_tracing_all_r_'+model_filename_south+'_south_cml360.txt'
                if (sslong != 0.0) & (fieldline_tracing == True) & (rj < 15.0):
                    filetext = directory_name + 'fieldline_tracing_2to15rj_'+model_filename_south+'_south_cml'+sslong_str+'.txt'
                if (sslong == 0.0) & (fieldline_tracing == True) & (rj < 15.0):
                    filetext = directory_name + 'fieldline_tracing_2to15rj_'+model_filename_south+'_south_cml360.txt'

                sslong = sslong1

                full_file = pandas.read_table(filetext,delim_whitespace=True, header=None)
                r1 = np.array(full_file.iloc[:,0].values.tolist())
                ltime1 = np.array(full_file.iloc[:,1].values.tolist())
                lat1 = np.array(full_file.iloc[:,2].values.tolist())
                long1 = np.array(full_file.iloc[:,3].values.tolist())
                r2 = np.array(full_file.iloc[:,4].values.tolist())
                ltime2 = np.array(full_file.iloc[:,5].values.tolist())
                lat2 = np.array(full_file.iloc[:,6].values.tolist())
                long2 = np.array(full_file.iloc[:,7].values.tolist())
                if fieldline_tracing == True:
                    r4 = np.array(full_file.iloc[:,8].values.tolist())
                    ltime4 = np.array(full_file.iloc[:,9].values.tolist())
                    lat4 = np.array(full_file.iloc[:,10].values.tolist())
                    long4 = np.array(full_file.iloc[:,11].values.tolist())
                    r3 = np.array(full_file.iloc[:,12].values.tolist())
                    ltime3 = np.array(full_file.iloc[:,13].values.tolist())
                    lat3 = np.array(full_file.iloc[:,14].values.tolist())
                    long3 = np.array(full_file.iloc[:,15].values.tolist())
                else:
                    r3 = np.array(full_file.iloc[:,8].values.tolist())
                    ltime3 = np.array(full_file.iloc[:,9].values.tolist())
                    lat3 = np.array(full_file.iloc[:,10].values.tolist())
                    long3 = np.array(full_file.iloc[:,11].values.tolist())
                    r4 = np.array(full_file.iloc[:,12].values.tolist())
                    ltime4 = np.array(full_file.iloc[:,13].values.tolist())
                    lat4 = np.array(full_file.iloc[:,14].values.tolist())
                    long4 = np.array(full_file.iloc[:,15].values.tolist())

                mag_to_ion_mapping_output = mag_to_ion_mapping('south',xpt,ypt,r1,ltime1,lat1,long1,r2,ltime2,lat2,long2,r3,ltime3,lat3,long3,r4,ltime4,lat4,long4)
                lat_interpol = mag_to_ion_mapping_output[0]
                long_interpol = mag_to_ion_mapping_output[1]
            #endif else begin if sslong % 10 != 0, mag to ion south
            elif sslong % 10 != 0:
                sslong1 = sslong
                sslong = sslong_old

                sslongmatch_min = np.arange(36)*10.0
                sslongmatch_max = np.arange(1,37)*10.0
                sslongmatch_file = np.arange(36)*10
                sslongmatch_file[0] = 360

                sslongmatch1 = np.array(np.where((sslong >= sslongmatch_min) & (sslong < sslongmatch_max)))
                sslong = sslong1

                if fieldline_tracing == True:
                    first_mapping = mapping_function_2025('mag_to_ion_south',sslongmatch_file[sslongmatch1[0]],rj_input,loctime_input,model,True)
                    second_mapping = mapping_function_2025('mag_to_ion_south',(sslongmatch_file[sslongmatch1[0]]+10) % 360,rj_input,loctime_input,model,True)
                else:
                    first_mapping = mapping_function_2025('mag_to_ion_south',sslongmatch_file[sslongmatch1[0]],rj_input,loctime_input,model,False)
                    second_mapping = mapping_function_2025('mag_to_ion_south',(sslongmatch_file[sslongmatch1[0]]+10) % 360,rj_input,loctime_input,model,False)


                sslong_first = float(sslongmatch_file[sslongmatch1[0]])
                if (sslong_first == 0) | (sslong_first == 360):
                    sslong_first = 0.0

                if abs(first_mapping[0]) >= 200.0:
                    lat_interpol = second_mapping[0]
                    long_interpol = second_mapping[1]
                elif abs(second_mapping[0]) >= 200.0:
                    lat_interpol = first_mapping[0]
                    long_interpol = first_mapping[1]
                else:
                    first_x = abs(math.sin((90.0 -first_mapping[0])*deg_to_rad))*math.cos(first_mapping[1]*deg_to_rad)
                    first_y = abs(math.sin((90.0 -first_mapping[0])*deg_to_rad))*math.sin(first_mapping[1]*deg_to_rad)
                    second_x = abs(math.sin((90.0 -second_mapping[0])*deg_to_rad))*math.cos(second_mapping[1]*deg_to_rad)
                    second_y = abs(math.sin((90.0 -second_mapping[0])*deg_to_rad))*math.sin(second_mapping[1]*deg_to_rad)
                    x_interpol = first_x + (second_x-first_x)*((360.0-sslong)-sslong_first)/10.0
                    y_interpol = first_y + (second_y-first_y)*((360.0-sslong)-sslong_first)/10.0
                    lat_interpol = math.asin(math.sqrt(x_interpol*x_interpol + y_interpol*y_interpol))*180.0/math.pi - 90.0
                    long_interpol = (math.atan2(y_interpol,x_interpol)*180.0/math.pi + 720.0) % 360.0


            #  end % 10 != 0
        #end mag to ion south mapping
    #end mapping
            
    sslong = sslong2

    if ((ion_to_mag == 1) & (abs(r_interpol) < 15.0) & (fieldline_tracing == False)) | ((ion_to_mag == 1) & (abs(r_interpol) < 2.0) & (fieldline_tracing == True)):
        r_interpol = 999.9
        ltime_interpol = 999.9
    if ion_to_mag == 1:
        return ([r_interpol,ltime_interpol])
    else:
        return ([lat_interpol,long_interpol])

