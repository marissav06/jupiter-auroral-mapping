# jupiter-auroral-mapping
Code for magnetosphere-ionosphere mapping at Jupiter

Jovian magnetosphere/ionosphere mapping program - README file
Written by Marissa Vogt, mvogt@psi.edu
Last updated July 2024

This mapping program allows a user to magnetically map a point between Jupiter's middle magnetosphere and the ionosphere. The program was originally intended to be used with the flux equivalence mapping model of Vogt et al. [2011, 2015] but also gives mappings obtained by tracing field lines from a global field model. 

Running the mapping program requires several files:
* mapping_function2024.pro: This is the function used to do the mapping.
* Two IDL routines (atan2, point_inside_polygon) called by mapping_function.pro
* A set of precalculated contour files (IDL .sav format) that can be downloaded at bit.ly/mappingfiles2024 (updated links will be posted on www.marissavogt.com/mapping and can be obtained by emailing mvogt@psi.edu). These files (sslong0_gam.sav, sslong10_jrm09.sav, etc.) are read by mapping_function2024.pro and should be placed in an appropriate directory - see line 233 of mapping_function2024.pro (replace ‘Documents/mapping_function_files/’ as appropriate)


See comments at the beginning of mapping_function2024.pro for detailed instructions on how to use the program.

Questions? Comments? Suggestions? Found an error? E-mail mvogt@psi.edu

Updates in 2024:
* Adds JRM33 as an option

Updates in 2020 version:
* Fieldline tracing option is now valid as close in as 2 Rj (flux mapping is valid beyond 15 Rj)
* I've recently made some minor changes to my CAN (Connerney et al. 1981) current sheet code, which affects all mapping options except those using the Khurana model. The expected changes should be very minor (a few RJ/fraction of an hour LT in mapping ionosphere->magnetosphere, or a fraction of a degree in lat/lon mapping magnetosphere->ionosphere). All mapping options except JRM09 use the CAN current sheet model with the VIP4 dipole values (tilt/longitude). The JRM09 flux mapping and field line tracing options use the CAN sheet with the JRM09 dipole values. 
* Speed improvements - program now loads in IDL .sav files instead of .txt files

Updates in June 2019 version:
* Now includes the a version of the Khurana model with JRM09 (instead of VIP4) as the internal field
* Fixed an error with local time mapping for the fieldline tracing option

Updates since 2017 version:
* Added the JRM09 field model (Connerney et al., 2018)
* Changed the way that the mapping function is called. 
* Changed the default mapping from GAM (north) and VIP4 (south) to JRM09 for both hemispheres.  
* Fixed an error that caused sharp discontinuities in mapping at certain subsolar longitudes. Thanks to Benjamin Palmaerts for finding this error. 
* The mapping is now done for Jupiter's oblate, flattened, surface rather than to a constant radial distance of 0.95 RJ.
* Other miscellaneous minor bug fixes.      


ACKNOWLEDGMENTS

Thanks to Bertrand Bonfond, Benjamin Palmaerts, and Andrew Steffl for helping to test early versions of this code and the online mapping tool. Very special thanks to Masafumi Imai for assistance implementing the JRM09 field model (2019).


REFERENCES

Vogt, M. F., M. G. Kivelson, K. K. Khurana, R. J. Walker, B. Bonfond, D. Grodent, and A. Radioti (2011), Improved mapping of Jupiter’s auroral features to magnetospheric sources, J. Geophys. Res., 116, A03220, doi:10.1029/2010JA016148, http://onlinelibrary.wiley.com/doi/10.1029/2010JA016148/abstract

Vogt, M. F., E. J. Bunce, M. G. Kivelson, K. K. Khurana, R. J. Walker, A. Radioti, B. Bonfond, and D. Grodent (2015), Magnetosphere-ionosphere mapping at Jupiter: Quantifying the effects of using different internal field models, J. Geophys. Res. Space Physics, 10.1002/2014JA020729, http://onlinelibrary.wiley.com/doi/10.1002/2014JA020729/abstract
