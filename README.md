# jupiter-auroral-mapping
Code for magnetosphere-ionosphere mapping at Jupiter

Jovian magnetosphere/ionosphere mapping program - README file

Written by Marissa Vogt, mvogt@psi.edu

Last updated September 2025

This mapping program allows a user to magnetically map a point between Jupiter's middle magnetosphere and the ionosphere. The program was originally intended to be used with the flux equivalence mapping model of Vogt et al. [2011, 2015] but also gives mappings obtained by tracing field lines from a global field model.

The program is available in both IDL and Python as of September 2025.

Running the IDL mapping program requires:
* mapping_function2025.pro: This is the function used to do the mapping.
* Two IDL routines (atan2.pro, point_inside_polygon.pro) called by mapping_function.pro
* A set of precalculated contour files (IDL .sav format) that can be downloaded at bit.ly/mappingfiles2025

Running the IDL mapping program requires:
* mapping_function2025.py: This is the function used to do the mapping.
* A set of precalculated contour files (.txt format) that can be downloaded at bit.ly/mappingfiles2025

All files can be downloaded from bit.ly/mappingfiles2025 (updated links will be posted on marissavogt.com/mapping). These files are read by the mapping_function2025 code and should be placed in an appropriate directory - replace ‘Documents/mapping_function_files/’ in the mapping function code as appropriate)

Questions? Comments? Suggestions? Found an error? E-mail mvogt@psi.edu

Updates in 2025:
* JRM33 is the default internal field option for both flux equivalence and fieldline tracing
* Code is now available in Python as well as IDL
* The IDL code has been significantly restructured and also runs about 40% faster thanks to several speed up suggestions by Rob Wilson   

Please see the comments at the beginning of the mapping function code for detailed instructions on how to use the program as well as a summary of previous revisions.
 

ACKNOWLEDGMENTS

Thanks to Bertrand Bonfond, Benjamin Palmaerts, and Andrew Steffl for helping to test early versions of this code and the online mapping tool. Very special thanks to Masafumi Imai for assistance implementing the JRM09 field model (2019). Thanks to Rob Wilson for help implementing several speed ups to the IDL version.


REFERENCES

Vogt, M. F., M. G. Kivelson, K. K. Khurana, R. J. Walker, B. Bonfond, D. Grodent, and A. Radioti (2011), Improved mapping of Jupiter’s auroral features to magnetospheric sources, J. Geophys. Res., 116, A03220, doi:10.1029/2010JA016148, http://onlinelibrary.wiley.com/doi/10.1029/2010JA016148/abstract

Vogt, M. F., E. J. Bunce, M. G. Kivelson, K. K. Khurana, R. J. Walker, A. Radioti, B. Bonfond, and D. Grodent (2015), Magnetosphere-ionosphere mapping at Jupiter: Quantifying the effects of using different internal field models, J. Geophys. Res. Space Physics, 10.1002/2014JA020729, http://onlinelibrary.wiley.com/doi/10.1002/2014JA020729/abstract
