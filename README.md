# scaTran
Exoplanet transit model with Monte Carlo light scattering.

Transit model with light scattering from Robinson, T.D., 2017. ''A Theory of Exoplanet
Transits with Light Scattering.'' ApJ, subm.

Please cite this paper if you use any aspects of this code.  All software authored by 
T. D. Robinson, except for xyinterp.f (which is authored by D. Crisp).  Questions, 
concerns, and requested updates can be sent to robinson.tyler.d@gmail.com.

Note that this version assumes a Henyey-Greenstein scattering phase function, and only 
has limb darkening parameters for the Sun.  Time permitting, generalizations of these 
treatments will be added in the future.

Compiling the Makefile, which uses ifort, generates an example executable called scatran.

Current setup is for a hot Jupiter-like exoplanet, with a strongly forward scattering 
(g=0.95) cloud at 1e-4 bar.  The cloud has slant scattering optical depth equal to 10, 
and is distributed about one pressure scale height.  This setup mimics the transit 
spectrum at 1.15 um in Figure 14 of Robinson (2017).
