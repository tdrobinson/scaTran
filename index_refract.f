      subroutine index_refract(ncomp,icomp,volmix,wl,p,t,nr)
c
ccccccccccccccccccccccccccc  index_refract  cccccccccccccccccccccccccccc
cc                                                                    cc
cc    p u r p o s e :                                                 cc
cc      computes the pressure- and temperature-dependent index of     cc
cc      refraction for the atmosphere, given information about the    cc
cc      main atmospheric constituents.  note that the assumed input   cc
cc      units are microns for wavelength, pascal for pressure, and    cc
cc      kelvin for temperature. index of refraction is determined     cc
cc      from analytic fits for individual gases, with data sources    cc
cc      given below.                                                  cc
cc                                                                    cc
cc    i n p u t s:                                                    cc
cc         ncomp --    number of primary atmospheric constituents     cc
cc         icomp --    index of primary atmospheric constituents:     cc
cc                     (1) air (2) co2 (3) n2 (4) o2 (5) h2 (6) he    cc
cc        volmix --    volume mixing ratio of major constituents      cc
cc            wl --    wavelength (in um) where index of refraction   cc
cc                     is desired                                     cc
cc             p --    pressure (in pa) where index is desired        cc
cc             t --    temperature (in K) at which index is desired   cc
cc                                                                    cc
cc    o u t p u t s:                                                  cc
cc            nr --    index of refraction at wl, p, and t            cc
cc                                                                    cc
cc                                                                    cc
cc    a u t h o r:                                                    cc
cc           Tyler D. Robinson (robinson.tyler.d@gmail.com)           cc
cc                                                                    cc
ccccccccccccccccccccccccccc  index_refract  cccccccccccccccccccccccccccc
c
      implicit none
c
      include 'param.inc'
c
      integer ncomp,icomp(6),i
      real p,t,wl,volmix(6)
      real p0,t0,nu,nu_i,sigma,nr
c
c*TDR inverse of wavelength (um)
c
      sigma = 1./wl
c
c*TDR determine refractivity of major constituents, 
c     and produce a volume mixing ratio average.  note 
c     that refractivities are given at a standard 
c     pressure (p0) and temperature (t0), and so are 
c     adjusted assuming the ideal gas law.
c
      nu = 0.
      do 1001 i=1,ncomp
c
c****   air [Birch (1994). Metrologia, 31:315]
c
        if (icomp(i) .eq. 1) then
          p0   = 101325.
          t0   = 298.
c
c*TDR     constant slope below 0.2 um
c
          if (wl .lt. 0.2) then
            nu_i = 0.000324092 + 0.000661155*(wl - 0.2)
          else
c
c*TDR       Birch (1994) expression
c
            nu_i = 1.e-8*( 8342.54 
     -              + 2406147./(130.-sigma**2)
     -              + 15998./(38.9-sigma**2) )
          endif
c
          nu_i = nu_i*(p/p0)/(t/t0)
        endif
c
c****   carbon dioxide [Bideau-Mehu et al. (1973). 
c       Opt. Commun., 9:432]
c
        if (icomp(i) .eq. 2) then
          p0   = 101325.
          t0   = 273.
c
c*TDR     constant slope below 0.2 um
c
          if (wl .lt. 0.2) then
            nu_i = 0.000526177 + 0.000949764*(wl - 0.2)
          else
c
c*TDR     Bideau-Mehu et al. (1973) expression, and 
c         assume constant above ~2 um
c
            if (wl .gt. 2.) sigma = 1./2.
            nu_i = 6.99100e-2/(166.175-sigma**2) 
     -               + 1.44720e-3/(79.609-sigma**2)
     -               + 6.42941e-5/(56.3064-sigma**2) 
     -               + 5.21306e-5/(46.0196-sigma**2)
     -               + 1.46847e-6/(0.0584738-sigma**2)
            if (wl .gt. 2.) sigma = 1./wl
          endif
c
          nu_i = nu_i*(p/p0)/(t/t0)
        endif
c
c****   nitrogen [Peck and Khanna (1963). J. Opt. Soc.
c       Am. 56:1059; 
c
        if (icomp(i) .eq. 3) then
          p0   = 101325.
          t0   = 273.
c
c*TDR     constant slope below 0.2 um
c
          if (wl .lt. 0.2) then
            nu_i = 0.000341086 + 0.000522109*(wl - 0.2)
          else
c
c*TDR     Peck and Khanna expression
c
            nu_i = 6.8552e-5 + 3.243157e-2/(144.-sigma**2)
          endif
c
          nu_i = nu_i*(p/p0)/(t/t0)
        endif
c
c****   oxygen [Zhang and Wang (2008). Appl. Opt. 47:3143]
c
c
        if (icomp(i) .eq. 4) then
          p0   = 101325.
          t0   = 293.
c
c*TDR     constant slope below 0.2 um
c
          if (wl .lt. 0.2) then
            nu_i = 0.000310787 + 0.000849213*(wl - 0.2)
          else
c
c*TDR     Zhang and Wang expression
c
            nu_i = 1.181494e-4 + 9.708931e-3/(75.4-sigma**2)
          endif
c
          nu_i = nu_i*(p/p0)/(t/t0)
        endif
c
c****   hydrogen [Leonard (1974).  Atomic Data and 
c       Nuclear Data Tables, 14:21]
c
c
        if (icomp(i) .eq. 5) then
          p0   = 101325.
          t0   = 273.
c
c*TDR     constant slope below 0.2 um
c
          if (wl .lt. 0.2) then
            nu_i = 0.000159383 + 0.000218933*(wl - 0.2)
          else
c
c*TDR       fitted expression
c
            nu_i = 1.49390e-05 + 0.0137745/(113.4-sigma**2)
          endif
c
          nu_i = nu_i*(p/p0)/(t/t0)
        endif
c
c****   helium [Mansfield and Peck (1969). J. Opt. 
c       Soc. Am. 59:199]
c
c
        if (icomp(i) .eq. 6) then
          p0   = 101000.
          t0   = 273.
c
c*TDR     constant slope below 0.2 um
c
          if (wl .lt. 0.2) then
            nu_i = 3.68462e-05 + 2.13404e-05*(wl - 0.2)
          else
c
c*TDR       Mansfield and Peck expression
c
            nu_i = 0.01470091/(423.98-sigma**2)
          endif
c
          nu_i = nu_i*(p/p0)/(t/t0)
        endif
c
c*TDR   volume mixing ratio-weighted average refractivity
c
        nu = nu + volmix(i)*nu_i
1001  continue
c
      nr = 1. + nu
c
      return
      end
