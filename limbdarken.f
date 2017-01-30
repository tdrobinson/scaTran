      subroutine limbdarken(iulimbd,ilimbd,logg,metal,
     -                      teff,wl,ldcoeff)
c
ccccccccccccccccccccccccccc   limbdarken   ccccccccccccccccccccccccccccc
cc                                                                    cc
cc    p u r p o s e :                                                 cc
cc      sets the coefficients for a wavelength-dependent limb         cc
cc      darkening model of the functional form:                       cc
cc        I(mu)/I0 = 1 - sum(i=1,4){C_i*[1 - mu**(i/2)]}              cc
cc      the coefficients (C_i = ldcoeff) are interpolated from the    cc
cc      extensive grids of Claret (2000).  entries are tabulated in   cc
cc      surface gravity, metallicity, and effective temperature, and  cc
cc      are given at 12 wavelength points.                            cc
cc                                                                    cc
cc        i n p u t :                                                 cc
cc          iulimbd  -- input/output unit for Claret (2000) table     cc
cc           ilimbd  -- model types: (1) Sun, (2) Phoenix,            cc
cc                                    (3) Atlas                       cc
cc             logg  -- for ilimbd of 2 or 3, the stellar log         cc
cc                       surface gravity (g in cm/s**2)               cc
cc            metal  -- for ilimbd of 2 or 3, the stellar metallicity cc
cc             teff  -- for ilimbd of 2 or 3, the stellar effective   cc
cc                        temperature (K)                             cc
cc               wl  -- the desired wavelength (in um) where the      cc
cc                      coefficients should be interpolated to        cc
cc                                                                    cc
cc        o u t p u t :                                               cc
cc          ldcoeff  -- vector of four (4) coefficients, interpolated cc
cc                      to wl, for use in the function given above    cc
cc                                                                    cc
cc    a u t h o r:                                                    cc
cc           Tyler D. Robinson (robinson.tyler.d@gmail.com)           cc
cc                                                                    cc
ccccccccccccccccccccccccccc   limbdarken   ccccccccccccccccccccccccccccc
c
      implicit none
c
      integer numa(19),nump(1)
      real teffa(61),teffp(40)
      real metala(19),metalp(1),mindiff
      real loggm(11)
c
      integer iulimbd,ilimbd
      integer i,iwl,nwl,map(12)
      real a(4,12),as(4,12),ai(12),ao,ldcoeff(4)
      real wlm(12),wl0(12),limbd_wl(12),limbd_wl0(12)
      real logg,metal,teff,wl
c
c*TDR wavelengths of model gridpoints, which are:
c           u v b y U B V R I J H K
c
      nwl = 12
      wlm = (/ 0.350, 0.411, 0.467, 0.547, 0.365,
     -         0.445, 0.551, 0.658, 0.806, 1.220,
     -         1.630, 2.190 /)
c
c*TDR mapping from wavelengths above to an order that
c     monotonically increases with wavelength
c
      map = (/ 1, 5, 2, 6, 3, 4, 7, 8, 9, 10, 11, 12 /)
c
c*TDR in order of increasing wavelength
c
      wl0 = (/ 0.350, 0.365, 0.411, 0.445, 0.467,
     -         0.547, 0.551, 0.658, 0.806, 1.220,
     -         1.630, 2.190 /)  
c
c*TDR coefficients for solar model
c
      as(1,1:12) = (/ 0.4927, 0.4614, 0.3899, 0.5258, 0.5928, 
     -                0.4767, 0.5311, 0.6071, 0.6785, 0.5239, 
     -                0.6705, 0.5613 /)
      as(2,1:12) = (/-0.5768,-0.2988, 0.2103,-0.0395,-0.8045,
     -               -0.1591,-0.0545,-0.2195,-0.4515, 0.1084,
     -                0.0845, 0.0325 /)
      as(3,1:12) = (/ 1.7919, 1.3589, 0.6405, 0.7261, 1.9579,
     -                1.0711, 0.7301, 0.7601, 0.8718,-0.0604,
     -               -0.4463,-0.2695 /)
      as(4,1:12) = (/-0.7834,-0.6212,-0.3800,-0.4058,-0.8251,
     -              -0.5154,-0.4053,-0.4058,-0.4292,-0.0241,
     -               0.2130, 0.1221 /)
c
c*TDR metallicity for atlas models
c
      metala = (/ -5.0, -4.5, -4.0, -3.5, -3.0,
     -            -2.5, -2.0, -1.5, -1.0, -0.5,
     -            -0.3, -0.2, -0.1,  0.0,  0.1,
     -             0.2,  0.3,  0.5,  1.0 /)
c
c*TDR metallicity for phoenix models
c
      metalp = (/ 0.0 /)
c
c*TDR number of atlas models for each metallicity
c
      numa = (/ 1544, 1568, 1568, 1588, 1612,
     -          1604, 1604, 1620, 1636, 1640,
     -          1640, 1640, 1640, 1904, 1564,
     -          1616, 1592, 1580, 1500 /)
c
c*TDR number of phoenix models for each metallicity
c
      nump = (/ 560 /)
c
c*TDR log gravity for models
c
      loggm = (/ 0.0, 0.5, 1.0, 1.5, 2.0, 2.5,
     -           3.0, 3.5, 4.0, 4.5, 5.0 /)
c
c*TDR temperature gridpoints for atlas models
c
      teffa = (/ 3500.,  3750.,  4000.,  4250.,  4500.,
     -           4750.,  5000.,  5250.,  5500.,  5750., 
     -           6000.,  6250.,  6500.,  6750.,  7000., 
     -           7250.,  7500.,  7750.,  8000.,  8250., 
     -           8500.,  8750.,  9000.,  9250.,  9500., 
     -           9750., 10000., 10500., 11000., 11500., 
     -          12000., 12500., 13000., 14000., 15000.,
     -          16000., 17000., 18000., 19000., 20000.,
     -          21000., 22000., 23000., 24000., 25000.,
     -          26000., 27000., 28000., 29000., 30000.,
     -          31000., 32000., 33000., 34000., 35000.,
     -          37500., 40000., 42500., 45000., 47500.,
     -          50000. /)
c
c*TDR temperature gridpoints for phoenix models
c
      teffp = (/ 2000.,  2200.,  2400.,  2600.,  2800.,
     -           3000.,  3200.,  3400.,  3600.,  3800.,
     -           4000.,  4200.,  4400.,  4600.,  4800.,
     -           5000.,  5200.,  5400.,  5600.,  5800.,
     -           6000.,  6200.,  6400.,  6600.,  6800.,
     -           7000.,  7200.,  7400.,  7600.,  7800.,
     -           8000.,  8200.,  8400.,  8600.,  8800.,
     -           9000.,  9200.,  9400.,  9600.,  9800./)
c
c*TDR   use solar model
c
      if( ilimbd .eq. 1 ) then
c
c*TDR   wavelength-dependent limb darkening, from Eqn. (7) 
c       of Claret (2000)
c
c        do 1055 iwl=1,nwl
c          limbd_wl(iwl) = 1.
c          do 1051 i=1,4
c            a(i,iwl) = as(i,iwl)
c            limbd_wl(iwl) = limbd_wl(iwl) - a(i,iwl)*(1 - mu**(i/2.))
c1051      continue
c1055    continue
c
        do 1055 iwl=1,nwl
          do 1051 i=1,4
            a(i,iwl) = as(i,map(iwl))
1051      continue
1055    continue
c
      else if( ilimbd .eq. 2 ) then
c
c*TDR   use Phoenix models
c
        write(*,*) 'WARNING: Phoenix limb darkening requested,'
        write(*,*) 'but not currently implemented.  Will assume'
        write(*,*) '*no* limb darkening.'
c
        do 2055 iwl=1,nwl
          do 2051 i=1,4
            a(i,iwl) = 0.
2051      continue
2055    continue
c        
      else
c
c*TDR   use Atlas models
c
        write(*,*) 'WARNING: Atlas limb darkening requested,'
        write(*,*) 'but not currently implemented.  Will assume'
        write(*,*) '*no* limb darkening.'
c
        do 2555 iwl=1,nwl
          do 2551 i=1,4
            a(i,iwl) = 0.
2551      continue
2555    continue
c
c*TDR   determine closest model metallicity to provide metallicity
c
c        mindiff = ABS(metal - metala(1))
c        imin    = 1
c        do 2021 im=1,19
c          if( ABS(metal - metala(im)) .lt. mindiff ) then
c            mindiff = ABS(metal - metala(im))
c            imin    = im
c          endif
c2021    continue
c
c*TDR   determine number of lines to skip
c
c        lskip = 0
c        do 2041 im=1,immin-1
c          lksip = lskip + numa(im)
c2041    continue
c
c*TDR   determine nearest two temperature gridpoints
c
c        mindiff = ABS(teff - teffa(1))
c        imin    = 1
c        do 2061 it=1,61
c          if( ABS(teff - teffa(it)) .lt. mindiff ) then
c            mindiff = ABS(teff - teffa(it))
c            imin    = im
c          endif
c2061    continue
c        imin1 = imin
c        if (imin1 .eq. 1) then
c          imin2 = imin1+1
c        else if (imin1 .eq. 61) then
c          imin2 = imin1-1
c        else
c          imin2 = imin1+1
c          if( ABS(teffa(imin1-1)-teff) .lt. 
c     -        ABS(teffa(imin1+1)-teff) ) imin2 = imin1-1
c        endif 
      endif 
c
c*TDR interpolate coefficients to given wavelength
c
      if( wl .lt. wl0(1) ) then
        do 3005 i=1,4
          ldcoeff(i) = a(i,1)
3005    continue
      else if( wl .gt. wl0(nwl) ) then
        do 3011 i=1,4
          ldcoeff(i) = a(i,12)
3011    continue
      else
        do 3025 i=1,4
          do 3021 iwl=1,nwl
            ai(iwl) = a(i,iwl)
3021      continue
c
          call xyinterp(wl0,ai,wl,ao,
     -                  nwl,1,nwl,1,1)
          ldcoeff(i) = ao
c
3025    continue
      endif
c
      return
      end