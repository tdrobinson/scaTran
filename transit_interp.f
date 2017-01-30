      subroutine transit_interp(ipi0,nlev,itrnst,ilimbd,irefract,
     -                          iulimbd,ncomp,icomp,volmix,nr0,p,t,alt,
     -                          grav,ratm,au,radius,radstar,logg,metal,
     -                          teff,bimpact,Ntheta,wn,copi,dtau_ex,
     -                          pathg,path,areag,area,do_mc,mu,trans_b,
     -                          mu_b,dtau_ab_b,dtrans_dtau_b,dmu_dtau_b,
     -                          zeff,tdepth)
c
cccccccccccccccccccccccccc  transit_interp  cccccccccccccccccccccccccccc
cc                                                                    cc
cc    p u r p o s e :                                                 cc
cc      Uses geometric path distribution, ray tracing, or scattering  cc
cc      Jacobians to compute a transit depth.                         cc
cc                                                                    cc
cc    i n p u t s :                                                   cc
cc      ipi0   - (0) non-scattering, (1) scattering                   cc
cc      nlev   - number of atmospheric levels                         cc
cc      itrnst - (1) geometric, (2) ray tracing, or (3) monte carlo   cc
cc      ilimbd - (0) no limb darkening, (1) limb darkening            cc
cc    irefract - (0) no refraction, (1) refraction                    cc
cc    iulimbd  - I/O unit for limb darkening parameters               cc
cc      ncomp  - number of background gases                           cc
cc      icomp  - background gas codes                                 cc
cc      volmix - volume mixing ratio of background gases              cc
cc       nr0   - index of refraction where path is computed           cc
cc     radstar - stellar radius (in solar radii)                      cc
cc     bimpact - distance from star center to planet center (in Rs)   cc
cc      radius - planetary radius (km)                                cc
cc      au     - planet-star separation (au)                          cc
cc      p      - atmospheric pressure profile (Pa)                    cc
cc      t      - atmospheric temperature profile (K)                  cc
cc      alt    - altitude profile (km)                                cc
cc      grav   - gravity profile (m/s/s)                              cc
cc      ratm   - atmospheric specific gas constant (J/K/kmole)        cc
cc      logg   - stellar log(gravity) (cgs)                           cc
cc      teff   - stellar effective temperature (K)                    cc
cc      metal  - stellar metalicity                                   cc
cc      wn     - wavenumber (cm**-1)                                  cc
cc      copi   - profile of layer single scattering co-albedos        cc
cc     dtau_ex - profile of layer extinction optical depths           cc
cc      areag  - area integration unit for geometric case (km**2)     cc
cc      area   - area integration unit for 2-D ray tracing (km**2)    cc
cc      pathg  - dimensionless geometric path distribution            cc
cc      path   - dimensionless path distribution for ray tracing      cc
cc      mu     - zenith angle on star for rays                        cc
cc      Ntheta - number of angular points in integration              cc
cc      do_mc  - profile of which altitudes have MC calculations      cc
cc     trans_b - altitude- and angle-dependent transmissions          cc
cc        mu_b - altitude and angle-dependent angles of stellar disk  cc
cc   dtau_ab_b - absorption optical depths at which transmissions and cc
cc               Jacobians are computed (i.e., most valid)            cc
cc dtrans_dtau_b - transmission Jacobian                              cc
cc  dmu_dtau_b - angle on stellar disk Jacobian                       cc
cc                                                                    cc
cc    o u t p u t s :                                                 cc
cc      zeff   - effective transit altitude (km)                      cc
cc    tdepth   - transit depth                                        cc
cc                                                                    cc
cc    a u t h o r:                                                    cc
cc           Tyler D. Robinson (robinson.tyler.d@gmail.com)           cc
cc                                                                    cc
cccccccccccccccccccccccccc  transit_interp  cccccccccccccccccccccccccccc
c
      implicit none
c
      include 'param.inc'
c
      integer nlev,nlyr,itrnst,Ntheta,ilimbd,irefract,iulimbd
      integer k,l,i,ith,do_mc(kp),ipi0,ncomp,icomp(6)
      real radstar,logg,metal,teff,bimpact,radius,p(kp),t(kp),alt(kp)
      real dtau_ex(kp),dtau(kp),grav(kp),ratm,rgas,a0,nr0,au
      real wl,limbd,mus,I_I0,I_I0s,ldcoeff(4),zbar(kp),volmix(6)
      real tdepth,tdzeff,zeff,copi(kp),dtau_ab
      double precision wn
      double precision path(kp,kp),pathg(kp,kp),areag(kp),area(kp,mxth)
      double precision pi,radsol,kmau,Rs,b,d1,d2,tau,trans,dtrans
      double precision Aatm,Arp,As,Asp,Az,Ap
      double precision areadsk(kp),tdisk(kp),ddtau_ab(kp)
      double precision theta,dtheta,r0s,mu0,dmu0
      double precision trans_b(kp,mxth),mu_b(kp,mxth)
      double precision dtrans_dtau_b(kp,kp,mxth),mu(kp,mxth)
      double precision dmu_dtau_b(kp,kp,mxth),dtau_ab_b(kp)
c
      pi  = ACOS(-1.d0)
c
c**** define avagadro's number (mks units: molecules/kmole)
c
      data a0/6.02297e26/
c
c**** define gas constant (mks units: J/K/kmole)
c
      data rgas/8.3144598e3/
c
      nlyr = nlev - 1
c
      wl = 1.e4/real(wn)
c
c*TDR solar radius and au in km
c
      radsol = 695700.
      kmau   = 1.49597870700e8
      Rs     = DBLE(radstar*radsol)
      b      = DBLE(bimpact*Rs)
c
c*TDR set limb darkening coefficients, if required
c
      if( ilimbd .ne. 0 ) then
        call limbdarken(iulimbd,ilimbd,logg,
     -                  metal,teff,wl,ldcoeff)
      endif
c
c*TDR update ray tracing path and geometry, if needed
c
      if( (itrnst .eq. 2 .and. irefract .eq. 1) .or. 
     -    (itrnst .eq. 3 .and. irefract .eq. 1 .and. ipi0 .eq. 0) ) then
c
        call transit_rayt(nlev,Ntheta,ncomp,icomp,volmix,nr0,
     -                    wn,bimpact,radstar,au,radius,p,t,
     -                    alt,path,mu)
c
      endif
c
c*TDR quantities needed for transit spectrum calculation for cases both
c     with and without atmospheric refraction and/or scattering
c
      do 1001 k=1,nlyr
c
c*TDR   layer mean altitude and thickness (km)
c
        zbar(k) = 0.5*(alt(k)+alt(k+1))
c
        if( itrnst .eq. 3 ) then
c
c*TDR     absorption optical depths
c
          dtau_ab = copi(k)*dtau_ex(k)
c
c*TDR     absorption optical depth change
c
          ddtau_ab(k) = dtau_ab - dtau_ab_b(k)
c
        endif
c
1001  continue
c
c*TDR for the geometric case, no efficiency is gained by interpolating,
c     so we simply use the geometric path distribution and optical 
c     depths (interpolated to wn) to set the transit spectrum.
c
      if( itrnst .eq. 1 ) then
c
        Aatm = 0.
        do 1505 k=1,nlyr
          tau = 0.
          do 1501 l=1,k
            tau = tau + dtau_ex(l)*pathg(k,l)
1501      continue
c
          if( tau .lt. 20. ) then
            Aatm = Aatm + (1.d0 - exp(-tau))*areag(k)
          else
            Aatm = Aatm + areag(k)
          endif
1505    continue
c
c*TDR otherwise we are doing ray tracing or monte carlo, which require 
c     integration over angle, and may require wavelength interpolation.
c
      else
c
c*TDR   initialize areas blocked by planet, and set angular integration
c       width
c
        Ap     = 0.0d0
        Asp    = 0.0d0
        Az     = 0.0d0
        dtheta = pi/REAL(Ntheta)
c
c*TDR   loop for integral over shells (i.e., tangent rays)
c
        do 3501 k=1,nlyr
c
c*TDR     area occulted by opaque disk, which we'll use to determine 
c         the effective transit height
c
          areadsk(k) = 0.d0
c
c*TDR     loop for integral over angle
c
          do 3001 ith=1,Ntheta
c
c*TDR       transmission is determined using different logic for monte 
c           carlo versus ray tracing.  begin with ray tracing, which is
c           more straightforward.  only need to do transmission
c           calculation for first angular element -- does not depend on 
c           theta
c
            if( (itrnst .eq. 2 .or. do_mc(k) .eq. 0) .or. 
     -          (itrnst .eq. 3 .and. ipi0 .eq. 0) ) then
c
              if( ith .eq. 1 ) then
c
c*TDR           initialize transmission
c
                trans = 1.0d0
c
c*TDR           transmission is zero if ray strikes the surface
c
                if( path(k,nlev) .gt. 0 ) trans = 0.d0
c
                do 2051 l=1,nlyr
c
                  trans = trans*EXP(-dtau_ex(l)*path(k,l))
c
2051            continue
c
              endif
c
c*TDR         set zenith angle on stellar disk to value for bin
c
              mu0 = mu(k,ith)
c
c*TDR       otherwise we must interpolate monte carlo results to 
c           this wavenumber (i.e., dtau_ab)
c
            else
c
c*TDR         initialize transmission and zenith angle to values for bin
c
              trans = trans_b(k,ith)
              mu0   = mu_b(k,ith)
c
c*TDR         loop over layers, adding first-order Jacobian contributions
c
              dtrans = 0.0d0
              dmu0   = 0.0d0
              do 2501 l=1,nlyr
                dtrans = dtrans + dtrans_dtau_b(k,l,ith)*ddtau_ab(l)
                dmu0   = dmu0   + dmu_dtau_b(k,l,ith)*ddtau_ab(l)
2501          continue
c
              trans = trans + dtrans
              mu0   = mu0 + dmu0
c
c*TDR         check and correct bounds on transmission and zenith
c
              if( trans .gt. 1.0d0 ) trans = 1.0d0
              if( trans .lt. 0.0d0 ) trans = 0.0d0
              if( mu0   .gt. 1.0d0 ) mu0   = 2.0d0 - mu0
c
            endif
c
c*TDR       determine where (in radius) on the stellar disk 
c           this r, theta sit.
c
            theta = REAL(ith)*dtheta - dtheta/2
            r0s   = (((radius + zbar(k))*sin(theta))**2 + 
     -              (b - (radius + zbar(k))*cos(theta))**2)
c
c*TDR       determine the limb darkening at this radius.  honestly, 
c           this is a real pain in the a**.  limb darkening is already
c           considered in standard transit observation reduction 
c           pipelines (see, e.g., the Mandel and Agol 2002 paper).  so, 
c           we're just concerned with deviations from that theory.  this
c           only occurs in the non-geometric limit, where rays map to 
c           different locations on the stellar disk than the geometric 
c           case would have.  either the area element is on the off-disk 
c           part of the planet, which could map back to the disk and has 
c           the limb darkening for that location on the disk.  or, the 
c           element is on the disk, and the ray maps to something other 
c           than the geometric, in which case we adjust the brightness of 
c           the area element by the ratio of the bent limb darkening to 
c           the straight line limb darkening.
c
            limbd = 1.
            if( mu0 .lt. 0 ) limbd = 0.
c
            if( ilimbd .ne. 0 ) then
              if( mu0 .ge. 0 ) then
                I_I0 = 1.0
                do i=1,4
                  I_I0 = I_I0 - 
     -                   ldcoeff(i)*(1. - REAL(mu0)**(i/2.))
                enddo
              else
                I_I0 = 0.0
              endif
c
c*TDR         if this location is on the stellar disk, adjust the
c             limb darkening
c
              if( r0s .le. Rs**2 ) then
                I_I0s = 1.0
                mus   = (1.0d0 - r0s/Rs**2)**0.5
                do i=1,4
                  I_I0s = I_I0s - ldcoeff(i)*(1. - mus**(i/2.))
                enddo
                limbd = I_I0/I_I0s
              else
                limbd = I_I0
              endif
c
            endif
c
c*TDR       determine planetary flux.  special case for monte 
c           carlo, as transmission must also account for rays 
c           that do not strike the stellar disk
c
            Ap = Ap + trans*limbd*area(k,ith)
c
c*TDR       stellar flux only contributed for r, theta locations 
c           that are on the stellar disk.  also, keep track of the
c           "positive" portion of the planetary flux, which comes 
c           from where the planet disk overlaps the stellar disk.  we 
c           use this to determine an effective transit altitude.
c
            if( r0s .le. Rs**2 ) then
              Az    = Az  + trans*limbd*area(k,ith)
              Asp   = Asp + area(k,ith)
              areadsk(k) = areadsk(k) + area(k,ith)
            endif
c
2751      continue
c
3001    continue
c
3501   continue
c
      endif
c
c*TDR also consider the area blocked by the "solid" (i.e., opaque) 
c     part of the planet, where r < Rp.  the first 'if' checks to 
c     see if the 'solid' part of the planet even overlaps the 
c     stellar disk.
c
      if( b .gt. Rs + radius ) then
        Arp = 0.d0
      else
c
c*TDR   remaining logic is just the geometry of overlapping circles.
c
        Arp = pi*radius**2
        if( b .gt. Rs - radius ) then
          d1 = (b**2 - radius**2 + Rs**2)/2/b
          d2 = b - d1
          Arp = Rs**2*ACOS(d1/Rs) - d1*(Rs**2 - d1**2)**0.5 + 
     -           radius**2*ACOS(d2/radius) - d2*(radius**2 - d2**2)**0.5
        endif
c
      endif
c
c*TDR   find total brightness (area) of stellar disk
c
      As = pi*Rs**2
c
c*TDR determine transit depth in the "simple" case, 
c     i.e., geometric
c   
      if( itrnst .eq. 1 ) then
        tdepth = REAL((Arp +  Aatm)/As)
        tdzeff = tdepth
      else
c
c*TDR   add area to stellar brightness for region occulted by
c       the "solid" part of the disk
c
        Asp = Asp + Arp
c
        tdepth = REAL(Asp - Ap)/As
        tdzeff = REAL(Asp - Az)/As
c
      endif
c
c*TDR effective transit height, only analytic if the whole area of 
c     pi(radius + zeff)**2 is on the stellar disk.  otherwise 
c     requires interpolation.
c
      if( tdzeff .le. 0 ) then
        zeff = 0.
      else
        zeff = Rs*tdzeff**0.5 - radius
c
        if( zeff .lt. 0 .and. itrnst .eq. 1 ) zeff = 0.
c
        if( b+radius+zeff .gt. Rs .or. 
     -      (zeff .lt. 0 .and. itrnst .ne. 1) ) then
          if( itrnst .eq. 1 ) then
            tdisk(nlyr) = (Arp + areag(nlyr))/As
          else
            tdisk(nlyr) = (Arp + areadsk(nlyr))/As
          endif
          do 6051 k=nlyr-1,1,-1
            if( itrnst .eq. 1 ) then
              tdisk(k) = tdisk(k+1) + areag(k)/As
            else
              tdisk(k) = tdisk(k+1) + areadsk(k)/As
            endif
6051      continue
          call xyinterp(REAL(tdisk),zbar,tdzeff,zeff,kp,1,nlyr,1,1)
        endif
      endif
c
      return
      end