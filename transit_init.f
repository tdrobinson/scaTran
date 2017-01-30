      subroutine transit_init(nlev,itrnst,irefract,iupath,ncomp,icomp,
     -                        volmix,radstar,bimpact,radius,au,p,t,alt,
     -                        Ntheta,nr0,g0,dtau_sc0,areag,area,
     -                        pathg,path,mu)
c
cccccccccccccccccccccccccc   transit_init   cccccccccccccccccccccccccccc
cc                                                                    cc
cc    p u r p o s e :                                                 cc
cc      Initializes any gray aspects of the transit calculation,      cc
cc      including the geometric path distribution, the ray tracing    cc
cc      path distribution, the integration areas, and aspects of the  cc
cc      ray tracing geometry.                                         cc          
cc                                                                    cc
cc    i n p u t s :                                                   cc
cc      nlev   - number of atmospheric levels                         cc
cc      itrnst - (1) geometric, (2) ray tracing, or (3) monte carlo   cc
cc    irefract - (0) no refraction, (1) refraction                    cc
cc      iupath - I/O scratch unit numbers for storing path distr.     cc
cc      ncomp  - number of background gases                           cc
cc      icomp  - background gas codes                                 cc
cc      volmix - volume mixing ratio of background gases              cc
cc     radstar - stellar radius (in solar radii)                      cc
cc     bimpact - distance from star center to planet center (in Rs)   cc
cc      radius - planetary radius (km)                                cc
cc      au     - planet-star separation (au)                          cc
cc      p      - atmospheric pressure profile (Pa)                    cc
cc      t      - atmospheric temperature profile (K)                  cc
cc      alt    - altitude profile (km)                                cc
cc                                                                    cc
cc    o u t p u t s :                                                 cc
cc      Ntheta - number of angular points in integration              cc
cc      nr0    - index of refraction at wn                            cc
cc      g0     - asymmetry parameter normalized to -1                 cc
cc    dtau_sc0 - scattering optical depth normalized to -1            cc
cc      areag  - area integration unit for geometric case (km**2)     cc
cc      area   - area integration unit for 2-D ray tracing (km**2)    cc
cc      pathg  - dimensionless geometric path distribution            cc
cc      path   - dimensionless path distribution for ray tracing      cc
cc      mu     - initialized bin zenith angle on star for rays        cc
cc                                                                    cc
cc    a u t h o r:                                                    cc
cc           Tyler D. Robinson (robinson.tyler.d@gmail.com)           cc
cc                                                                    cc
cccccccccccccccccccccccccc   transit_init   cccccccccccccccccccccccccccc
c
      implicit none
c
      include 'param.inc'
c
      integer nlev,nlyr,k,l,ith,ioerr,iupath(kp),itrnst,Ntheta
      integer irefract,ncomp,icomp(6)
      real p(kp),t(kp),alt(kp),volmix(6),zbar(kp),dz(kp),au
      real nr0,g0(kp),dtau_sc0(kp),radius,bimpact,radstar
      double precision pi,radsol,kmau,Rs,b,d
      double precision area(kp,mxth),areag(kp),dtheta
      double precision rangle(kp),r0(kp)
      double precision path0(kp),path(kp,kp),pathg(kp,kp),rangle0,rad0
      double precision mu(kp,mxth)
c
      pi  = ACOS(-1.d0)
c
      nlyr = nlev - 1
c
c*TDR solar radius and au in km
c
      radsol = 695700.
      kmau   = 1.49597870700e8
      Rs     = DBLE(radstar*radsol)
      b      = DBLE(bimpact*Rs)
c
c*TDR resolution in radius and angle for integrating the planetary disk
c     over the stellar disk when considering limb darkening
c
      if( (b-Rs) .ge. 0 ) then
        Ntheta = MIN(180+INT(mxth*(b-Rs)/radius),180)
      else
        Ntheta = MAX(180+INT(90*(b-Rs)/radius),30)
      endif
      if( Ntheta .gt. mxth ) Ntheta = mxth
c
c*TDR for perfectly symmetric case (b ~ 0), only need one angle
c
      if( b/Rs .lt. 1.e-3 ) Ntheta = 1
c
c*TDR   set geometric path distribution and area matrix
c
      do 1155 k=1,nlyr
c
        zbar(k) = 0.5*(alt(k)+alt(k+1))
        dz(k)   = alt(k)-alt(k+1)
c
        call path_geom(nlev,radius,alt,zbar(k),path0)
c
        do 1101 l=1,nlyr
          pathg(k,l) = path0(l)
1101    continue
        pathg(k,nlev) = 0.d0
c
        dtheta = pi
        if( b .gt. Rs - (radius + zbar(k)) ) then
          d = (b**2 - (radius + zbar(k))**2 + Rs**2)/2/b
          dtheta = ACOS((b-d)/(radius + zbar(k)))
        endif
        if( b .gt. Rs + (radius + zbar(k)) ) then
          dtheta = 0.
        endif
c
        areag(k) = 2*dtheta*(radius+zbar(k))*dz(k)
c
        dtheta = pi/REAL(Ntheta)
        do 1125 ith=1,Ntheta
          area(k,ith) = 2*dtheta*(radius + zbar(k))*dz(k)
1125    continue
c
1155  continue
c
c*TDR   if doing refraction, initialize refractive index 
c
      if( irefract .eq. 1 ) nr0 = -1.0
c
c*TDR   initialize ray tracing path to simple geometric limit
c
      do 1755 k=1,nlyr
        rangle(k) = 0.d0
        r0(k) = radius+zbar(k)
        do 1701 l=1,nlyr
          path(k,l) = pathg(k,l)
1701    continue
1755  continue
c
c*TDR     set zenith angle on star where ray traced paths intersect.
c
      call set_zenith(nlev,Ntheta,radius,Rs,
     -                b,au,alt,rangle,r0,mu)
c
c*TDR   if doing monte carlo, open scratch files and store asymmetry 
c       parameter and scattering optical depth.
c
      if( itrnst .eq. 3 ) then
        do 2001 k=1,nlyr
c
          open(UNIT=iupath(k), STATUS='scratch', 
     -         FORM='unformatted', IOSTAT=ioerr)
c
          g0(k) = -1.
          dtau_sc0(k) = -1.
c
2001    continue
      endif
c
      return
      end
