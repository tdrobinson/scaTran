      subroutine transit_rayt(nlev,Ntheta,ncomp,icomp,volmix,nr0,wn,
     -                        bimpact,radstar,au,radius,p,t,alt,path,mu)
c
cccccccccccccccccccccccccc   transit_rayt   cccccccccccccccccccccccccccc
cc                                                                    cc
cc    p u r p o s e :                                                 cc
cc      Initializes any gray aspects of the transit calculation,      cc
cc      including the geometric path distribution, the ray tracing    cc
cc      path distribution, the integration areas, and aspects of the  cc
cc      ray tracing geometry.                                         cc          
cc                                                                    cc
cc    i n p u t s :                                                   cc
cc      nlev   - number of atmospheric levels                         cc
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
cc      wn     - wavenumber to initialize at (cm**-1)                 cc
cc      Ntheta - number of angular points in integration              cc
cc      nr0    - index of refraction at last ray tracing              cc
cc                                                                    cc
cc    o u t p u t s :                                                 cc
cc      path   - dimensionless path distribution for ray tracing      cc
cc      mu     - initialized bin zenith angle on star for rays        cc
cc                                                                    cc
cc    a u t h o r:                                                    cc
cc           Tyler D. Robinson (robinson.tyler.d@gmail.com)           cc
cc                                                                    cc
cccccccccccccccccccccccccc   transit_rayt   cccccccccccccccccccccccccccc
c
      implicit none
c
      include 'param.inc'
c
      integer nlev,nlyr,k,l,Ntheta,ncomp,icomp(6),kmin,init
      real radius,nr0,nr(kp),wl,bimpact,tol,radstar
      real au,p(kp),t(kp),alt(kp),volmix(6),zbar(kp)
      double precision r0(kp),rangle(kp),rangle0,rad0,path0(kp)
      double precision hmin0,hmin
      double precision wn,mu(kp,mxth),path(kp,kp),dlnnrdz(kp)
      double precision b,pi,radsol,kmau,Rs
c
      pi  = ACOS(-1.d0)
c
      tol = 0.01
c
      nlyr = nlev - 1
c
      wl = 1.e4/REAL(wn)
c
      init = 1
c
c*TDR solar radius and au in km
c
      radsol = 695700.
      kmau   = 1.49597870700e8
      Rs     = DBLE(radstar*radsol)
      b      = DBLE(bimpact*Rs)
c
c*TDR check if index of refraction has changed by enough to warrant 
c     re-doing the ray tracing calculation
c
      call index_refract(ncomp,icomp,volmix,wl,
     -                   p(nlev),t(nlev),nr(nlev))
c
      if( ABS(nr(nlev)-nr0)/ABS(nr0-1.0) .gt. tol 
     -    .or. nr0 .le. 0 ) then
c
        if( nr0 .lt. 0 ) init = 0
c
        nr0 = nr(nlev)
c
        hmin0 = alt(1)
        kmin  = 1
c
c*TDR set derivative of index of refraction with respect to 
c     height. normalized to zero if irefract not equal to one.
c
        call index_refract(ncomp,icomp,volmix,wl,
     -                     p(1),t(1),nr(1))
c
        do 1501 k=1,nlyr
c
          call index_refract(ncomp,icomp,volmix,wl,
     -                       p(k+1),t(k+1),nr(k+1))
c
          dlnnrdz(k) = 
     -        dble((alog(nr(k))-alog(nr(k+1)))/(alt(k)-alt(k+1)))
c
1501    continue
c
c*TDR   with the derivative of the index of refraction, do the ray
c       tracing calculation
c
        do 1905 k=1,nlyr
c
          zbar(k) = 0.5*(alt(k)+alt(k+1))
c
          call path_rayt(nlev,radius,dlnnrdz,alt,
     -                   hmin,zbar(k),rangle0,rad0,path0)
c
          rangle(k) = rangle0
          r0(k)     = rad0
c
c*TDR     minimum altitude probed
c
          if( hmin .lt. hmin0 .and. hmin .gt. 0 ) then
            hmin0 = hmin
            kmin  = k
          endif
c
          do 1901 l=1,nlyr
            path(k,l) = path0(l)
1901      continue
c
1905    continue
c
        call set_zenith(nlev,Ntheta,radius,Rs,
     -                  b,au,alt,rangle,r0,mu)
c
c        if( init .eq. 0 ) then
c          write(*,'(/,A20,F4.2,A28,F4.2,A7,F6.2,A4)') 
c     -          'Minimum altitude of ', hmin0,
c     -          ' km for impact parameter at ', zbar(kmin), 
c     -          ' km at ', wl, ' um.'
c        endif
c
      endif
c
      return
      end
