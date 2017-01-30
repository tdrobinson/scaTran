      subroutine path_mc(iseed,nlev,irefract,Ntheta,Rs,b,au,
     -                   radius,dlnnrdz,a_sc,g,alt,h0,mu0,path)
c
ccccccccccccccccccccccccccc    path_mc    cccccccccccccccccccccccccccccc
cc                                                                    cc
cc    p u r p o s e :                                                 cc
cc      Computes the path distribution using ray tracing.             cc
cc                                                                    cc
cc    i n p u t s :                                                   cc
cc      nlev   - number of atmospheric levels                         cc
cc      radius - planetary radius (km)                                cc
cc      alt    - altitude profile (km)                                cc
cc                                                                    cc
cc    o u t p u t s :                                                 cc
cc      rangle - refraction angle for each tangent ray                cc
cc      x0     - radial distance at which ray exits the atmosphere    cc
cc      path   - dimensionless path distribution                      cc
cc                                                                    cc
cc    a u t h o r:                                                    cc
cc           Tyler D. Robinson (robinson.tyler.d@gmail.com)           cc
cc                                                                    cc
ccccccccccccccccccccccccccc    path_mc    cccccccccccccccccccccccccccccc
c
      implicit none
c
      include 'param.inc'
c
      integer irefract,Ntheta
      integer nlev,nlyr,k,l,iz,sflag,eflag,iseed,nsub,is,ith
      real zbar(kp),dz(kp),h0,alt(kp),radius,au,g(kp)
      double precision pi,a_sc(kp),dlnnrdz(kp)
      double precision x,y,z,r,h,ys,zs,mu,theta,dtheta,curv
      double precision mux,muy,muz,muxp,muyp,muzp,sinth,costh
      double precision rnum,tau0,tau,alpha,s,ds,mus,phis,sgn
      double precision r0,bq,cq,Rs,b,kmau,dspath,curvature
      double precision path(kp),mu0(mxth)
c
      nlyr = nlev - 1
c
      pi  = ACOS(-1.d0)
c
c*TDR au in km
c
      kmau   = 1.49597870700e8
c
      dtheta = pi/REAL(Ntheta)
c
c*TDR layer mid-heights and thickness
c
      do 901 k=1,nlyr
        zbar(k) = 0.5*(alt(k)+alt(k+1))
        dz(k)   = alt(k)-alt(k+1)
901   continue
c
c*TDR initialize matrixes.  the path distribution, for a tangent 
c     ray at altitude h0, has a total path in layer 'l' of
c     path(l)*dz(l).  having any pathlength at l = nlev 
c     indicates the ray struck the surface.  also initialize the 
c     area matrix, and certain matrixes that get averaged over 
c     in the monte carlo.
c
      do 1101 l=1,nlev
        path(l)  = 0.0
1101  continue
c
c*TDR initialize ray on cartesian grid, with planet center at
c     (0,0,0).  ray begins at far right (+z), at the top of 
c     the atmosphere for a tangent height of zbar.  factor of
c     dz(1)/1e3 ensure photon is "just" inside the atmosphere
c
      z     = dble(radius + h0)
      x     = ((dble(radius) + alt(1))**2 - z**2)**0.5 
     -        - dz(1)/1.d3
      y     = 0.
      r     = (x**2 + y**2 + z**2)**0.5
      h     = r - dble(radius)
c
c*TDR directs photon along the -x axis
c
      muz   = 0.d0
      mux   = (1.-muz**2)**0.5*cos(pi)
      muy   = (1.-muz**2)**0.5*sin(pi)
      iz    = 1
c
c*TDR initialize angle-dependent location struck on stellar disk
c
      do 2301 ith=1,Ntheta
        mu0(ith) = -1.d0
2301  continue
c
c*TDR flag that indicates the photon is absorbed, and flag for
c     photon exiting atmosphere
c
      sflag = 1
      eflag = 1
c
      do while ( sflag .eq. 1 .and. eflag .eq. 1 )
c
c*TDR   path traveled comes from random sampling of optical
c       depth
c
        rnum = ran(iseed)
        tau0 = -LOG(1 - rnum)
c
c*TDR   ray trace the photon through the atmosphere until it
c       has encountered tau0 of optical depths or leaves the
c       atmosphere
c
        tau = 0.d0
        do while ( tau .lt. tau0 .and. h .lt. alt(1) .and. h .gt. 0 )
c
c*TDR     absorption coefficient for current layer
c
          alpha = a_sc(iz)
c
c*TDR     set pathlength
c
          s = dspath(nlev,radius,alt,iz,x,y,z,r,mux,muy,muz)
c
c*TDR     if optical depth tau0 would be achieved before ds is
c         traveled, then adjust ds to just achieve tau0
c
          if( tau + alpha*s .ge. tau0 .and. alpha .gt. 0. ) then
            s = (tau0-tau)/alpha
          endif
c
c*TDR     trace path through ds.  path is sub-divided into 
c         smaller increments when doing refraction.
c
          if( irefract .eq. 0 ) then
c
            x = x + mux*s
            y = y + muy*s
            z = z + muz*s
            r = (x**2 + y**2 + z**2)**0.5
            h = r - dble(radius)
c
          else
c
c*TDR       sub-intervals set by current layer thickness and 
c           whether or not refraction matters much over the 
c           path 's'.
c
            mu    = mux*x/r + muy*y/r + muz*z/r
            theta = acos(mu)
            curv  = curvature(nlyr,zbar,dlnnrdz,theta,h)
            if( curv*s .lt. 1.e-5 ) then
              nsub = 1
            else
              nsub = nint(s/dz(iz))
              if( nsub .eq. 0 ) nsub = 1
            endif
            ds = s/dble(nsub)
            s  = 0.d0
c
c*TDR       loop over sub-intervals
c
            is = 1
            do while ( is .le. nsub .and. h .gt. 0 )
c
c*TDR         angle between current direction of travel and
c             radial unit vector
c
              mu     = mux*x/r + muy*y/r + muz*z/r
c
c*TDR         apply curvature through ds
c
              theta  = acos(mu)
              curv   = curvature(nlyr,zbar,dlnnrdz,theta,h)
              dtheta = -curv*ds
c
c*TDR         update direction cosines
c
              call mu_refract(x,y,z,r,mux,muy,muz,
     -                        dtheta,muxp,muyp,muzp)
c
              mux = muxp
              muy = muyp
              muz = muzp
c
c*TDR         increment through ds
c
              x = x + mux*ds
              y = y + muy*ds
              z = z + muz*ds
              r = (x**2 + y**2 + z**2)**0.5
              h = r - DBLE(radius)
c
              is = is + 1
              s  = s + ds
c
            enddo
c
          endif
c
c*TDR     add contributions to stored path and optical depth
c
          path(iz) = path(iz) + s/dz(iz)
          tau = tau + alpha*s
c
c
c*TDR     set layer for current ray height
c
          if( h .le. alt(iz+1) ) then
            iz = iz+1
          else if( h .ge. alt(iz) ) then
            iz = iz-1
          endif
c
        enddo
c
c*TDR   if the loop exits with the photon height being larger
c       than the top of atmosphere, then the photon passed 
c       through the remainder of the atmosphere without an 
c       extinction event. or if the photon strikes the 
c       surface.
c
        if( h .ge. alt(1) .or. h .le. 0 ) then
c
          if( h .le. 0 ) then
            sflag = 0
            path(nlev) = DBLE(1)
          else
c
c*TDR       ray has exited the atmosphere, do an integral over
c           angle to determine the area blocked by arcs of the 
c           atmosphere, but only if ray is traveling towards the
c           star (i.e., mux is less than zero).
c
            eflag = 0
c
            if( mux .lt. 0 ) then
c
              muxp = mux
              muyp = muy
              muzp = muz
c
              do 2401 ith=1,Ntheta
c
                theta = REAL(ith)*dtheta - dtheta/2
c
c*TDR           direction cosines updated for location on arc
c
                if( muxp**2 .lt. 1.0d0 ) then
                  costh = muzp/((1.0d0 - muxp**2)**0.5)
                else
                  costh = 1.0d0
                endif
                if( costh**2 .le. 1.0d0 ) then
                  sinth = (1 - costh**2)**0.5
                else
                  costh = 1.0d0
                  sinth = 0.0d0
                endif
c
                mux = -muxp
                if( muxp**2 .lt. 1.0d0 ) then
                  muy  = ((1 - mux**2)**0.5)*(sinth*cos(theta)
     -                                        + costh*sin(theta))
                  muz  = ((1 - mux**2)**0.5)*(costh*cos(theta)
     -                                        - sinth*sin(theta))
                else
                  muy  = 0.d0
                  muz  = 0.d0
                endif
c
c*TDR           terms for solving quadratic, same as above.
c
                r0 = (y**2 + z**2)**0.5
                bq = -2*( muz*(r0*cos(theta) - b) + muy*r0*sin(theta) + 
     -                    mux*(- ((radius + alt(1))**2 - r0**2)**0.5 ) +
     -                    (au*kmau)*(mux - 1.d0) )
                cq = (b**2 - 2*r0*b*cos(theta) + (radius+alt(1))**2) +
     -               2*(1.0d0-mux)*(au*kmau)**2 - 
     -               2*au*kmau*(muz*(-b+r0*cos(theta)) +
     -                 muy*r0*sin(theta) + 
     -                 (1.d0-mux)*((radius+alt(1))**2-r0**2)**0.5) - 
     -               (Rs)**2
c
c*TDR           if the quadratic has no solution, then the ray
c               never strikes the stellar disk.
c
                if( bq**2 - 4*cq .ge. 0 ) then
c
c*TDR             solve quadratic for path location where disk
c                 is struck
c
                  s  = 0.5*(-bq + (bq**2 - 4*cq)**0.5)
                  zs = -b + r0*cos(theta) - (s + au*kmau)*muz
                  ys = r0*sin(theta) - (s + au*kmau)*muy
c
                  mu0(ith) = (1 - (zs**2+ys**2)/Rs**2)**0.5
c
                endif
c
2401          continue
c
            endif
c
          endif
c
c*TDR     otherwise an extinction event occurs within the
c         atmosphere
c
        else            
c
c*TDR     photon is scattered, determine scattering angle 
c         -- currently using the H-G phase function.  the 
c         goto statement is due to tiny numerical errors in
c         computing mus, where a particular rnum causes mus 
c         to be a tiny smidge above unity.  in that case, 
c         we just re-do the sampling.
c
2501      continue
          rnum  = ran(iseed)
          mus   = (1 + g(iz)**2 - ((1 - g(iz)**2)/(1 - 
     -             g(iz) + 2*g(iz)*rnum))**2.)/(2*g(iz))
          if( g(iz) .lt. 0.001 ) mus = 1 - 2*ran(iseed)
          if( ABS(mus) .gt. 1. ) goto 2501
          phis  = 2*pi*ran(iseed)
c
c*TDR     special case for photons travelling straight down or
c         up, otherwise use standard logic
c
          if( ABS(muz) .gt. 0.99999 ) then
            sgn = 1.
          if( muz .lt. 0 ) sgn = -1.
            mux = cos(phis)*(1-mus**2)**0.5
            muy = sgn*sin(phis)*(1-mus**2)**0.5
            muz = sgn*mus
          else
            muxp = ((1-mus**2)**0.5/(1-muz**2)**0.5)*
     -             (mux*muz*cos(phis) - muy*sin(phis)) + mux*mus
            muyp = ((1-mus**2)**0.5/(1-muz**2)**0.5)*
     -             (muy*muz*cos(phis) + mux*sin(phis)) + muy*mus
            muzp = muz*mus - cos(phis)*(1-muz**2)**0.5*(1-mus**2)**0.5
c

            mux = muxp
            muy = muyp
            muz = muzp
c                  
          endif
c
        endif
c
      enddo
c
      return
      end
c
c*TDR function to determine pathlength to nearest atmospheric layer
c
      double precision function dspath(nlev,radius,alt,iz,
     -                                 x,y,z,r,mux,muy,muz)
c
        include 'param.inc'
c
        integer iz,nlev
        real radius,alt(kp)
        double precision bq,cq,cq1,sgn
        double precision x,y,z,r,dz
        double precision mux,muy,muz
c
        dz = alt(iz)-alt(iz+1)
c
c*TDR   given current trajectory, determine linear distance to level iz
c       or iz+1, i.e., where r(ds) = radius + alt(iz), which involves 
c       solving a quadratic eqn.
c
        bq  = 2*(mux*x + muy*y + muz*z)
        cq  = r**2 - (DBLE(radius) + alt(iz))**2
        cq1 = r**2 - (DBLE(radius) + alt(iz+1))**2
c
        sgn = 1.
        if( -bq - (bq**2 - 4*cq)**0.5 .gt. 0 ) sgn = -1.
        dspath  = 0.5*(-bq + sgn*(bq**2 - 4*cq)**0.5)
c
        if( bq**2 - 4*cq1 .gt. 0 ) then
          sgn = 1.
          if( -bq - (bq**2 - 4*cq1)**0.5 .gt. 0 ) sgn = -1.
          if( 0.5*(-bq + sgn*(bq**2 - 4*cq1)**0.5) .gt. 0 )
     -      dspath = MIN(ds,0.5*(-bq + sgn*(bq**2 - 4*cq1)**0.5))
        endif
c
c*TDR   special case when right at level boundaries, which gives ds=0
c
        if( dspath .lt. 1.e-5*dz ) then
c
c*TDR   implies at level iz
c
          if( cq .lt. 1.e-5 ) then
            if( bq .gt. 0 ) then
              if( iz .eq. 1 ) then
                dspath = (alt(iz)-alt(iz+1))/exp(1.)
              else
                cq  = r**2 - (DBLE(radius) + alt(iz-1))**2
                sgn = 1.
                if( -bq - (bq**2 - 4*cq)**0.5 .gt. 0 )sgn = -1.
                dspath  = 0.5*(-bq + sgn*(bq**2 - 4*cq)**0.5)
              endif
            else
              sgn = 1.
              if( -bq - (bq**2 - 4*cq1)**0.5 .gt. 0 )sgn = -1.
              dspath  = 0.5*(-bq + sgn*(bq**2 - 4*cq1)**0.5)
            endif
c
c*TDR     implies at level iz+1
c
          else if( cq1 .lt. 1.e-5 ) then
            if( bq .gt. 0 ) then
              sgn = 1.
              if( -bq - (bq**2 - 4*cq)**0.5 .gt. 0 ) sgn = -1.
              dspath  = 0.5*(-bq + sgn*(bq**2 - 4*cq)**0.5)
            else
              if( iz+1 .eq. nlev ) then
                dspath = (alt(iz)-alt(iz+1))/exp(1.)
              else
                cq1 = r**2 - (DBLE(radius) + alt(iz+2))**2
                sgn = 1.
                if( -bq - (bq**2 - 4*cq)**0.5 .gt. 0 ) sgn = -1.
                dspath  = 0.5*(-bq + sgn*(bq**2 - 4*cq)**0.5)
              endif
            endif
          endif
        endif
c
c*TDR   add small factor to prevent ray from being "stuck" at a boundary
c
        dspath = dspath + dz/1.d3
        return
      end
c
c*TDR subroutine to determine updated direction cosines given a
c     direction change due to refraction
c
      subroutine mu_refract(x,y,z,r,mux,muy,muz,dtheta,muxp,muyp,muzp)
c
        include 'param.inc'
c
        double precision x,y,z,r,mux,muy,muz,dtheta
        double precision muxp,muyp,muzp,muxp1,muyp1,muzp1
        double precision mur,mur1
        double precision a1,a2,a3,b1,b2,c1,c2
        double precision aq,bq,cq
c
        a1 = muy*z/r - muz*y/r
        a2 = muz*x/r - mux*z/r
        a3 = mux*y/r - muy*x/r
c
        b1 = muy - a2/a1*mux
        b2 = muz - a3/a1*mux
c
        c1 = a1*muy - a2*mux
        c2 = a1*c1
c
        aq = 1 + a3**2/a1**2 + b2**2/c1**2*(a1**2 + a2**2)
     -       - 2*b2*a2*a3/c1/a1
        bq = 2*COS(dtheta)/c2*(a2*a3 - a1*b2/c1*(a1**2 + a2**2))
        cq = COS(dtheta)**2/c1**2*(a1**2 + a2**2) - 1
c
        if( bq**2 - 4*aq*cq .ge. 0 ) then
c
c*TDR     different cases for +/- in quadratic
c
          muzp  = (-bq - (bq**2 - 4*aq*cq)**0.5)/2/aq
          muyp  = (COS(dtheta) - b2*muzp)/b1
          muxp  = -(a2*muyp + a3*muzp)/a1
c
          muzp1 = (-bq + (bq**2 - 4*aq*cq)**0.5)/2/aq
          muyp1 = (COS(dtheta) - b2*muzp)/b1
          muxp1 = -(a2*muyp + a3*muzp)/a1
c
c*TDR     select case that yields more "downward" ray, by computing
c         the angle between the new direction and the unit radial vector
c
          mur   = muxp*(x/r)  + muyp*(y/r)  + muzp*(z/r)
          mur1  = muxp1*(x/r) + muyp1*(y/r) + muzp1*(z/r)
          if( mur1 .lt. mur ) then
            muxp = muxp1
            muyp = muyp1
            muzp = muzp1
          endif
        else
          muzp = muz
          muyp = muy
          muzp = muz
        endif
c
        return
      end