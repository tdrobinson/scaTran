      subroutine path_rayt(nlev,radius,dlnnrdz,alt,
     -                     hmin,h0,rangle,rad,path)
c
ccccccccccccccccccccccccccc   path_rayt   cccccccccccccccccccccccccccccc
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
ccccccccccccccccccccccccccc   path_rayt   cccccccccccccccccccccccccccccc
c
      implicit none
c
      include 'param.inc'
c
      integer nlev,nlyr,k,l,sflag,iz
      real radius,alt(kp),zbar(kp),dz(kp),h0
      double precision s1,h1,theta1,phi1,xi1,
     -                 s2,h2,theta2,phi2,xi2
      double precision dlnnrdz(kp),x0,z0,r0,ds,hbar
      double precision path(kp),rangle,rad,hmin
c
      nlyr = nlev - 1
c
c*TDR layer mid-heights and thickness
c
      do 901 k=1,nlyr
        zbar(k) = 0.5*(alt(k)+alt(k+1))
        dz(k)   = alt(k)-alt(k+1)
901   continue
c
c*TDR initialize matrixes.  the path distribution, for a tangent 
c     ray at altitude 'h0', has a total path in layer 'l' of
c     path(l)*dz(l).  having any pathlength at l = nlev 
c     indicates the ray struck the surface.  also initialize the 
c     area matrix, and certain matrixes that get averaged over 
c     in the monte carlo.
c
      do 1101 l=1,nlev
        path(l)  = 0.0
1101  continue
c
c*TDR initialize ray on cartesian grid. at the top of the atmosphere, 
c     for a tangent height of zbar, traveling towards the stellar disk. 
c     a small factor of dz(1)/1e5 places ray "just" inside of the 
c     atmosphere.  note that we've shifted the z-coordinate to be 0 at 
c     the planet's center.
c
      z0    = DBLE(radius) + h0
      x0    = ((DBLE(radius) + alt(1))**2 - z0**2)**0.5 - dz(1)/1.d3
      r0    = (x0**2 + z0**2)**0.5
      phi1  = ATAN(x0/z0)
      s1    = 0.
      h1    = r0 - radius
      h2    = h1
      hmin  = h2
      theta1= ACOS(-x0/r0)
      xi1   = 0.
      iz    = 1
c
c*TDR flag to indicate surface has been struck
c
      sflag = 1
c
      do while ( h2 .lt. alt(1) .and. sflag .eq. 1 )
c
c*TDR   set path increment to be 1/e of current layer thickness,
c       to prevent multiple layer crossings.
c
        ds  = dz(iz)/exp(1.)
c
c*TDR   use ray tracing scheme to step ray forward through a 
c       pathlength of ds
c
        call raytrace(nlyr,radius,zbar,dlnnrdz,ds,
     -                s1,h1,theta1,phi1,xi1,
     -                s2,h2,theta2,phi2,xi2)
c
c*TDR   find nearest altitude layer for hbar, and add incremental
c       path length to stored total path
c
        hbar = 0.5*(h1 + h2)
c
        if( h2 .lt. hmin ) hmin = h2
c
        if( hbar .lt. alt(1) ) then
c
c*TDR     using path(k,nlev) to store rays that strike the surface
c
          if( h2 .le. 0 ) then
            path(nlev) = 1.
            sflag = 0
          else
c
c*TDR     otherwise add to path distribution
c
            if( iz .ne. nlyr .and. hbar .lt. zbar(iz+1) ) then
              path(iz+1) = path(iz+1) + ds/dz(iz)
            else if( iz .ne. 1 .and. hbar .gt. zbar(iz-1) ) then 
              path(iz-1) = path(iz-1) + ds/dz(iz-1)
            else
              path(iz) = path(iz) + ds/dz(iz)
            endif
c
          endif
c
        endif
c
c*TDR   update quantities for next increment
c
        s1     = s2
        h1     = h2
        theta1 = theta2
        phi1   = phi2
        xi1    = xi2
c
c*TDR   store refraction angle and x-coordinate at exit
c
        rangle = -xi2
        rad    = (h2 + radius)*cos(phi2)
c
c*TDR   set layer for current ray height
c
        if( h2 .lt. alt(iz+1) ) then
          iz = iz+1
        else if( h2 .gt. alt(iz) ) then
          iz = iz-1
        endif
c
      enddo
c
      return
      end