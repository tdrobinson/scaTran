      subroutine set_zenith(nlev,Ntheta,radius,Rs,
     -                      b,au,alt,rangle,r0,mu0)
c
ccccccccccccccccccccccccccc   set_zenith   ccccccccccccccccccccccccccccc
cc                                                                    cc
cc    p u r p o s e :                                                 cc
cc      For rays exiting atmosphere with coordinate r0 and refraction cc
cc      angle 'rangle', determines if and where they strike the       cc
cc      stellar disk.                                                 cc
cc                                                                    cc
cc    i n p u t s :                                                   cc
cc      nlev   - number of atmospheric levels                         cc
cc      Ntheta - number of angular integration points                 cc
cc      ng0    - bin index                                            cc
cc      radius - planetary radius (km)                                cc
cc      Rs     - stellar radius (km)                                  cc
cc      b      - impact parameter (km)                                cc
cc      au     - semi-major axis (au)                                 cc
cc      alt    - altitude profile (km)                                cc
cc      rangle - refraction angle for tangent rays                    cc
cc      r0     - x-coordinate where rays exit atmosphere (km)         cc
cc                                                                    cc
cc    o u t p u t s :                                                 cc
cc      mu0    - zenith angle on star where ray strikes disk (equal   cc
cc               to -1 if ray misses disk)                            cc
cc                                                                    cc
cc    a u t h o r:                                                    cc
cc           Tyler D. Robinson (robinson.tyler.d@gmail.com)           cc
cc                                                                    cc
ccccccccccccccccccccccccccc   set_zenith   ccccccccccccccccccccccccccccc
c
      implicit none
c
      include 'param.inc'
c
      integer nlev,nlyr,k,ith,Ntheta
      real radius,alt(kp),au
      double precision rangle(kp),r0(kp),b
      double precision pi,theta,dtheta,mux,muy,muz,bq,cq,kmau,Rs
      double precision s,z,y
      double precision mu0(kp,mxth)
c
      pi  = ACOS(-1.d0)
c
      nlyr = nlev - 1
c
      dtheta = pi/REAL(Ntheta)
c
c*TDR au in km
c
      kmau   = 1.49597870700e8
c
      do 1501 k=1,nlyr
c
        do 1001 ith=1,Ntheta
c
          mu0(k,ith)  = -1.
c
          theta = REAL(ith)*dtheta - dtheta/2
c
c*TDR     direction cosines
c
          muz = sin(rangle(k))*cos(theta)
          muy = sin(rangle(k))*sin(theta)
          mux = cos(rangle(k))
c
c*TDR     terms for solving quadratic, determining if ray ever
c         strikes the stellar disk. we've written the (x,y,z) path
c         of the ray as x(s) = r0 - s*mux, y(s) = y0 - s*muy, and 
c         z(s) = z0 - s*muz.  to (hopefully) minimize round-off 
c         errors, we expand r0, y0, and z0 around s = a, where 'a'
c         is the semi-major axis.  we then find if 
c         r(s) = sqrt(x**2 + y**2 + z**2) is ever less than the 
c         stellar radius.  note that: 
c         z0 = -b + r0*cos(theta) - a*muz, 
c         y0 = r0*sin(theta) - a*muy, and r0 = -h0 + (1-mux)*a, 
c         where h0**2 = (Rp + alt(1))**2 - r0**2, and r0 is the 
c         y-z coordinate for where the ray exited the atmosphere.
c
          bq = -2*(muz*(r0(k)*cos(theta) - b) + muy*r0(k)*sin(theta) + 
     -             mux*(- ( (radius + alt(1))**2 - r0(k)**2 )**0.5 )
     -             + (au*kmau)*(mux - 1.d0))
          cq = (b**2 - 2*r0(k)*b*COS(theta)+(radius+alt(1))**2) +
     -         2*(1.0d0-mux)*(au*kmau)**2 - 
     -         2*au*kmau*(muz*(-b + r0(k)*cos(theta)) + 
     -                    muy*r0(k)*sin(theta) + 
     -                    (1.d0-mux)*((radius + alt(1))**2 
     -                                - r0(k)**2)**0.5) -
     -         (Rs)**2
c
c*TDR     if the quadratic has no solution, then the ray never 
c         strikes the stellar disk.
c
          if( bq**2 - 4*cq .ge. 0 ) then
c
c*TDR       solve quadratic for path location where disk is struck
c
            s  = 0.5*(-bq + (bq**2 - 4*cq)**0.5)
            z  = -b + r0(k)*cos(theta) - (s + au*kmau)*muz
            y  = r0(k)*sin(theta) - (s + au*kmau)*muy
c
            mu0(k,ith) = (1 - (z**2 + y**2)/Rs**2)**0.5
c
          endif
c
1001    continue
c
1501  continue
c
      return
      end
