      subroutine raytrace(nlay,radius,zbar,dlnnrdz,ds,
     -                    s1,h1,theta1,phi1,xi1,
     -                    s2,h2,theta2,phi2,xi2)
c
ccccccccccccccccccccccccccc    raytrace    ccccccccccccccccccccccccccccc
cc                                                                    cc
cc    p u r p o s e :                                                 cc
cc      ray tracing scheme, which increments the ray through          cc
cc      a distance ds according to van der Werf (2008; Applied        cc
cc      Optics, 47(2):153-156).  using initial  quantities            cc
cc      (i.e., s1, h1, theta1, phi1, xi1), derivatives in a           cc
cc      fourth order Runge Kutta scheme are computed, and then        cc
cc      used to produce updated quantities (i.e., s2, h2, theta2,     cc
cc      phi2, xi2).  note that ds, alt, and radius should             cc
cc      use the same length unit, and dlnnrdz should have inverse     cc
cc      units of this length (e.g., km**-1).  note also that vdW's    cc
cc      'beta' (the tilt angle) has been swapped to 'theta' (the      cc
cc      local zenith angle) via beta = theta + pi/2.                  cc
cc                                                                    cc
cc    i n p u t s:                                                    cc
cc         nlay --    number of atmospheric layers                    cc
cc       radius --    planetary radius                                cc
cc         zbar --    altitudes of layer midpoints                    cc
cc      dlnnrdz --    d ln(index of refraction) / dz, defined         cc
cc                    at atmospheric layers                           cc
cc           ds --    path length increment                           cc
cc           s1 --    current path length                             cc
cc           h1 --    current altitude                                cc
cc       theta1 --    current local zenith angle (0 => up)            cc
cc         phi1 --    current polar angle measured from planet        cc
cc                    center                                          cc
cc         xi1 --    current "refraction integral"                    cc
cc                                                                    cc
cc    o u t p u t s:                                                  cc
cc           s2 --    updated path length                             cc
cc           h2 --    updated altitude                                cc
cc       theta2 --    updated local zenith angle                      cc
cc         phi2 --    updated polar angle measured from planet        cc
cc                    center                                          cc
cc         xi2 --    updated "refraction integral"                    cc
cc                                                                    cc
cc    a u t h o r:                                                    cc
cc           Tyler D. Robinson (robinson.tyler.d@gmail.com)           cc
cc                                                                    cc
ccccccccccccccccccccccccccc    raytrace    ccccccccccccccccccccccccccccc
c
      implicit none
c
      include 'param.inc'
c
      integer nlay
      logical planepar
      real radius,zbar(kp)
      double precision dlnnrdz(kp),ds
      double precision s1,s2,h1,h2,theta1,theta2,phi1,phi2,xi1,xi2
      double precision s_i,h_i,theta_i,phi_i,curv
      double precision k_h1,k_h2,k_h3,k_h4
      double precision k_theta1,k_theta2,k_theta3,k_theta4
      double precision k_phi1,k_phi2,k_phi3,k_phi4
      double precision k_xi1,k_xi2,k_xi3,k_xi4
      double precision curvature,dhds,dthetads,dphids,dxids
c
c*TDR allow for plane parallel case, where planetary 
c     radius goes to infinity. usually set to false,
c     meaning that the atmosphere is *not* plane par.
c
      planepar = .false.
c
c*TDR compute derivatives for Runge Kutta scheme
c
      s_i      = s1
      h_i      = h1
      theta_i  = theta1
      phi_i    = phi1
      curv     = curvature(nlay,zbar,dlnnrdz,theta_i,h_i)
      k_h1     = dhds(theta_i)
      k_theta1 = dthetads(radius,curv,h_i,theta_i,planepar)
      k_phi1   = dphids(radius,h_i,theta_i,planepar)
      k_xi1    = dxids(curv)
c
      s_i      = s1 + 0.5*ds
      h_i      = h1 + 0.5*k_h1*ds
      theta_i  = theta1 + 0.5*k_theta1*ds
      phi_i    = phi1 + 0.5*k_phi1*ds
      curv     = curvature(nlay,zbar,dlnnrdz,theta_i,h_i)
      k_h2     = dhds(theta_i)
      k_theta2 = dthetads(radius,curv,h_i,theta_i,planepar)
      k_phi2   = dphids(radius,h_i,theta_i,planepar)
      k_xi2    = dxids(curv)
c
      s_i      = s1 + 0.5*ds
      h_i      = h1 + 0.5*k_h2*ds
      theta_i  = theta1 + 0.5*k_theta2*ds
      phi_i    = phi1 + 0.5*k_phi2*ds
      curv     = curvature(nlay,zbar,dlnnrdz,theta_i,h_i)
      k_h3     = dhds(theta_i)
      k_theta3 = dthetads(radius,curv,h_i,theta_i,planepar)
      k_phi3   = dphids(radius,h_i,theta_i,planepar)
      k_xi3    = dxids(curv)
c
      s_i      = s1 + 0.5*ds
      h_i      = h1 + 0.5*k_h3*ds
      theta_i  = theta1 + 0.5*k_theta3*ds
      phi_i    = phi1 + 0.5*k_phi3*ds
      curv     = curvature(nlay,zbar,dlnnrdz,theta_i,h_i)
      k_h4     = dhds(theta_i)
      k_theta4 = dthetads(radius,curv,h_i,theta_i,planepar)
      k_phi4   = dphids(radius,h_i,theta_i,planepar)
      k_xi4    = dxids(curv)
c
c*TDR use derivatives to update quantities
c
      s2     = s1 + ds
      h2     = h1 + (k_h1 + 2*k_h2 + 2*k_h3 + k_h4)*ds/6
      theta2 = theta1 + (k_theta1 + 2*k_theta2 + 
     -                   2*k_theta3 + k_theta4)*ds/6
      phi2   = phi1  + (k_phi1 + 2*k_phi2 + 2*k_phi3 + k_phi4)*ds/6
      xi2    = xi1  + (k_xi1 + 2*k_xi2 + 2*k_xi3 + k_xi4)*ds/6
c
      return
      end
c
c*TDR local curvature, from vdW08 Eqn. 2a
c
      double precision function curvature(nlay,z,dlnnrdz,theta,h)
        include 'param.inc'
c
        integer nlay
        real z(kp),deriv
        double precision theta,h
        double precision dlnnrdz(kp)
c
c*TDR   if current height is above top of atmosphere, 
c       set value equal to that at top
c
        if (h .ge. z(1)) then
          deriv = 0.
        else
          call xyinterp(z,REAL(dlnnrdz),REAL(h),deriv,
     -                  kp,1,nlay,1,1)
        endif
c
        curvature = SIN(theta)*deriv
c
        return
      end
c
c*TDR derivative of h with respect to s, from
c     vdW08 Table 1, with s as integration variable
c
      double precision function dhds(theta)
        double precision theta
c
        dhds = COS(theta)
c
        return
      end
c
c*TDR derivative of theta with respect to s, from
c     vdW08 Table 1, with s as integration variable
c
      double precision function dthetads(radius,curv,h,theta,planepar)
        logical planepar
        real radius
        double precision theta,h,curv
c
        if( planepar ) then
          dthetads = curv
        else
          dthetads = -(SIN(theta)/(radius + h) + curv)
        endif
c
        return
      end
c
c*TDR derivative of phi with respect to s, from
c     vdW08 Table 1, with s as integration variable
c
      double precision function dphids(radius,h,theta,planepar)
        logical planepar
        real radius
        double precision theta,h
c
        if( planepar ) then
          dphids = 0.
        else
          dphids = SIN(theta)/(radius + h)
        endif
c
        return
      end
c
c*TDR derivative of xi with respect to s, from
c     vdW08 Table 1, with s as integration variable
c
      double precision function dxids(curv)
        double precision curv
c
        dxids = curv
c
        return
      end