      subroutine transit_optd(nlev,copi,dtau_ex,path,tau_ab,tau_sc)
c
cccccccccccccccccccccccccc   transit_optd   cccccccccccccccccccccccccccc
cc                                                                    cc
cc    p u r p o s e :                                                 cc
cc      Computes the slant absorption and scattering optical depths.  cc       
cc                                                                    cc
cc    i n p u t s :                                                   cc
cc      nlev   -  number of atmospheric levels                        cc
cc      copi   - single scattering co-albedo                          cc
cc    dtau_ex  - differential extinction optical depth                cc
cc      pathg  - dimensionless geometric path distribution            cc
cc                                                                    cc
cc    o u t p u t s :                                                 cc
cc     tau_ab  - slant absorption optical depth                       cc
cc     tau_sc  - slant scattering optical depth                       cc
cc                                                                    cc
cc    a u t h o r:                                                    cc
cc           Tyler D. Robinson (robinson.tyler.d@gmail.com)           cc
cc                                                                    cc
cccccccccccccccccccccccccc   transit_optd   cccccccccccccccccccccccccccc
c
      implicit none
c
      include 'param.inc'
c
      integer nlev,nlyr,k,l
      real dtau_sc(kp),dtau_ab(kp),copi(kp),dtau_ex(kp)
      real tau_ab(kp),tau_sc(kp)
      double precision path(kp,kp)
c
      nlyr = nlev - 1
c
c
c*TDR quantities needed for transit spectrum calculation for cases both
c     with and without atmospheric refraction and/or scattering
c
      do 1001 k=1,nlyr
c
c*TDR   scattering and absorption optical depths
c
        dtau_sc(k) = (1.0 - copi(k))*dtau_ex(k)
        dtau_ab(k) = copi(k)*dtau_ex(k)
c
1001  continue
c
c
c*TDR   compute absorption and scattering optical depths along each 
c       tangent ray.
c
      do 2101 k=1,nlyr
        tau_ab(k) = 0.0d0
        tau_sc(k) = 0.0d0
c
        do 2001 l=1,nlyr
          tau_ab(k) = tau_ab(k) + dtau_ab(l)*path(k,l)
          tau_sc(k) = tau_sc(k) + dtau_sc(l)*path(k,l)
2001    continue
c
2101  continue
c
      return
      end