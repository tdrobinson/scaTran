      subroutine path_geom(nlev,radius,alt,h0,path)
c
ccccccccccccccccccccccccccc   path_geom   cccccccccccccccccccccccccccccc
cc                                                                    cc
cc    p u r p o s e :                                                 cc
cc      Computes the path distribution in the geometric limit.        cc
cc                                                                    cc
cc    i n p u t s :                                                   cc
cc      nlev   - number of atmospheric levels                         cc
cc      radius - planetary radius (km)                                cc
cc      alt    - altitude profile (km)                                cc
cc                                                                    cc
cc    o u t p u t s :                                                 cc
cc      pathg  - dimensionless geometric path distribution            cc
cc                                                                    cc
cc    a u t h o r:                                                    cc
cc           Tyler D. Robinson (robinson.tyler.d@gmail.com)           cc
cc                                                                    cc
ccccccccccccccccccccccccccc   path_geom   cccccccccccccccccccccccccccccc
c
      implicit none
c
      include 'param.inc'
c
      integer nlev,nlyr,l
      real radius,alt(kp),h0,dz
      double precision path(kp),s_kl,s_kl1
c
      nlyr = nlev - 1
c
c*TDR loop over atmospheric layers
c
      do 1051 l=1,nlyr
c
        dz   = alt(l)-alt(l+1)
c
c*TDR   path distribution in geometric limit.  see equations on
c       Robinson and Fortney (2017)
c
        if( h0 .lt. alt(l) ) then
          s_kl = ( (radius+alt(l))**2 - (radius+h0)**2 )**0.5
          if ( h0 .gt. alt(l+1) ) then
            s_kl1 = 0.0
          else
            s_kl1 = ( (radius+alt(l+1))**2 - (radius+h0)**2 )**0.5
          endif
c
          path(l) = 2*(s_kl - s_kl1)/dz
c
        else
          path(l) = 0.d0
        endif
c
1051  continue
c
      return
      end
