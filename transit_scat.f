      subroutine transit_scat(nlev,iupath,irefract,itrnst,ncomp,
     -                        icomp,volmix,radstar,bimpact,radius,au,
     -                        p,t,alt,grav,ratm,wn,copi,dtau_ex,
     -                        gtot,pathg,Ntheta,do_mc,trans_b,mu_b,
     -                        dtau_ab_b,dtrans_dtau_b,dmu_dtau_b)
c
ccccccccccccccccccccccccccc  transit_scat  ccccccccccccccccccccccccccccc
cc                                                                    cc
cc    p u r p o s e :                                                 cc
cc      Runs a scattering Monte Carlo simulation for computing a      cc
cc      transit depth.  Atmospheric and optical properties are taken  cc
cc      as inputs.  Average transmission and stellar disk zenith      cc
cc      angles are computed, and the absorption optical depth for     cc
cc      this calculation is stored.  Jacobians for the transmission   cc
cc      and stellar disk zenith angle are output.                     cc
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
cc      grav   - gravity profile (m/s/s)                              cc
cc      ratm   - atmospheric specific gas constant (J/K/kmole)        cc
cc      wn     - wavenumber (cm**-1)                                  cc
cc      copi   - profile of layer single scattering co-albedos        cc
cc     dtau_ex - profile of layer extinction optical depths           cc
cc      gtot   - profile of layer asymmetry parameters                cc
cc      pathg  - geometric path distribution                          cc
cc      Ntheta - number of angular points in integration              cc
cc                                                                    cc
cc    o u t p u t s :                                                 cc
cc      do_mc  - profile of which altitudes have MC calculations      cc
cc     trans_b - altitude- and angle-dependent transmissions          cc
cc        mu_b - altitude and angle-dependent angles of stellar disk  cc
cc   dtau_ab_b - absorption optical depths at which transmissions and cc
cc               Jacobians are computed (i.e., most valid)            cc
cc dtrans_dtau_b - transmission Jacobian                              cc
cc  dmu_dtau_b - angle on stellar disk Jacobian                       cc
cc                                                                    cc
cc    a u t h o r:                                                    cc
cc           Tyler D. Robinson (robinson.tyler.d@gmail.com)           cc
cc                                                                    cc
ccccccccccccccccccccccccccc  transit_scat  ccccccccccccccccccccccccccccc
c
      implicit none
c
      include 'param.inc'
c
      integer nlev,nlyr,ng,iupath(kp),ncomp,icomp(6),irefract
      integer nmin,nmax,k,l,do_path,ith,nd(mxth),n,np,do_mc(kp),mc(kp)
      integer Ntheta,nmc,ioerr,itrnst,iseed
      real a0,rgas
      real radstar,bimpact,radius,volmix(6),alt(kp),grav(kp),p(kp)
      real ratm,copi(kp),dtau_ex(kp),gtot(kp),t(kp),au
      real wl,zbar(kp),dz(kp),g(kp),nr(kp),tol
      real nr0,g0(kp),dtau_sc0(kp)
      double precision wn
      double precision pi,radsol,kmau,Rs,b,d
      double precision mu0(mxth),mu(mxth)
      double precision dim(mxth),tmu(mxth),dtheta
      double precision dtau_sc(kp),dtau_ab(kp),a_sc(kp),tau_sc,tau_ab
      double precision dlnnrdz(kp),rangle(kp),r0(kp),trans
      double precision path0(kp),path(kp,kp),pathg(kp,kp),rangle0,rad0
      double precision ddim_dtau(kp,mxth),dmu_dtau(kp,mxth)
      double precision trans_b(kp,mxth),mu_b(kp,mxth)
      double precision dtrans_dtau_b(kp,kp,mxth)
      double precision dmu_dtau_b(kp,kp,mxth),dtau_ab_b(kp)
c
c*TDR saved variables
c
      save mc
c
c*TDR seed random number generator
c
      call system_clock(iseed)
c
c**** define avagadro's number (mks units: molecules/kmole)
c
      data a0/6.02297e26/
c
c**** define gas constant (mks units: J/K/kmole)
c
      data rgas/8.3144598e3/
c
      nlyr = nlev-1
c
      pi  = ACOS(-1.d0)
c
      wl = 1.e4/real(wn)
c
c*TDR tolerance for when to re-do ray tracing or when a monte carlo
c     simulation is favored over geometric or ray tracing approach.
c
      tol = 0.01
c
c*TDR minimum and maximum number of rays in monte carlo simulation
c
      nmin = 10000
      nmax = mxmc
c
c*TDR solar radius and au in km
c
      radsol = 695700.
      kmau   = 1.49597870700e8
      Rs     = DBLE(radstar*radsol)
      b      = DBLE(bimpact*Rs)
c
c*TDR quantities needed for transit spectrum calculation for cases both
c     with and without atmospheric refraction and/or scattering
c
      do 1001 k=1,nlyr
c
c*TDR   layer mean altitude and thickness (km)
c
        zbar(k) = 0.5*(alt(k)+alt(k+1))
        dz(k)   = alt(k)-alt(k+1)
c
c*TDR   scattering and absorption optical depths
c
        dtau_sc(k) = (1.0 - copi(k))*dtau_ex(k)
        dtau_ab(k) = copi(k)*dtau_ex(k)
c
c*TDR   store averages and Jacobians from monte carlo for stellar zenith
c
        dtau_ab_b(k) = dtau_ab(k)
c
c*TDR   asymmetry parameter
c
        if( dtau_sc(k) .gt. 0 ) then
          g(k) = gtot(k)
        else
          g(k) = 0.
        endif
c
c*TDR     single scattering albedo is one minus the single scattering 
c         co-albedo, now with Rayleigh portion removed
c
c          ssalb(k) = dtau_sc(k)/dtau_ex(k)
c
c*TDR     effective scattering coefficients (km**-1)
c
        a_sc(k) = dtau_sc(k)/dz(k)
c
1001  continue
c
c
c*TDR set derivative of index of refraction with respect to 
c     height. normalized to zero if irefract not equal to one.
c
      if( irefract .eq. 1 ) then
        call index_refract(ncomp,icomp,volmix,wl,
     -                     p(1),t(1),nr(1))
c
        do 2001 k=1,nlyr
c
          call index_refract(ncomp,icomp,volmix,wl,
     -                       p(k+1),t(k+1),nr(k+1))
c
          dlnnrdz(k) = 
     -        dble((alog(nr(k))-alog(nr(k+1)))/(alt(k)-alt(k+1)))
c
2001    continue
c
      else
c
        do 2051 k=1,nlyr
          dlnnrdz(k) = 0.d0
2051    continue
c
      endif
c
c*TDR loop over layers, seeing which need path calculations
c
      do 4001 k=1,nlyr
c
c*TDR   variable to keep track of whether or not we need
c       a path calculation
c
        do_path = 1
c
c*TDR   compare index of refraction at surface to see if the path 
c       calculation needs to be revisited.
c
        if( do_path .eq. 1 .and. irefract .eq. 1 ) then
          if( ABS(nr(nlev)-nr0)/nr0 .gt. tol ) then
            do_path = 0
            nr0 = nr(nlev)
          endif
        endif
c
c*TDR   compare the scattering optical depths and the asymmetry 
c       parameters to see if the path calculation needs to be 
c       revisited.
c
        if( do_path .eq. 1 ) then
          if( g0(k) .gt. 0 ) then
            if( ABS(g(k)-g0(k))/g0(k) .gt. tol ) then
              do_path = 0
            endif
          else
            if( g(k) .gt. 0 ) do_path = 0
          endif
c
          if( dtau_sc(k) .gt. 1.e-5 ) then
            if( dtau_sc0(k) .gt. 0 ) then
              if(ABS(dtau_sc(k)-dtau_sc0(k))/dtau_sc0(k) 
     -           .gt. tol) then
                do_path = 0
              endif
            else
              if( dtau_sc(k) .gt. 0 ) do_path = 0
            endif
          endif
c
c*TDR     store new scattering optical depths and asymmetry parameters,
c         if the path is being re-computed
c
          if( do_path .eq. 0 ) then
            g0(k) = g(k)
            dtau_sc0(k) = dtau_sc(k)
          endif
        endif
c
c*TDR   compare absorption and scattering optical depths along each 
c       tangent ray to decide which layers actually need the monte 
c       carlo calculation.
c
        tau_ab = 0.d0
        tau_sc = 0.d0
c
        do 2101 l=1,k
          tau_ab = tau_ab + dtau_ab(l)*pathg(k,l)
          tau_sc = tau_sc + dtau_sc(l)*pathg(k,l)
2101    continue
c
c*TDR   straight-line scattering transmission, which we will use later 
c       to decide how many photons are needed in the monte carlo 
c       simulation
c
        if( tau_sc .lt. 20. ) then
          trans = exp(-tau_sc)
        else
          trans = 0.
        endif
c
c*TDR   check if scattering optical depth is large
c
        if( tau_sc .gt. 1.d-5*tau_ab ) then
          do_mc(k) = 1
        else
          do_mc(k) = 0
        endif
c
c*TDR   if scattering calculation is needed, but we don't
c       yet have it stored away.
c
        if( do_mc(k) .eq. 1 .and. mc(k) .eq. 0 ) then
          do_path = 0
        else
          do_path = 1
        endif
c
c*TDR   do the path calculations, if warranted
c
        if( do_path .eq. 0 ) then
c
          mc(k) = 1
c
c*TDR     following de Kok and Stam (2012), the number of photons 
c         used in the monte carlo simulation is either 
c         ~nmin/(transmission) or nmax, whichever is smaller. note
c         that dK&S use nmin=10**5 and nmax=mcmx.
c
          if( trans .gt. real(nmin)/2.e9 ) then
            nmc = min(abs(int(nmin/trans)),nmax)
          else
            nmc = nmax
          endif
c
c*TDR     store number of photons in this simulation
c
c          write(*,*) 'mc:', wn, k
          write(iupath(k)) nmc
c
          do 3001 n=1,nmc
c
c*TDR       do monte carlo path for this photon
c
            call path_mc(iseed,nlev,irefract,Ntheta,Rs,b,au,
     -                   radius,dlnnrdz,a_sc,g,alt,zbar(k),
     -                   mu0,path0)
c
c*TDR       write path distribution and associated information to 
c           scratch file
c
            write(iupath(k)) (path0(l),l=1,nlev),
     -                       (mu0(ith),ith=1,Ntheta)
c
3001      continue
c
          rewind(iupath(k))
c
        endif
c
c*TDR   use monte carlo paths to determine transmission, zenith angle, and 
c       their Jacobians
c
        if( do_mc(k) .eq. 1 ) then
c
c*TDR     initialize before averaging over photons
c
          do 3205 ith=1,Ntheta
            dim(ith) = 0.0d0
            mu(ith)  = 0.0d0
            nd(ith)  = 0
            do 3201 l=1,nlyr
              ddim_dtau(l,ith) = 0.0d0
              dmu_dtau(l,ith)  = 0.0d0
3201        continue
3205      continue
c
c*TDR     average over photons
c
          read(iupath(k)) np
c
          do 3401 n=1,np
c
            read(iupath(k)) (path0(l),l=1,nlev),
     -                      (mu0(ith),ith=1,Ntheta)
c
c*TDR       path0(nlev) will be > 1 if photon 
c           struck surface, which is a case we ignore
c
            if( path0(nlev) .lt. 0.5 ) then
c
              trans = 1.d0
              do 3305 l=1,nlyr
                trans = trans*EXP(-dtau_ab(l)*path0(l))
3305          continue
c
              do 3351 ith=1,Ntheta
                if(mu0(ith) .ge. 0.d0 ) then
                  dim(ith) = dim(ith) + trans
                  nd(ith)  =  nd(ith) + 1
                  mu(ith)  =  mu(ith) + mu0(ith)*trans
c
                  do 3333 l=1,nlyr
                    ddim_dtau(l,ith) = ddim_dtau(l,ith) - 
     -                                 path0(l)*trans
                    dmu_dtau(l,ith)  = dmu_dtau(l,ith) - 
     -                                 mu0(ith)*path0(l)*trans
3333              continue
c
                endif
3351          continue
c
            endif
c
3401      continue
c
          do 3501 ith=1,Ntheta
            if( dim(ith) .gt. 1.e-6 .and. nd(ith) .gt. 0 ) then
              mu_b(k,ith) = mu(ith)/dim(ith)
              do 3451 l=1,nlyr
                dmu_dtau_b(k,l,ith) = (dim(ith)*dmu_dtau(l,ith) - 
     -                           mu(ith)*ddim_dtau(l,ith))/(dim(ith)**2)
3451          continue
            else
              mu_b(k,ith) = -1.
              do 3455 l=1,nlyr
                dmu_dtau_b(k,l,ith) = 0.0d0
3455          continue
            endif
c
c*TDR       store averages and Jacobians from monte carlo for transmission
c
            trans_b(k,ith) = dim(ith)/DBLE(np)
            do 3475 l=1,nlyr
              dtrans_dtau_b(k,l,ith) = ddim_dtau(l,ith)/DBLE(np)
3475        continue
c
3501      continue
c
          rewind(iupath(k))
c
        endif
c
4001  continue
c
      return
      end