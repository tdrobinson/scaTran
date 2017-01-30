      program scatran
c
      implicit none
      include 'param.inc'
c
c****   variable declarations
c
      integer ipi0,k
      integer nlev,iupath(kp),itrnst,irefract,ilimbd,iulimbd
      integer ncomp,icomp(6),do_mc(kp),Ntheta
      real p(kp),t(kp),alt(kp),grav(kp),volmix(6),ratm
      real au,radius,bimpact,radstar,logg,metal,teff
      real nr0,g0(kp),dtau_sc0(kp)
      real copi(kp),dtau_ex(kp),g(kp)
      real copi_i(kp),dtau_ex_i(kp)
      real tau_ab(kp),tau_sc(kp)
      real zeff,tdepth
      real tp,gravp,Hp,dtau_a,dtau_s,tautot
      real pmin,pmax,dlnp,altm
      double precision wn,wn_i
      double precision area(kp,mxth),areag(kp)
      double precision path(kp,kp),pathg(kp,kp)
      double precision mu(kp,mxth)
      double precision trans_b(kp,mxth),mu_b(kp,mxth)
      double precision dtrans_dtau_b(kp,kp,mxth)
      double precision dmu_dtau_b(kp,kp,mxth),dtau_ab_b(kp)
c
c****   set up model atmosphere and conditions
c
        nlev      = 126    ! number of levels
        itrnst    = 3      ! (1) geometric, (2) ray trace, or (3) monte carlo
        irefract  = 0      ! (0) no refraction, (1) with refraction
        ilimbd    = 0      ! (0) no limb darkening, (1) limb darkening
        
        ncomp     = 2      ! number of major background gases
        icomp(1)  = 5      ! (1) air (2) co2 (3) n2 (4) o2 (5) h2 (6) he
        icomp(2)  = 6
        volmix(1) = 0.85   ! volume mixing ratio of H2
        volmix(2) = 0.15   ! volume mixing ratio of He
        ratm      = 3.59e3 ! specific gas constant (J/K/kmole)

        radstar   = 0.78   ! stellar radius (in solar radii)
        au        = 0.031  ! orbital distance (au)
        bimpact   = 0.     ! impact parameter (in stellar radii)

        radius    = 81237. ! planet radius (km)
        tp        = 1500.  ! isothermal atmosphere temperature (K)
        gravp     = 21.8   ! planet surface gravity (m/s/s)

        pmin      = 1.e-4  ! minimum profile pressure (Pa)
        pmax      = 1.e6   ! maximum profile pressure
        dlnp      = (LOG(pmax)-LOG(pmin))/(REAL(nlev)-1.)

c
c****   do the hydrostatic calculation
c
        alt(nlev)  = 0.
        p(nlev)    = pmax
        t(nlev)    = tp
        grav(nlev) = gravp
        Hp         = ratm*tp/gravp/1.e3 ! scale height (km)
        do k=nlev-1,1,-1
          p(k)       = EXP(ALOG(p(k+1)) - dlnp)
          t(k)       = tp
          altm       = alt(k+1) + 
     -                 Hp*gravp/grav(k+1)*LOG(p(k+1)/0.5/(p(k+1)+p(k)))
          grav(k)    = gravp*radius**2/(radius + altm)**2
          alt(k)     = alt(k+1) + Hp*gravp/grav(k+1)*dlnp
          iupath(k)  = 100 + k          ! sets up scratch file I/O units
        enddo
c
c****   set up optical properties where Monte Carlo and Jacobians are computed
c
        wn = 8695.d0 ! wavenumber (cm**-1)
        tautot = 0.0
        do k=1,nlev-1
          dtau_a = 0.002*(p(k)/2.0e3) ! gives slant tau ~ 0.56 at 2e3 Pa
          dtau_s = 0.
          if( ABS(2.90e3-alt(k)) .lt. Hp/2. ) then
            dtau_s = 0.031*EXP(-ABS(2.90e3-alt(k))/Hp) ! distributes tau ~ 10 cloud around 10 Pa
          endif
          tautot = tautot + dtau_s
          dtau_ex(k) = dtau_a + dtau_s   ! extinction optical depth
          g(k)       = 0.95              ! asymmetry parameter
          copi(k)    = dtau_a/dtau_ex(k) ! single scattering co-albedo (i.e., 1 - single scattering albedo)
        enddo
c
c****   set up optical properties where transit depth is interpolated to, 
c       which can be the same as above
c
        wn_i = 8695.d0 ! wavenumber (cm**-1)
        do k=1,nlev-1
          dtau_ex_i(k) = dtau_ex(k)   ! extinction optical depth
          copi_i(k)    = copi(k)      ! single scattering co-albedo (i.e., 1 - single scattering albedo)
        enddo
c
c****   call routine to initialize variables and gray quantities
c
      call transit_init(nlev,itrnst,irefract,iupath,ncomp,icomp,
     -                  volmix,radstar,bimpact,radius,au,p,t,alt,
     -                  Ntheta,nr0,g0,dtau_sc0,areag,area,
     -                  pathg,path,mu)
c
c****   for monte carlo transit spectra, check if grey 
c       (i.e., geometric) calculation suffices by computing slant 
c       absorption and scattering optical depths
c
      ipi0 = 0
c
      if( itrnst .eq. 3 ) then
c
        call transit_optd(nlev,copi,dtau_ex,pathg,tau_ab,tau_sc)
c
        do 2001 k=1,nlev-1
          if( tau_sc(k) .gt. 1.e-5*tau_ab(k) ) ipi0 = 1
2001    continue
c
      endif
c
c****   if a scattering calculation is warranted, do the 
c       monte carlo simulation and write to scratch
c
      if( ipi0 .eq. 1 ) then
c
        call transit_scat(nlev,iupath,irefract,itrnst,
     -                    ncomp,icomp,volmix,radstar,bimpact,radius,au,
     -                    p,t,alt,grav,ratm,wn,copi,dtau_ex,
     -                    g,pathg,Ntheta,do_mc,trans_b,mu_b,
     -                    dtau_ab_b,dtrans_dtau_b,dmu_dtau_b)
c
      endif
c
c****   compute transit depth, interpolated to wni, with optical 
c       properties copi_i and dtau_ex_i.
c

      call transit_interp(ipi0,nlev,itrnst,ilimbd,irefract,iulimbd,
     -                    ncomp,icomp,volmix,nr0,p,t,alt,grav,
     -                    ratm,au,radius,radstar,logg,metal,teff,
     -                    bimpact,Ntheta,wn_i,copi_i,dtau_ex_i,pathg,
     -                    path,areag,area,do_mc,mu,trans_b,mu_b,
     -                    dtau_ab_b,dtrans_dtau_b,dmu_dtau_b,zeff,
     -                    tdepth)
c
c****   print output
c
      write(*,*) 'Transit depth: ', tdepth
c
      stop
      end