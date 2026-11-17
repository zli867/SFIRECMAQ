
      SUBROUTINE plumeRiseFreitas(m1,area_in,
     &                    up,vp,ta,pres,rv,zt_rams,zm_rams,
     &                    bflux_in,zplm,plm_bottom,plm_top,VMD_out)

      implicit none
!      REAL,PARAMETER :: g = 9.81  ! acceleration due to gravity (m {s}^-2)
      REAL,PARAMETER :: r_d= 287.       ! gas constant of dry air (Jdeg^-1 kg^-1)
      REAL,PARAMETER :: cp= 1004.6     ! Specific heat of dry air at constant pressure (J deg^-1 kg^-1)
      REAL,PARAMETER :: rcp= r_d/cp
!      real,parameter :: rgas=r_d
!      real,parameter :: cpor=1./rcp
      REAL,PARAMETER :: p1000mb=100000.    ! pressure at 1000 hPa (pa)
!      real,parameter :: p00=p1000mb
!  real, external:: esat_pr


      INTEGER, PARAMETER ::  iflam_frac=1, imean_frp=2, istd_frp=3,
     &                        imean_size=4, istd_size=5       ! RAR:

!--for cmaq
      real, intent (in) :: area_in, bflux_in
      INTEGER, intent (in) :: m1
      real, intent (in), dimension(m1) :: ta, pres, up, vp, rv,
     &                                     zt_rams, zm_rams
      real, intent (out) :: zplm, plm_bottom, plm_top
      real, intent (out), dimension(m1,2) :: VMD_out
!-fre
      real, dimension(5)  :: plume_fre


      integer :: i,j,k,iveg_ag,imm,ixx

      INTEGER ::  k1,k2

      integer :: ncall
      integer :: kmt

      real,    dimension(m1) :: theta,pp
      real                :: burnt_area,STD_burnt_area,dzi,FRP
      real,    dimension(2)  :: ztopmax


! From plumerise1.F routine
      integer, parameter :: nveg_agreg      = 4
      integer, parameter :: tropical_forest = 1
      integer, parameter :: boreal_forest   = 2
      integer, parameter :: savannah        = 3

      integer, parameter :: grassland       = 4
      real, dimension(nveg_agreg) :: firesize

      INTEGER, PARAMETER :: wind_eff = 1

      INTEGER, SAVE :: icall
      real   :: convert_smold_to_flam

      integer, parameter :: nkp = 200, ntime = 200
      real,dimension(nkp) ::  w,t,qv,qc,qh,qi,sc,    ! blob
     &                         vth,vti,rho,txs,
     &                         est,qsat,qpas,qtotal

      real,dimension(nkp) ::  wc,wt,tt,qvt,qct,qht,qit,sct
      real,dimension(nkp) ::  dzm,dzt,zm,zt,vctr1,vctr2
     &                        ,vt3dc,vt3df,vt3dk,vt3dg,scr1

      real,dimension(nkp) ::  pke,the,thve,thee,pe,
     &                        te,qvenv,rhe,dne,sce ! environment at plume grid
      real,dimension(nkp) ::  ucon,vcon,wcon,thtcon ,rvcon,picon
     &                        ,tmpcon,dncon,prcon
     &                        ,zcon,zzcon,scon,VISC ! environment at RAMS  grid
     &                        ,CVH,CVI,RADIUS
     &                        ,upe,vpe,vel_e,vel_p,rad_p,vel_t,rad_t

      real :: DZ,DQSDZ,VISCOSITY,TSTPF
      integer :: N,NM1,L
      real :: ADVW,ADVT,ADVV,ADVC,ADVH,ADVI,ADIABAT,
     &         WBAR,ALAST(10),VHREL,VIREL  ! advection
      real :: ZSURF,ZBASE,ZTOP
      real :: AREA,RSURF,ALPHA,FMOIST,BLOAD
      real,dimension(ntime) :: HEATING
      real :: DT,TIME,TDUR
      integer :: LBASE,MINTIME,MDUR,MAXTIME
      REAL,    DIMENSION(nkp,2)    :: W_VMD,VMD
      real :: ztop_(ntime)
!      print *,"plumeRiseFreitas"
! pre calculation
      do k=1,m1
        pp(k)=cp*(pres(k)/p1000mb)**rcp
        theta(k)=ta(k)/pp(k)*cp
      end do
      do i=1,nveg_agreg
        firesize(i)=area_in*4046.86  ! acres to m2
      end do
      plume_fre(1)=1.0
      plume_fre(2)=bflux_in
      plume_fre(3)=0.0
      plume_fre(4)=area_in*4046.86  ! acres to m2
      plume_fre(5)=0.0
!Fator de conversao de unidades
!!fcu=1.      !=> kg [gas/part] /kg [ar]
!!fcu =1.e+12   !=> ng [gas/part] /kg [ar]
!!real,parameter :: fcu =1.e+6 !=> mg [gas/part] /kg [ar]
!----------------------------------------------------------------------
!               indexacao para o array "plume(k,i,j)"
! k
! 1   => area media (m^2) dos focos  em biomas floresta dentro do gribox i,j
! 2   => area media (m^2) dos focos  em biomas savana dentro do gribox i,j
! 3   => area media (m^2) dos focos  em biomas pastagem dentro do gribox i,j
! 4   => desvio padrao da area media (m^2) dos focos : floresta
! 5   => desvio padrao da area media (m^2) dos focos : savana
! 6   => desvio padrao da area media (m^2) dos focos : pastagem
! 7 a 9 =>  sem uso
!10(=k_CO_smold) => parte da emissao total de CO correspondente a fase smoldering
!11, 12 e 13 => este array guarda a relacao entre
!               qCO( flaming, floresta) e a quantidade total emitida
!               na fase smoldering, isto e;
!               qCO( flaming, floresta) =  plume(11,i,j)*plume(10,i,j)
!               qCO( flaming, savana  ) =  plume(12,i,j)*plume(10,i,j)
!               qCO( flaming, pastagem) =  plume(13,i,j)*plume(10,i,j)
!20(=k_PM25_smold),21,22 e 23 o mesmo para PM25
!
!24-n1 =>  sem uso
!----------------------------------------------------------------------
      ncall=0
      if (ncall == 0) then
        ncall = 1
        w=0.0;t=0.0;qv=0.0;qc=0.0;qh=0.0;qi=0.0;sc=0.0
        vth=0.0;vti=0.0;rho=0.0;txs=0.0
        est=0.0;qsat=0.0;qpas=0.0;qtotal=0.0
        wc=0.0;wt=0.0;tt=0.0;qvt=0.0;qct=0.0;qht=0.0;qit=0.0;sct=0.0
        dzm=0.0;dzt=0.0;zm=0.0;zt=0.0;vctr1=0.0;vctr2=0.0
        vt3dc=0.0;vt3df=0.0;vt3dk=0.0;vt3dg=0.0;scr1=0.0
        pke=0.0;the=0.0;thve=0.0;thee=0.0;pe=0.0;te=0.0
        qvenv=0.0;rhe=0.0;dne=0.0;sce=0.0
        ucon=0.0;vcon=0.0;wcon=0.0;thtcon =0.0;rvcon=0.0
        picon=0.0;tmpcon=0.0;dncon=0.0;prcon=0.0
        zcon=0.0;zzcon=0.0;scon=0.0
        dz=0.0;dqsdz=0.0;visc=0.0;viscosity=0.0;tstpf=0.0
        advw=0.0;advt=0.0;advv=0.0;advc=0.0;advh=0.0;advi=0.0
        cvh=0.0;cvi=0.0;adiabat=0.0
        wbar=0.0;alast=0.0;vhrel=0.0;virel=0.0
        zsurf=0.0;zbase=0.0;ztop=0.0;area=0.0;rsurf=0.0
        alpha=0.0;radius=0.0;heating=0.0
        fmoist=0.0;bload=0.0;dt=0.0;time=0.0;tdur=0.0
        ztop_=0.0
        upe =0.0
        vpe =0.0
        vel_e =0.0
        vel_p =0.0
        rad_p =0.0
        vel_t =0.0
        rad_t =0.0
        W_VMD=0.0
        VMD=0.0
        n=0;nm1=0;l=0;lbase=0;mintime=0;mdur=0;maxtime=0
      endif
      j=1
      i=1
!- if the max value of flaming is close to zero => there is not emission with
!- plume rise => cycle
!later add calculate theta,pi
      do k = 1,m1 ! loop over vertical grid
         ucon  (k)=up(k)          ! u wind
         vcon  (k)=vp(k)          ! v wind
         thtcon(k)=theta(k)       ! pot temperature
         picon (k)=pp(k)                ! exner function
         rvcon (k)=rv(k)          ! water vapor mixing ratio
         zcon  (k)=zt_rams(k)               ! termod-point height
         zzcon (k)=zm_rams(k)               ! W-point height
      enddo
!- get envinronmental state (temp, water vapor mix ratio, ...)
      call get_env_condition(1,m1,kmt,wind_eff,
     &                       qvenv,zcon,zt,ucon,upe,pke,
     &                       vcon,vpe,thtcon,the,rvcon,picon,
     &                       thve,pe,dne,vel_e,te,dz,dzm,dzt,zm,
     &                       zsurf)
!      print *,'get_env_condition zzcon',zzcon
!-from the amount emitted during the smoldering phase
      convert_smold_to_flam=plume_fre(iflam_frac)

!- check if there is only one fire in a given grid box (=> std =0.)
      if(plume_fre(istd_frp ) < 1.0e-6) then

!- if yes, we will set it as a 20% of the mean frp as a gross estimation
!- of the retrieval uncertainty by the sensors.
!- (we are not taking care about the fire size retrieval)
         plume_fre(istd_frp )=0.2*plume_fre(imean_frp )
      endif

!- loop over the minimum and maximum heat fluxes/FRP
      lp_minmax: do imm=1,2
        if(imm==1 ) then
!for imm = 1 => lower injection height
          burnt_area = max(1.0e4,plume_fre(imean_size) -
     &                      0.5*plume_fre(istd_size))
          FRP = max(1000.,plume_fre(imean_frp ) -
     &               0.5*plume_fre(istd_frp ))
        elseif(imm==2 ) then
!for imm = 2 => higher injection height
          burnt_area = max(1.0e4,plume_fre(imean_size) +
     &                      0.5*plume_fre(istd_size))
          FRP = max(1000.,plume_fre(imean_frp ) +
     &               0.5*plume_fre(istd_frp ))
        endif
!- get fire properties (burned area, plume radius, heating rates ...)
        call get_fire_properties(imm,iveg_ag,burnt_area,
     &                            STD_burnt_area,FRP,
     &                                AREA,BLOAD,alpha,RSURF,
     &                       FMOIST,TDUR,MDUR,maxtime,HEATING)
!------  generates the plume rise    ------
        call makeplume (kmt,ztopmax(imm),ixx,imm,
     &        tstpf,viscosity,ztop,TIME,
     &        dt,wbar,adiabat,alpha,DQSDZ,
     &        L,nm1,mintime,maxtime,
     &        zm,wc,wt,dzm,w,tt,t,te,
     &        qvt,qv,qvenv,qct,qc,qht,qh,qit,qi,
     &        vel_e,vel_p,vel_t,rad_p,rad_t,pe,zt,
     &        qsat,est,SCR1,txs,
     &        ztop_,W_VMD,Visc,vth,vti,RSURF,tdur,
     &        FMOIST,heating,cvi,area,sc,cvh,VHREL,VIREL,
     &        Vt3dc,vt3df,vt3dk,vt3dg,
     &        vctr1,vctr2,dzt,sct)
      enddo lp_minmax

!- define o dominio vertical onde a emissao flaming ira ser colocada
      call set_flam_vert(ztopmax,k1,k2,nkp,zzcon,W_VMD,VMD)
!- thickness of the vertical layer between k1 and k2 eta levels (lower and upper bounds for the injection height )
      dzi= 1./(zzcon(k2)-zzcon(k1-1))
      zplm = 0.5*(zzcon(k2)+zzcon(k1-1))
      plm_bottom = zzcon(k1-1)
      plm_top = zzcon(k2)
          !-- map plume-grid VMD to RAMS vertical grid
      VMD_out(:,:) = 0.0

      do imm = 1, 2
         do k = 1, m1
           VMD_out(k,imm) = VMD(k,imm)
         enddo
      enddo
      end subroutine plumeRiseFreitas
!-------------------------------------------------------------------------

      subroutine get_env_condition(k1,k2,kmt,wind_eff,
     &                       qvenv,zcon,zt,ucon,upe,pke,
     &                       vcon,vpe,thtcon,the,rvcon,picon,
     &                       thve,pe,dne,vel_e,te,dz,dzm,dzt,zm,
     &                       ZSURF)

!use module_zero_plumegen_coms
!use rconstants
      implicit none
      integer :: k1,k2,k,kmt,nk
      real :: znz
      integer :: n_setgrid
      integer :: wind_eff
      integer, parameter :: nkp = 200
      real,dimension(nkp) :: qvenv,zcon,zt,ucon,upe,pke,
     &                       vcon,vpe,thtcon,the,rvcon,picon,
     &                       thve,pe,dne,vel_e,te
      REAL,PARAMETER :: g = 9.81  ! acceleration due to gravity (m {s}^-2)
      REAL,PARAMETER :: r_d= 287.       ! gas constant of dry air (Jdeg^-1 kg^-1)
      REAL,PARAMETER :: cp= 1004.6     ! Specific heat of dry air at constant pressure (J deg^-1 kg^-1)
      REAL,PARAMETER :: rcp= r_d/cp
      real,parameter :: rgas=r_d
      real,parameter :: cpor=1./rcp
      REAL,PARAMETER :: p1000mb=100000.    ! pressure at 1000 hPa (pa)
      real,parameter :: p00=p1000mb
      real :: dz,ZSURF
      real,dimension(nkp) :: zm,dzm,dzt
!      print *,"get_env_condition",n_setgrid
      n_setgrid =0
      if( n_setgrid == 0) then
        n_setgrid = 1
        call set_grid(dz,zsurf,zt,dzt,zm,dzm) ! define vertical grid of plume model
                      ! zt(k) =  thermo and water levels
                      ! zm(k) =  dynamical levels
      endif

      znz=zcon(k2)
      do k=nkp,1,-1
        if(zt(k).lt.znz)go to 13
      enddo
      stop ' envir stop 12'
13    continue
!-srf-mb
      kmt=min(k,nkp-1)

      nk=k2-k1+1
      call htint(nk,  ucon,zcon,kmt,upe,zt)
      call htint(nk,  vcon,zcon,kmt,vpe,zt)
      call htint(nk,thtcon,zcon,kmt,the  ,zt)
      call htint(nk, rvcon,zcon,kmt,qvenv,zt)
      do k=1,kmt
        qvenv(k)=max(qvenv(k),1e-8)
      enddo

      pke(1)=picon(1)
      do k=1,kmt
        thve(k)=the(k)*(1.+.61*qvenv(k)) ! virtual pot temperature
      enddo
      do k=2,kmt
        pke(k)=pke(k-1)-g*2.*(zt(k)-zt(k-1))
     &         /(thve(k)+thve(k-1))
      enddo
      do k=1,kmt
        te(k)  = the(k)*pke(k)/cp         ! temperature (K)
        pe(k)  = (pke(k)/cp)**cpor*p00    ! pressure (Pa)
        dne(k)= pe(k)/(rgas*te(k)*(1.+.61*qvenv(k))) !  dry air density (kg/m3)
        vel_e(k) = sqrt(upe(k)**2+vpe(k)**2)         !-env wind (m/s)
      enddo

!-ewe - env wind effect
      if(wind_eff < 1)  then
        vel_e(1:kmt) = 0.
      endif
!--------- converte press de Pa para kPa para uso modelo de plumerise
      do k=1,kmt
        pe(k) = pe(k)*1.e-3
      enddo
      return
      end subroutine get_env_condition

!-------------------------------------------------------------------------

      subroutine set_grid(dz,zsurf,zt,dzt,zm,dzm)
      implicit none
      integer, parameter :: nkp = 200
      integer :: k,mzp
      real :: DZ,ZSURF
      real,dimension(nkp) :: zt,zm,dzm,dzt
!      print *,"set_grid"
      dz=100. ! set constant grid spacing of plume grid model(meters)

      mzp=nkp
      zt(1) = zsurf
      zm(1) = zsurf
      zt(2) = zt(1) + 0.5*dz
      zm(2) = zm(1) + dz
      do k=3,mzp
        zt(k) = zt(k-1) + dz ! thermo and water levels
        zm(k) = zm(k-1) + dz ! dynamical levels
      enddo
      do k = 1,mzp-1
        dzm(k) = 1. / (zt(k+1) - zt(k))
      enddo
      dzm(mzp)=dzm(mzp-1)

      do k = 2,mzp
        dzt(k) = 1. / (zm(k) - zm(k-1))
      enddo
      dzt(1) = dzt(2) * dzt(2) / dzt(3)

!   dzm(1) = 0.5/dz
!   dzm(2:mzp) = 1./dz
      return
      end subroutine set_grid
!-------------------------------------------------------------------------

      SUBROUTINE set_flam_vert(ztopmax,k1,k2,nkp,zzcon,W_VMD,VMD)
      implicit none
      REAL    , INTENT(IN)  :: ztopmax(2)
      INTEGER , INTENT(OUT) :: k1
      INTEGER , INTENT(OUT) :: k2

! plumegen_coms
      INTEGER , INTENT(IN)  :: nkp
      REAL    , INTENT(IN)  :: zzcon(nkp)

      INTEGER imm,k
      INTEGER, DIMENSION(2)  :: k_lim

!- version 2
      REAL    , INTENT(IN)  :: W_VMD(nkp,2)
      REAL    , INTENT(OUT) ::   VMD(nkp,2)
      real   w_thresold,xxx
      integer k_initial,k_final,ko,kk4,kl
!      print *,'set_flam_vert'
!- version 1
      DO imm=1,2
         ! checar
         !    do k=1,m1-1
        DO k=1,nkp-1
          IF(zzcon(k) > ztopmax(imm) ) EXIT
        ENDDO
        k_lim(imm) = k
      ENDDO
      k1=MAX(3,k_lim(1))
      k2=MAX(3,k_lim(2))

      IF(k2 < k1) THEN
        k2=k1
        !stop 1234
      ENDIF

      !- version 2
      !- vertical mass distribution
      !-
      w_thresold = 1.
      DO imm=1,2
        VMD(1:nkp,imm)= 0.
        xxx=0.
        k_initial= 0
        k_final  = 0

       !- define range of the upper detrainemnt layer
        do ko=nkp-10,2,-1
          if(w_vmd(ko,imm) < w_thresold) cycle
          if(k_final==0) then
             k_final=ko
          endif
          if(w_vmd(ko,imm)-1. > w_vmd(ko-1,imm)) then
            k_initial=ko
          exit
          endif

        enddo
       !- if there is a non zero depth layer, make the mass vertical distribution
        if(k_final > 0 .and. k_initial > 0) then
          k_initial=int((k_final+k_initial)*0.5)

        !- parabolic vertical distribution between k_initial and k_final
          kk4 = k_final-k_initial+2
          do ko=1,kk4-1
            kl=ko+k_initial-1
            VMD(kl,imm) = 6.* float(ko)/float(kk4)**2 *
     &                    (1. - float(ko)/float(kk4))
          enddo
          if(sum(VMD(1:NKP,imm)) .ne. 1.) then
            xxx= ( 1.- sum(VMD(1:NKP,imm)) )/float(k_final-k_initial+1)
            do ko=k_initial,k_final
              VMD(ko,imm) = VMD(ko,imm)+ xxx !- values between 0 and 1.
            enddo
           ! print*,'new mass=',sum(mass)*100.,xxx
           !pause
          endif
        endif !k_final > 0 .and. k_initial >

      ENDDO

      END SUBROUTINE set_flam_vert
!-------------------------------------------------------------------------

      subroutine get_fire_properties(imm,iveg_ag,burnt_area,
     &                                STD_burnt_area,FRP,
     &                                AREA,BLOAD,alpha,RSURF,
     &                       FMOIST,TDUR,MDUR,maxtime,HEATING)
      implicit none
      integer ::  moist,  i,  icount,imm,iveg_ag
      real::   bfract,  effload,  heat,  hinc ,burnt_area,
     &          STD_burnt_area,heat_fluxW,FRP
      real,    dimension(2,4) :: heat_flux
      INTEGER, parameter :: use_last = 0
!real, parameter :: beta = 5.0   !ref.: Wooster et al., 2005
      REAL, parameter :: beta = 0.88  !ref.: Paugam et al., 2015

      real :: AREA,BLOAD,alpha,RSURF, FMOIST,TDUR
      integer :: MDUR,maxtime
      integer, parameter :: ntime = 200
      real,    dimension(ntime) :: HEATING

      data heat_flux/
!---------------------------------------------------------------------
!  heat flux      !IGBP Land Cover      !
! min  ! max      !Legend and           ! reference
!    kW/m^2       !description          !
!--------------------------------------------------------------------
     & 30.0,     80.0,   ! Tropical Forest         ! igbp 2 & 4
     & 30.0,   80.0,      ! Boreal forest           ! igbp 1 & 3
     & 4.4,  23.0,   ! cerrado/woody savanna   | igbp  5 thru 9
     & 3.3,   3.3    /! Grassland/cropland      ! igbp 10 thru 17
!--------------------------------------------------------------------
!-- fire at surface
!
!area = 20.e+4   ! area of burn, m^2
!      print *, 'get_fire_properties'
      area = burnt_area! area of burn, m^2

! "beta" factor converts FRP to convective energy
      heat_fluxW = beta*(FRP/area)/0.55 ! in W/m^2
      ! Zongrun
      PRINT '(A,F10.2)', 'heat_fluxW = ', FRP / area

      mdur = 53        ! duration of burn, minutes
      bload = 10.      ! total loading, kg/m**2
      moist = 10       ! fuel moisture, %. average fuel moisture,percent dry
      maxtime =mdur+2  ! model time, min
!heat = 21.e6    !- joules per kg of fuel consumed
!heat = 15.5e6   !joules/kg - cerrado
      heat = 19.3e6    !joules/kg - floresta em alta floresta (mt)
!alpha = 0.1      !- entrainment constant
      alpha = 0.05      !- entrainment constant

!IF (MOD (MAXTIME, 2) .NE.0) MAXTIME = MAXTIME+1  !make maxtime even

      MAXTIME = MAXTIME * 60  ! and put in seconds
!
      RSURF = SQRT (AREA / 3.14159) !- entrainment surface radius (m)

      FMOIST   = MOIST / 100.       !- fuel moisture fraction
!
!
! calculate the energy flux and water content at lboundary.
! fills heating() on a minute basis. could ask for a file at this po
! in the program. whatever is input has to be adjusted to a one
! minute timescale.
!

      DO I = 1, ntime         !- make sure of energy release
        HEATING (I) = 0.0001  !- avoid possible divide by 0
      enddo
!
      TDUR = MDUR * 60.       !- number of seconds in the burn

      bfract = 1.             !- combustion factor

      EFFLOAD = BLOAD * BFRACT  !- patchy burning

!     spread the burning evenly over the interval
!     except for the first few minutes for stability
      ICOUNT = 1
!
      if(MDUR > NTIME) then
        STOP
        print *,'Increase time duration (ntime) in min'
      endif

      DO WHILE (ICOUNT.LE.MDUR)
!  HEATING (ICOUNT) = HEAT * EFFLOAD / TDUR  ! W/m**2
!  HEATING (ICOUNT) = 80000.  * 0.55         ! W/m**2

        HEATING (ICOUNT) = heat_fluxW  * 0.55     ! W/m**2 (0.55 converte para energia convectiva)
        ICOUNT = ICOUNT + 1
      ENDDO
!     ramp for 5 minutes
      IF(use_last /= 1) THEN

        HINC = HEATING (1) / 4.
        HEATING (1) = 0.1
        HEATING (2) = HINC
        HEATING (3) = 2. * HINC
        HEATING (4) = 3. * HINC
      ELSE
        IF(imm==1) THEN
          HINC = HEATING (1) / 4.
          HEATING (1) = 0.1
          HEATING (2) = HINC
          HEATING (3) = 2. * HINC
          HEATING (4) = 3. * HINC
        ELSE
          HINC = (HEATING (1) - heat_flux(imm-1,iveg_ag)
     &             * 1000. *0.55)/ 4.
          HEATING (1) = heat_flux(imm-1,iveg_ag) * 1000. *0.55 + 0.1
          HEATING (2) = HEATING (1)+ HINC
          HEATING (3) = HEATING (2)+ HINC
          HEATING (4) = HEATING (3)+ HINC
        ENDIF
      ENDIF

      return
      end subroutine get_fire_properties
!-------------------------------------------------------------------------------
!
      SUBROUTINE MAKEPLUME ( kmt,ztopmax,ixx,imm,
     &        tstpf,viscosity,ztop,TIME,
     &        dt,wbar,adiabat,alpha,DQSDZ,
     &        L,nm1,mintime,maxtime,
     &        zm,wc,wt,dzm,w,tt,t,te,
     &        qvt,qv,qvenv,qct,qc,qht,qh,qit,qi,
     &        vel_e,vel_p,vel_t,rad_p,rad_t,pe,zt,
     &        qsat,est,SCR1,txs,
     &        ztop_,W_VMD,Visc,vth,vti,RSURF,tdur,
     &        FMOIST,heating,cvi,area,sc,cvh,VHREL,VIREL,
     &        Vt3dc,vt3df,vt3dk,vt3dg,
     &        vctr1,vctr2,dzt,sct)
!
! *********************************************************************
!
!    EQUATION SOURCE--Kessler Met.Monograph No. 32 V.10 (K)
!    Alan Weinstein, JAS V.27 pp 246-255. (W),
!    Ogura and Takahashi, Monthly Weather Review V.99,pp895-911 (OT)
!    Roger Pielke,Mesoscale Meteorological Modeling,Academic Press,1984
!    Originally developed by: Don Latham (USFS)
!
!
! ************************ VARIABLE ID ********************************
!
!     DT=COMPUTING TIME INCREMENT (SEC)
!     DZ=VERTICAL INCREMENT (M)
!     LBASE=LEVEL ,CLOUD BASE
!
!     CONSTANTS:
!       G = GRAVITATIONAL ACCELERATION 9.80796 (M/SEC/SEC).
!       R = DRY AIR GAS CONSTANT (287.04E6 JOULE/KG/DEG K)
!       CP = SPECIFIC HT. (1004 JOULE/KG/DEG K)
!       HEATCOND = HEAT OF CONDENSATION (2.5E6 JOULE/KG)
!       HEATFUS = HEAT OF FUSION (3.336E5 JOULE/KG)
!       HEATSUBL = HEAT OF SUBLIMATION (2.83396E6 JOULE/KG)
!       EPS = RATIO OF MOL.WT. OF WATER VAPOR TO THAT OF DRY AIR (0.622)
!       DES = DIFFERENCE BETWEEN VAPOR PRESSURE OVER WATER AND ICE (MB)
!       TFREEZE = FREEZING TEMPERATURE (K)
!
!
!     PARCEL VALUES:
!       T = TEMPERATURE (K)
!       TXS = TEMPERATURE EXCESS (K)
!       QH = HYDROMETEOR WATER CONTENT (G/G DRY AIR)
!       QHI = HYDROMETEOR ICE CONTENT (G/G DRY AIR)
!       QC = WATER CONTENT (G/G DRY AIR)
!       QVAP = WATER VAPOR MIXING RATIO (G/G DRY AIR)
!       QSAT = SATURATION MIXING RATIO (G/G DRY AIR)
!       RHO = DRY AIR DENSITY (G/M**3) MASSES = RHO*Q'S IN G/M**3
!       ES = SATURATION VAPOR PRESSURE (kPa)
!
!     ENVIRONMENT VALUES:
!       TE = TEMPERATURE (K)
!       PE = PRESSURE (kPa)
!       QVENV = WATER VAPOR (G/G)
!       RHE = RELATIVE HUMIDITY FRACTION (e/esat)
!       DNE = dry air density (kg/m^3)
!
!     HEAT VALUES:
!       HEATING = HEAT OUTPUT OF FIRE (WATTS/M**2)
!       MDUR = DURATION OF BURN, MINUTES
!
!       W = VERTICAL VELOCITY (M/S)
!       RADIUS=ENTRAINMENT RADIUS (FCN OF Z)
!   RSURF = ENTRAINMENT RADIUS AT GROUND (SIMPLE PLUME, TURNER)
!   ALPHA = ENTRAINMENT CONSTANT
!       MAXTIME = TERMINATION TIME (MIN)
!
!
!**********************************************************************
!**********************************************************************
      implicit none
!logical :: endspace
      integer, parameter :: nkp = 200, ntime = 200
      character (len=10) :: varn
      integer ::  izprint, iconv,  itime, k, kk,
     &             kkmax, deltak,ilastprint,kmt
     &            ,ixx,nrectotal,i_micro,n_sub_step
      real ::  vc, g,  r,  cp,  eps,
     &          tmelt,  heatsubl,  heatfus,  heatcond, tfreeze,
     &          ztopmax, wmax, rmaxtime, es, esat, heat,dt_save !ESAT_PR,
      character (len=2) :: cixx
! Set threshold to be the same as dz=100., the constant grid spacing of plume grid model(meters) found in set_grid()
      REAL :: DELZ_THRESOLD = 100.

      INTEGER     :: imm

!  real, external:: esat_pr!
!
! ******************* SOME CONSTANTS **********************************
!
!      XNO=10.0E06 median volume diameter raindrop (K table 4)
!      VC = 38.3/(XNO**.125) mean volume fallspeed eqn. (K)
!
      parameter (vc = 5.107387)
      parameter (g = 9.80796, r = 287.04, cp = 1004.)
      parameter (eps = 0.622,  tmelt = 273.3)
      parameter (heatsubl = 2.834e6)
      parameter (heatfus = 3.34e5, heatcond = 2.501e6)
      parameter (tfreeze = 269.3)
!
      real,dimension(nkp) :: rho,radius

      REAL :: tstpf,viscosity,ztop,TIME
      REAL :: dt,wbar,adiabat,alpha,DQSDZ
      integer :: L,nm1,mintime,maxtime
      real,dimension(nkp) ::  zm,wc,wt,dzm,w,tt,t,te,
     &         qvt,qv,qvenv,qct,qc,qht,qh,qit,qi,
     &         vel_e,vel_p,vel_t,rad_p,rad_t,pe,zt,
     &         qsat,est,SCR1,txs
      real,    DIMENSION(ntime):: ztop_
      REAL,    DIMENSION(nkp,2)    :: W_VMD

      real,dimension(nkp) ::  cvi,cvh
      real,dimension(nkp) ::  VISC,vth,vti,sc
      real :: RSURF,TDUR,FMOIST,area
      real,dimension(ntime) :: HEATING
      real,dimension(nkp) ::  Vt3dc,vt3df,vt3dk,vt3dg,
     &                       vctr1,vctr2,dzt,sct
      real :: VHREL,VIREL
      REAL :: ESAT_PR
!      print *, 'MAKEPLUME'
      tstpf = 2.0   !- timestep factor
      viscosity = 500.!- viscosity constant (original value: 0.001)

      nrectotal=150
!
!*************** PROBLEM SETUP AND INITIAL CONDITIONS *****************
      mintime = 1
      ztopmax = 0.
      ztop    = 0.
         time = 0.
           dt = 1.
         wmax = 1.
      kkmax   = 10
      deltaK  = 20
      ilastprint=0
      L       = 1   ! L initialization
!--- initialization
      CALL INITIAL(kmt,W,T,WC,WT,TXS,qv,qc,qh,qi,VISC,
     &   TE,vth,vti,QVENV,PE,QSAT,RHO,VEL_P,RAD_P,RADIUS,zt,
     &   est,RSURF,alpha,VISCOSITY,Time,mintime,tdur,FMOIST,
     &   heating,cvi,wbar,dqsdz,L,area,dt)
!--- initial print fields:
      izprint  = 0          ! if = 0 => no printout
      if (izprint.ne.0) then
        write(cixx(1:2),'(i2.2)') ixx
        open(2, file = 'debug.'//cixx//'.dat')
        open(19,file='plumegen9.'//cixx//'.gra',
     &      form='unformatted',access='direct',status='unknown',
     &      recl=4*nrectotal)  !PC
!     recl=1*nrectotal) !sx6 e tupay
        call printout (izprint,nrectotal,
     &                        DT,TIME,ZTOP,MINTIME,
     &                        pe,t,te,QV,QSAT,QC,QH,
     &                        qvenv,QI,zt,W,VTH,SC)
        ilastprint=2
      endif
! ******************* model evolution ******************************
      rmaxtime = float(maxtime)
!
      DO WHILE (TIME.LE.RMAXTIME)  !beginning of time loop

!   do itime=1,120

!-- set model top integration
        nm1 = min(kmt, kkmax + deltak)

!-- set timestep
!dt = (zm(2)-zm(1)) / (tstpf * wmax)
        dt = min(5.,(zm(2)-zm(1)) / (tstpf * wmax))
!- - elapsed time, sec
        time = time+dt
!-- elapsed time, minutes
        mintime = 1 + int (time) / 60
        wmax = 1.  !no zeroes allowed.
!************************** BEGIN SPACE LOOP **************************

!-- zerout all model tendencies
        call tend0_plumerise(nm1,wt,tt,qvt,qct,qht,
     &                           qit,vel_t,rad_t)
!-- bounday conditions (k=1)
        L=1
        call lboundf(AREA,RSURF,ALPHA,DT,VISCOSITY,
     &            w,t,qv,qc,qh,qi,rho,wc,
     &            QVENV,PE,te,vel_p,rad_p,VISC,txs,est,
     &            qsat,vth,vti,Time,mintime,tdur,FMOIST,heating,
     &            cvi,wbar,dqsdz,L)
!-- dynamics for the level k>1
!-- W advection
!   call vel_advectc_plumerise(NM1,WC,WT,DNE,DZM)
        call vel_advectc_plumerise(NM1,WC,WT,RHO,DZM)
!-- scalars advection 1
        call scl_advectc_plumerise('SC',NM1,dt,
     &                       rho,vt3dc,vt3df,w,wc,dzm,
     &                       vctr1,vctr2,zt,zm,vt3dk,dzt,
     &                       scr1,vt3dg,t,tt,qv,qc,qi,qh,
     &                       vel_p,rad_p,sc,qvt,qct,qit,qht,
     &                       VEL_T,rad_t,SCT)
!-- scalars advection 2
        !call scl_advectc_plumerise2('SC',NM1)

!-- scalars entrainment, adiabatic
        call scl_misc(NM1,WBAR,ADIABAT,ALPHA,
     &                w,RADIUS,tt,t,te,qvt,qv,qvenv,qct,qc,
     &                qht,qh,qit,qi,VEL_T,VEL_P,VEL_E,rad_t)
!-- scalars dinamic entrainment
        call  scl_dyn_entrain(NM1,nkp,wbar,w,adiabat,alpha,radius,
     &               tt,t,te,qvt,qv,qvenv,qct,qc,qht,qh,qit,qi,
     &               vel_e,vel_p,vel_t,rad_p,rad_t)
!-- gravity wave damping using Rayleigh friction layer fot T
        call damp_grav_wave(1,nm1,deltak,dt,zt,zm,w,t,tt,
     &                       qv,qh,qi,qc,te,pe,qvenv)

!-- microphysics
!   goto 101 ! bypass microphysics
        dt_save=dt
        n_sub_step=3
        dt=dt/float(n_sub_step)

        do i_micro=1,n_sub_step
!-- sedim ?
          call fallpart(NM1,Rho,VTH,w,CVH,CVI,QH,QI,zm,qht,qit,
     &                    vti,VHREL,VIREL)
!-- microphysics
          do L=2,nm1-1
            WBAR    = 0.5*(W(L)+W(L-1))
            ES      = ESAT_PR (T(L))            !BLOB SATURATION VAPOR PRESSURE, EM KPA
            QSAT(L) = (EPS * ES) / (PE(L) - ES)  !BLOB SATURATION LWC G/G DRY AIR
            EST (L) = ES
            RHO (L) = 3483.8 * PE (L) / T (L) ! AIR PARCEL DENSITY , G/M**3
!srf18jun2005
!   IF (W(L) .ge. 0.) DQSDZ = (QSAT(L  ) - QSAT(L-1)) / (ZT(L  ) -ZT(L-1))
!   IF (W(L) .lt. 0.) DQSDZ = (QSAT(L+1) - QSAT(L  )) / (ZT(L+1) -ZT(L  ))
            IF (W(L) .ge. 0.) then
              DQSDZ = (QSAT(L+1) - QSAT(L-1)) / (ZT(L+1 )-ZT(L-1))
            ELSE
              DQSDZ = (QSAT(L+1) - QSAT(L-1)) / (ZT(L+1) -ZT(L-1))
            ENDIF

            call waterbal(qc,qh,qi,Qsat,qv,t,RHO,
     &                    cvi,est,WBAR,dqsdz,dt,L)
          enddo
        enddo
        dt=dt_save
!
101     continue
!
!-- W-viscosity for stability
        call visc_W(nm1,deltak,kmt,
     &                  zt,zm,w,t,qv,qh,qc,qi,vel_p,rad_p,wt,tt,
     &                  visc,qvt,qct,qht,qit,vel_t,rad_t)
!-- update scalars
        call update_plumerise(nm1,'S',w,wt,t,tt,qv,qvt,qc,
     &           qct,qh,qht,qi,qit,VEL_P,VEL_T,rad_p,rad_t,dt)
        call hadvance_plumerise(1,nm1,dt,WC,W,mintime)
!-- Buoyancy
        call buoyancy_plumerise(NM1, T, TE, QV, QVENV,
     &                           QH, QI, QC, WT, SCR1)
!-- Entrainment
        call entrainment(NM1,W,WT,RADIUS,ALPHA,
     &                       WBAR,vel_p,VEL_E)
!-- update W
        call update_plumerise(nm1,'W',w,wt,t,tt,qv,qvt,qc,
     &           qct,qh,qht,qi,qit,VEL_P,VEL_T,rad_p,rad_t,dt)

        call hadvance_plumerise(2,nm1,dt,WC,W,mintime)

!-- misc
        do k=2,nm1
!    pe esta em kpa  - esat do rams esta em mbar = 100 Pa = 0.1 kpa
!    es       = 0.1*esat (t(k)) !blob saturation vapor pressure, em kPa
!    rotina do plumegen calcula em kPa
          es       = esat_pr (t(k))  !blob saturation vapor pressure, em kPa
          qsat(k) = (eps * es) / (pe(k) - es)  !blob saturation lwc g/g dry air
          est (k) = es
          txs (k) = t(k) - te(k)
          rho (k) = 3483.8 * pe (k) / t (k) ! air parcel density , g/m**3
                                       ! no pressure diff with radius
          if ((abs(wc(k))).gt.wmax) then
            wmax = abs(wc(k)) ! keep wmax largest w
          endif
        enddo

! Gravity wave damping using Rayleigh friction layer for W
        call damp_grav_wave(2,nm1,deltak,dt,zt,zm,w,t,
     &                       tt,qv,qh,qi,qc,te,pe,qvenv)
!---
      !- update radius
        do k=2,nm1
          radius(k) = rad_p(k)
        enddo
      !-- try to find the plume top (above surface height)
        kk = 1
        DO WHILE (w (kk) .GT. 1.)
          kk = kk + 1
          ztop =  zm(kk)
        ENDDO
       !
        ztop_(mintime) = ztop
        ztopmax = MAX (ztop, ztopmax)
        kkmax   = MAX (kk  , kkmax  )
       !
       ! if the solution is going to a stationary phase, exit
        IF(mintime > 10) THEN
         !   if(mintime > 20) then
         !    if( abs(ztop_(mintime)-ztop_(mintime-10)) < DZ ) exit
          IF(ABS(ztop_(mintime)-ztop_(mintime-10))< DELZ_THRESOLD) then
            do k=2,nm1
              W_VMD(k,imm) = w(k)
            enddo
            EXIT ! finish the integration
          ENDIF
        ENDIF

        if(ilastprint == mintime) then
          call printout(izprint,nrectotal,
     &                        DT,TIME,ZTOP,MINTIME,
     &                        pe,t,te,QV,QSAT,QC,QH,
     &                        qvenv,QI,zt,W,VTH,SC)
          ilastprint = mintime+1
        endif

      ENDDO   !do next timestep
!
!the last printout
      if (izprint.ne.0) then
        call printout(izprint,nrectotal,
     &                        DT,TIME,ZTOP,MINTIME,
     &                        pe,t,te,QV,QSAT,QC,QH,
     &                        qvenv,QI,zt,W,VTH,SC)
        close (2)
        close (19)
      endif

      RETURN
      END SUBROUTINE MAKEPLUME
!-------------------------------------------------------------------------------
!
      SUBROUTINE BURN(EFLUX,WATER,DT,TIME,TDUR,
     &                FMOIST,HEATING,mintime)
!
!- calculates the energy flux and water content at lboundary
!real, parameter :: HEAT = 21.E6 !Joules/kg
!real, parameter :: HEAT = 15.5E6 !Joules/kg - cerrado
      implicit none
      real, parameter :: HEAT = 19.3E6 !Joules/kg - floresta em Alta Floresta (MT)
      real :: eflux,water
      integer, parameter :: ntime = 200
      real :: DT,TIME,TDUR,FMOIST
      real :: HEATING(ntime)
      integer :: mintime
!
! The emission factor for water is 0.5. The water produced, in kg,
! is then  fuel mass*0.5 + (moist/100)*mass per square meter.
! The fire burns for DT out of TDUR seconds, the total amount of
! fuel burned is AREA*BLOAD*(DT/TDUR) kg. this amount of fuel is
! considered to be spread over area AREA and so the mass burned per
! unit area is BLOAD*(DT/TDUR), and the rate is BLOAD/TDUR.
!
!      print *, 'BURN'
      IF (TIME.GT.TDUR) THEN !is the burn over?
        EFLUX = 0.000001    !prevent a potential divide by zero
        WATER = 0.
        RETURN
      ELSE
!
        EFLUX = HEATING (MINTIME)                          ! Watts/m**2
!  WATER = EFLUX * (DT / HEAT) * (0.5 + FMOIST)       ! kg/m**2
        WATER = EFLUX * (DT / HEAT) * (0.5 + FMOIST) /0.55 ! kg/m**2
        WATER = WATER * 1000.                              ! g/m**2
!
!        print*,'BURN:',time,EFLUX/1.e+9
      ENDIF
!
      RETURN
      END SUBROUTINE BURN
!-------------------------------------------------------------------------------
!
      SUBROUTINE LBOUNDF(AREA,RSURF,ALPHA,DT,VISCOSITY,
     &            w,t,qv,qc,qh,qi,rho,wc,
     &            QVENV,PE,te,vel_p,rad_p,VISC,txs,est,
     &            qsat,vth,vti,Time,mintime,tdur,FMOIST,heating,
     &            cvi,wbar,dqsdz,L)
!
! ********** BOUNDARY CONDITIONS AT ZSURF FOR PLUME AND CLOUD ********
!
! source of equations: J.S. Turner Buoyancy Effects in Fluids
!                      Cambridge U.P. 1973 p.172,
!                      G.A. Briggs Plume Rise, USAtomic Energy Commissio
!                      TID-25075, 1969, P.28
!
! fundamentally a point source below ground. at surface, this produces
! a velocity w(1) and temperature T(1) which vary with time. There is
! also a water load which will first saturate, then remainder go into
! QC(1).
! EFLUX = energy flux at ground,watt/m**2 for the last DT
!
      implicit none
      integer, parameter :: nkp = 200, ntime = 200
      real, parameter :: g = 9.80796, r = 287.04, cp = 1004.6
      real, parameter :: eps = 0.622,tmelt = 273.3
      real, parameter :: pi = 3.14159
      real, parameter :: e1 = 1./3., e2 = 5./3.
      real :: es,  esat, eflux, water,  pres, c1,  c2, f, zv,
     &        denscor, xwater !,ESAT_PR
!  real, external:: esat_pr!
      real :: AREA,RSURF,ALPHA,DT,VISCOSITY
      real,dimension(nkp) ::  w,t,qv,qc,qh,qi,rho,wc,
     &        QVENV,PE,te,vel_p,rad_p,VISC,txs,est,
     &        qsat,vth,vti
      real :: TIME,TDUR,FMOIST
      integer :: MINTIME
      real,dimension(ntime) :: HEATING
      real,dimension(nkp) ::  cvi
      real :: WBAR,dqsdz
      integer :: L
      REAL :: ESAT_PR
!
!      print *,'LBOUNDF'
      QH (1) = QH (2)   !soak up hydrometeors
      QI (1) = QI (2)
      QC (1) = 0.       !no cloud here
!
!
      CALL BURN (EFLUX, WATER,DT,TIME,TDUR,
     &                FMOIST,HEATING,mintime)
!
!  calculate parameters at boundary from a virtual buoyancy point source
!
      PRES = PE (1) * 1000.   !need pressure in N/m**2

      C1 = 5. / (6. * ALPHA)  !alpha is entrainment constant

      C2 = 0.9 * ALPHA

      F = EFLUX / (PRES * CP * PI)

      F = G * R * F * AREA  !buoyancy flux

      ZV = C1 * RSURF  !virtual boundary height

      W (1) = C1 * ( (C2 * F) **E1) / ZV**E1  !boundary velocity

      DENSCOR = C1 * F / G / (C2 * F) **E1 / ZV**E2   !density correction

      T (1) = TE (1) / (1. - DENSCOR)    !temperature of virtual plume at zsurf

!
      WC(1) = W(1)
      VEL_P(1) = 0.
      rad_p(1) = rsurf

!SC(1) = SCE(1)+F/1000.*dt  ! gas/particle (g/g)

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!     match dw/dz,dt/dz at the boundary. F is conserved.
!
!WBAR = W (1) * (1. - 1. / (6. * ZV) )
         !ADVW = WBAR * W (1) / (3. * ZV)
         !ADVT = WBAR * (5. / (3. * ZV) ) * (DENSCOR / (1. - DENSCOR) )
         !ADVC = 0.
         !ADVH = 0.
         !ADVI = 0.
         !ADIABAT = - WBAR * G / CP
      VTH (1) = - 4.
      VTI (1) = - 3.
      TXS (1) = T (1) - TE (1)

      VISC (1) = VISCOSITY

      RHO (1) = 3483.8 * PE (1) / T (1)   !air density at level 1, g/m**3

      XWATER = WATER / (W (1) * DT * RHO (1) )   !firewater mixing ratio

      QV (1) = XWATER + QVENV (1)  !plus what's already there


!  PE esta em kPa  - ESAT do RAMS esta em mbar = 100 Pa = 0.1 kPa
!  ES       = 0.1*ESAT (T(1)) !blob saturation vapor pressure, em kPa
!  rotina do plumegen ja calcula em kPa
      ES       = ESAT_PR (T(1))  !blob saturation vapor pressure, em kPa

      EST  (1)  = ES
      QSAT (1) = (EPS * ES) / (PE (1) - ES)   !blob saturation lwc g/g dry air

      IF (QV (1) .gt. QSAT (1) ) THEN
        QC (1) = QV   (1) - QSAT (1) + QC (1)  !remainder goes into cloud drops
        QV (1) = QSAT (1)
      ENDIF
!
      CALL WATERBAL(qc,qh,qi,Qsat,qv,t,RHO,
     &                    cvi,est,WBAR,dqsdz,dt,L)
!
      RETURN
      END SUBROUTINE LBOUNDF
!-------------------------------------------------------------------------------
!
      SUBROUTINE INITIAL ( kmt,W,T,WC,WT,TXS,qv,qc,qh,qi,VISC,
     &   TE,vth,vti,QVENV,PE,QSAT,RHO,VEL_P,RAD_P,RADIUS,zt,
     &   est,RSURF,alpha,VISCOSITY,Time,mintime,tdur,FMOIST,
     &   heating,cvi,wbar,dqsdz,L,AREA,DT)
!
! ************* SETS UP INITIAL CONDITIONS FOR THE PROBLEM ************
      implicit none
      integer ::  kmt
      real ::     es !,ESAT_PR
            integer, parameter :: nkp = 200, ntime = 200
      integer :: N,k
      real,dimension(nkp) ::  W,T,WC,WT,TXS,qv,qc,qh,qi,VISC,
     &   TE,vth,vti,QVENV,PE,QSAT,RHO,VEL_P,RAD_P,RADIUS,zt,
     &   est
      real :: RSURF,alpha,VISCOSITY
      REAL :: ESAT_PR
      real :: area,dt
      REAL :: TIME,TDUR,FMOIST,wbar,dqsdz
      integer :: mintime,L
      real,dimension(ntime) :: HEATING
      real,dimension(nkp) ::  cvi
!
!      print *,'INITIAL'
      N=kmt
! initialize temperature structure,to the end of equal spaced sounding,
      do k = 1, N
        TXS (k) = 0.0
          W (k) = 0.0
          T (k) = TE(k)       !blob set to environment
          WC(k) = 0.0
          WT(k) = 0.0
          QV(k) = QVENV (k)   !blob set to environment
         VTH(k) = 0.        !initial rain velocity = 0
         VTI(k) = 0.        !initial ice  velocity = 0
          QH(k) = 0.        !no rain
          QI(k) = 0.        !no ice
          QC(k) = 0.        !no cloud drops
!  PE esta em kPa  - ESAT do RAMS esta em mbar = 100 Pa = 0.1 kPa
!  ES       = 0.1*ESAT (T(k)) !blob saturation vapor pressure, em kPa
!  rotina do plumegen calcula em kPa
         ES       = ESAT_PR (T(k))  !blob saturation vapor pressure, em kPa
         EST  (k) = ES
         QSAT (k) = (.622 * ES) / (PE (k) - ES) !saturation lwc g/g
         RHO  (k) = 3483.8 * PE (k) / T (k)     !dry air density g/m**3
         VEL_P(k) = 0.
         rad_p(k) = 0.
      enddo

! Initialize the entrainment radius, Turner-style plume
      radius(1) = rsurf
      do k=2,N
        radius(k) = radius(k-1)+(6./5.)*alpha*(zt(k)-zt(k-1))
      enddo
! Initialize the entrainment radius, Turner-style plume
      radius(1) = rsurf
      rad_p(1)  = rsurf
      DO k=2,N
        radius(k) = radius(k-1)+(6./5.)*alpha*(zt(k)-zt(k-1))
        rad_p(k)  = radius(k)
      ENDDO

!  Initialize the viscosity
      VISC (1) = VISCOSITY
      do k=2,N
      !VISC (k) = VISCOSITY!max(1.e-3,visc(k-1) - 1.* VISCOSITY/float(nkp))
        VISC (k) = max(1.e-3,visc(k-1) - 1.* VISCOSITY/float(nkp))
      enddo
!--   Initialize gas/concentration
!DO k =10,20
!   SC(k) = 20.
!ENDDO
!stop 333

      CALL LBOUNDF(AREA,RSURF,ALPHA,DT,VISCOSITY,
     &            w,t,qv,qc,qh,qi,rho,wc,
     &            QVENV,PE,te,vel_p,rad_p,VISC,txs,est,
     &            qsat,vth,vti,Time,mintime,tdur,FMOIST,heating,
     &            cvi,wbar,dqsdz,L)

      RETURN
      END SUBROUTINE INITIAL
!-------------------------------------------------------------------------------
!
      subroutine damp_grav_wave(ifrom,nm1,deltak,dt,zt,zm,w,t,tt,qv,
     &                          qh,qi,qc,te,pe,qvenv)
      implicit none
      integer nm1,ifrom,deltak
      real dt
      real, dimension(nm1) :: w,t,tt,qv,qh,qi,
     &                  qc,te,pe,qvenv,dummy,zt,zm
!      print *, 'damp_grav_wave'
      if(ifrom==1) then
        call friction(ifrom,nm1,deltak,dt,zt,zm,t,tt,te)
!call friction(ifrom,nm1,dt,zt,zm,qv,qvt,qvenv)
        return
      endif

      dummy(:) = 0.
      if(ifrom==2) then
        call friction(ifrom,nm1,deltak,dt,zt,zm,w,dummy,
     &                           dummy)
      endif
!call friction(ifrom,nm1,dt,zt,zm,qi,qit ,dummy)
!call friction(ifrom,nm1,dt,zt,zm,qh,qht ,dummy)
!call friction(ifrom,nm1,dt,zt,zm,qc,qct ,dummy)
      return
      end subroutine damp_grav_wave
!-------------------------------------------------------------------------------
!
      subroutine friction(ifrom,nm1,deltak,dt,zt,zm,var1,vart,var2)
      implicit none
      real, dimension(nm1) :: var1,var2,vart,zt,zm
      integer k,nfpt,kf,nm1,ifrom,deltak
      real zmkf,ztop,distim,c1,c2,dt
!      print *, 'friction'
!nfpt=50
!kf = nm1 - nfpt
!kf = nm1 - int(deltak/2)
      kf = nm1 - int(deltak)

      zmkf = zm(kf) !old: float(kf )*dz
      ztop = zm(nm1)
!distim = min(4.*dt,200.)
!distim = 60.
      distim = min(3.*dt,60.)

      c1 = 1. / (distim * (ztop - zmkf))
      c2 = dt * c1

      if(ifrom == 1) then
        do k = nm1,2,-1
         if (zt(k) .le. zmkf) cycle
         vart(k) = vart(k)   + c1 * (zt(k) - zmkf)*(var2(k) - var1(k))
        enddo
      elseif(ifrom == 2) then
        do k = nm1,2,-1
         if (zt(k) .le. zmkf) cycle
         var1(k) =  var1(k) + c2 * (zt(k) - zmkf)*(var2(k) - var1(k))
        enddo
      endif
      return
      end subroutine friction
!-------------------------------------------------------------------------------
!
      subroutine vel_advectc_plumerise(m1,wc,wt,rho,dzm)

      implicit none
      integer, parameter :: nkp = 200
      integer :: k,m1
!      real, dimension(m1) :: wc,wt,flxw,dzm,rho
!      real, dimension(m1) :: dn0 ! var local
      real, dimension(nkp) :: wc,wt,flxw,dzm,rho
      real, dimension(nkp) :: dn0 ! var local
      real :: c1z
!      print *,'vel_advectc_plumerise'
!dzm(:)= 1./dz

      dn0(1:m1)=rho(1:m1)*1.e-3 ! converte de cgs para mks

      flxw(1) = wc(1) * dn0(1)

      do k = 2,m1-1
        flxw(k) = wc(k) * .5 * (dn0(k) + dn0(k+1))
      enddo

! Compute advection contribution to W tendency

      c1z = .5

      do k = 2,m1-2
         wt(k) = wt(k)
     &         + c1z*dzm(k)/(dn0(k)+dn0(k+1))*
     &       ((flxw(k)+flxw(k-1))*(wc(k)+wc(k-1))
     &       -(flxw(k)+flxw(k+1))*(wc(k)+wc(k+1))
     &       +(flxw(k+1)-flxw(k-1))*2.*wc(k))
      enddo

      return
      end subroutine vel_advectc_plumerise
!-------------------------------------------------------------------------------
!
      subroutine hadvance_plumerise(iac,m1,dt,wc,wp,mintime)

      implicit none
      integer :: k,iac
      integer :: m1,mintime
      integer, parameter :: nkp = 200
!      real, dimension(m1) :: dummy, wc,wt,wp
      real, dimension(nkp),intent(inout)  :: wc,wp
!      real, dimension(nkp) :: dummy
      real eps,dt
!     It is here that the Asselin filter is applied.  For the velocities
!     and pressure, this must be done in two stages, the first when
!     IAC=1 and the second when IAC=2.

!      print *, 'hadvance_plumerise'
      eps = .2
      if(mintime == 1) eps=0.5

!     For both IAC=1 and IAC=2, call PREDICT for U, V, W, and P.
!
      call predict_plumerise(m1,wc,wp,iac,2.*dt,eps)
!do k=1,m1
!   print*,'W-HAD',k,wc(k),wp(k),wt(k)
!enddo
      return
      end subroutine hadvance_plumerise
!-------------------------------------------------------------------------------
!
      subroutine predict_plumerise(npts,ac,ap,iac,dtlp,epsu)
      implicit none
      integer :: npts,iac,m
      real :: epsu,dtlp
      integer, parameter :: nkp = 200
      real, dimension(nkp),intent(inout) :: ac,ap
      real, dimension(nkp) ::  af

!     For IAC=3, this routine moves the arrays AC and AP forward by
!     1 time level by adding in the prescribed tendency. It also
!     applies the Asselin filter given by:

!              {AC} = AC + EPS * (AP - 2 * AC + AF)

!     where AP,AC,AF are the past, current and future time levels of A.
!     All IAC=1 does is to perform the {AC} calculation without the AF
!     term present.  IAC=2 completes the calculation of {AC} by adding
!     the AF term only, and advances AC by filling it with input AP
!     values which were already updated in ACOUSTC.
!
!      print *, 'predict_plumerise'
      if (iac .eq. 1) then
         do m = 1,npts
            ac(m) = ac(m) + epsu * (ap(m) - 2. * ac(m))
         enddo
         return
      elseif (iac .eq. 2) then
         do m = 1,npts
            af(m) = ap(m)
            ap(m) = ac(m) + epsu * af(m)
         enddo
!elseif (iac .eq. 3) then
!   do m = 1,npts
!      af(m) = ap(m) + dtlp * fa(m)
!   enddo
!   if (ngrid .eq. 1 .and. ipara .eq. 0) call cyclic(nzp,nxp,nyp,af,'T')
!   do m = 1,npts
!      ap(m) = ac(m) + epsu * (ap(m) - 2. * ac(m) + af(m))
!   enddo
      endif

      do m = 1,npts
        ac(m) = af(m)
      enddo
      return
      end subroutine predict_plumerise
!-------------------------------------------------------------------------------
!
      subroutine  buoyancy_plumerise(m1, T, TE, QV, QVENV, QH, QI, QC,
     &                                WT, scr1)
      implicit none
      integer :: k,m1
      real, parameter :: g = 9.8, eps = 0.622, gama = 0.5 ! mass virtual coeff.
      integer, parameter :: nkp = 200
      real, dimension(nkp),intent(inout) :: T, TE, QV, QVENV,
     &            QH, QI, QC, WT, scr1
!      real, dimension(m1) :: T, TE, QV, QVENV, QH, QI, QC, WT, scr1
      real :: TV,TVE,QWTOTL,umgamai
      real, parameter :: mu = 0.15
!      print *, 'buoyancy_plumerise'
!- orig
      umgamai = 1./(1.+gama) ! compensa a falta do termo de aceleracao associado `as
                             ! das pertubacoes nao-hidrostaticas no campo de pressao

!- new                 ! Siesbema et al, 2004
!umgamai = 1./(1.-2.*mu)

      do k = 2,m1-1

          TV =   T(k) * (1. + (QV(k)   /EPS))/(1. + QV(k)   )  !blob virtual temp.
          TVE = TE(k) * (1. + (QVENV(k)/EPS))/(1. + QVENV(k))  !and environment

          QWTOTL = QH(k) + QI(k) + QC(k)                       ! QWTOTL*G is drag
!- orig
!scr1(k)= G*( umgamai*(  TV - TVE) / TVE   - QWTOTL)
          scr1(k)= G*  umgamai*( (TV - TVE) / TVE   - QWTOTL)
!if(k .lt. 10)print*,'BT',k,TV,TVE,TVE,QWTOTL
      enddo

      do k = 2,m1-2
          wt(k) = wt(k)+0.5*(scr1(k)+scr1(k+1))
!   print*,'W-BUO',k,wt(k),scr1(k),scr1(k+1)
      enddo

      end subroutine  buoyancy_plumerise
!-------------------------------------------------------------------------------
!
      subroutine ENTRAINMENT(m1,w,wt,radius,ALPHA,
     &                       WBAR,vel_p,VEL_E)
      implicit none
      integer :: k,m1
      integer, parameter :: nkp = 200
!      real, dimension(m1) :: w,wt
!      real, dimension(m1) :: radius
      real, dimension(nkp) :: w,wt
      real, dimension(nkp) :: radius
      real :: DMDTM,WBAR,RADIUSBAR,umgamai,DYN_ENTR
      real :: ALPHA
      real, parameter :: mu = 0.15 ,gama = 0.5 ! mass virtual coeff.

      real,dimension(nkp) ::  vel_p,VEL_E
!      print *,'ENTRAINMENT'
!- new - Siesbema et al, 2004
!umgamai = 1./(1.-2.*mu)

!- orig
!umgamai = 1
      umgamai = 1./(1.+gama) ! compensa a falta do termo de aceleracao associado as
                       ! das pertubacoes nao-hidrostaticas no campo de pressao

!
!-- ALPHA/RADIUS(L) = (1/M)DM/DZ  (W 14a)
      do k=2,m1-1
        WBAR=W(k)
        RADIUSBAR = 0.5*(RADIUS(k) + RADIUS(k-1))
! orig
!DMDTM =           2. * ALPHA * ABS (WBAR) / RADIUS_BAR  != (1/M)DM/DT
        DMDTM = umgamai * 2. * ALPHA * ABS (WBAR) / RADIUSBAR  != (1/M)DM/DT

!--  DMDTM*W(L) entrainment,
        wt(k) = wt(k)  - DMDTM*ABS (WBAR)

      !if(VEL_P (k) - VEL_E (k) > 0.) cycle

       !-   dynamic entrainment
        DYN_ENTR =  (2./3.1416)*0.5*ABS (VEL_P(k)-VEL_E(k)+VEL_P(k-1)-
     &               VEL_E(k-1)) /RADIUSBAR

        wt(k) = wt(k)  - DYN_ENTR*ABS (WBAR)
       !- entraiment acceleration for output only
       !dwdt_entr(k) =  - DMDTM*ABS (WBAR)- DYN_ENTR*ABS (WBAR)
      enddo
      end subroutine  ENTRAINMENT
!-------------------------------------------------------------------------------
!
      subroutine scl_advectc_plumerise(varn,mzp,dt,
     &                       rho,vt3dc,vt3df,w,wc,dzm,
     &                       vctr1,vctr2,zt,zm,vt3dk,dzt,
     &                       scr1,vt3dg,t,tt,qv,qc,qi,qh,
     &                       vel_p,rad_p,sc,qvt,qct,qit,qht,
     &                       VEL_T,rad_t,SCT)
      implicit none
      integer, parameter :: nkp = 200
      integer :: mzp
      character(len=*) :: varn
      real :: dtlto2,dt
      integer :: k
      real,dimension(nkp),intent(inout) :: rho,vt3dc,vt3df,w,wc,dzm,
     &                       vctr1,vctr2,zt,zm,vt3dk,dzt,
     &                       scr1,vt3dg,t,tt,qv,qc,qi,qh,
     &                       vel_p,rad_p,sc,qvt,qct,qit,qht,
     &                       VEL_T,rad_t,SCT
!      print *,' scl_advectc_plumerise'
!  wp => w
!- Advect  scalars
      dtlto2   = .5 * dt
!  vt3dc(1) =      (w(1) + wc(1)) * dtlto2 * dne(1)
      vt3dc(1) =      (w(1) + wc(1)) * dtlto2 * rho(1)*1.e-3!converte de CGS p/ MKS
      vt3df(1) = .5 * (w(1) + wc(1)) * dtlto2 * dzm(1)

      do k = 2,mzp
!     vt3dc(k) =  (w(k) + wc(k)) * dtlto2 *.5 * (dne(k) + dne(k+1))
        vt3dc(k) =  (w(k) + wc(k)) * dtlto2 *0.5 * (rho(k) +
     &              rho(k+1))*1.e-3
        vt3df(k) =  (w(k) + wc(k)) * dtlto2 *0.5 *  dzm(k)
      enddo


!-srf-24082005
!  do k = 1,mzp-1
      do k = 1,mzp
        vctr1(k) = (zt(k+1) - zm(k)) * dzm(k)
        vctr2(k) = (zm(k)   - zt(k)) * dzm(k)
!      vt3dk(k) = dzt(k) / dne(k)
        vt3dk(k) = dzt(k) /(rho(k)*1.e-3)
      enddo

!      scalarp => scalar_tab(n,ngrid)%var_p
!      scalart => scalar_tab(n,ngrid)%var_t

!- temp advection tendency (TT)
      scr1=T
!      call fa_zc_plumerise(mzp
!     &                    ,T        ,scr1  (1)
!     &                    ,vt3dc (1) ,vt3df (1)
!     &                    ,vt3dg (1) ,vt3dk (1)
!     &                    ,vctr1,vctr2  )
      call fa_zc_plumerise(mzp
     &                    ,T     ,scr1
     &                    ,vt3dc ,vt3df
     &                    ,vt3dg ,vt3dk
     &                    ,vctr1,vctr2  )
      call advtndc_plumerise(mzp,T,scr1,TT,dt)
!- water vapor advection tendency (QVT)
      scr1=QV
      call fa_zc_plumerise(mzp
     &                     ,QV       ,scr1
     &                     ,vt3dc  ,vt3df
     &                     ,vt3dg  ,vt3dk
     &                     ,vctr1,vctr2    )
      call advtndc_plumerise(mzp,QV,scr1,QVT,dt)

!- liquid advection tendency (QCT)
      scr1=QC
      call fa_zc_plumerise(mzp
     &                     ,QC       ,scr1
     &                     ,vt3dc ,vt3df
     &                     ,vt3dg ,vt3dk
     &                     ,vctr1,vctr2    )

      call advtndc_plumerise(mzp,QC,scr1,QCT,dt)

!- ice advection tendency (QIT)
      scr1=QI
      call fa_zc_plumerise(mzp
     &                  ,QI        ,scr1
     &                  ,vt3dc  ,vt3df
     &                  ,vt3dg  ,vt3dk
     &                  ,vctr1,vctr2    )

      call advtndc_plumerise(mzp,QI,scr1,QIT,dt)

!- hail/rain advection tendency (QHT)
!   if(ak1 > 0. .or. ak2 > 0.) then

      scr1=QH
      call fa_zc_plumerise(mzp
     &                    ,QH       ,scr1
     &                    ,vt3dc  ,vt3df
     &                    ,vt3dg  ,vt3dk
     &                    ,vctr1,vctr2         )

      call advtndc_plumerise(mzp,QH,scr1,QHT,dt)
!   endif
!- horizontal wind advection tendency (VEL_T)
      scr1=VEL_P
      call fa_zc_plumerise(mzp
     &              ,VEL_P     ,scr1
     &          ,vt3dc  ,vt3df
     &          ,vt3dg  ,vt3dk
     &          ,vctr1,vctr2         )

      call advtndc_plumerise(mzp,VEL_P,scr1,VEL_T,dt)

!- vertical radius transport

      scr1=rad_p
      call fa_zc_plumerise(mzp
     &                  ,rad_p     ,scr1
     &                  ,vt3dc  ,vt3df
     &                  ,vt3dg  ,vt3dk
     &                  ,vctr1,vctr2         )

      call advtndc_plumerise(mzp,rad_p,scr1,rad_t,dt)


      return
!
!- gas/particle advection tendency (SCT)
!    if(varn == 'SC')return
      scr1=SC
      call fa_zc_plumerise(mzp
     &                 ,SC   ,scr1
     &                 ,vt3dc  ,vt3df
     &                 ,vt3dg  ,vt3dk
     &                 ,vctr1,vctr2      )

      call advtndc_plumerise(mzp,SC,scr1,SCT,dt)

      return
      end subroutine scl_advectc_plumerise
!-------------------------------------------------------------------------------
!
      subroutine fa_zc_plumerise(m1,scp,scr1,vt3dc,vt3df,vt3dg,
     &             vt3dk,vctr1,vctr2)

      implicit none
      integer :: m1,k
      real :: dfact
      integer, parameter :: nkp = 200
      real, dimension(nkp) :: scp,scr1,vt3dc,vt3df,vt3dg,vt3dk
      real, dimension(nkp) :: vctr1,vctr2
!      real, dimension(m1) :: scp,scr1,vt3dc,vt3df,vt3dg,vt3dk
!      real, dimension(m1) :: vctr1,vctr2
!      print *, 'fa_zc_plumerise'
      dfact = .5

! Compute scalar flux VT3DG
      do k = 1,m1-1
         vt3dg(k) = vt3dc(k)
     &             * (vctr1(k) * scr1(k)
     &             +  vctr2(k) * scr1(k+1)
     &             +  vt3df(k) * (scr1(k) - scr1(k+1)))
      enddo

! Modify fluxes to retain positive-definiteness on scalar quantities.
!    If a flux will remove 1/2 quantity during a timestep,
!    reduce to first order flux. This will remain positive-definite
!    under the assumption that ABS(CFL(i)) + ABS(CFL(i-1)) < 1.0 if
!    both fluxes are evacuating the box.

      do k = 1,m1-1
        if (vt3dc(k) .gt. 0.) then
          if (vt3dg(k) * vt3dk(k)    .gt. dfact * scr1(k)) then
            vt3dg(k) = vt3dc(k) * scr1(k)
          endif
        elseif (vt3dc(k) .lt. 0.) then
          if (-vt3dg(k) * vt3dk(k+1) .gt. dfact * scr1(k+1)) then
            vt3dg(k) = vt3dc(k) * scr1(k+1)
          endif
        endif
      enddo

! Compute flux divergence

      do k = 2,m1-1
          scr1(k) = scr1(k)
     &             + vt3dk(k) * ( vt3dg(k-1) - vt3dg(k)
     &             + scp  (k) * ( vt3dc(k)   - vt3dc(k-1)))
      enddo
      return
      end subroutine fa_zc_plumerise
!-------------------------------------------------------------------------------
!
      subroutine advtndc_plumerise(m1,scp,sca,sct,dtl)
      implicit none
      integer :: m1,k
      real :: dtl,dtli
      integer, parameter :: nkp = 200
!      real, dimension(m1) :: scp,sca,sct
      real, dimension(nkp),intent(in) :: scp,sca
      real, dimension(nkp),intent(inout) :: sct
!      print *,'advtndc_plumerise'
      dtli = 1. / dtl
      do k = 2,m1-1
        sct(k) = sct(k) + (sca(k)-scp(k)) * dtli
      enddo
      return
      end subroutine advtndc_plumerise
!-------------------------------------------------------------------------------
!
      subroutine tend0_plumerise(nm1,wt,tt,qvt,qct,qht,
     &                           qit,vel_t,rad_t)
      implicit none
      integer, parameter :: nkp = 200
      integer :: nm1
      real,dimension(nkp),intent(inout) ::
     &        wt,tt,qvt,qct,qht,qit,vel_t,rad_t
!      print *,'tend0_plumerise'
       wt(1:nm1)  = 0.
       tt(1:nm1)  = 0.
      qvt(1:nm1)  = 0.
      qct(1:nm1)  = 0.
      qht(1:nm1)  = 0.
      qit(1:nm1)  = 0.
      vel_t(1:nm1)  = 0.
      rad_t(1:nm1)  = 0.
!sct(1:nm1)  = 0.
      end subroutine tend0_plumerise

!     ****************************************************************

      subroutine scl_misc(m1,WBAR,ADIABAT,ALPHA,
     &                    w,RADIUS,tt,t,te,qvt,qv,qvenv,qct,qc,
     &                    qht,qh,qit,qi,VEL_T,VEL_P,VEL_E,rad_t)
      implicit none
      real, parameter :: g = 9.81, cp=1004.
      integer m1,k
      real dmdtm
      integer, parameter :: nkp = 200
      real :: WBAR,ADIABAT,ALPHA
      real,dimension(nkp) ::  w,RADIUS,tt,t,te,qvt,qv,qvenv,qct,qc,
     &                      qht,qh,qit,qi,VEL_T,VEL_P,VEL_E,rad_t
!      print *, 'scl_misc'
      do k=2,m1-1
        WBAR    = 0.5*(W(k)+W(k-1))
!-- dry adiabat
        ADIABAT = - WBAR * G / CP
!
!-- entrainment
        DMDTM = 2. * ALPHA * ABS (WBAR) / RADIUS (k)  != (1/M)DM/DT

!-- tendency temperature = adv + adiab + entrainment
        TT(k) = TT(K) + ADIABAT - DMDTM * ( T  (k) -    TE (k) )

!-- tendency water vapor = adv  + entrainment
        QVT(K) = QVT(K)         - DMDTM * ( QV (k) - QVENV (k) )

        QCT(K) = QCT(K)       - DMDTM * ( QC (k)  )
        QHT(K) = QHT(K)       - DMDTM * ( QH (k)  )
        QIT(K) = QIT(K)       - DMDTM * ( QI (k)  )

      !-- tendency horizontal speed = adv  + entrainment
        VEL_T(K) = VEL_T(K)     - DMDTM * ( VEL_P (k) - VEL_E (k) )

      !-- tendency horizontal speed = adv  + entrainment
        rad_t(K) = rad_t(K)     + 0.5*DMDTM*(6./5.)*RADIUS (k)
!-- tendency gas/particle = adv  + entrainment
!      SCT(K) = SCT(K)         - DMDTM * ( SC (k) -   SCE (k) )

      enddo
      end subroutine scl_misc
!     ****************************************************************

      SUBROUTINE scl_dyn_entrain(m1,nkp,wbar,w,adiabat,alpha,radius,
     &                     tt,t,te,qvt,qv,qvenv,qct,qc,qht,qh,qit,qi,
     &                     vel_e,vel_p,vel_t,rad_p,rad_t)
      implicit none
      INTEGER , INTENT(IN)    :: m1

! plumegen_coms
      INTEGER , INTENT(IN)    :: nkp
      REAL    , INTENT(INOUT) :: wbar
      REAL    , INTENT(IN)    :: w(nkp)
      REAL    , INTENT(INOUT) :: adiabat
      REAL    , INTENT(IN)    :: alpha
      REAL    , INTENT(IN)    :: radius(nkp)
      REAL    , INTENT(INOUT) :: tt(nkp)
      REAL    , INTENT(IN)    :: t(nkp)
      REAL    , INTENT(IN)    :: te(nkp)
      REAL    , INTENT(INOUT) :: qvt(nkp)
      REAL    , INTENT(IN)    :: qv(nkp)
      REAL    , INTENT(IN)    :: qvenv(nkp)
      REAL    , INTENT(INOUT) :: qct(nkp)
      REAL    , INTENT(IN)    :: qc(nkp)
      REAL    , INTENT(INOUT) :: qht(nkp)
      REAL    , INTENT(IN)    :: qh(nkp)
      REAL    , INTENT(INOUT) :: qit(nkp)
      REAL    , INTENT(IN)    :: qi(nkp)

      REAL    , INTENT(IN)    :: vel_e(nkp)
      REAL    , INTENT(IN)    :: vel_p(nkp)
      REAL    , INTENT(INOUT) :: vel_t(nkp)
      REAL    , INTENT(INOUT) :: rad_T(nkp)
      REAL    , INTENT(IN)    :: rad_p(nkp)

      real, parameter :: g = 9.81, cp=1004., pi=3.1416

      integer k
      real dmdtm
!      print *,'scl_dyn_entrain'
      DO k=2,m1-1
      !
      !-- tendency horizontal radius from dyn entrainment
       !rad_t(K) = rad_t(K)   + (vel_e(k)-vel_p(k)) /pi
         rad_t(K) = rad_t(K)   + ABS((vel_e(k)-vel_p(k)))/pi

      !-- entrainment
       !DMDTM = (2./3.1416)  *     (VEL_E (k) - VEL_P (k)) / RADIUS (k)
         DMDTM = (2./3.1416)  *  ABS(VEL_E (k) - VEL_P (k)) /RADIUS (k)

      !-- tendency horizontal speed  from dyn entrainment
         VEL_T(K) = VEL_T(K)     - DMDTM * ( VEL_P (k) - VEL_E (k) )

      !     if(VEL_P (k) - VEL_E (k) > 0.) cycle

      !-- tendency temperature  from dyn entrainment
         TT(k) = TT(K)    - DMDTM * ( T (k) - TE  (k) )

      !-- tendency water vapor  from dyn entrainment
         QVT(K) = QVT(K)        - DMDTM * ( QV (k) - QVENV (k) )

         QCT(K) = QCT(K)        - DMDTM * ( QC (k)  )
         QHT(K) = QHT(K)        - DMDTM * ( QH (k)  )
         QIT(K) = QIT(K)        - DMDTM * ( QI (k)  )

      !-- tendency gas/particle  from dyn entrainment
      !  SCT(K) = SCT(K)     - DMDTM * ( SC (k) - SCE (k) )

      ENDDO
      END SUBROUTINE scl_dyn_entrain

!     ****************************************************************

      subroutine visc_W(m1,deltak,kmt,
     &                  zt,zm,w,t,qv,qh,qc,qi,vel_p,rad_p,wt,tt,
     &                  visc,qvt,qct,qht,qit,vel_t,rad_t)
      implicit none
      integer, parameter :: nkp = 200
      integer m1,k,deltak,kmt,m2
      real dz1t,dz1m,dz2t,dz2m,d2wdz,d2tdz  ,d2qvdz ,d2qhdz ,d2qcdz ,
     &      d2qidz ,d2scdz, d2vel_pdz,d2rad_dz
      real,dimension(nkp),intent(in) ::  zt,zm,w,t,qv,qh,qc,qi,
     &                        vel_p,rad_p,visc
      real,dimension(nkp),intent(inout) :: wt,tt,
     &                        qvt,qct,qht,qit,vel_t,rad_t
!      print *,'visc_W'
!srf--- 17/08/2005
!m2=min(m1+deltak,kmt)
      m2=min(m1,kmt)
!do k=2,m1-1
      do k=2,m2-1
        DZ1T   = 0.5*(ZT(K+1)-ZT(K-1))
        DZ2T   = VISC (k) / (DZ1T * DZ1T)
        DZ1M   = 0.5*(ZM(K+1)-ZM(K-1))
        DZ2M   = VISC (k) / (DZ1M * DZ1M)
        D2WDZ  = (W  (k + 1) - 2 * W  (k) + W  (k - 1) ) * DZ2M
        D2TDZ  = (T  (k + 1) - 2 * T  (k) + T  (k - 1) ) * DZ2T
        D2QVDZ = (QV (k + 1) - 2 * QV (k) + QV (k - 1) ) * DZ2T
        D2QHDZ = (QH (k + 1) - 2 * QH (k) + QH (k - 1) ) * DZ2T
        D2QCDZ = (QC (k + 1) - 2 * QC (k) + QC (k - 1) ) * DZ2T
        D2QIDZ = (QI (k + 1) - 2 * QI (k) + QI (k - 1) ) * DZ2T
        !D2SCDZ = (SC (k + 1) - 2 * SC (k) + SC (k - 1) ) * DZ2T
        d2vel_pdz=(vel_P(k+1)- 2 * vel_P  (k) + vel_P  (k - 1) ) * DZ2T
        d2rad_dz =(rad_p(k+1)- 2 * rad_p  (k) + rad_p  (k - 1) ) * DZ2T

        WT(k) =   WT(k) + D2WDZ
        TT(k) =   TT(k) + D2TDZ
        QVT(k) =  QVT(k) + D2QVDZ
        QCT(k) =  QCT(k) + D2QCDZ
        QHT(k) =  QHT(k) + D2QHDZ
        QIT(k) =  QIT(k) + D2QIDZ
        vel_t(k) =   vel_t(k) + d2vel_pdz
        rad_t(k) =   rad_t(k) + d2rad_dz
        !SCT(k) =  SCT(k) + D2SCDZ
      enddo
      return
      end subroutine visc_W

!     ****************************************************************

      subroutine update_plumerise(m1,varn,w,wt,t,tt,qv,qvt,qc,
     &          qct,qh,qht,qi,qit,VEL_P,VEL_T,rad_p,rad_t,dt)
      implicit none
      integer m1,k
      integer, parameter :: nkp = 200
      character(len=*) :: varn
      real,dimension(nkp),intent(in) ::  wt,tt,qvt,qct,qht,
     &                        qit,VEL_T,rad_t
      real,dimension(nkp),intent(inout) :: w,t,qv,qc,qh,
     &                        qi,VEL_P,rad_p
      real,intent(inout) :: DT
!      print *,'update_plumerise'
      if(varn == 'W') then

        do k=2,m1-1
          W(k) =  W(k) +  WT(k) * DT
        enddo
        return

      else
        do k=2,m1-1
           T(k) =  T(k) +  TT(k) * DT
          QV(k) = QV(k) + QVT(k) * DT
          QC(k) = QC(k) + QCT(k) * DT !cloud drops travel with air
          QH(k) = QH(k) + QHT(k) * DT
          QI(k) = QI(k) + QIT(k) * DT
! SC(k) = SC(k) + SCT(k) * DT

!srf---18jun2005
          QV(k) = max(0., QV(k))
          QC(k) = max(0., QC(k))
          QH(k) = max(0., QH(k))
          QI(k) = max(0., QI(k))

          VEL_P(k) =  VEL_P(k) + VEL_T(k) * DT
          rad_p(k) =  rad_p(k) + rad_t(k) * DT
! SC(k) = max(0., SC(k))
        enddo
        return
      endif
      end subroutine update_plumerise
!-------------------------------------------------------------------------------
!
      subroutine fallpart(m1,Rho,VTH,w,CVH,CVI,QH,QI,zm,qht,qit,
     &                    vti,VHREL,VIREL)
      implicit none
      integer m1,k
      real vtc, dfhz,dfiz,dz1
!srf==================================
!   verificar se o gradiente esta correto
!
!srf==================================
!
!     XNO=1.E7  [m**-4] median volume diameter raindrop,Kessler
!     VC = 38.3/(XNO**.125), median volume fallspeed eqn., Kessler
!     for ice, see (OT18), use F0=0.75 per argument there. rho*q
!     values are in g/m**3, velocities in m/s

      real, PARAMETER :: VCONST = 5.107387, EPS = 0.622, F0 = 0.75
      real, PARAMETER :: G = 9.81, CP = 1004.
!
      integer, parameter :: nkp = 200
      real,dimension(nkp) ::  Rho,VTH,w,CVH,CVI,QH,QI,zm,
     &                        qht,qit,vti
      real :: VHREL,VIREL
!      print *,'fallpart'
      do k=2,m1-1

         VTC = VCONST * RHO (k) **.125   ! median volume fallspeed (KTable4)

!  hydrometeor assembly velocity calculations (K Table4)
!  VTH(k)=-VTC*QH(k)**.125  !median volume fallspeed, water
         VTH (k) = - 4.     !small variation with qh

         VHREL = W (k) + VTH (k)  !relative to surrounding cloud

!  rain ventilation coefficient for evaporation
         CVH(k) = 1.6 + 0.57E-3 * (ABS (VHREL) ) **1.5
!
!  VTI(k)=-VTC*F0*QI(k)**.125    !median volume fallspeed,ice
         VTI (k) = - 3.                !small variation with qi

         VIREL = W (k) + VTI (k)       !relative to surrounding cloud
!
!  ice ventilation coefficient for sublimation
         CVI(k) = 1.6 + 0.57E-3 * (ABS (VIREL) ) **1.5 / F0
!
!
         IF (VHREL.GE.0.0) THEN
           DFHZ=QH(k)*(RHO(k  )*VTH(k  )-RHO(k-1)*VTH(k-1))/RHO(k-1)
         ELSE
           DFHZ=QH(k)*(RHO(k+1)*VTH(k+1)-RHO(k  )*VTH(k  ))/RHO(k)
         ENDIF
         IF (VIREL.GE.0.0) THEN
           DFIZ=QI(k)*(RHO(k  )*VTI(k  )-RHO(k-1)*VTI(k-1))/RHO(k-1)
         ELSE
           DFIZ=QI(k)*(RHO(k+1)*VTI(k+1)-RHO(k  )*VTI(k  ))/RHO(k)
         ENDIF

         DZ1=ZM(K)-ZM(K-1)

         qht(k) = qht(k) - DFHZ / DZ1 !hydrometeors don't

         qit(k) = qit(k) - DFIZ / DZ1  !nor does ice? hail, what about

      enddo
      end subroutine fallpart
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!
      subroutine printout (izprint,nrectotal,
     &                        DT,TIME,ZTOP,MINTIME,
     &                        pe,t,te,QV,QSAT,QC,QH,
     &                        qvenv,QI,zt,W,VTH,SC)
      real, parameter :: tmelt = 273.3
      integer, save :: nrec
      data nrec/0/
      integer :: ko,izprint,interval,nrectotal
      real :: pea, btmp,etmp,vap1,vap2,gpkc,gpkh,gpki,deficit
      integer, parameter :: nkp = 200
      real :: DT,TIME,ZTOP
      INTEGER :: MINTIME
      real,dimension(nkp) ::  pe,t,te,QV,QSAT,QC,QH,
     &                        qvenv,QI,zt,W,VTH,SC
!      print *,'printout'
      interval = 1              !debug time interval,min
!
      IF (IZPRINT.EQ.0) RETURN
      IF(MINTIME == 1) then
        nrec = 0
      endif
!
      WRITE (2, 430) MINTIME, DT, TIME
      WRITE (2, 431) ZTOP
      WRITE (2, 380)
!
! do the print
!
      DO 390 KO = 1, nrectotal, interval
         PEA = PE (KO) * 10.       !pressure is stored in decibars(kPa),print in mb;
         BTMP = T (KO) - TMELT     !temps in Celsius
         ETMP = T (KO) - TE (KO)   !temperature excess
         VAP1 = QV (KO)   * 1000.  !printout in g/kg for all water,
         VAP2 = QSAT (KO) * 1000.  !vapor (internal storage is in g/g)
         GPKC = QC (KO)   * 1000.  !cloud water
         GPKH = QH (KO)   * 1000.  !raindrops
         GPKI = QI (KO)   * 1000.  !ice particles
         DEFICIT = VAP2 - VAP1     !vapor deficit
!
         WRITE (2, 400) zt(KO)/1000., PEA, W (KO), BTMP, ETMP, VAP1,
     &             VAP2, GPKC, GPKH, GPKI, VTH (KO), SC(KO)
!
!
!                                    !end of printout

390   CONTINUE

      nrec=nrec+1
      write (19,rec=nrec) (W (KO), KO=1,nrectotal)
      nrec=nrec+1
      write (19,rec=nrec) (T (KO), KO=1,nrectotal)
      nrec=nrec+1
      write (19,rec=nrec) (TE(KO), KO=1,nrectotal)
      nrec=nrec+1
      write (19,rec=nrec) (QV(KO)*1000., KO=1,nrectotal)
      nrec=nrec+1
      write (19,rec=nrec) (QC(KO)*1000., KO=1,nrectotal)
      nrec=nrec+1
      write (19,rec=nrec) (QH(KO)*1000., KO=1,nrectotal)
      nrec=nrec+1
      write (19,rec=nrec) (QI(KO)*1000., KO=1,nrectotal)
      nrec=nrec+1
!   write (19,rec=nrec) (SC(KO), KO=1,nrectotal)
      write (19,rec=nrec) (QSAT(KO)*1000., KO=1,nrectotal)
      nrec=nrec+1
      write (19,rec=nrec) (QVENV(KO)*1000., KO=1,nrectotal)
!
      RETURN
!
! ************** FORMATS *********************************************
!
380   FORMAT(/,' Z(KM) P(MB) W(MPS) T(C)  T-TE   VAP   SAT   QC    QH'
     &               '     QI    VTH(MPS) SCAL'/)
!
400   FORMAT(1H , F4.1,F7.2,F7.2,F6.1,6F6.2,F7.2,1X,F6.2)
!
430   FORMAT(1H ,//I5,' MINUTES       DT= ',F6.2,' SECONDS   TIME= '
     &             ,F8.2,' SECONDS')
431   FORMAT(' ZTOP= ',F10.2)
!
      end subroutine printout
!
! *********************************************************************
      SUBROUTINE WATERBAL(qc,qh,qi,Qsat,qv,t,RHO,
     &                    cvi,est,WBAR,dqsdz,dt,L)
!
      implicit none
      integer :: L
      integer, parameter :: nkp = 200
      real,dimension(nkp) :: qc,qh,qi
      real,dimension(nkp) ::  Qsat,qv,t,RHO,cvi,est
      real :: WBAR,dqsdz,dt

!      print *,'WATERBAL'
      IF (QC (L) .LE.1.0E-10)then
        QC (L) = 0.  !DEFEAT UNDERFLOW PROBLEM
      endif
      IF (QH (L) .LE.1.0E-10) then
        QH (L) = 0.
      endif
      IF (QI (L) .LE.1.0E-10) then
        QI (L) = 0.
      endif
!
      CALL EVAPORATE(qsat,qv,qc,t,RHO,QH,qi,cvi,est,
     &           WBAR,dqsdz,dt,l)    !vapor to cloud,cloud to vapor
!
      CALL SUBLIMATE(T,qv,qsat,rho,qi,est,dt,l)    !vapor to ice
!
      CALL GLACIATE(Qh,qv,qsat,t,qi,dt,l)     !rain to ice

      CALL MELT(Qi,t,rho,cvi,qh,dt,L)         !ice to rain
!
!if(ak1 > 0. .or. ak2 > 0.)
      CALL CONVERT(T,rho,qh,qc,DT,L)  !(auto)conversion and accretion
!

      RETURN
      END SUBROUTINE WATERBAL
! *********************************************************************
      SUBROUTINE EVAPORATE(qsat,qv,qc,t,RHO,QH,qi,cvi,est,
     &           WBAR,dqsdz,dt,l)
!
!- evaporates cloud,rain and ice to saturation
!
      implicit none
!
!     XNO=10.0E06
!     HERC = 1.93*1.E-6*XN035        !evaporation constant
!
      real, PARAMETER :: HERC = 5.44E-4, CP = 1.004, HEATCOND = 2.5E3
      real, PARAMETER :: HEATSUBL = 2834., TMELT = 273., TFREEZE = 269.3

      real, PARAMETER :: FRC = HEATCOND / CP, SRC = HEATSUBL / CP

      real :: evhdt, evidt, evrate, evap, sd,
     &        quant, dividend, divisor, devidt

      integer, parameter :: nkp = 200
      real,dimension(nkp) ::  Qsat,qv,qc,t,RHO,QH,qi,cvi,est
      real :: WBAR,dqsdz,dt
      integer :: L
!
!
!      print *,'EVAPORATE'
      SD = QSAT (L) - QV (L)  !vapor deficit
      IF (SD.EQ.0.0)  RETURN
!IF (abs(SD).lt.1.e-7)  RETURN

      EVHDT = 0.
      EVIDT = 0.
!evrate =0.; evap=0.; sd=0.0; quant=0.0; dividend=0.0; divisor=0.0; devidt=0.0

      EVRATE = ABS (WBAR * DQSDZ)   !evaporation rate (Kessler 8.32)
      EVAP = EVRATE * DT            !what we can get in DT
      IF (SD.LE.0.0) THEN  !     condense. SD is negative
         IF (EVAP.GE.ABS (SD) ) THEN    !we get it all
            QC (L) = QC  (L) - SD  !deficit,remember?
            QV (L) = QSAT(L)       !set the vapor to saturation
            T  (L) = T   (L) - SD * FRC  !heat gained through condensation
                                !per gram of dry air
             RETURN
         ELSE
            QC (L) = QC (L) + EVAP         !get what we can in DT
            QV (L) = QV (L) - EVAP         !remove it from the vapor
            T  (L) = T  (L) + EVAP * FRC   !get some heat
            RETURN
         ENDIF
!
      ELSE                                !SD is positive, need some water
!
! not saturated. saturate if possible. use everything in order
! cloud, rain, ice. SD is positive

         IF (EVAP.LE.QC (L) ) THEN        !enough cloud to last DT
!
            IF (SD.LE.EVAP) THEN          !enough time to saturate
               QC (L) = QC (L) - SD       !remove cloud
               QV (L) = QSAT (L)          !saturate
               T (L) = T (L) - SD * FRC   !cool the parcel
               RETURN  !done
!
            ELSE   !not enough time

               SD = SD-EVAP               !use what there is
               QV (L) = QV (L) + EVAP     !add vapor
               T (L) = T (L) - EVAP * FRC !lose heat
               QC (L) = QC (L) - EVAP     !lose cloud
               !go on to rain.
            ENDIF
!
         ELSE                !not enough cloud to last DT
!
            IF (SD.LE.QC (L) ) THEN   !but there is enough to sat

               QV (L) = QSAT (L)  !use it
               QC (L) = QC (L) - SD
               T  (L) = T (L) - SD * FRC
               RETURN

            ELSE            !not enough to sat
               SD = SD-QC (L)
               QV (L) = QV (L) + QC (L)
               T  (L) = T (L) - QC (L) * FRC
               QC (L) = 0.0  !all gone

            ENDIF       !on to rain
         ENDIF          !finished with cloud
!
!  but still not saturated, so try to use some rain
!  this is tricky, because we only have time DT to evaporate. if there
!  is enough rain, we can evaporate it for dt. ice can also sublimate
!  at the same time. there is a compromise here.....use rain first, then
!  ice. saturation may not be possible in one DT time.
!  rain evaporation rate (W12),(OT25),(K Table 4). evaporate rain first
!  sd is still positive or we wouldn't be here.


         IF (QH (L) .LE.1.E-10) GOTO 33

!srf-25082005
!  QUANT = ( QC (L)  + QV (L) - QSAT (L) ) * RHO (L)   !g/m**3
         QUANT = ( QSAT (L)- QC (L) - QV (L)   ) * RHO (L)   !g/m**3
!
         EVHDT = (DT*HERC*(QUANT)*(QH (L)* RHO (L) ) **.65) / RHO (L)
!             rain evaporation in time DT

         IF (EVHDT.LE.QH (L) ) THEN           !enough rain to last DT

            IF (SD.LE.EVHDT) THEN       !enough time to saturate
               QH (L) = QH (L) - SD     !remove rain
               QV (L) = QSAT (L)        !saturate
               T (L) = T (L) - SD * FRC     !cool the parcel

               RETURN           !done
!
            ELSE                               !not enough time
               SD = SD-EVHDT         !use what there is
               QV (L) = QV (L) + EVHDT       !add vapor
               T (L) = T (L) - EVHDT * FRC       !lose heat
               QH (L) = QH (L) - EVHDT       !lose rain

            ENDIF                 !go on to ice.
!
         ELSE  !not enough rain to last DT
!
            IF (SD.LE.QH (L) ) THEN             !but there is enough to sat
               QV (L) = QSAT (L)                !use it
               QH (L) = QH (L) - SD
               T (L) = T (L) - SD * FRC
               RETURN
!
            ELSE                              !not enough to sat
               SD = SD-QH (L)
               QV (L) = QV (L) + QH (L)
               T (L) = T (L) - QH (L) * FRC
               QH (L) = 0.0                   !all gone

            ENDIF                             !on to ice
!

         ENDIF                                !finished with rain
!
!
!  now for ice
!  equation from (OT); correction factors for units applied
!
33       continue
         IF (QI (L) .LE.1.E-10) RETURN            !no ice there
!
         DIVIDEND = ( (1.E6 / RHO (L) ) **0.475) * (SD / QSAT (L)
     &            - 1) * (QI (L) **0.525) * 1.13
         DIVISOR = 7.E5 + 4.1E6 / (10. * EST (L) )

         DEVIDT = - CVI(L) * DIVIDEND / DIVISOR   !rate of change

         EVIDT = DEVIDT * DT                      !what we could get
!
! logic here is identical to rain. could get fancy and make subroutine
! but duplication of code is easier. God bless the screen editor.
!

         IF (EVIDT.LE.QI (L) ) THEN             !enough ice to last DT
!

            IF (SD.LE.EVIDT) THEN         !enough time to saturate
               QI (L) = QI (L) - SD       !remove ice
               QV (L) = QSAT (L)          !saturate
               T (L) = T (L) - SD * SRC       !cool the parcel

               RETURN             !done
!

            ELSE                                !not enough time

               SD = SD-EVIDT          !use what there is
               QV (L) = QV (L) + EVIDT        !add vapor
                T (L) =  T (L) - EVIDT * SRC      !lose heat
               QI (L) = QI (L) - EVIDT        !lose ice

            ENDIF                 !go on,unsatisfied
!
         ELSE                                   !not enough ice to last DT
!
            IF (SD.LE.QI (L) ) THEN             !but there is enough to sat

               QV (L) = QSAT (L)                !use it
               QI (L) = QI   (L) - SD
                T (L) =  T   (L) - SD * SRC

               RETURN
!
            ELSE                                 !not enough to sat
               SD = SD-QI (L)
               QV (L) = QV (L) + QI (L)
               T (L) = T (L) - QI (L) * SRC
               QI (L) = 0.0                      !all gone

            ENDIF                                !on to better things
                                           !finished with ice
         ENDIF
!
      ENDIF                                      !finished with the SD decision
!
      RETURN
!
      END SUBROUTINE EVAPORATE
!
! *********************************************************************
      SUBROUTINE CONVERT(T,rho,qh,qc,DT,L)
!
!- ACCRETION AND AUTOCONVERSION
!
!
      implicit none
      real,      PARAMETER ::  AK1 = 0.001    !conversion rate constant
      real,      PARAMETER ::  AK2 = 0.0052   !collection (accretion) rate
      real,      PARAMETER ::  TH  = 0.5      !Kessler threshold
      integer,   PARAMETER ::iconv = 1        !- Kessler conversion (=0)

!real, parameter :: ANBASE =  50.!*1.e+6 !Berry-number at cloud base #/m^3(maritime)
      real, parameter :: ANBASE =100000.!*1.e+6 !Berry-number at cloud base #/m^3(continental)
!real, parameter :: BDISP = 0.366       !Berry--size dispersion (maritime)
      real, parameter :: BDISP = 0.146       !Berry--size dispersion (continental)
      real, parameter :: TFREEZE = 269.3  !ice formation temperature
!
      real ::   accrete, con, q, h, bc1,   bc2,  total
      integer, parameter :: nkp = 200
      real,dimension(nkp) :: T,rho,qh,qc
      real :: DT
      integer :: L
!      print *,'CONVERT'
      IF (T (L)  .LE. TFREEZE) RETURN  !process not allowed above ice
!
      IF (QC (L) .EQ. 0.     ) RETURN

      ACCRETE = 0.
CON = 0.
      Q = RHO (L) * QC (L)
      H = RHO (L) * QH (L)
!
!     selection rules
!
!
      IF (QH (L) .GT. 0.     ) then
        ACCRETE = AK2 * Q * (H**.875)  !accretion, Kessler
      endif
!
      IF (ICONV.NE.0) THEN   !select Berry or Kessler
!
!old   BC1 = 120.
!old   BC2 = .0266 * ANBASE * 60.
!old   CON = BDISP * Q * Q * Q / (BC1 * Q * BDISP + BC2)

         CON = Q*Q*Q*BDISP/(60.*(5.*Q*BDISP+0.0366*ANBASE))
!
      ELSE
!
!   CON = AK1 * (Q - TH)   !Kessler autoconversion rate
!
!   IF (CON.LT.0.0) CON = 0.0   !havent reached threshold

         CON = max(0.,AK1 * (Q - TH)) ! versao otimizada
!
      ENDIF
!
!
      TOTAL = (CON + ACCRETE) * DT / RHO (L)

!
      IF (TOTAL.LT.QC (L) ) THEN
!
         QC (L) = QC (L) - TOTAL
         QH (L) = QH (L) + TOTAL    !no phase change involved
         RETURN
!
      ELSE
!
         QH (L) = QH (L) + QC (L)    !uses all there is
         QC (L) = 0.0
!
      ENDIF
!
      RETURN
!
      END SUBROUTINE CONVERT
!
!**********************************************************************
!
! ice - effect on temperature
!      TTD = 0.0
!      TTE = 0.0
!       CALL ICE(QSATW,QSATE,Y(1),Y(2),Y(3), &
!               TTA,TTB,TTC,DZ,ROH,D,C,TTD,TTE)
!       DYDX(1) = DYDX(1) + TTD  + TTE ! DT/DZ on Temp
!
!**********************************************************************
!
      SUBROUTINE SUBLIMATE( T,qv,qsat,rho,qi,est,dt,l)
!
! ********************* VAPOR TO ICE (USE EQUATION OT22)***************
!
      implicit none
      real, PARAMETER :: EPS = 0.622, HEATFUS = 334.,
     &                         HEATSUBL = 2834., CP = 1.004
      real, PARAMETER :: SRC = HEATSUBL / CP, FRC = HEATFUS / CP,
     &                         TMELT = 273.3
      real, PARAMETER :: TFREEZE = 269.3

      real ::dtsubh,  dividend,divisor, subl
      integer, parameter :: nkp = 200
      real,dimension(nkp) ::  T,qv,qsat,rho,qi,est
      real :: dt
      integer :: L
!
!      print *,'SUBLIMATE'
      DTSUBH = 0.
!
!selection criteria for sublimation
      IF (T (L)  .GT. TFREEZE  ) RETURN
      IF (QV (L) .LE. QSAT (L) ) RETURN
!
!     from (OT); correction factors for units applied
!
      DIVIDEND = ( (1.E6 / RHO (L) ) **0.475) * (QV (L) / QSAT (L)
     &           - 1) * (QI (L) **0.525) * 1.13
      DIVISOR = 7.E5 + 4.1E6 / (10. * EST (L) )
!

      DTSUBH = ABS (DIVIDEND / DIVISOR)   !sublimation rate
      SUBL = DTSUBH * DT                  !and amount possible
!
!     again check the possibilities
!
      IF (SUBL.LT.QV (L) ) THEN
!
         QV (L) = QV (L) - SUBL             !lose vapor
         QI (L) = QI (L) + SUBL           !gain ice
         T (L) = T (L) + SUBL * SRC         !energy change, warms air

         RETURN
!
      ELSE
!
         QI (L) = QV (L)                    !use what there is
         T  (L) = T (L) + QV (L) * SRC      !warm the air
         QV (L) = 0.0
!
      ENDIF
!
      RETURN
      END SUBROUTINE SUBLIMATE
!
! *********************************************************************
!
      SUBROUTINE GLACIATE(Qh,qv,qsat,t,qi,dt,l)
!
! *********************** CONVERSION OF RAIN TO ICE *******************
!     uses equation OT 16, simplest. correction from W not applied, but
!     vapor pressure differences are supplied.
!
!
      implicit none
      real, PARAMETER :: HEATFUS = 334., CP = 1.004,
     &                   EPS = 0.622, HEATSUBL = 2834.
      real, PARAMETER :: FRC = HEATFUS / CP, FRS = HEATSUBL / CP,
     &                   TFREEZE =  269.3
      real, PARAMETER :: GLCONST = 0.025   !glaciation time constant, 1/sec
      real dfrzh
      integer, parameter :: nkp = 200
      real,dimension(nkp) ::  Qh,qv,qsat,t,qi
      real :: DT
      integer :: L
!
!      print *,'GLACIATE'
      DFRZH = 0.    !rate of mass gain in ice
!
!selection rules for glaciation
      IF (QH (L) .LE. 0.       ) RETURN
      IF (QV (L) .LT. QSAT (L) ) RETURN
      IF (T  (L) .GT. TFREEZE  ) RETURN
!
!      NT=TMELT-T(L)
!      IF (NT.GT.50) NT=50
!

      DFRZH = DT * GLCONST * QH (L)    ! from OT(16)
!
      IF (DFRZH.LT.QH (L) ) THEN
!
         QI (L) = QI (L) + DFRZH
         QH (L) = QH (L) - DFRZH
         T (L) = T (L) + FRC * DFRZH  !warms air



         RETURN
!
      ELSE
!
         QI (L) = QI (L) + QH (L)
         T  (L) = T  (L) + FRC * QH (L)
         QH (L) = 0.0

!
      ENDIF
!
      RETURN
!
      END SUBROUTINE GLACIATE
!
!
! *********************************************************************
      SUBROUTINE MELT(Qi,t,rho,cvi,qh,dt,l)
!
! ******************* MAKES WATER OUT OF ICE **************************
!
      real, PARAMETER :: FRC = 332.27, TMELT = 273., F0 = 0.75   !ice velocity factor
      real DTMELT
      integer, parameter :: nkp = 200
      real,dimension(nkp) ::  Qi,t,rho,cvi,qh
      real :: DT
      integer :: L
!
!      print *, 'MELT'
      DTMELT = 0.   !conversion,ice to rain
!
!selection rules
      IF (QI (L) .LE. 0.0  ) RETURN
      IF (T (L)  .LT. TMELT) RETURN
!
                                                      !OT(23,24)
      DTMELT = DT * (2.27 / RHO (L) ) * CVI(L) * (T (L) -
     &               TMELT) * ( (RHO(L)
     &          * QI (L) * 1.E-6) **0.525) * (F0** ( - 0.42) )
                                                      !after Mason,1956
!
!     check the possibilities
!
      IF (DTMELT.LT.QI (L) ) THEN
!
         QH (L) = QH (L) + DTMELT
         QI (L) = QI (L) - DTMELT
         T  (L) = T (L) - FRC * DTMELT     !cools air


         RETURN
!
      ELSE
!
         QH (L) = QH (L) + QI (L)   !get all there is to get
         T  (L) = T (L) - FRC * QI (L)
         QI (L) = 0.0
!
      ENDIF
!
      RETURN
!
      END SUBROUTINE MELT

      SUBROUTINE htint (nzz1, vctra, eleva, nzz2, vctrb, elevb)
        IMPLICIT NONE
        INTEGER, INTENT(IN ) :: nzz1
        INTEGER, INTENT(IN ) :: nzz2
        REAL,    INTENT(IN ) :: vctra(nzz1)
        REAL,    INTENT(OUT) :: vctrb(nzz2)
        REAL,    INTENT(IN ) :: eleva(nzz1)
        REAL,    INTENT(IN ) :: elevb(nzz2)

        INTEGER :: l
        INTEGER :: k
        INTEGER :: kk
        REAL    :: wt
!        print *,'hint'
        l=1

        DO k=1,nzz2
           DO
              IF ( (elevb(k) <  eleva(1)) .OR.
     &             ((elevb(k) >= eleva(l)) .AND.
     &             (elevb(k) <= eleva(l+1))) ) THEN
                 wt = (elevb(k)-eleva(l))/(eleva(l+1)-eleva(l))
                 vctrb(k) = vctra(l)+(vctra(l+1)-vctra(l))*wt
                 EXIT
              ELSE IF ( elevb(k) >  eleva(nzz1))  THEN
                 wt= (elevb(k)-eleva(nzz1))/(eleva(nzz1-1)-eleva(nzz1))
                 vctrb(k) = vctra(nzz1)+(vctra(nzz1-1)-vctra(nzz1))*wt
                 EXIT
              END IF

              l=l+1
              IF(l == nzz1) THEN
                 PRINT *,'htint:nzz1',nzz1
                 DO kk=1,l
                    PRINT*,'kk,eleva(kk),elevb(kk)',eleva(kk),elevb(kk)
                 END DO
                 STOP 'htint'
              END IF
           END DO
        END DO
      END SUBROUTINE htint
!-----------------------------------------------------------------------------
      FUNCTION ESAT_PR (TEM)
!
! ******* Vapor Pressure  A.L. Buck JAM V.20 p.1527. (1981) ***********
!
      IMPLICIT NONE
      real, PARAMETER :: CI1 = 6.1115, CI2 = 22.542, CI3 = 273.48
      real, PARAMETER :: CW1 = 6.1121, CW2 = 18.729,
     &                   CW3 = 257.87, CW4 = 227.3
      real, PARAMETER :: TMELT = 273.3

      real ESAT_PR
      real temc , tem,esatm
!
!     formulae from Buck, A.L., JAM 20,1527-1532
!     custom takes esat wrt water always. formula for h2o only
!     good to -40C so:
!
!
!      print *,'ESAT_PR'
      TEMC = TEM - TMELT
      IF (TEMC.GT. - 40.0) GOTO 230
      ESATM = CI1 * EXP (CI2 * TEMC / (TEMC + CI3) )  !ice, millibars
      ESAT_PR = ESATM / 10. !kPa

      RETURN
!
230   ESATM = CW1*EXP ( ( (CW2 - (TEMC / CW4) ) * TEMC) / (TEMC + CW3))

      ESAT_PR = ESATM / 10. !kPa
      RETURN
      END function ESAT_PR
!     ******************************************************************