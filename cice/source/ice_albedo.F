c $Id: $
c=======================================================================
!---! the albedo parameterization
!---!
!---! authors:  Bruce P. Briegleb, NCAR 
!---!           Elizabeth C. Hunke, LANL
c=======================================================================

      module ice_albedo

      use ice_kinds_mod
      use ice_domain

      implicit none

      real (kind=dbl_kind), parameter ::
     &   albocn = 0.06_dbl_kind  ! ocean albedo

      ! weights for albedos match those for isccp shortwave forcing
      real (kind=dbl_kind), parameter ::             ! currently used only
     &   awtvdr = 0.29_dbl_kind  ! visible, direct   ! for history and
     &,  awtidr = 0.31_dbl_kind  ! near IR, direct   ! diagnostics
     &,  awtvdf = 0.24_dbl_kind  ! visible, diffuse
     &,  awtidf = 0.16_dbl_kind  ! near IR, diffuse

      ! parameter for fractional snow area 
      real (kind=dbl_kind), parameter ::
     &   snowpatch = 0.04_dbl_kind  

      ! albedos for ice in each category
      real (kind=dbl_kind) ::
     &   alvdrn (ilo:ihi,jlo:jhi,ncat) ! visible, direct   (fraction)
     &,  alidrn (ilo:ihi,jlo:jhi,ncat) ! near-ir, direct   (fraction)
     &,  alvdfn (ilo:ihi,jlo:jhi,ncat) ! visible, diffuse  (fraction)
     &,  alidfn (ilo:ihi,jlo:jhi,ncat) ! near-ir, diffuse  (fraction)

      ! albedos aggregated over categories
      real (kind=dbl_kind) ::
     &   alvdr (ilo:ihi,jlo:jhi)  ! visible, direct   (fraction)
     &,  alidr (ilo:ihi,jlo:jhi)  ! near-ir, direct   (fraction)
     &,  alvdf (ilo:ihi,jlo:jhi)  ! visible, diffuse  (fraction)
     &,  alidf (ilo:ihi,jlo:jhi)  ! near-ir, diffuse  (fraction)

c=======================================================================

      contains

c=======================================================================

      subroutine albedos

!---!-------------------------------------------------------------------
!---! compute albedos and aggregate them
!---! note: ice albedo is zero if no ice present
!---!-------------------------------------------------------------------

      use ice_constants
      use ice_grid
      use ice_state

      real (kind=dbl_kind), parameter :: 
     &   ahmax = p5    ! thickness above which ice albedo is constant, m
     &,  albicev   = 0.78_dbl_kind  ! visible ice albedo for h > ahmax
     &,  albicei   = 0.36_dbl_kind  ! near-ir ice albedo for h > ahmax
     &,  albsnowv = 0.98_dbl_kind   ! cold snow albedo, visible
     &,  albsnowi = 0.70_dbl_kind   ! cold snow albedo, near IR
     &,  dT_mlt = 5._dbl_kind  ! change in temp to give dalb_mlt albedo change
     &,  dalb_mlt = -0.075_dbl_kind  ! albedo change per dT_mlt change in temp
         ! albedo vis change per dT_mlt change in temp for snow
     &,  dalb_mltv = -0.100_dbl_kind
         ! albedo nir change per dT_mlt change in temp for snow 
     &,  dalb_mlti = -0.150_dbl_kind

      integer (kind=int_kind) :: i, j, nc
      real (kind=dbl_kind) :: 
     &   hi      ! ice thickness  (m)
     &,  hs      ! snow thickness (m)
     &,  albo    ! effective ocean albedo, function of ice thickness
     &,  asnow   ! snow-covered area fraction
     &,  asnwv   ! snow albedo, visible 
     &,  asnwi   ! snow albedo, near IR
     &,  fh      ! piecewise linear function of thickness 
     &,  fT      ! piecewise linear function of surface temperature
     &,  dTs     ! difference of Tsfc and Timelt

      do j=jlo,jhi
        do i=ilo,ihi
          if (tmask(i,j)) then
            do nc = 1,ncat
              if (aicen(i,j,nc).gt.puny) then
                hi = vicen(i,j,nc) / aicen(i,j,nc)
                hs = vsnon(i,j,nc) / aicen(i,j,nc)

                ! bare ice, thickness dependence
                fh = min(hi/ahmax,c1)
                albo = albocn*(c1-fh)
                alvdfn(i,j,nc) = albicev*fh + albo
                alidfn(i,j,nc) = albicei*fh + albo

                ! bare ice, temperature dependence
                dTs = Timelt - Tsfcn(i,j,nc)
                fT = min(dTs/dT_mlt-c1,c0)
                alvdfn(i,j,nc) = alvdfn(i,j,nc) - dalb_mlt*fT
                alidfn(i,j,nc) = alidfn(i,j,nc) - dalb_mlt*fT

                if( hs .gt. puny ) then
                  ! fractional area of snow on ice (thickness dependent)
                  asnow = hs / ( hs + snowpatch ) 
                  asnwv = albsnowv
                  asnwi = albsnowi
                  ! snow on ice, temperature dependence
                  asnwv = asnwv - dalb_mltv*fT
                  asnwi = asnwi - dalb_mlti*fT

                  ! combine ice and snow albedos
                  alvdfn(i,j,nc) = alvdfn(i,j,nc)*(c1-asnow) +
     $                             asnwv*asnow
                  alidfn(i,j,nc) = alidfn(i,j,nc)*(c1-asnow) +
     $                             asnwi*asnow
                endif
                alvdrn(i,j,nc) = alvdfn(i,j,nc)
                alidrn(i,j,nc) = alidfn(i,j,nc)

              else    ! no ice
                alvdfn(i,j,nc) = albocn
                alidfn(i,j,nc) = albocn
                alvdrn(i,j,nc) = albocn
                alidrn(i,j,nc) = albocn
              endif  ! aicen > puny
            enddo  ! nc

            ! aggregate
            alvdf(i,j) = c0
            alidf(i,j) = c0
            alvdr(i,j) = c0
            alidr(i,j) = c0
            if (aice(i,j) .gt. puny) then
              do nc = 1,ncat
                alvdf(i,j) = alvdf(i,j) + alvdfn(i,j,nc)*aicen(i,j,nc)
                alidf(i,j) = alidf(i,j) + alidfn(i,j,nc)*aicen(i,j,nc)
                alvdr(i,j) = alvdr(i,j) + alvdrn(i,j,nc)*aicen(i,j,nc)
                alidr(i,j) = alidr(i,j) + alidrn(i,j,nc)*aicen(i,j,nc)
              enddo
              alvdf(i,j) = alvdf(i,j) / aice(i,j)
              alidf(i,j) = alidf(i,j) / aice(i,j)
              alvdr(i,j) = alvdr(i,j) / aice(i,j)
              alidr(i,j) = alidr(i,j) / aice(i,j)
            endif  ! aicen > puny
          endif ! tmask
        enddo  
      enddo

      end subroutine albedos

c=======================================================================

      end module ice_albedo

c=======================================================================
