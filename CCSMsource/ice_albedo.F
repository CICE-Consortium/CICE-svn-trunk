c $Id: ice_albedo.F,v 1.11 2006/06/08 22:16:07 eclare Exp $
!=======================================================================
!BOP
!
! !MODULE: ice_albedo - snow and ice albedo parameterization
!
! !DESCRIPTION:
!
! The albedo parameterization
!
! !REVISION HISTORY:
!
! authors:  Bruce P. Briegleb, NCAR 
!           Elizabeth C. Hunke, LANL
!
! Vectorized by Clifford Chen (Fujitsu) and William H. Lipscomb (LANL)
!
! !INTERFACE:
!
      module ice_albedo
!
! !USES:
!
      use ice_kinds_mod
      use ice_domain
!
!EOP
!
      implicit none
      save

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
     &   snowpatch = 0.02_dbl_kind

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

      ! baseline albedos
      real (kind=dbl_kind) ::
     &   albicev   ! visible ice albedo for h > ahmax
     &,  albicei   ! near-ir ice albedo for h > ahmax
     &,  albsnowv  ! cold snow albedo, visible
     &,  albsnowi  ! cold snow albedo, near IR

!=======================================================================

      contains

!=======================================================================
!BOP
!
! !IROUTINE: albedos - compute snow/ice albedos and aggregate
!
! !INTERFACE:
!
      subroutine albedos
!
! !DESCRIPTION:
!
! Compute albedos and aggregate them \\
! note: ice albedo is zero if no ice present
!
! !REVISION HISTORY:
!
! authors:  Bruce P. Briegleb, NCAR 
!           Elizabeth C. Hunke, LANL
!
! !USES:
!
      use ice_constants
      use ice_grid
      use ice_state
!
! !INPUT/OUTPUT PARAMETERS:
!
!EOP
!
      real (kind=dbl_kind), parameter :: 
     &   ahmax     = 0.5_dbl_kind    ! thickness above which ice albedo constant (m)
     &,  dT_mlt    = 1._dbl_kind     ! change in temp to give dalb_mlt albedo change
     &,  dalb_mlt  = -0.075_dbl_kind ! albedo change per dT_mlt change 
                                     ! in temp for ice
     &,  dalb_mltv = -0.100_dbl_kind ! albedo vis change per dT_mlt change 
                                     ! in temp for snow
     &,  dalb_mlti = -0.150_dbl_kind ! albedo nir change per dT_mlt change 
                                     ! in temp for snow

      integer (kind=int_kind) :: i, j, n

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
     &,  fhtan   ! factor used in albedo dependence on ice thickness

      integer (kind=int_kind) ::
     &   icells  ! number of ice/ocean grid cells
     &,  ij      ! horizontal index, combines i and j loops

      integer (kind=int_kind), dimension (1:(ihi-ilo+1)*(jhi-jlo+1)) ::
     &   indxi   ! compressed indices for ice/ocean cells
     &,  indxj

      fhtan = atan(ahmax*c4)

      icells = 0
      do j = jlo, jhi
      do i = ilo, ihi
         if (tmask(i,j)) then
            icells = icells + 1
            indxi(icells) = i
            indxj(icells) = j
         endif                  ! tmask
      enddo
      enddo

      !-----------------------------------------------------------------
      ! albedo for each thickness category
      !-----------------------------------------------------------------

      do n = 1, ncat
!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
         do ij = 1, icells
            i = indxi(ij)
            j = indxj(ij)

            if (aicen(i,j,n) > puny) then
               hi = vicen(i,j,n) / aicen(i,j,n)
               hs = vsnon(i,j,n) / aicen(i,j,n)

               ! bare ice, thickness dependence
               fh = min(atan(hi*c4)/fhtan,c1)
               albo = albocn*(c1-fh)
               alvdfn(i,j,n) = albicev*fh + albo
               alidfn(i,j,n) = albicei*fh + albo

               ! bare ice, temperature dependence
               dTs = Timelt - Tsfcn(i,j,n)
               fT = min(dTs/dT_mlt-c1,c0)
               alvdfn(i,j,n) = alvdfn(i,j,n) - dalb_mlt*fT
               alidfn(i,j,n) = alidfn(i,j,n) - dalb_mlt*fT

               ! avoid negative albedos for thin, bare, melting ice
               alvdfn(i,j,n) = max (alvdfn(i,j,n), albocn)
               alidfn(i,j,n) = max (alidfn(i,j,n), albocn)

               if( hs > puny ) then

                  ! fractional area of snow on ice (thickness dependent)
                  asnow = hs / ( hs + snowpatch ) 
                  asnwv = albsnowv
                  asnwi = albsnowi

                  ! snow on ice, temperature dependence
                  asnwv = asnwv - dalb_mltv*fT
                  asnwi = asnwi - dalb_mlti*fT

                  ! combine ice and snow albedos
                  alvdfn(i,j,n) = alvdfn(i,j,n)*(c1-asnow) +
     &                            asnwv*asnow
                  alidfn(i,j,n) = alidfn(i,j,n)*(c1-asnow) +
     &                            asnwi*asnow
               endif            ! hs > puny

               alvdrn(i,j,n) = alvdfn(i,j,n)
               alidrn(i,j,n) = alidfn(i,j,n)

            else                ! no ice
               alvdfn(i,j,n) = albocn
               alidfn(i,j,n) = albocn
               alvdrn(i,j,n) = albocn
               alidrn(i,j,n) = albocn
            endif               ! aicen > puny
         enddo                  ! ij
      enddo                     ! ncat

      !-----------------------------------------------------------------
      ! aggregate
      !-----------------------------------------------------------------

      do j = jlo, jhi
      do i = ilo, ihi
         alvdf(i,j) = c0
         alidf(i,j) = c0
         alvdr(i,j) = c0
         alidr(i,j) = c0
      enddo
      enddo

      do n = 1, ncat
!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
         do ij = 1, icells
            i = indxi(ij)
            j = indxj(ij)
            if (aice(i,j) > puny) then
               alvdf(i,j) = alvdf(i,j) + alvdfn(i,j,n)*aicen(i,j,n)
               alidf(i,j) = alidf(i,j) + alidfn(i,j,n)*aicen(i,j,n)
               alvdr(i,j) = alvdr(i,j) + alvdrn(i,j,n)*aicen(i,j,n)
               alidr(i,j) = alidr(i,j) + alidrn(i,j,n)*aicen(i,j,n)
            endif
         enddo                  ! ij
      enddo                     ! ncat

!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
      do ij = 1, icells
         i = indxi(ij)
         j = indxj(ij)
         if (aice(i,j) > puny) then
            alvdf(i,j) = alvdf(i,j) / aice(i,j)
            alidf(i,j) = alidf(i,j) / aice(i,j)
            alvdr(i,j) = alvdr(i,j) / aice(i,j)
            alidr(i,j) = alidr(i,j) / aice(i,j)
         endif                  ! aicen > puny
      enddo                     ! ij

      end subroutine albedos

!=======================================================================

      end module ice_albedo

!=======================================================================
