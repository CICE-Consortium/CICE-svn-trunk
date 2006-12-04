!=======================================================================
!BOP
!
! !MODULE: ice_ocean - ocean mixed layer internal to sea ice model
!
! !DESCRIPTION:
!
! Ocean mixed layer calculation (internal to sea ice model).
! Allows heat storage in ocean for uncoupled runs.
!
! !REVISION HISTORY:
!  SVN:$Id$
!
! authors:   John Weatherly, CRREL
!            C.M. Bitz, UW
!            Elizabeth C. Hunke, LANL
!            Bruce P. Briegleb, NCAR
!            William H. Lipscomb, LANL
!
! 2004: Block structure added by William Lipscomb
! 2005: Ocean-to-atmosphere fluxes added as 3D arrays, William Lipscomb
! 2006: Converted to free source form (F90) by Elizabeth Hunke
!
! !INTERFACE:
!
      module ice_ocean
!
! !USES:
!
      use ice_kinds_mod
      use ice_constants
!
!EOP
!
      implicit none
      save

      logical (kind=log_kind) :: &
         oceanmixed_ice           ! if true, use ocean mixed layer

      real (kind=dbl_kind), parameter :: &
         cprho = cp_ocn*rhow

!=======================================================================

      contains

!=======================================================================
!BOP
!
! !ROUTINE: ocean_mixed_layer - compute SST and freeze/melt potential
!
! !DESCRIPTION:
!
! Compute the mixed layer heat balance and update the SST.
! Compute the energy available to freeze or melt ice.
! NOTE: SST changes due to fluxes through the ice are computed in
!       ice_therm_vertical.
!
! !REVISION HISTORY: same as module
!
! !INTERFACE:
!
      subroutine ocean_mixed_layer (dt)
!
! !USES:
!
      use ice_blocks
      use ice_domain
      use ice_state, only: aice
      use ice_flux
      use ice_grid, only: tmask
      use ice_atmo, only: atmo_boundary_layer
!
! !INPUT/OUTPUT PARAMETERS:
!
      real (kind=dbl_kind), intent(in) :: &
         dt      ! time step
!
!EOP
!
      integer (kind=int_kind) :: &
         i, j           , & ! horizontal indices
         iblk           , & ! block index
         ilo,ihi,jlo,jhi    ! beginning and end of physical domain

      real (kind=dbl_kind), dimension(nx_block,ny_block) :: &
         delt  , & ! potential temperature difference   (K)
         delq  , & ! specific humidity difference   (kg/kg)
         shcoef, & ! transfer coefficient for sensible heat
         lhcoef    ! transfer coefficient for latent heat

      integer (kind=int_kind) :: &
         ij        ! combined ij index

      integer (kind=int_kind), save :: &
         icells    ! number of ocean cells

      integer (kind=int_kind), dimension(nx_block*ny_block), save :: &
         indxi, indxj    ! compressed indices for ocean cells

      type (block) :: &
         this_block           ! block information for current block

      do iblk = 1, nblocks

      !-----------------------------------------------------------------
      ! Identify ocean cells.
      ! Set fluxes to zero in land cells.
      !-----------------------------------------------------------------

         icells = 0
         do j = 1, ny_block
         do i = 1, nx_block
            if (tmask(i,j,iblk)) then
               icells = icells + 1
               indxi(icells) = i
               indxj(icells) = j
            else
               sst       (i,j,iblk) = c0
               frzmlt    (i,j,iblk) = c0
               flwout_ocn(i,j,iblk) = c0
               fsens_ocn (i,j,iblk) = c0
               flat_ocn  (i,j,iblk) = c0
               evap_ocn  (i,j,iblk) = c0
            endif
         enddo                  ! i
         enddo                  ! j

      !-----------------------------------------------------------------
      ! Compute boundary layer quantities
      !-----------------------------------------------------------------
         call atmo_boundary_layer (nx_block,  ny_block,   &
                                   'ocn',     icells,     &
                                   indxi,     indxj,      &
                                   sst        (:,:,iblk), &    
                                   potT       (:,:,iblk), &
                                   uatm       (:,:,iblk), &   
                                   vatm       (:,:,iblk), &   
                                   wind       (:,:,iblk), &   
                                   zlvl       (:,:,iblk), &   
                                   Qa         (:,:,iblk), &     
                                   rhoa       (:,:,iblk), &
                                   strairx_ocn(:,:,iblk), & 
                                   strairy_ocn(:,:,iblk), & 
                                   Tref_ocn   (:,:,iblk), & 
                                   Qref_ocn   (:,:,iblk), & 
                                   delt       (:,:),      &    
                                   delq       (:,:),      &
                                   lhcoef     (:,:),      &
                                   shcoef     (:,:) )

      !-----------------------------------------------------------------
      ! Ocean  albedo
      ! For now, assume albedo = albocn in each spectral band.
      !-----------------------------------------------------------------

         alvdr_ocn(:,:,iblk) = albocn
         alidr_ocn(:,:,iblk) = albocn
         alvdf_ocn(:,:,iblk) = albocn
         alidf_ocn(:,:,iblk) = albocn

      !-----------------------------------------------------------------
      ! Compute ocean fluxes and update SST
      !-----------------------------------------------------------------


         call ocean_energy_budget (nx_block, ny_block,   &
                                   dt,       icells,     &
                                   indxi,    indxj,      &
                                   delt      (:,:),      &   
                                   delq      (:,:),      &
                                   lhcoef    (:,:),      &
                                   shcoef    (:,:),      &
                                   aice      (:,:,iblk), &
                                   Tf        (:,:,iblk), &
                                   swvdr     (:,:,iblk), &
                                   swidr     (:,:,iblk), &
                                   swvdf     (:,:,iblk), &
                                   swidf     (:,:,iblk), &    
                                   alvdr_ocn (:,:,iblk), &
                                   alidr_ocn (:,:,iblk), &
                                   alvdf_ocn (:,:,iblk), &
                                   alidf_ocn (:,:,iblk), &
                                   flw       (:,:,iblk), &
                                   qdp       (:,:,iblk), &
                                   hmix      (:,:,iblk), &
                                   flwout_ocn(:,:,iblk), &
                                   fsens_ocn (:,:,iblk), &
                                   flat_ocn  (:,:,iblk), &
                                   evap_ocn  (:,:,iblk), &
                                   sst       (:,:,iblk), &
                                   frzmlt    (:,:,iblk))

         evap_ocn(:,:,iblk) = flat_ocn(:,:,iblk) / Lvap

      enddo                     ! iblk

      end subroutine ocean_mixed_layer

!=======================================================================
!BOP
!
! !IROUTINE: ocean_energy_budget - sfc energy balance, sst change, frzmlt
!
! !INTERFACE:
!
      subroutine ocean_energy_budget (nx_block,   ny_block,  &
                                      dt,         icells,    &
                                      indxi,      indxj,     &
                                      delt,       delq,      &
                                      lhcoef,     shcoef,    &
                                      aice,       Tf,        &
                                      swvdr,      swidr,     &
                                      swvdf,      swidf,     &
                                      alvdr_ocn,  alidr_ocn, &
                                      alvdf_ocn,  alidf_ocn, &
                                      flw,                   &
                                      qdp,        hmix,      &
                                      flwout_ocn, fsens_ocn, &
                                      flat_ocn,   evap_ocn,  &
                                      sst,        frzmlt)

! !DESCRIPTION:
!
! Compute ocean energy budget and update SST accordingly.
! Compute freeze-melt potential.
! 
! !REVISION HISTORY: same as module
!
! !USES:
!
! !INPUT/OUTPUT PARAMETERS:
!
      integer (kind=int_kind), intent(in) :: &
         nx_block, ny_block, & ! block dimensions
         icells                ! number of cells that require atmo fluxes

      integer (kind=int_kind), dimension(nx_block*ny_block), &
         intent(in) :: &
         indxi, indxj    ! compressed i and j indices

      real (kind=dbl_kind), intent(in) :: &
         dt              ! time step

      real (kind=dbl_kind), dimension(nx_block,ny_block), &
         intent(in) :: &
         delt         , & ! potential T difference   (K)
         delq         , & ! humidity difference      (kg/kg)
         shcoef       , & ! transfer coefficient for sensible heat
         lhcoef       , & ! transfer coefficient for latent heat
         aice         , & ! fractional ice area
         Tf           , & ! ocean freezing temperature (C)
         swvdr        , & ! incoming shortwave, visible direct (W/m^2)
         swidr        , & ! incoming shortwave, near IR direct (W/m^2)
         swvdf        , & ! incoming shortwave, visible diffuse (W/m^2)
         swidf        , & ! incoming shortwave, near IR diffuse (W/m^2)
         alvdr_ocn    , & ! visible albedo, direct   (fraction)
         alidr_ocn    , & ! near-ir albedo, direct   (fraction)
         alvdf_ocn    , & ! visible albedo, diffuse  (fraction)
         alidf_ocn    , & ! near-ir albedo, diffuse  (fraction)
         flw          , & ! incoming longwave (W/m^2)
         hmix             ! ocean mixed layer depth (m)

      real (kind=dbl_kind), dimension(nx_block,ny_block), &
         intent(inout) :: &
         sst          , & ! sea surface temperature (C)
         qdp          , & ! deep ocean heat flux (W/m^2)
         frzmlt       , & ! freeze-melt potential (W/m^2)
         fsens_ocn    , & ! sensible heat flux (W/m^2)
         flat_ocn     , & ! latent heat flux   (W/m^2)
         flwout_ocn   , & ! outgoing longwave  (W/m^2)
         evap_ocn         ! evaporative vapor flux (kg/m^2/s)
!
!EOP
!
      real (kind=dbl_kind) :: &
         TsfK , & ! surface temperature (K)
         swabs    ! surface absorbed shortwave heat flux (W/m^2)

      real (kind=dbl_kind), parameter :: &
         frzmlt_max = c1000   ! max magnitude of frzmlt (W/m^2)

      integer (kind=int_kind) :: i, j, ij

!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
      do ij = 1, icells
         i = indxi(ij)
         j = indxj(ij)

         ! shortwave radiative flux
         swabs = (c1-alvdr_ocn(i,j)) * swvdr(i,j) &
               + (c1-alidr_ocn(i,j)) * swidr(i,j) &
               + (c1-alvdf_ocn(i,j)) * swvdf(i,j) &
               + (c1-alidf_ocn(i,j)) * swidf(i,j) 

         ! ocean surface temperature in Kelvin
         TsfK = sst(i,j) + Tffresh

         ! longwave radiative flux
         flwout_ocn(i,j) = -stefan_boltzmann * TsfK**4

         ! downward latent and sensible heat fluxes
         fsens_ocn(i,j) =  shcoef(i,j) * delt(i,j)
         flat_ocn (i,j) =  lhcoef(i,j) * delq(i,j)
         evap_ocn (i,j) = -flat_ocn(i,j) / Lvap

         ! Compute sst change due to exchange with atm/ice above
         ! Note: fhnet, fswthru are added in ice_therm_vertical.F
         sst(i,j) = sst(i,j) + &
              (fsens_ocn(i,j) + flat_ocn(i,j) + flwout_ocn(i,j) &
             + flw(i,j) + swabs) * (c1-aice(i,j)) * dt &
             / (cprho*hmix(i,j))

         ! adjust qdp if cooling of mixed layer would occur when sst <= Tf
         if (sst(i,j) <= Tf(i,j) .and. qdp(i,j) > c0) qdp(i,j) = c0

         ! computed T change due to exchange with deep layers:
         sst(i,j) = sst(i,j) - qdp(i,j)*dt/(cprho*hmix(i,j))

         ! compute potential to freeze or melt ice
         frzmlt(i,j) = (Tf(i,j)-sst(i,j))*cprho*hmix(i,j)/dt
         frzmlt(i,j) = min(max(frzmlt(i,j),-frzmlt_max),frzmlt_max)

         ! if sst is below freezing, reset sst to Tf
         if (sst(i,j) <= Tf(i,j)) sst(i,j) = Tf(i,j)

      enddo                     ! ij


      end subroutine ocean_energy_budget

!=======================================================================

      end module ice_ocean

!=======================================================================
