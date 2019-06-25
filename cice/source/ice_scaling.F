c $Id: ice_scaling.F,v 1.2 2001/02/23 21:06:56 schramm Exp $
c=======================================================================

      module ice_scaling

      use ice_domain
      use ice_kinds_mod
      use ice_constants
      use ice_state

      implicit none

c=======================================================================

      contains

c=======================================================================

      subroutine scale_fluxes

      use ice_flux
      use ice_albedo
      use ice_history
      use ice_grid

      real (kind=dbl_kind) :: ar

      integer (kind=int_kind) :: i, j

      do j=jlo,jhi
       do i=ilo,ihi
        if (tmask(i,j).and.(aice(i,j).gt.0)) then
          ar = c1/aice(i,j)
          Fsensible(i,j) = Fsensible(i,j) *ar
          Flatent(i,j)   = Flatent(i,j)   *ar
          Flwout(i,j)    = Flwout(i,j)    *ar
          sabs(i,j)      = sabs(i,j)      *ar ! not needed for coupling
          Fswthru(i,j)   = Fswthru(i,j)   *ar
          evap(i,j)      = evap(i,j)      *ar
          Fhnet(i,j)     = Fhnet(i,j)     *ar
          Fsalt(i,j)     = Fsalt(i,j)     *ar
          Fresh(i,j)     = Fresh(i,j)     *ar
         else
          Fsensible(i,j) = c0
          Flatent(i,j)   = c0
          Flwout(i,j)    = -stefan_boltzmann*(Tffresh - Tf(i,j))**4 
                                  ! to make longwave up over ocean
                                  ! reasonable for history file
          Fswthru(i,j)   = c0     
          sabs(i,j)      = c0
          alvdf(i,j)     = c0
          alidf(i,j)     = c0
          alvdr(i,j)     = c0
          alidr(i,j)     = c0
          evap(i,j)      = c0
          Fhnet(i,j)     = c0
          Fsalt(i,j)     = c0
          Fresh(i,j)     = c0
         endif 
       enddo  
      enddo

      end subroutine scale_fluxes

c=======================================================================

      subroutine scale_hist_fluxes

      use ice_flux
      use ice_albedo
      use ice_history
      use ice_grid

      real (kind=dbl_kind) :: ar

      integer (kind=int_kind) :: i, j

      do j=jlo,jhi
       do i=ilo,ihi
	! use aice_init for correct diagnostics 
        if (tmask(i,j).and.(aice_init(i,j).gt.0)) then 	
          ar = c1/aice_init(i,j)	
          Frsh_hist(i,j)     = Frsh_hist(i,j)*ar
          Foht_hist(i,j)     = Foht_hist(i,j)*ar
          Fswthru_hist(i,j)  = Fswthru_hist(i,j)*ar
         else
          Foht_hist(i,j)     = c0
          Frsh_hist(i,j)     = c0
          Fswthru_hist(i,j)  = c0
         endif 
       enddo  
      enddo

      end subroutine scale_hist_fluxes

c=======================================================================

      end module ice_scaling

c=======================================================================
