!=======================================================================
!BOP
!
! !MODULE: ice_transport_upwind - upwind horizontal advection
!
! !DESCRIPTION:
!
! Calculates horizontal advection using an upwind algorithm
!
! !REVISION HISTORY:
!  SVN:$Id$
!
! authors: Elizabeth C. Hunke and William H. Lipscomb, LANL 
!
! 2004: Revised by William Lipscomb from ice_transport_mpdata.
!       Stripped out mpdata, retained upwind, and added block structure.
!   
!
! !INTERFACE:

      module ice_transport_upwind
!
! !USES:
!
      use ice_domain_size
      use ice_constants
      use ice_fileunits
      use ice_work, only:  worka, workb
!
!EOP
!
      implicit none

!=======================================================================

      contains

!=======================================================================
!BOP
!
! !IROUTINE: transport - computes transport equations for one timestep
!
! !INTERFACE:
!
      subroutine transport_upwind (dt)
!
! !DESCRIPTION:
!
! Computes the transport equations for one timestep using upwind. Sets
! several fields into a work array and passes it to upwind routine.
!
! !REVISION HISTORY:
!
! same as module
!
! !USES:
!
      use ice_boundary
      use ice_blocks
      use ice_domain
      use ice_state
      use ice_grid, only: HTE, HTN, tarea
      use ice_timers
!
! !INPUT/OUTPUT PARAMETERS:
!
      real (kind=dbl_kind), intent(in) ::
     &   dt      ! time step
!
!EOP
!
      integer (kind=int_kind), parameter ::
     &   narr = 1 + ncat*(3+ntrcr)   ! number of state variable arrays
                                     ! not including eicen, esnon

      integer (kind=int_kind) ::
     &   i, j, iblk         ! horizontal indices
     &,  ilo,ihi,jlo,jhi    ! beginning and end of physical domain

      real (kind=dbl_kind), dimension (nx_block,ny_block,nblocks) ::
     &   uee, vnn           ! cell edge velocities

      real (kind=dbl_kind),
     &   dimension (nx_block,ny_block,narr,max_blocks) :: 
     &   works              ! work array

      type (block) ::
     &   this_block           ! block information for current block

      call ice_timer_start(timer_advect)  ! advection 

    !-------------------------------------------------------------------
    ! Get ghost cell values of state variables.
    ! (Assume velocities are already known for ghost cells.)
    !-------------------------------------------------------------------

      call bound_state (aicen, trcrn,
     &                  vicen, vsnon, 
     &                  eicen, esnon)

    !-------------------------------------------------------------------
    ! Average corner velocities to edges.
    !-------------------------------------------------------------------
      
      do iblk = 1, nblocks
         this_block = get_block(blocks(iblk),iblk)         
         ilo = this_block%ilo
         ihi = this_block%ihi
         jlo = this_block%jlo
         jhi = this_block%jhi

         do j = jlo, jhi
         do i = ilo, ihi
            uee(i,j,iblk) = p5*(uvel(i,j,iblk) + uvel(i,j-1,iblk))
            vnn(i,j,iblk) = p5*(vvel(i,j,iblk) + vvel(i-1,j,iblk))
         enddo
         enddo
      enddo

      call ice_timer_start(timer_bound)
      call update_ghost_cells (uee,             bndy_info,
     &                         field_loc_Eface, field_type_scalar)
      call update_ghost_cells (vnn,             bndy_info,
     &                         field_loc_Nface, field_type_scalar)
      call ice_timer_stop(timer_bound)

      do iblk = 1, nblocks

      !-----------------------------------------------------------------
      ! fill work arrays with fields to be advected
      !-----------------------------------------------------------------

         call state_to_work (nx_block,             ny_block,
     &                       narr,                 trcr_depend,
     &                       aicen (:,:,  :,iblk), trcrn (:,:,:,:,iblk),
     &                       vicen (:,:,  :,iblk), vsnon (:,:,  :,iblk),
     &                       aice0 (:,:,    iblk), works (:,:,  :,iblk))

      !-----------------------------------------------------------------
      ! advect
      !-----------------------------------------------------------------

         call upwind_field (nx_block,       ny_block,
     &                      nghost,         dt,
     &                      narr,           works(:,:,:,iblk),  
     &                      uee(:,:,iblk),  vnn    (:,:,iblk), 
     &                      HTE(:,:,iblk),  HTN    (:,:,iblk), 
     &                      tarea(:,:,iblk))

         call upwind_field (nx_block,       ny_block,
     &                      nghost,         dt,
     &                      ntilyr,         eicen(:,:,:,iblk), 
     &                      uee(:,:,iblk),  vnn    (:,:,iblk), 
     &                      HTE(:,:,iblk),  HTN    (:,:,iblk), 
     &                      tarea(:,:,iblk))

         call upwind_field (nx_block,       ny_block,
     &                      nghost,         dt,
     &                      ntslyr,         esnon(:,:,:,iblk), 
     &                      uee(:,:,iblk),  vnn    (:,:,iblk), 
     &                      HTE(:,:,iblk),  HTN    (:,:,iblk), 
     &                      tarea(:,:,iblk))

      !-----------------------------------------------------------------
      ! convert work arrays back to state variables
      !-----------------------------------------------------------------

         call work_to_state (nx_block,            ny_block,
     &                       narr,                trcr_depend,
     &                       aicen(:,:,  :,iblk), trcrn (:,:,:,:,iblk),
     &                       vicen(:,:,  :,iblk), vsnon (:,:,  :,iblk),
     &                       aice0(:,:,    iblk), works (:,:,  :,iblk)) 

      enddo                     ! iblk
 
      call ice_timer_stop(timer_advect)  ! advection 

      end subroutine transport_upwind

!=======================================================================
!BOP
!
! !IROUTINE: state_to_work - fill work arrays with state variables
!
! !INTERFACE:
!
      subroutine state_to_work (nx_block, ny_block,
     &                          narr,     trcr_depend,
     &                          aicen,    trcrn,
     &                          vicen,    vsnon,
     &                          aice0,    works)

!
! !DESCRIPTION:
!
! Fill work array with state variables in preparation for upwind transport
!
! !REVISION HISTORY:
!
! same as module
!
! !USES:
!
      use ice_itd, only: ilyr1, slyr1
!
! !INPUT/OUTPUT PARAMETERS:
!
      integer (kind=int_kind), intent(in) ::
     &   nx_block, ny_block   ! block dimensions
     &,  narr        ! number of 2D state variable arrays in works array

      integer (kind=int_kind), dimension (ntrcr), intent(in) ::
     &   trcr_depend ! = 0 for aicen tracers, 1 for vicen, 2 for vsnon

      real (kind=dbl_kind), dimension (nx_block,ny_block,ncat),
     &   intent(in) ::
     &   aicen     ! concentration of ice
     &,  vicen     ! volume per unit area of ice          (m)
     &,  vsnon     ! volume per unit area of snow         (m)

      real (kind=dbl_kind), dimension (nx_block,ny_block,ntrcr,ncat),
     &   intent(in) ::
     &   trcrn     ! ice tracers

      real (kind=dbl_kind), dimension (nx_block,ny_block),
     &   intent(in) ::
     &   aice0     ! concentration of open water

      real (kind=dbl_kind), dimension(nx_block,ny_block,narr),
     &   intent (out) :: 
     &   works     ! work array
!
!EOP
!
      integer (kind=int_kind) :: 
     &   i, j, k, n, it   ! counting indices
     &,  narrays          ! counter for number of state variable arrays

      !-----------------------------------------------------------------
      ! This array is used for performance (balance memory/cache vs
      ! number of bound calls);  a different number of arrays may perform
      ! better depending on the machine used, number of processors, etc.
      ! --tested on SGI R2000, using 4 pes for the ice model under MPI
      !-----------------------------------------------------------------

      do j = 1, ny_block
      do i = 1, nx_block
         works(i,j,1) = aice0(i,j)
      enddo
      enddo
      narrays = 1

      do n=1, ncat

         do j = 1, ny_block
         do i = 1, nx_block
            works(i,j,narrays+1) = aicen(i,j,n)
            works(i,j,narrays+2) = vicen(i,j,n)
            works(i,j,narrays+3) = vsnon(i,j,n)
         enddo                  ! i
         enddo                  ! j
         narrays = narrays + 3

         do it = 1, ntrcr
            if (trcr_depend(it) == 0) then
               do j = 1, ny_block
               do i = 1, nx_block
                  works(i,j,narrays+it) = aicen(i,j,n)*trcrn(i,j,it,n)
               enddo
               enddo
            elseif (trcr_depend(it) ==1) then
               do j = 1, ny_block
               do i = 1, nx_block
                  works(i,j,narrays+it) = vicen(i,j,n)*trcrn(i,j,it,n)
               enddo
               enddo
            elseif (trcr_depend(it) ==1) then
               do j = 1, ny_block
               do i = 1, nx_block
                  works(i,j,narrays+it) = vsnon(i,j,n)*trcrn(i,j,it,n)
               enddo
               enddo
            endif
         enddo
         narrays = narrays + ntrcr

      enddo                     ! n

      if (narr /= narrays) write(nu_diag,*) 
     &     "Wrong number of arrays in transport bound call"

      end subroutine state_to_work

!=======================================================================
!BOP
!
! !IROUTINE: work_to_state - convert work arrays back to state variables
!
! !INTERFACE:
!
      subroutine work_to_state (nx_block, ny_block,
     &                          narr,     trcr_depend,
     &                          aicen,    trcrn,
     &                          vicen,    vsnon,
     &                          aice0,    works)

!
! !DESCRIPTION:
!
! Convert work array back to state variables
!
! !REVISION HISTORY:
!
! same as module
!
! !USES:
!
      use ice_itd, only: ilyr1, slyr1, compute_tracers
!
! !INPUT/OUTPUT PARAMETERS:
!
      integer (kind=int_kind), intent (in) ::
     &   nx_block, ny_block   ! block dimensions
     &,  narr        ! number of 2D state variable arrays in works array

      integer (kind=int_kind), dimension (ntrcr), intent(in) ::
     &   trcr_depend ! = 0 for aicen tracers, 1 for vicen, 2 for vsnon

      real (kind=dbl_kind), intent (in) :: 
     &   works (nx_block,ny_block,narr)

      real (kind=dbl_kind), dimension (nx_block,ny_block,ncat),
     &   intent(out) ::
     &   aicen     ! concentration of ice
     &,  vicen     ! volume per unit area of ice          (m)
     &,  vsnon     ! volume per unit area of snow         (m)

      real (kind=dbl_kind), dimension (nx_block,ny_block,ntrcr,ncat),
     &   intent(out) ::
     &   trcrn     ! ice tracers

      real (kind=dbl_kind), dimension (nx_block,ny_block),
     &   intent(out) ::
     &   aice0     ! concentration of open water
!
!EOP
!
      integer (kind=int_kind) :: 
     &   i, j, k, n , it  ! counting indices
     &,  narrays          ! counter for number of state variable arrays
     &,  icells           ! number of ocean/ice cells

      integer (kind=int_kind), dimension (nx_block*ny_block) ::
     &  indxi, indxj

      real (kind=dbl_kind), dimension (nx_block*ny_block,narr) :: 
     &   work 

      ! for call to compute_tracers
      icells = 0
      do j = 1, ny_block
      do i = 1, nx_block
         icells = icells + 1
         indxi(icells) = i
         indxj(icells) = j
         work (icells,:) = works(i,j,:)
      enddo
      enddo

      do j=1,ny_block
      do i=1,nx_block
         aice0(i,j) = works(i,j,1)
      enddo
      enddo
      narrays = 1               ! aice0 is first array

      do n=1,ncat

         do j=1,ny_block
         do i=1,nx_block
            aicen(i,j,n) = works(i,j,narrays+1)
            vicen(i,j,n) = works(i,j,narrays+2)
            vsnon(i,j,n) = works(i,j,narrays+3)
         enddo
         enddo
         narrays = narrays + 3

         call compute_tracers (nx_block,     ny_block,
     &                         icells,       indxi,   indxj,
     &                         trcr_depend,
     &                         work (:,narrays+1:narrays+ntrcr),
     &                         aicen(:,:,n),
     &                         vicen(:,:,n), vsnon(:,:,n),
     &                         trcrn(:,:,:,n))

         narrays = narrays + ntrcr

      enddo                     ! ncat

      end subroutine work_to_state

!=======================================================================
!BOP
!
! !IROUTINE: upwind_field - advection according to upwind
!
! !INTERFACE:
!
      subroutine upwind_field (nx_block, ny_block,
     &                         nghost,   dt,
     &                         narrays,  phi, 
     &                         uee,      vnn, 
     &                         HTE,      HTN, 
     &                         tarea)
!
! !DESCRIPTION:
!
! upwind transport algorithm
!
! !REVISION HISTORY:
!
! same as module
!
! !USES:
!
! !INPUT/OUTPUT PARAMETERS:
!
      integer (kind=int_kind), intent (in) ::
     &   nx_block, ny_block   ! block dimensions
     &,  nghost               ! number of ghost cells
     &,  narrays              ! number of 2D arrays to be transported

      real (kind=dbl_kind), intent(in) ::
     &   dt                   ! time step

      real (kind=dbl_kind), dimension(nx_block,ny_block,narrays),
     &   intent(inout) ::
     &   phi                  ! scalar field

      real (kind=dbl_kind), dimension(nx_block,ny_block), 
     &   intent(in)::
     &   uee, vnn             ! cell edge velocities

      real (kind=dbl_kind), dimension(nx_block,ny_block), 
     &   intent(in)::
     &   HTE,                 ! length of east cell edge 
     &   HTN,                 ! length of north cell edge 
     &   tarea                ! grid cell area
!
!EOP
!
      integer (kind=int_kind) ::
     &   i, j, k, n
     &,  ilo,ihi,jlo,jhi      ! beginning and end of physical domain

      real (kind=dbl_kind) ::
     &   upwind, y1, y2, a, h   ! function

    !-------------------------------------------------------------------
    ! Define upwind function
    !-------------------------------------------------------------------

      upwind(y1,y2,a,h)=p5*dt*h*((a+abs(a))*y1+(a-abs(a))*y2)

    !-------------------------------------------------------------------
    ! upwind transport
    !-------------------------------------------------------------------

      ilo = 1 + nghost
      ihi = nx_block - nghost
      jlo = 1 + nghost
      jhi = ny_block - nghost

      do n = 1, narrays

         do j = 1, jhi
         do i = 1, ihi
            worka(i,j)=
     &         upwind(phi(i,j,n),phi(i+1,j,n),uee(i,j),HTE(i,j))
            workb(i,j)=
     &         upwind(phi(i,j,n),phi(i,j+1,n),vnn(i,j),HTN(i,j))
         enddo
         enddo

         do j = jlo, jhi
         do i = ilo, ihi
            phi(i,j,n) = phi(i,j,n) - ( worka(i,j)-worka(i-1,j)
     &                                + workb(i,j)-workb(i,j-1) )
     &                                / tarea(i,j)
         enddo
         enddo

      enddo                     ! narrays

      end subroutine upwind_field

!=======================================================================

      end module ice_transport_upwind

!=======================================================================
