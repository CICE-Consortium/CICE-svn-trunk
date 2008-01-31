!=======================================================================
!BOP
!
! !MODULE: ice_transport_driver - drivers for ice transport
!
! !DESCRIPTION:
!
! Drivers for remapping and upwind ice transport
!
! !REVISION HISTORY:
!  SVN:$Id: ice_transport_upwind.F 28 2006-11-03 20:32:53Z eclare $
!
! authors: Elizabeth C. Hunke and William H. Lipscomb, LANL 
!
! 2004: Revised by William Lipscomb from ice_transport_mpdata.
!       Stripped out mpdata, retained upwind, and added block structure.
! 2006: Incorporated remap transport driver and renamed from
!       ice_transport_upwind.  
!
! !INTERFACE:

      module ice_transport_driver
!
! !USES:
      use ice_kinds_mod
      use ice_communicate, only: my_task, master_task
      use ice_domain_size
      use ice_constants
      use ice_fileunits
!
!EOP
!
      implicit none
      save

      character (len=char_len) ::     &
         advection   ! type of advection scheme used
                     ! 'upwind' => 1st order donor cell scheme
                     ! 'remap' => remapping scheme

!lipscomb - Remove when there is again just one remapping routine
      logical, parameter ::    &
!         newremap = .false.    ! if true, call new remapping scheme
         newremap = .true.    ! if true, call new remapping scheme
                               ! if false, call old (CICE 3.14) scheme

      logical, parameter :: & ! if true, prescribe area flux across each edge  
         l_fixed_area = .false.

! NOTE: For remapping, hice, hsno, qice, and qsno are considered tracers.
!       ntrace is not equal to ntrcr!

      integer (kind=int_kind), parameter ::                      &
         ntrace = 2+ntrcr+nilyr+nslyr  ! hice,hsno,qice,qsno,trcr
                          
      integer (kind=int_kind), dimension (ntrace) ::             &
         tracer_type       ,&! = 1, 2, or 3 (see comments below)
         depend              ! tracer dependencies (see below)

      logical (kind=log_kind), dimension (ntrace) ::             &
         has_dependents      ! true if a tracer has dependent tracers

      integer (kind=int_kind), parameter ::                      &
         integral_order = 3   ! polynomial order of quadrature integrals
                              ! linear=1, quadratic=2, cubic=3

      logical (kind=log_kind), parameter ::     &
         l_dp_midpt = .true.  ! if true, find departure points using
                              ! corrected midpoint velocity
                          
!=======================================================================

      contains

!=======================================================================

!BOP
!
! !IROUTINE: init_transport - initializations for horizontal transport
!
! !INTERFACE:
!
      subroutine init_transport
!
! !DESCRIPTION:
!
! This subroutine is a wrapper for init_remap, which initializes the
! remapping transport scheme.  If the model is run with upwind
! transport, no initializations are necessary.
!
! !REVISION HISTORY:
!
! authors William H. Lipscomb, LANL
!
! !USES:
!
      use ice_state, only: trcr_depend
      use ice_exit
      use ice_timers
      use ice_transport_remap, only: init_remap

!lipscomb - delete later
      use ice_transport_remap, only: init_remap_old
!
!EOP
!
      integer (kind=int_kind) ::       &
         k, nt, nt1     ! tracer indices

      call ice_timer_start(timer_advect)  ! advection 

      if (trim(advection)=='remap') then

!lipscomb - two branches for now; consolidate later
       if (newremap) then

         ! define tracer dependency arrays
         ! see comments in remapping routine

          depend(1:2)         = 0 ! hice, hsno
          tracer_type(1:2)    = 1 ! no dependency
      
          k = 2

          do nt = 1, ntrcr
             depend(k+nt) = trcr_depend(nt) ! 0 for ice area tracers
                                            ! 1 for ice volume tracers
                                            ! 2 for snow volume tracers
             if (trcr_depend(nt) == 0) then
                tracer_type(k+nt) = 1
             else               ! trcr_depend = 1 or 2
                tracer_type(k+nt) = 2
             endif
          enddo

          k = k + ntrcr
          
          depend(k+1:k+nilyr) = 1 ! qice depends on hice
          tracer_type(k+1:k+nilyr) = 2 

          k = k + nilyr

          depend(k+1:k+nslyr) = 2 ! qsno depends on hsno
          tracer_type(k+1:k+nslyr) = 2 

          has_dependents = .false.
          do nt = 1, ntrace
             if (depend(nt) > 0) then
                nt1 = depend(nt)
                has_dependents(nt1) = .true.
                if (nt1 > nt) then
                   write(nu_diag,*)     &
                      'Tracer nt2 =',nt,' depends on tracer nt1 =',nt1
                   call abort_ice       &
                      ('ice: remap transport: Must have nt2 > nt1')
                endif
             endif
          enddo                 ! ntrace

          call init_remap    ! grid quantities

       else      ! old remapping

         call init_remap_old (ntrace,    tracer_type,  &
                              depend,    has_dependents)

       endif  ! newremap

      else   ! upwind

         continue

      endif

      call ice_timer_stop(timer_advect)  ! advection 

      end subroutine init_transport

!=======================================================================
!BOP
!
! !IROUTINE: transport_remap - wrapper for remapping transport scheme
!
! !INTERFACE:
!
      subroutine transport_remap (dt)
!
! !DESCRIPTION:
!
! This subroutine solves the transport equations for one timestep
! using the conservative remapping scheme developed by John Dukowicz
! and John Baumgardner (DB) and modified for sea ice by William
! Lipscomb and Elizabeth Hunke.
!
! This scheme preserves monotonicity of ice area and tracers.  That is,
! it does not produce new extrema.  It is second-order accurate in space,
! except where gradients are limited to preserve monotonicity. 
!
! !REVISION HISTORY:
!
! authors William H. Lipscomb, LANL
!
! !USES:
!
      use ice_boundary
      use ice_global_reductions
      use ice_domain
      use ice_blocks
      use ice_state
      use ice_grid, only: tarea, HTE, HTN
      use ice_exit
      use ice_work, only: work1
      use ice_calendar, only: istep1
      use ice_timers
      use ice_transport_remap, only: horizontal_remap, make_masks

!lipscomb - delete later
      use ice_transport_remap, only: horizontal_remap_old

!
! !INPUT/OUTPUT PARAMETERS:
!
      real (kind=dbl_kind), intent(in) ::     &
         dt      ! time step
!
!EOP
!
      ! local variables

      integer (kind=int_kind) ::     &
         i, j           ,&! horizontal indices
         iblk           ,&! block index
         ilo,ihi,jlo,jhi,&! beginning and end of physical domain
         n              ,&! ice category index
         nt, nt1, nt2     ! tracer indices

      real (kind=dbl_kind),      &
         dimension (nx_block,ny_block,0:ncat,max_blocks) ::     &
         aim            ,&! mean ice category areas in each grid cell
         aimask           ! = 1. if ice is present, = 0. otherwise

      real (kind=dbl_kind),      &
         dimension (nx_block,ny_block,ntrace,ncat,max_blocks) ::     &
         trm            ,&! mean tracer values in each grid cell
         trmask           ! = 1. if tracer is present, = 0. otherwise

      logical (kind=log_kind) ::     &
         l_stop           ! if true, abort the model

      integer (kind=int_kind) ::     &
         istop, jstop     ! indices of grid cell where model aborts 

      integer (kind=int_kind), dimension(0:ncat,max_blocks) ::     &
         icellsnc         ! number of cells with ice

      integer (kind=int_kind),      &
         dimension(nx_block*ny_block,0:ncat,max_blocks) ::     &
         indxinc, indxjnc   ! compressed i/j indices

    !-------------------------------------------------------------------
    ! If l_fixed_area is true, the area of each departure region is
    !  computed in advance (e.g., by taking the divergence of the 
    !  velocity field and passed to locate_triangles.  The departure 
    !  regions are adjusted to obtain the desired area.
    ! If false, edgearea is computed in locate_triangles and passed out.
    !-------------------------------------------------------------------

      real (kind=dbl_kind), dimension(nx_block,ny_block,max_blocks) ::   &
         edgearea_e     ,&! area of departure regions for east edges
         edgearea_n       ! area of departure regions for north edges

      ! variables related to optional bug checks

      logical (kind=log_kind), parameter ::     &
         l_conservation_check = .false. ,&! if true, check conservation
         l_monotonicity_check = .false.   ! if true, check monotonicity

      real (kind=dbl_kind), dimension(0:ncat) ::     &
         asum_init      ,&! initial global ice area
         asum_final       ! final global ice area

      real (kind=dbl_kind), dimension(ntrace,ncat) ::     &
         atsum_init     ,&! initial global ice area*tracer
         atsum_final      ! final global ice area*tracer

      real (kind=dbl_kind), dimension (:,:,:,:,:), allocatable ::     &
         tmin         ,&! local min tracer
         tmax           ! local max tracer

      integer (kind=int_kind) :: alloc_error

      call ice_timer_start(timer_advect)  ! advection 

!---!-------------------------------------------------------------------
!---! Prepare for remapping.
!---! Initialize, update ghost cells, fill tracer arrays.
!---!-------------------------------------------------------------------

      l_stop = .false.
      istop = 0
      jstop = 0

      ilo = 1 + nghost
      ihi = nx_block - 1
      jlo = 1 + nghost
      jhi = ny_block - 1

    !-------------------------------------------------------------------
    ! Compute open water area in each grid cell.
    ! Note: An aggregate_area call is needed only if the open
    !       water area has changed since the previous call.
    !       Here we assume that aice0 is up to date.
    !-------------------------------------------------------------------

!      do iblk = 1, nblocks
!         call aggregate_area (nx_block, ny_block,
!                              iblk,     &
!                              aicen(:,:,:,iblk),     &
!                              aice (:,:,  iblk),     &
!                              aice0(:,:,  iblk)) 
!      enddo

    !-------------------------------------------------------------------
    ! Ghost cell updates for state variables.
    ! Commented out because ghost cells are updated after cleanup_itd.
    !-------------------------------------------------------------------
!      call ice_timer_start(timer_bound)

!      call ice_HaloUpdate (aice0,            halo_info,     &
!                           field_loc_center, field_type_scalar)

!      call bound_state (aicen, trcrn,     &
!                        vicen, vsnon,      &
!                        eicen, esnon)

!      call ice_timer_stop(timer_bound)

    !-------------------------------------------------------------------
    ! Ghost cell updates for ice velocity.
    ! Commented out because ghost cell velocities are computed
    !  in ice_dyn_evp.
    !-------------------------------------------------------------------

!      call ice_timer_start(timer_bound)
!      call ice_HaloUpdate (uvel,               halo_info,     &
!                           field_loc_NEcorner, field_type_vector)
!      call ice_HaloUpdate (vvel,               halo_info,     &
!                           field_loc_NEcorner, field_type_vector)
!      call ice_timer_stop(timer_bound)


      do iblk = 1, nblocks

    !-------------------------------------------------------------------
    ! Fill arrays with fields to be remapped.
    !-------------------------------------------------------------------

         call state_to_tracers(nx_block,          ny_block,             &
                               aice0(:,:,  iblk),                       &
                               aicen(:,:,:,iblk), trcrn(:,:,:,:,iblk),  &
                               vicen(:,:,:,iblk), vsnon(:,:,  :,iblk),  &
                               eicen(:,:,:,iblk), esnon(:,:,  :,iblk),  &
                               aim  (:,:,:,iblk), trm  (:,:,:,:,iblk))

      enddo

!---!-------------------------------------------------------------------
!---! Optional conservation and monotonicity checks.
!---!-------------------------------------------------------------------

      if (l_conservation_check) then

    !-------------------------------------------------------------------
    ! Compute initial values of globally conserved quantities.
    !-------------------------------------------------------------------

         do n = 0, ncat
            asum_init(n) = global_sum(aim(:,:,n,:),     distrb_info,       &
                                      field_loc_center, tarea)
         enddo

         do n = 1, ncat
            do nt = 1, ntrace
               if (tracer_type(nt)==1) then ! does not depend on another tracer
                  atsum_init(nt,n) =      &
                      global_sum_prod(trm(:,:,nt,n,:), aim(:,:,n,:),       &
                                      distrb_info,     field_loc_center,   &
                                      tarea)
               elseif (tracer_type(nt)==2) then ! depends on another tracer
                  nt1 = depend(nt)
                  work1(:,:,:) = trm(:,:,nt,n,:)*trm(:,:,nt1,n,:)
                  atsum_init(nt,n) =     &
                      global_sum_prod(work1(:,:,:), aim(:,:,n,:),          &
                                      distrb_info,  field_loc_center,      &
                                      tarea)
               elseif (tracer_type(nt)==3) then ! depends on two tracers
                  nt1 = depend(nt)
                  nt2 = depend(nt1)
                  work1(:,:,:) = trm(:,:,nt,n,:)*trm(:,:,nt1,n,:)          &
                                                *trm(:,:,nt2,n,:)
                  atsum_init(nt,n) =     &
                      global_sum_prod(work1(:,:,:), aim(:,:,n,:),          &
                                      distrb_info,  field_loc_center,      &
                                      tarea)
               endif            ! tracer_type
            enddo               ! nt
         enddo                  ! n

      endif                     ! l_conservation_check
      
      if (l_monotonicity_check) then

         allocate(tmin(nx_block,ny_block,ntrace,ncat,max_blocks),     &
                  tmax(nx_block,ny_block,ntrace,ncat,max_blocks),     &
                  STAT=alloc_error)

         if (alloc_error /= 0)      &
              call abort_ice ('ice: allocation error')

         tmin(:,:,:,:,:) = c0
         tmax(:,:,:,:,:) = c0

         do iblk = 1, nblocks

    !------------------------------------------------------------------- 
    ! Compute masks.
    ! Masks are used to prevent tracer values in cells without ice
    !  from being used in the monotonicity check.
    !------------------------------------------------------------------- 

            call make_masks (nx_block,          ny_block,              &
                             nghost,            has_dependents,        &
                             icellsnc(:,iblk),                         &
                             indxinc(:,:,iblk), indxjnc(:,:,iblk),     &
                             aim(:,:,:,iblk),   aimask(:,:,:,iblk),    &
                             trm(:,:,:,:,iblk), trmask(:,:,:,:,iblk))

    !-------------------------------------------------------------------
    ! Compute local max and min of tracer fields.
    !-------------------------------------------------------------------

            do n = 1, ncat
               call local_max_min                                      &  
                            (nx_block,           ny_block,             &
                             nghost,                                   &
                             trm (:,:,:,n,iblk),                       &
                             tmin(:,:,:,n,iblk), tmax  (:,:,:,n,iblk), &
                             aimask(:,:,n,iblk), trmask(:,:,:,n,iblk))
            enddo
         enddo

         call ice_timer_start(timer_bound)
         call ice_HaloUpdate (tmin,             halo_info,     &
                              field_loc_center, field_type_scalar)
         call ice_HaloUpdate (tmax,             halo_info,     &
                              field_loc_center, field_type_scalar)
         call ice_timer_stop(timer_bound)

         do iblk = 1, nblocks
            do n = 1, ncat
               call quasilocal_max_min (nx_block, ny_block,     &
                                        nghost,     &
                                        tmin(:,:,:,n,iblk),      &
                                        tmax(:,:,:,n,iblk))
            enddo
         enddo

      endif                     ! l_monotonicity_check

    !-------------------------------------------------------------------
    ! Main remapping routine: Step ice area and tracers forward in time.
    ! Two choices: (1) Old code, based on Lipscomb and Hunke (2004)
    !              (2) New code, incorporating Mats Bentsen's ideas
    !                  for prescribing area fluxes across each edge
    !-------------------------------------------------------------------

      if (newremap) then

    !-------------------------------------------------------------------
    ! If l_fixed_area is true, compute edgearea by taking the divergence
    !  of the velocity field.  Otherwise, initialize edgearea.
    !-------------------------------------------------------------------

         do iblk = 1, nblocks
         do j = 1, ny_block
         do i = 1, nx_block
            edgearea_e(i,j,iblk) = c0
            edgearea_n(i,j,iblk) = c0
         enddo
         enddo
         enddo

         if (l_fixed_area) then

            do iblk = 1, nblocks
               do j = jlo, jhi
               do i = ilo-1, ihi
                  edgearea_e(i,j,iblk) = (uvel(i,j,iblk) + uvel(i,j-1,iblk)) &
                                        * p5 * HTE(i,j,iblk) * dt
               enddo
               enddo

               do j = jlo-1, jhi
               do i = ilo, ihi
                  edgearea_n(i,j,iblk) = (vvel(i,j,iblk) + vvel(i-1,j,iblk)) &
                                        * p5 * HTN(i,j,iblk) * dt
               enddo
               enddo

            enddo  ! iblk

         endif

         call horizontal_remap (dt,                                    &
                                uvel      (:,:,:), vvel      (:,:,:),  &
                                aim     (:,:,:,:), trm   (:,:,:,:,:),  &
                                l_fixed_area,                          &
                                edgearea_e(:,:,:), edgearea_n(:,:,:),  &
                                tracer_type,       depend,             &
                                has_dependents,    integral_order,     &
                                l_dp_midpt)

      else  ! standard remapping

         call horizontal_remap_old (dt,                                &
                                uvel   (:,:,:),    vvel      (:,:,:),  &
                                aim  (:,:,:,:),    trm   (:,:,:,:,:),  &
                                edgearea_e(:,:,:), edgearea_n(:,:,:),  &
                                tracer_type,       depend,             &
                                has_dependents,    integral_order,     &
                                l_dp_midpt)

      endif  ! newremap

    !-------------------------------------------------------------------
    ! Given new fields, recompute state variables.
    !-------------------------------------------------------------------

      do iblk = 1, nblocks

         call tracers_to_state (nx_block,          ny_block,            &
                                aim  (:,:,:,iblk), trm  (:,:,:,:,iblk), &
                                aice0(:,:,  iblk),                      &
                                aicen(:,:,:,iblk), trcrn(:,:,:,:,iblk), &
                                vicen(:,:,:,iblk), vsnon(:,:,  :,iblk), &
                                eicen(:,:,:,iblk), esnon(:,:,  :,iblk)) 

      enddo                     ! iblk

    !-------------------------------------------------------------------
    ! Ghost cell updates for state variables.
    !-------------------------------------------------------------------

      call ice_timer_start(timer_bound)

      call bound_state (aicen, trcrn,     &
                        vicen, vsnon,      &
                        eicen, esnon)

      call ice_timer_stop(timer_bound)

!---!-------------------------------------------------------------------
!---! Optional conservation and monotonicity checks
!---!-------------------------------------------------------------------

    !-------------------------------------------------------------------
    ! Compute final values of globally conserved quantities.
    ! Check global conservation of area and area*tracers.  (Optional)
    !-------------------------------------------------------------------

      if (l_conservation_check) then

         do n = 0, ncat
            asum_final(n) = global_sum(aim(:,:,n,:),     distrb_info,      &
                                       field_loc_center, tarea)
         enddo

         do n = 1, ncat
            do nt = 1, ntrace
               if (tracer_type(nt)==1) then ! does not depend on another tracer
                  atsum_final(nt,n) =      &
                      global_sum_prod(trm(:,:,nt,n,:), aim(:,:,n,:),       &
                                      distrb_info,     field_loc_center,   &
                                      tarea)
               elseif (tracer_type(nt)==2) then ! depends on another tracer
                  nt1 = depend(nt)
                  work1(:,:,:) = trm(:,:,nt,n,:)*trm(:,:,nt1,n,:)
                  atsum_final(nt,n) =     &
                      global_sum_prod(work1(:,:,:), aim(:,:,n,:),          &
                                      distrb_info,  field_loc_center,      &
                                      tarea)
               elseif (tracer_type(nt)==3) then ! depends on two tracers
                  nt1 = depend(nt)
                  nt2 = depend(nt1)
                  work1(:,:,:) = trm(:,:,nt,n,:)*trm(:,:,nt1,n,:)          &
                                                *trm(:,:,nt2,n,:)
                  atsum_final(nt,n) =     &
                      global_sum_prod(work1(:,:,:), aim(:,:,n,:),          &
                                      distrb_info,  field_loc_center,      &
                                      tarea)
               endif            ! tracer_type
            enddo               ! nt
         enddo                  ! n


         if (my_task == master_task) then
            call global_conservation (l_stop,     &
                                      asum_init(0), asum_final(0))

            if (l_stop) then
               write (nu_diag,*) 'istep1, my_task, iblk =',     &
                                  istep1, my_task, iblk
               write (nu_diag,*) 'transport: conservation error, cat 0'
               call abort_ice('ice remap transport: conservation error')
            endif

            do n = 1, ncat               
               call global_conservation                                 &
                                     (l_stop,                           &
                                      asum_init(n),    asum_final(n),   &
                                      atsum_init(:,n), atsum_final(:,n))

               if (l_stop) then
                  write (nu_diag,*) 'istep1, my_task, iblk, cat =',     &
                                     istep1, my_task, iblk, n
                  write (nu_diag,*) 'transport: conservation error, cat ',n
                  call abort_ice     &
                       ('ice remap transport: conservation error')
               endif
            enddo               ! n

         endif                  ! my_task = master_task

      endif                     ! l_conservation_check

    !-------------------------------------------------------------------
    ! Check tracer monotonicity.  (Optional)
    !-------------------------------------------------------------------

      if (l_monotonicity_check) then
         do iblk = 1, nblocks
            do n = 1, ncat
               call check_monotonicity      &
                               (nx_block,           ny_block,     &
                                nghost,     &
                                iblk,     &
                                tmin(:,:,:,n,iblk), tmax(:,:,:,n,iblk),  &
                                aim (:,:,  n,iblk), trm (:,:,:,n,iblk),  &
                                l_stop,     &
                                istop,              jstop)
            enddo               ! n

            if (l_stop) then
               write (nu_diag,*) 'istep1, my_task, iblk, cat =',     &
                                  istep1, my_task, iblk, n
               call abort_ice('ice remap transport: monotonicity error')
            endif

         enddo                  ! iblk

         deallocate(tmin, tmax, STAT=alloc_error)
         if (alloc_error /= 0) call abort_ice ('deallocation error')

      endif                     ! l_monotonicity_check

      call ice_timer_stop(timer_advect)  ! advection 

      end subroutine transport_remap

!=======================================================================
!BOP
!
! !IROUTINE: transport_upwind - upwind transport
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
      real (kind=dbl_kind), intent(in) ::     &
         dt      ! time step
!
!EOP
!
      integer (kind=int_kind), parameter ::     &
         narr = 1 + ncat*(3+ntrcr)   ! number of state variable arrays
                                     ! not including eicen, esnon

      integer (kind=int_kind) ::     &
         i, j, iblk       ,&! horizontal indices
         ilo,ihi,jlo,jhi    ! beginning and end of physical domain

      real (kind=dbl_kind), dimension (nx_block,ny_block,nblocks) ::     &
         uee, vnn           ! cell edge velocities

      real (kind=dbl_kind),     &
         dimension (nx_block,ny_block,narr,max_blocks) ::      &
         works              ! work array

      type (block) ::     &
         this_block           ! block information for current block

      call ice_timer_start(timer_advect)  ! advection 

    !-------------------------------------------------------------------
    ! Get ghost cell values of state variables.
    ! (Assume velocities are already known for ghost cells, also.)
    !-------------------------------------------------------------------
!      call bound_state (aicen, trcrn,     &
!                        vicen, vsnon,     &
!                        eicen, esnon)

    !-------------------------------------------------------------------
    ! Average corner velocities to edges.
    !-------------------------------------------------------------------
      
      do iblk = 1, nblocks
         this_block = get_block(blocks_ice(iblk),iblk)         
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
      call ice_HaloUpdate (uee,             halo_info,     &
                           field_loc_Eface, field_type_scalar)
      call ice_HaloUpdate (vnn,             halo_info,     &
                           field_loc_Nface, field_type_scalar)
      call ice_timer_stop(timer_bound)

      do iblk = 1, nblocks

      !-----------------------------------------------------------------
      ! fill work arrays with fields to be advected
      !-----------------------------------------------------------------

         call state_to_work (nx_block,             ny_block,             &
                             narr,                 trcr_depend,          &
                             aicen (:,:,  :,iblk), trcrn (:,:,:,:,iblk), &
                             vicen (:,:,  :,iblk), vsnon (:,:,  :,iblk), &
                             aice0 (:,:,    iblk), works (:,:,  :,iblk))

      !-----------------------------------------------------------------
      ! advect
      !-----------------------------------------------------------------

         call upwind_field (nx_block,       ny_block,               &
                            nghost,         dt,                     &
                            narr,           works(:,:,:,iblk),      &
                            uee(:,:,iblk),  vnn    (:,:,iblk),      &
                            HTE(:,:,iblk),  HTN    (:,:,iblk),      &
                            tarea(:,:,iblk))

         call upwind_field (nx_block,       ny_block,               &
                            nghost,         dt,                     &
                            ntilyr,         eicen(:,:,:,iblk),      &
                            uee(:,:,iblk),  vnn    (:,:,iblk),      &
                            HTE(:,:,iblk),  HTN    (:,:,iblk),      &
                            tarea(:,:,iblk))

         call upwind_field (nx_block,       ny_block,               &
                            nghost,         dt,                     &
                            ntslyr,         esnon(:,:,:,iblk),      &
                            uee(:,:,iblk),  vnn    (:,:,iblk),      &
                            HTE(:,:,iblk),  HTN    (:,:,iblk),      &
                            tarea(:,:,iblk))

      !-----------------------------------------------------------------
      ! convert work arrays back to state variables
      !-----------------------------------------------------------------

         call work_to_state (nx_block,            ny_block,     &
                             narr,                trcr_depend,     &
                             aicen(:,:,  :,iblk), trcrn (:,:,:,:,iblk),     &
                             vicen(:,:,  :,iblk), vsnon (:,:,  :,iblk),     &
                             aice0(:,:,    iblk), works (:,:,  :,iblk)) 

      enddo                     ! iblk
 
    !-------------------------------------------------------------------
    ! Ghost cell updates for state variables.
    !-------------------------------------------------------------------

      call ice_timer_start(timer_bound)

      call bound_state (aicen, trcrn,     &
                        vicen, vsnon,      &
                        eicen, esnon)

      call ice_timer_stop(timer_bound)

      call ice_timer_stop(timer_advect)  ! advection 

      end subroutine transport_upwind

!=======================================================================
! The next few subroutines (through check_monotonicity) are called
! by transport_remap.
!=======================================================================
!
!BOP
!
! !IROUTINE: state_to_tracers -fill ice area and tracer arrays
!
! !INTERFACE:
!
      subroutine state_to_tracers (nx_block, ny_block,   &
                                   aice0,                &
                                   aicen,    trcrn,      &
                                   vicen,    vsnon,      &
                                   eicen,    esnon,      &
                                   aim,      trm)
!
! !DESCRIPTION:
!
! Fill ice area and tracer arrays.
! Assume that the advected tracers are hicen, hsnon, trcrn, 
!  qicen(1:nilyr), and qsnon(1:nslyr).
! This subroutine must be modified if a different set of tracers
!   is to be transported.  The rule for ordering tracers
!   is that a dependent tracer (such as qice) must have a larger
!   tracer index than the tracer it depends on (i.e., hice).
!
! !REVISION HISTORY:
!
! author William H. Lipscomb, LANL
!
! !USES:
!
      use ice_work, only: worka, workb
      use ice_itd, only: ilyr1, slyr1
!
! !INPUT/OUTPUT PARAMETERS:
!
      integer (kind=int_kind), intent(in) ::     &
           nx_block, ny_block  ! block dimensions

      real (kind=dbl_kind), dimension (nx_block,ny_block),     &
           intent(in) ::     &
           aice0     ! fractional open water area

      real (kind=dbl_kind), dimension (nx_block,ny_block,ncat),     &
           intent(in) ::     &
           aicen   ,&! fractional ice area
           vicen   ,&! volume per unit area of ice          (m)
           vsnon     ! volume per unit area of snow         (m)

      real (kind=dbl_kind), dimension (nx_block,ny_block,ntrcr,ncat),     &
           intent(in) ::     &
           trcrn     ! ice area tracers

      real (kind=dbl_kind), dimension (nx_block,ny_block,ntilyr),     &
           intent(in) ::     &
           eicen     ! energy of melting for each ice layer (J/m^2)

      real (kind=dbl_kind), dimension (nx_block,ny_block,ntslyr),     &
           intent(in) ::     &
           esnon     ! energy of melting for each snow layer (J/m^2)

      real (kind=dbl_kind), dimension (nx_block,ny_block,0:ncat),     &
            intent(out)::     &
           aim       ! mean ice area in each grid cell

      real (kind=dbl_kind), dimension (nx_block,ny_block,ntrace,ncat),  &
           intent(out) ::     &
           trm       ! mean tracer values in each grid cell
!
!EOP
!
      integer (kind=int_kind) ::     &
           i, j, k, n   ,&! standard indices
           it, kt       ,&! tracer indices
           ij             ! combined i/j index

      real (kind=dbl_kind) ::     &
           w1             ! work variable

      integer (kind=int_kind), dimension(nx_block*ny_block,0:ncat) ::  &
           indxi        ,&! compressed i/j indices
           indxj

      integer (kind=int_kind), dimension(0:ncat) ::     &
           icells         ! number of cells with ice

      aim(:,:,0) = aice0(:,:)

      do n = 1, ncat

         trm(:,:,:,n) = c0

    !-------------------------------------------------------------------
    ! Find grid cells where ice is present and fill area array.
    !-------------------------------------------------------------------

         icells(n) = 0
         do j = 1, ny_block
         do i = 1, nx_block
            aim(i,j,n) = aicen(i,j,n)
            if (aim(i,j,n) > puny) then
               icells(n) = icells(n) + 1
               ij = icells(n)
               indxi(ij,n) = i
               indxj(ij,n) = j
            endif               ! aim > puny
         enddo
         enddo
      
    !-------------------------------------------------------------------
    ! Fill tracer array
    ! Note: If aice > 0, then hice > 0, but we can have hsno = 0.
    ! Alse note: We transport qice*nilyr rather than qice, so as to
    !  avoid extra operations here and in tracers_to_state.
    !-------------------------------------------------------------------

         do ij = 1, icells(n)
            i = indxi(ij,n)
            j = indxj(ij,n)
            w1 = c1 / aim(i,j,n)
            worka(i,j) = c1 / vicen(i,j,n)
            trm(i,j,1,n) = vicen(i,j,n) * w1 ! hice
            trm(i,j,2,n) = vsnon(i,j,n) * w1 ! hsno
            if (trm(i,j,2,n) > puny) then
               workb(i,j) = c1 / vsnon(i,j,n)
            else
               workb(i,j) = c0
            endif
         enddo
         kt = 2

         do it = 1, ntrcr
            do ij = 1, icells(n)
               i = indxi(ij,n)
               j = indxj(ij,n)
               trm(i,j,kt+it,n) = trcrn(i,j,it,n) ! ice area tracers
            enddo
         enddo
         kt = kt + ntrcr

         do k =1, nilyr
            do ij = 1, icells(n)
               i = indxi(ij,n)
               j = indxj(ij,n)
               trm(i,j,kt+k,n) = eicen(i,j,ilyr1(n)+k-1)*worka(i,j) ! qice
            enddo               ! ij
         enddo                  ! ilyr
         kt = kt + nilyr

         do k = 1, nslyr
            do ij = 1, icells(n)
               i = indxi(ij,n)
               j = indxj(ij,n)
               if (trm(i,j,2,n) > puny)    &    ! hsno > puny
                 trm(i,j,kt+k,n) = esnon(i,j,slyr1(n)+k-1)*workb(i,j) & ! qsno
                                 + rhos*Lfresh
            enddo               ! ij
         enddo                  ! nslyr

      enddo                     ! ncat
 
      end subroutine state_to_tracers

!=======================================================================
!BOP
!
! !IROUTINE: tracers_to_state - convert tracer array to state variables
!
! !INTERFACE:
!
      subroutine tracers_to_state (nx_block, ny_block,   &
                                   aim,      trm,        &
                                   aice0,                &
                                   aicen,    trcrn,      &
                                   vicen,    vsnon,      &
                                   eicen,    esnon) 
!
! !DESCRIPTION:
!
! Convert area and tracer arrays back to state variables.
!
! !REVISION HISTORY:
!
! author William H. Lipscomb, LANL
!
! !USES:
!
      use ice_itd, only: ilyr1, slyr1
!
! !INPUT/OUTPUT PARAMETERS:
!
      integer (kind=int_kind), intent(in) ::     &
           nx_block, ny_block  ! block dimensions

      real (kind=dbl_kind), dimension (nx_block,ny_block,0:ncat),     &
           intent(in) ::     &
           aim       ! fractional ice area

      real (kind=dbl_kind), dimension (nx_block,ny_block,ntrace,ncat),  &
           intent(in) ::     &
           trm       ! mean tracer values in each grid cell

      real (kind=dbl_kind), dimension (nx_block,ny_block),     &
           intent(inout) ::     &
           aice0     ! fractional ice area

      real (kind=dbl_kind), dimension (nx_block,ny_block,ncat),     &
           intent(inout) ::     &
           aicen   ,&! fractional ice area
           vicen   ,&! volume per unit area of ice          (m)
           vsnon     ! volume per unit area of snow         (m)

      real (kind=dbl_kind), dimension (nx_block,ny_block,ntrcr,ncat),  &
           intent(inout) ::     &
           trcrn     ! tracers

      real (kind=dbl_kind), dimension (nx_block,ny_block,ntilyr),     &
           intent(inout) ::     &
           eicen ! energy of melting for each ice layer (J/m^2)

      real (kind=dbl_kind), dimension (nx_block,ny_block,ntslyr),     &
           intent(inout) ::     &
           esnon ! energy of melting for each snow layer (J/m^2)
!
!EOP
!
      integer (kind=int_kind) ::     &
           i, j, k, n      ,&! standard indices
           it, kt          ,&! tracer indices
           icells          ,&! number of cells with ice
           ij

      integer (kind=int_kind), dimension (nx_block*ny_block) ::     &
           indxi, indxj      ! compressed indices

      aice0(:,:) = aim(:,:,0)

      do n = 1, ncat

      icells = 0
      do j = 1, ny_block
      do i = 1, nx_block
         if (aim(i,j,n) > c0) then
            icells = icells + 1
            indxi(icells) = i
            indxj(icells) = j
         endif
      enddo
      enddo

    !-------------------------------------------------------------------
    ! Compute state variables.
    !-------------------------------------------------------------------

         do ij = 1, icells
            i = indxi(ij)
            j = indxj(ij)
            aicen(i,j,n) = aim(i,j,n)
            vicen(i,j,n) = aim(i,j,n)*trm(i,j,1,n) ! aice*hice
            vsnon(i,j,n) = aim(i,j,n)*trm(i,j,2,n) ! aice*hsno
         enddo                  ! ij
         kt = 2

         do it = 1, ntrcr
         do ij = 1, icells
            i = indxi(ij)
            j = indxj(ij)
               trcrn(i,j,it,n) = trm(i,j,kt+it,n)  ! ice tracers
            enddo               ! ij
         enddo                  ! ntrcr
         kt = kt + ntrcr

         do k = 1, nilyr
         do ij = 1, icells
            i = indxi(ij)
            j = indxj(ij)
               eicen(i,j,ilyr1(n)+k-1) = vicen(i,j,n)*trm(i,j,kt+k,n) 
            enddo               ! ij
         enddo                  ! nilyr
         kt = kt + nilyr

         do k = 1, nslyr
         do ij = 1, icells
            i = indxi(ij)
            j = indxj(ij)
               esnon(i,j,slyr1(n)+k-1) = (trm(i,j,kt+k,n) - rhos*Lfresh) &
                                         * vsnon(i,j,n)
            enddo               ! ij
         enddo                  ! nslyr

      enddo                     ! ncat

      end subroutine tracers_to_state

!=======================================================================
!
!BOP
!
! !IROUTINE: global_conservation - check for changes in conserved quantities
!
! !INTERFACE:
!
      subroutine global_conservation (l_stop,                     &
                                      asum_init,  asum_final,     &
                                      atsum_init, atsum_final)
!
! !DESCRIPTION:
!
! Check whether values of conserved quantities have changed.
! An error probably means that ghost cells are treated incorrectly.
!
! !REVISION HISTORY:
!
! author William H. Lipscomb, LANL
!
! !USES:
!
! !INPUT/OUTPUT PARAMETERS:
!
      real (kind=dbl_kind), intent(in) ::     &
         asum_init   ,&! initial global ice area
         asum_final    ! final global ice area

      real (kind=dbl_kind), dimension(ntrace), intent(in), optional :: &
         atsum_init  ,&! initial global ice area*tracer
         atsum_final   ! final global ice area*tracer

      logical (kind=log_kind), intent(inout) ::     &
         l_stop    ! if true, abort on return
!
!EOP
!
      integer (kind=int_kind) ::     &
           nt            ! tracer index

      real (kind=dbl_kind) ::     &
           diff          ! difference between initial and final values


      if (asum_init > puny) then
         diff = asum_final - asum_init
         if (abs(diff/asum_init) > puny) then
            l_stop = .true.
            write (nu_diag,*)
            write (nu_diag,*) 'Ice area conserv error'
            write (nu_diag,*) 'Initial global area =', asum_init
            write (nu_diag,*) 'Final global area =', asum_final
            write (nu_diag,*) 'Fractional error =', abs(diff)/asum_init
            write (nu_diag,*) 'asum_final-asum_init =', diff
         endif
      endif

      if (present(atsum_init)) then
       do nt = 1, ntrace
         if (abs(atsum_init(nt)) > puny) then
            diff = atsum_final(nt) - atsum_init(nt)
            if (abs(diff/atsum_init(nt)) > puny) then
               l_stop = .true.
               write (nu_diag,*)
               write (nu_diag,*) 'area*tracer conserv error'
               write (nu_diag,*) 'tracer index =', nt
               write (nu_diag,*) 'Initial global area*tracer =',   &
                                  atsum_init(nt)
               write (nu_diag,*) 'Final global area*tracer =',     &
                                  atsum_final(nt)
               write (nu_diag,*) 'Fractional error =',             &
                                  abs(diff)/atsum_init(nt)
               write (nu_diag,*) 'atsum_final-atsum_init =', diff
            endif
         endif
       enddo
      endif                     ! present(atsum_init)

      end subroutine global_conservation

!=======================================================================
!BOP
!
! !IROUTINE: local_max_min - compute local max and min of a scalar field
!
! !INTERFACE:
!
      subroutine local_max_min (nx_block, ny_block,     &
                                nghost,                 &
                                trm,                    &
                                tmin,     tmax,         &
                                aimask,   trmask)
!
! !DESCRIPTION:
!
! At each grid point, compute the local max and min of a scalar
! field phi: i.e., the max and min values in the nine-cell region
! consisting of the home cell and its eight neighbors.
! 
! To extend to the neighbors of the neighbors (25 cells in all),
! follow this call with a call to quasilocal_max_min.
!
! !REVISION HISTORY:
!
! author William H. Lipscomb, LANL
!
! !USES:
!
! !INPUT/OUTPUT PARAMETERS:
!
      integer (kind=int_kind), intent(in) ::     &
           nx_block, ny_block,&! block dimensions
           nghost              ! number of ghost cells

      real (kind=dbl_kind), intent(in),        &
           dimension(nx_block,ny_block) ::     &
           aimask         ! ice area mask

      real (kind=dbl_kind), intent(in),               &
           dimension (nx_block,ny_block,ntrace) ::    &
           trm          ,&! tracer fields
           trmask         ! tracer mask

      real (kind=dbl_kind), intent(out),              &
           dimension (nx_block,ny_block,ntrace) ::    &
           tmin         ,&! local min tracer
           tmax           ! local max tracer
!
!EOP
!
      integer (kind=int_kind) ::     &
           i, j         ,&! horizontal indices
           nt, nt1      ,&! tracer indices
           ilo,ihi,jlo,jhi     ! beginning and end of physical domain

      real (kind=dbl_kind), dimension(nx_block,ny_block) ::     &
           phimask        ! aimask or trmask, as appropriate

      real (kind=dbl_kind) ::     &
           phi_nw, phi_n, phi_ne ,&! field values in 8 neighbor cells
           phi_w, phi_e          ,&
           phi_sw, phi_s, phi_se

      ilo = 1 + nghost
      ihi = nx_block - nghost
      jlo = 1 + nghost
      jhi = ny_block - nghost

      do nt = 1, ntrace

         if (tracer_type(nt)==1) then  ! does not depend on another tracer

            do j = 1, ny_block
            do i = 1, nx_block
               phimask(i,j) = aimask(i,j)
            enddo
            enddo

         else   ! depends on another tracer

            nt1 = depend(nt)
            do j = 1, ny_block
            do i = 1, nx_block
               phimask(i,j) = trmask(i,j,nt1)
            enddo
            enddo

         endif

!-----------------------------------------------------------------------
!  Store values of trm in the 8 neighbor cells.
!  If aimask = 1, use the true value; otherwise use the home cell value
!  so that non-physical values of phi do not contribute to the gradient.
!-----------------------------------------------------------------------

         do j = jlo, jhi
            do i = ilo, ihi

               phi_nw = phimask(i-1,j+1) * trm(i-1,j+1,nt)     &
                  + (c1-phimask(i-1,j+1))* trm(i,  j,  nt)
               phi_n  = phimask(i,  j+1) * trm(i,  j+1,nt)     &
                  + (c1-phimask(i,  j+1))* trm(i,  j,  nt)
               phi_ne = phimask(i+1,j+1) * trm(i+1,j+1,nt)     &
                  + (c1-phimask(i+1,j+1))* trm(i,  j,  nt)
               phi_w  = phimask(i-1,j)   * trm(i-1,j,  nt)     &
                  + (c1-phimask(i-1,j))  * trm(i,  j,  nt)
               phi_e  = phimask(i+1,j)   * trm(i+1,j,  nt)     &
                  + (c1-phimask(i+1,j))  * trm(i,  j,  nt)
               phi_sw = phimask(i-1,j-1) * trm(i-1,j-1,nt)     &
                  + (c1-phimask(i-1,j-1))* trm(i,  j,  nt)
               phi_s  = phimask(i,  j-1) * trm(i,  j-1,nt)     &
                  + (c1-phimask(i,  j-1))* trm(i,  j,  nt)
               phi_se = phimask(i+1,j-1) * trm(i+1,j-1,nt)     &
                  + (c1-phimask(i+1,j-1))* trm(i,  j,  nt)

!-----------------------------------------------------------------------
!     Compute the minimum and maximum among the nine local cells.
!-----------------------------------------------------------------------

               tmax(i,j,nt) = max (phi_nw, phi_n,  phi_ne, phi_w,     &
                      trm(i,j,nt), phi_e,  phi_sw, phi_s,  phi_se)

               tmin(i,j,nt) = min (phi_nw, phi_n,  phi_ne, phi_w,     &
                      trm(i,j,nt), phi_e,  phi_sw, phi_s,  phi_se)

            enddo               ! i
         enddo                  ! j

      enddo                     ! nt

      end subroutine local_max_min

!=======================================================================
!BOP
!
! !IROUTINE: quasilocal_max_min - look one grid cell farther away
!
! !INTERFACE:
!
      subroutine quasilocal_max_min (nx_block, ny_block,     &
                                     nghost,     &
                                     tmin,     tmax)
!
! !DESCRIPTION:
!
! Extend the local max and min by one grid cell in each direction.
! Incremental remapping is monotone for the "quasilocal" max and min,
! but in rare cases may violate monotonicity for the local max and min.
!
! !REVISION HISTORY:
!
! author William H. Lipscomb, LANL
!
! !USES:
!
! !INPUT/OUTPUT PARAMETERS:
!
      integer (kind=int_kind), intent(in) ::     &
         nx_block, ny_block,&! block dimensions
         nghost              ! number of ghost cells

      real (kind=dbl_kind), intent(inout),     &
           dimension (nx_block,ny_block,ntrace) ::     &
           tmin         ,&! local min tracer
           tmax           ! local max tracer
!
!EOP
!
      integer (kind=int_kind) ::     &
           i, j          ,&! horizontal indices
           nt            ,&! tracer index
           ilo,ihi,jlo,jhi ! beginning and end of physical domain

      ilo = 1 + nghost
      ihi = nx_block - nghost
      jlo = 1 + nghost
      jhi = ny_block - nghost

      do nt = 1, ntrace

         do j = jlo, jhi
         do i = ilo, ihi

            tmax(i,j,nt) =     &
              max (tmax(i-1,j+1,nt), tmax(i,j+1,nt), tmax(i+1,j+1,nt),     &
                   tmax(i-1,j,  nt), tmax(i,j,  nt), tmax(i+1,j,  nt),     &
                   tmax(i-1,j-1,nt), tmax(i,j-1,nt), tmax(1+1,j-1,nt))

            tmin(i,j,nt) =     &
              min (tmin(i-1,j+1,nt), tmin(i,j+1,nt), tmin(i+1,j+1,nt),     &
                   tmin(i-1,j,  nt), tmin(i,j,  nt), tmin(i+1,j,  nt),     &
                   tmin(i-1,j-1,nt), tmin(i,j-1,nt), tmin(1+1,j-1,nt))

         enddo                  ! i
         enddo                  ! j

      enddo

      end subroutine quasilocal_max_min

!======================================================================
!
!BOP
!
! !IROUTINE: check_monotonicity - check bounds on new tracer values
!
! !INTERFACE:
!
      subroutine check_monotonicity (nx_block, ny_block,     &
                                     nghost,                 &
                                     iblk,                   &
                                     tmin,     tmax,         &
                                     aim,      trm,          &
                                     l_stop,                 &
                                     istop,    jstop)
!
! !DESCRIPTION:
!
! At each grid point, make sure that the new tracer values
! fall between the local max and min values before transport.
!
! !REVISION HISTORY:
!
! author William H. Lipscomb, LANL
!
! !USES:
!
! !INPUT/OUTPUT PARAMETERS:
!
      integer (kind=int_kind), intent(in) ::     &
           nx_block, ny_block,&! block dimensions
           nghost            ,&! number of ghost cells
           iblk                ! block index (diagnostic only)

      real (kind=dbl_kind), intent(in),         &
           dimension (nx_block,ny_block) ::     &
           aim            ! new ice area

      real (kind=dbl_kind), intent(in),                &
           dimension (nx_block,ny_block,ntrace) ::     &
           trm            ! new tracers

      real (kind=dbl_kind), intent(in),                &
           dimension (nx_block,ny_block,ntrace) ::     &
           tmin         ,&! local min tracer
           tmax           ! local max tracer

      logical (kind=log_kind), intent(inout) ::     &
         l_stop    ! if true, abort on return

      integer (kind=int_kind), intent(inout) ::     &
         istop, jstop     ! indices of grid cell where model aborts 
!
!EOP
!
      integer (kind=int_kind) ::     &
           i, j           ,&! horizontal indices
           nt, nt1, nt2   ,&! tracer indices
           ilo,ihi,jlo,jhi  ! beginning and end of physical domain

      real (kind=dbl_kind) ::     &
           w1, w2         ! work variables

      logical (kind=log_kind), dimension (nx_block, ny_block) ::   &
           l_check        ! if true, check monotonicity

    !-------------------------------------------------------------------
    ! Initialize
    !-------------------------------------------------------------------

      ilo = 1 + nghost
      ihi = nx_block - nghost
      jlo = 1 + nghost
      jhi = ny_block - nghost

      do nt = 1, ntrace

    !-------------------------------------------------------------------
    ! Load logical array to identify tracers that need checking.
    !-------------------------------------------------------------------

         if (tracer_type(nt)==1) then ! does not depend on another tracer

            do j = jlo, jhi
            do i = ilo, ihi
               if (aim(i,j) > puny) then 
                  l_check(i,j) = .true.
               else
                  l_check(i,j) = .false.
               endif
            enddo
            enddo

         elseif (tracer_type(nt)==2) then ! depends on another tracer

            nt1 = depend(nt)
            do j = jlo, jhi
            do i = ilo, ihi
               if (abs(trm(i,j,nt1)) > puny) then
                  l_check(i,j) = .true.
               else
                  l_check(i,j) = .false.
               endif
            enddo
            enddo

         elseif (tracer_type(nt)==2) then ! depends on another tracer

            nt1 = depend(nt)
            nt2 = depend(nt1)
            do j = jlo, jhi
            do i = ilo, ihi
               if (abs(trm(i,j,nt1)) > puny .and.     &
                   abs(trm(i,j,nt2)) > puny) then
                  l_check(i,j) = .true.
               else
                  l_check(i,j) = .false.
               endif
            enddo
            enddo
         endif

    !-------------------------------------------------------------------
    ! Make sure new values lie between tmin and tmax
    !-------------------------------------------------------------------

         do j = jlo, jhi
         do i = ilo, ihi

            if (l_check(i,j)) then
               ! w1 and w2 allow for roundoff error when abs(trm) is big
               w1 = max(c1, abs(tmin(i,j,nt)))
               w2 = max(c1, abs(tmax(i,j,nt)))
               if (trm(i,j,nt) < tmin(i,j,nt)-w1*puny) then
                  l_stop = .true.
                  istop = i
                  jstop = j
                  write (nu_diag,*) ' '
                  write (nu_diag,*) 'new tracer < tmin'
                  write (nu_diag,*) 'i, j, nt =', i, j, nt
                  write (nu_diag,*) 'new tracer =', trm (i,j,nt)
                  write (nu_diag,*) 'tmin ='      , tmin(i,j,nt)
                  write (nu_diag,*) 'ice area ='  , aim(i,j)
               elseif (trm(i,j,nt) > tmax(i,j,nt)+w2*puny) then
                  l_stop = .true.
                  istop = i
                  jstop = j
                  write (nu_diag,*) ' '
                  write (nu_diag,*) 'new tracer > tmax'
                  write (nu_diag,*) 'i, j, nt =', i, j, nt
                  write (nu_diag,*) 'new tracer =', trm (i,j,nt)
                  write (nu_diag,*) 'tmax ='      , tmax(i,j,nt)
                  write (nu_diag,*) 'ice area ='  , aim(i,j)
               endif
            endif

         enddo                  ! i
         enddo                  ! j

      enddo                     ! nt

      end subroutine check_monotonicity

!=======================================================================
! The remaining subroutines are called by transport_upwind.
!=======================================================================
!BOP
!
! !IROUTINE: state_to_work - fill work arrays with state variables
!
! !INTERFACE:
!
      subroutine state_to_work (nx_block, ny_block,        &
                                narr,     trcr_depend,     &
                                aicen,    trcrn,           &
                                vicen,    vsnon,           &
                                aice0,    works)

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
      integer (kind=int_kind), intent(in) ::     &
         nx_block, ny_block ,&! block dimensions
         narr        ! number of 2D state variable arrays in works array

      integer (kind=int_kind), dimension (ntrcr), intent(in) ::     &
         trcr_depend ! = 0 for aicen tracers, 1 for vicen, 2 for vsnon

      real (kind=dbl_kind), dimension (nx_block,ny_block,ncat),     &
         intent(in) ::     &
         aicen   ,&! concentration of ice
         vicen   ,&! volume per unit area of ice          (m)
         vsnon     ! volume per unit area of snow         (m)

      real (kind=dbl_kind), dimension (nx_block,ny_block,ntrcr,ncat),     &
         intent(in) ::     &
         trcrn     ! ice tracers

      real (kind=dbl_kind), dimension (nx_block,ny_block),         &
         intent(in) ::        &
         aice0     ! concentration of open water

      real (kind=dbl_kind), dimension(nx_block,ny_block,narr),     &
         intent (out) ::      &
         works     ! work array
!
!EOP
!
      integer (kind=int_kind) ::      &
         i, j, k, n, it ,&! counting indices
         narrays          ! counter for number of state variable arrays

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

      if (narr /= narrays) write(nu_diag,*)      &
           "Wrong number of arrays in transport bound call"

      end subroutine state_to_work

!=======================================================================
!BOP
!
! !IROUTINE: work_to_state - convert work arrays back to state variables
!
! !INTERFACE:
!
      subroutine work_to_state (nx_block, ny_block,        &
                                narr,     trcr_depend,     &
                                aicen,    trcrn,           &
                                vicen,    vsnon,           &
                                aice0,    works)

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
      integer (kind=int_kind), intent (in) ::                       &
         nx_block, ny_block ,&! block dimensions
         narr        ! number of 2D state variable arrays in works array

      integer (kind=int_kind), dimension (ntrcr), intent(in) ::     &
         trcr_depend ! = 0 for aicen tracers, 1 for vicen, 2 for vsnon

      real (kind=dbl_kind), intent (in) ::                          &
         works (nx_block,ny_block,narr)

      real (kind=dbl_kind), dimension (nx_block,ny_block,ncat),     &
         intent(out) ::     &
         aicen   ,&! concentration of ice
         vicen   ,&! volume per unit area of ice          (m)
         vsnon     ! volume per unit area of snow         (m)

      real (kind=dbl_kind), dimension (nx_block,ny_block,ntrcr,ncat), &
         intent(out) ::     &
         trcrn     ! ice tracers

      real (kind=dbl_kind), dimension (nx_block,ny_block),          &
         intent(out) ::     &
         aice0     ! concentration of open water
!
!EOP
!
      integer (kind=int_kind) ::      &
         i, j, k, n , it,&! counting indices
         narrays        ,&! counter for number of state variable arrays
         icells           ! number of ocean/ice cells

      integer (kind=int_kind), dimension (nx_block*ny_block) ::        &
        indxi, indxj

      real (kind=dbl_kind), dimension (nx_block*ny_block,narr) ::      &
         work 

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

         call compute_tracers (nx_block,     ny_block,               &
                               icells,       indxi,   indxj,         &
                               trcr_depend,                          &
                               work (:,narrays+1:narrays+ntrcr),     &
                               aicen(:,:,n),                         &
                               vicen(:,:,n), vsnon(:,:,n),           &
                               trcrn(:,:,:,n))

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
      subroutine upwind_field (nx_block, ny_block,   &
                               nghost,   dt,         &
                               narrays,  phi,        &
                               uee,      vnn,        &
                               HTE,      HTN,        &
                               tarea)
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
      use ice_work, only:  worka, workb
!
! !INPUT/OUTPUT PARAMETERS:
!
      integer (kind=int_kind), intent (in) ::     &
         nx_block, ny_block ,&! block dimensions
         nghost             ,&! number of ghost cells
         narrays              ! number of 2D arrays to be transported

      real (kind=dbl_kind), intent(in) ::         &
         dt                   ! time step

      real (kind=dbl_kind), dimension(nx_block,ny_block,narrays), &
         intent(inout) ::                                         &
         phi                  ! scalar field

      real (kind=dbl_kind), dimension(nx_block,ny_block),         &
         intent(in)::     &
         uee, vnn             ! cell edge velocities

      real (kind=dbl_kind), dimension(nx_block,ny_block), intent(in) :: &
         HTE                ,&! length of east cell edge 
         HTN                ,&! length of north cell edge
         tarea                ! grid cell area
!
!EOP
!
      integer (kind=int_kind) ::     &
         i, j, k, n         ,&! standard indices
         ilo,ihi,jlo,jhi      ! beginning and end of physical domain

      real (kind=dbl_kind) ::        &
         upwind, y1, y2, a, h   ! function

    !-------------------------------------------------------------------
    ! Define upwind function
    !-------------------------------------------------------------------

      upwind(y1,y2,a,h) = p5*dt*h*((a+abs(a))*y1+(a-abs(a))*y2)

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
            worka(i,j)=     &
               upwind(phi(i,j,n),phi(i+1,j,n),uee(i,j),HTE(i,j))
            workb(i,j)=     &
               upwind(phi(i,j,n),phi(i,j+1,n),vnn(i,j),HTN(i,j))
         enddo
         enddo

         do j = jlo, jhi
         do i = ilo, ihi
            phi(i,j,n) = phi(i,j,n) - ( worka(i,j)-worka(i-1,j)      &
                                      + workb(i,j)-workb(i,j-1) )    &
                                      / tarea(i,j)
         enddo
         enddo

      enddo                     ! narrays

      end subroutine upwind_field

!=======================================================================

      end module ice_transport_driver

!=======================================================================
