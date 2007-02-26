!=======================================================================
!BOP
!
! !MODULE: ice_transport_remap - horizontal transport via incremental remapping
!
! !DESCRIPTION:
!
! Transports quantities using the second-order conservative remapping
! scheme developed by John Dukowicz and John Baumgardner (DB) and modified
! for sea ice by William Lipscomb and Elizabeth Hunke.
!
! References:
!
! Dukowicz, J. K., and J. R. Baumgardner, 2000: Incremental
!  remapping as a transport/advection algorithm, J. Comput. Phys.,
!  160, 318-335.
!
! Lipscomb, W. H., and E. C. Hunke, 2004: Modeling sea ice
!  transport using incremental remapping, Mon. Wea. Rev., 132,
!  1341-1354.
!
! !REVISION HISTORY:
!  SVN:$Id: ice_transport_remap.F 33 2006-11-13 19:51:14Z eclare $
!
! authors William H. Lipscomb, LANL
!         John Baumgardner, LANL
!
! 2003: Vectorized by Clifford Chen (Fujitsu) and William Lipscomb
! 2004-05: Block structure added by William Lipscomb
!
! !INTERFACE:
!
      module ice_transport_remap
!
! !USES:
!
      use ice_kinds_mod
      use ice_communicate, only: my_task, master_task
      use ice_domain_size
      use ice_constants
      use ice_fileunits, only: nu_diag
!
!EOP
!
      implicit none
      save
      private
      public :: init_remap, horizontal_remap, make_masks

! NOTE: For remapping, hice, hsno, qice, and qsno are considered tracers.
!       ntrace is not equal to ntrcr!

      integer (kind=int_kind), parameter, public ::         &
         ntrace = 2+ntrcr+nilyr+nslyr  ! hice,hsno,qice,qsno,trcr
                          
      integer (kind=int_kind), dimension (ntrace), public ::     &
         tracer_type       ,&! = 1, 2, or 3 (see comments below)
         depend              ! tracer dependencies (see below)

      logical (kind=log_kind), dimension (ntrace) ::     &
         has_dependents      ! true if a tracer has dependent tracers

      integer (kind=int_kind), parameter ::     &
         ngroups  = 5        ! number of groups of triangles that
                             ! contribute transports across each edge

      integer (kind=int_kind), parameter ::     &
         integral_order = 3   ! polynomial order of quadrature integrals
                              ! linear=1, quadratic=2, cubic=3

      ! for triangle integral formulas
      real (kind=dbl_kind), parameter ::  & 
         p5625m = -9._dbl_kind/16._dbl_kind    ,&
         p52083 = 25._dbl_kind/48._dbl_kind

      logical (kind=log_kind), parameter ::     &
         l_dp_midpt = .true.  ! if true, find departure points using
                              ! corrected midpoint velocity

!=======================================================================
!
! The following comments are intended for users who want to 
! transport additional variables in CICE, or who want to apply 
! incremental remapping in another model.  For further details, 
! see Lipscomb and Hunke (2004). 
! 
! The incremental remapping algorithm is designed to solve the horizontal 
! transport equations in CICE.  These equations may be written as follows: 
! 
!               da/dt = del*(u*a)        (1) 
! dv/dt =   d(a*h)/dt = del*(u*a*h)      (2) 
! de/dt = d(a*h*q)/dt = del*(u*a*h*q)    (3) 
!            d(aT)/dt = del*(u*a*t)      (4) 
! 
! where d is a partial derivative, del is the 2D divergence operator, 
! a = fractional ice area, v = ice or snow volume, h = v/a = thickness, 
! e = ice or snow internal energy (J/m^2), q = e/v = internal energy per 
! unit volume (J/m^3), and T is a generic tracer such as surface 
! temperature.  These equations represent conservation of ice area, 
! volume, internal energy, and area-weighted tracer, respectively. 
! 
! In CICE, a, v and e are prognostic quantities, whereas h and q are 
! diagnosed.  The remapping routine works with tracers, which means 
! that h and q must be derived from a, v, and e at the beginning 
! of the routine (in state_to_tracers). 
! 
! Previous versions of CICE have assumed fixed ice and snow density. 
! Beginning with CICE 4.0, the ice and snow density can be variable. 
! In this case, equations (2) and (3) are replaced by 
! 
! dv/dt =        d(a*h)/dt = del*(u*a*h)        (5) 
! dm/dt =    d(a*h*rho)/dt = del*(u*a*h*rho)    (6) 
! de/dt = d(a*h*rho*qm)/dt = del*(u*a*h*rho*qm) (7) 
! 
! where rho = density and qm = internal energy per unit mass (J/kg). 
! Eq. (6) expresses mass conservation, which in the variable-density 
! case is no longer equivalent to volume conservation. 
! 
! We could have introduced optional state variables mi = a*hi*rhoi 
!  and ms = a*hs*rhos for ice and snow, but then we would need 
!  a lot of extra logic. 
! It is simpler to introduce optional tracer arrays rhoi and/or rhos. 
!  These are held in the trcrn array, whose first (and only required) 
!  element is Tsf, the surface temperature. 
!
! The remapping routine is written to transport a generic mass-like 
! field (in CICE, the ice area) along with an arbitrary number of 
! tracers in two dimensions.  The velocity components are assumed 
! to lie at grid cell corners and the transported scalars at cell centers. 
! Incremental remapping has the following nice properties: 
! 
! (1) Tracer monotonicity is preserved.  That is, no new local 
!     extrema are produced in fields like ice thickness or internal 
!     energy. 
! (2) The reconstucted mass and tracer fields vary linearly in x and y. 
!     This means that remapping is 2nd-order accurate in space, 
!     except where horizontal gradients are limited to preserve 
!     monotonicity. 
! (3) There are economies of scale.  Transporting a single field 
!     is rather expensive, but additional fields have a relatively 
!     low marginal cost. 
! 
! The following generic conservation equations may be solve: 
! 
!            da/dt = del*(u*a) 
!       d(a*T1)/dt = del*(u*a*T1)          (8) 
!    d(a*T1*T2)/dt = del*(u*a*T1*T2)       (9) 
! d(a*T1*T2*T3)/dt = del*(u*a*T1*T2*T3)    (10) 
! 
! which are of the same form as eq. (1)-(7). 
! Tracers satisfying equations of the form (1) are called "type 1." 
! The paradigmatic type 1 tracers are hi and hs. 
! 
! Tracers satisfying equations of the form (2) are called "type 2". 
! The paradigmatic type 2 tracers are qi and qs (or rhoi and rhos 
!  in the variable-density case). 
! 
! Tracers satisfying equations of the form (3) are called (you 
! guessed it) "type 3."  The paradigmatic type 3 tracers are qmi 
! and qms in the variable-density case.  There are no such tracers 
! in the constant-density case. 
! 
! The fields a, T1, and T2 are reconstructed in each grid cell with 
! 2nd-order accuracy.  T3 is reconstructed with 1st-order accuracy 
! (i.e., it is transported in upwind fashion) in order to avoid 
! additional mathematical complexity. 
! 
! The mass-like field lives in the array "aim" (shorthand for mean 
! ice area) and the tracers fields in the array "trm" (mean tracers). 
! In order to transport tracers correctly, the remapping routine 
! needs to know the tracers types and relationships.  This is done 
! as follows: 
! 
! Each field in the "trm" array is assigned an index, 1:ntrace. 
! (Note: ntrace is not the same as ntrcr, the number of tracers 
! in the trcrn state variable array.  For remapping purposes we 
! have additional tracers hi, hs, qi and qs.) 
! For standard CICE with ntrcr = 1, nilyr = 4, and nslyr = 1, 
! the indexing is as follows: 
! 1   = hi 
! 2   = hs 
! 3   = Ts 
! 4-7 = qi 
! 8   = qs 
! 
! The tracer types (1,2,3) are contained in the "tracer_type" array. 
! For standard CICE: 
! 
!     tracer_type = (1 1 1 2 2 2 2 2) 
! 
! Type 2 and type 3 tracers are said to depend on type 1 tracers. 
! For instance, qi depends on hi, which is to say that 
! there is a conservation equation of the form (3). 
! Thus we define a "depend" array.  For standard CICE: 
! 
!          depend = (0 0 0 1 1 1 1 2) 
! 
! which implies that elements 1-3 (hi, hs, Ts) are type 1, 
! elements 4-7 (qi) depend on element 1 (hi), and element 8 (qs) 
! depends on element 2 (hs). 
!
! We also define a logical array "has_dependents".  In standard CICE: 
! 
!  has_dependents = (T T F F F F F F), 
! 
! which means that only elements 1 and 2 (hi and hs) have dependent 
! tracers. 
! 
! For the variable-density case, things are a bit more complicated. 
! Suppose we have 4 variable-density ice layers and one variable- 
! density snow layers.  Then the indexing is as follows: 
! 1    = hi 
! 2    = hs 
! 3    = Ts 
! 4-7  = rhoi 
! 8    = rhos 
! 9-12 = qmi 
! 13   = qms 
! 
! The key arrays are: 
! 
!    tracer_type = (1 1 1 2 2 2 2 2 3 3 3 3 3) 
! 
!         depend = (0 0 0 1 1 1 1 2 4 5 6 7 8) 
! 
! has_dependents = (T T F T T T T T F F F F F) 
! 
! which imply that hi and hs are type 1 with dependents rhoi and rhos, 
! while rhoi and rhos are type 2 with dependents qmi and qms. 
! 
! Tracers added to the ntrcr array are handled automatically 
! by the remapping with little extra coding.  It is necessary 
! only to provide the correct type and dependency information 
! in subroutine init_remap. 
!
!=======================================================================

      contains

!=======================================================================
!
!BOP
!
! !IROUTINE: init_remap - initialize grid quantities used by remapping
!
! !INTERFACE:
!
      subroutine init_remap
!
! !DESCRIPTION:
!
! Grid quantities used by the remapping transport scheme
!
! !REVISION HISTORY:
!
! author William H. Lipscomb, LANL
!
! !USES:
!
      use ice_boundary
      use ice_domain
      use ice_blocks
      use ice_grid, only: HTE, HTN, dxt, dyt, dxhy, dyhx, tarear, &
                          mne, mnw, mse, msw,                     &
                          xav, yav, xxav, xyav, yyav,             &
                          xxxav, xxyav, xyyav, yyyav
      use ice_exit
      use ice_state, only: trcr_depend
      use ice_timers
!
! !INPUT/OUTPUT PARAMETERS:
!
!EOP
!
      integer (kind=int_kind) ::     &
           i, j, k, iblk, nt, nt1  ! standard indices

      character (char_len) :: &
         bc                   ! boundary condition type (Dirichlet, Neumann)

    !-------------------------------------------------------------------
    ! Compute tracer type and dependency vectors
    ! NOTE: May need to change these if transporting
    !       a different set of tracers.  See comments above.
    !-------------------------------------------------------------------

      depend(1:2)         = 0   ! hice, hsno
      tracer_type(1:2)    = 1   ! no dependency
      
      k = 2

      do nt = 1, ntrcr
         depend(k+nt) = trcr_depend(nt)  ! 0 for ice area tracers
                                         ! 1 for ice volume tracers
                                         ! 2 for snow volume tracers
         if (trcr_depend(nt) == 0) then
            tracer_type(k+nt) = 1
         else                   ! trcr_depend = 1 or 2
            tracer_type(k+nt) = 2
         endif
      enddo

      k = k + ntrcr

      depend(k+1:k+nilyr) = 1        ! qice depends on hice
      tracer_type(k+1:k+nilyr) = 2 

      k = k + nilyr

      depend(k+1:k+nslyr) = 2         ! qsno depends on hsno
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
      enddo


      do iblk = 1, nblocks

      ! Compute some geometric quantities needed for non-rectangular grids.

         call remapping_matrix_elements                                   &
                               (nx_block,          ny_block,              &
                                nghost,                                   &
                                HTE    (:,:,iblk), HTN    (:,:,iblk),     &
                                dxhy   (:,:,iblk), dyhx   (:,:,iblk),     &
                                mne(:,:,:,:,iblk), mnw(:,:,:,:,iblk),     &
                                mse(:,:,:,:,iblk), msw(:,:,:,:,iblk))

      ! Construct mean values of geometric quantities over each grid cell.

         call geometric_means(nx_block,        ny_block,             &
                              nghost,          tarear(:,:,iblk),     &
                              HTE  (:,:,iblk), HTN   (:,:,iblk),     &
                              dxt  (:,:,iblk), dyt   (:,:,iblk),     &
                              xav  (:,:,iblk), yav   (:,:,iblk),     &
                              xxav (:,:,iblk), xyav  (:,:,iblk),     &
                              yyav (:,:,iblk),                       &
                              xxxav(:,:,iblk), xxyav (:,:,iblk),     &
                              xyyav(:,:,iblk), yyyav (:,:,iblk))
      enddo

      ! Compute ghost cell values

      call ice_timer_start(timer_bound)
      bc = 'Neumann'
      call update_ghost_cells (xav,              bndy_info,     &
                               field_loc_center, field_type_scalar, bc)
      call update_ghost_cells (yav,              bndy_info,     &
                               field_loc_center, field_type_scalar, bc)
      call update_ghost_cells (xxav,             bndy_info,     &
                               field_loc_center, field_type_scalar, bc)
      call update_ghost_cells (xyav,             bndy_info,     &
                               field_loc_center, field_type_scalar, bc)
      call update_ghost_cells (yyav,             bndy_info,     &
                               field_loc_center, field_type_scalar, bc)
      call update_ghost_cells (xxxav,            bndy_info,     &
                               field_loc_center, field_type_scalar, bc)
      call update_ghost_cells (xxyav,            bndy_info,     &
                               field_loc_center, field_type_scalar, bc)
      call update_ghost_cells (xyyav,            bndy_info,     &
                               field_loc_center, field_type_scalar, bc)
      call update_ghost_cells (yyyav,            bndy_info,     &
                               field_loc_center, field_type_scalar, bc)
      call ice_timer_stop(timer_bound)

      end subroutine init_remap

!=======================================================================
!BOP
!
! !IROUTINE: geometric_means - geometric means used for remapping
!
! !INTERFACE:
!
      subroutine geometric_means(nx_block,  ny_block,   &
                                 nghost,    tarear,     &
                                 HTE,       HTN,        &
                                 dxt,       dyt,        &
                                 xav,       yav,        &
                                 xxav,      xyav,       &
                                 yyav,                  &
                                 xxxav,     xxyav,      &
                                 xyyav,     yyyav)
!
! !DESCRIPTION:
!
! In each grid cell, compute various geometric means used for remapping.
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
         nx_block, ny_block ,&! block dimensions
         nghost               ! number of ghost cells

      real (kind=dbl_kind), dimension(nx_block,ny_block),     &
         intent(in) ::     &
         tarear         ,&! 1 / grid cell area
         HTN            ,&! length of northern edge of T-cell (m)
         HTE            ,&! length of eastern edge of T-cell (m)
         dxt            ,&! width of T-cell through the middle (m)
         dyt              ! height of T-cell through the middle (m)

      real (kind=dbl_kind), dimension(nx_block, ny_block),     &
         intent(out) ::         &
         xav,  yav             ,&! mean T-cell values of x, y
         xxav, xyav, yyav      ,&! mean T-cell values of xx, xy, yy
         xxxav,xxyav,xyyav,yyyav ! mean T-cell values of xxx, xxy, xyy, yyy

! Local variables

      integer (kind=int_kind) ::     &
         i, j, m        ,&! counting indices
         ilo,ihi,jlo,jhi  ! beginning and end of physical domain

      real (kind=dbl_kind) ::             &
           ar                            ,&! triangle area
           xv1, xv2, xv3, yv1, yv2, yv3  ,&! vertices
           x0, x1, x2, x3, y0, y1, y2, y3,&! interior points
           x0sq, x1sq, x2sq, x3sq        ,&! x0^2, etc.
           y0sq, y1sq, y2sq, y3sq

      !-----------------------------------------------------------------
      ! Construct mean values of geometric quantities over each grid cell,
      ! relative to the origin (0,0) at the geometric center of the cell.
      ! (The geometric center is located at the intersection of the
      ! line joining the midpoints of the north and south edges with
      ! the line joining the midpoints of the east and west edges.
      ! The intersection is assumed to form a right angle.)
      ! These mean values are used to compute the cell centroid, center of
      ! ice area, and centers of ice and snow volume.
      ! The calculation is done by summing contributions from each of
      ! four triangles, labeled N, E, S, and W.  These triangles are
      ! formed by connecting the geometric center to the four cell
      ! corners.  Integrals are computed using the method described
      ! in subroutine transport_integrals.
      !------------------------------------------------------------------

      ilo = 1 + nghost
      ihi = nx_block - nghost
      jlo = 1 + nghost
      jhi = ny_block - nghost

      ! Initialize
      do j = 1, ny_block
      do i = 1, nx_block
         xav  (i,j) = c0
         yav  (i,j) = c0
         xxav (i,j) = c0
         xyav (i,j) = c0
         yyav (i,j) = c0
         xxxav(i,j) = c0
         xxyav(i,j) = c0
         xyyav(i,j) = c0
         yyyav(i,j) = c0
      enddo
      enddo

      ! Compute vertices and area of each triangle.

      do m = 1, 4
            
         do j = jlo,jhi
         do i = ilo,ihi

            ! xv1, yv1 = 0 for each triangle

            if (m==1) then
               xv2 =  p5*HTN(i,j) ! east triangle
               yv2 =  p5*HTE(i,j)
               xv3 =  p5*HTN(i,j-1)
               yv3 = -yv2
               ar  = p25*HTE(i,j)*dxt(i,j)

            elseif (m==2) then
               xv2 = -p5*HTN(i,j) ! west triangle
               yv2 =  p5*HTE(i-1,j)
               xv3 = -p5*HTN(i,j-1)
               yv3 = -p5*HTE(i-1,j)
               ar  = p25*HTE(i-1,j)*dxt(i,j)
                     
            elseif (m==3) then
               xv2 = -p5*HTN(i,j) ! north triangle
               yv2 =  p5*HTE(i-1,j)
               xv3 =  p5*HTN(i,j)
               yv3 =  p5*HTE(i,j)
               ar  = p25*HTN(i,j)*dyt(i,j)

            elseif (m==4) then
               xv2 = -p5*HTN(i,j-1) ! south triangle
               yv2 = -p5*HTE(i-1,j)
               xv3 =  p5*HTN(i,j-1)
               yv3 = -p5*HTE(i,j)
               ar  = p25*HTN(i,j-1)*dyt(i,j)
            endif

            ! Contribution to means from each triangle (E, W, N, S)

            x0 = (xv2 + xv3) / c3 ! midpoint
            y0 = (yv2 + yv3) / c3
            x1 = p6*x0          ! other 3 points needed for integral
            x2 = p6*x0 + p4*xv2
            x3 = p6*x0 + p4*xv3
            y1 = p6*y0
            y2 = p6*y0 + p4*yv2
            y3 = p6*y0 + p4*yv3

            x0sq = x0*x0
            y0sq = y0*y0
            x1sq = x1*x1
            y1sq = y1*y1
            x2sq = x2*x2
            y2sq = y2*y2
            x3sq = x3*x3
            y3sq = y3*y3
            
            xav(i,j) = xav(i,j) + ar*x0
            yav(i,j) = yav(i,j) + ar*y0
            xxav(i,j) = xxav(i,j) + ar *               &
                                  ( p5625m * x0*x0     &
                                  + p52083 * (x1sq + x2sq + x3sq) )
            xyav(i,j) = xyav(i,j) + ar *               &
                                  ( p5625m * x0*y0     &
                                  + p52083 * (x1*y1 + x2*y2 + x3*y3) )
            yyav(i,j) = yyav(i,j) + ar *               &
                                  ( p5625m * y0*y0     &
                                  + p52083 * (y1sq + y2sq + y3sq) )
            xxxav(i,j) = xxxav(i,j) + ar *             &
                            ( p5625m * x0*x0sq         &
                            + p52083 * (x1*x1sq + x2*x2sq + x3*x3sq) )
            xxyav(i,j) = xxyav(i,j) + ar *             &
                            ( p5625m * x0sq*y0         &
                            + p52083 * (x1sq*y1 + x2sq*y2 + x3sq*y3) )
            xyyav(i,j) = xyyav(i,j) + ar *             &
                            ( p5625m * x0*y0sq         &
                            + p52083 * (x1*y1sq + x2*y2sq + x3*y3sq) )
            yyyav(i,j) = yyyav(i,j) + ar *             &
                            ( p5625m * y0*y0sq         &
                            + p52083 * (y1*y1sq + y2*y2sq + y3*y3sq) )

         enddo                  ! i
         enddo                  ! j
      enddo                     ! m (loop over 4 triangles)

      ! Divide by grid cell area

      do j = jlo, jhi
      do i = ilo, ihi
         xav  (i,j) = xav  (i,j) * tarear(i,j)
         yav  (i,j) = yav  (i,j) * tarear(i,j)
         xxav (i,j) = xxav (i,j) * tarear(i,j)
         xyav (i,j) = xyav (i,j) * tarear(i,j)
         yyav (i,j) = yyav (i,j) * tarear(i,j)
         xxxav(i,j) = xxxav(i,j) * tarear(i,j)
         xxyav(i,j) = xxyav(i,j) * tarear(i,j)
         xyyav(i,j) = xyyav(i,j) * tarear(i,j)
         yyyav(i,j) = yyyav(i,j) * tarear(i,j)
      enddo                     ! i
      enddo                     ! j

      end subroutine geometric_means

!=======================================================================
!BOP
!
! !IROUTINE: remapping_matrix_elements - matrix elements used for remapping
!
! !INTERFACE:
!
      subroutine remapping_matrix_elements(nx_block,  ny_block,    &
                                           nghost,                 &
                                           HTE,       HTN,         &
                                           dxhy,      dyhx,        &
                                           mne,       mnw,         &
                                           mse,       msw)

!
! !DESCRIPTION:
!
! Compute matrix elements needed to transform from a reference frame
! whose origin is at a cell corner (a U-cell) to a reference frame
! whose origin is at a neighboring cell center (a T-cell).
! Transformations are needed because the axes of these two reference
! frames are not parallel on a curved grid.
!
! !REVISION HISTORY:
!
! author William H. Lipscomb
!
! !USES:
!
      use ice_work, only: worka, workb, workc, workd

! !INPUT/OUTPUT PARAMETERS:
!
      integer (kind=int_kind), intent(in) ::    &                
         nx_block, ny_block, &! block dimensions
         nghost               ! number of ghost cells

      real (kind=dbl_kind), dimension(nx_block,ny_block),     &
         intent(in) ::     &
         HTN         ,&! length of northern edge of T-cell (m)
         HTE         ,&! length of eastern edge of T-cell (m)
         dxhy        ,&! 0.5*(HTE(i,j) - HTE(i-1,j))
         dyhx          ! 0.5*(HTN(i,j) - HTN(i,j-1))

      real (kind=dbl_kind), dimension(2,2,nx_block,ny_block),   &
         intent(out) ::     &
         mne, mnw    ,&! matrices used for coordinate transformations
         mse, msw      ! ne = northeast corner, nw = northwest, etc.

! Local variables

      integer (kind=int_kind) ::     &
         i, j
                           
      real (kind=dbl_kind) ::     &
           theta1, theta2, theta3, theta4  ! angles

      !-----------------------------------------------------------------
      ! Compute angles between the U-cell coordinate axes and T-cell
      ! coordinate axes.  The U-cell axes lie along the cell edges,
      ! which do not meet at right angles in the T-cell reference frame.
      ! The angles are defined as positive if the U-cell axes are obtained
      ! by a counterclockwise rotation from the T-cell axes.
      ! Four independent angles are needed.
      ! NOTE: Where theta2 and theta3 are set to zero, they are not
      !         needed because they lie outside the transport domain.
      !       OK to have HTN = 0 or HTE = 0 for ghost cells at the edge
      !         of the global domain if the boundary is closed.
      !         Ice velocities should not be needed at the edges of a
      !         closed (i.e., land-locked) domain.
      !
      ! You might ask: Why not compute matrix elements on physical cells
      ! and then do a ghost-cell update?  The problem is that for tripole
      ! updates, these matrix elements transform across the tripole cut
      ! in a complicated way for which there is no update routine.
      !-----------------------------------------------------------------

      do j = 1, ny_block
      do i = 1, nx_block
         worka(i,j) = c0           ! theta1, NE corner
         workb(i,j) = c0           ! theta2, NW corner
         workc(i,j) = c0           ! theta3, SW corner
         workd(i,j) = c0           ! theta4, SW corner

         if (HTN(i,j) > puny) worka(i,j) = atan( dxhy(i,j)/HTN(i,j))
         if (HTE(i,j) > puny) workd(i,j) = atan(-dyhx(i,j)/HTE(i,j)) 
      enddo
      enddo

      do j = 1, ny_block
      do i = nghost+1, nx_block
         workb(i,j) = atan(dyhx(i,j)/HTE(i-1,j)) 
      enddo
      enddo
        
      do j = nghost+1, ny_block
      do i = 1, nx_block
         workc(i,j) = atan(-dxhy(i,j)/HTN(i,j-1))
      enddo
      enddo

      !-----------------------------------------------------------------
      ! Compute matrix elements.
      !-----------------------------------------------------------------

      do j = 1, ny_block
      do i = 1, nx_block

         theta1 = worka(i,j)
         theta2 = workb(i,j)
         theta3 = workc(i,j)
         theta4 = workd(i,j)

         mne(1,1,i,j) =  cos(theta1)
         mne(2,1,i,j) =  sin(theta1)
         mne(1,2,i,j) = -sin(theta4)
         mne(2,2,i,j) =  cos(theta4)
         
         mnw(1,1,i,j) =  cos(theta1)
         mnw(2,1,i,j) =  sin(theta1)
         mnw(1,2,i,j) = -sin(theta2)
         mnw(2,2,i,j) =  cos(theta2)
         
         msw(1,1,i,j) =  cos(theta3)
         msw(2,1,i,j) =  sin(theta3)
         msw(1,2,i,j) = -sin(theta2)
         msw(2,2,i,j) =  cos(theta2)
         
         mse(1,1,i,j) =  cos(theta3)
         mse(2,1,i,j) =  sin(theta3)
         mse(1,2,i,j) = -sin(theta4)
         mse(2,2,i,j) =  cos(theta4)
         
      enddo                     ! i
      enddo                     ! j

      end subroutine remapping_matrix_elements

!=======================================================================
!BOP
!
! !IROUTINE: horizontal_remap - incremental remapping transport scheme
!
! !INTERFACE:
!
      subroutine horizontal_remap (dt,                    &
                                   uvel,        vvel,     &
                                   aim,         trm)
!
! !DESCRIPTION:

! Solve the transport equations for one timestep using the incremental
! remapping scheme developed by John Dukowicz and John Baumgardner (DB)
! and modified for sea ice by William Lipscomb and Elizabeth Hunke.
!
! This scheme preserves monotonicity of ice area and tracers.  That is,
! it does not produce new extrema.  It is second-order accurate in space,
! except where gradients are limited to preserve monotonicity. 
!
! !REVISION HISTORY:
!
! author William H. Lipscomb, LANL
! 2006: Moved driver (subroutine transport_remap) into separate module. 
!       Developed new remapping scheme that will supersede this one.
!       For now, keep this module for testing and back compatibility.
!
! !USES:
!
      use ice_boundary
      use ice_global_reductions
      use ice_domain
      use ice_blocks
      use ice_grid, only: HTE, HTN, dxt, dyt, dxhy, dyhx,     &
                          tarea, tarear, hm,                  &
                          mne, mnw, mse, msw,                 &
                          xav, yav, xxav, xyav, yyav,         &
                          xxxav, xxyav, xyyav, yyyav
      use ice_exit
      use ice_work, only: worka
      use ice_calendar, only: istep1
      use ice_timers
!
! !INPUT/OUTPUT PARAMETERS:
!
      real (kind=dbl_kind), intent(in) ::     &
         dt      ! time step

      real (kind=dbl_kind), intent(in),       &
                dimension(nx_block,ny_block,max_blocks) ::           &
         uvel       ,&! x-component of velocity (m/s)
         vvel         ! y-component of velocity (m/s)

      real (kind=dbl_kind), intent(inout),     &
         dimension (nx_block,ny_block,0:ncat,max_blocks) ::          &
         aim          ! mean ice category areas in each grid cell

      real (kind=dbl_kind), intent(inout),     &
         dimension (nx_block,ny_block,ntrace,ncat,max_blocks) ::     &
         trm              ! mean tracer values in each grid cell
!
!EOP
!

      ! local variables

      integer (kind=int_kind) ::     &
         i, j           ,&! horizontal indices
         iblk           ,&! block indices
         n                ! ice category index

      integer (kind=int_kind), dimension(0:ncat,max_blocks) ::     &
         icellsnc         ! number of cells with ice

      integer (kind=int_kind),     &
         dimension(nx_block*ny_block,0:ncat,max_blocks) ::     &
         indxinc, indxjnc   ! compressed i/j indices

      type (block) ::     &
         this_block       ! block information for current block

      real (kind=dbl_kind), dimension (nx_block,ny_block,max_blocks) ::     &
         dpx            ,&! x coordinates of departure points at cell corners
         dpy            ,&! y coordinates of departure points at cell corners
         extmask          ! extended ice mask

      real (kind=dbl_kind),      &
         dimension (nx_block,ny_block,0:ncat,max_blocks) ::     &
         aic            ,&! ice area at geometric center of cell
         aix, aiy      ,&! limited derivative of ice area wrt x and y
         aimask           ! = 1. if ice is present, = 0. otherwise

      real (kind=dbl_kind),      &
         dimension (nx_block,ny_block,ntrace,ncat,max_blocks) ::     &
         trc            ,&! tracer values at geometric center of cell
         trx, try       ,&! limited derivative of tracer wrt x and y
         trmask           ! = 1. if tracer is present, = 0. otherwise

      real (kind=dbl_kind), dimension (nx_block,ny_block,0:ncat) ::     &
         aiflxe, aiflxn   ! ice area transports across E and N cell edges

      real (kind=dbl_kind), dimension (nx_block,ny_block,ntrace,ncat) ::     &
         atflxe, atflxn   ! area*tracer transports across E and N cell edges

      real (kind=dbl_kind), dimension (nx_block,ny_block,ngroups) ::     &
         triarea        ,&! area of east-edge departure triangle
         xp0, yp0       ,&! x and y coordinates of special triangle points
         xp1, yp1       ,&! (need 4 points for triangle integrals)
         xp2, yp2       ,&
         xp3, yp3

      integer (kind=int_kind),     &
         dimension (nx_block,ny_block,ngroups) ::     &
         iflux          ,&! i index of cell contributing transport
         jflux            ! j index of cell contributing transport

      integer (kind=int_kind), dimension(ngroups,max_blocks) ::     &
         icellsng         ! number of cells with ice

      integer (kind=int_kind),     &
         dimension(nx_block*ny_block,ngroups,max_blocks) ::     &
         indxing, indxjng ! compressed i/j indices

      logical (kind=log_kind) ::     &
         l_stop           ! if true, abort the model

      integer (kind=int_kind) ::     &
         istop, jstop     ! indices of grid cell where model aborts

      l_stop = .false.
      istop = 0
      jstop = 0

!---!-------------------------------------------------------------------
!---! Remap the ice area and associated tracers.
!---! Remap the open water area (without tracers).
!---!-------------------------------------------------------------------

      do iblk = 1, nblocks

    !------------------------------------------------------------------- 
    ! Compute masks and count ice cells.
    ! Masks are used to prevent tracer values in cells without ice from
    !  being used to compute tracer gradients.
    !------------------------------------------------------------------- 

         call make_masks (nx_block,          ny_block,     &
                          nghost,     &
                          icellsnc(:,iblk),     &
                          indxinc(:,:,iblk), indxjnc(:,:,iblk),      &
                          aim(:,:,:,iblk),   aimask(:,:,:,iblk),     &
                          trm(:,:,:,:,iblk), trmask(:,:,:,:,iblk))

    !-------------------------------------------------------------------
    ! Construct linear fields, limiting gradients to preserve monotonicity.
    !-------------------------------------------------------------------

         ! open water

!echmod         if (icellsnc(0,iblk) > 0) &
         call construct_fields(nx_block,            ny_block,             &
                               nghost,              icellsnc (0,iblk),    &
                               indxinc  (:,0,iblk), indxjnc(:,0,iblk),    &
                               HTN      (:,:,iblk), HTE  (:,:,iblk),      &
                               hm       (:,:,iblk), xav  (:,:,iblk),      &
                               yav      (:,:,iblk), xxav (:,:,iblk),      &
                               xyav     (:,:,iblk), yyav (:,:,iblk),      &
                               xxxav    (:,:,iblk), xxyav(:,:,iblk),      &
                               xyyav    (:,:,iblk), yyyav(:,:,iblk),      &
                               dxt      (:,:,iblk), dyt  (:,:,iblk),      &
                               aim    (:,:,0,iblk), aic(:,:,0,iblk),      &
                               aix    (:,:,0,iblk), aiy(:,:,0,iblk),      &
                               aimask (:,:,0,iblk) )

         ! ice categories

         do n = 1, ncat

!echmod            if (icellsnc(n,iblk) > 0) &
            call construct_fields                                         &
                                (nx_block,            ny_block,           &
                                 nghost,              icellsnc (n,iblk),  &
                                 indxinc  (:,n,iblk), indxjnc(:,n,iblk),  &
                                 HTN      (:,:,iblk), HTE    (:,:,iblk),  &
                                 hm       (:,:,iblk), xav    (:,:,iblk),  &
                                 yav      (:,:,iblk), xxav   (:,:,iblk),  &
                                 xyav     (:,:,iblk), yyav   (:,:,iblk),  &
                                 xxxav    (:,:,iblk), xxyav  (:,:,iblk),  &
                                 xyyav    (:,:,iblk), yyyav  (:,:,iblk),  &
                                 dxt      (:,:,iblk), dyt    (:,:,iblk),  &
                                 aim    (:,:,n,iblk), aic  (:,:,n,iblk),  &
                                 aix    (:,:,n,iblk), aiy  (:,:,n,iblk),  &
                                 aimask (:,:,n,iblk),                     &
                                 trm  (:,:,:,n,iblk), trc(:,:,:,n,iblk),  &
                                 trx  (:,:,:,n,iblk), try(:,:,:,n,iblk),  &
                                 trmask(:,:,:,n,iblk) )

         enddo                  ! n

    !-------------------------------------------------------------------
    ! Given velocity field at cell corners, compute departure points
    ! of trajectories.
    !-------------------------------------------------------------------

         call departure_points(nx_block,          ny_block,               &
                               nghost,            dt,                     &
                               uvel   (:,:,iblk), vvel   (:,:,iblk),      &
                               dpx    (:,:,iblk), dpy    (:,:,iblk),      &
                               HTN    (:,:,iblk), HTE    (:,:,iblk),      &
                               dxt    (:,:,iblk), dyt    (:,:,iblk),      &
                               dxhy   (:,:,iblk), dyhx   (:,:,iblk),      &
                               mne(:,:,:,:,iblk), mnw(:,:,:,:,iblk),      &
                               mse(:,:,:,:,iblk), msw(:,:,:,:,iblk),      &
                               l_dp_midpt,        l_stop,                 &
                               istop,             jstop)

         if (l_stop) then
            this_block = get_block(blocks_ice(iblk),iblk)         
            write(nu_diag,*) 'istep1, my_task, iblk =',     &
                              istep1, my_task, iblk
            write (nu_diag,*) 'Global block:', this_block%block_id
            if (istop > 0 .and. jstop > 0)     &
                 write(nu_diag,*) 'Global i and j:',     &
                                  this_block%i_glob(istop),     &
                                  this_block%j_glob(jstop) 
            call abort_ice('remap transport: bad departure points')
         endif

      enddo                     ! iblk

      call ice_timer_start(timer_bound)
      call update_ghost_cells (dpx,                bndy_info,     &
                               field_loc_NEcorner, field_type_vector)
      call update_ghost_cells (dpy,                bndy_info,     &
                               field_loc_NEcorner, field_type_vector)
      call ice_timer_stop(timer_bound)

    !-------------------------------------------------------------------
    ! define extended ice mask
    !-------------------------------------------------------------------

      do iblk = 1, nblocks
         do j = 1, ny_block
         do i = 1, nx_block
            if (dpx(i,j,iblk)/=c0 .or. dpy(i,j,iblk)/=c0) then
               worka(i,j) = c1
            else
               worka(i,j) = c0
            endif
         enddo
         enddo

         do j = 1+nghost, ny_block-nghost
         do i = 1+nghost, nx_block-nghost
         if (worka(i-1,j+1)==c1 .or. worka(i,j+1)==c1 .or.      &
             worka(i+1,j+1)==c1 .or.     &
             worka(i-1,j)  ==c1 .or. worka(i,j)  ==c1 .or.      &
             worka(i+1,j)  ==c1 .or.     &
             worka(i-1,j-1)==c1 .or. worka(i,j-1)==c1 .or.      &
             worka(i+1,j-1)==c1 ) then
             extmask(i,j,iblk) = c1
         endif
         enddo
         enddo
      enddo

      call ice_timer_start(timer_bound)
      call update_ghost_cells (extmask,            bndy_info,   &
                               field_loc_NEcorner, field_type_vector)

    !-------------------------------------------------------------------
    ! Ghost cell updates for centroids and gradients
    !-------------------------------------------------------------------

      call update_ghost_cells (aic,              bndy_info,     &
                               field_loc_center, field_type_scalar)
      call update_ghost_cells (aix,              bndy_info,     &
                               field_loc_center, field_type_vector)
      call update_ghost_cells (aiy,              bndy_info,     &
                               field_loc_center, field_type_vector)

      call update_ghost_cells (trc,              bndy_info,     &
                               field_loc_center, field_type_scalar)
      call update_ghost_cells (trx,              bndy_info,     &
                               field_loc_center, field_type_vector)
      call update_ghost_cells (try,              bndy_info,     &
                               field_loc_center, field_type_vector)
      call ice_timer_stop(timer_bound)

      do iblk = 1, nblocks

    !-------------------------------------------------------------------
    ! Transports for east cell edges.
    !-------------------------------------------------------------------

    !-------------------------------------------------------------------
    ! Compute areas and vertices of departure triangles.
    !-------------------------------------------------------------------

         call locate_triangles_east     &
                                 (nx_block,           ny_block,              &
                                  nghost,                                    &
                                  extmask(:,:,iblk),  icellsng (:,iblk),     &
                                  indxing(:,:,iblk),  indxjng(:,:,iblk),     &
                                  dpx    (:,:,iblk),  dpy    (:,:,iblk),     &
                                  HTN    (:,:,iblk),  HTE    (:,:,iblk),     &
                                  mne(:,:,:,:,iblk),  mnw(:,:,:,:,iblk),     &
                                  mse(:,:,:,:,iblk),  msw(:,:,:,:,iblk),     &
                                  xp0,                xp1,     &
                                  xp2,                xp3,     &
                                  yp0,                yp1,     &
                                  yp2,                yp3,     &
                                  iflux,              jflux,   &
                                  triarea)

    !-------------------------------------------------------------------
    ! Given triangle vertices, compute coordinates of triangle points
    !  needed for transport integrals.
    !-------------------------------------------------------------------

         call triangle_coordinates (nx_block, ny_block,     &
                                    integral_order,     &
                                    icellsng (:,iblk),     &
                                    indxing(:,:,iblk),       &
                                    indxjng(:,:,iblk),     &
                                    xp0,      xp1,     &
                                    xp2,      xp3,     &
                                    yp0,      yp1,     &
                                    yp2,      yp3)

    !-------------------------------------------------------------------
    ! Compute the transport across east cell edges by summing contributions
    ! from each triangle.
    !-------------------------------------------------------------------

         ! open water
         call transport_integrals(nx_block,        ny_block,     &
                                  integral_order,     &
                                  triarea,         icellsng (:,iblk),      &
                                  indxing(:,:,iblk),indxjng(:,:,iblk),     &
                                  iflux,           jflux,   &
                                  xp0,             xp1,     &
                                  xp2,             xp3,     &
                                  yp0,             yp1,     &
                                  yp2,             yp3,     &
                                  aic(:,:,0,iblk), aix   (:,:,0,iblk),     &
                                  aiy(:,:,0,iblk), aiflxe(:,:,0))

         ! ice categories
         do n = 1, ncat
            call transport_integrals     &
                               (nx_block,          ny_block,              &
                                integral_order,                           &
                                triarea,           icellsng (:,iblk),     &
                                indxing(:,:,iblk), indxjng(:,:,iblk),     &
                                iflux,             jflux,   &
                                xp0,               xp1,     &
                                xp2,               xp3,     &
                                yp0,               yp1,     &
                                yp2,               yp3,     &
                                aic(:,:,  n,iblk), aix   (:,:,  n,iblk),     &
                                aiy(:,:,  n,iblk), aiflxe(:,:,  n),          &
                                trc(:,:,:,n,iblk), trx   (:,:,:,n,iblk),     &
                                try(:,:,:,n,iblk), atflxe(:,:,:,n))
         enddo

    !-------------------------------------------------------------------
    ! Repeat for north edges
    !-------------------------------------------------------------------

         call locate_triangles_north     &
                                 (nx_block,           ny_block,     &
                                  nghost,     &
                                  extmask(:,:,iblk),  icellsng (:,iblk),     &
                                  indxing(:,:,iblk),  indxjng(:,:,iblk),     &
                                  dpx    (:,:,iblk),  dpy    (:,:,iblk),     &
                                  HTN    (:,:,iblk),  HTE    (:,:,iblk),     &
                                  mne(:,:,:,:,iblk),  mnw(:,:,:,:,iblk),     &
                                  mse(:,:,:,:,iblk),  msw(:,:,:,:,iblk),     &
                                  xp0,                xp1,     &
                                  xp2,                xp3,     &
                                  yp0,                yp1,     &
                                  yp2,                yp3,     &
                                  iflux,              jflux,   &
                                  triarea)

         call triangle_coordinates (nx_block, ny_block,    &
                                    integral_order,        &
                                    icellsng (:,iblk),     &
                                    indxing(:,:,iblk),     &
                                    indxjng(:,:,iblk),     &
                                    xp0,      xp1,         &
                                    xp2,      xp3,         &
                                    yp0,      yp1,         &
                                    yp2,      yp3)

         ! open water
         call transport_integrals(nx_block,         ny_block,     &
                                  integral_order,                          &
                                  triarea,          icellsng (:,iblk),     &
                                  indxing(:,:,iblk),indxjng(:,:,iblk),     &
                                  iflux,            jflux,   &
                                  xp0,              xp1,     &
                                  xp2,              xp3,     &
                                  yp0,              yp1,     &
                                  yp2,              yp3,     &
                                  aic(:,:,0,iblk),  aix(:,:,0,iblk),     &
                                  aiy(:,:,0,iblk),  aiflxn(:,:,0))

         call update_fields (nx_block,           ny_block,     &
                             nghost,             tarear(:,:,iblk),           &
                             l_stop,     &
                             istop,              jstop,     &
                             aiflxe(:,:,0),      aiflxn(:,:,0),     &
                             aim   (:,:,0,iblk))

         if (l_stop) then
            this_block = get_block(blocks_ice(iblk),iblk)         
            write (nu_diag,*) 'istep1, my_task, iblk, cat =',     &
                               istep1, my_task, iblk, n
            write (nu_diag,*) 'Global block:', this_block%block_id
            if (istop > 0 .and. jstop > 0)     &
                 write(nu_diag,*) 'Global i and j:',     &
                                  this_block%i_glob(istop),     &
                                  this_block%j_glob(jstop) 
            call abort_ice ('ice remap_transport: negative area, open water')
         endif

         ! ice categories
         do n = 1, ncat
            call transport_integrals     &
                               (nx_block,           ny_block,     &
                                integral_order,     &
                                triarea,            icellsng (:,iblk),     &
                                indxing(:,:,iblk),  indxjng(:,:,iblk),     &
                                iflux,              jflux,   &
                                xp0,                xp1,     &
                                xp2,                xp3,     &
                                yp0,                yp1,     &
                                yp2,                yp3,     &
                                aic  (:,:,n,iblk),  aix  (:,:,n,iblk),     &
                                aiy  (:,:,n,iblk),  aiflxn(:,:,n),         &
                                trc(:,:,:,n,iblk),  trx(:,:,:,n,iblk),     &
                                try(:,:,:,n,iblk),  atflxn(:,:,:,n))

            call update_fields(nx_block,             ny_block,             &
                               nghost,               tarear(:,:,iblk),     &
                               l_stop,                                     &
                               istop,                jstop,                &
                               aiflxe(:,:,  n),      aiflxn(:,:,  n),      &
                               aim   (:,:,  n,iblk),                       &
                               atflxe(:,:,:,n),      atflxn(:,:,:,n),      &
                               trm   (:,:,:,n,iblk))

            if (l_stop) then
               this_block = get_block(blocks_ice(iblk),iblk)         
               write (nu_diag,*) 'istep1, my_task, iblk, cat =',     &
                                  istep1, my_task, iblk, n
               write (nu_diag,*) 'Global block:', this_block%block_id
               if (istop > 0 .and. jstop > 0)     &
                    write(nu_diag,*) 'Global i and j:',     &
                                     this_block%i_glob(istop),     &
                                     this_block%j_glob(jstop) 
               call abort_ice ('ice remap_transport: negative area, ice')
            endif
         enddo                  ! n

      enddo                     ! iblk

      end subroutine horizontal_remap

!=======================================================================
!
!BOP
!
! !IROUTINE: make_masks - make area and tracer masks
!
! !INTERFACE:
!
      subroutine make_masks (nx_block, ny_block,            &
                             nghost,                        &
                             icells,   indxi,    indxj,     &
                             aim,      aimask,              &
                             trm,      trmask)
!
! !DESCRIPTION:
!
! Make area and tracer masks.
!
! If an area is masked out (aim < puny), then the values of tracers
!  in that grid cell are assumed to have no physical meaning.
!
! Similarly, if a tracer with dependents is masked out
!  (abs(trm) < puny), then the values of its dependent tracers in that
!  grid cell are assumed to have no physical meaning.
! For example, the enthalpy value has no meaning if the thickness
!  is zero.
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
           nx_block, ny_block  ,&! block dimensions
           nghost                ! number of ghost cells

      integer (kind=int_kind), dimension(0:ncat),     &
           intent(out) ::     &
           icells         ! number of cells with ice

      integer (kind=int_kind), dimension(nx_block*ny_block,0:ncat),     &
           intent(out) ::     &
           indxi        ,&! compressed i/j indices
           indxj

      real (kind=dbl_kind), dimension (nx_block,ny_block,0:ncat),     &
           intent(in) ::     &
           aim            ! mean ice area in each grid cell

      real (kind=dbl_kind), dimension (nx_block,ny_block,0:ncat),     &
           intent(out) ::     &
           aimask         ! = 1. if ice is present, else = 0.

      real (kind=dbl_kind),      &
           dimension (nx_block, ny_block, ntrace, ncat),     &
           intent(in), optional ::     &
           trm            ! mean tracer values in each grid cell

      real (kind=dbl_kind),      &
           dimension (nx_block, ny_block, ntrace, ncat),     &
           intent(out), optional ::     &
           trmask         ! = 1. if tracer is present, else = 0.
!
!EOP
!
      integer (kind=int_kind) ::     &
           i, j, ij       ,&! horizontal indices
           n              ,&! ice category index
           nt             ,&! tracer index
           ilo,ihi,jlo,jhi  ! beginning and end of physical domain

      ilo = 1 + nghost
      ihi = nx_block - nghost
      jlo = 1 + nghost
      jhi = ny_block - nghost

      do n = 0, ncat
         do ij = 1, nx_block*ny_block
            indxi(ij,n) = 0
            indxj(ij,n) = 0
         enddo
      enddo

    !-------------------------------------------------------------------
    ! open water mask
    !-------------------------------------------------------------------

      icells(0) = 0
      do j = 1, ny_block
      do i = 1, nx_block
         if (aim(i,j,0) > puny) then
            aimask(i,j,0) = c1
            icells(0) = icells(0) + 1
            ij = icells(0)
            indxi(ij,0) = i
            indxj(ij,0) = j
         else
            aimask(i,j,0) = c0
         endif
      enddo
      enddo


      do n = 1, ncat

    !-------------------------------------------------------------------
    ! Find grid cells where ice is present.
    !-------------------------------------------------------------------

         icells(n) = 0
         do j = 1, ny_block
         do i = 1, nx_block
            if (aim(i,j,n) > puny) then
               icells(n) = icells(n) + 1
               ij = icells(n)
               indxi(ij,n) = i
               indxj(ij,n) = j
            endif               ! aim > puny
         enddo
         enddo

    !-------------------------------------------------------------------
    ! ice area mask
    !-------------------------------------------------------------------

         aimask(:,:,n) = c0
         do ij = 1, icells(n)
            i = indxi(ij,n)
            j = indxj(ij,n)
            aimask(i,j,n) = c1
         enddo

    !-------------------------------------------------------------------
    ! tracer masks
    !-------------------------------------------------------------------

         if (present(trm)) then

            trmask(:,:,:,n) = c0
            do nt = 1, ntrace
               if (has_dependents(nt)) then
                  do ij = 1, icells(n)
                     i = indxi(ij,n)
                     j = indxj(ij,n)
                     if (abs(trm(i,j,nt,n)) > puny) then
                        trmask(i,j,nt,n) = c1
                     endif
                  enddo
               endif
            enddo

         endif                     ! present(trm)

    !-------------------------------------------------------------------
    ! Redefine icells (without ghost cells)
    !-------------------------------------------------------------------

         icells(n) = 0
         do j = jlo, jhi
         do i = ilo, ihi
            if (aim(i,j,n) > puny) then
               icells(n) = icells(n) + 1
               ij = icells(n)
               indxi(ij,n) = i
               indxj(ij,n) = j
            endif               ! aim > puny
         enddo
         enddo
      
      enddo ! n

      end subroutine make_masks

!=======================================================================
!BOP
!
! !IROUTINE: departure_points - compute departure points of trajectories
!
! !INTERFACE:
!
      subroutine departure_points (nx_block,   ny_block,    &
                                   nghost,     dt,          &
                                   uvel,       vvel,        &
                                   dpx,        dpy,         &
                                   HTN,        HTE,         &
                                   dxt,        dyt,         &
                                   dxhy,       dyhx,        &
                                   mne,        mnw,         &
                                   mse,        msw,         &
                                   l_dp_midpt, l_stop,      &
                                   istop,      jstop)
!
! !DESCRIPTION:
!
! Given velocity fields on cell corners, compute departure points
! of trajectories using a midpoint approximation.
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
         nx_block, ny_block ,&! block dimensions
         nghost               ! number of ghost cells

      real (kind=dbl_kind), intent(in) ::     &
         dt               ! time step (s)

      real (kind=dbl_kind), dimension (nx_block,ny_block), intent(in) ::     &
         uvel           ,&! x-component of velocity (m/s)
         vvel           ,&! y-component of velocity (m/s)
         HTN            ,&! length of northern edge of T-cell (m)
         HTE            ,&! length of eastern edge of T-cell (m)
         dxt            ,&! width of T-cell through the middle (m)
         dyt            ,&! height of T-cell through the middle (m)
         dxhy           ,&! 0.5*(HTE(i,j) - HTE(i-1,j))
         dyhx             ! 0.5*(HTN(i,j) - HTN(i,j-1))

      real (kind=dbl_kind), dimension(2,2,nx_block,ny_block),     &
         intent(in) ::     &
         mne, mnw       ,&! matrices used for coordinate transformations
         mse, msw         ! ne = northeast corner, nw = northwest, etc.

      real (kind=dbl_kind), dimension (nx_block,ny_block),     &
         intent(out) ::     &
         dpx            ,&! x coordinates of departure points at cell corners
         dpy              ! y coordinates of departure points at cell corners

      logical (kind=log_kind), intent(in) ::     &
         l_dp_midpt          ! if true, find departure points using
                             ! corrected midpoint velocity

      logical (kind=log_kind), intent(inout) ::     &
         l_stop       ! if true, abort on return

      integer (kind=int_kind), intent(inout) ::     &
         istop, jstop     ! indices of grid cell where model aborts 

!
!EOP
!
      integer (kind=int_kind) ::     &
         i, j           ,&! horizontal indices
         i2, j2         ,&! horizontal indices
         niter          ,&! iteration counter
         ilo,ihi,jlo,jhi  ! beginning and end of physical domain

      real (kind=dbl_kind) ::     &
         mpx, mpy       ,&! coordinates of midpoint of back trajectory
         u1t, v1t       ,&! transformed velocity, SW corner
         u2t, v2t       ,&! transformed velocity, SE corner
         u3t, v3t       ,&! transformed velocity, NE corner
         u4t, v4t       ,&! transformed velocity, NW corner
         umpt, vmpt     ,&! midpoint velocity in transformed coordinates
         ump, vmp       ,&! midpoint velocity in original corner coordinates
         w1,w2,w3,w4      ! work variables

      real (kind=dbl_kind), dimension (nx_block,ny_block) ::     &
         cxt, cyt       ,&! transformed cell corner coordinates
         mpxt, mpyt     ,&! midpoint transformed to cell-ctr coordinates
         mpxs, mpys     ,&! midpoint in stretched coordinates
         mat11, mat12   ,&! transformation matrix, cell corner to center
         mat21, mat22     ! transformation matrix, cell corner to center

      integer (kind=int_kind), dimension (nx_block,ny_block) ::     &
         hindi, hindj     ! horizontal indices array

    !-------------------------------------------------------------------
    ! Initialize
    !-------------------------------------------------------------------

      ilo = 1 + nghost
      ihi = nx_block - nghost
      jlo = 1 + nghost
      jhi = ny_block - nghost

    !-------------------------------------------------------------------
    ! Estimate departure points.
    ! This estimate is 1st-order accurate; improve accuracy by
    !  using midpoint approximation below.
    !-------------------------------------------------------------------

      do j = 1, ny_block
      do i = 1, nx_block
         dpx(i,j) = -dt*uvel(i,j)
         dpy(i,j) = -dt*vvel(i,j)
      enddo
      enddo

    !-------------------------------------------------------------------
    ! Check for values out of bounds (more than one grid cell away).
    !-------------------------------------------------------------------


      do j = jlo, jhi
      do i = ilo, ihi
         if ((dpx(i,j) < -HTN(i,j)) .or.      &
             (dpx(i,j) >  HTN(i+1,j)) .or.     &
             (dpy(i,j) < -HTE(i,j)) .or.     &
             (dpy(i,j) >  HTE(i,j+1))) then
            l_stop = .true.
            istop = i
            jstop = j
         endif
      enddo
      enddo

      if (l_stop) then
         i = istop
         j = jstop
         write (nu_diag,*) ' '
         write (nu_diag,*)     &
                    'Warning: Departure points out of bounds in remap'
         write (nu_diag,*) 'i, j =', i, j
         write (nu_diag,*) 'dpx, dpy =', dpx(i,j), dpy(i,j)
         return
      endif

      if (l_dp_midpt) then ! find dep pts using corrected midpt velocity 

         do j = jlo, jhi
         do i = ilo, ihi
            if (uvel(i,j)/=c0 .or. vvel(i,j)/=c0) then

    !-------------------------------------------------------------------
    ! Estimate midpoint of backward trajectory relative to corner (i,j).
    !-------------------------------------------------------------------
            mpx = p5 * dpx(i,j)
            mpy = p5 * dpy(i,j)

    !-------------------------------------------------------------------
    ! Determine the indices (i2,j2) of the cell where the trajectory lies
    ! and compute for that cell:
    ! (1) the matrix 'mat' needed to transform vectors from the reference
    !     frame of cell corner (i,j) to that of cell center (i2,j2)
    ! (2) the coordinates (cxt,cyt) of corner (i,j) relative to center
    !     (i2,j2) in the (i2,j2) reference frame
    ! (3) a rough guess for the midpoint location in stretched coordinates
    !     (mpxs, mpys), used in the first pass of the iteration below
    ! Note: Coordinates in the (i2,j2) reference frame have a 't' at
    !       the end.
    !-------------------------------------------------------------------

            if (mpx >= c0 .and. mpy >= c0) then ! cell (i+1,j+1)
               i2 = i+1
               j2 = j+1
               hindi(i,j) = i2
               hindj(i,j) = j2
               mat11(i,j) = msw(1,1,i2,j2)
               mat21(i,j) = msw(2,1,i2,j2)
               mat12(i,j) = msw(1,2,i2,j2)
               mat22(i,j) = msw(2,2,i2,j2)
               cxt(i,j)   = -p5*HTN(i2,j2-1)
               cyt(i,j)   = -p5*HTE(i2-1,j2)
               mpxs(i,j)  = -c1
               mpys(i,j)  = -c1
            elseif (mpx < c0 .and. mpy < c0) then ! cell (i,j)
               i2 = i
               j2 = j
               hindi(i,j) = i2
               hindj(i,j) = j2
               mat11(i,j) = mne(1,1,i2,j2)
               mat21(i,j) = mne(2,1,i2,j2)
               mat12(i,j) = mne(1,2,i2,j2)
               mat22(i,j) = mne(2,2,i2,j2)
               cxt(i,j)   = p5*HTN(i2,j2)
               cyt(i,j)   = p5*HTE(i2,j2)
               mpxs(i,j)  = c1
               mpys(i,j)  = c1
            elseif (mpx >= c0 .and. mpy < c0) then ! cell (i+1,j)
               i2 = i+1
               j2 = j
               hindi(i,j) = i2
               hindj(i,j) = j2
               mat11(i,j) = mnw(1,1,i2,j2)
               mat21(i,j) = mnw(2,1,i2,j2)
               mat12(i,j) = mnw(1,2,i2,j2)
               mat22(i,j) = mnw(2,2,i2,j2)
               cxt(i,j)   = -p5*HTN(i2,j2)
               cyt(i,j)   =  p5*HTE(i2-1,j2)
               mpxs(i,j)  = -c1
               mpys(i,j)  =  c1
            elseif (mpx < c0 .and. mpy >= c0) then ! cell (i,j+1)
               i2 = i
               j2 = j+1
               hindi(i,j) = i2
               hindj(i,j) = j2
               mat11(i,j) = mse(1,1,i2,j2)
               mat21(i,j) = mse(2,1,i2,j2)
               mat12(i,j) = mse(1,2,i2,j2)
               mat22(i,j) = mse(2,2,i2,j2)
               cxt(i,j)   =  p5*HTN(i2,j2-1)
               cyt(i,j)   = -p5*HTE(i2,j2)
               mpxs(i,j)  =  c1
               mpys(i,j)  = -c1
            endif
            
    !-------------------------------------------------------------------
    ! Transform coordinates of the trajectory midpoint to the (i2,j2)
    ! reference frame.
    !-------------------------------------------------------------------

            mpxt(i,j) = cxt(i,j) + mat11(i,j)*mpx + mat12(i,j)*mpy
            mpyt(i,j) = cyt(i,j) + mat21(i,j)*mpx + mat22(i,j)*mpy

    !-------------------------------------------------------------------
    ! Transform (mpxt,mpyt) to a stretched coordinate system in which
    ! the coordinates of the four corners relative to the center are
    ! (-1,-1), (1,-1), (1,1), and (-1,1).
    !
    ! Iterate a couple of times for accuracy.
    ! (Occasionally abs(mpxs) or abs(mpys) > 1 after first iteration.)
    !-------------------------------------------------------------------

               i2 = hindi(i,j)
               j2 = hindj(i,j)
               w1 = c2*mpxt(i,j) * dyt(i2,j2)
               w3 = c2*mpyt(i,j) * dxt(i2,j2)
               ! 1st iteration
               w2 =   dxt(i2,j2) * dyt(i2,j2)     &
                   + dyhx(i2,j2) * (c2*mpyt(i,j)     &
                               - mpxs(i,j)*mpys(i,j)*dxhy(i2,j2))
               w4 =   dxt(i2,j2) * dyt(i2,j2)     &
                   + dxhy(i2,j2) * (c2*mpxt(i,j)     &
                               - mpxs(i,j)*mpys(i,j)*dyhx(i2,j2))
               mpxs(i,j) = w1/w2
               mpys(i,j) = w3/w4
               ! 2nd iteration
               w2 =   dxt(i2,j2) * dyt(i2,j2)     &
                   + dyhx(i2,j2) * (c2*mpyt(i,j)     &
                               - mpxs(i,j)*mpys(i,j)*dxhy(i2,j2))
               w4 =   dxt(i2,j2) * dyt(i2,j2)     &
                   + dxhy(i2,j2) * (c2*mpxt(i,j)     &
                               - mpxs(i,j)*mpys(i,j)*dyhx(i2,j2))
               mpxs(i,j) = w1/w2
               mpys(i,j) = w3/w4
               ! 3rd iteration
               w2 =   dxt(i2,j2) * dyt(i2,j2)     &
                   + dyhx(i2,j2) * (c2*mpyt(i,j)     &
                               - mpxs(i,j)*mpys(i,j)*dxhy(i2,j2))
               w4 =   dxt(i2,j2) * dyt(i2,j2)     &
                   + dxhy(i2,j2) * (c2*mpxt(i,j)     &
                               - mpxs(i,j)*mpys(i,j)*dyhx(i2,j2))
               mpxs(i,j) = w1/w2
               mpys(i,j) = w3/w4

    !-------------------------------------------------------------------
    ! Transform the four corner velocities to the (i2,j2) reference frame.
    !-------------------------------------------------------------------

            i2   = hindi(i,j)
            j2   = hindj(i,j)

            u1t = msw(1,1,i2,j2)*uvel(i2-1,j2-1)     &
                + msw(1,2,i2,j2)*vvel(i2-1,j2-1)
            v1t = msw(2,1,i2,j2)*uvel(i2-1,j2-1)     &
                + msw(2,2,i2,j2)*vvel(i2-1,j2-1)

            u2t = mse(1,1,i2,j2)*uvel(i2,j2-1)     &
                + mse(1,2,i2,j2)*vvel(i2,j2-1)
            v2t = mse(2,1,i2,j2)*uvel(i2,j2-1)     &
                + mse(2,2,i2,j2)*vvel(i2,j2-1)

            u3t = mne(1,1,i2,j2)*uvel(i2,j2)     &
                + mne(1,2,i2,j2)*vvel(i2,j2)
            v3t = mne(2,1,i2,j2)*uvel(i2,j2)     &
                + mne(2,2,i2,j2)*vvel(i2,j2)

            u4t = mnw(1,1,i2,j2)*uvel(i2-1,j2)     &
                + mnw(1,2,i2,j2)*vvel(i2-1,j2)
            v4t = mnw(2,1,i2,j2)*uvel(i2-1,j2)     &
                + mnw(2,2,i2,j2)*vvel(i2-1,j2)

    !-------------------------------------------------------------------
    ! Using a bilinear approximation, estimate the velocity at the
    ! trajectory midpoint in the (i2,j2) reference frame.
    !-------------------------------------------------------------------

            umpt = p25 * ( u1t*(mpxs(i,j)-c1)*(mpys(i,j)-c1)     &
                         - u2t*(mpxs(i,j)+c1)*(mpys(i,j)-c1)     &
                         + u3t*(mpxs(i,j)+c1)*(mpys(i,j)+c1)     &
                         - u4t*(mpxs(i,j)-c1)*(mpys(i,j)+c1) )

            vmpt = p25 * ( v1t*(mpxs(i,j)-c1)*(mpys(i,j)-c1)     &
                         - v2t*(mpxs(i,j)+c1)*(mpys(i,j)-c1)     &
                         + v3t*(mpxs(i,j)+c1)*(mpys(i,j)+c1)     &
                         - v4t*(mpxs(i,j)-c1)*(mpys(i,j)+c1) )

    !-------------------------------------------------------------------
    ! Transform the velocity back to the cell corner reference frame,
    ! using the inverse of matrix 'mat'.
    !-------------------------------------------------------------------

            w1 = c1 / (mat11(i,j)*mat22(i,j) - mat12(i,j)*mat21(i,j))
            ump = w1 * ( mat22(i,j)*umpt - mat12(i,j)*vmpt)
            vmp = w1 * (-mat21(i,j)*umpt + mat11(i,j)*vmpt)

    !-------------------------------------------------------------------
    ! Use the midpoint velocity to estimate the coordinates of the
    ! departure point relative to corner (i,j).
    !-------------------------------------------------------------------

            dpx(i,j) = -dt * ump
            dpy(i,j) = -dt * vmp

            endif               ! nonzero velocity
         enddo                  ! i
         enddo                  ! j

      endif                     ! l_dp_mipdt

      end subroutine departure_points

!=======================================================================
!
!BOP
!
! !IROUTINE: locate_triangles_east - triangle info for east edges
!
! !INTERFACE:
!
      subroutine locate_triangles_east (nx_block,  ny_block,   &
                                        nghost,                &
                                        extmask,  icells,      &
                                        indxi,  indxj,         &
                                        dpx,       dpy,        &
                                        HTN,       HTE,        &
                                        mne,       mnw,        &
                                        mse,       msw,        &
                                        xp0,       xp1,        &
                                        xp2,       xp3,        &
                                        yp0,       yp1,        &
                                        yp2,       yp3,        &
                                        iflux,     jflux,      &
                                        triarea)
!
! !DESCRIPTION:
!
! Compute areas and vertices of transport triangles for east cell edges.
!
! !REVISION HISTORY:
!
! authors William H. Lipscomb, LANL
!         John R. Baumgardner, LANL
!
! !USES:
!
! !INPUT/OUTPUT PARAMETERS:
!
      integer (kind=int_kind), intent(in) ::     &
         nx_block, ny_block ,&! block dimensions
         nghost               ! number of ghost cells

      real (kind=dbl_kind), dimension (nx_block,ny_block),     &
         intent(in) ::     &
         extmask          ! extended ice extent mask

      real (kind=dbl_kind), dimension (nx_block,ny_block),     &
           intent(in) ::     &
         dpx            ,&! x coordinates of departure points at cell corners
         dpy            ,&! y coordinates of departure points at cell corners
         HTN            ,&! length of northern edge of T-cell (m)
         HTE              ! length of southern edge of T-cell

      real (kind=dbl_kind), dimension(2,2,nx_block,ny_block),      &
            intent(in) ::     &
         mne, mnw       ,&! matrices used for coordinate transformations
         mse, msw         ! ne = northeast corner, nw = northwest, etc.

      real (kind=dbl_kind), dimension (nx_block, ny_block, ngroups),     &
           intent(out) ::     &
         triarea        ,&! area of east-edge departure triangle
         xp0, yp0       ,&! coordinates of special triangle points
         xp1, yp1       ,&
         xp2, yp2       ,&
         xp3, yp3

      integer (kind=int_kind), intent(out),     &
           dimension (nx_block, ny_block, ngroups) ::     &
         iflux          ,&! i index of cell contributing east transport
         jflux            ! j index of cell contributing east transport

      integer (kind=int_kind), dimension (ngroups), intent(out) ::     &
         icells              ! number of cells where triarea > puny

      integer (kind=int_kind), dimension (nx_block*ny_block,ngroups),     &
         intent(out) ::     &
         indxi ,&! compressed index in i-direction
         indxj   ! compressed index in j-direction
!
!EOP
!
      integer (kind=int_kind) ::     &
         i, j           ,&! horizontal indices of cell edge
         i2, j2         ,&! horizontal indices of cell contributing transport
         ng             ,&! triangle group index
         ilo,ihi,jlo,jhi,&! beginning and end of physical domain
         ij               ! compressed horizontal index

      real (kind=dbl_kind) ::     &
         x1, y1, x2, y2 ,&! x,y coordinates of departure points, as in DB
         xa, ya, xb, yb ,&! x,y coordinates of points where the lines joining
                          ! (x1,y1) and (x2,y2) cross cell edges, as in DB
         xca, yca       ,&! transformed coordinates of corner point a
         xda, yda       ,&! transformed coordinates of departure point a
         xxa, yxa       ,&! transformed coordinates of intersection point xa
         xya, yya       ,&! transformed coordinates of intersection point ya
         xcb, ycb       ,&! transformed coordinates of corner point b
         xdb, ydb       ,&! transformed coordinates of departure point b
         xxb, yxb       ,&! transformed coordinates of intersection point xb
         xyb, yyb       ,&! transformed coordinates of intersection point yb
         xic, yic       ,&! transformed coordinates of point where the
                          ! line joining dep pts intersects the face
         w1               ! work variable

      integer (kind=int_kind), dimension (nx_block,ny_block,ngroups) ::  &
         fluxsign         ! = 1 for positive flux, -1 for negative

      logical :: cnd1, cnd2, cnd3   ! conditionals

    !-------------------------------------------------------------------
    ! Triangle notation:
    ! For each edge, there are 20 triangles that can contribute,
    ! but many of these are mutually exclusive.  It turns out that
    ! at most 5 triangles can contribute to transport integrals at once.
    !
    ! For the east edge, these triangles are referred to as:
    ! (1) NE, NE1, N, N2
    ! (2) SE, SE1, S, S2
    ! (3) NE2, N1, SE2, S1
    ! (4) H1a, H1b, E1a, E2b
    ! (5) H2a, H2b, N2a, N2b
    !
    ! See Figure 3 in DB for pictures of these triangles.
    ! See Table 1 in DB for logical conditions.
    !
    ! Many triangle vertices lie at points whose coordinates are
    ! (x1,y1), (x2,y2), (xa,0), (xb,0), (0,ya), and (0,yb) in a
    ! reference frame whose origin is the cell corner.  These
    ! coordinates must be transformed to the reference frame whose
    ! origin is the geometric center of the T-cell in which the triangle
    ! is located.  The transformation is carried out using pre-computed
    ! 2x2 matrices.  There are 4 matrices (one for each corner)
    ! associated with each grid cell.  They do not describe pure
    ! rotations, because they do not preserve length.
    !-------------------------------------------------------------------

    !-------------------------------------------------------------------
    ! Initialize
    !-------------------------------------------------------------------

      ilo = 1 + nghost
      ihi = nx_block - nghost
      jlo = 1 + nghost
      jhi = ny_block - nghost

      do ng = 1, ngroups
      do j = 1, ny_block
      do i = 1, nx_block
         fluxsign(i,j,ng) = 0
         iflux   (i,j,ng) = 0
         jflux   (i,j,ng) = 0
         xp0     (i,j,ng) = c0
         xp1     (i,j,ng) = c0
         xp2     (i,j,ng) = c0
         xp3     (i,j,ng) = c0
         yp0     (i,j,ng) = c0
         yp1     (i,j,ng) = c0
         yp2     (i,j,ng) = c0
         yp3     (i,j,ng) = c0
      enddo
      enddo
      enddo

    !-------------------------------------------------------------------
    ! Main loop
    !-------------------------------------------------------------------

      do j = jlo, jhi
      do i = ilo-1, ihi      ! includes W edge of cells with index ilo

         if (extmask(i,j) > puny) then

    !-------------------------------------------------------------------
    ! coordinates of departure points
    !-------------------------------------------------------------------
         x1 = dpx(i,j)
         y1 = dpy(i,j)
         x2 = dpx(i,j-1)
         y2 = dpy(i,j-1)
         w1 =  c1 / (y2 - HTE(i,j)  - y1)
         xa = (x1*(y2 - HTE(i,j)) - x2*y1) * w1
         xb = (x1*y2 - x2*(y1 + HTE(i,j))) * w1
         if (abs(xb-xa) > puny) then
            ya = xa * HTE(i,j) / (xb - xa)
            yb = ya + HTE(i,j)
         else
            ya = c0
            yb = c0
         endif

    !-------------------------------------------------------------------
    ! contribution from triangles in NE cell
    ! (All corner cells follow the same pattern.)
    !-------------------------------------------------------------------
         ! load horizontal indices of NE cell
         i2 = i+1
         j2 = j+1

         ! find vertex coordinates relative to center of NE cell

         xca = -p5*HTN(i2,j2-1)                         ! corner pt
         yca = -p5*HTE(i2-1,j2)

         xda = xca + msw(1,1,i2,j2)*x1 + msw(1,2,i2,j2)*y1  ! departure pt
         yda = yca + msw(2,1,i2,j2)*x1 + msw(2,2,i2,j2)*y1

         xxa = xca + msw(1,1,i2,j2)*xa                ! xa
         yxa = yca + msw(2,1,i2,j2)*xa

         xya = xca + msw(1,2,i2,j2)*ya                ! ya
         yya = yca + msw(2,2,i2,j2)*ya

         ! vertices of 2 potential group 1 triangles
         ng = 1

         cnd1 = xa > c0 .and. y1 > c0 .and. x1 >= c0   ! NE (group 1)
         if (cnd1) then
            xp1     (i,j,ng) = xca
            yp1     (i,j,ng) = yca
            xp2     (i,j,ng) = xxa
            yp2     (i,j,ng) = yxa
            xp3     (i,j,ng) = xda
            yp3     (i,j,ng) = yda
            iflux   (i,j,ng) = i2
            jflux   (i,j,ng) = j2
            fluxsign(i,j,ng) = -1
         endif

         cnd2 = xa < c0 .and. y1 > c0 .and. x1 >= c0   ! NE1 (group 1)
         if (cnd2) then
            xp1     (i,j,ng) = xca
            yp1     (i,j,ng) = yca
            xp2     (i,j,ng) = xya
            yp2     (i,j,ng) = yya
            xp3     (i,j,ng) = xda
            yp3     (i,j,ng) = yda
            iflux   (i,j,ng) = i2
            jflux   (i,j,ng) = j2
            fluxsign(i,j,ng) = 1
         endif

         ! vertices of potential group 3 triangle
         ng = 3

         cnd3 = xa > c0 .and. y1 > c0 .and. x1 < c0    ! NE2 (group 3)
         if (cnd3) then
            xp1     (i,j,ng) = xca
            yp1     (i,j,ng) = yca
            xp2     (i,j,ng) = xxa
            yp2     (i,j,ng) = yxa
            xp3     (i,j,ng) = xya
            yp3     (i,j,ng) = yya
            iflux   (i,j,ng) = i2
            jflux   (i,j,ng) = j2
            fluxsign(i,j,ng) = -1
         endif

    !-------------------------------------------------------------------
    ! contribution from triangles in N cell
    !-------------------------------------------------------------------
         i2 = i
         j2 = j+1

         xca =  p5*HTN(i2,j2-1)
         yca = -p5*HTE(i2,j2)

         xda = xca + mse(1,1,i2,j2)*x1 + mse(1,2,i2,j2)*y1
         yda = yca + mse(2,1,i2,j2)*x1 + mse(2,2,i2,j2)*y1

         xxa = xca + mse(1,1,i2,j2)*xa
         yxa = yca + mse(2,1,i2,j2)*xa

         xya = xca + mse(1,2,i2,j2)*ya
         yya = yca + mse(2,2,i2,j2)*ya

         ng = 1

         cnd1 =  xa < c0 .and. y1 > c0 .and. x1 < c0    ! N (group 1)
         if (cnd1) then
            xp1     (i,j,ng) = xca
            yp1     (i,j,ng) = yca
            xp2     (i,j,ng) = xxa
            yp2     (i,j,ng) = yxa
            xp3     (i,j,ng) = xda
            yp3     (i,j,ng) = yda
            iflux   (i,j,ng) = i2
            jflux   (i,j,ng) = j2
            fluxsign(i,j,ng) = 1
         endif

         cnd2 = xa > c0 .and. y1 > c0 .and. x1 < c0     ! N2 (group 1)
         if (cnd2) then
            xp1     (i,j,ng) = xca
            yp1     (i,j,ng) = yca
            xp2     (i,j,ng) = xya
            yp2     (i,j,ng) = yya
            xp3     (i,j,ng) = xda
            yp3     (i,j,ng) = yda
            iflux   (i,j,ng) = i2
            jflux   (i,j,ng) = j2
            fluxsign(i,j,ng) = -1
         endif

         ng = 3

         cnd3 =  xa < c0 .and. y1 > c0 .and. x1 >= c0   ! N1 (group 3)
         if (cnd3) then
            xp1     (i,j,ng) = xca
            yp1     (i,j,ng) = yca
            xp2     (i,j,ng) = xxa
            yp2     (i,j,ng) = yxa
            xp3     (i,j,ng) = xya
            yp3     (i,j,ng) = yya
            iflux   (i,j,ng) = i2
            jflux   (i,j,ng) = j2
            fluxsign(i,j,ng) = 1
         endif

    !-------------------------------------------------------------------
    ! contribution from triangles in SE cell
    !-------------------------------------------------------------------
         i2 = i+1
         j2 = j-1

         xcb = -p5*HTN(i2,j2)
         ycb =  p5*HTE(i2-1,j2)

         xdb = xcb + mnw(1,1,i2,j2)*x2 + mnw(1,2,i2,j2)*y2
         ydb = ycb + mnw(2,1,i2,j2)*x2 + mnw(2,2,i2,j2)*y2

         xxb = xcb + mnw(1,1,i2,j2)*xb
         yxb = ycb + mnw(2,1,i2,j2)*xb

         xyb = xcb + mnw(1,2,i2,j2)*yb
         yyb = ycb + mnw(2,2,i2,j2)*yb

         ng = 2

         cnd1 = xb > c0 .and. y2 < c0 .and. x2 >= c0    ! SE (group 2)
         if (cnd1) then
            xp1     (i,j,ng) = xcb
            yp1     (i,j,ng) = ycb
            xp2     (i,j,ng) = xxb
            yp2     (i,j,ng) = yxb
            xp3     (i,j,ng) = xdb
            yp3     (i,j,ng) = ydb
            iflux   (i,j,ng) = i2
            jflux   (i,j,ng) = j2
            fluxsign(i,j,ng) = -1
         endif

         cnd2 = xb < c0 .and. y2 < c0 .and. x2 >= c0    ! SE1 (group 2)
         if (cnd2) then
            xp1     (i,j,ng) = xcb
            yp1     (i,j,ng) = ycb
            xp2     (i,j,ng) = xyb
            yp2     (i,j,ng) = yyb
            xp3     (i,j,ng) = xdb
            yp3     (i,j,ng) = ydb
            iflux   (i,j,ng) = i2
            jflux   (i,j,ng) = j2
            fluxsign(i,j,ng) = 1
         endif

         ng = 3

         cnd3 = xb > c0 .and. y2 < c0 .and. x2 < c0     ! SE2 (group 3)
         if (cnd3) then
            xp1     (i,j,ng) = xcb
            yp1     (i,j,ng) = ycb
            xp2     (i,j,ng) = xxb
            yp2     (i,j,ng) = yxb
            xp3     (i,j,ng) = xyb
            yp3     (i,j,ng) = yyb
            iflux   (i,j,ng) = i2
            jflux   (i,j,ng) = j2
            fluxsign(i,j,ng) = -1
         endif

    !-------------------------------------------------------------------
    ! contribution from triangles in S cell
    !-------------------------------------------------------------------
         i2 = i
         j2 = j-1

         xcb = p5*HTN(i2,j2)
         ycb = p5*HTE(i2,j2)

         xdb = xcb + mne(1,1,i2,j2)*x2 + mne(1,2,i2,j2)*y2
         ydb = ycb + mne(2,1,i2,j2)*x2 + mne(2,2,i2,j2)*y2

         xxb = xcb + mne(1,1,i2,j2)*xb
         yxb = ycb + mne(2,1,i2,j2)*xb

         xyb = xcb + mne(1,2,i2,j2)*yb
         yyb = ycb + mne(2,2,i2,j2)*yb

         ng = 2

         cnd1 = xb < c0 .and. y2 < c0 .and. x2 < c0     ! S (group 2)
         if (cnd1) then
            xp1     (i,j,ng) = xcb
            yp1     (i,j,ng) = ycb
            xp2     (i,j,ng) = xxb
            yp2     (i,j,ng) = yxb
            xp3     (i,j,ng) = xdb
            yp3     (i,j,ng) = ydb
            iflux   (i,j,ng) = i2
            jflux   (i,j,ng) = j2
            fluxsign(i,j,ng) = 1
         endif

         cnd2 = xb > c0 .and. y2 < c0 .and. x2 < c0     ! S2 (group 2)
         if (cnd2) then
            xp1     (i,j,ng) = xcb
            yp1     (i,j,ng) = ycb
            xp2     (i,j,ng) = xyb
            yp2     (i,j,ng) = yyb
            xp3     (i,j,ng) = xdb
            yp3     (i,j,ng) = ydb
            iflux   (i,j,ng) = i2
            jflux   (i,j,ng) = j2
            fluxsign(i,j,ng) = -1
         endif

         ng = 3

         cnd3 = xb < c0 .and. y2 < c0 .and. x2 >= c0    ! S1 (group 3)
         if (cnd3) then
            xp1     (i,j,ng) = xcb
            yp1     (i,j,ng) = ycb
            xp2     (i,j,ng) = xxb
            yp2     (i,j,ng) = yxb
            xp3     (i,j,ng) = xyb
            yp3     (i,j,ng) = yyb
            iflux   (i,j,ng) = i2
            jflux   (i,j,ng) = j2
            fluxsign(i,j,ng) = 1
         endif

    !-------------------------------------------------------------------
    ! redefine departure points if not in home or east cell
    !-------------------------------------------------------------------

         if (y1 > c0) then
            x1 = xa
            y1 = c0
         endif

         if (y2 < c0) then
            x2 = xb
            y2 = c0
         endif

         ! quantity used to compute intersection point

         if (abs(xb-xa) > puny) then
            w1 = min (c1, max(c0, xb/(xb-xa)))
         else
            w1 = c0
         endif

    !-------------------------------------------------------------------
    ! contribution from triangles inside home cell
    ! Note that home and facing cells follow the same pattern.
    !-------------------------------------------------------------------

         ! load horizontal indices
         i2 = i
         j2 = j

         ! triangle vertices relative to center of home cell

         xca =  p5*HTN(i2,j2)
         yca =  p5*HTE(i2,j2)

         xcb =  p5*HTN(i2,j2-1)
         ycb = -p5*HTE(i2,j2)

         xda = xca + mne(1,1,i2,j2)*x1 + mne(1,2,i2,j2)*y1
         yda = yca + mne(2,1,i2,j2)*x1 + mne(2,2,i2,j2)*y1

         xdb = xcb + mse(1,1,i2,j2)*x2 + mse(1,2,i2,j2)*y2
         ydb = ycb + mse(2,1,i2,j2)*x2 + mse(2,2,i2,j2)*y2

         xic = p5 * (w1*(HTN(i2,j2)-HTN(i2,j2-1)) + HTN(i2,j2-1))
         yic = (w1 - p5) * HTE(i2,j2)

         ! contribution from triangle that includes the
         ! E cell edge (part of convex quadrilateral inside home cell)

         ng = 4

         cnd1 = xa*xb >= c0 .and. xa+xb < c0            ! H1a (group 4)
         if (cnd1) then
            xp1     (i,j,ng) = xca
            yp1     (i,j,ng) = yca
            xp2     (i,j,ng) = xcb
            yp2     (i,j,ng) = ycb
            xp3     (i,j,ng) = xda
            yp3     (i,j,ng) = yda
            iflux   (i,j,ng) = i2
            jflux   (i,j,ng) = j2
            fluxsign(i,j,ng) = 1
         endif

        ! contribution from triangle lying along the upper part
        ! of E edge for case of line xa-xb intersecting the edge

         cnd2 = xa*xb < c0 .and. x1 < c0                ! H1b (group 4)
         if (cnd2) then
            xp1     (i,j,ng) = xca
            yp1     (i,j,ng) = yca
            xp2     (i,j,ng) = xic
            yp2     (i,j,ng) = yic
            xp3     (i,j,ng) = xda
            yp3     (i,j,ng) = yda
            iflux   (i,j,ng) = i2
            jflux   (i,j,ng) = j2
            fluxsign(i,j,ng) = 1
         endif

         ! contribution from triangle touching but not
         ! lying along the E edge (other part of convex quadrilateral)

         ng = 5

         cnd1 = xa*xb >= c0 .and. xa+xb < c0            ! H2a (group 5)
         if (cnd1) then
            xp1     (i,j,ng) = xcb
            yp1     (i,j,ng) = ycb
            xp2     (i,j,ng) = xda
            yp2     (i,j,ng) = yda
            xp3     (i,j,ng) = xdb
            yp3     (i,j,ng) = ydb
            iflux   (i,j,ng) = i2
            jflux   (i,j,ng) = j2
            fluxsign(i,j,ng) = 1
         endif

         ! contribution from triangle lying along the lower part
         ! of E edge for case of line xa-xb intersecting the edge

         cnd2 = xa*xb < c0 .and. x2 < c0                ! H2b (group 5)
         if (cnd2) then
            xp1     (i,j,ng) = xcb
            yp1     (i,j,ng) = ycb
            xp2     (i,j,ng) = xic
            yp2     (i,j,ng) = yic
            xp3     (i,j,ng) = xdb
            yp3     (i,j,ng) = ydb
            iflux  (i,j,ng)  = i2
            jflux  (i,j,ng)  = j2
            fluxsign(i,j,ng) = 1
         endif

    !-------------------------------------------------------------------
    ! contribution from triangles in E cell
    !-------------------------------------------------------------------

         i2 = i+1
         j2 = j

         xca = -p5*HTN(i2,j2)
         yca =  p5*HTE(i2-1,j2)

         xcb = -p5*HTN(i2,j2-1)
         ycb = -p5*HTE(i2-1,j2)

         xda = xca + mnw(1,1,i2,j2)*x1 + mnw(1,2,i2,j2)*y1
         yda = yca + mnw(2,1,i2,j2)*x1 + mnw(2,2,i2,j2)*y1

         xdb = xcb + msw(1,1,i2,j2)*x2 + msw(1,2,i2,j2)*y2
         ydb = ycb + msw(2,1,i2,j2)*x2 + msw(2,2,i2,j2)*y2

         xic = -p5 * (w1*(HTN(i2,j2)-HTN(i2,j2-1)) + HTN(i2,j2-1))
         yic = (w1 - p5) * HTE(i2-1,j2)

         ! contribution from triangle that includes the
         ! W cell edge (part of convex quadrilateral inside E cell)

         ng = 4

         cnd1 = xa*xb >= c0 .and. xa+xb > c0            ! E1a (group 4)
         if (cnd1) then
            xp1     (i,j,ng) = xca
            yp1     (i,j,ng) = yca
            xp2     (i,j,ng) = xcb
            yp2     (i,j,ng) = ycb
            xp3     (i,j,ng) = xda
            yp3     (i,j,ng) = yda
            iflux   (i,j,ng) = i2
            jflux   (i,j,ng) = j2
            fluxsign(i,j,ng) = -1
         endif

         ! contribution from triangle lying along the upper part
         ! of W edge for case of line xa-xb intersecting the edge

         cnd2 = xa*xb < c0 .and. x1 > c0                ! E1b (group 4)
         if (cnd2) then
            xp1     (i,j,ng) = xca
            yp1     (i,j,ng) = yca
            xp2     (i,j,ng) = xic
            yp2     (i,j,ng) = yic
            xp3     (i,j,ng) = xda
            yp3     (i,j,ng) = yda
            iflux   (i,j,ng) = i2
            jflux   (i,j,ng) = j2
            fluxsign(i,j,ng) = -1
         endif

         ! contribution from triangle touching but not
         ! lying along the W edge (other part of convex quadrilateral)

         ng = 5

         cnd1 = xa*xb >= c0 .and. xa+xb > c0            ! E2a (group 5)
         if (cnd1) then
            xp1     (i,j,ng) = xcb
            yp1     (i,j,ng) = ycb
            xp2     (i,j,ng) = xda
            yp2     (i,j,ng) = yda
            xp3     (i,j,ng) = xdb
            yp3     (i,j,ng) = ydb
            iflux   (i,j,ng) = i2
            jflux   (i,j,ng) = j2
            fluxsign(i,j,ng) = -1
         endif

         ! contribution from triangle lying along the lower part
         ! of W edge for case of line xa-xb intersecting the edge

         cnd2 = xa*xb < c0 .and. x2 > c0                ! E2b (group 5)
         if (cnd2) then
            xp1     (i,j,ng) = xcb
            yp1     (i,j,ng) = ycb
            xp2     (i,j,ng) = xic
            yp2     (i,j,ng) = yic
            xp3     (i,j,ng) = xdb
            yp3     (i,j,ng) = ydb
            iflux   (i,j,ng) = i2
            jflux   (i,j,ng) = j2
            fluxsign(i,j,ng) = -1
         endif

         endif                  ! extmask

      enddo                     ! i
      enddo                     ! j

    !-------------------------------------------------------------------
    ! compute triangle areas with appropriate sign
    !-------------------------------------------------------------------

      do ng = 1, ngroups

         icells(ng) = 0
         do j = 1, ny_block
         do i = 1, nx_block

            w1 = p5 * abs( (xp2(i,j,ng)-xp1(i,j,ng)) *     &
                           (yp3(i,j,ng)-yp1(i,j,ng))     &
                         - (yp2(i,j,ng)-yp1(i,j,ng)) *     &
                           (xp3(i,j,ng)-xp1(i,j,ng)) )

            triarea(i,j,ng) = fluxsign(i,j,ng) * w1

            if (abs(triarea(i,j,ng)) <= puny) then
               triarea(i,j,ng) = c0 
            else
               icells(ng) = icells(ng) + 1 
               ij = icells(ng)
               indxi(ij,ng) = i
               indxj(ij,ng) = j
            endif

         enddo                  ! i
         enddo                  ! j
      enddo                     ! ng

      end subroutine locate_triangles_east

!=======================================================================
!
!BOP
!
! !IROUTINE: locate_triangles_north - triangle info for north edges
!
! !INTERFACE:
!
      subroutine locate_triangles_north(nx_block,  ny_block,   &
                                        nghost,                &
                                        extmask,   icells,     &
                                        indxi,     indxj,      &
                                        dpx,       dpy,        &
                                        HTN,       HTE,        &
                                        mne,       mnw,        &
                                        mse,       msw,        &
                                        xp0,       xp1,        &
                                        xp2,       xp3,        &
                                        yp0,       yp1,        &
                                        yp2,       yp3,        &
                                        iflux,     jflux,      &
                                        triarea)
!
! !DESCRIPTION:
!
! Compute areas and vertices of transport triangles for north cell edges.
! Note: The north edges follow the same pattern as the east edges.
! With some work, it would be possible to have a single generic
!  subroutine for both east and north edges.
!
! !REVISION HISTORY:
!
! authors William H. Lipscomb, LANL
!         John R. Baumgardner, LANL
!
! !USES:
!
! !INPUT/OUTPUT PARAMETERS:
!
      integer (kind=int_kind), intent(in) ::     &
         nx_block, ny_block ,&! block dimensions
         nghost               ! number of ghost cells

      real (kind=dbl_kind), dimension (nx_block,ny_block),     &
         intent(in) ::     &
         extmask          ! extended ice extent mask

      real (kind=dbl_kind), dimension (nx_block,ny_block),     &
           intent(in) ::     &
         dpx            ,&! x coordinates of departure points at cell corners
         dpy            ,&! y coordinates of departure points at cell corners
         HTN            ,&! length of northern edge of T-cell (m)
         HTE              ! length of southern edge of T-cell

      real (kind=dbl_kind), dimension(2,2,nx_block,ny_block),      &
            intent(in) ::     &
         mne, mnw       ,&! matrices used for coordinate transformations
         mse, msw         ! ne = northeast corner, nw = northwest, etc.

      real (kind=dbl_kind), dimension (nx_block, ny_block, ngroups),     &
           intent(out) ::     &
         triarea        ,&! area of north-edge departure triangle
         xp0,   yp0     ,&! coordinates of special triangle points
         xp1,   yp1     ,&
         xp2,   yp2     ,&
         xp3,   yp3

      integer (kind=int_kind), intent(out),     &
           dimension (nx_block, ny_block, ngroups) ::     &
         iflux          ,&! i index of cell contributing north transport
         jflux            ! j index of cell contributing north transport

      integer (kind=int_kind), dimension (ngroups), intent(out) ::     &
         icells              ! number of cells where triarea > puny

      integer (kind=int_kind), dimension (nx_block*ny_block,ngroups),     &
         intent(out) ::     &
         indxi ,&! compressed index in i-direction
         indxj   ! compressed index in j-direction
!
!EOP
!
      integer (kind=int_kind) ::     &
         i, j           ,&! horizontal indices of cell edge
         i2, j2         ,&! horizontal indices of cell contributing transport
         ng             ,&! triangle group index
         ilo,ihi,jlo,jhi,&! beginning and end of physical domain
         ij               ! compressed horizontal index

      real (kind=dbl_kind) ::     &
         x1, y1, x2, y2 ,&! x,y coordinates of departure points, as in DB
         xa, ya, xb, yb ,&! x,y coordinates of points where the lines joining
                          ! (x1,y1) and (x2,y2) cross cell edges, as in DB
         xca, yca       ,&! transformed coordinates of corner point a
         xda, yda       ,&! transformed coordinates of departure point a
         xxa, yxa       ,&! transformed coordinates of intersection point xa
         xya, yya       ,&! transformed coordinates of intersection point ya
         xcb, ycb       ,&! transformed coordinates of corner point b
         xdb, ydb       ,&! transformed coordinates of departure point b
         xxb, yxb       ,&! transformed coordinates of intersection point xb
         xyb, yyb       ,&! transformed coordinates of intersection point yb
         xic, yic       ,&! transformed coordinates of point where the
                          ! line joining dep pts intersects the face
         w1               ! work variable

      integer (kind=int_kind), dimension (nx_block,ny_block,ngroups) ::     &
         fluxsign         ! = 1 for positive flux, -1 for negative

      logical :: cnd1, cnd2, cnd3   ! conditionals

    !-------------------------------------------------------------------
    ! Triangle notation:
    ! For each edge, there are 20 triangles that can contribute,
    ! but many of these are mutually exclusive.  It turns out that
    ! at most 5 triangles can contribute to transport integrals at once.
    !
    ! For the north edge, these triangles are referred to as:
    ! (1) NW, NW1, W, W2
    ! (2) NE, NE1, E, E2
    ! (3) NW2, W1, NE2, E1
    ! (4) H1a, H1b, N1a, N1b
    ! (5) H2a, H2b, N2a, N2b
    !
    ! See Figure 3 in DB for pictures of these triangles.
    ! See Table 1 in DB for logical conditions.
    !
    ! Many triangle vertices lie at points whose coordinates are
    ! (x1,y1), (x2,y2), (xa,0), (xb,0), (0,ya), and (0,yb) in a
    ! reference frame whose origin is the cell corner.  These
    ! coordinates must be transformed to the reference frame whose
    ! origin is the geometric center of the T-cell in which the triangle
    ! is located.  The transformation is carried out using pre-computed
    ! 2x2 matrices.  There are 4 matrices (one for each corner)
    ! associated with each grid cell.  They do not describe pure
    ! rotations, because they do not preserve length.
    !-------------------------------------------------------------------

    !-------------------------------------------------------------------
    ! Initialize
    !-------------------------------------------------------------------

      ilo = 1 + nghost
      ihi = nx_block - nghost
      jlo = 1 + nghost
      jhi = ny_block - nghost

      do ng = 1, ngroups
      do j = 1, ny_block
      do i = 1, nx_block
         fluxsign(i,j,ng) = 0
         iflux   (i,j,ng) = 0
         jflux   (i,j,ng) = 0
         xp0     (i,j,ng) = c0
         xp1     (i,j,ng) = c0
         xp2     (i,j,ng) = c0
         xp3     (i,j,ng) = c0
         yp0     (i,j,ng) = c0
         yp1     (i,j,ng) = c0
         yp2     (i,j,ng) = c0
         yp3     (i,j,ng) = c0
      enddo
      enddo
      enddo

    !-------------------------------------------------------------------
    ! Main loop
    !-------------------------------------------------------------------

      do j = jlo-1, jhi  ! includes S edge of cells with index jlo
      do i = ilo, ihi

         if (extmask(i,j) > puny) then

    !-------------------------------------------------------------------
    ! coordinates of departure points
    !-------------------------------------------------------------------
         x2 = dpx(i,j)
         y2 = dpy(i,j)
         x1 = dpx(i-1,j)
         y1 = dpy(i-1,j)
         w1 =  c1 / (x1 - HTN(i,j)  - x2)
         ya = (x1*y2 - y1*(HTN(i,j) + x2)) * w1
         yb = (y2*(x1 - HTN(i,j)) - x2*y1) * w1
         if (abs(ya-yb) > puny) then
            xa = ya*HTN(i,j) / (ya - yb)
            xb = xa - HTN(i,j)
         else
            xa = c0
            xb = c0
         endif

    !-------------------------------------------------------------------
    ! contribution from triangles in NW cell
    !-------------------------------------------------------------------
         i2 = i-1
         j2 = j+1

         xca =  p5*HTN(i2,j2-1)                         ! corner pt
         yca = -p5*HTE(i2,j2)

         xda = xca + mse(1,1,i2,j2)*x1 + mse(1,2,i2,j2)*y1  ! departure pt
         yda = yca + mse(2,1,i2,j2)*x1 + mse(2,2,i2,j2)*y1

         xya = xca + mse(1,2,i2,j2)*ya                  ! ya
         yya = yca + mse(2,2,i2,j2)*ya

         xxa = xca + mse(1,1,i2,j2)*xa                  ! xa
         yxa = yca + mse(2,1,i2,j2)*xa

         ng = 1

         cnd1 = ya > c0 .and. x1 < c0 .and. y1 >= c0    ! NW (group 1)
         if (cnd1) then
            xp1     (i,j,ng) = xca
            yp1     (i,j,ng) = yca
            xp2     (i,j,ng) = xya
            yp2     (i,j,ng) = yya
            xp3     (i,j,ng) = xda
            yp3     (i,j,ng) = yda
            iflux   (i,j,ng) = i2
            jflux   (i,j,ng) = j2
            fluxsign(i,j,ng) = -1
         endif

         cnd2 = ya < c0 .and. x1 < c0 .and. y1 >= c0    ! NW1 (group 1)
         if (cnd2) then
            xp1     (i,j,ng) = xca
            yp1     (i,j,ng) = yca
            xp2     (i,j,ng) = xxa
            yp2     (i,j,ng) = yxa
            xp3     (i,j,ng) = xda
            yp3     (i,j,ng) = yda
            iflux   (i,j,ng) = i2
            jflux   (i,j,ng) = j2
            fluxsign(i,j,ng) = 1
         endif

         ng = 3

         cnd3 = ya > c0 .and. x1 < c0 .and. y1 < c0     ! NW2 (group 3)
         if (cnd3) then
            xp1     (i,j,ng) = xca
            yp1     (i,j,ng) = yca
            xp2     (i,j,ng) = xya
            yp2     (i,j,ng) = yya
            xp3     (i,j,ng) = xxa
            yp3     (i,j,ng) = yxa
            iflux   (i,j,ng) = i2
            jflux   (i,j,ng) = j2
            fluxsign(i,j,ng) = -1
         endif

    !-------------------------------------------------------------------
    ! contribution from triangles in W cell
    !-------------------------------------------------------------------
         i2 = i-1
         j2 = j

         xca = p5*HTN(i2,j2)
         yca = p5*HTE(i2,j2)

         xda = xca + mne(1,1,i2,j2)*x1 + mne(1,2,i2,j2)*y1
         yda = yca + mne(2,1,i2,j2)*x1 + mne(2,2,i2,j2)*y1

         xya = xca + mne(1,2,i2,j2)*ya
         yya = yca + mne(2,2,i2,j2)*ya

         xxa = xca + mne(1,1,i2,j2)*xa
         yxa = yca + mne(2,1,i2,j2)*xa

         ng = 1

         cnd1 = ya < c0 .and. x1 < c0 .and. y1 < c0     ! W (group 1)
         if (cnd1) then
            xp1     (i,j,ng) = xca
            yp1     (i,j,ng) = yca
            xp2     (i,j,ng) = xya
            yp2     (i,j,ng) = yya
            xp3     (i,j,ng) = xda
            yp3     (i,j,ng) = yda
            iflux   (i,j,ng) = i2
            jflux   (i,j,ng) = j2
            fluxsign(i,j,ng) = 1
         endif

         cnd2 = ya > c0 .and. x1 < c0 .and. y1 < c0     ! W2 (group 1)
         if (cnd2) then
            xp1     (i,j,ng) = xca
            yp1     (i,j,ng) = yca
            xp2     (i,j,ng) = xxa
            yp2     (i,j,ng) = yxa
            xp3     (i,j,ng) = xda
            yp3     (i,j,ng) = yda
            iflux   (i,j,ng) = i2
            jflux   (i,j,ng) = j2
            fluxsign(i,j,ng) = -1
         endif

         ng = 3

         cnd3 = ya < c0 .and. x1 < c0 .and. y1 >= c0    ! W1 (group 3)
         if (cnd3) then
            xp1     (i,j,ng) = xca
            yp1     (i,j,ng) = yca
            xp2     (i,j,ng) = xya
            yp2     (i,j,ng) = yya
            xp3     (i,j,ng) = xxa
            yp3     (i,j,ng) = yxa
            iflux   (i,j,ng) = i2
            jflux   (i,j,ng) = j2
            fluxsign(i,j,ng) = 1
         endif

    !-------------------------------------------------------------------
    ! contribution from triangles in NE cell
    !-------------------------------------------------------------------
         i2 = i+1
         j2 = j+1

         xcb = -p5*HTN(i2,j2-1)
         ycb = -p5*HTE(i2-1,j2)

         xdb = xcb + msw(1,1,i2,j2)*x2 + msw(1,2,i2,j2)*y2
         ydb = ycb + msw(2,1,i2,j2)*x2 + msw(2,2,i2,j2)*y2

         xyb = xcb + msw(1,2,i2,j2)*yb
         yyb = ycb + msw(2,2,i2,j2)*yb

         xxb = xcb + msw(1,1,i2,j2)*xb
         yxb = ycb + msw(2,1,i2,j2)*xb

         ng = 2

         cnd1 = yb > c0 .and. x2 > c0 .and. y2 >= c0    ! NE (group 2)
         if (cnd1) then
            xp1     (i,j,ng) = xcb
            yp1     (i,j,ng) = ycb
            xp2     (i,j,ng) = xyb
            yp2     (i,j,ng) = yyb
            xp3     (i,j,ng) = xdb
            yp3     (i,j,ng) = ydb
            iflux   (i,j,ng) = i2
            jflux   (i,j,ng) = j2
            fluxsign(i,j,ng) = -1
         endif

         cnd2 = yb < c0 .and. x2 > c0  .and. y2 >= c0   ! NE1 (group 2)
         if (cnd2) then
            xp1     (i,j,ng) = xcb
            yp1     (i,j,ng) = ycb
            xp2     (i,j,ng) = xxb
            yp2     (i,j,ng) = yxb
            xp3     (i,j,ng) = xdb
            yp3     (i,j,ng) = ydb
            iflux   (i,j,ng) = i2
            jflux   (i,j,ng) = j2
            fluxsign(i,j,ng) = 1
         endif

         ng = 3

         cnd3 = yb > c0 .and. x2 > c0 .and. y2 < c0     ! NE2 (group 3)
         if (cnd3) then
            xp1     (i,j,ng) = xcb
            yp1     (i,j,ng) = ycb
            xp2     (i,j,ng) = xyb
            yp2     (i,j,ng) = yyb
            xp3     (i,j,ng) = xxb
            yp3     (i,j,ng) = yxb
            iflux   (i,j,ng) = i2
            jflux   (i,j,ng) = j2
            fluxsign(i,j,ng) = -1
         endif

    !-------------------------------------------------------------------
    ! contribution from triangles in E cell
    !-------------------------------------------------------------------
         i2 = i+1
         j2 = j

         xcb = -p5*HTN(i2,j2)
         ycb =  p5*HTE(i2-1,j2)

         xdb = xcb + mnw(1,1,i2,j2)*x2 + mnw(1,2,i2,j2)*y2
         ydb = ycb + mnw(2,1,i2,j2)*x2 + mnw(2,2,i2,j2)*y2

         xyb = xcb + mnw(1,2,i2,j2)*yb
         yyb = ycb + mnw(2,2,i2,j2)*yb

         xxb = xcb + mnw(1,1,i2,j2)*xb
         yxb = ycb + mnw(2,1,i2,j2)*xb

         ng = 2

         cnd1 = yb < c0 .and. x2 > c0 .and. y2 < c0     ! E (group 2)
         if (cnd1) then
            xp1     (i,j,ng) = xcb
            yp1     (i,j,ng) = ycb
            xp2     (i,j,ng) = xyb
            yp2     (i,j,ng) = yyb
            xp3     (i,j,ng) = xdb
            yp3     (i,j,ng) = ydb
            iflux   (i,j,ng) = i2
            jflux   (i,j,ng) = j2
            fluxsign(i,j,ng) = 1
         endif

         cnd2 = yb > c0 .and. x2 > c0 .and. y2 < c0     ! E2 (group 2)
         if (cnd2) then
            xp1     (i,j,ng) = xcb
            yp1     (i,j,ng) = ycb
            xp2     (i,j,ng) = xxb
            yp2     (i,j,ng) = yxb
            xp3     (i,j,ng) = xdb
            yp3     (i,j,ng) = ydb
            iflux   (i,j,ng) = i2
            jflux   (i,j,ng) = j2
            fluxsign(i,j,ng) = -1
         endif

         ng = 3

         cnd3 = yb < c0 .and. x2 > c0 .and. y2 >= c0    ! E1 (group 3)
         if (cnd3) then
            xp1     (i,j,ng) = xcb
            yp1     (i,j,ng) = ycb
            xp2     (i,j,ng) = xyb
            yp2     (i,j,ng) = yyb
            xp3     (i,j,ng) = xxb
            yp3     (i,j,ng) = yxb
            iflux   (i,j,ng) = i2
            jflux   (i,j,ng) = j2
            fluxsign(i,j,ng) = 1
         endif

    !-------------------------------------------------------------------
    ! redefine departure points if not in home or north cell
    !-------------------------------------------------------------------
         if (x1 < c0) then
            x1 = c0
            y1 = ya
         endif

         if (x2 > c0) then
            x2 = c0
            y2 = yb
         endif

         ! quantity used to compute intersection point

         if (abs(yb-ya) > puny) then
            w1 = min (c1, max(c0, yb/(yb-ya)))
         else
            w1 = c0
         endif

    !-------------------------------------------------------------------
    ! contribution from triangles inside home cell
    !-------------------------------------------------------------------
         i2 = i
         j2 = j

         xca = -p5*HTN(i2,j2)
         yca =  p5*HTE(i2-1,j2)

         xcb =  p5*HTN(i2,j2)
         ycb =  p5*HTE(i2,j2)

         xda = xca + mnw(1,1,i2,j2)*x1 + mnw(1,2,i2,j2)*y1
         yda = yca + mnw(2,1,i2,j2)*x1 + mnw(2,2,i2,j2)*y1

         xdb = xcb + mne(1,1,i2,j2)*x2 + mne(1,2,i2,j2)*y2
         ydb = ycb + mne(2,1,i2,j2)*x2 + mne(2,2,i2,j2)*y2

         xic = (p5 - w1) * HTN(i2,j2)          ! intersection w/ N edge
         yic = p5 * (w1*(HTE(i2-1,j2)-HTE(i2,j2)) + HTE(i2,j2))

         ! contribution from triangle that includes the
         ! N cell edge (part of convex quadrilateral inside home cell)

         ng = 4

         cnd1 = ya*yb >= c0 .and. ya+yb < c0            ! H1a (group 4)
         if (cnd1) then
            xp1     (i,j,ng) = xca
            yp1     (i,j,ng) = yca
            xp2     (i,j,ng) = xcb
            yp2     (i,j,ng) = ycb
            xp3     (i,j,ng) = xda
            yp3     (i,j,ng) = yda
            iflux   (i,j,ng) = i2
            jflux   (i,j,ng) = j2
            fluxsign(i,j,ng) = 1
         endif

         ! contribution from triangle lying along the left part of
         ! the N edge for case of line ya-yb intersecting the edge

         cnd2 = ya*yb < c0 .and. y1 < c0                ! H1b (group 4)
         if (cnd2) then
            xp1     (i,j,ng) = xca
            yp1     (i,j,ng) = yca
            xp2     (i,j,ng) = xic
            yp2     (i,j,ng) = yic
            xp3     (i,j,ng) = xda
            yp3     (i,j,ng) = yda
            iflux   (i,j,ng) = i2
            jflux   (i,j,ng) = j2
            fluxsign(i,j,ng) = 1
         endif

         ! contribution from triangle touching but not lying
         ! along the N edge (other part of convex quadrilateral)

         ng = 5

         cnd1 = ya*yb >= c0 .and. ya+yb < c0            ! H2a (group 5)
         if (cnd1) then
            xp1     (i,j,ng) = xcb
            yp1     (i,j,ng) = ycb
            xp2     (i,j,ng) = xda
            yp2     (i,j,ng) = yda
            xp3     (i,j,ng) = xdb
            yp3     (i,j,ng) = ydb
            iflux   (i,j,ng) = i2
            jflux   (i,j,ng) = j2
            fluxsign(i,j,ng) = 1
         endif

         ! contribution from triangle lying along the right part
         ! of the N edge for case of line ya-yb intersecting the edge

         cnd2 = ya*yb < c0 .and. y2 < c0                ! H2b (group 5)
         if (cnd2) then
            xp1     (i,j,ng) = xcb
            yp1     (i,j,ng) = ycb
            xp2     (i,j,ng) = xic
            yp2     (i,j,ng) = yic
            xp3     (i,j,ng) = xdb
            yp3     (i,j,ng) = ydb
            iflux   (i,j,ng) = i2
            jflux   (i,j,ng) = j2
            fluxsign(i,j,ng) = 1
         endif

    !-------------------------------------------------------------------
    ! contribution from triangles in N cell
    !-------------------------------------------------------------------
         i2 = i
         j2 = j+1

         xca = -p5*HTN(i2,j2-1)
         yca = -p5*HTE(i2-1,j2)

         xcb =  p5*HTN(i2,j2-1)
         ycb = -p5*HTE(i2,j2)

         xda = xca + msw(1,1,i2,j2)*x1 + msw(1,2,i2,j2)*y1
         yda = yca + msw(2,1,i2,j2)*x1 + msw(2,2,i2,j2)*y1

         xdb = xcb + mse(1,1,i2,j2)*x2 + mse(1,2,i2,j2)*y2
         ydb = ycb + mse(2,1,i2,j2)*x2 + mse(2,2,i2,j2)*y2

         xic = (p5 - w1)*HTN(i2,j2-1)
         yic = -p5 * (w1*(HTE(i2-1,j2)-HTE(i2,j2)) + HTE(i2,j2))

         ! contribution from triangle that includes the
         ! S cell edge (part of convex quadrilateral inside home cell)

         ng = 4

         cnd1 = ya*yb >= c0 .and. ya+yb > c0            ! N1a (group 4)
         if (cnd1) then
            xp1     (i,j,ng) = xca
            yp1     (i,j,ng) = yca
            xp2     (i,j,ng) = xcb
            yp2     (i,j,ng) = ycb
            xp3     (i,j,ng) = xda
            yp3     (i,j,ng) = yda
            iflux   (i,j,ng) = i2
            jflux   (i,j,ng) = j2
            fluxsign(i,j,ng) = -1
         endif

         ! contribution from triangle lying along the left part
         ! of the S edge for case of line ya-yb intersecting the edge

         cnd2 = ya*yb < c0 .and. y1 > c0                ! N1b (group 4)
         if (cnd2) then
            xp1     (i,j,ng) = xca
            yp1     (i,j,ng) = yca
            xp2     (i,j,ng) = xic
            yp2     (i,j,ng) = yic
            xp3     (i,j,ng) = xda
            yp3     (i,j,ng) = yda
            iflux   (i,j,ng) = i2
            jflux   (i,j,ng) = j2
            fluxsign(i,j,ng) = -1
         endif

         ! contribution from triangle touching but not
         ! lying along the S edge (other part of convex quadrilateral)

         ng = 5

         cnd1 = ya*yb >= c0 .and. ya+yb > c0            ! N2a (group 5)
         if (cnd1) then
            xp1     (i,j,ng) = xcb
            yp1     (i,j,ng) = ycb
            xp2     (i,j,ng) = xda
            yp2     (i,j,ng) = yda
            xp3     (i,j,ng) = xdb
            yp3     (i,j,ng) = ydb
            iflux   (i,j,ng) = i2
            jflux   (i,j,ng) = j2
            fluxsign(i,j,ng) = -1
         endif

         ! contribution from triangle lying along the right part
         ! of the S edge for case of line ya-yb intersecting the edge

         cnd2 = ya*yb < c0 .and. y2 > c0                ! N2b (group 5)
         if (cnd2) then
            xp1     (i,j,ng) = xcb
            yp1     (i,j,ng) = ycb
            xp2     (i,j,ng) = xic
            yp2     (i,j,ng) = yic
            xp3     (i,j,ng) = xdb
            yp3     (i,j,ng) = ydb
            iflux   (i,j,ng) = i2
            jflux   (i,j,ng) = j2
            fluxsign(i,j,ng) = -1
         endif

         endif                  ! extmask

      enddo                     ! i loop
      enddo                     ! j loop

    !-------------------------------------------------------------------
    ! compute triangle areas with appropriate sign
    !-------------------------------------------------------------------

      do ng = 1, ngroups

         icells(ng) = 0
         do j = 1, ny_block
         do i = 1, nx_block

            w1 = p5 * abs( (xp2(i,j,ng)-xp1(i,j,ng)) *     &
                           (yp3(i,j,ng)-yp1(i,j,ng))     &
                         - (yp2(i,j,ng)-yp1(i,j,ng)) *     &
                           (xp3(i,j,ng)-xp1(i,j,ng)) )

            triarea(i,j,ng) = fluxsign(i,j,ng) * w1

            if (abs(triarea(i,j,ng)) <= puny) then
               triarea(i,j,ng) = c0 
            else
               icells(ng) = icells(ng) + 1 
               ij = icells(ng)
               indxi(ij,ng) = i
               indxj(ij,ng) = j
            endif

         enddo                  ! i
         enddo                  ! j
      enddo                     ! ng

      end subroutine locate_triangles_north

!=======================================================================
!
!BOP
!
! !IROUTINE: triangle_coordinates - find coordinates of quadrature points
!
! !INTERFACE:
!
      subroutine triangle_coordinates (nx_block, ny_block,  &
                                       integral_order,      &
                                       icells,              &
                                       indxi,  indxj,       &
                                       xp0,      xp1,       &
                                       xp2,      xp3,       &
                                       yp0,      yp1,       &
                                       yp2,      yp3)
!
! !DESCRIPTION:
!
! For each triangle, find the coordinates of the quadrature points needed
!  to compute integrals of linear, quadratic, or cubic polynomials,
!  using formulas from A.H. Stroud, Approximate Calculation of Multiple
!  Integrals, Prentice-Hall, 1971.  (Section 8.8, formula 3.1.)
! Linear functions can be integrated exactly by evaluating the function 
!  at just one point (the midpoint).  Quadratic functions require
!  3 points, and cubics require 4 points.
! The default is cubic, but the code can be sped up slightly using 
!  linear or quadratic integrals, usually with little loss of accuracy.
!
! The formulas are as follows:
!
! I1 = integral of f(x,y)*dA
!    = A * f(x0,y0)
! where A is the traingle area and (x0,y0) is the midpoint.
!
! I2 = A * (f(x1,y1) + f(x2,y2) + f(x3,y3))
! where these three points are located halfway between the midpoint
! and the three vertics of the triangle.
!
! I3 = A * [ -9/16 *  f(x0,y0)
!           + 25/48 * (f(x1,y1) + f(x2,y2) + f(x3,y3))]
! where (x0,y0) is the midpoint, and the other three points are
! located 2/5 of the way from the midpoint to the three vertices.
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
           nx_block, ny_block ,&! block dimensions
           integral_order       ! polynomial order for quadrature integrals 

      integer (kind=int_kind), dimension (ngroups), intent(in) ::     &
           icells              ! number of cells where triarea > puny

      integer (kind=int_kind), dimension (nx_block*ny_block,ngroups),     &
           intent(in) ::     &
           indxi ,&! compressed index in i-direction
           indxj   ! compressed index in j-direction

      real (kind=dbl_kind), intent(inout),     &
           dimension (nx_block, ny_block, ngroups) ::     &
           xp0, yp0          ,&! coordinates of triangle points
           xp1, yp1          ,&
           xp2, yp2          ,&
           xp3, yp3
!
!EOP
!
      integer (kind=int_kind) ::     &
           i, j              ,&! horizontal indices
           ng                ,&! triangle index
           ij                  ! compressed horizontal index


      if (integral_order == 1) then ! linear (1-point formula)

         do ng = 1, ngroups
         do ij = 1, icells(ng)
            i = indxi(ij,ng)
            j = indxj(ij,ng)

            ! coordinates of midpoint
               xp0(i,j,ng) = p333     &
                           * (xp1(i,j,ng) + xp2(i,j,ng) + xp3(i,j,ng))
               yp0(i,j,ng) = p333     &
                           * (yp1(i,j,ng) + yp2(i,j,ng) + yp3(i,j,ng))

         enddo                  ! ij
         enddo                  ! ng

      elseif (integral_order == 2) then ! quadratic (3-point formula)

!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
         do ng = 1, ngroups
         do ij = 1, icells(ng)
            i = indxi(ij,ng)
            j = indxj(ij,ng)

            ! coordinates of midpoint
               xp0(i,j,ng) = p333     &
                           * (xp1(i,j,ng) + xp2(i,j,ng) + xp3(i,j,ng))
               yp0(i,j,ng) = p333     &
                           * (yp1(i,j,ng) + yp2(i,j,ng) + yp3(i,j,ng))

            ! coordinates of the 3 points needed for integrals

               xp1(i,j,ng) = p5*xp1(i,j,ng) + p5*xp0(i,j,ng)
               yp1(i,j,ng) = p5*yp1(i,j,ng) + p5*yp0(i,j,ng)

               xp2(i,j,ng) = p5*xp2(i,j,ng) + p5*xp0(i,j,ng)
               yp2(i,j,ng) = p5*yp2(i,j,ng) + p5*yp0(i,j,ng)

               xp3(i,j,ng) = p5*xp3(i,j,ng) + p5*xp0(i,j,ng)
               yp3(i,j,ng) = p5*yp3(i,j,ng) + p5*yp0(i,j,ng)

         enddo                  ! ij
         enddo                  ! ng

      else                      ! cubic (4-point formula)

!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
         do ng = 1, ngroups
         do ij = 1, icells(ng)
            i = indxi(ij,ng)
            j = indxj(ij,ng)

            ! coordinates of midpoint
               xp0(i,j,ng) = p333     &
                           * (xp1(i,j,ng) + xp2(i,j,ng) + xp3(i,j,ng))
               yp0(i,j,ng) = p333     &
                           * (yp1(i,j,ng) + yp2(i,j,ng) + yp3(i,j,ng))

            ! coordinates of the other 3 points needed for integrals

               xp1(i,j,ng) = p4*xp1(i,j,ng) + p6*xp0(i,j,ng)
               yp1(i,j,ng) = p4*yp1(i,j,ng) + p6*yp0(i,j,ng)

               xp2(i,j,ng) = p4*xp2(i,j,ng) + p6*xp0(i,j,ng)
               yp2(i,j,ng) = p4*yp2(i,j,ng) + p6*yp0(i,j,ng)

               xp3(i,j,ng) = p4*xp3(i,j,ng) + p6*xp0(i,j,ng)
               yp3(i,j,ng) = p4*yp3(i,j,ng) + p6*yp0(i,j,ng)

         enddo                  ! ij
         enddo                  ! ng

      endif

      end subroutine triangle_coordinates

!=======================================================================
!
!BOP
!
! !IROUTINE: construct_fields - construct fields of ice area and tracers
!
! !INTERFACE:
!
      subroutine construct_fields (nx_block, ny_block,   &
                                   nghost,               &
                                   icells,               &
                                   indxi,   indxj,       &
                                   HTN,      HTE,        &
                                   hm,       xav,        &
                                   yav,      xxav,       &
                                   xyav,     yyav,       &
                                   xxxav,    xxyav,      &
                                   xyyav,    yyyav,      &
                                   dxt,      dyt,        &
                                   aim,      aic,        &
                                   aix,      aiy,        &
                                   aimask,               &
                                   trm,      trc,        &
                                   trx,      try,        &
                                   trmask)
!
! !DESCRIPTION:
!
! Construct fields of ice area and tracers.
!
! !REVISION HISTORY:
!
! authors William H. Lipscomb, LANL
!         John R. Baumgardner, LANL
!
! !USES:
!
! !INPUT/OUTPUT PARAMETERS:
!
      integer (kind=int_kind), intent(in) ::     &
         nx_block, ny_block ,&! block dimensions
         nghost              ! number of ghost cells

      integer (kind=int_kind), dimension(nx_block*ny_block),     &
         intent(in) ::     &
         indxi        ,&! compressed i/j indices
         indxj

      integer (kind=int_kind), intent(in) ::     &
         icells         ! number of cells with ice

      real (kind=dbl_kind), dimension (nx_block,ny_block),     &
         intent(in) ::     &
         hm             ,&! land/boundary mask, thickness (T-cell)
         HTN            ,&! length of northern edge of T-cell (m)
         HTE            ,&! length of eastern edge of T-cell (m)
         xav,  yav               ,&! mean T-cell values of x, y
         xxav, xyav, yyav        ,&! mean values of xx, xy, yy
         xxxav,xxyav,xyyav,yyyav ,&! mean values of xxx, xxy, xyy, yyy
         dxt            ,&! grid cell width (m)
         dyt              ! grid cell height (m)

      real (kind=dbl_kind), dimension (nx_block,ny_block),     &
         intent(in) ::     &
         aim            ,&! mean ice area
         aimask           ! = 1. if ice is present, = 0. otherwise

      real (kind=dbl_kind), dimension (nx_block,ny_block,ntrace),     &
         intent(in), optional ::     &
         trm            ,&! mean tracer
         trmask           ! = 1. if tracer is present, = 0. otherwise

      real (kind=dbl_kind), dimension (nx_block,ny_block),     &
         intent(out) ::     &
         aic            ,&! ice area at geometric center of cell
         aix, aiy         ! limited derivative of ice area wrt x and y

      real (kind=dbl_kind), dimension (nx_block,ny_block,ntrace),     &
         intent(out), optional ::     &
         trc            ,&! tracer at geometric center of cell
         trx, try         ! limited derivative of tracer wrt x and y
!
!EOP
!
      integer (kind=int_kind) ::     &
         i, j           ,&! horizontal indices
         nt, nt1        ,&! tracer indices
         ilo,ihi,jlo,jhi,&! beginning and end of physical domain
         ij, m          ,&! combined i/j horizontal indices
         ilim

      integer (kind=int_kind), dimension(nx_block*ny_block) ::     &
         indxii, indxjj ,&! combined i/j horizontal indices
         indxij

      real (kind=dbl_kind), dimension (icells) ::     &
         axav, ayav       ! x,y coordinates of center of ice area

      real (kind=dbl_kind), dimension (icells,ntrace) ::     &
         atxav, atyav     ! x,y coordinates of center of area*tracer

      real (kind=dbl_kind), dimension(nx_block*ny_block) ::     &
         work1, work2     ! work arrays

      real (kind=dbl_kind) ::     &
         w1, w2, w3, w4, w5, w6, w7   ! work variables

    !-------------------------------------------------------------------
    ! Compute field values at the geometric center of each grid cell,
    ! and compute limited gradients in the x and y directions.
    !
    ! For second order accuracy, each state variable is approximated as
    ! a field varying linearly over x and y within each cell.  For each
    ! category, the integrated value of aicen(x,y) over the cell must
    ! equal aicen(i,j,n)*tarea(i,j), where tarea(i,j) is the cell area.
    ! Similarly, the integrated value of aicen(x,y)*hicen(x,y) must equal
    ! the total ice volume, aicen(i,j,n)*hicen(i,j,n)*tarea(i,j).
    ! And for a given layer, the integrated value of
    ! aicen(x,y)*hice(x,y)*qice(x,y) must equal the total ice layer
    ! energy, aicen(i,j,n)*hice(i,j)*qice(i,j,k)*tarea(i,j).
    !
    ! These integral conditions are satisfied for linear fields if we
    ! stipulate the following:
    ! (1) The mean ice area, aicen(i,j,n), is equal to the area at
    ! the cell centroid: the point where an equal grid cell area
    ! (not ice area!) lies to the north and south, and to the east
    ! and west.
    ! (2) The mean ice thickness, hice(i,j), is equal to the
    ! thickness at the center of ice area: the point where an equal
    ! ice area lies to the north and south, and to the east and west.
    ! (And similarly for hsno(i,j) and trcrn(i,j,n).)
    ! (3) The mean ice enthalpy, qice(i,j,k) is equal to the enthalpy
    ! at the center of ice volume: the point where an equal
    ! ice volume lies to the north and south, and to the east and west.
    ! (And similarly for qsno(i,j), which is the snow enthalpy
    ! at the center of snow volume.)
    !
    ! We want to find the value of each state variable at a standard
    ! reference point, which we choose to be the geometric center of
    ! the cell.  The geometric center is located at the intersection
    ! of the line joining the midpoints of the north and south edges
    ! with the line joining the midpoints of the east and west edges.
    ! To find the value at the geometric center, we must know the
    ! location of the cell centroid/center of ice area/center of ice
    ! or snow volume relative to the geometric center, along with
    ! the field gradients with respect to x and y.
    !
    ! The cell gradients are first computed from the difference between
    ! values in the neighboring cells, then limited by requiring that
    ! no new extrema are created within the cell.
    !-------------------------------------------------------------------

      ilo = 1 + nghost
      ihi = nx_block - nghost
      jlo = 1 + nghost
      jhi = ny_block - nghost

      ilim = 0
!      do j = 1, ny_block
!      do i = 1, nx_block
      do j = jlo, jhi
      do i = ilo, ihi
         ilim = ilim + 1
         indxii(ilim) = i
         indxjj(ilim) = j
         indxij(ilim) = ilim

         work1(ilim) = xav(i,j)
         work2(ilim) = yav(i,j)

         aic(i,j) = c0
      enddo
      enddo

      call limited_gradient (nx_block, ny_block,          &
                             ilim,     ilim,              &
                             indxii,   indxjj,   indxij,  &
                             aim,      hm,                &
                             work1,    work2,             &
                             HTN,      HTE,               &
                             dxt,      dyt,               &
                             aix,      aiy)

      ! ice area

      if (icells > 0) then

      do ij = 1,icells   ! ice is present
         i = indxi(ij)
         j = indxj(ij)

         ! ice area at geometric center
         aic(i,j) = aim(i,j) - xav(i,j)*aix(i,j)     &
                             - yav(i,j)*aiy(i,j)
         ! center of ice area (axav,ayav) for each cell
         axav(ij) = (aix(i,j)*xxav(i,j)      &
                    + aiy(i,j)*xyav(i,j)     &
                    + aic(i,j)*xav (i,j)) / aim(i,j)
         ayav(ij) = (aix(i,j)*xyav(i,j)      &
                    + aiy(i,j)*yyav(i,j)     &
                    + aic(i,j)*yav(i,j)) / aim(i,j)
      enddo                     ! ij

      ! tracers

      if (present(trm)) then

       do nt = 1, ntrace

         trc  (:,:,nt) = c0

         if (tracer_type(nt)==1) then ! independent of other tracers
                                 ! (surface temp, ice and snow thickness)

            do ij = 1, icells
               indxij(ij) = ij
            enddo

            call limited_gradient(nx_block,     ny_block,   &
                                  icells,       icells,     &
                                  indxi, indxj, indxij,     &
                                  trm(:,:,nt),  aimask,     &
                                  axav,         ayav,       &
                                  HTN,          HTE,        &
                                  dxt,          dyt,        &
                                  trx(:,:,nt),  try(:,:,nt)) 

            if (has_dependents(nt)) then   ! need center of area*tracer

!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
               do ij = 1, icells  ! Note: no trx or try in ghost cells
                                  ! (bound calls are later)
                  i = indxi(ij)
                  j = indxj(ij)

                  ! tracer value at geometric center
                  trc(i,j,nt) = trm(i,j,nt) - trx(i,j,nt)*axav(ij)     &
                                            - try(i,j,nt)*ayav(ij)

                  if (trmask(i,j,nt) > puny) then

                     ! center of area*tracer
                     w1 = aic(i,j)*trc(i,j,nt)
                     w2 = aic(i,j)*trx(i,j,nt)     &
                        + aix(i,j)*trc(i,j,nt)
                     w3 = aic(i,j)*try(i,j,nt)     &
                        + aiy(i,j)*trc(i,j,nt)
                     w4 = aix(i,j)*trx(i,j,nt)
                     w5 = aix(i,j)*try(i,j,nt)     &
                        + aiy(i,j)*trx(i,j,nt)
                     w6 = aiy(i,j)*try(i,j,nt)
                     w7 = c1 / (aim(i,j)*trm(i,j,nt))
                     atxav(ij,nt) = (w1*xav (i,j)  + w2*xxav (i,j)       &
                                    + w3*xyav (i,j) + w4*xxxav(i,j)      &
                                    + w5*xxyav(i,j) + w6*xyyav(i,j))     &
                                    * w7
                     atyav(ij,nt) = (w1*yav(i,j)   + w2*xyav (i,j)       &
                                    + w3*yyav(i,j)  + w4*xxyav(i,j)      &
                                    + w5*xyyav(i,j) + w6*yyyav(i,j))     &
                                    * w7
                  else
                     atxav(ij,nt) = c0
                     atyav(ij,nt) = c0
                  endif         ! trmask

               enddo            ! ij

            else                ! no dependents

               do ij = 1, icells      ! ice is present
                  i = indxi(ij)
                  j = indxj(ij)

                  ! tracer value at geometric center
                  trc(i,j,nt) = trm(i,j,nt) - trx(i,j,nt)*axav(ij)     &
                                            - try(i,j,nt)*ayav(ij)
               enddo            ! ij

            endif               ! has_dependents

         elseif (tracer_type(nt)==2) then   ! tracer nt depends on nt1
                                             ! (ice and snow enthalpy)
            nt1 = depend(nt)

            ilim = 0
            do ij = 1, icells
               i =  indxi(ij)
               j =  indxj(ij)
               if (trmask(i,j,nt1) > puny) then
                  ilim = ilim + 1
                  indxii(ilim) = i
                  indxjj(ilim) = j
                  indxij(ilim) = ij
               endif                  ! phimask > puny
            enddo

            call limited_gradient(nx_block,     ny_block,         &
                                  icells,       ilim,             &
                                  indxii, indxjj, indxij,         &
                                  trm(:,:,nt),  trmask(:,:,nt1),  &
                                  atxav(:,nt1), atyav(:,nt1),     &
                                  HTN,          HTE,              &
                                  dxt,          dyt,              &
                                  trx(:,:,nt),  try(:,:,nt))    

               trc(:,:,nt) = c0
               do ij = 1, ilim      ! ice is present
                  i = indxii(ij)
                  j = indxjj(ij)
                  m = indxij(ij)
                  trc(i,j,nt) = trm(i,j,nt)     &
                              - trx(i,j,nt) * atxav(m,nt1)     &
                              - try(i,j,nt) * atyav(m,nt1)
            enddo               ! ij

         elseif (tracer_type(nt) ==3) then  ! upwind approx; gradient = 0

            do j = jlo, jhi
            do i = ilo, ihi
               trc(i,j,nt) = trm(i,j,nt)
               trx(i,j,nt) = c0
               try(i,j,nt) = c0
            enddo               ! i
            enddo               ! j

         endif                  ! tracer_type
       enddo                    ! ntrace

      endif                     ! present (trm)

      endif                     ! icells > 0

      end subroutine construct_fields

!=======================================================================
!
!BOP
!
! !IROUTINE: transport_integrals - compute transports across each edge
!
! !INTERFACE:
!
      subroutine transport_integrals (nx_block, ny_block,   &
                                      integral_order,       &
                                      triarea,  icells,     &
                                      indxi,    indxj,      &
                                      iflux,    jflux,      &
                                      xp0,      xp1,        &
                                      xp2,      xp3,        &
                                      yp0,      yp1,        &
                                      yp2,      yp3,        &
                                      aic,      aix,        &
                                      aiy,      aiflx,      &
                                      trc,      trx,        &
                                      try,      atflx)
!
! !DESCRIPTION:
!
! Compute the transports across each edge by integrating the ice area
! and tracers over each departure triangle.
! Input variables have the same meanings as in the main subroutine.
! Repeated use of certain sums makes the calculation more efficient.
! Integral formulas are described in triangle_coordinates subroutine.
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
           nx_block, ny_block ,&! block dimensions
           integral_order   ! polynomial order for quadrature integrals 

      integer (kind=int_kind), dimension (ngroups), intent(in) ::     &
           icells              ! number of cells where triarea > puny

      integer (kind=int_kind), dimension (nx_block*ny_block,ngroups),     &
           intent(in) ::     &
           indxi ,&! compressed index in i-direction
           indxj   ! compressed index in j-direction

      real (kind=dbl_kind), intent(in),     &
           dimension (nx_block, ny_block, ngroups) ::     &
           triarea     ,&
           xp0, yp0    ,&
           xp1, yp1    ,&
           xp2, yp2    ,&
           xp3, yp3

      integer (kind=int_kind), intent(in),     &
           dimension (nx_block, ny_block, ngroups) ::     &
           iflux  ,&
           jflux

      real (kind=dbl_kind), intent(in),     &
           dimension (nx_block, ny_block) ::     &
           aic, aix, aiy

      real (kind=dbl_kind), intent(out),     &
           dimension (nx_block, ny_block) ::     &
           aiflx

      real (kind=dbl_kind), intent(in),     &
           dimension (nx_block, ny_block, ntrace), optional ::     &
           trc, trx, try

      real (kind=dbl_kind), intent(out),     &
           dimension (nx_block, ny_block, ntrace), optional ::     &
           atflx
!
!EOP
!
      integer (kind=int_kind) ::     &
           i, j          ,&! horizontal indices of edge
           i2, j2        ,&! horizontal indices of cell contributing transport
           ng            ,&! triangle index
           nt, nt1       ,&! tracer indices
           ij              ! combined i/j index

      real (kind=dbl_kind) ::     &
           a0, a1, a2, a3         ,&! ice area at internal points
           w0, w1, w2, w3           ! work variables

      real (kind=dbl_kind), dimension (nx_block, ny_block) ::     &
           asum, axsum, aysum     ,&! sum of area, area*x, and area*y
           axxsum, axysum, ayysum   ! sum of area*x*x, area*x*y, area*y*y

      real (kind=dbl_kind), dimension (nx_block, ny_block, ntrace) ::     &
           atsum            ,&! sum of area*tracer
           atxsum           ,&! sum of area*tracer*x
           atysum             ! sum of area*tracer*y

    !-------------------------------------------------------------------
    ! Initialize
    !-------------------------------------------------------------------

      aiflx(:,:) = c0
      if (present(atflx)) then
         do nt = 1, ntrace
            atflx(:,:,nt) = c0
         enddo
      endif

    !-------------------------------------------------------------------
    ! Main loop
    !-------------------------------------------------------------------

      do ng = 1, ngroups

         if (integral_order == 1) then  ! linear (1-point formula)

!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
            do ij = 1, icells(ng)
               i = indxi(ij,ng)
               j = indxj(ij,ng)

               i2 = iflux(i,j,ng)
               j2 = jflux(i,j,ng)

               ! area transports

               a0 = aic(i2,j2) + xp0(i,j,ng)*aix(i2,j2)     &
                               + yp0(i,j,ng)*aiy(i2,j2)
               asum(i,j) = a0

               aiflx(i,j) = aiflx(i,j) + triarea(i,j,ng)*asum(i,j)

               ! quantities needed for tracer transports
               axsum(i,j)  =         a0*xp0(i,j,ng) 
               axxsum(i,j) = axsum(i,j)*xp0(i,j,ng) 
               axysum(i,j) = axsum(i,j)*yp0(i,j,ng) 
               aysum(i,j)  =         a0*yp0(i,j,ng) 
               ayysum(i,j) = aysum(i,j)*yp0(i,j,ng) 
            enddo               ! ij

         elseif (integral_order == 2) then  ! quadratic (3-point formula)

!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
            do ij = 1, icells(ng)
               i = indxi(ij,ng)
               j = indxj(ij,ng)

               i2 = iflux(i,j,ng)
               j2 = jflux(i,j,ng)

               ! area transports
               ! Weighting factor of 1/3 is incorporated into the ice
               ! area terms a1, a2, and a3.
               a1 = p333 * (aic(i2,j2) + xp1(i,j,ng)*aix(i2,j2)     &
                                       + yp1(i,j,ng)*aiy(i2,j2))
               a2 = p333 * (aic(i2,j2) + xp2(i,j,ng)*aix(i2,j2)     &
                                       + yp2(i,j,ng)*aiy(i2,j2))
               a3 = p333 * (aic(i2,j2) + xp3(i,j,ng)*aix(i2,j2)     &
                                       + yp3(i,j,ng)*aiy(i2,j2))
               asum(i,j) = a1 + a2 + a3
               aiflx(i,j) = aiflx(i,j) + triarea(i,j,ng)*asum(i,j)

               ! quantities needed for volume transports
               w1 = a1 * xp1(i,j,ng)
               w2 = a2 * xp2(i,j,ng)
               w3 = a3 * xp3(i,j,ng)

               axsum(i,j) = w1 + w2 + w3

               axxsum(i,j) = w1*xp1(i,j,ng) + w2*xp2(i,j,ng)     &
                           + w3*xp3(i,j,ng) 

               axysum(i,j) = w1*yp1(i,j,ng) + w2*yp2(i,j,ng)     &
                           + w3*yp3(i,j,ng)

               w1 = a1 * yp1(i,j,ng)
               w2 = a2 * yp2(i,j,ng)
               w3 = a3 * yp3(i,j,ng)

               aysum(i,j) = w1 + w2 + w3

               ayysum(i,j) = w1*yp1(i,j,ng) + w2*yp2(i,j,ng)     &
                           + w3*yp3(i,j,ng)
            enddo               ! ij

         else                   ! cubic (4-point formula)

!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
            do ij = 1, icells(ng)
               i = indxi(ij,ng)
               j = indxj(ij,ng)

               i2 = iflux(i,j,ng)
               j2 = jflux(i,j,ng)

               ! area transports

               ! Weighting factors are incorporated into the ice
               ! area terms a0, a1, a2, and a3.
               a0 = p5625m * (aic(i2,j2) + xp0(i,j,ng)*aix(i2,j2)     &
                                         + yp0(i,j,ng)*aiy(i2,j2))
               a1 = p52083 * (aic(i2,j2) + xp1(i,j,ng)*aix(i2,j2)     &
                                         + yp1(i,j,ng)*aiy(i2,j2))
               a2 = p52083 * (aic(i2,j2) + xp2(i,j,ng)*aix(i2,j2)     &
                                         + yp2(i,j,ng)*aiy(i2,j2))
               a3 = p52083 * (aic(i2,j2) + xp3(i,j,ng)*aix(i2,j2)     &
                                         + yp3(i,j,ng)*aiy(i2,j2))
               asum(i,j) = a0 + a1 + a2 + a3
               aiflx(i,j) = aiflx(i,j) + triarea(i,j,ng)*asum(i,j)

               ! quantities needed for tracer transports
               w0 = a0 * xp0(i,j,ng)
               w1 = a1 * xp1(i,j,ng)
               w2 = a2 * xp2(i,j,ng)
               w3 = a3 * xp3(i,j,ng)

               axsum(i,j) = w0 + w1 + w2 + w3

               axxsum(i,j) = w0*xp0(i,j,ng) + w1*xp1(i,j,ng)     &
                           + w2*xp2(i,j,ng) + w3*xp3(i,j,ng)

               axysum(i,j) = w0*yp0(i,j,ng) + w1*yp1(i,j,ng)     &
                           + w2*yp2(i,j,ng) + w3*yp3(i,j,ng)

               w0 = a0 * yp0(i,j,ng)
               w1 = a1 * yp1(i,j,ng)
               w2 = a2 * yp2(i,j,ng)
               w3 = a3 * yp3(i,j,ng)

               aysum(i,j) = w0 + w1 + w2 + w3

               ayysum(i,j) = w0*yp0(i,j,ng) + w1*yp1(i,j,ng)     &
                           + w2*yp2(i,j,ng) + w3*yp3(i,j,ng)

            enddo               ! ij


         endif                  ! integral_order

         ! area * tracer transports

         if (present(atflx)) then

            do nt = 1, ntrace
               if (tracer_type(nt)==1) then ! does not depend on another tracer

!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
                  do ij = 1, icells(ng)
                     i = indxi(ij,ng)
                     j = indxj(ij,ng)

                     i2 = iflux(i,j,ng)
                     j2 = jflux(i,j,ng)

                     atsum(i,j,nt) =  asum(i,j) * trc(i2,j2,nt)     &
                                   + axsum(i,j) * trx(i2,j2,nt)     &
                                   + aysum(i,j) * try(i2,j2,nt)

                     atflx(i,j,nt) = atflx(i,j,nt)     &
                                 + triarea(i,j,ng) * atsum(i,j,nt)

                     ! quantities needed for dependent tracers

                     atxsum(i,j,nt) =  axsum(i,j) * trc(i2,j2,nt)     &
                                    + axxsum(i,j) * trx(i2,j2,nt)     &
                                    + axysum(i,j) * try(i2,j2,nt)

                     atysum(i,j,nt) =  aysum(i,j) * trc(i2,j2,nt)     &
                                    + axysum(i,j) * trx(i2,j2,nt)     &
                                    + ayysum(i,j) * try(i2,j2,nt)
                  enddo         ! ij

               elseif (tracer_type(nt)==2) then ! depends on another tracer
                  nt1 = depend(nt)

!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
                  do ij = 1, icells(ng)
                     i = indxi(ij,ng)
                     j = indxj(ij,ng)

                     i2 = iflux(i,j,ng)
                     j2 = jflux(i,j,ng)

                     atsum(i,j,nt) =  atsum(i,j,nt1) * trc(i2,j2,nt)     &
                                   + atxsum(i,j,nt1) * trx(i2,j2,nt)     &
                                   + atysum(i,j,nt1) * try(i2,j2,nt)

                     atflx(i,j,nt) = atflx(i,j,nt)     &
                                   + triarea(i,j,ng) * atsum(i,j,nt)
                  enddo         ! ij


               elseif (tracer_type(nt)==3) then ! depends on two tracers
                  nt1 = depend(nt)

!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
                  do ij = 1, icells(ng)
                     i = indxi(ij,ng)
                     j = indxj(ij,ng)

                     i2 = iflux(i,j,ng)
                     j2 = jflux(i,j,ng)

                     ! upwind approx (trx=try=0) for type 3 tracers
                     atsum(i,j,nt) =  atsum(i,j,nt1) * trc(i2,j2,nt)

                     atflx(i,j,nt) = atflx(i,j,nt)     &
                                   + triarea(i,j,ng) * atsum(i,j,nt)
                  enddo         ! ij

               endif            ! tracer type
            enddo               ! ntrace
         endif                  ! present(atflx)
      enddo                     ! ng

      end subroutine transport_integrals

!=======================================================================
!
!BOP
!
! !IROUTINE: update_fields - compute new area and tracers
!
! !INTERFACE:
!
      subroutine update_fields (nx_block, ny_block,     &
                                nghost,   tarear,     &
                                l_stop,   istop,    jstop,     &
                                aiflxe,   aiflxn,   aim,     &
                                atflxe,   atflxn,   trm)
!
! !DESCRIPTION:
!
! Given transports through cell edges, compute new area and tracers.
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

      real (kind=dbl_kind), dimension (nx_block, ny_block),     &
         intent(in) ::     &
         aiflxe, aiflxn ,&! aice transport across east and north cell edges
         tarear           ! 1/tarea

      real (kind=dbl_kind), dimension (nx_block, ny_block),     &
         intent(inout) ::     &
         aim              ! mean ice area

      real (kind=dbl_kind), dimension (nx_block, ny_block, ntrace),     &
         intent(in), optional ::     &
         atflxe, atflxn   ! aice*tracer transport across E and N cell edges

      real (kind=dbl_kind), dimension (nx_block, ny_block, ntrace),     &
         intent(inout), optional ::     &
         trm              ! mean tracers

      logical (kind=log_kind), intent(inout) ::     &
         l_stop       ! if true, abort on return

      integer (kind=int_kind), intent(inout) ::     &
         istop, jstop     ! indices of grid cell where model aborts 

!
!EOP
!
      integer (kind=int_kind) ::     &
         i, j           ,&! horizontal indices
         nt, nt1, nt2   ,&! tracer indices
         ilo,ihi,jlo,jhi  ! beginning and end of physical domain

      real (kind=dbl_kind), dimension(nx_block,ny_block,ntrace) ::     &
         atold            ! starting area*tracer

      real (kind=dbl_kind) ::     &
         w1, w2           ! work variables

      integer (kind=int_kind), dimension(nx_block*ny_block) ::     &
         indxi          ,&! compressed indices in i and j directions
         indxj

      integer (kind=int_kind) ::     &
         icells         ,&! number of cells with aim > 0.
         ij               ! combined i/j horizontal index

    !-------------------------------------------------------------------
    ! Initialize
    !-------------------------------------------------------------------

      ilo = 1 + nghost
      ihi = nx_block - nghost
      jlo = 1 + nghost
      jhi = ny_block - nghost

    !-------------------------------------------------------------------
    ! Save starting values of area*tracer
    !-------------------------------------------------------------------

      if (present(trm)) then
         do nt = 1, ntrace
            if (tracer_type(nt)==1) then ! does not depend on other tracers
               do j = jlo, jhi
               do i = ilo, ihi
                  atold(i,j,nt) = aim(i,j) * trm(i,j,nt)
               enddo            ! i
               enddo              ! j
            elseif (tracer_type(nt)==2) then  ! depends on another tracer
               nt1 = depend(nt)
               do j = jlo, jhi
               do i = ilo, ihi
                  atold(i,j,nt) = aim(i,j) * trm(i,j,nt1) * trm(i,j,nt)
               enddo            ! i
               enddo            ! j
            elseif (tracer_type(nt)==3) then  ! depends on two tracers
               nt1 = depend(nt)
               nt2 = depend(nt1)
               do j = jlo, jhi
               do i = ilo, ihi
                  atold(i,j,nt) = aim(i,j)      &
                            * trm(i,j,nt2) * trm(i,j,nt1) * trm(i,j,nt)
               enddo            ! i
               enddo            ! j

 
            endif               ! depend(nt) = 0
         enddo                  ! nt
      endif                     ! present(trm)

    !-------------------------------------------------------------------
    ! Update ice area
    !-------------------------------------------------------------------

      do j = jlo, jhi
      do i = ilo, ihi

         w1 = aiflxe(i,j) - aiflxe(i-1,j)     &
            + aiflxn(i,j) - aiflxn(i,j-1)
         aim(i,j) = aim(i,j) - w1*tarear(i,j)

         if (aim(i,j) < -puny) then    ! abort with negative value
            l_stop = .true.
            istop = i
            jstop = j
         elseif (aim(i,j) < c0) then   ! set to zero
            aim(i,j) = c0
         endif

      enddo
      enddo

      if (l_stop) then
         i = istop
         j = jstop
         w1 = aiflxe(i,j) - aiflxe(i-1,j)     &
            + aiflxn(i,j) - aiflxn(i,j-1)
         write (nu_diag,*) ' '
         write (nu_diag,*) 'New area < 0, i, j =', i, j
         write (nu_diag,*) 'Old area =', aim(i,j) + w1*tarear(i,j)
         write (nu_diag,*) 'New area =', aim(i,j)
         write (nu_diag,*) 'Net transport =', -w1*tarear(i,j)
         return
      endif

    !-------------------------------------------------------------------
    ! Update tracers
    !-------------------------------------------------------------------

      if (present(trm)) then

         icells = 0
         do j = jlo, jhi
         do i = ilo, ihi
            if (aim(i,j) > c0) then ! grid cells with positive areas
               icells = icells + 1
               indxi(icells) = i
               indxj(icells) = j
            endif
         enddo                  ! i
         enddo                  ! j

         do nt = 1, ntrace

            trm(:,:,nt) = c0

            if (tracer_type(nt)==1) then ! does not depend on other tracers

               do ij = 1, icells
                  i = indxi(ij)
                  j = indxj(ij)

                  w1  = atflxe(i,j,nt) - atflxe(i-1,j,nt)     &
                      + atflxn(i,j,nt) - atflxn(i,j-1,nt)
                  trm(i,j,nt) = (atold(i,j,nt) - w1*tarear(i,j))     &
                                / aim(i,j)
               enddo            ! ij

            elseif (tracer_type(nt)==2) then ! depends on another tracer
               nt1 = depend(nt)

!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
               do ij = 1, icells
                  i = indxi(ij)
                  j = indxj(ij)

                  if (abs(trm(i,j,nt1)) > puny) then
                     w1  = atflxe(i,j,nt) - atflxe(i-1,j,nt)     &
                         + atflxn(i,j,nt) - atflxn(i,j-1,nt)
                     trm(i,j,nt) = (atold(i,j,nt) - w1*tarear(i,j))     &
                                 / (aim(i,j) * trm(i,j,nt1))
                  endif

               enddo            ! ij

            elseif (tracer_type(nt)==3) then ! depends on two tracers
               nt1 = depend(nt)
               nt2 = depend(nt1)

!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
               do ij = 1, icells
                  i = indxi(ij)
                  j = indxj(ij)

                  if (abs(trm(i,j,nt1)) > puny .and.     &
                      abs(trm(i,j,nt2)) > puny) then
                     w1  = atflxe(i,j,nt) - atflxe(i-1,j,nt)     &
                         + atflxn(i,j,nt) - atflxn(i,j-1,nt)
                     trm(i,j,nt) = (atold(i,j,nt) - w1*tarear(i,j))     &
                              / (aim(i,j) * trm(i,j,nt2) * trm(i,j,nt1))
                  endif
               enddo            ! ij

            endif               ! tracer_type
         enddo                  ! nt
      endif                     ! present(trm)

      end subroutine update_fields

!=======================================================================
!
!BOP
!
! !IROUTINE: limited_gradient - limited gradient of a scalar field
!
! !INTERFACE:
!
      subroutine limited_gradient (nx_block, ny_block,       &
                                   icells,   ilim,           &
                                   indxii, indxjj, indxij,   &
                                   phi,      phimask,        &
                                   cnx,      cny,            &
                                   HTN,      HTE,            &
                                   dxt,      dyt,            &
                                   gx,       gy)
!
! !DESCRIPTION:
!
! Compute a limited gradient of the scalar field phi.
! "Limited" means that we do not create new extrema in phi.  For
! instance, field values at the cell corners can neither exceed the
! maximum of phi(i,j) in the cell and its eight neighbors, nor fall
! below the minimum.
!
! !REVISION HISTORY:
!
! authors William H. Lipscomb, LANL
!         John R. Baumgardner, LANL
!
! !USES:
!
! !INPUT/OUTPUT PARAMETERS:
!
      integer (kind=int_kind), intent(in) ::     &
          nx_block, ny_block,&! block dimensions
          ilim              ,&
          icells

      integer (kind=int_kind), dimension (nx_block*ny_block),     &
          intent(in) ::     &
          indxii, indxjj ,&! combined i/j horizontal indices
          indxij

      real (kind=dbl_kind), dimension (icells),     &
          intent (in) ::     &
          cnx    ,&! x-coordinate of phi relative to geometric center of cell
          cny      ! y-coordinate of phi relative to geometric center of cell

      real (kind=dbl_kind), dimension (nx_block,ny_block),     &
          intent (in) ::     &
          phi    ,&! input tracer field (mean values in each grid cell)
          dxt    ,&! grid cell width (m)
          dyt    ,&! grid cell height (m)
          phimask,&
          ! phimask(i,j) = 1 if phi(i,j) has physical meaning, = 0 otherwise.
          ! For instance, aice has no physical meaning on land points,
          ! and hice no physical meaning where aice = 0.
          HTN    ,&! length of northern edge of T-cell (m)
          HTE      ! length of eastern edge of T-cell (m)

      real (kind=dbl_kind), dimension (nx_block,ny_block),     &
          intent(out) ::     &
          gx     ,&! limited x-direction gradient
          gy       ! limited y-direction gradient
!
!EOP
!
      integer (kind=int_kind) ::     &
          i, j

      real (kind=dbl_kind) ::     &
          phi_nw, phi_n, phi_ne ,&! values of phi in 8 neighbor cells
          phi_w,         phi_e  ,&
          phi_sw, phi_s, phi_se ,& 
          qmn, qmx     ,&! min and max value of phi within grid cell
          pmn, pmx     ,&! min and max value of phi among neighbor cells
          w1, w2, w3, w4 ! work variables

      integer (kind=int_kind) ::     &
          ij, m          ! combined i/j horizontal indices

      real (kind=dbl_kind) ::     &
          gxtmp, gytmp   ! temporary term for x- and y- limited gradient

      do j = 1, ny_block
      do i = 1, nx_block
         gx(i,j) = c0
         gy(i,j) = c0
      enddo
      enddo

      do ij = 1, ilim
         i = indxii(ij)
         j = indxjj(ij)
         m = indxij(ij)

         ! Store values of phi in the 8 neighbor cells.
         ! Note: phimask = 1. or 0.  If phimask = 1., use the true value;
         !  if phimask = 0., use the home cell value so that non-physical
         !  values of phi do not contribute to the gradient.
         phi_nw = phimask(i-1,j+1) * phi(i-1,j+1)     &
            + (c1-phimask(i-1,j+1))* phi(i,j)
         phi_n  = phimask(i,j+1)   * phi(i,j+1)       &
            + (c1-phimask(i,j+1))  * phi(i,j)
         phi_ne = phimask(i+1,j+1) * phi(i+1,j+1)     &
            + (c1-phimask(i+1,j+1))* phi(i,j)
         phi_w  = phimask(i-1,j)   * phi(i-1,j)       &
            + (c1-phimask(i-1,j))  * phi(i,j)
         phi_e  = phimask(i+1,j)   * phi(i+1,j)       &
            + (c1-phimask(i+1,j))  * phi(i,j)
         phi_sw = phimask(i-1,j-1) * phi(i-1,j-1)     &
            + (c1-phimask(i-1,j-1))* phi(i,j)
         phi_s  = phimask(i,j-1)   * phi(i,j-1)       &
            + (c1-phimask(i,j-1))  * phi(i,j)
         phi_se = phimask(i+1,j-1) * phi(i+1,j-1)     &
            + (c1-phimask(i+1,j-1))* phi(i,j)

         ! unlimited gradient components
         ! (factors of two cancel out)
         gxtmp = (phi_e - phi(i,j)) / (dxt(i,j)   + dxt(i+1,j))     &
               + (phi(i,j) - phi_w) / (dxt(i-1,j) + dxt(i,j)  )
         gytmp = (phi_n - phi(i,j)) / (dyt(i,j)   + dyt(i,j+1))     &
               + (phi(i,j) - phi_s) / (dyt(i,j-1) + dyt(i,j)  )

         ! minimum and maximum among the nine local cells
         pmn = min (phi_nw, phi_n,  phi_ne, phi_w, phi(i,j),     &
                    phi_e,  phi_sw, phi_s,  phi_se)
         pmx = max (phi_nw, phi_n,  phi_ne, phi_w, phi(i,j),     &
                    phi_e,  phi_sw, phi_s,  phi_se)

         pmn = pmn - phi(i,j)
         pmx = pmx - phi(i,j)

         ! minimum and maximum deviation of phi within the cell
         w1  =  (p5*HTN(i,j)   - cnx(m)) * gxtmp     &
              + (p5*HTE(i,j)   - cny(m)) * gytmp
         w2  =  (p5*HTN(i,j-1) - cnx(m)) * gxtmp     &
              - (p5*HTE(i,j)   + cny(m)) * gytmp
         w3  = -(p5*HTN(i,j-1) + cnx(m)) * gxtmp     &
              - (p5*HTE(i-1,j) + cny(m)) * gytmp
         w4  =  (p5*HTE(i-1,j) - cny(m)) * gytmp     &
              - (p5*HTN(i,j)   + cnx(m)) * gxtmp

         qmn = min (w1, w2, w3, w4)
         qmx = max (w1, w2, w3, w4)

         ! the limiting coefficient
         if (abs(qmn) > c0) then ! 'abs(qmn) > puny' not sufficient
            w1 = max(c0, pmn/qmn)
         else
            w1 = c1
         endif

         if (abs(qmx) > c0) then
            w2 = max(c0, pmx/qmx)
         else
            w2 = c1
         endif

         w1 = min(c1, w1, w2)

         ! Limit the gradient components
         gx(i,j) = w1 * gxtmp
         gy(i,j) = w1 * gytmp

      enddo                     ! ij

      end subroutine limited_gradient

!=======================================================================

      end module ice_transport_remap

!=======================================================================
