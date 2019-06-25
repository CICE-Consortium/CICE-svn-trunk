c $Id: ice_transport_remap.F,v 1.19 2004/03/01 16:56:24 eclare Exp $
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
! !REVISION HISTORY:
!
! authors William H. Lipscomb, LANL
!         John Baumgardner, LANL
!
! Vectorized by Clifford Chen (Fujitsu) and William Lipscomb (LANL)
!
! !INTERFACE:
!
      module ice_transport_remap
!
! !USES:
!
      use ice_model_size
      use ice_kinds_mod
      use ice_domain
      use ice_constants
      use ice_grid
      use ice_fileunits
      use ice_calendar
      use ice_state
      use ice_timers
      use ice_itd
      use ice_work, only:  work_l1
!
!EOP
!
      implicit none

!=======================================================================
!
! INSTRUCTIONS FOR CHANGING THE SET OF TRANSPORTED VARIABLES:
!
! The default version of remapping assumes that the transported 
! variables are ice area, ice volume, snow volume, surface temperature, 
! ice energy (in several vertical layers), and snow energy.
!
! Thus the transported tracers are hice, hsno, Tsfc, qice(1:nilyr),
! and qsno.  The number of tracers is 4 + nilyr. 
!
! To add or remove tracers from this list:
!
! (1) Change ntracer immediately below these instructions.
!
! (2) At the beginning of the initialization subroutine, init_remap,
!     change the definition of the depend vector.  This is defined
!     as follows:
!
!     Any tracer T1 satisfying d(a*T1)/dt + del*(u*a*T1) = 0
!     is said to be independent of other tracers and has depend = 0.
!     Such tracers are given indices nt = 1, 2, 3, ...
!     For example, in the default case, hice, hsno, and Tsfc have
!      nt = 1, 2, and 3, respectively.  So depend(1:3) = 0.
!
!     Any tracer T2 satisfying d(a*T1*T2)/dt + del*(u*a*T1*T2) = 0
!     is said to be dependent on tracer T1.  It has depend = nt1,
!     where nt1 is the index of tracer T1.
!     For example, in the default case, we have depend(nt) = 1 for
!      nt = 4 to (3+nilyr), corresponding to qice(1:nilyr).
!      We have depend(nt) = 2 for nt = 4+nilyr, corresponding to qsno.
!
! (3) In subroutine load\_tracers, change the loading of tracer
!     array trm from state variables.  Make sure the tracer indices
!     in this array are consistent with the depend vector.
!
! (4) In subroutine unload\_tracers, change the conversion of tracer
!     array trm back to state variables. 
!    
!=======================================================================

      integer (kind=int_kind), parameter ::      
     &   ntracer = nilyr + 4    ! qice + hice, hsno, Tsfc, qsno

      integer (kind=int_kind), dimension (ntracer) ::
     &   depend                 ! describes tracer dependencies
                                ! (see instructions above)

      logical (kind=log_kind), dimension (ntracer) ::
     &   has_dependents         ! true if a tracer has dependent tracers

      real (kind=dbl_kind), dimension (imt_local,jmt_local) ::
     &   dxt2r      ! 1. / (dxt(i+1,j) + dxt(i,j))
     &,  dyt2r      ! 1. / (dyt(i,j+1) + dyt(i,j))
     &,  xav, yav   ! mean T-cell values of x, y  
     &,  xxav, xyav, yyav  ! mean T-cell values of xx, xy, yy  
     &,  xxxav,xxyav,xyyav,yyyav ! mean T-cell values of xxx, xxy, xyy, yyy 

      real (kind=dbl_kind), dimension (2,2,imt_local,jmt_local) ::
     &   mne, mnw   ! matrices used for coordinate transformations
     &,  msw, mse   ! ne = northeast corner, nw = northwest corner, etc.

      integer (kind=int_kind), parameter ::      
     &   ngroups  = 5        ! number of groups of triangles that
                             ! contribute fluxes across each edge

      real (kind=dbl_kind), parameter :: ! for triangle integral formulas      
     &   p5625m = -c9/c16
     &,  p52083 = c25/48._dbl_kind

      character (len=char_len) :: 
     &   advection   ! type of advection algorithm used
                     ! 'upwind'  => 1st order mpdata scheme (donor cell)
                     ! 'mpdata'  => 2nd order mpdata scheme
                     ! 'remap' (or anything else) => remapping scheme      

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
! authors William H. Lipscomb, LANL
!         John R. Baumgardner, LANL
!
! !USES:
!
      use ice_exit
!
! !INPUT/OUTPUT PARAMETERS:
!
!EOP
!
      integer (kind=int_kind) ::
     &     i, j                            ! horizontal indices
     &,    m                               ! triangle index
     &,    nt, nt1                         ! tracer index

      real (kind=dbl_kind), dimension (4) ::  
     &     ar,                             ! triangle area
     &     xv1, xv2, xv3, yv1, yv2, yv3    ! vertices

      real (kind=dbl_kind) ::
     &     x0, x1, x2, x3, y0, y1, y2, y3  ! interior points
     &,    theta1, theta2, theta3, theta4  ! angles

    !-------------------------------------------------------------------
    ! Compute tracer dependency vectors
    ! NOTE: Must change determination of 'depend' if transporting 
    !       a different set of tracers.  See instructions above.
    !-------------------------------------------------------------------

      depend(1:3)       = 0     ! hice, hsno, Tsfc
      depend(4:3+nilyr) = 1     ! qice depends on hice
      depend(nilyr+4)   = 2     ! qsno depends on hsno

      has_dependents = .false.
      do nt = 1, ntracer
         if (depend(nt) > 0) then
            nt1 = depend(nt)
            has_dependents(nt1) = .true.
            if (nt1 > nt) then
               write(nu_diag,*) 
     &              'Tracer nt2 =',nt,' depends on tracer nt1 =',nt1
               call abort_ice('remap transport: Must have nt2 > nt1')
            endif
         endif
      enddo

      !-----------------------------------------------------------------
      ! Compute cell widths used in gradient-limiting subroutine.
      !-----------------------------------------------------------------

      do j = jlo,jhi
      do i = ilo,ihi
         dxt2r(i,j) = c1 / (dxt(i+1,j) + dxt(i,j))
         dyt2r(i,j) = c1 / (dyt(i,j+1) + dyt(i,j))
      enddo
      enddo
      call bound(dxt2r)
      call bound(dyt2r)

      !-----------------------------------------------------------------
      ! Construct mean values of geometic quantities over the T-cells,
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
      ! in subroutine flux_integrals.
      !
      ! xav  = mean value of x
      ! yav  = mean value of y
      ! xxav = mean value of x^2
      ! xyav = mean value of xy
      ! yyav = mean value of y^2
      ! xxxav = mean value of x^3
      ! xxyav = mean value of x^2 * y
      ! xyyav = mean value of x * y^2
      ! yyyav = mean value of y^3
      !------------------------------------------------------------------

      do j = jlo,jhi
      do i = ilo,ihi
         ! Initialize
         xav(i,j) = c0
         yav(i,j) = c0
         xxav(i,j) = c0
         xyav(i,j) = c0
         yyav(i,j) = c0
         xxxav(i,j) = c0
         xxyav(i,j) = c0
         xyyav(i,j) = c0
         yyyav(i,j) = c0
         
         ! Compute vertices and area of each triangle.
         xv1(1) = c0            ! East triangle
         yv1(1) = c0
         xv2(1) = p5*HTN(i,j)
         yv2(1) = p5*HTE(i,j)
         xv3(1) = p5*HTS(i,j)
         yv3(1) = -p5*HTE(i,j)
         ar(1)  = p25*HTE(i,j)*dxt(i,j)
         
         xv1(2) = c0            ! West triangle
         yv1(2) = c0
         xv2(2) = -p5*HTN(i,j)
         yv2(2) = p5*HTW(i,j)
         xv3(2) = -p5*HTS(i,j)
         yv3(2) = -p5*HTW(i,j)
         ar(2)  = p25*HTW(i,j)*dxt(i,j)

         xv1(3) = c0            ! North triangle
         yv1(3) = c0
         xv2(3) = -p5*HTN(i,j)
         yv2(3) = p5*HTW(i,j)
         xv3(3) = p5*HTN(i,j)
         yv3(3) = p5*HTE(i,j)
         ar(3)  = p25*HTN(i,j)*dyt(i,j)

         xv1(4) = c0            ! South triangle
         yv1(4) = c0
         xv2(4) = -p5*HTS(i,j)
         yv2(4) = -p5*HTW(i,j)
         xv3(4) = p5*HTS(i,j)
         yv3(4) = -p5*HTE(i,j)
         ar(4)  = p25*HTS(i,j)*dyt(i,j)

         ! Compute contribution to means from each triangle (E, W, N, S)
         do m = 1, 4

            x0 = (xv1(m) + xv2(m) + xv3(m)) / c3 ! midpoint
            y0 = (yv1(m) + yv2(m) + yv3(m)) / c3
            x1 = p6*x0 + p4*xv1(m) ! other 3 points needed for integral 
            x2 = p6*x0 + p4*xv2(m)
            x3 = p6*x0 + p4*xv3(m)
            y1 = p6*y0 + p4*yv1(m)
            y2 = p6*y0 + p4*yv2(m)
            y3 = p6*y0 + p4*yv3(m)

            xav(i,j) = xav(i,j) + ar(m)*x0
            yav(i,j) = yav(i,j) + ar(m)*y0
            xxav(i,j) = xxav(i,j) + ar(m) * 
     &                  ( p5625m * x0*x0 
     &                  + p52083 * (x1*x1 + x2*x2 + x3*x3) )
            xyav(i,j) = xyav(i,j) + ar(m) * 
     &                  ( p5625m * x0*y0 
     &                  + p52083 * (x1*y1 + x2*y2 + x3*y3) )
            yyav(i,j) = yyav(i,j) + ar(m) * 
     &                  ( p5625m * y0*y0 
     &                  + p52083 * (y1*y1 + y2*y2 + y3*y3) )
            xxxav(i,j) = xxxav(i,j) + ar(m) * 
     &                  ( p5625m * x0*x0*x0
     &                  + p52083 * (x1*x1*x1 + x2*x2*x2 + x3*x3*x3) )
            xxyav(i,j) = xxyav(i,j) + ar(m) * 
     &                  ( p5625m * x0*x0*y0
     &                  + p52083 * (x1*x1*y1 + x2*x2*y2 + x3*x3*y3) )
            xyyav(i,j) = xyyav(i,j) + ar(m) * 
     &                  ( p5625m * x0*y0*y0
     &                  + p52083 * (x1*y1*y1 + x2*y2*y2 + x3*y3*y3) )
            yyyav(i,j) = yyyav(i,j) + ar(m) * 
     &                  ( p5625m * y0*y0*y0
     &                  + p52083 * (y1*y1*y1 + y2*y2*y2 + y3*y3*y3) )

         enddo                  ! m (loop over 4 triangles)

         ! Divide by grid cell area
         xav(i,j) = xav(i,j) * tarear(i,j) 
         yav(i,j) = yav(i,j) * tarear(i,j)
         xxav(i,j) = xxav(i,j) * tarear(i,j)
         xyav(i,j) = xyav(i,j) * tarear(i,j)
         yyav(i,j) = yyav(i,j) * tarear(i,j)
         xxxav(i,j) = xxxav(i,j) * tarear(i,j)
         xxyav(i,j) = xxyav(i,j) * tarear(i,j)
         xyyav(i,j) = xyyav(i,j) * tarear(i,j)
         yyyav(i,j) = yyyav(i,j) * tarear(i,j)
         
      enddo                     ! i
      enddo                     ! j

      ! Compute ghost cell values
      call bound(xav)
      call bound(yav)
      call bound(xxav)
      call bound(xyav)
      call bound(yyav)
      call bound(xxxav)
      call bound(xxyav)
      call bound(xyyav)
      call bound(yyyav)

      !-----------------------------------------------------------------
      ! Compute matrix elements needed to transform from a reference frame
      ! whose origin is at a cell corner (a U-cell frame) to a reference 
      ! frame whose origin is at a cell center (a T-cell frame).
      ! Transformations are needed because the axes of these two reference
      ! frames are not parallel on a curved grid.
      !----------------------------------------------------------------- 

      do j = 1, jmt_local
      do i = 1, imt_local

      !----------------------------------------------------------------- 
      ! Compute angles between the U-cell coordinate axes and T-cell
      ! coordinate axes.  The U-cell axes lie along the cell edges,
      ! which do not meet at right angles in the T-cell reference frame.
      ! The angles are defined as positive if the U-cell axes are obtained 
      ! by a counterclockwise rotation from the T-cell axes.  
      ! Four independent angles are needed.
      !-----------------------------------------------------------------

         theta1 = atan(p5*(HTE(i,j)-HTW(i,j))/HTN(i,j))
         theta2 = atan(p5*(HTN(i,j)-HTS(i,j))/HTW(i,j))
         theta3 = atan(p5*(HTW(i,j)-HTE(i,j))/HTS(i,j))
         theta4 = atan(p5*(HTS(i,j)-HTN(i,j))/HTE(i,j))

      !-----------------------------------------------------------------
      ! Compute matrix elements.
      !-----------------------------------------------------------------

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

      enddo
      enddo

      end subroutine init_remap

!=======================================================================
!
!BOP
!
! !IROUTINE: transport_remap - remapping transport scheme
!
! !INTERFACE:
!
      subroutine transport_remap
!
! !DESCRIPTION:
! 
! This subroutine solves the transport equations for one timestep 
! using the conservative remapping scheme developed by John Dukowicz 
! and John Baumgardner (DB) and modified for sea ice by William 
! Lipscomb and Elizabeth Hunke.  
!
! This scheme is second-order accurate, except where gradients are
! limited to preserve monotonicity.  It is compatible for tracers; 
! that is, it does not produce new extrema in thickness or enthalpy.  
!
! References:
! 
! Dukowicz, J. K., and J. R. Baumgardner, 2000: Incremental 
!  remapping as a transport/advection algorithm, J. Comput. Phys., 
!  160, 318-335.
!
! Lipscomb, W. H., and E. C. Hunke, 2003: Modeling sea ice
!  transport using incremental remapping, accepted by Mon. Wea. Rev.
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
!EOP
!
      integer (kind=int_kind) ::
     &   i, j             ! horizontal indices
     &,  k                ! layer index
     &,  n                ! ice category index

      real (kind=dbl_kind), dimension (imt_local,jmt_local) ::
     &   dpx              ! x coordinates of departure points at cell corners
     &,  dpy              ! y coordinates of departure points at cell corners

      real (kind=dbl_kind), dimension (imt_local, jmt_local) ::
     &   aimask           ! = 1. if ice is present, = 0. otherwise
     &,  aic              ! ice area at geometric center of cell
     &,  aix, aiy         ! limited derivative of ice area wrt x and y
     &,  aiflxe, aiflxn   ! flux of ice area thru east and north cell edges

      real (kind=dbl_kind), dimension (imt_local, jmt_local, ntracer) ::
     &   trm              ! mean tracer values in each grid cell
     &,  trmask           ! = 1. if tracer is present, = 0. otherwise
     &,  trc              ! tracer values at geometric center of cell
     &,  trx, try         ! limited derivative of tracer wrt x and y
     &,  atflxe, atflxn   ! flux of area*tracer thru E and N cell edges

      real (kind=dbl_kind), dimension (imt_local, jmt_local, ngroups) ::
     &   triarea_e        ! area of east-edge flux triangle
     &,  triarea_n        ! area of north-edge flux triangle
     &,  xp0_e, yp0_e     ! x and y coordinates of special triangle points
     &,  xp1_e, yp1_e     ! (need 4 points for triangle integrals)
     &,  xp2_e, yp2_e   
     &,  xp3_e, yp3_e     ! e for east edges, n for north edges
     &,  xp0_n, yp0_n
     &,  xp1_n, yp1_n
     &,  xp2_n, yp2_n
     &,  xp3_n, yp3_n

      integer (kind=int_kind), 
     &   dimension (imt_local,jmt_local,ngroups) ::
     &   iflux_e          ! i index of cell contributing east flux
     &,  jflux_e          ! j index of cell contributing east flux
     &,  iflux_n          ! i index of cell contributing north flux
     &,  jflux_n          ! j index of cell contributing north flux

      ! variables used for optional bug checks

      logical (kind=log_kind), parameter ::
     &   l_conservation_check = .false. ! if true, check conservation 
     &,  l_monotonicity_check = .false. ! if true, check monotonicity 
     
      real (kind=dbl_kind) ::
     &   asum_init        ! initial global ice area
     &,  asum_final       ! final global ice area

      real (kind=dbl_kind), dimension(ntracer) ::
     &   atsum_init       ! initial global ice area*tracer
     &,  atsum_final      ! final global ice area*tracer

      real (kind=dbl_kind), dimension (:,:,:), allocatable ::
     &   tmin           ! local min tracer
     &,  tmax           ! local max tracer

      call ice_timer_start(3)   ! advection

    !-------------------------------------------------------------------
    ! Get ghost cell velocities.
    !-------------------------------------------------------------------

      call bound(uvel)
      call bound(vvel)

    !-------------------------------------------------------------------
    ! Given velocity field at cell corners, compute departure points
    ! of trajectories.
    !-------------------------------------------------------------------

      call departure_points (dpx, dpy)

    !-------------------------------------------------------------------
    ! Fill in values along south and west boundaries of physical domain.
    !-------------------------------------------------------------------

      call bound_sw (dpx)
      call bound_sw (dpy)

    !-------------------------------------------------------------------
    ! Compute areas and vertices of flux triangles for east 
    !  and north cell edges.
    !-------------------------------------------------------------------

      call locate_triangles (dpx, dpy,
     &        xp0_e, xp1_e, xp2_e, xp3_e, yp0_e, yp1_e, yp2_e, yp3_e, 
     &        xp0_n, xp1_n, xp2_n, xp3_n, yp0_n, yp1_n, yp2_n, yp3_n,
     &        triarea_e, iflux_e, jflux_e,
     &        triarea_n, iflux_n, jflux_n)

    !-------------------------------------------------------------------
    ! Given triangle vertices, compute coordinates of triangle points 
    !  needed for flux integrals.
    !-------------------------------------------------------------------

      call triangle_coordinates (triarea_e,
     &                           xp0_e, xp1_e, xp2_e, xp3_e,
     &                           yp0_e, yp1_e, yp2_e, yp3_e)

      call triangle_coordinates (triarea_n,
     &                           xp0_n, xp1_n, xp2_n, xp3_n,
     &                           yp0_n, yp1_n, yp2_n, yp3_n)

    !-------------------------------------------------------------------
    ! Compute open water area in each grid cell.
    ! Note: An aggregate_area call is needed only if the open
    !       water area has changed since the previous call.
    !       Here we assume that aice(i,j) is up to date.
    !-------------------------------------------------------------------
c      call aggregate_area
c      call bound (aice0)

    !-------------------------------------------------------------------
    ! Get ghost cell values of state variables.
    !-------------------------------------------------------------------

      call bound_state

!---!-------------------------------------------------------------------
!---! Transport the open water area.
!---!-------------------------------------------------------------------

    !-------------------------------------------------------------------
    ! Compute mask for open water area.
    !-------------------------------------------------------------------

      call make_masks (aice0, aimask)

    !-------------------------------------------------------------------
    ! Compute initial value of global open water area.
    !-------------------------------------------------------------------

      if (l_conservation_check) then
         call conserved_sums (aice0, asum_init)
      endif

    !-------------------------------------------------------------------
    ! Construct linear fields of open water area, limiting gradients
    ! to preserve monotonicity.
    !-------------------------------------------------------------------

      call construct_fields (aice0, aic,  aix,  aiy,  aimask)

      call bound (aix)
      call bound (aiy)
      call bound (aic)

    !-------------------------------------------------------------------
    ! Compute the flux of open water area across the cell edges by 
    ! summing contributions from each triangle.
    !-------------------------------------------------------------------

      ! east edges
      call flux_integrals (triarea_e, iflux_e,   jflux_e,
     &                     xp0_e,     xp1_e,     xp2_e,    xp3_e,
     &                     yp0_e,     yp1_e,     yp2_e,    yp3_e,
     &                     aic,       aix,       aiy,      aiflxe)

      ! north edges
      call flux_integrals (triarea_n, iflux_n,   jflux_n,
     &                     xp0_n,     xp1_n,     xp2_n,    xp3_n,
     &                     yp0_n,     yp1_n,     yp2_n,    yp3_n,
     &                     aic,       aix,       aiy,      aiflxn)

    !-------------------------------------------------------------------
    ! Transfer fluxes across cell edges and update the open water area.
    !-------------------------------------------------------------------

      call update_fields (0, aiflxe, aiflxn, aice0)

    !-------------------------------------------------------------------
    ! Compute final value of global open water area and check conservation.
    ! (Optional)
    !-------------------------------------------------------------------

      if (l_conservation_check) then
         call conserved_sums (aice0, asum_final)
         if (my_task.eq.master_task) 
     &        call global_conservation (asum_init,  asum_final)
      endif

!---!-------------------------------------------------------------------
!---! Loop over ice categories
!---!-------------------------------------------------------------------

      do n = 1, ncat

    !-------------------------------------------------------------------
    ! Load tracer array and get ghost cell values.
    ! Assume that the transported tracers are hice, hsno, Tsfc,
    !  qice(1:nilyr), and qsno.  Must change this subroutine if
    !  transporting a different set of tracers.
    !-------------------------------------------------------------------

         call load_tracers (trm, n)

    !-------------------------------------------------------------------
    ! Compute masks for ice area and tracers.
    !-------------------------------------------------------------------

         call make_masks (aicen(:,:,n), aimask, trm, trmask)
         
    !-------------------------------------------------------------------
    ! Compute initial values of globally conserved quantities.
    !-------------------------------------------------------------------

         if (l_conservation_check) then
            call conserved_sums (aicen(:,:,n), asum_init, 
     &                           trm,          atsum_init)
         endif

    !-------------------------------------------------------------------
    ! Compute local max and min of tracer fields.
    !-------------------------------------------------------------------
      
         if (l_monotonicity_check) then
            allocate(tmin(imt_local,jmt_local,ntracer),
     &               tmax(imt_local,jmt_local,ntracer))
            call local_max_min (trm, tmin, tmax, aimask, trmask)
         endif

    !-------------------------------------------------------------------
    ! Construct linear fields of ice area and tracers, limiting gradients
    ! to preserve monotonicity.
    !-------------------------------------------------------------------

         call construct_fields
     &        (aicen(:,:,n), aic,  aix,  aiy,  aimask,
     &         trm,          trc,  trx,  try,  trmask)

         call bound (aix)
         call bound (aiy)
         call bound (aic)

         call bound_narr (ntracer, trx)
         call bound_narr (ntracer, try)
         call bound_narr (ntracer, trc)
         
    !-------------------------------------------------------------------
    ! Compute the flux of ice area and tracers across the cell edges by 
    ! summing contributions from each triangle.
    !-------------------------------------------------------------------

         ! east edges
         call flux_integrals 
     &                    (triarea_e, iflux_e,   jflux_e,
     &                     xp0_e,     xp1_e,     xp2_e,    xp3_e,
     &                     yp0_e,     yp1_e,     yp2_e,    yp3_e,
     &                     aic,       aix,       aiy,      aiflxe, 
     &                     trc,       trx,       try,      atflxe)

         ! north edges
         call flux_integrals
     &                    (triarea_n, iflux_n,   jflux_n,
     &                     xp0_n,     xp1_n,     xp2_n,    xp3_n,
     &                     yp0_n,     yp1_n,     yp2_n,    yp3_n,
     &                     aic,       aix,       aiy,      aiflxn, 
     &                     trc,       trx,       try,      atflxn)

    !-------------------------------------------------------------------
    ! Transfer fluxes across cell edges and update ice area and tracers.
    !-------------------------------------------------------------------

         call update_fields (n, aiflxe, aiflxn, aicen(:,:,n),
     &                          atflxe, atflxn, trm)
         
    !-------------------------------------------------------------------
    ! Compute final values of globally conserved quantities.
    ! Check global conservation of area and area*tracers.
    ! (Optional)
    !-------------------------------------------------------------------

         if (l_conservation_check) then
            call conserved_sums (aicen(:,:,n), asum_final, 
     &                           trm,          atsum_final)
            if (my_task.eq.master_task) 
     &         call global_conservation (asum_init,  asum_final,
     &                                   atsum_init, atsum_final)
         endif

    !-------------------------------------------------------------------
    ! Check tracer monotonicity.  (Optional)
    !-------------------------------------------------------------------

         if (l_monotonicity_check) then
            call check_monotonicity (aicen(:,:,n), trm, tmin, tmax)
            deallocate(tmin,tmax)
         endif

    !-------------------------------------------------------------------
    ! Given new tracers, recompute state variables.
    !-------------------------------------------------------------------

         call unload_tracers (trm, n)

      enddo                     ! n loop

      call ice_timer_stop(3)    ! advection

    !-------------------------------------------------------------------
    ! mask
    !-------------------------------------------------------------------

      do n=1,ncat
       do j=1,jmt_local
        do i=1,imt_local
         if (.not.tmask(i,j)) then
          aicen(i,j,n) = c0
          vicen(i,j,n) = c0
          vsnon(i,j,n) = c0
          esnon(i,j,n) = c0
          Tsfcn(i,j,n) = c0
         endif
        enddo
       enddo
      enddo

      do k=1,ntilay
       do j=1,jmt_local
        do i=1,imt_local
         if (.not.tmask(i,j)) then
          eicen(i,j,k) = c0
         endif
        enddo
       enddo
      enddo

      end subroutine transport_remap

!=======================================================================
!
!BOP
!
! !IROUTINE: departure points - compute departure points of trajectories
!
! !INTERFACE:
!
      subroutine departure_points (dpx, dpy)
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
      real (kind=dbl_kind), dimension (imt_local, jmt_local), 
     &     intent(out) ::
     &   dpx              ! x coordinates of departure points at cell corners
     &,  dpy              ! y coordinates of departure points at cell corners
!
!EOP
!
      integer (kind=int_kind) ::
     &   i, j             ! horizontal indices
     &,  i2, j2           ! horizontal indices
     &,  niter            ! iteration counter
 
      real (kind=dbl_kind) ::
     &   mpx, mpy         ! coordinates of midpoint of back trajectory
     &,  u1t, v1t         ! transformed velocity, SW corner
     &,  u2t, v2t         ! transformed velocity, SE corner
     &,  u3t, v3t         ! transformed velocity, NE corner
     &,  u4t, v4t         ! transformed velocity, NW corner         
     &,  umpt, vmpt       ! midpoint velocity in transformed coordinates
     &,  ump, vmp         ! midpoint velocity in original corner coordinates
     &,  w1,w2,w3,w4      ! work variables
     
      real (kind=dbl_kind), dimension (ilo:ihi, jlo:jhi) ::
     &   cxt, cyt         ! transformed cell corner coordinates
     &,  mpxt, mpyt       ! midpoint transformed to cell-ctr coordinates
     &,  mpxs, mpys       ! midpoint in stretched coordinates
     &,  mat11, mat12     ! transformation matrix, cell corner to center
     &,  mat21, mat22     ! transformation matrix, cell corner to center

      integer (kind=int_kind), dimension (ilo:ihi, jlo:jhi) ::
     &   hindi, hindj     ! horizontal indices array

      logical (kind=log_kind) ::
     &   dept_out_bound   ! flag for departure point out of range (-1,1)

    !-------------------------------------------------------------------
    ! Initialize departure points
    !-------------------------------------------------------------------
      dpx(:,:) = c0
      dpy(:,:) = c0

      do j = jlo, jhi
      do i = ilo, ihi

    !-------------------------------------------------------------------
    ! Estimate midpoint of backward trajectory relative to corner (i,j).
    !-------------------------------------------------------------------
         mpx = -p5 * dt * uvel(i,j)
         mpy = -p5 * dt * vvel(i,j)

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

         if (mpx >= c0 .and. mpy >= c0) then    ! cell (i+1,j+1)
            i2 = i+1
            j2 = j+1
            hindi(i,j) = i2
            hindj(i,j) = j2
            mat11(i,j) = msw(1,1,i2,j2)
            mat21(i,j) = msw(2,1,i2,j2)
            mat12(i,j) = msw(1,2,i2,j2)
            mat22(i,j) = msw(2,2,i2,j2)
            cxt(i,j)   = -p5*hts(i2,j2)
            cyt(i,j)   = -p5*htw(i2,j2)
            mpxs(i,j)  = -c1
            mpys(i,j)  = -c1
         elseif (mpx < c0 .and. mpy < c0) then  ! cell (i,j)
            i2 = i
            j2 = j
            hindi(i,j) = i2
            hindj(i,j) = j2
            mat11(i,j) = mne(1,1,i2,j2)
            mat21(i,j) = mne(2,1,i2,j2)
            mat12(i,j) = mne(1,2,i2,j2)
            mat22(i,j) = mne(2,2,i2,j2)
            cxt(i,j)   = p5*htn(i2,j2)
            cyt(i,j)   = p5*hte(i2,j2)
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
            cxt(i,j)   = -p5*htn(i2,j2)
            cyt(i,j)   =  p5*htw(i2,j2)
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
            cxt(i,j)   =  p5*hts(i2,j2)
            cyt(i,j)   = -p5*hte(i2,j2)
            mpxs(i,j)  =  c1
            mpys(i,j)  = -c1
         endif

    !-------------------------------------------------------------------
    ! Transform coordinates of the trajectory midpoint to the (i2,j2)
    ! reference frame.
    !-------------------------------------------------------------------

         mpxt(i,j) = cxt(i,j) + mat11(i,j)*mpx + mat12(i,j)*mpy
         mpyt(i,j) = cyt(i,j) + mat21(i,j)*mpx + mat22(i,j)*mpy

      enddo                     ! i
      enddo                     ! j

    !-------------------------------------------------------------------
    ! Transform (mpxt,mpyt) to a stretched coordinate system in which
    ! the coordinates of the four corners relative to the center are
    ! (-1,-1), (1,-1), (1,1), and (-1,1).
    !
    ! Iterate a couple of times for accuracy.
    ! (Occasionally abs(mpxs) or abs(mpys) > 1 after first iteration.)
    !-------------------------------------------------------------------

      do niter = 1, 3
         do j = jlo, jhi
         do i = ilo, ihi
            i2 = hindi(i,j)
            j2 = hindj(i,j)
            w1 = c2*mpxt(i,j) * dyt(i2,j2)
            w2 =   dxt(i2,j2) * dyt(i2,j2)
     &          + dyhx(i2,j2) * (c2*mpyt(i,j)
     &                         - mpxs(i,j)*mpys(i,j)*dxhy(i2,j2))
            w3 = c2*mpyt(i,j) * dxt(i2,j2)
            w4 =   dxt(i2,j2) * dyt(i2,j2)
     &          + dxhy(i2,j2) * (c2*mpxt(i,j)
     &                         - mpxs(i,j)*mpys(i,j)*dyhx(i2,j2))
            mpxs(i,j) = w1/w2
            mpys(i,j) = w3/w4
         enddo                  ! i
         enddo                  ! j
      enddo                     ! niter

    !-------------------------------------------------------------------
    ! Check for departure points out of bounds
    !-------------------------------------------------------------------
      dept_out_bound = .false.

      do j = jlo, jhi
      do i = ilo, ihi
         if (abs(mpxs(i,j)) > c1+puny .or. 
     &       abs(mpys(i,j)) > c1+puny) then
           dept_out_bound = .true.
         endif
      enddo
      enddo

      if (dept_out_bound) then
         do j = jlo, jhi
         do i = ilo, ihi
            if (abs(mpxs(i,j)) > c1+puny .or. 
     &          abs(mpys(i,j)) > c1+puny) then
               write (nu_diag,*) ''
               write (nu_diag,*)
     &              'Warning: Departure points out of bounds in remap'
               write (nu_diag,*) 'Should be in range (-1,1)'
               write (nu_diag,*) 'istep1, my_task, i, j =',
     &                            istep1, my_task, i, j
               write (nu_diag,*) 'x, y =', mpxs(i,j), mpys(i,j)
            endif
         enddo
         enddo
      endif

      do j = jlo, jhi
      do i = ilo, ihi

    !-------------------------------------------------------------------
    ! Transform the four corner velocities to the (i2,j2) reference frame.
    !-------------------------------------------------------------------

         i2   = hindi(i,j)
         j2   = hindj(i,j)

         u1t = msw(1,1,i2,j2)*uvel(i2-1,j2-1)
     &       + msw(1,2,i2,j2)*vvel(i2-1,j2-1)
         v1t = msw(2,1,i2,j2)*uvel(i2-1,j2-1)
     &       + msw(2,2,i2,j2)*vvel(i2-1,j2-1)

         u2t = mse(1,1,i2,j2)*uvel(i2,j2-1)
     &       + mse(1,2,i2,j2)*vvel(i2,j2-1)
         v2t = mse(2,1,i2,j2)*uvel(i2,j2-1)
     &       + mse(2,2,i2,j2)*vvel(i2,j2-1)

         u3t = mne(1,1,i2,j2)*uvel(i2,j2)
     &       + mne(1,2,i2,j2)*vvel(i2,j2)
         v3t = mne(2,1,i2,j2)*uvel(i2,j2)
     &       + mne(2,2,i2,j2)*vvel(i2,j2)

         u4t = mnw(1,1,i2,j2)*uvel(i2-1,j2)
     &       + mnw(1,2,i2,j2)*vvel(i2-1,j2)
         v4t = mnw(2,1,i2,j2)*uvel(i2-1,j2)
     &       + mnw(2,2,i2,j2)*vvel(i2-1,j2)

    !-------------------------------------------------------------------
    ! Using a bilinear approximation, estimate the velocity at the
    ! trajectory midpoint in the (i2,j2) reference frame.
    !-------------------------------------------------------------------

         umpt = p25 * ( u1t*(mpxs(i,j)-c1)*(mpys(i,j)-c1)
     &                - u2t*(mpxs(i,j)+c1)*(mpys(i,j)-c1)
     &                + u3t*(mpxs(i,j)+c1)*(mpys(i,j)+c1)
     &                - u4t*(mpxs(i,j)-c1)*(mpys(i,j)+c1) )

         vmpt = p25 * ( v1t*(mpxs(i,j)-c1)*(mpys(i,j)-c1)
     &                - v2t*(mpxs(i,j)+c1)*(mpys(i,j)-c1)
     &                + v3t*(mpxs(i,j)+c1)*(mpys(i,j)+c1)
     &                - v4t*(mpxs(i,j)-c1)*(mpys(i,j)+c1) )

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

      enddo                     ! i
      enddo                     ! j

      end subroutine departure_points

!=======================================================================
!
!BOP
!
! !IROUTINE: locate_triangles - compute areas and vertices of flux triangles
!
! !INTERFACE:
!
      subroutine locate_triangles (dpx, dpy,
     &        xp0_e, xp1_e, xp2_e, xp3_e, yp0_e, yp1_e, yp2_e, yp3_e, 
     &        xp0_n, xp1_n, xp2_n, xp3_n, yp0_n, yp1_n, yp2_n, yp3_n,
     &        triarea_e, iflux_e, jflux_e,
     &        triarea_n, iflux_n, jflux_n)
!
! !DESCRIPTION:
!
! Compute areas and vertices of flux triangles for east 
! and north cell edges.
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
      real (kind=dbl_kind), dimension (imt_local,jmt_local), 
     &     intent(in) ::
     &   dpx              ! x coordinates of departure points at cell corners
     &,  dpy              ! y coordinates of departure points at cell corners

      real (kind=dbl_kind), dimension (imt_local, jmt_local, ngroups),
     &     intent(out) ::
     &   triarea_e        ! area of east-edge flux triangle
     &,  triarea_n        ! area of north-edge flux triangle
     &,  xp0_e, yp0_e     ! coordinates of special triangle points
     &,  xp1_e, yp1_e     ! e for east edges, n for north edges
     &,  xp2_e, yp2_e   
     &,  xp3_e, yp3_e     
     &,  xp0_n, yp0_n   
     &,  xp1_n, yp1_n
     &,  xp2_n, yp2_n
     &,  xp3_n, yp3_n
 
      integer (kind=int_kind), intent(out),
     &     dimension (imt_local, jmt_local, ngroups) ::
     &   iflux_e          ! i index of cell contributing east flux
     &,  jflux_e          ! j index of cell contributing east flux
     &,  iflux_n          ! i index of cell contributing north flux
     &,  jflux_n          ! j index of cell contributing north flux
!
!EOP
!
      integer (kind=int_kind) ::
     &   i, j             ! horizontal indices of cell edge
     &,  i2, j2           ! horizontal indices of cell contributing flux
     &,  ng               ! triangle group index

      real (kind=dbl_kind) ::
     &   x1, y1, x2, y2   ! x,y coordinates of departure points from
                          ! east and north edges, as in DB
     &,  xa, ya, xb, yb   ! x,y coordinates of points where the lines joining
                          ! (x1,y1) and (x2,y2) cross cell edges, as in DB
     &,  xca, yca         ! transformed coordinates of corner point a
     &,  xda, yda         ! transformed coordinates of departure point a
     &,  xxa, yxa         ! transformed coordinates of intersection point xa
     &,  xya, yya         ! transformed coordinates of intersection point ya
     &,  xcb, ycb         ! transformed coordinates of corner point b
     &,  xdb, ydb         ! transformed coordinates of departure point b
     &,  xxb, yxb         ! transformed coordinates of intersection point xb
     &,  xyb, yyb         ! transformed coordinates of intersection point yb
     &,  xic, yic         ! transformed coordinates of point where the
                          ! line joining dep pts intersects the face
     &,  w1               ! work variable

      real (kind=dbl_kind), dimension (imt_local,jmt_local,ngroups) ::
     &   fluxsign         ! = 1. for positive flux, -1. for negative

      logical :: cnd1, cnd2, cnd3   ! conditionals

    !-------------------------------------------------------------------
    ! Triangle notation:
    ! For each edge, there are 20 triangles that can contribute a flux,
    ! but many of these are mutually exclusive.  It turns out that
    ! at most 5 triangles can contribute to flux integrals at one time.
    !
    ! For the east edge, these triangles are referred to as:
    ! (1) NE, NE1, N, N2
    ! (2) SE, SE1, S, S2
    ! (3) NE2, N1, SE2, S1
    ! (4) H1a, H1b, E1a, E2b
    ! (5) H2a, H2b, N2a, N2b
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
    ! East edges
    !-------------------------------------------------------------------

    !-------------------------------------------------------------------
    ! Initialize
    !-------------------------------------------------------------------

      do ng = 1, ngroups
         do j = 1, jmt_local
         do i = 1, imt_local
            iflux_e (i,j,ng) = 0
            jflux_e (i,j,ng) = 0
            xp0_e   (i,j,ng) = c0
            xp1_e   (i,j,ng) = c0
            xp2_e   (i,j,ng) = c0
            xp3_e   (i,j,ng) = c0
            yp0_e   (i,j,ng) = c0
            yp1_e   (i,j,ng) = c0
            yp2_e   (i,j,ng) = c0
            yp3_e   (i,j,ng) = c0
            fluxsign(i,j,ng) = c0
         enddo
         enddo
      enddo

    !-------------------------------------------------------------------
    ! Main loop
    !-------------------------------------------------------------------

      do j = jlo, jhi
      do i = ilo-1, ihi      ! includes W edge of cells with index ilo

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
    ! flux contribution from triangles in NE cell
    ! (All corner cells follow the same pattern.)
    !-------------------------------------------------------------------
         ! load horizontal indices of NE cell
         i2 = i+1
         j2 = j+1

         ! find vertex coordinates relative to center of NE cell

         xca = -p5*HTS(i2,j2)                         ! corner pt
         yca = -p5*HTW(i2,j2)

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
            xp1_e   (i,j,ng) = xca
            yp1_e   (i,j,ng) = yca
            xp2_e   (i,j,ng) = xxa
            yp2_e   (i,j,ng) = yxa
            xp3_e   (i,j,ng) = xda
            yp3_e   (i,j,ng) = yda
            iflux_e (i,j,ng) = i2
            jflux_e (i,j,ng) = j2
            fluxsign(i,j,ng) = -c1
         endif

         cnd2 = xa < c0 .and. y1 > c0 .and. x1 >= c0   ! NE1 (group 1)
         if (cnd2) then
            xp1_e   (i,j,ng) = xca
            yp1_e   (i,j,ng) = yca
            xp2_e   (i,j,ng) = xya
            yp2_e   (i,j,ng) = yya
            xp3_e   (i,j,ng) = xda
            yp3_e   (i,j,ng) = yda
            iflux_e (i,j,ng) = i2
            jflux_e (i,j,ng) = j2
            fluxsign(i,j,ng) = c1
         endif
            
         ! vertices of potential group 3 triangle
         ng = 3
         
         cnd3 = xa > c0 .and. y1 > c0 .and. x1 < c0    ! NE2 (group 3)
         if (cnd3) then
            xp1_e   (i,j,ng) = xca
            yp1_e   (i,j,ng) = yca
            xp2_e   (i,j,ng) = xxa
            yp2_e   (i,j,ng) = yxa
            xp3_e   (i,j,ng) = xya
            yp3_e   (i,j,ng) = yya
            iflux_e (i,j,ng) = i2
            jflux_e (i,j,ng) = j2
            fluxsign(i,j,ng) = -c1
         endif

    !-------------------------------------------------------------------
    ! flux contribution from triangles in N cell
    !-------------------------------------------------------------------
         i2 = i
         j2 = j+1

         xca =  p5*HTS(i2,j2) 
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
            xp1_e   (i,j,ng) = xca
            yp1_e   (i,j,ng) = yca
            xp2_e   (i,j,ng) = xxa
            yp2_e   (i,j,ng) = yxa
            xp3_e   (i,j,ng) = xda
            yp3_e   (i,j,ng) = yda
            iflux_e (i,j,ng) = i2
            jflux_e (i,j,ng) = j2
            fluxsign(i,j,ng) = c1
         endif

         cnd2 = xa > c0 .and. y1 > c0 .and. x1 < c0     ! N2 (group 1)
         if (cnd2) then
            xp1_e   (i,j,ng) = xca
            yp1_e   (i,j,ng) = yca
            xp2_e   (i,j,ng) = xya
            yp2_e   (i,j,ng) = yya
            xp3_e   (i,j,ng) = xda
            yp3_e   (i,j,ng) = yda
            iflux_e (i,j,ng) = i2
            jflux_e (i,j,ng) = j2
            fluxsign(i,j,ng) = -c1
         endif

         ng = 3
         
         cnd3 =  xa < c0 .and. y1 > c0 .and. x1 >= c0   ! N1 (group 3)
         if (cnd3) then
            xp1_e   (i,j,ng) = xca
            yp1_e   (i,j,ng) = yca
            xp2_e   (i,j,ng) = xxa
            yp2_e   (i,j,ng) = yxa
            xp3_e   (i,j,ng) = xya
            yp3_e   (i,j,ng) = yya
            iflux_e (i,j,ng) = i2
            jflux_e (i,j,ng) = j2
            fluxsign(i,j,ng) = c1
         endif

    !-------------------------------------------------------------------
    ! flux contribution from triangles in SE cell
    !-------------------------------------------------------------------
         i2 = i+1
         j2 = j-1

         xcb = -p5*HTN(i2,j2)
         ycb =  p5*HTW(i2,j2)

         xdb = xcb + mnw(1,1,i2,j2)*x2 + mnw(1,2,i2,j2)*y2
         ydb = ycb + mnw(2,1,i2,j2)*x2 + mnw(2,2,i2,j2)*y2

         xxb = xcb + mnw(1,1,i2,j2)*xb
         yxb = ycb + mnw(2,1,i2,j2)*xb

         xyb = xcb + mnw(1,2,i2,j2)*yb
         yyb = ycb + mnw(2,2,i2,j2)*yb

         ng = 2

         cnd1 = xb > c0 .and. y2 < c0 .and. x2 >= c0    ! SE (group 2)
         if (cnd1) then
            xp1_e   (i,j,ng) = xcb
            yp1_e   (i,j,ng) = ycb
            xp2_e   (i,j,ng) = xxb
            yp2_e   (i,j,ng) = yxb
            xp3_e   (i,j,ng) = xdb
            yp3_e   (i,j,ng) = ydb
            iflux_e (i,j,ng) = i2
            jflux_e (i,j,ng) = j2
            fluxsign(i,j,ng) = -c1
         endif

         cnd2 = xb < c0 .and. y2 < c0 .and. x2 >= c0    ! SE1 (group 2)
         if (cnd2) then
            xp1_e   (i,j,ng) = xcb
            yp1_e   (i,j,ng) = ycb
            xp2_e   (i,j,ng) = xyb
            yp2_e   (i,j,ng) = yyb
            xp3_e   (i,j,ng) = xdb
            yp3_e   (i,j,ng) = ydb
            iflux_e (i,j,ng) = i2
            jflux_e (i,j,ng) = j2
            fluxsign(i,j,ng) = c1
         endif

         ng = 3

         cnd3 = xb > c0 .and. y2 < c0 .and. x2 < c0     ! SE2 (group 3)
         if (cnd3) then
            xp1_e   (i,j,ng) = xcb
            yp1_e   (i,j,ng) = ycb
            xp2_e   (i,j,ng) = xxb
            yp2_e   (i,j,ng) = yxb
            xp3_e   (i,j,ng) = xyb
            yp3_e   (i,j,ng) = yyb
            iflux_e (i,j,ng) = i2
            jflux_e (i,j,ng) = j2
            fluxsign(i,j,ng) = -c1
         endif

    !-------------------------------------------------------------------
    ! flux contribution from triangles in S cell
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
            xp1_e   (i,j,ng) = xcb
            yp1_e   (i,j,ng) = ycb
            xp2_e   (i,j,ng) = xxb
            yp2_e   (i,j,ng) = yxb
            xp3_e   (i,j,ng) = xdb
            yp3_e   (i,j,ng) = ydb
            iflux_e (i,j,ng) = i2
            jflux_e (i,j,ng) = j2
            fluxsign(i,j,ng) = c1
         endif

         cnd2 = xb > c0 .and. y2 < c0 .and. x2 < c0     ! S2 (group 2)
         if (cnd2) then
            xp1_e   (i,j,ng) = xcb
            yp1_e   (i,j,ng) = ycb
            xp2_e   (i,j,ng) = xyb
            yp2_e   (i,j,ng) = yyb
            xp3_e   (i,j,ng) = xdb
            yp3_e   (i,j,ng) = ydb
            iflux_e (i,j,ng) = i2
            jflux_e (i,j,ng) = j2
            fluxsign(i,j,ng) = -c1
         endif

         ng = 3

         cnd3 = xb < c0 .and. y2 < c0 .and. x2 >= c0    ! S1 (group 3)
         if (cnd3) then
            xp1_e   (i,j,ng) = xcb
            yp1_e   (i,j,ng) = ycb
            xp2_e   (i,j,ng) = xxb
            yp2_e   (i,j,ng) = yxb
            xp3_e   (i,j,ng) = xyb
            yp3_e   (i,j,ng) = yyb
            iflux_e (i,j,ng) = i2
            jflux_e (i,j,ng) = j2
            fluxsign(i,j,ng) = c1
         endif

    !-------------------------------------------------------------------
    ! redefine departure points if not in home or north cells
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
    ! flux from triangles inside home cell
    ! Note that home and facing cells follow the same pattern.
    !-------------------------------------------------------------------

         ! load horizontal indices
         i2 = i
         j2 = j

         ! triangle vertices relative to center of home cell

         xca =  p5*HTN(i2,j2)
         yca =  p5*HTE(i2,j2)

         xcb =  p5*HTS(i2,j2)  
         ycb = -p5*HTE(i2,j2)

         xda = xca + mne(1,1,i2,j2)*x1 + mne(1,2,i2,j2)*y1
         yda = yca + mne(2,1,i2,j2)*x1 + mne(2,2,i2,j2)*y1

         xdb = xcb + mse(1,1,i2,j2)*x2 + mse(1,2,i2,j2)*y2
         ydb = ycb + mse(2,1,i2,j2)*x2 + mse(2,2,i2,j2)*y2

         xic = p5 * (w1*(HTN(i2,j2)-HTS(i2,j2)) + HTS(i2,j2))
         yic = (w1 - p5) * HTE(i2,j2)

         ! flux contribution from triangle that includes the
         ! E cell edge (part of convex quadrilateral inside home cell)

         ng = 4
         
         cnd1 = xa*xb >= c0 .and. xa+xb < c0            ! H1a (group 4)
         if (cnd1) then
            xp1_e   (i,j,ng) = xca
            yp1_e   (i,j,ng) = yca
            xp2_e   (i,j,ng) = xcb
            yp2_e   (i,j,ng) = ycb
            xp3_e   (i,j,ng) = xda
            yp3_e   (i,j,ng) = yda
            iflux_e (i,j,ng) = i2
            jflux_e (i,j,ng) = j2
            fluxsign(i,j,ng) = c1
         endif

        ! flux contribution from triangle lying along the upper part 
        ! of E edge for case of line xa-xb intersecting the edge

         cnd2 = xa*xb < c0 .and. x1 < c0                ! H1b (group 4)
         if (cnd2) then
            xp1_e   (i,j,ng) = xca
            yp1_e   (i,j,ng) = yca
            xp2_e   (i,j,ng) = xic
            yp2_e   (i,j,ng) = yic
            xp3_e   (i,j,ng) = xda
            yp3_e   (i,j,ng) = yda
            iflux_e (i,j,ng) = i2
            jflux_e (i,j,ng) = j2
            fluxsign(i,j,ng) = c1
         endif

         ! flux contribution from triangle touching but not
         ! lying along the E edge (other part of convex quadrilateral)    

         ng = 5

         cnd1 = xa*xb >= c0 .and. xa+xb < c0            ! H2a (group 5)
         if (cnd1) then
            xp1_e   (i,j,ng) = xcb
            yp1_e   (i,j,ng) = ycb
            xp2_e   (i,j,ng) = xda
            yp2_e   (i,j,ng) = yda
            xp3_e   (i,j,ng) = xdb
            yp3_e   (i,j,ng) = ydb
            iflux_e (i,j,ng) = i2
            jflux_e (i,j,ng) = j2
            fluxsign(i,j,ng) = c1
         endif

         ! flux contribution from triangle lying along the lower part 
         ! of E edge for case of line xa-xb intersecting the edge

         cnd2 = xa*xb < c0 .and. x2 < c0                ! H2b (group 5)
         if (cnd2) then
            xp1_e   (i,j,ng) = xcb
            yp1_e   (i,j,ng) = ycb
            xp2_e   (i,j,ng) = xic
            yp2_e   (i,j,ng) = yic
            xp3_e   (i,j,ng) = xdb
            yp3_e   (i,j,ng) = ydb
            iflux_e(i,j,ng)  = i2
            jflux_e(i,j,ng)  = j2
            fluxsign(i,j,ng) = c1
         endif
    
    !-------------------------------------------------------------------
    ! flux from triangles in E cell
    !-------------------------------------------------------------------

         i2 = i+1
         j2 = j

         xca = -p5*HTN(i2,j2)
         yca =  p5*HTW(i2,j2)

         xcb = -p5*HTS(i2,j2)                    
         ycb = -p5*HTW(i2,j2)

         xda = xca + mnw(1,1,i2,j2)*x1 + mnw(1,2,i2,j2)*y1
         yda = yca + mnw(2,1,i2,j2)*x1 + mnw(2,2,i2,j2)*y1 
 
         xdb = xcb + msw(1,1,i2,j2)*x2 + msw(1,2,i2,j2)*y2
         ydb = ycb + msw(2,1,i2,j2)*x2 + msw(2,2,i2,j2)*y2 

         xic = -p5 * (w1*(HTN(i2,j2)-HTS(i2,j2)) + HTS(i2,j2))
         yic = (w1 - p5) * HTW(i2,j2)

         ! flux contribution from triangle that includes the
         ! W cell edge (part of convex quadrilateral inside E cell)

         ng = 4

         cnd1 = xa*xb >= c0 .and. xa+xb > c0            ! E1a (group 4)
         if (cnd1) then
            xp1_e   (i,j,ng) = xca
            yp1_e   (i,j,ng) = yca
            xp2_e   (i,j,ng) = xcb
            yp2_e   (i,j,ng) = ycb
            xp3_e   (i,j,ng) = xda
            yp3_e   (i,j,ng) = yda
            iflux_e (i,j,ng) = i2
            jflux_e (i,j,ng) = j2
            fluxsign(i,j,ng) = -c1
         endif

         ! flux contribution from triangle lying along the upper part 
         ! of W edge for case of line xa-xb intersecting the edge

         cnd2 = xa*xb < c0 .and. x1 > c0                ! E1b (group 4)
         if (cnd2) then
            xp1_e   (i,j,ng) = xca
            yp1_e   (i,j,ng) = yca
            xp2_e   (i,j,ng) = xic
            yp2_e   (i,j,ng) = yic
            xp3_e   (i,j,ng) = xda
            yp3_e   (i,j,ng) = yda
            iflux_e (i,j,ng) = i2
            jflux_e (i,j,ng) = j2
            fluxsign(i,j,ng) = -c1
         endif

         ! flux contribution from triangle touching but not
         ! lying along the W edge (other part of convex quadrilateral)

         ng = 5

         cnd1 = xa*xb >= c0 .and. xa+xb > c0            ! E2a (group 5)
         if (cnd1) then
            xp1_e   (i,j,ng) = xcb
            yp1_e   (i,j,ng) = ycb
            xp2_e   (i,j,ng) = xda
            yp2_e   (i,j,ng) = yda
            xp3_e   (i,j,ng) = xdb
            yp3_e   (i,j,ng) = ydb
            iflux_e (i,j,ng) = i2
            jflux_e (i,j,ng) = j2
            fluxsign(i,j,ng) = -c1
         endif

         ! flux contribution from triangle lying along the lower part 
         ! of W edge for case of line xa-xb intersecting the edge

         cnd2 = xa*xb < c0 .and. x2 > c0                ! E2b (group 5)
         if (cnd2) then
            xp1_e   (i,j,ng) = xcb
            yp1_e   (i,j,ng) = ycb
            xp2_e   (i,j,ng) = xic
            yp2_e   (i,j,ng) = yic
            xp3_e   (i,j,ng) = xdb
            yp3_e   (i,j,ng) = ydb
            iflux_e (i,j,ng) = i2
            jflux_e (i,j,ng) = j2
            fluxsign(i,j,ng) = -c1
         endif

      enddo                     ! i
      enddo                     ! j

    !-------------------------------------------------------------------
    ! compute triangle areas with appropriate sign
    !-------------------------------------------------------------------

      do ng = 1, ngroups
         do j = 1, jmt_local
         do i = 1, imt_local

            w1 = p5 * abs( (xp2_e(i,j,ng)-xp1_e(i,j,ng)) *
     &                     (yp3_e(i,j,ng)-yp1_e(i,j,ng))
     &                   - (yp2_e(i,j,ng)-yp1_e(i,j,ng)) *
     &                     (xp3_e(i,j,ng)-xp1_e(i,j,ng)) )

            triarea_e(i,j,ng) = fluxsign(i,j,ng) * w1

         enddo                  ! i
         enddo                  ! j
      enddo                     ! ng

    !-------------------------------------------------------------------
    ! North edges
    ! (Follow the same pattern as east edges)
    !-------------------------------------------------------------------

    !-------------------------------------------------------------------
    ! Initialize
    !-------------------------------------------------------------------

      do ng = 1, ngroups
         do j = 1, jmt_local
         do i = 1, imt_local
            iflux_n (i,j,ng) = 0
            jflux_n (i,j,ng) = 0
            xp0_n   (i,j,ng) = c0
            xp1_n   (i,j,ng) = c0
            xp2_n   (i,j,ng) = c0
            xp3_n   (i,j,ng) = c0
            yp0_n   (i,j,ng) = c0
            yp1_n   (i,j,ng) = c0
            yp2_n   (i,j,ng) = c0
            yp3_n   (i,j,ng) = c0
            fluxsign(i,j,ng) = c0
         enddo
         enddo
      enddo

    !-------------------------------------------------------------------
    ! Main loop
    !-------------------------------------------------------------------

      do j = jlo-1, jhi  ! includes S edge of cells with index jlo
      do i = ilo, ihi

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
    ! flux contribution from triangles in NW cell
    !-------------------------------------------------------------------
         i2 = i-1
         j2 = j+1

         xca =  p5*HTS(i2,j2)                           ! corner pt
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
            xp1_n   (i,j,ng) = xca
            yp1_n   (i,j,ng) = yca
            xp2_n   (i,j,ng) = xya
            yp2_n   (i,j,ng) = yya
            xp3_n   (i,j,ng) = xda
            yp3_n   (i,j,ng) = yda
            iflux_n (i,j,ng) = i2
            jflux_n (i,j,ng) = j2
            fluxsign(i,j,ng) = -c1
         endif

         cnd2 = ya < c0 .and. x1 < c0 .and. y1 >= c0    ! NW1 (group 1)
         if (cnd2) then
            xp1_n   (i,j,ng) = xca
            yp1_n   (i,j,ng) = yca
            xp2_n   (i,j,ng) = xxa
            yp2_n   (i,j,ng) = yxa
            xp3_n   (i,j,ng) = xda
            yp3_n   (i,j,ng) = yda
            iflux_n (i,j,ng) = i2
            jflux_n (i,j,ng) = j2
            fluxsign(i,j,ng) = c1
         endif

         ng = 3
         
         cnd3 = ya > c0 .and. x1 < c0 .and. y1 < c0     ! NW2 (group 3)
         if (cnd3) then 
            xp1_n   (i,j,ng) = xca
            yp1_n   (i,j,ng) = yca
            xp2_n   (i,j,ng) = xya
            yp2_n   (i,j,ng) = yya
            xp3_n   (i,j,ng) = xxa
            yp3_n   (i,j,ng) = yxa
            iflux_n (i,j,ng) = i2
            jflux_n (i,j,ng) = j2
            fluxsign(i,j,ng) = -c1
         endif

    !-------------------------------------------------------------------
    ! flux contribution from triangles in W cell
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
            xp1_n   (i,j,ng) = xca
            yp1_n   (i,j,ng) = yca
            xp2_n   (i,j,ng) = xya
            yp2_n   (i,j,ng) = yya
            xp3_n   (i,j,ng) = xda
            yp3_n   (i,j,ng) = yda
            iflux_n (i,j,ng) = i2
            jflux_n (i,j,ng) = j2
            fluxsign(i,j,ng) = c1
         endif

         cnd2 = ya > c0 .and. x1 < c0 .and. y1 < c0     ! W2 (group 1)
         if (cnd2) then
            xp1_n   (i,j,ng) = xca
            yp1_n   (i,j,ng) = yca
            xp2_n   (i,j,ng) = xxa
            yp2_n   (i,j,ng) = yxa
            xp3_n   (i,j,ng) = xda
            yp3_n   (i,j,ng) = yda
            iflux_n (i,j,ng) = i2
            jflux_n (i,j,ng) = j2
            fluxsign(i,j,ng) = -c1
         endif

         ng = 3

         cnd3 = ya < c0 .and. x1 < c0 .and. y1 >= c0    ! W1 (group 3)
         if (cnd3) then  
            xp1_n   (i,j,ng) = xca
            yp1_n   (i,j,ng) = yca
            xp2_n   (i,j,ng) = xya
            yp2_n   (i,j,ng) = yya
            xp3_n   (i,j,ng) = xxa
            yp3_n   (i,j,ng) = yxa
            iflux_n (i,j,ng) = i2
            jflux_n (i,j,ng) = j2
            fluxsign(i,j,ng) = c1
         endif

    !-------------------------------------------------------------------
    ! flux contribution from triangles in NE cell
    !-------------------------------------------------------------------
         i2 = i+1
         j2 = j+1

         xcb = -p5*HTS(i2,j2) 
         ycb = -p5*HTW(i2,j2)
         
         xdb = xcb + msw(1,1,i2,j2)*x2 + msw(1,2,i2,j2)*y2
         ydb = ycb + msw(2,1,i2,j2)*x2 + msw(2,2,i2,j2)*y2

         xyb = xcb + msw(1,2,i2,j2)*yb
         yyb = ycb + msw(2,2,i2,j2)*yb

         xxb = xcb + msw(1,1,i2,j2)*xb   
         yxb = ycb + msw(2,1,i2,j2)*xb

         ng = 2

         cnd1 = yb > c0 .and. x2 > c0 .and. y2 >= c0    ! NE (group 2)
         if (cnd1) then
            xp1_n   (i,j,ng) = xcb
            yp1_n   (i,j,ng) = ycb
            xp2_n   (i,j,ng) = xyb
            yp2_n   (i,j,ng) = yyb
            xp3_n   (i,j,ng) = xdb
            yp3_n   (i,j,ng) = ydb
            iflux_n (i,j,ng) = i2
            jflux_n (i,j,ng) = j2
            fluxsign(i,j,ng) = -c1
         endif

         cnd2 = yb < c0 .and. x2 > c0  .and. y2 >= c0   ! NE1 (group 2)
         if (cnd2) then
            xp1_n   (i,j,ng) = xcb
            yp1_n   (i,j,ng) = ycb
            xp2_n   (i,j,ng) = xxb
            yp2_n   (i,j,ng) = yxb
            xp3_n   (i,j,ng) = xdb
            yp3_n   (i,j,ng) = ydb
            iflux_n (i,j,ng) = i2
            jflux_n (i,j,ng) = j2
            fluxsign(i,j,ng) = c1
         endif
            
         ng = 3

         cnd3 = yb > c0 .and. x2 > c0 .and. y2 < c0     ! NE2 (group 3)
         if (cnd3) then
            xp1_n   (i,j,ng) = xcb
            yp1_n   (i,j,ng) = ycb
            xp2_n   (i,j,ng) = xyb
            yp2_n   (i,j,ng) = yyb
            xp3_n   (i,j,ng) = xxb
            yp3_n   (i,j,ng) = yxb
            iflux_n (i,j,ng) = i2
            jflux_n (i,j,ng) = j2
            fluxsign(i,j,ng) = -c1
         endif

    !-------------------------------------------------------------------
    ! flux contribution from triangles in E cell
    !-------------------------------------------------------------------
         i2 = i+1
         j2 = j

         xcb = -p5*HTN(i2,j2)
         ycb =  p5*HTW(i2,j2)

         xdb = xcb + mnw(1,1,i2,j2)*x2 + mnw(1,2,i2,j2)*y2
         ydb = ycb + mnw(2,1,i2,j2)*x2 + mnw(2,2,i2,j2)*y2

         xyb = xcb + mnw(1,2,i2,j2)*yb
         yyb = ycb + mnw(2,2,i2,j2)*yb

         xxb = xcb + mnw(1,1,i2,j2)*xb
         yxb = ycb + mnw(2,1,i2,j2)*xb

         ng = 2

         cnd1 = yb < c0 .and. x2 > c0 .and. y2 < c0     ! E (group 2)
         if (cnd1) then
            xp1_n   (i,j,ng) = xcb
            yp1_n   (i,j,ng) = ycb
            xp2_n   (i,j,ng) = xyb
            yp2_n   (i,j,ng) = yyb
            xp3_n   (i,j,ng) = xdb
            yp3_n   (i,j,ng) = ydb
            iflux_n (i,j,ng) = i2
            jflux_n (i,j,ng) = j2
            fluxsign(i,j,ng) = c1
         endif

         cnd2 = yb > c0 .and. x2 > c0 .and. y2 < c0     ! E2 (group 2)
         if (cnd2) then
            xp1_n   (i,j,ng) = xcb
            yp1_n   (i,j,ng) = ycb
            xp2_n   (i,j,ng) = xxb
            yp2_n   (i,j,ng) = yxb
            xp3_n   (i,j,ng) = xdb
            yp3_n   (i,j,ng) = ydb
            iflux_n (i,j,ng) = i2
            jflux_n (i,j,ng) = j2
            fluxsign(i,j,ng) = -c1
         endif

         ng = 3
            
         cnd3 = yb < c0 .and. x2 > c0 .and. y2 >= c0    ! E1 (group 3)
         if (cnd3) then                   
            xp1_n   (i,j,ng) = xcb
            yp1_n   (i,j,ng) = ycb
            xp2_n   (i,j,ng) = xyb
            yp2_n   (i,j,ng) = yyb
            xp3_n   (i,j,ng) = xxb
            yp3_n   (i,j,ng) = yxb
            iflux_n (i,j,ng) = i2
            jflux_n (i,j,ng) = j2
            fluxsign(i,j,ng) = c1
         endif

    !-------------------------------------------------------------------
    ! redefine departure points if not in home or north cells
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
    ! flux from triangles inside home cell
    !-------------------------------------------------------------------
         i2 = i
         j2 = j

         xca = -p5*HTN(i2,j2)                 
         yca =  p5*HTW(i2,j2)

         xcb =  p5*HTN(i2,j2)                 
         ycb =  p5*HTE(i2,j2)

         xda = xca + mnw(1,1,i2,j2)*x1 + mnw(1,2,i2,j2)*y1
         yda = yca + mnw(2,1,i2,j2)*x1 + mnw(2,2,i2,j2)*y1

         xdb = xcb + mne(1,1,i2,j2)*x2 + mne(1,2,i2,j2)*y2
         ydb = ycb + mne(2,1,i2,j2)*x2 + mne(2,2,i2,j2)*y2

         xic = (p5 - w1) * HTN(i2,j2)          ! intersection w/ N edge
         yic = p5 * (w1*(HTW(i2,j2)-HTE(i2,j2)) + HTE(i2,j2))
         
         ! flux contribution from triangle that includes the
         ! N cell edge (part of convex quadrilateral inside home cell)

         ng = 4

         cnd1 = ya*yb >= c0 .and. ya+yb < c0            ! H1a (group 4)
         if (cnd1) then
            xp1_n   (i,j,ng) = xca
            yp1_n   (i,j,ng) = yca
            xp2_n   (i,j,ng) = xcb
            yp2_n   (i,j,ng) = ycb
            xp3_n   (i,j,ng) = xda
            yp3_n   (i,j,ng) = yda
            iflux_n (i,j,ng) = i2
            jflux_n (i,j,ng) = j2
            fluxsign(i,j,ng) = c1
         endif
         
         ! flux contribution from triangle lying along the left part of 
         ! the N edge for case of line ya-yb intersecting the edge

         cnd2 = ya*yb < c0 .and. y1 < c0                ! H1b (group 4)
         if (cnd2) then
            xp1_n   (i,j,ng) = xca
            yp1_n   (i,j,ng) = yca
            xp2_n   (i,j,ng) = xic
            yp2_n   (i,j,ng) = yic
            xp3_n   (i,j,ng) = xda
            yp3_n   (i,j,ng) = yda
            iflux_n (i,j,ng) = i2
            jflux_n (i,j,ng) = j2
            fluxsign(i,j,ng) = c1
         endif

         ! flux contribution from triangle touching but not lying 
         ! along the N edge (other part of convex quadrilateral)

         ng = 5

         iflux_n(i,j,ng) = i2
         jflux_n(i,j,ng) = j2
            
         cnd1 = ya*yb >= c0 .and. ya+yb < c0            ! H2a (group 5)
         if (cnd1) then
            xp1_n   (i,j,ng) = xcb
            yp1_n   (i,j,ng) = ycb
            xp2_n   (i,j,ng) = xda
            yp2_n   (i,j,ng) = yda
            xp3_n   (i,j,ng) = xdb
            yp3_n   (i,j,ng) = ydb
            iflux_n (i,j,ng) = i2
            jflux_n (i,j,ng) = j2
            fluxsign(i,j,ng) = c1
         endif

         ! flux contribution from triangle lying along the right part 
         ! of the N edge for case of line ya-yb intersecting the edge

         cnd2 = ya*yb < c0 .and. y2 < c0                ! H2b (group 5)
         if (cnd2) then
            xp1_n   (i,j,ng) = xcb
            yp1_n   (i,j,ng) = ycb
            xp2_n   (i,j,ng) = xic
            yp2_n   (i,j,ng) = yic
            xp3_n   (i,j,ng) = xdb
            yp3_n   (i,j,ng) = ydb
            iflux_n (i,j,ng) = i2
            jflux_n (i,j,ng) = j2
            fluxsign(i,j,ng) = c1
         endif
            
    !-------------------------------------------------------------------
    ! flux from triangles in N cell
    !-------------------------------------------------------------------
         i2 = i
         j2 = j+1
         
         xca = -p5*HTS(i2,j2)
         yca = -p5*HTW(i2,j2)

         xcb =  p5*HTS(i2,j2)
         ycb = -p5*HTE(i2,j2)

         xda = xca + msw(1,1,i2,j2)*x1 + msw(1,2,i2,j2)*y1
         yda = yca + msw(2,1,i2,j2)*x1 + msw(2,2,i2,j2)*y1

         xdb = xcb + mse(1,1,i2,j2)*x2 + mse(1,2,i2,j2)*y2
         ydb = ycb + mse(2,1,i2,j2)*x2 + mse(2,2,i2,j2)*y2

         xic = (p5 - w1)*HTS(i2,j2)
         yic = -p5 * (w1*(HTW(i2,j2)-HTE(i2,j2)) + HTE(i2,j2)) 

         ! flux contribution from triangle that includes the
         ! S cell edge (part of convex quadrilateral inside home cell)

         ng = 4

         cnd1 = ya*yb >= c0 .and. ya+yb > c0            ! N1a (group 4)
         if (cnd1) then
            xp1_n   (i,j,ng) = xca
            yp1_n   (i,j,ng) = yca
            xp2_n   (i,j,ng) = xcb
            yp2_n   (i,j,ng) = ycb
            xp3_n   (i,j,ng) = xda
            yp3_n   (i,j,ng) = yda
            iflux_n (i,j,ng) = i2
            jflux_n (i,j,ng) = j2
            fluxsign(i,j,ng) = -c1
         endif

         ! flux contribution from triangle lying along the left part 
         ! of the S edge for case of line ya-yb intersecting the edge

         cnd2 = ya*yb < c0 .and. y1 > c0                ! N1b (group 4)
         if (cnd2) then
            xp1_n   (i,j,ng) = xca
            yp1_n   (i,j,ng) = yca
            xp2_n   (i,j,ng) = xic
            yp2_n   (i,j,ng) = yic
            xp3_n   (i,j,ng) = xda
            yp3_n   (i,j,ng) = yda
            iflux_n (i,j,ng) = i2
            jflux_n (i,j,ng) = j2
            fluxsign(i,j,ng) = -c1
         endif

         ! flux contribution from triangle touching but not
         ! lying along the S edge (other part of convex quadrilateral)

         ng = 5
 
         cnd1 = ya*yb >= c0 .and. ya+yb > c0            ! N2a (group 5)
         if (cnd1) then
            xp1_n   (i,j,ng) = xcb
            yp1_n   (i,j,ng) = ycb
            xp2_n   (i,j,ng) = xda
            yp2_n   (i,j,ng) = yda
            xp3_n   (i,j,ng) = xdb
            yp3_n   (i,j,ng) = ydb
            iflux_n (i,j,ng) = i2
            jflux_n (i,j,ng) = j2
            fluxsign(i,j,ng) = -c1
         endif

         ! flux contribution from triangle lying along the right part
         ! of the S edge for case of line ya-yb intersecting the edge

         cnd2 = ya*yb < c0 .and. y2 > c0                ! N2b (group 5)
         if (cnd2) then
            xp1_n   (i,j,ng) = xcb
            yp1_n   (i,j,ng) = ycb
            xp2_n   (i,j,ng) = xic
            yp2_n   (i,j,ng) = yic
            xp3_n   (i,j,ng) = xdb
            yp3_n   (i,j,ng) = ydb
            iflux_n (i,j,ng) = i2
            jflux_n (i,j,ng) = j2
            fluxsign(i,j,ng) = -c1
         endif

      enddo                     ! i loop
      enddo                     ! j loop

    !-------------------------------------------------------------------
    ! compute triangle areas with appropriate sign
    !-------------------------------------------------------------------

      do ng = 1, ngroups
         do j = 1, jmt_local
         do i = 1, imt_local

            w1 = p5 * abs( (xp2_n(i,j,ng)-xp1_n(i,j,ng)) *
     &                     (yp3_n(i,j,ng)-yp1_n(i,j,ng))
     &                   - (yp2_n(i,j,ng)-yp1_n(i,j,ng)) *
     &                     (xp3_n(i,j,ng)-xp1_n(i,j,ng)) )

            triarea_n(i,j,ng) = fluxsign(i,j,ng) * w1

         enddo                  ! i
         enddo                  ! j
      enddo                     ! ng

    !-------------------------------------------------------------------
    ! Zero out tiny triangle areas.
    !-------------------------------------------------------------------

      do ng = 1, ngroups
      do j = 1, jmt_local
      do i = 1, imt_local
         if (abs(triarea_e(i,j,ng)) < puny) triarea_e(i,j,ng) = c0
         if (abs(triarea_n(i,j,ng)) < puny) triarea_n(i,j,ng) = c0
      enddo
      enddo
      enddo

      end subroutine locate_triangles

!=======================================================================
!
!BOP
!
! !IROUTINE: triangle_coordinates - find coordinates of 4 points
!
! !INTERFACE:
!
      subroutine triangle_coordinates (triarea, 
     &                                 xp0, xp1, xp2, xp3,
     &                                 yp0, yp1, yp2, yp3)
!
! !DESCRIPTION:
!
! For each triangle, find the coordinates of the 4 points needed
! to compute integrals of cubic polynomials, using a formula from 
! A.H. Stroud, Approximate Calculation of Multiple Integrals, 
! Prentice-Hall, 1971.  (Section 8.8, formula 3.1.)  
! Quadratic functions can be integrated using 3-point formulas, but 
! it is more efficient to use a single formula for both quadratics 
! and cubics.
!
! The formula is as follows:
! I3 = integral of f(x,y)*dA  
!    = AR * [ -9/16 *  f(x0,y0) 
!           + 25/48 * (f(x1,y1) + f(x2,y2) + f(x3,y3))]
! where I3 is the integral of a polynomial of 3rd order or lower,
! AR is the area of the triangle,  (x0,y0) is the midpoint,
! and the other three points are located 2/5 of the way from 
! the midpoint to each of the three vertices.
!
! !REVISION HISTORY:
!
! author William H. Lipscomb, LANL
!
! !USES:
!
! !INPUT/OUTPUT PARAMETERS:
!

      real (kind=dbl_kind), intent(in),
     &     dimension (imt_local, jmt_local, ngroups) ::
     &     triarea              ! areas of flux triangles 

      real (kind=dbl_kind), intent(inout), 
     &     dimension (imt_local, jmt_local, ngroups) ::
     &     xp0, yp0             ! coordinates of triangle points
     &,    xp1, yp1
     &,    xp2, yp2
     &,    xp3, yp3
!
!EOP
!
      integer (kind=int_kind) :: 
     &     i, j                 ! horizontal indices 
     &,    ng                   ! triangle index

      do ng = 1, ngroups
      do j = 1, jmt_local
      do i = 1, imt_local

         if (abs(triarea(i,j,ng)) > puny) then

         ! coordinates of midpoint
            xp0(i,j,ng) = p333 
     &                  * (xp1(i,j,ng) + xp2(i,j,ng) + xp3(i,j,ng))
            yp0(i,j,ng) = p333 
     &                  * (yp1(i,j,ng) + yp2(i,j,ng) + yp3(i,j,ng))

         ! coordinates of the other 3 points needed for integrals

            xp1(i,j,ng) = p4*xp1(i,j,ng) + p6*xp0(i,j,ng)
            yp1(i,j,ng) = p4*yp1(i,j,ng) + p6*yp0(i,j,ng)

            xp2(i,j,ng) = p4*xp2(i,j,ng) + p6*xp0(i,j,ng)
            yp2(i,j,ng) = p4*yp2(i,j,ng) + p6*yp0(i,j,ng)

            xp3(i,j,ng) = p4*xp3(i,j,ng) + p6*xp0(i,j,ng)
            yp3(i,j,ng) = p4*yp3(i,j,ng) + p6*yp0(i,j,ng)

         endif                  ! triarea > puny

      enddo                     ! i
      enddo                     ! j 
      enddo                     ! ng 

      end subroutine triangle_coordinates

!=======================================================================
!
!BOP
!
! !IROUTINE: load_tracers - load tracer array and compute tracer dependency vectors
!
! !INTERFACE:
!
      subroutine load_tracers (trm, n)
!
! !DESCRIPTION:
!
! Load tracer array and compute tracer dependency vectors.
! Default version assumes that the advected tracers are hice, hsno, 
!    Tsfc, qice(1:nilyr), and qsno.
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
! !INPUT/OUTPUT PARAMETERS:
!
      integer (kind=int_kind), intent(in) ::
     &     n              ! ice category index

      real (kind=dbl_kind), dimension (imt_local, jmt_local, ntracer),
     &     intent(out) ::
     &     trm            ! mean tracer values in each grid cell
!
!EOP
!
      integer (kind=int_kind) ::
     &     i, j           ! horizontal indices
     &,    k              ! ice layer index
     &,    nt             ! tracer index

      real (kind=dbl_kind) ::
     &     w1             ! work variables

      integer (kind=int_kind), dimension(imt_local*jmt_local) ::
     &     indxi          ! compressed i/j indices
     &,    indxj

      integer (kind=int_kind) ::
     &     icells         ! number of cells with ice
     &,    ij             ! combined i/j index

      icells = 0
      do j = 1, jmt_local
      do i = 1, imt_local
         if (aicen(i,j,n) > puny) then
            icells = icells + 1
            indxi(icells) = i
            indxj(icells) = j
         endif                  ! aicen > puny
      enddo
      enddo

      trm = c0

    !-------------------------------------------------------------------
    ! Load tracer array
    ! Note: If aice > 0, then hice > 0, but we can have hsno = 0.
    ! Alse note: We transport qice*nilyr rather than qice, so as to
    !  avoid extra operations here and in unload_tracer.
    !-------------------------------------------------------------------

      do ij = 1, icells
         i = indxi(ij)
         j = indxj(ij)
         w1 = c1 / aicen(i,j,n)
         work_l1(i,j) = c1 / vicen(i,j,n)
         trm(i,j,1) = vicen(i,j,n) * w1 ! hice
         trm(i,j,2) = vsnon(i,j,n) * w1 ! hsno
         trm(i,j,3) = Tsfcn(i,j,n)      ! Tsfc
      enddo

      do k = 1, nilyr
         do ij = 1, icells
           i = indxi(ij)
           j = indxj(ij)
           trm(i,j,3+k) = eicen(i,j,(n-1)*nilyr+k) * work_l1(i,j)  ! qice
         enddo      !ij
      enddo         !k

cdir$ ivdep      !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
      do ij = 1, icells
         i = indxi(ij)
         j = indxj(ij)
         if (trm(i,j,2) > puny)                            ! hsno > puny
     &        trm(i,j,4+nilyr) = esnon(i,j,n)/vsnon(i,j,n) ! qsno
      enddo

      end subroutine load_tracers

!=======================================================================
!
!BOP
!
! !IROUTINE: make_mask - make area and tracer masks
!
! !INTERFACE:
!
      subroutine make_masks (aim, aimask, trm, trmask)
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
      real (kind=dbl_kind), dimension (imt_local, jmt_local),
     &     intent(inout) ::
     &     aim            ! mean ice area in each grid cell

      real (kind=dbl_kind), dimension (imt_local, jmt_local),
     &     intent(out) ::
     &     aimask         ! = 1. if ice is present, else = 0.

      real (kind=dbl_kind), dimension (imt_local, jmt_local, ntracer),
     &     intent(inout), optional ::
     &     trm            ! mean tracer values in each grid cell

      real (kind=dbl_kind), dimension (imt_local, jmt_local, ntracer),
     &     intent(out), optional ::
     &     trmask         ! = 1. if tracer is present, else = 0.
!
!EOP
!
      integer (kind=int_kind) ::
     &     i, j           ! horizontal indices
     &,    nt             ! tracer index

    !-------------------------------------------------------------------
    ! ice area mask
    !-------------------------------------------------------------------

      do j = 1, jmt_local
      do i = 1, imt_local
         if (aim(i,j) > puny) then
            aimask(i,j) = c1
         else
            aimask(i,j) = c0
         endif
      enddo
      enddo

    !-------------------------------------------------------------------
    ! tracer masks
    !-------------------------------------------------------------------

      if (present(trm)) then

         do nt = 1, ntracer
            if (has_dependents(nt)) then
               do j = 1, jmt_local
               do i = 1, imt_local
                  if (abs(trm(i,j,nt)) > puny) then
                     trmask(i,j,nt) = c1
                  else
                     trmask(i,j,nt) = c0
                  endif
               enddo
               enddo
            else                ! mask is not needed, so set to zero
               do j = 1, jmt_local
               do i = 1, imt_local
                  trmask(i,j,nt) = c0
               enddo
               enddo
            endif
         enddo

      endif                     ! present(trm)

      end subroutine make_masks

!=======================================================================
!BOP
!
! !IROUTINE: unload_tracers - convert tracer array back to state variables
!
! !INTERFACE:
!
      subroutine unload_tracers (trm, n)
!
! !DESCRIPTION:
!
! Convert from tracer array back to state variables.
! Like load_tracers, this subroutine must be modified if a different 
! set of tracers is to be transported.  
!
! !REVISION HISTORY:
!
! author William H. Lipscomb, LANL
!
! !USES:
!
! !INPUT/OUTPUT PARAMETERS:
!
      integer (kind=int_kind), intent(in) ::
     &     n              ! ice category index

      real (kind=dbl_kind), dimension (imt_local, jmt_local, ntracer),
     &     intent(in) ::
     &     trm            ! mean tracer values in each grid cell
!
!EOP
!
      integer (kind=int_kind) ::
     &     i, j           ! horizontal indices
     &,    k              ! ice layer index

    !-------------------------------------------------------------------
    ! Compute state variables.
    !-------------------------------------------------------------------

      do j = jlo,jhi
      do i = ilo,ihi
         vicen(i,j,n) = aicen(i,j,n)*trm(i,j,1)       ! aice*hice
         vsnon(i,j,n) = aicen(i,j,n)*trm(i,j,2)       ! aice*hsno
         Tsfcn(i,j,n) = trm(i,j,3)                    ! Tsfc
         esnon(i,j,n) = vsnon(i,j,n)*trm(i,j,4+nilyr) ! vsno*qsno
      enddo                     ! i
      enddo                     ! j

      do k = 1, nilyr           
         do j = jlo,jhi
         do i = ilo,ihi                               ! vice*qice  
            eicen(i,j,ilyr1(n)+k-1) = vicen(i,j,n)*trm(i,j,3+k)
         enddo                  ! i
         enddo                  ! j
      enddo                     ! k

      end subroutine unload_tracers

!=======================================================================
!
!BOP
!
! !IROUTINE: construct_fields - construct fields of ice area and tracers
!
! !INTERFACE:
!
      subroutine construct_fields
     &        (aim,   aic,    aix,    aiy,   aimask,
     &         trm,   trc,    trx,    try,   trmask)
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
      real (kind=dbl_kind), dimension (imt_local, jmt_local),
     &   intent(in) ::
     &   aim              ! mean ice area
     &,  aimask           ! = 1. if ice is present, = 0. otherwise

      real (kind=dbl_kind), dimension (imt_local, jmt_local, ntracer),
     &   intent(in), optional ::
     &   trm              ! mean tracer
     &,  trmask           ! = 1. if tracer is present, = 0. otherwise

      real (kind=dbl_kind), dimension (imt_local, jmt_local),
     &   intent(out) ::
     &   aic              ! ice area at geometric center of cell
     &,  aix, aiy         ! limited derivative of ice area wrt x and y

      real (kind=dbl_kind), dimension (imt_local, jmt_local, ntracer),
     &   intent(out), optional ::
     &   trc              ! tracer at geometric center of cell
     &,  trx, try         ! limited derivative of tracer wrt x and y
!
!EOP
!
      integer (kind=int_kind) ::
     &   i, j             ! horizontal indices
     &,  nt, nt1          ! tracer indices

      real (kind=dbl_kind), dimension (imt_local, jmt_local) ::
     &   axav, ayav       ! x,y coordinates of center of ice area

      real (kind=dbl_kind), dimension (imt_local, jmt_local, ntracer) ::
     &   atxav, atyav     ! x,y coordinates of center of area*tracer

      real (kind=dbl_kind) ::
     &   w1, w2, w3, w4, w5, w6, w7   ! work variables

      integer (kind=int_kind), dimension(1:(ihi-ilo+1)*(jhi-jlo+1)) ::
     &   indxi            ! compressed indices in i/j directions         
     &,  indxj

      integer (kind=int_kind) ::
     &   icells            ! number of cells with ice present
     &,  ij               ! combined i/j horizontal index

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
    ! (And similarly for hsno(i,j) and Tsfcn(i,j,n).)
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

      call limited_gradient (aix, aiy, aim, 
     &                       xav, yav, hm)

      ! ice area

      icells = 0
      do j = jlo,jhi
      do i = ilo,ihi
         if (aimask(i,j) > puny) then ! ice is present
            icells = icells + 1
            indxi(icells) = i
            indxj(icells) = j
         else
            aic(i,j)  = c0
            axav(i,j) = c0
            ayav(i,j) = c0
        endif      ! aimask > puny
      enddo
      enddo

      do ij = 1,icells   ! ice is present
         i = indxi(ij)
         j = indxj(ij)

         ! ice area at geometric center
         aic(i,j) = aim(i,j) - xav(i,j)*aix(i,j) 
     &                       - yav(i,j)*aiy(i,j)
         ! center of ice area (axav,ayav) for each cell
         axav(i,j) = (aix(i,j)*xxav(i,j) + aiy(i,j)*xyav(i,j) 
     &              + aic(i,j)*xav(i,j)) / aim(i,j)
         ayav(i,j) = (aix(i,j)*xyav(i,j) + aiy(i,j)*yyav(i,j) 
     &              + aic(i,j)*yav(i,j)) / aim(i,j)
      enddo                     ! ij

      ! tracers

      if (present(trm)) then

       do nt = 1, ntracer

         do j = jlo, jhi
         do i = ilo, ihi
            trc(i,j,nt) = c0
         enddo
         enddo

         if (depend(nt)==0) then ! independent of other tracers
                                 ! (surface temp, ice and snow thickness)
 
            call limited_gradient 
     &           (trx(:,:,nt), try(:,:,nt), trm(:,:,nt),
     &            axav,        ayav,        aimask)

            if (has_dependents(nt)) then   ! need center of area*tracer

               do j = jlo, jhi
               do i = ilo, ihi
                  atxav(i,j,nt) = c0
                  atyav(i,j,nt) = c0
               enddo
               enddo

cdir$ ivdep      !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
               do ij = 1, icells  ! Note: no trx or try in ghost cells
                                  ! (bound calls are below)
                  i = indxi(ij)
                  j = indxj(ij)

                  ! tracer value at geometric center
                  trc(i,j,nt) = trm(i,j,nt) - trx(i,j,nt)*axav(i,j) 
     &                                      - try(i,j,nt)*ayav(i,j)

                  if (trmask(i,j,nt) > puny) then 

                     ! center of area*tracer
                     w1 = aic(i,j)*trc(i,j,nt)
                     w2 = aic(i,j)*trx(i,j,nt) 
     &                  + aix(i,j)*trc(i,j,nt)
                     w3 = aic(i,j)*try(i,j,nt) 
     &                  + aiy(i,j)*trc(i,j,nt)
                     w4 = aix(i,j)*trx(i,j,nt)
                     w5 = aix(i,j)*try(i,j,nt) 
     &                  + aiy(i,j)*trx(i,j,nt)
                     w6 = aiy(i,j)*try(i,j,nt)
                     w7 = c1 / (aim(i,j)*trm(i,j,nt))
                     atxav(i,j,nt) = (w1*xav(i,j) + w2*xxav(i,j)  
     &                             + w3*xyav(i,j)  + w4*xxxav(i,j) 
     &                             + w5*xxyav(i,j) + w6*xyyav(i,j)) * w7
                     atyav(i,j,nt) = (w1*yav(i,j)   + w2*xyav(i,j)  
     &                             + w3*yyav(i,j)  + w4*xxyav(i,j) 
     &                             + w5*xyyav(i,j) + w6*yyyav(i,j)) * w7
                  endif         ! trmask

               enddo            ! ij

            else                ! no dependents

               do ij = 1, icells      ! ice is present
                  i = indxi(ij)
                  j = indxj(ij)

                  ! tracer value at geometric center
                  trc(i,j,nt) = trm(i,j,nt) - trx(i,j,nt)*axav(i,j) 
     &                                      - try(i,j,nt)*ayav(i,j)
               enddo            ! ij

            endif               ! has_dependents

         else                   ! tracer nt depends on tracer nt1
                                ! (ice and snow enthalpy)
            nt1 = depend(nt)

            call limited_gradient 
     &           (trx(:,:,nt),    try(:,:,nt),    trm(:,:,nt), 
     &            atxav(:,:,nt1), atyav(:,:,nt1), trmask(:,:,nt1) )

            do j = jlo, jhi
            do i = ilo, ihi

               if (trmask(i,j,nt1) > puny) then
                  trc(i,j,nt) = trm(i,j,nt) 
     &                        - trx(i,j,nt) * atxav(i,j,nt1) 
     &                        - try(i,j,nt) * atyav(i,j,nt1)
               endif            ! trmask(nt1)

            enddo               ! i
            enddo               ! j

         endif                  ! depend(nt) = 0
       enddo                    ! ntracer

      endif                     ! present (trm)

      end subroutine construct_fields

!=======================================================================
!
!BOP
!
! !IROUTINE: flux_integrals - compute the fluxes across each edge
!
! !INTERFACE:
!
      subroutine flux_integrals (triarea, iflux,   jflux,
     &                           xp0,     xp1,     xp2,   xp3,
     &                           yp0,     yp1,     yp2,   yp3,
     &                           aic,     aix,     aiy,   aiflx, 
     &                           trc,     trx,     try,   atflx)
!
! !DESCRIPTION:
!
! Compute the fluxes across each edge by integrating the ice area 
! and tracers over each flux triangle.
! Input variables have the same meanings as in the main subroutine.
! Repeated use of certain sums makes the calculation more efficient.
!
! !REVISION HISTORY:
!
! author William H. Lipscomb, LANL
!
! !USES:
!
! !INPUT/OUTPUT PARAMETERS:
!
      real (kind=dbl_kind), intent(in),
     &     dimension (imt_local, jmt_local, ngroups) ::
     &     triarea
     &,    xp0, yp0
     &,    xp1, yp1
     &,    xp2, yp2
     &,    xp3, yp3

      integer (kind=int_kind), intent(in),
     &     dimension (imt_local, jmt_local, ngroups) ::
     &     iflux
     &,    jflux

      real (kind=dbl_kind), intent(in),
     &     dimension (imt_local, jmt_local) :: 
     &     aic, aix, aiy

      real (kind=dbl_kind), intent(out),
     &     dimension (imt_local, jmt_local) :: 
     &     aiflx

      real (kind=dbl_kind), intent(in),
     &     dimension (imt_local, jmt_local, ntracer), optional ::  
     &     trc, trx, try

      real (kind=dbl_kind), intent(out),
     &     dimension (imt_local, jmt_local, ntracer), optional :: 
     &     atflx
!
!EOP
!
      integer (kind=int_kind) :: 
     &     i, j               ! horizontal indices of edge
     &,    i2, j2             ! horizontal indices of cell contributing flux
     &,    ng                 ! triangle index
     &,    nt, nt1            ! tracer indices 

      real (kind=dbl_kind) ::
     &     a0, a1, a2, a3           ! ice area at internal points
     &,    w0, w1, w2, w3           ! work variables

      real (kind=dbl_kind), dimension (imt_local, jmt_local) ::
     &     asum, axsum, aysum       ! sum of area, area*x, and area*y
     &,    axxsum, axysum, ayysum   ! sum of area*x*x, area*x*y, area*y*y

      real (kind=dbl_kind), dimension (imt_local, jmt_local, ntracer) ::
     &     atsum              ! sum of area*tracer
     &,    atxsum             ! sum of area*tracer*x
     &,    atysum             ! sum of area*tracer*y
 
      integer (kind=int_kind), dimension (1:(ihi-ilo)*(jhi-jlo)) ::
     &   indxi   ! compressed index in i-direction
     &,  indxj   ! compressed index in j-direction

      integer (kind=int_kind) ::
     &   icells  ! number of cells where triarea > puny
     &,  ij      ! combined i/j index

      ! initialize fluxes
      aiflx = c0
      if (present(atflx)) atflx = c0

    !-------------------------------------------------------------------
    ! Main loop
    !-------------------------------------------------------------------

      do ng = 1, ngroups

         icells = 0
         do j = jlo-1, jhi
         do i = ilo-1, ihi
            if (abs(triarea(i,j,ng)) > puny) then
               icells = icells+1
               indxi(icells) = i
               indxj(icells) = j
            endif
         enddo                  ! i
         enddo                  ! j

cdir$ ivdep      !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu

         do ij = 1, icells
            i = indxi(ij)
            j = indxj(ij)

            i2 = iflux(i,j,ng)
            j2 = jflux(i,j,ng)

            ! area fluxes  

            ! Weighting factors are incorporated into the ice
            ! area terms a0, a1, a2, and a3.
            a0 = p5625m * (aic(i2,j2) + xp0(i,j,ng)*aix(i2,j2) 
     &                                + yp0(i,j,ng)*aiy(i2,j2))
            a1 = p52083 * (aic(i2,j2) + xp1(i,j,ng)*aix(i2,j2) 
     &                                + yp1(i,j,ng)*aiy(i2,j2))   
            a2 = p52083 * (aic(i2,j2) + xp2(i,j,ng)*aix(i2,j2) 
     &                                + yp2(i,j,ng)*aiy(i2,j2))
            a3 = p52083 * (aic(i2,j2) + xp3(i,j,ng)*aix(i2,j2) 
     &                                + yp3(i,j,ng)*aiy(i2,j2)) 
            asum(i,j) = a0 + a1 + a2 + a3
            aiflx(i,j) = aiflx(i,j) + triarea(i,j,ng)*asum(i,j)

            ! quantities needed for volume fluxes
            w0 = a0 * xp0(i,j,ng)
            w1 = a1 * xp1(i,j,ng)
            w2 = a2 * xp2(i,j,ng)
            w3 = a3 * xp3(i,j,ng)
               
            axsum(i,j) = w0 + w1 + w2 + w3

            axxsum(i,j) = w0*xp0(i,j,ng) + w1*xp1(i,j,ng)
     &                  + w2*xp2(i,j,ng) + w3*xp3(i,j,ng)

            axysum(i,j) = w0*yp0(i,j,ng) + w1*yp1(i,j,ng)
     &                  + w2*yp2(i,j,ng) + w3*yp3(i,j,ng)

            w0 = a0 * yp0(i,j,ng)
            w1 = a1 * yp1(i,j,ng)
            w2 = a2 * yp2(i,j,ng)
            w3 = a3 * yp3(i,j,ng)

            aysum(i,j) = w0 + w1 + w2 + w3

            ayysum(i,j) = w0*yp0(i,j,ng) + w1*yp1(i,j,ng)
     &                  + w2*yp2(i,j,ng) + w3*yp3(i,j,ng)

         enddo                  ! ij

         ! area * tracer fluxes

         if (present(atflx)) then

          do nt = 1, ntracer
            if (depend(nt)==0) then ! does not depend on another tracer

cdir$ ivdep      !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
               do ij = 1, icells
                  i = indxi(ij)
                  j = indxj(ij)

                  i2 = iflux(i,j,ng)
                  j2 = jflux(i,j,ng)

                  atsum(i,j,nt) =  asum(i,j) * trc(i2,j2,nt) 
     &                          + axsum(i,j) * trx(i2,j2,nt) 
     &                          + aysum(i,j) * try(i2,j2,nt)

                  atflx(i,j,nt) = atflx(i,j,nt)
     &                          + triarea(i,j,ng) * atsum(i,j,nt)

                  ! quantities needed for dependent tracers

                  atxsum(i,j,nt) =  axsum(i,j) * trc(i2,j2,nt) 
     &                           + axxsum(i,j) * trx(i2,j2,nt) 
     &                           + axysum(i,j) * try(i2,j2,nt)

                  atysum(i,j,nt) =  aysum(i,j) * trc(i2,j2,nt) 
     &                           + axysum(i,j) * trx(i2,j2,nt) 
     &                           + ayysum(i,j) * try(i2,j2,nt)
               enddo            ! ij

            else                ! depends on another tracer
               nt1 = depend(nt)

cdir$ ivdep      !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
               do ij = 1, icells
                  i = indxi(ij)
                  j = indxj(ij)

                  i2 = iflux(i,j,ng)
                  j2 = jflux(i,j,ng)

                  atsum(i,j,nt) =  atsum(i,j,nt1) * trc(i2,j2,nt) 
     &                          + atxsum(i,j,nt1) * trx(i2,j2,nt) 
     &                          + atysum(i,j,nt1) * try(i2,j2,nt)

                  atflx(i,j,nt) = atflx(i,j,nt) 
     &                          + triarea(i,j,ng) * atsum(i,j,nt)
                  
               enddo            ! ij

            endif               ! depend(nt) = 0
         enddo                  ! ntracer
        endif                   ! present(atflx)
      enddo                     ! ng

      end subroutine flux_integrals

!=======================================================================
!
!BOP
!
! !IROUTINE: update_fields - compute new area and tracers
!
! !INTERFACE:
!
      subroutine update_fields (n, aiflxe, aiflxn, aim,
     &                             atflxe, atflxn, trm)
!
! !DESCRIPTION:
!
! Given fluxes through cell edges, compute new area and tracers.
!
! !REVISION HISTORY:
!
! author William H. Lipscomb, LANL
!
! !USES:
!
      use ice_exit
!
! !INPUT/OUTPUT PARAMETERS:
!
      real (kind=dbl_kind), dimension (imt_local, jmt_local),
     &   intent(in) ::
     &   aiflxe, aiflxn   ! flux of aice thru east and north cell edges
     
      real (kind=dbl_kind), dimension (imt_local, jmt_local),
     &   intent(inout) ::
     &   aim              ! mean ice area

      real (kind=dbl_kind), dimension (imt_local, jmt_local, ntracer),
     &   intent(in), optional ::
     &   atflxe, atflxn   ! flux of aice*tracer thru E and N cell edges

      real (kind=dbl_kind), dimension (imt_local, jmt_local, ntracer),
     &   intent(inout), optional ::
     &   trm              ! mean tracers

      integer (kind=int_kind), intent(in) ::
     &   n                ! ice category index (for error diagnostics)
!
!EOP
!
      integer (kind=int_kind) ::
     &   i, j             ! horizontal indices
     &,  nt, nt1          ! tracer indices

      real (kind=dbl_kind), dimension(imt_local,jmt_local,ntracer) ::
     &   atold            ! starting area*tracer

      real (kind=dbl_kind) ::
     &   w1, w2           ! work variables

      logical (kind=log_kind) ::
     &   aim_negative         ! flag for aim < -puny

      integer (kind=int_kind), dimension(1:(ihi-ilo+1)*(jhi-jlo+1)) ::
     &   indxi            ! compressed indices in i and j directions
     &,  indxj

      integer (kind=int_kind) ::
     &   icells           ! number of cells with aim > 0. 
     &,  ij               ! combined i/j horizontal index

    !-------------------------------------------------------------------
    ! Save starting values of area*tracer
    !-------------------------------------------------------------------

      if (present(trm)) then
         do nt = 1, ntracer
            if (depend(nt)==0) then ! does not depend on other tracers
               do j = jlo, jhi
               do i = ilo, ihi
                  atold(i,j,nt) = aim(i,j) * trm(i,j,nt)
               enddo            ! i
               enddo              ! j
            else                ! depends on another tracer
               nt1 = depend(nt)
               do j = jlo, jhi
               do i = ilo, ihi
                  atold(i,j,nt) = aim(i,j) * trm(i,j,nt1) * trm(i,j,nt)   
               enddo            ! i
               enddo            ! j
            endif               ! depend(nt) = 0
         enddo                  ! nt
      endif                     ! present(trm)

    !-------------------------------------------------------------------
    ! Update ice area
    !-------------------------------------------------------------------

      aim_negative = .false.        ! initialization

      do j = jlo, jhi
      do i = ilo, ihi

         w1 = aiflxe(i,j) - aiflxe(i-1,j)
     &      + aiflxn(i,j) - aiflxn(i,j-1)
         aim(i,j) = aim(i,j) - w1*tarear(i,j)

         if (aim(i,j) < -puny) then    ! there is a bug
            aim_negative = .true.
         elseif (aim(i,j) < c0) then   ! set to zero
            aim(i,j) = c0
         endif

      enddo
      enddo

      if (aim_negative) then
         do j = jlo, jhi
         do i = ilo, ihi
            if (aim(i,j) < -puny) then
               write (nu_diag,*) ''
               write (nu_diag,*) 'New area < 0, istep =', istep1
               write (nu_diag,*) '(my_task,i,j,n) =', my_task,i,j,n
               write (nu_diag,*) 'Old area =', aim(i,j) + w1*tarear(i,j)
               write (nu_diag,*) 'New area =', aim(i,j)
               write (nu_diag,*) 'Net flux =', -w1*tarear(i,j)
               call abort_ice('remap transport: negative area')
            endif
         enddo                  ! i
         enddo                  ! j
      endif                     ! aim_negative

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
      
         do nt = 1, ntracer

            trm(:,:,nt) = c0

            if (depend(nt)==0) then ! does not depend on other tracers

               do ij = 1, icells
                  i = indxi(ij)
                  j = indxj(ij)

                  w1  = atflxe(i,j,nt) - atflxe(i-1,j,nt)
     &                + atflxn(i,j,nt) - atflxn(i,j-1,nt)
                  trm(i,j,nt) = (atold(i,j,nt) - w1*tarear(i,j)) 
     &                          / aim(i,j)
               enddo            ! ij

            else                ! depends on another tracer
               nt1 = depend(nt)

cdir$ ivdep      !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
               do ij = 1, icells
                  i = indxi(ij)
                  j = indxj(ij)

                  if (abs(trm(i,j,nt1)) > puny) then
                     w1  = atflxe(i,j,nt) - atflxe(i-1,j,nt)
     &                   + atflxn(i,j,nt) - atflxn(i,j-1,nt)
                     trm(i,j,nt) = (atold(i,j,nt) - w1*tarear(i,j)) 
     &                           / (aim(i,j) * trm(i,j,nt1))
                  endif

               enddo            ! ij

            endif               ! depend(nt) = 0
         enddo                  ! nt
      endif                     ! present(trm)

      end subroutine update_fields

!=======================================================================
!
!BOP
!
! !IROUTINE: limited_gradient - compute limited gradient the scalar field phi
!
! !INTERFACE:
!
      subroutine limited_gradient (gx, gy, phi, cnx, cny, phimask)
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
      real (kind=dbl_kind), dimension (imt_local,jmt_local), 
     &     intent (in) ::
     &    phi     ! input tracer field (mean values in each grid cell)
     &,   cnx     ! x-coordinate of phi relative to geometric center of cell
     &,   cny     ! y-coordinate of phi relative to geometric center of cell
     &,   phimask
          ! phimask(i,j) = 1 if phi(i,j) has physical meaning, = 0 otherwise.
          ! For instance, aice has no physical meaning on land points,
          ! and hice no physical meaning where aice = 0.

      real (kind=dbl_kind), intent(out), 
     &   dimension (imt_local,jmt_local) :: 
     &    gx      ! limited x-direction gradient
     &,   gy      ! limited y-direction gradient
!
!EOP
!
      integer (kind=int_kind) :: i, j

      real (kind=dbl_kind) ::
     &    phi_nw, phi_n, phi_ne   ! values of phi in 8 neighbor cells 
     &,   phi_w,         phi_e
     &,   phi_sw, phi_s, phi_se
     &,   qmn, qmx       ! min and max value of phi within grid cell
     &,   pmn, pmx       ! min and max value of phi among neighbor cells
     &,   w1, w2, w3, w4 ! work variables

      integer (kind=int_kind), dimension(1:(ihi-ilo+1)*(jhi-jlo+1)) ::
     &    indxi          ! compressed indices for i/j direction 
     &,   indxj

      integer (kind=int_kind) ::
     &    icells         ! number of unmasked grid cells
     &,   ij             ! combined i/j horizontal index

      real (kind=dbl_kind) ::
     &    gxtmp, gytmp   ! temporary term for x- and y- limited gradient

      icells = 0
      do j = jlo,jhi
      do i = ilo,ihi
         if (phimask(i,j) > puny) then
            icells = icells + 1
            indxi(icells) = i
            indxj(icells) = j
         endif                  ! phimask > puny
         gx(i,j) = c0           ! initialization
         gy(i,j) = c0
      enddo
      enddo

      do ij = 1, icells
         i = indxi(ij)
         j = indxj(ij)

         ! Store values of phi in the 8 neighbor cells.
         ! Note: phimask = 1. or 0.  If phimask = 1., use the true value;
         !  if phimask = 0., use the home cell value so that non-physical 
         !  values of phi do not contribute to the gradient.
         phi_nw = phimask(i-1,j+1) * phi(i-1,j+1) 
     &      + (c1-phimask(i-1,j+1))* phi(i,j)
         phi_n  = phimask(i,j+1)   * phi(i,j+1) 
     &      + (c1-phimask(i,j+1))  * phi(i,j)
         phi_ne = phimask(i+1,j+1) * phi(i+1,j+1) 
     &      + (c1-phimask(i+1,j+1))* phi(i,j)
         phi_w  = phimask(i-1,j)   * phi(i-1,j) 
     &      + (c1-phimask(i-1,j))  * phi(i,j)
         phi_e  = phimask(i+1,j)   * phi(i+1,j) 
     &      + (c1-phimask(i+1,j))  * phi(i,j)
         phi_sw = phimask(i-1,j-1) * phi(i-1,j-1) 
     &      + (c1-phimask(i-1,j-1))* phi(i,j)
         phi_s  = phimask(i,j-1)   * phi(i,j-1) 
     &      + (c1-phimask(i,j-1))  * phi(i,j)
         phi_se = phimask(i+1,j-1) * phi(i+1,j-1) 
     &      + (c1-phimask(i+1,j-1))* phi(i,j)
               
         ! unlimited gradient components
         ! grid lengths computed in init_remap for efficiency:
         !       dxt2r(i,j) = 1 / (dxt(i+1,j) + dxt(i,j))
         !       dyt2r(i,j) = 1 / (dyt(i,j+1) + dyt(i,j))
         gxtmp = (phi_e - phi(i,j)) * dxt2r(i,j) 
     &         + (phi(i,j) - phi_w) * dxt2r(i-1,j)
         gytmp = (phi_n - phi(i,j)) * dyt2r(i,j) 
     &         + (phi(i,j) - phi_s) * dyt2r(i,j-1)

         ! minimum and maximum among the nine local cells
         pmn = min (phi_nw, phi_n,  phi_ne, phi_w, phi(i,j),
     &              phi_e,  phi_sw, phi_s,  phi_se)
         pmx = max (phi_nw, phi_n,  phi_ne, phi_w, phi(i,j),
     &              phi_e,  phi_sw, phi_s,  phi_se)
            
         pmn = pmn - phi(i,j)
         pmx = pmx - phi(i,j)

         ! minimum and maximum deviation of phi within the cell
         w1  =  (p5*HTN(i,j) - cnx(i,j)) * gxtmp 
     &        + (p5*HTE(i,j) - cny(i,j)) * gytmp
         w2  =  (p5*HTS(i,j) - cnx(i,j)) * gxtmp
     &        - (p5*HTE(i,j) + cny(i,j)) * gytmp
         w3  = -(p5*HTS(i,j) + cnx(i,j)) * gxtmp
     &        - (p5*HTW(i,j) + cny(i,j)) * gytmp
         w4  =  (p5*HTW(i,j) - cny(i,j)) * gytmp
     &        - (p5*HTN(i,j) + cnx(i,j)) * gxtmp
         
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
!
!BOP
!
! !IROUTINE: conserved_sums - computes global sums
!
! !INTERFACE:
!
      subroutine conserved_sums (aim, asum, trm, atsum)
!
! !DESCRIPTION:
!
! Computes global sums of conserved variables over the physical grid
!
! !REVISION HISTORY:
!
! author William H. Lipscomb, LANL
!
! !USES:
!
      use ice_mpi_internal
!
! !INPUT/OUTPUT PARAMETERS:
!
      real (kind=dbl_kind), dimension (imt_local,jmt_local),
     &     intent(in) ::
     &     aim             ! mean ice area

      real (kind=dbl_kind), intent(out) ::
     &     asum            ! global sum of area

      real (kind=dbl_kind), dimension (imt_local,jmt_local,ntracer),
     &     intent(in), optional ::
     &     trm             ! mean tracer

      real (kind=dbl_kind), dimension(ntracer), intent(out), optional ::
     &     atsum           ! global sum of area*tracer
!
!EOP
!
      integer (kind=int_kind) ::
     &     i,j             ! horizontal indices
     &,    nt, nt1         ! tracer indices
     &,    nwork           ! number of points in work array
 
      real (kind=dbl_kind), dimension(ilo:ihi,jlo:jhi,0:ntracer) ::
     &     work            ! work array

      nwork = (ihi-ilo+1)*(jhi-jlo+1)

      do j = jlo, jhi
      do i = ilo, ihi
         work(i,j,0) = tarea(i,j)*aim(i,j)
      enddo
      enddo

      asum = ice_global_real_sum(nwork,work(:,:,0))

      if (present(trm)) then 
         do j = jlo, jhi
         do i = ilo, ihi
            do nt = 1, ntracer
               if (depend(nt)==0) then ! does not depend on another tracer
                  work(i,j,nt) = tarea(i,j)*aim(i,j)*trm(i,j,nt)
               else
                  nt1 = depend(nt)
                  work(i,j,nt) = 
     &                 tarea(i,j)*aim(i,j)*trm(i,j,nt1)*trm(i,j,nt)
               endif
            enddo               ! nt
         enddo                  ! i
         enddo                  ! j

         do nt = 1, ntracer
            atsum(nt) = ice_global_real_sum(nwork,work(:,:,nt))
         enddo
      endif                     ! present(trm)
      
      end subroutine conserved_sums

!=======================================================================
!
!BOP
!
! !IROUTINE: global_conservation - check for changes in conserved quantities
!
! !INTERFACE:
!
      subroutine global_conservation (asum_init,  asum_final, 
     &                                atsum_init, atsum_final)
!
! !DESCRIPTION:
!
! Check whether values of conserved quantities have changed
! An error probably means that ghost cells are treated incorrectly.
!
! !REVISION HISTORY:
!
! author William H. Lipscomb, LANL
!
! !USES:
!
      use ice_exit
!
! !INPUT/OUTPUT PARAMETERS:
!
      real (kind=dbl_kind), intent(in) ::
     &     asum_init     ! initial global ice area
     &,    asum_final    ! final global ice area

      real (kind=dbl_kind), dimension(ntracer), intent(in), optional ::
     &     atsum_init    ! initial global ice area*tracer
     &,    atsum_final   ! final global ice area*tracer
!
!EOP
!
      integer (kind=int_kind) ::
     &     nt            ! tracer index

      real (kind=dbl_kind) ::
     &     diff          ! difference between initial and final values

      if (asum_init > puny) then
         diff = asum_init - asum_final
         if (abs(diff/asum_init) > puny) then
            write (nu_diag,*)
            write (nu_diag,*) 'Ice area conserv error, istep =',istep1
            write (nu_diag,*) 'Initial global area =', asum_init
            write (nu_diag,*) 'Final global area =', asum_final
            write (nu_diag,*) 'Fractional error =', abs(diff)/asum_init
            write (nu_diag,*) 'asum_final-asum_init =', diff 
            call abort_ice ('remap transport: conservation error')
         endif
      endif

      if (present(atsum_init)) then
       do nt = 1, ntracer
         if (atsum_init(nt) > puny) then
            diff = atsum_init(nt) - atsum_final(nt)
            if (abs(diff/atsum_init(nt)) > puny) then
               write (nu_diag,*)
               write (nu_diag,*) 'area*tracer conserv error, istep =',
     &                             istep1
               write (nu_diag,*) 'tracer index =', nt
               write (nu_diag,*) 'Initial global area*tracer =',
     &                            atsum_init(nt)
               write (nu_diag,*) 'Final global area*tracer =', 
     &                            atsum_final(nt)
               write (nu_diag,*) 'Fractional error =',
     &                            abs(diff)/atsum_init(nt)
               write (nu_diag,*) 'atsum_final-atsum_init =', diff 
               call abort_ice ('remap transport: conservation error')
            endif
         endif
       enddo
      endif                     ! present(atsum_init)

      end subroutine global_conservation

!=======================================================================
!
!BOP
!
! !IROUTINE: local_max_min - compute local max and min of a scalar field
!
! !INTERFACE:
!
      subroutine local_max_min (trm, tmin, tmax, aimask, trmask)
!
! !DESCRIPTION:
!
! At each grid point, compute the local max and min of a scalar 
! field phi: i.e., the max and min values in the nine-cell region 
! consisting of the home cell and its eight neighbors, plus the
! neighbors of the neighbors (25 cells in all). 
!
! !REVISION HISTORY:
!
! author William H. Lipscomb, LANL
!
! !USES:
!
! !INPUT/OUTPUT PARAMETERS:
!
      real (kind=dbl_kind), intent(in), 
     &     dimension(imt_local,jmt_local) ::
     &     aimask

      real (kind=dbl_kind), intent(in),
     &     dimension (imt_local,jmt_local,ntracer) ::
     &     trm
     &,    trmask

      real (kind=dbl_kind), intent(out),
     &     dimension (imt_local,jmt_local,ntracer) ::
     &     tmin           ! local min tracer
     &,    tmax           ! local max tracer
!
!EOP
!
      integer (kind=int_kind) ::
     &     i, j           ! horizontal indices
     &,    nt, nt1        ! tracer indices

      real (kind=dbl_kind), dimension(imt_local,jmt_local) ::
     &     phimask        ! aimask or trmask, as appropriate

      real (kind=dbl_kind) ::
     &     phi_nw, phi_n, phi_ne   ! field values in 8 neighbor cells 
     &,    phi_w, phi_e
     &,    phi_sw, phi_s, phi_se

      do nt = 1, ntracer

         if (depend(nt)==c0) then  ! does not depend on another tracer

            do j = 1, jmt_local
            do i = 1, imt_local
               phimask(i,j) = aimask(i,j)
            enddo
            enddo

         else   ! depends on another tracer

            nt1 = depend(nt)
            do j = 1, jmt_local
            do i = 1, imt_local
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

               phi_nw = phimask(i-1,j+1) * trm(i-1,j+1,nt) 
     &            + (c1-phimask(i-1,j+1))* trm(i,  j,  nt)
               phi_n  = phimask(i,  j+1) * trm(i,  j+1,nt) 
     &            + (c1-phimask(i,  j+1))* trm(i,  j,  nt)
               phi_ne = phimask(i+1,j+1) * trm(i+1,j+1,nt) 
     &            + (c1-phimask(i+1,j+1))* trm(i,  j,  nt)
               phi_w  = phimask(i-1,j)   * trm(i-1,j,  nt) 
     &            + (c1-phimask(i-1,j))  * trm(i,  j,  nt)
               phi_e  = phimask(i+1,j)   * trm(i+1,j,  nt) 
     &            + (c1-phimask(i+1,j))  * trm(i,  j,  nt)
               phi_sw = phimask(i-1,j-1) * trm(i-1,j-1,nt) 
     &            + (c1-phimask(i-1,j-1))* trm(i,  j,  nt)
               phi_s  = phimask(i,  j-1) * trm(i,  j-1,nt) 
     &            + (c1-phimask(i,  j-1))* trm(i,  j,  nt)
               phi_se = phimask(i+1,j-1) * trm(i+1,j-1,nt) 
     &            + (c1-phimask(i+1,j-1))* trm(i,  j,  nt)

!-----------------------------------------------------------------------
!     Compute the minimum and maximum among the nine local cells.
!-----------------------------------------------------------------------

               tmax(i,j,nt) = max (phi_nw, phi_n,  phi_ne, phi_w, 
     &                trm(i,j,nt), phi_e,  phi_sw, phi_s,  phi_se)   

               tmin(i,j,nt) = min (phi_nw, phi_n,  phi_ne, phi_w, 
     &                trm(i,j,nt), phi_e,  phi_sw, phi_s,  phi_se)

            enddo               ! i
         enddo                  ! j

! Now extend by one more cell, to include the neighbors of the neighbors.

         call bound (tmax(:,:,nt))
         call bound (tmin(:,:,nt))

         do j = jlo, jhi
         do i = ilo, ihi

            tmax(i,j,nt) = 
     &        max (tmax(i-1,j+1,nt), tmax(i,j+1,nt), tmax(i+1,j+1,nt),
     &             tmax(i-1,j,  nt), tmax(i,j,  nt), tmax(i+1,j,  nt),
     &             tmax(i-1,j-1,nt), tmax(i,j-1,nt), tmax(1+1,j-1,nt))

            tmin(i,j,nt) =
     &        min (tmin(i-1,j+1,nt), tmin(i,j+1,nt), tmin(i+1,j+1,nt),
     &             tmin(i-1,j,  nt), tmin(i,j,  nt), tmin(i+1,j,  nt),
     &             tmin(i-1,j-1,nt), tmin(i,j-1,nt), tmin(1+1,j-1,nt))

            enddo               ! i
         enddo                  ! j

      enddo                     ! nt

      end subroutine local_max_min

!======================================================================
!
!BOP
!
! !IROUTINE: check_monotonicity - check bounds on new tracer values
!
! !INTERFACE:
!
      subroutine check_monotonicity (aim, trm, tmin, tmax)
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
      use ice_exit
!
! !INPUT/OUTPUT PARAMETERS:
!
      real (kind=dbl_kind), intent(in),
     &     dimension (imt_local,jmt_local) ::
     &     aim            ! new ice area

      real (kind=dbl_kind), intent(in),
     &     dimension (imt_local,jmt_local,ntracer) ::
     &     trm            ! new tracers
     &,    tmin           ! local min tracer
     &,    tmax           ! local max tracer
!
!EOP
!
      integer (kind=int_kind) ::
     &     i, j           ! horizontal indices
     &,    nt, nt1        ! tracer indices

      real (kind=dbl_kind) ::
     &     w1, w2         ! work variables

      logical (kind=log_kind), dimension(imt_local,jmt_local) ::
     &     l_check        ! if true, check monotonicity

      do nt = 1, ntracer

    !-------------------------------------------------------------------
    ! Load logical array to identify tracers that need checking.   
    !-------------------------------------------------------------------

         if (depend(nt)==0) then ! does not depend on another tracer

            do j = jlo, jhi
            do i = ilo, ihi
               if (aim(i,j) > puny) then
                  l_check(i,j) = .true.
               else
                  l_check(i,j) = .false.
               endif
            enddo
            enddo

         else                   ! depends on another tracer

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
                  write (nu_diag,*) ''
                  write (nu_diag,*) 'new tracer < tmin, istep =', istep1
                  write (nu_diag,*) 'my_task, i, j, nt =',my_task,i,j,nt
                  write (nu_diag,*) 'new tracer =', trm (i,j,nt)
                  write (nu_diag,*) 'tmin ='      , tmin(i,j,nt)
                  write (nu_diag,*) 'ice area =', aim(i,j)
                  call abort_ice('remap transport: monotonicity error')
               elseif (trm(i,j,nt) > tmax(i,j,nt)+w2*puny) then
                  write (nu_diag,*) ''
                  write (nu_diag,*) 'new tracer > tmax, istep =', istep1
                  write (nu_diag,*) 'my_task, i, j, nt =',my_task,i,j,nt
                  write (nu_diag,*) 'new tracer =', trm (i,j,nt)
                  write (nu_diag,*) 'tmax ='      , tmax(i,j,nt)
                  write (nu_diag,*) 'ice area =', aim(i,j)
                  call abort_ice('remap transport: monotonicity error')
               endif
            endif

         enddo                  ! i
         enddo                  ! j
         
      enddo                     ! nt

      end subroutine check_monotonicity

!=======================================================================

      end module ice_transport_remap

!=======================================================================






