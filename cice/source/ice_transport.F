c $Id: $
c=======================================================================
!---! horizontal advection 
!---!
!---! authors Elizabeth C. Hunke
!---!         William H. Lipscomb
!---!         Fluid Dynamics Group, Los Alamos National Laboratory
!---!
!---! Copyright, 2001.  The Regents of the University of California.
!---! This software was produced under a U.S. Government contract 
!---! (W-7405-ENG-36) by Los Alamos National Laboratory, which is 
!---! operated by the University of California for the U.S. Department 
!---! of Energy.  The U.S. Government is licensed to use, reproduce, and 
!---! distribute this software.  Permission is granted to the public to 
!---! copy and use this software without charge, provided that this 
!---! Notice and any statement of authorship are reproduced on all 
!---! copies.  Neither the Government nor the University makes any 
!---! warranty, express or implied, or assumes any liability or 
!---! responsibility for the use of this software.
c=======================================================================

      module ice_transport

      use ice_model_size
      use ice_domain
      use ice_constants
      use ice_grid
      use ice_fileunits

      implicit none

      character (len=char_len) :: 
     &   advection   ! type of advection algorithm used
                     ! upwind  => 1st order mpdata scheme (donor cell)
                     ! mpdata2 => 2nd order mpdata scheme
                     ! remap   => remapping scheme

c=======================================================================

      contains

c=======================================================================

      subroutine transport_mpdata

!---!-------------------------------------------------------------------
!---! computes the transport equations for one timestep
!---!-------------------------------------------------------------------

      use ice_flux
      use ice_timers
      use ice_state
      use ice_calendar
      use ice_mechred_cice

      integer (kind=int_kind), parameter :: narr = 1 + 5*ncat ! number
                                            ! of state variable arrays
      integer (kind=int_kind) :: 
     &   i, j, layer, nc  ! standard indices
     &,  narrays          ! counter for number of state variable arrays

      real (kind=dbl_kind) :: works(imt_local,jmt_local,narr)
      real (kind=dbl_kind) :: worke(imt_local,jmt_local,ntilay)
      real (kind=dbl_kind) :: asum  ! frac area of ice plus open water 

      call ice_timer_start(3)  ! advection

      !-----------------------------------------------------------------
      ! fill work arrays with fields to be advected
      !-----------------------------------------------------------------
      ! two arrays are used for performance (balance memory/cache vs 
      ! number of bound calls);  one array or more than two may perform 
      ! better depending on the machine used, number of processors, etc.
      ! --tested on SGI R2000, using 4 pes for the ice model under MPI
      !-----------------------------------------------------------------
      do j=jlo,jhi
       do i=ilo,ihi
        works(i,j,1) = aice0(i,j)
        narrays = 1
        do nc=1,ncat
         works(i,j,narrays+1) = aicen(i,j,nc)
         works(i,j,narrays+2) = vicen(i,j,nc)
         works(i,j,narrays+3) = vsnon(i,j,nc)
         works(i,j,narrays+4) = -aicen(i,j,nc)*Tsfcn(i,j,nc)
         works(i,j,narrays+5) = -(esnon(i,j,nc)
     &                          +rhos*Lfresh*vsnon(i,j,nc)) ! for Tsnow<0
         narrays = narrays + c5    ! if there were 3 arrays in this 
        enddo                      ! loop, use c3 instead and reset narr
        do layer=1,ntilay
         worke(i,j,layer) = -eicen(i,j,layer)
        enddo
       enddo
      enddo
      if (narr .ne. narrays) 
     & write(nu_diag,*) "Wrong number of arrays in transport bound call"

      !-----------------------------------------------------------------
      ! update local domain boundaries
      !-----------------------------------------------------------------

      call bound_narr(narr,  works) 
      call bound_narr(ntilay,worke)

      !-----------------------------------------------------------------
      ! advect
      !-----------------------------------------------------------------

      call mpdata(narr,  works)
      call mpdata(ntilay,worke)

      !-----------------------------------------------------------------
      ! update local domain boundaries again
      !-----------------------------------------------------------------

      call bound_narr(narr,  works) 
      call bound_narr(ntilay,worke)

      !-----------------------------------------------------------------
      ! retrieve advected fields from work array
      !-----------------------------------------------------------------
       do j=1,jmt_local
        do i=1,imt_local
        aice0(i,j) = works(i,j,1)
        narrays = 1
        do nc=1,ncat
         aicen(i,j,nc) = works(i,j,narrays+1)
         vicen(i,j,nc) = works(i,j,narrays+2)
         vsnon(i,j,nc) = works(i,j,narrays+3)
         esnon(i,j,nc) = -works(i,j,narrays+5)
     &                          -rhos*Lfresh*vsnon(i,j,nc) ! for Tsnow<0 
        narrays = narrays + c5     ! if there were 3 arrays in this 
        enddo                      ! loop, use c3 instead and reset narr
        do layer=1,ntilay
         eicen(i,j,layer) = -worke(i,j,layer)
        enddo
       enddo
      enddo

      do j=jlo,jhi
       do i=ilo,ihi
        narrays = 1
        do nc=1,ncat
         if (aicen(i,j,nc).gt.puny) then     ! advected surface temperature
          Tsfcn(i,j,nc) = -works(i,j,narrays+4)/aicen(i,j,nc)
         else
          Tsfcn(i,j,nc) = Tf(i,j)
         endif
         narrays = narrays + c5    ! if there were 3 arrays in this 
        enddo                      ! loop, use c3 instead and reset narr
       enddo
      enddo

      !-----------------------------------------------------------------
      ! compute divergence rate due to advection
      !-----------------------------------------------------------------

      do j=jlo,jhi
       do i=ilo,ihi
        asum = aice0(i,j)
        do nc = 1, ncat
         asum = asum + aicen(i,j,nc)
        enddo
        divu_adv(i,j) = (c1-asum)/dt
       enddo
      enddo

      call ice_timer_stop(3) ! advection

      !-----------------------------------------------------------------
      ! mask
      !-----------------------------------------------------------------
      do nc=1,ncat
       do j=1,jmt_local
        do i=1,imt_local
         if (.not.tmask(i,j)) then
          aicen(i,j,nc)    = c0
          vicen(i,j,nc)    = c0
          vsnon(i,j,nc)    = c0
          Tsfcn(i,j,nc)    = c0
          esnon(i,j,nc)    = c0
         endif
        enddo
       enddo
      enddo

      do layer=1,ntilay
       do j=1,jmt_local
        do i=1,imt_local
         if (.not.tmask(i,j)) then
          eicen(i,j,layer) = c0
         endif
        enddo
       enddo
      enddo

      end subroutine transport_mpdata

c=======================================================================

      subroutine mpdata(narrays,phi)

!---!-------------------------------------------------------------------
!---! Smolarkiewicz, P. K., 1984:  A fully multidimensional positive 
!---! definite advection transport algorithm with small implicit 
!---! diffusion, J. Comput. Phys., 54, 325-362.
!---!-------------------------------------------------------------------

      use ice_calendar
      use ice_dyn_evp

      integer (kind=int_kind), intent(in) :: narrays

      real (kind=dbl_kind), intent(inout) ::
     &   phi(imt_local,jmt_local,narrays)

      integer (kind=int_kind) :: i, j, k, ix, iy, n
      integer (kind=int_kind), parameter :: 
     &   iord=3  ! 2nd order MPDATA paramter

      real (kind=dbl_kind) ::
     &   dive(ilo:ihi,jlo:jhi)
     &,  divn(ilo:ihi,jlo:jhi)
     &,  phiavg(imt_local,jmt_local)
     &,  uee(imt_local,jmt_local,narrays)
     &,  vnn(imt_local,jmt_local,narrays)
     &,  fx(imt_local,jmt_local) 
     &,  fy(imt_local,jmt_local)
     &,  eps
     &,  upwind, y1, y2, a, h   ! function

      upwind(y1,y2,a,h)=p5*dt*h*((a+abs(a))*y1+(a-abs(a))*y2)

      eps = eps15

      do n = 1,narrays
       do j=jlo,jhi
        do i=ilo,ihi
         uee(i,j,n)=p5*(u(i,j)+u(i,j-1))
         vnn(i,j,n)=p5*(v(i,j)+v(i-1,j))
        enddo
       enddo
      enddo ! narrays

      call bound_narr(narrays,uee)
      call bound_narr(narrays,vnn)
 
      ! upwind 
      do n = 1,narrays
       do j=1,jhi
        do i=1,ihi
         fx(i,j)=upwind(phi(i,j,n),phi(i+1,j,n),uee(i,j,n),HTE(i,j))
         fy(i,j)=upwind(phi(i,j,n),phi(i,j+1,n),vnn(i,j,n),HTN(i,j))
        enddo
       enddo
       do j=jlo,jhi
        do i=ilo,ihi
         phi(i,j,n)=phi(i,j,n)-(fx(i,j)-fx(i-1,j)+fy(i,j)-fy(i,j-1))
     &                          /tarea(i,j)
        enddo
       enddo
      enddo ! narrays

      call bound_narr(narrays,phi)

      if (advection .eq. 'upwind') goto 555   ! Upwind (Donor Cell)

      do k=1,iord                             ! 2nd order MPDATA

      do n = 1,narrays
       do j=1,jhi
        do i=1,ihi
         phiavg(i,j) = p25*(phi(i,j,n) + phi(i+1,j,n)
     &               + phi(i+1,j+1,n) + phi(i,j+1,n))
        enddo
       enddo

       do j=jlo,jhi
        do i=ilo,ihi
         dive(i,j) = 
     1  ((dyt(i+1,j)*(uee(i+1,j,n)+uee(i,j,n))*phi(i+1,j,n)
     1  - dyt(i,j)*(uee(i-1,j,n)+uee(i,j,n))*phi(i,j,n))
     1  /(phi(i+1,j,n)+phi(i,j,n)+eps)
     2  + (dxu(i,j)*(vnn(i+1,j,n)+vnn(i,j,n))*phiavg(i,j)
     2  - dxu(i,j-1)*(vnn(i+1,j-1,n)+vnn(i,j-1,n))*phiavg(i,j-1))
     2  /(phiavg(i,j)+phiavg(i,j-1)+eps))
     3         /(HTE(i,j)*(dxu(i,j)+dxu(i,j-1)))

         divn(i,j) = 
     1  ((dxt(i,j+1)*(vnn(i,j+1,n)+vnn(i,j,n))*phi(i,j+1,n)
     1  - dxt(i,j)*(vnn(i,j-1,n)+vnn(i,j,n))*phi(i,j,n))
     1  /(phi(i,j+1,n)+phi(i,j,n)+eps)
     2  + (dyu(i,j)*(uee(i,j+1,n)+uee(i,j,n))*phiavg(i,j)
     2  - dyu(i-1,j)*(uee(i-1,j,n)+uee(i-1,j+1,n))*phiavg(i-1,j))
     2  /(phiavg(i,j)+phiavg(i-1,j)+eps))
     3         /(HTN(i,j)*(dyu(i,j)+dyu(i-1,j)))
        enddo
       enddo

       ! antidiffusive velocities
       do j=jlo,jhi
        do i=ilo,ihi
         uee(i,j,n) = abs(uee(i,j,n))*(phi(i+1,j,n)-phi(i,j,n))
     1   /(phi(i+1,j,n)+phi(i,j,n)+eps) - dt*uee(i,j,n)*dive(i,j)

         vnn(i,j,n) = abs(vnn(i,j,n))*(phi(i,j+1,n)-phi(i,j,n))
     1   /(phi(i,j+1,n)+phi(i,j,n)+eps) - dt*vnn(i,j,n)*divn(i,j)
        enddo
       enddo
      enddo ! narrays

      call bound_narr(narrays,uee)
      call bound_narr(narrays,vnn)

      ! upwind with antidiffusive velocities
      do n = 1,narrays
       do j=1,jhi
        do i=1,ihi
         fx(i,j)=upwind(phi(i,j,n),phi(i+1,j,n),uee(i,j,n),HTE(i,j))
         fy(i,j)=upwind(phi(i,j,n),phi(i,j+1,n),vnn(i,j,n),HTN(i,j))
        enddo
       enddo

       ix=-1
       do j=jlo,jhi
        do i=ilo,ihi
         phi(i,j,n)=phi(i,j,n)-(fx(i,j)-fx(i-1,j)+fy(i,j)-fy(i,j-1))
     &                          /tarea(i,j)

         if (phi(i,j,n).lt.-eps12) then 
          ix=i
          iy=j
         elseif (phi(i,j,n).lt.0.) then
          phi(i,j,n) = c0 
         endif
        enddo
       enddo
       if (ix.ge.0) then
        write (nu_diag,*)  my_task,ix,iy,' transport unstable ',istep,k
        write (nu_diag,*)  'mpdata phi = ',phi(ix,iy,n),' n = ',n
        stop
       endif
      enddo ! narrays

      call bound_narr(narrays,phi)

      enddo ! iord
 555  continue

      end subroutine mpdata

c=======================================================================

      subroutine transport_remap

!---!-------------------------------------------------------------------
!---! Remapping transport scheme
!---! 
!---! authors William H. Lipscomb and John R. Baumgardner, LANL
!---!
!---! This subroutine solves the transport equations for one timestep 
!---! using the second-order conservative remapping scheme developed by 
!---! John Dukowicz and John Baumgardner (DB) and modified for sea ice 
!---! by William Lipscomb and Elizabeth Hunke.  
!---!
!---! This scheme is second-order accurate and conservative.  It is
!---! compatible for tracers; i.e, it does not produce new extrema 
!---! in thickness or enthalpy.  It is monotone for ice area if the 
!---! velocity is non-divergent.
!---!
!---! References:
!---! 
!---! Dukowicz, J. K., and J. R. Baumgardner, 2000: Incremental 
!---!  remapping as a transport/advection algorithm, J. Comput. Phys., 
!---!  160, 318-335.
!---!
!---! Lipscomb, W. H., and E. C. Hunke, 2001: Modeling sea ice
!---!  transport with incremental remapping, in preparation.
!---!-------------------------------------------------------------------
    !
    ! An outline of the scheme:
    !
    ! (1) Given velocities at the cell corners, identify the departure 
    !     regions for fluxes across each cell face.
    ! 
    !     (a) Compute areas and vertices of departure triangles for 
    !         east cell faces.  From the triangle vertices, compute 
    !         the coordinates of points used in triangle integration 
    !         formulas.
    !     (b) Repeat for north cell faces.
    !
    ! (2) Compute fluxes and update state variables (aicen, Tsfcn,
    !       vicen, vsnon, eicen, esnon).
    !   
    !     (a) Construct 2D scalar fields for each variable in each cell.
    !         This involves finding the field value in the geometric 
    !         center of the cell, along with gradients in the x and y 
    !         directions.  The gradient is limited to preserve 
    !         monotonicity.
    !     (b) Given these field reconstructions, integrate the fluxes
    !         over each departure triangle for each cell face.
    !     (c) Transfer the fluxes across the cell faces to obtain new
    !         values of the state variables.
    !
    ! Note that part (1) is done only once, whereas part (2) is 
    ! repeated for each ice category. 
    !-------------------------------------------------------------------

      use ice_timers
      use ice_state
      use ice_calendar
      use ice_dyn_evp
      use ice_mpi_internal
      use ice_itd
      use ice_mechred_cice

      integer (kind=int_kind) ::
     &   i, j             ! horizontal indices
     &,  k                ! layer index
     &,  n                ! ice category index
     &,  m                ! triangle index


      integer (kind=int_kind), dimension(8), parameter :: 
           ! i and j values of a grid cell compared to the home cell
     &   ishift_e = (/ 1, 0, 1, 0, 0, 0, 1, 1 /)
     &,  jshift_e = (/ 1, 1,-1,-1, 0, 0, 0, 0 /)
                    ! NE  N SE  S  H  H  E  E

     &,  ishift_n = (/-1,-1, 1, 1, 0, 0, 0, 0 /)
     &,  jshift_n = (/ 1, 0, 1, 0, 0, 0, 1, 1 /)
                    ! NW  W NE  E  H  H  N  N
           

      real (kind=dbl_kind) ::
     &   x1, y1, x2, y2   ! x,y coordinates of backward trajectories from
                          ! corners of east and north faces, as in DB
     &,  xa, ya, xb, yb   ! x,y coordinates of points where the lines joining
                          ! (x1,y1) and (x2,y2) cross cell faces, as in DB
     &,  x1t,y1t,x2t,y2t  ! x1, y1, x2, y2 transformed to T-cell ref frame
     &,  xat,yat,xbt,ybt  ! xa, ya, xb, yb transformed to T-cell ref frame
     &,  aicesum          ! sum of aicen over categories
     &,  aixsum           ! sum of aixn over categories
     &,  aiysum           ! sum of aiyn over categories
     &,  aTsfold          ! old value of aicen*Tsfcn
     &,  w1, w2, w3, w4, w5, w6, w7  ! work variables


      real (kind=dbl_kind), dimension (imt_local, jmt_local) ::
     &   iceflag          ! = 1. if ice is present, = 0. otherwise
     &,  snoflag          ! = 1. if snow is present, = 0. otherwise
     &,  hice             ! ice thickness (m)
     &,  hsno             ! snow thickness (m)
     &,  qsno             ! snow enthalpy/unit volume (J/m^3)
     &,  cax, cay         ! x,y coordinates of center of ice area
     &,  cvix, cviy       ! x,y coordinates of center of ice volume
     &,  cvsx, cvsy       ! x,y coordinates of center of snow volume
     &,  aigc             ! ice area at geometric center of cell
     &,  aiflxe, aiflxn   ! flux of aice thru east and north cell faces
     &,  Tsfgc            ! surface temperature at geometric center of cell
     &,  Tsfx,Tsfy        ! limited derivative of sfc temp wrt x and y
     &,  aTflxe,aTflxn    ! flux of aice*Tsfc thru E and N cell faces
     &,  higc             ! ice thickness at geometric center of cell
     &,  hix, hiy         ! limited derivative of ice thickness wrt x and y
     &,  viflxe, viflxn   ! flux of ice volume thru E and N cell faces
     &,  hsgc             ! snow thickness at geometric center of cell
     &,  hsx, hsy         ! limited derivative of snow thickness wrt x and y
     &,  vsflxe, vsflxn   ! flux of snow volume thru E and N cell faces
     &,  qsgc             ! snow enthalpy at geometric center of cell
     &,  qsx, qsy         ! limited derivative of snow enthalpy wrt x and y
     &,  esflxe, esflxn   ! flux of snow energy (esno) thru E and N cell faces
     &,  areaflxe, areaflxn  ! flux of total area thru E and N cell faces


      real (kind=dbl_kind), dimension (imt_local, jmt_local, ncat) ::
     &   aixn, aiyn       ! limited derivative of ice area wrt x and y


      real (kind=dbl_kind), dimension (imt_local, jmt_local, nilyr) ::
     &   qice             ! nilyr * ice enthalpy/unit volume (J/m^3)
                          ! (saves time not to divide by nilyr)
     &,  qigc             ! ice enthalpy at geometric center of cell
     &,  qix, qiy         ! limited derivative of ice enthalpy wrt x and y
     &,  eiflxe, eiflxn   ! flux of ice energy (eice) thru E and N cell faces


      real (kind=dbl_kind), dimension (8, imt_local, jmt_local) ::
     &   area_e           ! area of east-face flux triangle
     &,  area_n           ! area of north-face flux triangle
                          ! (there are 8 kinds of triangles, as in DB)

      real (kind=dbl_kind), dimension (0:3, 8, imt_local, jmt_local) ::
     &	 xpe, ype         ! x and y coordinates of special triangle points
                          ! (need 4 points for triangle integrals)
     &,	 xpn, ypn         ! e for east faces, n for north faces
                          

      real (kind=dbl_kind) ::  ! globally conserved quantities
     &   aicegsum1, aicegsum2  ! aice at start and end of time step
     &,  vicegsum1, vicegsum2  ! similarly for vice, etc
     &,  vsnogsum1, vsnogsum2
     &,  eicegsum1, eicegsum2
     &,  esnogsum1, esnogsum2

      logical, dimension (8, imt_local, jmt_local) :: 
     &    cnd_e, cnd_n    ! true if there is a flux contribution from
                          ! a given triangle (1-8) in a given cell
      logical :: cnd2     ! same as cnd_e or cnd_n

      call ice_timer_start(3)   ! advection

    !-------------------------------------------------------------------
    ! Compute global sums of conserved properties.  
    ! (This is a bug check; it is not part of the numerical scheme.)
    !-------------------------------------------------------------------

c      call conserved_sums (aicegsum1, vicegsum1, vsnogsum1,
c     &	     eicegsum1, esnogsum1)

    !-------------------------------------------------------------------
    ! Part 1: Compute areas and vertices of flux triangles for east 
    !         and north cell faces.
    !
    ! Triangle notation:
    ! For each face, there are 20 triangles that can contribute a flux,
    ! but many of these are mutually exclusive.  It turns out that
    ! at most 8 triangles must be evaluated.
    !
    ! For the east face, these triangles are referred to as:
    ! (1) NE, NE1, NE2
    ! (2) N, N1, N2
    ! (3) SE, SE1, SE2
    ! (4) S, S1, S2
    ! (5) H1a, H1b
    ! (6) H2a, H2b
    ! (7) E1a, E1b
    ! (8) E2a, E2b
    !
    ! For the north face, these triangles are referred to as:
    ! (1) NW, NW1, NW2
    ! (2) W, W1, W2
    ! (3) NE, NE1, NE2
    ! (4) E, E1, E2
    ! (5) H1a, H1b
    ! (6) H2a, H2b
    ! (7) N1a, N1b
    ! (8) N2a, N2b
    !
    ! See Figure 3 in DB for pictures of these triangles.
    ! See Table 1 in DB for logical conditions.
    !
    ! The vertices of many triangles lie at points whose coordinates are
    ! (x1,y1), (x2,y2), (xa,0), (xb,0), (0,ya), and (0,yb) in a
    ! reference frame whose origin is the cell corner.  These
    ! coordinates must be transformed to the reference frame whose
    ! origin is the geometric center of the T-cell in which the triangle
    ! is located.  The transformation is carried out using pre-computed
    ! 2x2 matrices.  There are 4 matrices (one for each corner)
    ! associated with each grid cell.  They do not describe pure
    ! rotations, because they do not preserve length.
    !
    ! The function 'triarea' computes the triangle area, which is always
    ! positive, but the variable 'area' is positive or negative
    ! depending on the flow direction across the face.
    !-------------------------------------------------------------------

    !-------------------------------------------------------------------
    ! East faces
    !-------------------------------------------------------------------

      ! initialize cnd
      do j = 1, jmt_local
         do i = 1, imt_local
            do m = 1, 8
               cnd_e(m,i,j) = .false.
            enddo
         enddo
      enddo

      do j = jlo, jhi
         do i = ilo-1, ihi  ! includes W face of cells with index ilo

            ! velocity coordinates
            x1 = -dt*u(i,j)
            y1 = -dt*v(i,j)
            x2 = -dt*u(i,j-1)
            y2 = -dt*v(i,j-1)
            w1 =  c1 / (y2 - HTE(i,j)  - y1)
            xa = (x1*(y2 - HTE(i,j)) - x2*y1) * w1
            xb = (x1*y2 - x2*(y1 + HTE(i,j))) * w1
            if (abs(xb-xa) .ge. puny) then
               ya = xa * HTE(i,j) / (xb - xa)
               yb = ya + HTE(i,j)
            else
               ya = c0
               yb = c0
            endif

            ! flux contribution from triangles in NE cell
            cnd_e(1,i+1,j+1) =  xa.gt.c0 .and. y1.gt.c0         ! NE
            x1t = msw(1,1,i+1,j+1)*x1 + msw(1,2,i+1,j+1)*y1
            y1t = msw(2,1,i+1,j+1)*x1 + msw(2,2,i+1,j+1)*y1
            xat = msw(1,1,i+1,j+1)*xa
            yat = msw(2,1,i+1,j+1)*xa
            xpe(1,1,i+1,j+1) = -p5*HTS(i+1,j+1)
            xpe(2,1,i+1,j+1) =  xpe(1,1,i+1,j+1) + xat
            xpe(3,1,i+1,j+1) =  xpe(1,1,i+1,j+1) + x1t
            ype(1,1,i+1,j+1) = -p5*HTW(i+1,j+1)
            ype(2,1,i+1,j+1) =  ype(1,1,i+1,j+1) + yat
            ype(3,1,i+1,j+1) =  ype(1,1,i+1,j+1) + y1t

            if (cnd_e(1,i+1,j+1) .and. x1.lt.c0) then           ! NE2
                          ! If x1 .ge. 0, we have NE triangle.
                          ! Similarly below for E, SE, and S.
               xat = msw(1,2,i+1,j+1)*ya
               yat = msw(2,2,i+1,j+1)*ya
               xpe(3,1,i+1,j+1) = xpe(1,1,i+1,j+1) + xat
               ype(3,1,i+1,j+1) = ype(1,1,i+1,j+1) + yat
            endif
            
            area_e(1,i+1,j+1) = -triarea(xpe(1:3,1,i+1,j+1),
     &                                   ype(1:3,1,i+1,j+1))

            cnd2 = xa.lt.c0 .and. x1.gt.c0 .and. y1.gt.c0       ! NE1 
            if (cnd2) then
               xat = msw(1,2,i+1,j+1)*ya
               yat = msw(2,2,i+1,j+1)*ya
               xpe(2,1,i+1,j+1) = xpe(1,1,i+1,j+1) + xat
               ype(2,1,i+1,j+1) = ype(1,1,i+1,j+1) + yat
               area_e(1,i+1,j+1) = triarea(xpe(1:3,1,i+1,j+1),
     &                                     ype(1:3,1,i+1,j+1))
            endif
            
            cnd_e(1,i+1,j+1) = cnd_e(1,i+1,j+1) .or. cnd2

            ! flux contribution from triangles in N cell
            cnd_e(2,i,j+1) =  xa.lt.c0 .and. y1.gt.c0           ! N
            x1t = mse(1,1,i,j+1)*x1 + mse(1,2,i,j+1)*y1
            y1t = mse(2,1,i,j+1)*x1 + mse(2,2,i,j+1)*y1
            xat = mse(1,1,i,j+1)*xa
            yat = mse(2,1,i,j+1)*xa
            xpe(1,2,i,j+1) =  p5*HTS(i,j+1) 
            xpe(2,2,i,j+1) =  xpe(1,2,i,j+1) + xat 
            xpe(3,2,i,j+1) =  xpe(1,2,i,j+1) + x1t
            ype(1,2,i,j+1) = -p5*HTE(i,j+1) 
            ype(2,2,i,j+1) =  ype(1,2,i,j+1) + yat
            ype(3,2,i,j+1) =  ype(1,2,i,j+1) + y1t 

            if (cnd_e(2,i,j+1) .and. x1.gt.c0) then             ! N1
               xat = mse(1,2,i,j+1)*ya
               yat = mse(2,2,i,j+1)*ya
               xpe(3,2,i,j+1) = xpe(1,2,i,j+1) + xat
               ype(3,2,i,j+1) = ype(1,2,i,j+1) + yat 
            endif

            area_e(2,i,j+1) = triarea(xpe(1:3,2,i,j+1),
     &                                ype(1:3,2,i,j+1))

            cnd2 = xa.gt.c0 .and. x1.lt.c0 .and. y1.gt.c0       ! N2
            if (cnd2) then
               xat = mse(1,2,i,j+1)*ya
               yat = mse(2,2,i,j+1)*ya
               xpe(2,2,i,j+1) = xpe(1,2,i,j+1) + xat
               ype(2,2,i,j+1) = ype(1,2,i,j+1) + yat
               area_e(2,i,j+1) = -triarea(xpe(1:3,2,i,j+1),
     &                                    ype(1:3,2,i,j+1))
            endif

            cnd_e(2,i,j+1) = cnd_e(2,i,j+1) .or. cnd2

            ! flux contribution from triangles in SE cell
            cnd_e(3,i+1,j-1) = xb.gt.c0 .and. y2.lt.c0          ! SE

            x2t = mnw(1,1,i+1,j-1)*x2 + mnw(1,2,i+1,j-1)*y2
            y2t = mnw(2,1,i+1,j-1)*x2 + mnw(2,2,i+1,j-1)*y2
            xbt = mnw(1,1,i+1,j-1)*xb
            ybt = mnw(2,1,i+1,j-1)*xb
            xpe(1,3,i+1,j-1) = -p5*HTN(i+1,j-1)
            xpe(2,3,i+1,j-1) =  xpe(1,3,i+1,j-1) + xbt
            xpe(3,3,i+1,j-1) =  xpe(1,3,i+1,j-1) + x2t
            ype(1,3,i+1,j-1) =  p5*HTW(i+1,j-1)
            ype(2,3,i+1,j-1) =  ype(1,3,i+1,j-1) + ybt
            ype(3,3,i+1,j-1) =  ype(1,3,i+1,j-1) + y2t

            if (cnd_e(3,i+1,j-1) .and. x2.lt.c0) then           ! SE2
               xbt = mnw(1,2,i+1,j-1)*yb
               ybt = mnw(2,2,i+1,j-1)*yb
               xpe(3,3,i+1,j-1) = xpe(1,3,i+1,j-1) + xbt
               ype(3,3,i+1,j-1) = ype(1,3,i+1,j-1) + ybt
            endif

            area_e(3,i+1,j-1) = -triarea(xpe(1:3,3,i+1,j-1),
     &                                   ype(1:3,3,i+1,j-1))

            cnd2 = xb.lt.c0 .and. x2.gt.c0 .and. y2.lt.c0       ! SE1
            if (cnd2) then
               xbt = mnw(1,2,i+1,j-1)*yb
               ybt = mnw(2,2,i+1,j-1)*yb
               xpe(2,3,i+1,j-1) = xpe(1,3,i+1,j-1) + xbt
               ype(2,3,i+1,j-1) = ype(1,3,i+1,j-1) + ybt
               area_e(3,i+1,j-1) = triarea(xpe(1:3,3,i+1,j-1),
     &                                     ype(1:3,3,i+1,j-1))
            endif

            cnd_e(3,i+1,j-1) = cnd_e(3,i+1,j-1) .or. cnd2

            ! flux contribution from triangles in S cell
            cnd_e(4,i,j-1) = xb.lt.c0 .and. y2.lt.c0            ! S

            x2t = mne(1,1,i,j-1)*x2 + mne(1,2,i,j-1)*y2
            y2t = mne(2,1,i,j-1)*x2 + mne(2,2,i,j-1)*y2
            xbt = mne(1,1,i,j-1)*xb
            ybt = mne(2,1,i,j-1)*xb
            xpe(1,4,i,j-1) =  p5*HTN(i,j-1)
            xpe(2,4,i,j-1) =  xpe(1,4,i,j-1) + xbt
            xpe(3,4,i,j-1) =  xpe(1,4,i,j-1) + x2t
            ype(1,4,i,j-1) =  p5*HTE(i,j-1)
            ype(2,4,i,j-1) =  ype(1,4,i,j-1) + ybt
            ype(3,4,i,j-1) =  ype(1,4,i,j-1) + y2t

            if (cnd_e(4,i,j-1) .and. x2.gt.c0) then             ! S1
               xbt = mne(1,2,i,j-1)*yb
               ybt = mne(2,2,i,j-1)*yb
               xpe(3,4,i,j-1) = xpe(1,4,i,j-1) + xbt
               ype(3,4,i,j-1) = ype(1,4,i,j-1) + ybt
            endif

            area_e(4,i,j-1) = triarea(xpe(1:3,4,i,j-1),
     &                                ype(1:3,4,i,j-1))

            cnd2 = xb.gt.c0 .and. x2.lt.c0 .and. y2.lt.c0       ! S2
            if (cnd2) then
               xbt = mne(1,2,i,j-1)*yb
               ybt = mne(2,2,i,j-1)*yb
               xpe(2,4,i,j-1) = xpe(1,4,i,j-1) + xbt
               ype(2,4,i,j-1) = ype(1,4,i,j-1) + ybt
               area_e(4,i,j-1) = -triarea(xpe(1:3,4,i,j-1),
     &                                    ype(1:3,4,i,j-1))
            endif

            cnd_e(4,i,j-1) = cnd_e(4,i,j-1) .or. cnd2

            ! flux from triangles inside home cell
            if (y1 .gt. c0) then
               x1 = xa
               y1 = c0
            endif
            
            if (y2 .lt. c0) then
               x2 = xb
               y2 = c0
            endif
 
            if (abs(xb-xa) .ge. puny) then
               w1 = min (c1, max(c0, xb/(xb-xa)))
            else
               w1 = c0
            endif
 
            x1t = mne(1,1,i,j)*x1 + mne(1,2,i,j)*y1
            y1t = mne(2,1,i,j)*x1 + mne(2,2,i,j)*y1
            x2t = mse(1,1,i,j)*x2 + mse(1,2,i,j)*y2
            y2t = mse(2,1,i,j)*x2 + mse(2,2,i,j)*y2

            ! flux contribution from triangle that includes the
            ! E cell face (part of convex quadrilateral inside home cell)
            xpe(1,5,i,j) =  p5*HTS(i,j)                         ! H1a
            xpe(2,5,i,j) =  p5*HTN(i,j)
            xpe(3,5,i,j) =  p5*HTN(i,j) + x1t
            ype(1,5,i,j) = -p5*HTE(i,j)
            ype(2,5,i,j) =  p5*HTE(i,j)
            ype(3,5,i,j) =  p5*HTE(i,j) + y1t
 
            ! flux contribution from triangle lying along the upper part 
            ! of E face for case of line xa-xb intersecting the face
            cnd_e(5,i,j) = xa*xb.lt.c0 .and. x1.lt.c0           ! H1b
            if (cnd_e(5,i,j)) then
               xpe(1,5,i,j) =  p5 * 
     &              (w1*(HTN(i,j)-HTS(i,j)) + HTS(i,j))
               ype(1,5,i,j) = (w1 - p5) * HTE(i,j)
            endif
            
            cnd_e(5,i,j) = cnd_e(5,i,j) .or. 
     &           (xa*xb.ge.c0 .and. xa+xb.lt.c0)

            area_e(5,i,j) = triarea(xpe(1:3,5,i,j),
     &                              ype(1:3,5,i,j))

            ! flux contribution from triangle touching but not
            ! lying along the E face (other part of convex quadrilateral)
            xpe(1,6,i,j) =  p5*HTS(i,j)                         ! H2a
            xpe(2,6,i,j) =  p5*HTN(i,j) + x1t
            xpe(3,6,i,j) =  p5*HTS(i,j) + x2t
            ype(1,6,i,j) = -p5*HTE(i,j)
            ype(2,6,i,j) =  p5*HTE(i,j) + y1t
            ype(3,6,i,j) = -p5*HTE(i,j) + y2t

            ! flux contribution from triangle lying along the lower part 
            ! of E face for case of line xa-xb intersecting the face
            cnd_e(6,i,j) = xa*xb.lt.c0 .and. x2.lt.c0           ! H2b
            if (cnd_e(6,i,j)) then
               xpe(2,6,i,j) =  p5 *
     &              (w1*(HTN(i,j)-HTS(i,j)) + HTS(i,j))
               ype(2,6,i,j) = (w1 - p5)*HTE(i,j)
            endif
            
            cnd_e(6,i,j) = cnd_e(6,i,j) .or.
     &           (xa*xb.ge.c0 .and. xa+xb.lt.c0)

            area_e(6,i,j) = triarea(xpe(1:3,6,i,j),
     &                              ype(1:3,6,i,j))

            ! flux from triangles in E cell
            x1t = mnw(1,1,i+1,j)*x1 + mnw(1,2,i+1,j)*y1
            y1t = mnw(2,1,i+1,j)*x1 + mnw(2,2,i+1,j)*y1
            x2t = msw(1,1,i+1,j)*x2 + msw(1,2,i+1,j)*y2
            y2t = msw(2,1,i+1,j)*x2 + msw(2,2,i+1,j)*y2

            ! flux contribution from triangle that includes the
            ! W cell face (part of convex quadrilateral inside E cell)
            xpe(1,7,i+1,j) = -p5*HTS(i+1,j)                     ! E1a
            xpe(2,7,i+1,j) = -p5*HTN(i+1,j)
            xpe(3,7,i+1,j) = -p5*HTN(i+1,j) + x1t
            ype(1,7,i+1,j) = -p5*HTW(i+1,j)
            ype(2,7,i+1,j) =  p5*HTW(i+1,j)
            ype(3,7,i+1,j) =  p5*HTW(i+1,j) + y1t
            
            ! flux contribution from triangle lying along the upper part 
            ! of W face for case of line xa-xb intersecting the face
            cnd_e(7,i+1,j) = xa*xb.lt.c0 .and. x1.gt.c0         ! E1b
            if (cnd_e(7,i+1,j)) then
               xpe(1,7,i+1,j) =  -p5 *
     &              (w1*(HTN(i+1,j)-HTS(i+1,j)) + HTS(i+1,j))
               ype(1,7,i+1,j) = (w1 - p5) * HTW(i+1,j)
            endif
            
            cnd_e(7,i+1,j) = cnd_e(7,i+1,j) .or. 
     &           (xa*xb.ge.c0 .and. xa+xb.gt.c0)

            area_e(7,i+1,j) = -triarea(xpe(1:3,7,i+1,j),
     &                                 ype(1:3,7,i+1,j))
            
            ! flux contribution from triangle touching but not
            ! lying along the W face (other part of convex quadrilateral)
            xpe(1,8,i+1,j) = -p5*HTS(i+1,j)
            xpe(2,8,i+1,j) = -p5*HTN(i+1,j) + x1t
            xpe(3,8,i+1,j) = -p5*HTS(i+1,j) + x2t
            ype(1,8,i+1,j) = -p5*HTW(i+1,j)
            ype(2,8,i+1,j) =  p5*HTW(i+1,j) + y1t
            ype(3,8,i+1,j) = -p5*HTW(i+1,j) + y2t
            
            ! flux contribution from triangle lying along the lower part 
            ! of W face for case of line xa-xb intersecting the face
            cnd_e(8,i+1,j) = xa*xb.lt.c0 .and. x2.gt.c0         ! E2b
            if (cnd_e(8,i+1,j)) then
               xpe(2,8,i+1,j) = -p5 *
     &              (w1*(HTN(i+1,j)-HTS(i+1,j)) + HTS(i+1,j))
               ype(2,8,i+1,j) = (w1 - p5) * HTW(i+1,j)
            endif

            cnd_e(8,i+1,j) = cnd_e(8,i+1,j) .or.
     &           (xa*xb.ge.c0 .and. xa+xb.gt.c0)

            area_e(8,i+1,j) = -triarea(xpe(1:3,8,i+1,j),
     &                                 ype(1:3,8,i+1,j))            

         enddo                  ! i loop
      enddo                     ! j loop

      ! coordinates of triangle points needed for flux integrals
      call triangle_coordinates (xpe, ype, cnd_e)

    !-------------------------------------------------------------------
    ! North faces
    !-------------------------------------------------------------------

      ! Initialize cnd
      do j = 1, jmt_local
         do i = 1, imt_local
            do m = 1, 8
               cnd_n(m,i,j) = .false.
            enddo
         enddo
      enddo

      do j = jlo-1, jhi  ! includes S face of cells with index jlo
         do i = ilo, ihi

            ! velocity coordinates
            x2 = -dt*u(i,j)
            y2 = -dt*v(i,j)
            x1 = -dt*u(i-1,j)
            y1 = -dt*v(i-1,j)
            w1 =  c1 / (x1 - HTN(i,j)  - x2)
            ya = (x1*y2 - y1*(HTN(i,j) + x2)) * w1
            yb = (y2*(x1 - HTN(i,j)) - x2*y1) * w1
            if (abs(ya-yb) .ge. puny) then
               xa = ya*HTN(i,j) / (ya - yb)
               xb = xa - HTN(i,j)
            else
               xa = c0
               xb = c0
            endif

            ! flux contribution from triangles in NW cell
            cnd_n(1,i-1,j+1) = ya.gt.c0 .and. x1.lt.c0          ! NW
            x1t = mse(1,1,i-1,j+1)*x1 + mse(1,2,i-1,j+1)*y1
            y1t = mse(2,1,i-1,j+1)*x1 + mse(2,2,i-1,j+1)*y1
            xat = mse(1,2,i-1,j+1)*ya
            yat = mse(2,2,i-1,j+1)*ya
            xpn(1,1,i-1,j+1) =  p5*HTS(i-1,j+1)
            xpn(2,1,i-1,j+1) =  xpn(1,1,i-1,j+1) + xat
            xpn(3,1,i-1,j+1) =  xpn(1,1,i-1,j+1) + x1t
            ypn(1,1,i-1,j+1) = -p5*HTE(i-1,j+1)
            ypn(2,1,i-1,j+1) =  ypn(1,1,i-1,j+1) + yat
            ypn(3,1,i-1,j+1) =  ypn(1,1,i-1,j+1) + y1t

            if (cnd_n(1,i-1,j+1) .and. y1.lt.c0) then           ! NW2
                                ! If y1 .ge. 0, we have NW triangle.
                                ! Similarly below for W, NE, and E.
               xat = mse(1,1,i-1,j+1)*xa
               yat = mse(2,1,i-1,j+1)*xa
               xpn(3,1,i-1,j+1) = xpn(1,1,i-1,j+1) + xat
               ypn(3,1,i-1,j+1) = ypn(1,1,i-1,j+1) + yat
            endif

            area_n(1,i-1,j+1) = -triarea(xpn(1:3,1,i-1,j+1),
     &                                   ypn(1:3,1,i-1,j+1))
            
            cnd2 = ya.lt.c0 .and. x1.lt.c0 .and. y1.gt.c0       ! NW1
            if (cnd2) then
               xat = mse(1,1,i-1,j+1)*xa
               yat = mse(2,1,i-1,j+1)*xa
               xpn(2,1,i-1,j+1) = xpn(1,1,i-1,j+1) + xat
               ypn(2,1,i-1,j+1) = ypn(1,1,i-1,j+1) + yat
               area_n(1,i-1,j+1) = triarea(xpn(1:3,1,i-1,j+1),
     &                                     ypn(1:3,1,i-1,j+1))
            endif

            cnd_n(1,i-1,j+1) = cnd_n(1,i-1,j+1) .or. cnd2

            ! flux contribution from triangles in W cell
            cnd_n(2,i-1,j) = ya.lt.c0 .and. x1.lt.c0            ! W
            x1t = mne(1,1,i-1,j)*x1 + mne(1,2,i-1,j)*y1
            y1t = mne(2,1,i-1,j)*x1 + mne(2,2,i-1,j)*y1
            xat = mne(1,2,i-1,j)*ya
            yat = mne(2,2,i-1,j)*ya
            xpn(1,2,i-1,j) =  p5*HTN(i-1,j)
            xpn(2,2,i-1,j) =  xpn(1,2,i-1,j) + xat
            xpn(3,2,i-1,j) =  xpn(1,2,i-1,j) + x1t
            ypn(1,2,i-1,j) =  p5*HTE(i-1,j)
            ypn(2,2,i-1,j) =  ypn(1,2,i-1,j) + yat
            ypn(3,2,i-1,j) =  ypn(1,2,i-1,j) + y1t

            if (cnd_n(2,i-1,j) .and. y1.gt.c0) then             ! W1
               xat = mne(1,1,i-1,j)*xa
               yat = mne(2,1,i-1,j)*xa
               xpn(3,2,i-1,j) = xpn(1,2,i-1,j) + xat
               ypn(3,2,i-1,j) = ypn(1,2,i-1,j) + yat
            endif
            
            area_n(2,i-1,j) = triarea(xpn(1:3,2,i-1,j),
     &                                ypn(1:3,2,i-1,j))

            cnd2 = ya.gt.c0 .and. x1.lt.c0 .and. y1.lt.c0       ! W2
            if (cnd2) then
               xat = mne(1,1,i-1,j)*xa
               yat = mne(2,1,i-1,j)*xa
               xpn(2,2,i-1,j) = xpn(1,2,i-1,j) + xat
               ypn(2,2,i-1,j) = ypn(1,2,i-1,j) + yat
               area_n(2,i-1,j) = -triarea(xpn(1:3,2,i-1,j),
     &                                    ypn(1:3,2,i-1,j))
            endif

            cnd_n(2,i-1,j) = cnd_n(2,i-1,j) .or. cnd2

            ! flux contribution from triangles in NE cell
            cnd_n(3,i+1,j+1) =  yb.gt.c0 .and. x2.gt.c0         ! NE
            x2t = msw(1,1,i+1,j+1)*x2 + msw(1,2,i+1,j+1)*y2
            y2t = msw(2,1,i+1,j+1)*x2 + msw(2,2,i+1,j+1)*y2
            xbt = msw(1,2,i+1,j+1)*yb
            ybt = msw(2,2,i+1,j+1)*yb
            xpn(1,3,i+1,j+1) = -p5*HTS(i+1,j+1)
            xpn(2,3,i+1,j+1) =  xpn(1,3,i+1,j+1) + xbt
            xpn(3,3,i+1,j+1) =  xpn(1,3,i+1,j+1) + x2t
            ypn(1,3,i+1,j+1) = -p5*HTW(i+1,j+1)
            ypn(2,3,i+1,j+1) =  ypn(1,3,i+1,j+1) + ybt
            ypn(3,3,i+1,j+1) =  ypn(1,3,i+1,j+1) + y2t

            if (cnd_n(3,i+1,j+1) .and. y2.lt.c0) then           ! NE2
               xbt = msw(1,1,i+1,j+1)*xb
               ybt = msw(2,1,i+1,j+1)*xb
               xpn(3,3,i+1,j+1) = xpn(1,3,i+1,j+1) + xbt
               ypn(3,3,i+1,j+1) = ypn(1,3,i+1,j+1) + ybt
            endif
            
            area_n(3,i+1,j+1) = -triarea(xpn(1:3,3,i+1,j+1),
     &                                   ypn(1:3,3,i+1,j+1))

            cnd2 = yb.lt.c0 .and. x2.gt.c0  .and. y2.gt.c0      ! NE1
            if (cnd2) then
               xbt = msw(1,1,i+1,j+1)*xb
               ybt = msw(2,1,i+1,j+1)*xb
               xpn(2,3,i+1,j+1) = xpn(1,3,i+1,j+1) + xbt
               ypn(2,3,i+1,j+1) = ypn(1,3,i+1,j+1) + ybt
               area_n(3,i+1,j+1) = triarea(xpn(1:3,3,i+1,j+1),
     &                                     ypn(1:3,3,i+1,j+1))
            endif
            
            cnd_n(3,i+1,j+1) = cnd_n(3,i+1,j+1) .or. cnd2

            ! flux contribution from triangles in E cell
            cnd_n(4,i+1,j) =  yb.lt.c0 .and. x2.gt.c0           ! E
            x2t = mnw(1,1,i+1,j)*x2 + mnw(1,2,i+1,j)*y2
            y2t = mnw(2,1,i+1,j)*x2 + mnw(2,2,i+1,j)*y2
            xbt = mnw(1,2,i+1,j)*yb
            ybt = mnw(2,2,i+1,j)*yb
            xpn(1,4,i+1,j) = -p5*HTN(i+1,j)
            xpn(2,4,i+1,j) =  xpn(1,4,i+1,j) + xbt
            xpn(3,4,i+1,j) =  xpn(1,4,i+1,j) + x2t
            ypn(1,4,i+1,j) =  p5*HTW(i+1,j)
            ypn(2,4,i+1,j) =  ypn(1,4,i+1,j) + ybt
            ypn(3,4,i+1,j) =  ypn(1,4,i+1,j) + y2t

            if (cnd_n(4,i+1,j) .and. y2.gt.c0) then             ! E1
               xbt = mnw(1,1,i+1,j)*xb
               ybt = mnw(2,1,i+1,j)*xb
               xpn(3,4,i+1,j) = xpn(1,4,i+1,j) + xbt
               ypn(3,4,i+1,j) = ypn(1,4,i+1,j) + ybt
            endif
 
            area_n(4,i+1,j) = triarea(xpn(1:3,4,i+1,j),
     &                                ypn(1:3,4,i+1,j))

            cnd2 = yb.gt.c0 .and. x2.gt.c0 .and. y2.lt.c0       ! E2
            if (cnd2) then
               xbt = mnw(1,1,i+1,j)*xb
               ybt = mnw(2,1,i+1,j)*xb
               xpn(2,4,i+1,j) = xpn(1,4,i+1,j) + xbt
               ypn(2,4,i+1,j) = ypn(1,4,i+1,j) + ybt
               area_n(4,i+1,j) = -triarea(xpn(1:3,4,i+1,j),
     &                                    ypn(1:3,4,i+1,j))
            endif
            
            cnd_n(4,i+1,j) = cnd_n(4,i+1,j) .or. cnd2

            ! flux from triangles inside home cell
            if (x1 .lt. c0) then
               x1 = c0
               y1 = ya
            endif
            
            if (x2 .gt. c0) then
               x2 = c0
               y2 = yb
            endif
            
            if (abs(yb-ya) .ge. puny) then
               w1 = min (c1, max(c0, yb/(yb-ya)))
            else
               w1 = c0
            endif
            
            x1t = mnw(1,1,i,j)*x1 + mnw(1,2,i,j)*y1
            y1t = mnw(2,1,i,j)*x1 + mnw(2,2,i,j)*y1
            x2t = mne(1,1,i,j)*x2 + mne(1,2,i,j)*y2
            y2t = mne(2,1,i,j)*x2 + mne(2,2,i,j)*y2

            ! flux contribution from triangle that includes the
            ! N cell face (part of convex quadrilateral inside home cell)
            xpn(1,5,i,j) =  p5*HTN(i,j)                         ! H1a
            xpn(2,5,i,j) = -p5*HTN(i,j)
            xpn(3,5,i,j) = -p5*HTN(i,j) + x1t
            ypn(1,5,i,j) =  p5*HTE(i,j)
            ypn(2,5,i,j) =  p5*HTW(i,j)
            ypn(3,5,i,j) =  p5*HTW(i,j) + y1t
            
            ! flux contribution from triangle lying along the left part of 
            ! the N face for case of line ya-yb intersecting the face
            cnd_n(5,i,j) = ya*yb.lt.c0 .and. y1.lt.c0           ! H1b
            if (cnd_n(5,i,j)) then
               xpn(1,5,i,j) = (p5 - w1)*HTN(i,j)
               ypn(1,5,i,j) = p5 *
     &              (w1*(HTW(i,j)-HTE(i,j)) + HTE(i,j))
            endif
            
            cnd_n(5,i,j) = cnd_n(5,i,j) .or. 
     &           (ya*yb.ge.c0 .and. ya+yb.lt.c0)

            area_n(5,i,j) = triarea(xpn(1:3,5,i,j),
     &                              ypn(1:3,5,i,j))

            ! flux contribution from triangle touching but not lying 
            ! along the N face (other part of convex quadrilateral)
            xpn(1,6,i,j) =  p5*HTN(i,j)                         ! H2a
            xpn(2,6,i,j) = -p5*HTN(i,j) + x1t
            xpn(3,6,i,j) =  p5*HTN(i,j) + x2t
            ypn(1,6,i,j) =  p5*HTE(i,j)
            ypn(2,6,i,j) =  p5*HTW(i,j) + y1t
            ypn(3,6,i,j) =  p5*HTE(i,j) + y2t
            
            ! flux contribution from triangle lying along the right part 
            ! of the N face for case of line ya-yb intersecting the face
            cnd_n(6,i,j) = ya*yb.lt.c0 .and. y2.lt.c0           ! H2b
            if (cnd_n(6,i,j)) then
               xpn(2,6,i,j) = (p5 - w1)*HTN(i,j)
               ypn(2,6,i,j) = p5 *
     &              (w1*(HTW(i,j)-HTE(i,j)) + HTE(i,j))
            endif
            
            cnd_n(6,i,j) = cnd_n(6,i,j) .or. 
     &           (ya*yb.ge.c0 .and. ya+yb.lt.c0)
      
            area_n(6,i,j) = triarea(xpn(1:3,6,i,j),
     &                              ypn(1:3,6,i,j))

            ! flux from triangles in N cell
            x1t = msw(1,1,i,j+1)*x1 + msw(1,2,i,j+1)*y1
            y1t = msw(2,1,i,j+1)*x1 + msw(2,2,i,j+1)*y1
            x2t = mse(1,1,i,j+1)*x2 + mse(1,2,i,j+1)*y2
            y2t = mse(2,1,i,j+1)*x2 + mse(2,2,i,j+1)*y2

            ! flux contribution from triangle that includes the
            ! S cell face (part of convex quadrilateral inside home cell)
            xpn(1,7,i,j+1) =  p5*HTS(i,j+1)                     ! N1a
            xpn(2,7,i,j+1) = -p5*HTS(i,j+1)
            xpn(3,7,i,j+1) = -p5*HTS(i,j+1) + x1t
            ypn(1,7,i,j+1) = -p5*HTE(i,j+1)
            ypn(2,7,i,j+1) = -p5*HTW(i,j+1)
            ypn(3,7,i,j+1) = -p5*HTW(i,j+1) + y1t
            
            ! flux contribution from triangle lying along the left part 
            ! of the S face for case of line ya-yb intersecting the face
            cnd_n(7,i,j+1) = ya*yb.lt.c0 .and. y1.gt.c0         ! N1b
            if (cnd_n(7,i,j+1)) then
               xpn(1,7,i,j+1) = (p5 - w1)*HTS(i,j+1)
               ypn(1,7,i,j+1) = -p5 *
     &              (w1*(HTW(i,j+1)-HTE(i,j+1)) + HTE(i,j+1)) 
            endif
            
            cnd_n(7,i,j+1) = cnd_n(7,i,j+1) .or. 
     &           (ya*yb.ge.c0 .and. ya+yb.gt.c0)

            area_n(7,i,j+1) = -triarea(xpn(1:3,7,i,j+1),
     &                                 ypn(1:3,7,i,j+1))

            ! flux contribution from triangle touching but not
            ! lying along the S face (other part of convex quadrilateral)
            xpn(1,8,i,j+1) =  p5*HTS(i,j+1)
            xpn(2,8,i,j+1) = -p5*HTS(i,j+1) + x1t
            xpn(3,8,i,j+1) =  p5*HTS(i,j+1) + x2t
            ypn(1,8,i,j+1) = -p5*HTE(i,j+1)
            ypn(2,8,i,j+1) = -p5*HTW(i,j+1) + y1t
            ypn(3,8,i,j+1) = -p5*HTE(i,j+1) + y2t
 
            ! flux contribution from triangle lying along the right part
            ! of the S face for case of line ya-yb intersecting the face
            cnd_n(8,i,j+1) = ya*yb.lt.c0 .and. y2.gt.c0         ! N2b
            if (cnd_n(8,i,j+1)) then
               xpn(2,8,i,j+1) = (p5 - w1)*HTS(i,j+1)
               ypn(2,8,i,j+1) = -p5 *
     &              (w1*(HTW(i,j+1)-HTE(i,j+1)) + HTE(i,j+1))
            endif
            
            cnd_n(8,i,j+1) = cnd_n(8,i,j+1) .or. 
     &           (ya*yb.ge.c0 .and. ya+yb.gt.c0)

            area_n(8,i,j+1) = -triarea(xpn(1:3,8,i,j+1),
     &                                 ypn(1:3,8,i,j+1))

         enddo                  ! i loop
      enddo                     ! j loop

    !-------------------------------------------------------------------
    ! Compute velocity divergence implied by remapping (used for
    ! mechanical redistribution).
    !-------------------------------------------------------------------

      ! east face fluxes of total area (ice plus open water)
      call area_flux (area_e, cnd_e, ishift_e, jshift_e, areaflxe)

      ! north face fluxes
      call area_flux (area_n, cnd_n, ishift_n, jshift_n, areaflxn)

      do j = jlo, jhi
         do i = ilo, ihi
            w1 = areaflxe(i,j) - areaflxe(i-1,j) 
     &         + areaflxn(i,j) - areaflxn(i,j-1)  
            divu_adv(i,j) = w1*tarear(i,j)/dt
         enddo
      enddo

      ! coordinates of triangle points needed for flux integrals
      call triangle_coordinates (xpn, ypn, cnd_n)

    !-------------------------------------------------------------------
    ! Part 2: Transfer fluxes and update state variables.
    !     This is done for one thickness category at a time.
    !-------------------------------------------------------------------

    !-------------------------------------------------------------------
    ! Compute limited gradients of ice area.  (See comments below.)
    !-------------------------------------------------------------------

      do n = 1, ncat

         ! boundary values
         call bound (aicen(:,:,n))

         ! limited gradients
         call limited_gradient (aixn(:,:,n), aiyn(:,:,n), aicen(:,:,n), 
     &        cx, cy, hm)      

      enddo

    !-------------------------------------------------------------------
    ! Limit the gradients further to ensure aice(x,y) <= 1 everywhere.
    ! (aice is the total fractional ice area)
    !-------------------------------------------------------------------

      do j = 1, jmt_local
         do i = 1, imt_local

            if (tmask(i,j)) then
               ! total ice area and gradients
               aixsum = c0
               aiysum = c0
               aicesum = c0
               do n = 1, ncat
                  aicesum = aicesum + aicen(i,j,n) 
                  aixsum  = aixsum + aixn(i,j,n)
                  aiysum  = aiysum + aiyn(i,j,n)
               enddo

              if (aicesum .gt. puny) then
                 ! total ice area in cell corners
                 w1 = aicesum + ( p5*HTN(i,j) - cx(i,j))*aixsum ! NE 
     &                        + ( p5*HTE(i,j) - cy(i,j))*aiysum
                 w2 = aicesum + (-p5*HTN(i,j) - cx(i,j))*aixsum ! NW
     &                        + ( p5*HTW(i,j) - cy(i,j))*aiysum
                 w3 = aicesum + (-p5*HTS(i,j) - cx(i,j))*aixsum ! SW
     &                        + (-p5*HTW(i,j) - cy(i,j))*aiysum
                 w4 = aicesum + ( p5*HTS(i,j) - cx(i,j))*aixsum ! SE
     &                        + (-p5*HTE(i,j) - cy(i,j))*aiysum
                 w1 = max (w1, w2, w3, w4)

                 if (w1 .gt. c1+puny) then     ! limit aixn and aiyn
                    w2 = (c1 - aicesum) / (w1 - aicesum)
                    do n = 1, ncat
                       aixn(i,j,n) = aixn(i,j,n)*w2
                       aiyn(i,j,n) = aiyn(i,j,n)*w2
                    enddo
                 endif          ! aimax > 1
              endif             ! aicesum > puny
           endif                ! tmask
         enddo                  ! i
      enddo                     ! j

    !-------------------------------------------------------------------
    ! Main loop over categories
    !-------------------------------------------------------------------

      do n = 1, ncat

    !-------------------------------------------------------------------
    ! boundary values of state variables
    !-------------------------------------------------------------------

         call bound(Tsfcn(:,:,n))
         call bound(vicen(:,:,n))
         call bound(vsnon(:,:,n))
         do k = (n-1)*nilyr+1, n*nilyr
            call bound(eicen(:,:,k))
         enddo
         call bound(esnon(:,:,n))

    !-------------------------------------------------------------------
    ! Compute hice, hsno, qice, qsno in each grid cell.
    ! Note that hice(i,j) is the ice thickness in category n; it is
    ! not averaged over categories like vice(i,j).  The same
    ! convention applies to hsno, qice, qsno, and other variables
    ! inside the n loop.  The n index is omitted to save memory.
    ! Also compute ice and snow masks.
    !-------------------------------------------------------------------

         do j = 1, jmt_local
            do i = 1, imt_local
 
               if (aicen(i,j,n).gt.puny) then

                  iceflag(i,j) = c1
                  hice(i,j) = vicen(i,j,n)/aicen(i,j,n)
                  do k = 1, nilyr
                     qice(i,j,k) = eicen(i,j,(n-1)*nilyr+k)/vicen(i,j,n)
                  enddo
                  if (vsnon(i,j,n)/aicen(i,j,n).gt.puny) then
                     snoflag(i,j) = c1
                     hsno(i,j) = vsnon(i,j,n)/aicen(i,j,n)		
                     qsno(i,j) = esnon(i,j,n)/vsnon(i,j,n)
                  else
                     snoflag(i,j) = c0
                     hsno(i,j) = c0
                     qsno(i,j) = c0
                  endif

               else             ! aicen < puny

                  iceflag(i,j) = c0
                  snoflag(i,j) = c0
                  hice(i,j) = c0
                  do k = 1, nilyr
                     qice(i,j,k) = c0
                  enddo
                  hsno(i,j) = c0
                  qsno(i,j) = c0

               endif            ! aicen > puny

            enddo               ! i
         enddo                  ! j

    !-------------------------------------------------------------------
    ! Compute field values (aicen, Tsfcn, hice, hsno, qice, qsno) at the 
    ! geometric center of each grid cell, and compute limited gradients 
    ! in the x and y directions.
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
    ! of the line joining the midpoints of the north and south faces
    ! with the line joining the midpoints of the east and west faces.  
    ! To find the value at the geometric center, we must know the 
    ! location of the cell centroid/center of ice area/center of ice 
    ! or snow volume relative to the geometric center, along with 
    ! the field gradients with respect to x and y. 
    !
    ! The cell gradients are first computed from the difference between 
    ! values in the neighboring cells, then limited by requiring that 
    ! no new extrema are created within the cell.
    !-------------------------------------------------------------------

         ! ice area 

         do j = 1, jmt_local
            do i = 1, imt_local

               if (iceflag(i,j).gt.p5) then  ! ice is present
                 ! ice area at geometric center
                 aigc(i,j) = aicen(i,j,n) - cx(i,j)*aixn(i,j,n) 
     &                                    - cy(i,j)*aiyn(i,j,n)
                 ! center of ice area (cax,cay) for each cell
                 cax(i,j) = (aixn(i,j,n)*cxx(i,j) + aiyn(i,j,n)*cxy(i,j) 
     &                    + aigc(i,j)*cx(i,j)) / aicen(i,j,n)
                 cay(i,j) = (aixn(i,j,n)*cxy(i,j) + aiyn(i,j,n)*cyy(i,j) 
     &                    + aigc(i,j)*cy(i,j)) / aicen(i,j,n)
               else
                 aigc(i,j) = c0
                 cax(i,j)  = c0
                 cay(i,j)  = c0
               endif
            enddo               ! i
         enddo                  ! j

         ! surface temperature, ice thickness, and snow thickness
         call limited_gradient (Tsfx, Tsfy, Tsfcn(:,:,n), cax, cay, 
     &        iceflag)

         call limited_gradient (hix, hiy, hice, cax, cay, iceflag)

         call limited_gradient (hsx, hsy, hsno, cax, cay, iceflag)

         do j = 1, jmt_local
          do i = 1, imt_local

             if (iceflag(i,j).gt.p5) then  ! ice is present

                ! surface temperature at geometric center
                Tsfgc(i,j) = Tsfcn(i,j,n) - Tsfx(i,j)*cax(i,j) 
     &                                    - Tsfy(i,j)*cay(i,j)

                ! ice thickness
                ! geometric center
                higc(i,j) = hice(i,j) - hix(i,j)*cax(i,j) 
     &                                - hiy(i,j)*cay(i,j)
                ! center of ice thickness
                w1 = aigc(i,j)*higc(i,j)
                w2 = aigc(i,j)*hix(i,j) + aixn(i,j,n)*higc(i,j)
                w3 = aigc(i,j)*hiy(i,j) + aiyn(i,j,n)*higc(i,j)
                w4 = aixn(i,j,n)*hix(i,j)
                w5 = aixn(i,j,n)*hiy(i,j) + aiyn(i,j,n)*hix(i,j)
                w6 = aiyn(i,j,n)*hiy(i,j)
                w7 = c1 / vicen(i,j,n)
                cvix(i,j) = (w1*cx(i,j) + w2*cxx(i,j)  + w3*cxy(i,j)
     &                 + w4*cxxx(i,j) + w5*cxxy(i,j) + w6*cxyy(i,j))
     &                 * w7
                cviy(i,j) = (w1*cy(i,j) + w2*cxy(i,j)  + w3*cyy(i,j)
     &                 + w4*cxxy(i,j) + w5*cxyy(i,j) + w6*cyyy(i,j))
     &                 * w7

                ! snow thickness
                if (snoflag(i,j).gt.p5) then  ! snow is present
                   ! geometric center
                   hsgc(i,j) = hsno(i,j) - hsx(i,j)*cax(i,j) 
     &                                   - hsy(i,j)*cay(i,j)
                   ! center of snow thickness
                   w1 = aigc(i,j)*hsgc(i,j)
                   w2 = aigc(i,j)*hsx(i,j) + aixn(i,j,n)*hsgc(i,j)
                   w3 = aigc(i,j)*hsy(i,j) + aiyn(i,j,n)*hsgc(i,j)
                   w4 = aixn(i,j,n)*hsx(i,j)
                   w5 = aixn(i,j,n)*hsy(i,j) + aiyn(i,j,n)*hsx(i,j)
                   w6 = aiyn(i,j,n)*hsy(i,j)
                   w7 = c1 / vsnon(i,j,n)
                   cvsx(i,j) = (w1*cx(i,j) + w2*cxx(i,j)  + w3*cxy(i,j)
     &                    + w4*cxxx(i,j) + w5*cxxy(i,j) + w6*cxyy(i,j))
     &                    * w7
                   cvsy(i,j) = (w1*cy(i,j) + w2*cxy(i,j)  + w3*cyy(i,j)
     &                    + w4*cxxy(i,j) + w5*cxyy(i,j) + w6*cyyy(i,j))
     &                    * w7
                else            ! no snow
                   hsgc(i,j) = c0
                   cvsx(i,j) = c0
                   cvsy(i,j) = c0
                endif           ! snoflag

             else               ! no ice
                Tsfgc(i,j) = c0
                higc(i,j)  = c0
                cvix(i,j)  = c0
                cviy(i,j)  = c0
                hsgc(i,j)  = c0
                cvsx(i,j)  = c0
                cvsy(i,j)  = c0
             endif              ! iceflag

          enddo
         enddo

         ! ice and snow enthalpy
         do k = 1, nilyr
            call limited_gradient (qix(:,:,k), qiy(:,:,k), qice(:,:,k), 
     &           cvix, cviy, iceflag)
         enddo

         call limited_gradient (qsx, qsy, qsno, cvsx, cvsy, snoflag)

         do j = 1, jmt_local
            do i = 1, imt_local

               if (iceflag(i,j).gt.p5) then  ! ice is present
                  ! ice enthalpy
                  do k = 1, nilyr
                     qigc(i,j,k) = qice(i,j,k) - qix(i,j,k)*cvix(i,j) 
     &                                         - qiy(i,j,k)*cviy(i,j)
                  enddo
                  if (snoflag(i,j).gt.p5) then  ! snow is present
                  ! snow enthalpy
                     qsgc(i,j) = qsno(i,j) - qsx(i,j)*cvsx(i,j) 
     &                                     - qsy(i,j)*cvsy(i,j)
                  else
                     qsgc(i,j) = c0
                  endif

               else             ! no ice
                  do k = 1, nilyr
                     qigc(i,j,k) = c0
                  enddo
                  qsgc(i,j) = c0
               endif

            enddo               ! i
         enddo                  ! j

    !-------------------------------------------------------------------
    ! Compute the flux of ice area and tracers across the cell faces by 
    ! summing contributions from each triangle.
    !-------------------------------------------------------------------

         ! east faces
         call sum_fluxes (iceflag, 
     &     xpe, ype, area_e, cnd_e, ishift_e, jshift_e,
     &     aigc, aixn(:,:,n), aiyn(:,:,n), aiflxe, 
     &     Tsfgc, Tsfx, Tsfy, aTflxe,  
     &     higc, hix, hiy, viflxe, hsgc, hsx, hsy, vsflxe, 
     &     qigc, qix, qiy, eiflxe, qsgc, qsx, qsy, esflxe)     

         ! north faces
         call sum_fluxes (iceflag, 
     &     xpn, ypn, area_n, cnd_n, ishift_n, jshift_n,
     &     aigc, aixn(:,:,n), aiyn(:,:,n), aiflxn, 
     &     Tsfgc, Tsfx, Tsfy, aTflxn,
     &     higc, hix, hiy, viflxn, hsgc, hsx, hsy, vsflxn, 
     &     qigc, qix, qiy, eiflxn, qsgc, qsx, qsy, esflxn) 

    !-------------------------------------------------------------------
    ! update the state variables
    ! aicen, Tsfcn, vicen, vsnon, eicen, esnon
    !-------------------------------------------------------------------
         do j=jlo,jhi
            do i=ilo,ihi

               aTsfold = aicen(i,j,n)*Tsfcn(i,j,n)  ! save starting value

               ! ice area
               w1 = aiflxe(i,j) - aiflxe(i-1,j)
     &            + aiflxn(i,j) - aiflxn(i,j-1)
               aicen(i,j,n) = aicen(i,j,n) - w1*tarear(i,j)

               if (aicen(i,j,n) .le. -puny) then  ! there's a bug
                  write (nu_diag,*) ''
                  write (nu_diag,*) 'New area < 0, istep =', istep1
                  write (nu_diag,*) '(my_task,i,j,n) =', my_task,i,j,n
                  write (nu_diag,*) 'New area =', aicen(i,j,n)
                  write (nu_diag,*) 'Net flux =', -w1*tarear(i,j)
                  stop

               elseif (aicen(i,j,n) .le. puny) then
                  aicen(i,j,n) = c0
                  Tsfcn(i,j,n) = c0
                  vicen(i,j,n) = c0
                  vsnon(i,j,n) = c0
                  do k = 1, nilyr
                     eicen(i,j,(n-1)*nilyr+k) = c0
                  enddo
                  esnon(i,j,n) = c0

               else  ! aicen(i,j,n) > puny

                  ! surface temperature
                  w1  = aTflxe(i,j) - aTflxe(i-1,j)
     &                + aTflxn(i,j) - aTflxn(i,j-1)
                  Tsfcn(i,j,n) = (aTsfold - w1*tarear(i,j))
     &                          / aicen(i,j,n)
                        
                  ! ice volume
                  w1  = viflxe(i,j) - viflxe(i-1,j)
     &                + viflxn(i,j) - viflxn(i,j-1)
                  vicen(i,j,n) = vicen(i,j,n) - w1*tarear(i,j)

                  ! snow volume
                  w1  = vsflxe(i,j) - vsflxe(i-1,j)
     &                + vsflxn(i,j) - vsflxn(i,j-1)
                  vsnon(i,j,n) = vsnon(i,j,n) - w1*tarear(i,j)

                  ! ice energy
                  do k = 1, nilyr
                     w1 = eiflxe(i,j,k) - eiflxe(i-1,j,k)
     &                  + eiflxn(i,j,k) - eiflxn(i,j-1,k)
                     eicen(i,j,(n-1)*nilyr+k) = 
     &            	eicen(i,j,(n-1)*nilyr+k)	- w1*tarear(i,j)
                  enddo

                  ! snow energy
                  w1 = esflxe(i,j) - esflxe(i-1,j)
     &               + esflxn(i,j) - esflxn(i,j-1)
                  esnon(i,j,n) = esnon(i,j,n) - w1*tarear(i,j)

               endif            ! aicen(i,j,n)

            enddo               ! i loop
         enddo                  ! j loop

      enddo                     ! n loop

    !-------------------------------------------------------------------
    ! Check global conservation of area, volume, and energy
    !-------------------------------------------------------------------

c      call conserved_sums (aicegsum2, vicegsum2, vsnogsum2,
c     &     eicegsum2, esnogsum2)

c      if (my_task.eq.master_task) then
c         call global_conservation (aicegsum1, aicegsum2,
c     &        vicegsum1, vicegsum2, vsnogsum1, vsnogsum2, 
c     &        eicegsum1, eicegsum2, esnogsum1, esnogsum2)
c      endif

      call ice_timer_stop(3)  ! advection

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

c=======================================================================

      subroutine triangle_coordinates (xp, yp, cnd)

!---!-------------------------------------------------------------------
!---! For each triangle, find the coordinates of the 4 points needed
!---! to compute integrals of cubic polynomials, using a formula from 
!---! A.H. Stroud, Approximate Calculation of Multiple Integrals, 
!---! Prentice-Hall, 1971.  (Section 8.8, formula 3.1.)  
!---! Quadratic functions can be integrated using 3-point formulas, but 
!---! it is more efficient to use a single formula for both quadratics 
!---! and cubics.
!---!
!---! The formula is as follows:
!---! I3 = integral of f(x,y)*dA  
!---!    = AR * [ -9/16 *  f(x0,y0) 
!---!           + 25/48 * (f(x1,y1) + f(x2,y2) + f(x3,y3))]
!---! where I3 is the integral of a polynomial of 3rd order or lower,
!---! AR is the area of the triangle,  (x0,y0) is the midpoint,
!---! and the other three points are located 2/5 of the way from 
!---! the midpoint to each of the three vertices.
!---!-------------------------------------------------------------------

      real (kind=dbl_kind), intent(inout), 
     & dimension (0:3, 8, imt_local, jmt_local) ::
     &     xp, yp

      logical (kind=log_kind), intent(in) ::
     &     cnd (8, imt_local, jmt_local)

      integer (kind=int_kind) :: 
     &     i, j, m            ! loop indices

      do j = 1, jmt_local
         do i = 1, imt_local
            do m = 1, 8
               if (cnd(m,i,j)) then

                  ! coordinates of midpoint
                  xp(0,m,i,j) = p333 * 
     &                 (xp(1,m,i,j) + xp(2,m,i,j) + xp(3,m,i,j))
                  yp(0,m,i,j) = p333 * 
     &                 (yp(1,m,i,j) + yp(2,m,i,j) + yp(3,m,i,j))

                  ! coordinates of the other 3 points needed for integrals
                  ! (2/5 of the way from the midpoint to the 3 vertices)
                  xp(1,m,i,j) = p4*xp(1,m,i,j) + p6*xp(0,m,i,j)
                  xp(2,m,i,j) = p4*xp(2,m,i,j) + p6*xp(0,m,i,j)
                  xp(3,m,i,j) = p4*xp(3,m,i,j) + p6*xp(0,m,i,j)
                  yp(1,m,i,j) = p4*yp(1,m,i,j) + p6*yp(0,m,i,j)
                  yp(2,m,i,j) = p4*yp(2,m,i,j) + p6*yp(0,m,i,j)
                  yp(3,m,i,j) = p4*yp(3,m,i,j) + p6*yp(0,m,i,j)

               endif            ! cnd
            enddo               ! m (triangles 
         enddo                  ! i 
      enddo                     ! j 

      end subroutine triangle_coordinates


c=======================================================================

      subroutine area_flux (area, cnd, ishift, jshift, areaflx)

!---!-------------------------------------------------------------------
!---! Compute fluxes of total area (ice plus open water).
!---! Input variables have the same meanings as in the main subroutine.
!---!-------------------------------------------------------------------
 
      integer (kind=int_kind), intent(in) :: 
     &     ishift(8), jshift(8)

      real (kind=dbl_kind), intent(in) ::
     &     area (8, imt_local, jmt_local)

      logical (kind=log_kind), intent(in) ::
     &     cnd (8, imt_local, jmt_local)

      real (kind=dbl_kind), intent(out),
     & dimension (imt_local, jmt_local) :: 
     &     areaflx                     

      integer (kind=int_kind) :: 
     &     i, j, m            ! loop indices

      ! initialize fluxes
      do j = 1, jhi           ! don't need fluxes for jmt_local, imt_local,
         do i = 1, ihi        ! but need E flux for i=1, N flux for j=1
            areaflx(i,j) = c0
         enddo
      enddo   

    !-------------------------------------------------------------------
    ! Main loop
    !-------------------------------------------------------------------

      do j = 1, jmt_local
         do i = 1, imt_local
            do m = 1, 8
               if (cnd(m,i,j)) then
                  areaflx(i-ishift(m),j-jshift(m)) =
     &                 areaflx(i-ishift(m),j-jshift(m)) + area(m,i,j)
               endif
            enddo
         enddo
      enddo

      end subroutine area_flux

c=======================================================================

      subroutine sum_fluxes 
     &     (iceflag, xp, yp, area, cnd, ishift, jshift,
     &     aigc, aix, aiy, aiflx, Tsfgc, Tsfx, Tsfy, aTflx,  
     &     higc, hix, hiy, viflx, hsgc, hsx, hsy, vsflx, 
     &     qigc, qix, qiy, eiflx, qsgc, qsx, qsy, esflx) 
 
!---!-------------------------------------------------------------------
!---! Compute the fluxes across each face by integrating the ice area 
!---! and tracers over each flux triangle.
!---! Input variables have the same meanings as in the main subroutine.
!---! Repeated use of certain sums makes the calculation more efficient.
!---!-------------------------------------------------------------------

      integer (kind=int_kind), intent(in) :: 
     &     ishift(8), jshift(8)

      real (kind=dbl_kind), intent(inout), 
     & dimension (0:3, 8, imt_local, jmt_local) ::
     &     xp, yp          

      real (kind=dbl_kind), intent(in) ::
     &     area (8, imt_local, jmt_local)

      real (kind=dbl_kind), intent(in),
     & dimension (imt_local, jmt_local) :: 
     &     iceflag,
     &     aigc, aix, aiy, 
     &     Tsfgc,Tsfx,Tsfy,
     &     higc, hix, hiy, 
     &     hsgc, hsx, hsy,
     &     qsgc, qsx, qsy

      real (kind=dbl_kind), intent(in),
     & dimension (imt_local, jmt_local, nilyr) :: 
     &     qigc, qix, qiy

      real (kind=dbl_kind), intent(out),
     & dimension (imt_local, jmt_local) :: 
     &     aiflx, aTflx, viflx, vsflx, esflx                     
     
      real (kind=dbl_kind), intent(out),
     & dimension (imt_local, jmt_local, nilyr) :: 
     &     eiflx

      logical (kind=log_kind), intent(in) ::
     &     cnd (8, imt_local, jmt_local)

      integer (kind=int_kind) :: 
     &      i, j, k, m               ! loop indices

      real (kind=dbl_kind) ::
     &      a0, a1, a2, a3           ! ice area at internal points
     &,     asum, axsum, aysum       ! sums of a, a*x, and a*y
     &,     aTsum, ahisum, ahssum    ! sums of a*Tsfc, a*hice, a*hsno
     &,     axxsum, axysum, ayysum   ! sums of a*x*x, a*x*y, a*y*y
     &,     axhisum, ayhisum         ! sums of a*x*hice, a*y*hice
     &,     axhssum, ayhssum         ! sums of a*x*hsno, a*y*hsno
     &,     ahqisum, ahqssum         ! sums of a*hice*qice, a*hsno*qsno
 
      ! initialize fluxes
      do j = 1, jhi           ! don't need fluxes for jmt_local, imt_local,
         do i = 1, ihi        ! but need E flux for i=1, N flux for j=1
            aiflx(i,j) = c0
            aTflx(i,j) = c0
            viflx(i,j) = c0
            vsflx(i,j) = c0
            esflx(i,j) = c0
         enddo
      enddo   

      do k = 1, nilyr
         do j = 1, jhi        ! don't need fluxes for jmt_local, imt_local,
            do i = 1, ihi     ! but need E flux for i=1, N flux for j=1
               eiflx(i,j,k) = c0
            enddo
         enddo
      enddo

    !-------------------------------------------------------------------
    ! Main loop
    !-------------------------------------------------------------------

      do j = 1, jmt_local
         do i = 1, imt_local
            do m = 1, 8
               if (cnd(m,i,j) .and. iceflag(i,j).gt.p5) then

                  ! area fluxes  
                  ! find the ice area at the midpoint
                  a0 = aigc(i,j) + xp(0,m,i,j)*aix(i,j) 
     &                           + yp(0,m,i,j)*aiy(i,j)
                  aiflx(i-ishift(m),j-jshift(m)) =
     &                 aiflx(i-ishift(m),j-jshift(m)) 
     &                 + area(m,i,j) * a0

                  ! quantities needed for volume fluxes
                  ! Weighting factors are incorporated into the ice
                  ! area terms a0, a1, a2, and a3.
                  a0 = p5625m * a0 
                  a1 = p52083 * (aigc(i,j) + xp(1,m,i,j)*aix(i,j) 
     &                                     + yp(1,m,i,j)*aiy(i,j))   
                  a2 = p52083 * (aigc(i,j) + xp(2,m,i,j)*aix(i,j) 
     &                                     + yp(2,m,i,j)*aiy(i,j))
                  a3 = p52083 * (aigc(i,j) + xp(3,m,i,j)*aix(i,j) 
     &                                     + yp(3,m,i,j)*aiy(i,j)) 
                  asum = a0 + a1 + a2 + a3
                  axsum = a0*xp(0,m,i,j) + a1*xp(1,m,i,j) 
     &                  + a2*xp(2,m,i,j) + a3*xp(3,m,i,j)
                  aysum = a0*yp(0,m,i,j) + a1*yp(1,m,i,j) 
     &                  + a2*yp(2,m,i,j) + a3*yp(3,m,i,j)

                  ! surface temperature fluxes
                  aTsum = asum*Tsfgc(i,j) + axsum*Tsfx(i,j) 
     &                                    + aysum*Tsfy(i,j)
                  aTflx(i-ishift(m),j-jshift(m)) = 
     &                 aTflx(i-ishift(m),j-jshift(m))
     &                 + area(m,i,j) * aTsum

                  ! ice volume fluxes
                  ahisum = asum*higc(i,j) + axsum*hix(i,j) 
     &                                    + aysum*hiy(i,j)
                  viflx(i-ishift(m),j-jshift(m)) = 
     &                 viflx(i-ishift(m),j-jshift(m)) 
     &                 + area(m,i,j) * ahisum

                  ! snow volume fluxes
                  ahssum = asum*hsgc(i,j) + axsum*hsx(i,j) 
     &                                    + aysum*hsy(i,j)
                  vsflx(i-ishift(m),j-jshift(m)) = 
     &                 vsflx(i-ishift(m),j-jshift(m))
     &                 + area(m,i,j) * ahssum

                  ! quantities needed for enthalpy fluxes
                  axxsum = a0*xp(0,m,i,j)*xp(0,m,i,j)
     &                   + a1*xp(1,m,i,j)*xp(1,m,i,j)
     &                   + a2*xp(2,m,i,j)*xp(2,m,i,j)
     &                   + a3*xp(3,m,i,j)*xp(3,m,i,j) 
                  axysum = a0*xp(0,m,i,j)*yp(0,m,i,j)
     &                   + a1*xp(1,m,i,j)*yp(1,m,i,j)
     &                   + a2*xp(2,m,i,j)*yp(2,m,i,j)
     &                   + a3*xp(3,m,i,j)*yp(3,m,i,j) 
                  ayysum = a0*yp(0,m,i,j)*yp(0,m,i,j)
     &                   + a1*yp(1,m,i,j)*yp(1,m,i,j)
     &                   + a2*yp(2,m,i,j)*yp(2,m,i,j)
     &                   + a3*yp(3,m,i,j)*yp(3,m,i,j)

                  ! ice enthalpy fluxes
                  axhisum = axsum*higc(i,j) + axxsum*hix(i,j) 
     &                                      + axysum*hiy(i,j)
                  ayhisum = aysum*higc(i,j) + axysum*hix(i,j) 
     &                                      + ayysum*hiy(i,j)

                  do k = 1,nilyr
                     ahqisum = ahisum*qigc(i,j,k) + axhisum*qix(i,j,k) 
     &                                            + ayhisum*qiy(i,j,k)
                     eiflx(i-ishift(m),j-jshift(m),k) =
     &                    eiflx(i-ishift(m),j-jshift(m),k) +
     &                    area(m,i,j) * ahqisum
                  enddo

                  ! snow enthalpy fluxes
                  axhssum = axsum*hsgc(i,j) + axxsum*hsx(i,j) 
     &                                      + axysum*hsy(i,j)
                  ayhssum = aysum*hsgc(i,j) + axysum*hsx(i,j) 
     &                                      + ayysum*hsy(i,j)
                  
                  ahqssum = ahssum*qsgc(i,j) + axhssum*qsx(i,j) 
     &                                       + ayhssum*qsy(i,j)
                  esflx(i-ishift(m),j-jshift(m)) =
     &                 esflx(i-ishift(m),j-jshift(m)) +
     &                 area(m,i,j) * ahqssum
                  
               endif            ! triangle contributes to flux
            enddo               ! m (triangles)
         enddo                  ! i
      enddo                     ! j

      end subroutine sum_fluxes

c=======================================================================

      subroutine limited_gradient (gx, gy, phi, cnx, cny, iflag)
 
!---!-------------------------------------------------------------------
!---! Compute a limited gradient of the scalar field phi.
!---!
!---! "Limited" means that we do not create new extrema in phi.  For 
!---! instance, field values at the cell corners can neither exceed the 
!---! maximum of phi(i,j) in the cell and its eight neighbors, nor fall
!---! below the minimum.
!---!-------------------------------------------------------------------
 

      real (kind=dbl_kind), intent(in),
     & dimension (imt_local,jmt_local) ::
     &    phi     ! input tracer field (mean values in each grid cell)
     &,   cnx     ! x-coordinate of phi relative to geometric center of cell
     &,   cny     ! y-coordinate of phi relative to geometric center of cell
     &,   iflag
          ! iflag(i,j) = 1 if phi(i,j) has physical meaning, = 0 otherwise.
          ! For instance, aice has no physical meaning on land points,
          ! and hice no physical meaning where aice = 0.

      real (kind=dbl_kind), intent(out), 
     & dimension (imt_local,jmt_local) :: 
     &    gx      ! limited x-direction gradient
     &,   gy      ! limited y-direction gradient

      integer (kind=int_kind) :: i, j

      real (kind=dbl_kind) ::
     &    phi_nw, phi_n, phi_ne   ! values of phi in 8 neighbor cells 
     &,   phi_w,         phi_e
     &,   phi_sw, phi_s, phi_se
     &,   qmn, qmx       ! min and max value of phi within grid cell
     &,   pmn, pmx       ! min and max value of phi among neighbor cells
     &,   w1, w2, w3, w4 ! work variables

      do j = jlo, jhi 
         do i = ilo, ihi 

            if (iflag(i,j).gt.p5) then

               ! Store values of phi in the 8 neighbor cells.
               ! Note: iflag = c1 or c0.  If iflag = c1, use the true value;
               ! if iflag = c0, use the home cell value so that 
               ! non-physical values of phi do not contribute to the 
               ! gradient.
               phi_nw = iflag(i-1,j+1) * phi(i-1,j+1) 
     &             + (1-iflag(i-1,j+1))* phi(i,j)
               phi_n  = iflag(i,j+1)   * phi(i,j+1) 
     &             + (1-iflag(i,j+1))  * phi(i,j)
               phi_ne = iflag(i+1,j+1) * phi(i+1,j+1) 
     &             + (1-iflag(i+1,j+1))* phi(i,j)
               phi_w  = iflag(i-1,j)   * phi(i-1,j) 
     &             + (1-iflag(i-1,j))  * phi(i,j)
               phi_e  = iflag(i+1,j)   * phi(i+1,j) 
     &             + (1-iflag(i+1,j))  * phi(i,j)
               phi_sw = iflag(i-1,j-1) * phi(i-1,j-1) 
     &             + (1-iflag(i-1,j-1))* phi(i,j)
               phi_s  = iflag(i,j-1)   * phi(i,j-1) 
     &             + (1-iflag(i,j-1))  * phi(i,j)
               phi_se = iflag(i+1,j-1) * phi(i+1,j-1) 
     &             + (1-iflag(i+1,j-1))* phi(i,j)
               
               ! the unlimited gradient components
               ! grid lengths computed in ice_grid.F for efficiency:
               !       dxt2r(i,j) = 1 / (dxt(i+1,j) + dxt(i,j))
               !       dyt2r(i,j) = 1 / (dyt(i,j+1) + dyt(i,j))
               gx(i,j) = (phi_e - phi(i,j)) * dxt2r(i,j) 
     &                 + (phi(i,j) - phi_w) * dxt2r(i-1,j)
               gy(i,j) = (phi_n - phi(i,j)) * dyt2r(i,j) 
     &                 + (phi(i,j) - phi_s) * dyt2r(i,j-1)
               
               ! minimum and maximum among the nine local cells
               pmn = min (phi_nw, phi_n,  phi_ne, phi_w, phi(i,j),
     &                    phi_e,  phi_sw, phi_s,  phi_se)
               pmx = max (phi_nw, phi_n,  phi_ne, phi_w, phi(i,j),
     &                    phi_e,  phi_sw, phi_s,  phi_se)
            
               pmn = pmn - phi(i,j)
               pmx = pmx - phi(i,j)

               ! minimum and maximum deviation of phi within the cell
               w1  =  (p5*HTN(i,j) - cnx(i,j)) * gx(i,j) 
     &              + (p5*HTE(i,j) - cny(i,j)) * gy(i,j)
               w2  =  (p5*HTS(i,j) - cnx(i,j)) * gx(i,j) 
     &              - (p5*HTE(i,j) + cny(i,j)) * gy(i,j)
               w3  = -(p5*HTS(i,j) + cnx(i,j)) * gx(i,j) 
     &              - (p5*HTW(i,j) + cny(i,j)) * gy(i,j)
               w4  =  (p5*HTW(i,j) - cny(i,j)) * gy(i,j)
     &              - (p5*HTN(i,j) + cnx(i,j)) * gx(i,j)  

               qmn = min (w1, w2, w3, w4)
               qmx = max (w1, w2, w3, w4)
            
               ! the limiting coefficient
               if (abs(qmn) .ge. puny) then
                  w1 = max(c0, pmn/qmn)
               else
                  w1 = c1
               endif

               if (abs(qmx) .ge. puny) then
                  w2 = max(c0, pmx/qmx)
               else
                  w2 = c1
               endif
            
               w1 = min(c1, w1, w2)
 
               ! Limit the gradient components
               gx(i,j) = w1*gx(i,j)
               gy(i,j) = w1*gy(i,j)
               
            else                ! iflag(i,j) = 0.

               gx(i,j) = c0
               gy(i,j) = c0

            endif               ! iflag

         enddo                  ! j
      enddo                     ! i   

      ! update boundaries
      call bound (gx)
      call bound (gy)

      end subroutine limited_gradient

c=======================================================================

      subroutine conserved_sums (aicegsum, vicegsum, vsnogsum, 
     &     eicegsum, esnogsum)

!---!-------------------------------------------------------------------
!---! global sums of conserved variables over the physical grid
!---!-------------------------------------------------------------------

      use ice_state
      use ice_itd
      use ice_mpi_internal

      real (kind=dbl_kind), intent(out) ::
     &     aicegsum        ! global sum of ice area
     &,    vicegsum        ! global sum of ice volume
     &,    vsnogsum        ! global sum of snow volume
     &,    eicegsum        ! global sum of ice energy
     &,    esnogsum        ! global sum of snow energy

      real (kind=dbl_kind), dimension(ilo:ihi,jlo:jhi) :: 
     &     work1, work2, work3, work4, work5

      integer (kind=int_kind) :: i,j,n

      call aggregate

      do j = jlo, jhi
         do i = ilo, ihi
            work1(i,j) = tarea(i,j)*aice(i,j)
            work2(i,j) = tarea(i,j)*vice(i,j)
            work3(i,j) = tarea(i,j)*vsno(i,j)
            work4(i,j) = tarea(i,j)*eice(i,j)
            work5(i,j) = tarea(i,j)*esno(i,j)
         enddo
      enddo

      n = (ihi-ilo+1)*(jhi-jlo+1)
      aicegsum = ice_global_real_sum(n,work1)
      vicegsum = ice_global_real_sum(n,work2)
      vsnogsum = ice_global_real_sum(n,work3)
      eicegsum = ice_global_real_sum(n,work4)
      esnogsum = ice_global_real_sum(n,work5)
      
      end subroutine conserved_sums

c=======================================================================

      subroutine global_conservation (ai1, ai2, vi1, vi2, vs1, vs2,
     &     ei1, ei2, es1, es2)

!---!-------------------------------------------------------------------
!---! Check whether values of conserved quantities have changed
!---!-------------------------------------------------------------------

      use ice_calendar

      real (kind=dbl_kind), intent(in) ::
     &     ai1, ai2     ! initial and final global ice area
     &,    vi1, vi2     ! initial and final global ice volume
     &,    vs1, vs2     ! initial and final global snow volume
     &,    ei1, ei2     ! initial and final global ice energy
     &,    es1, es2     ! initial and final global snow energy

      if (ai1.gt.puny 
     &     .and. abs(ai1-ai2)/ai1 .gt. puny) then
         write (nu_diag,*)
         write (nu_diag,*) 'Ice area conserv error, istep =',istep1
         write (nu_diag,*) 'Initial global area =', ai1
         write (nu_diag,*) 'Final global area =', ai2
         write (nu_diag,*) 'Fractional error =', abs((ai1-ai2)/ai1)
         write (nu_diag,*) 'aice2-aice1 =', ai2 - ai1
         stop
      endif

      if (vi1.gt.puny 
     &     .and. abs(vi1-vi2)/vi1 .gt. puny) then
         write (nu_diag,*)
         write (nu_diag,*) 'Ice volume conserv error, istep =', istep1
         write (nu_diag,*) 'Initial global volume =', vi1
         write (nu_diag,*) 'Final global volume =', vi2
         write (nu_diag,*) 'Fractional error =', abs((vi1-vi2)/vi1)
         write (nu_diag,*) 'vice2-vice1 =', vi2 - vi1
         stop
      endif
         
      if (vs1.gt.puny 
     &     .and. abs(vs1-vs2)/vs1 .gt. puny) then
         write (nu_diag,*)
         write (nu_diag,*) 'Snow volume conserv error, istep =', istep1
         write (nu_diag,*) 'Initial global volume =', vs1
         write (nu_diag,*) 'Final global volume =', vs2
         write (nu_diag,*) 'Fractional error =', abs((vs1-vs2)/vs1)
         write (nu_diag,*) 'vsno2-vsno1 =', vs2 - vs1
         stop
      endif

      if (abs(ei1).gt.puny 
     &     .and. abs((ei1-ei2)/ei1) .gt. puny) then
         write (nu_diag,*)
         write (nu_diag,*) 'Ice energy conserv error, istep =', istep1
         write (nu_diag,*) 'Initial global energy =', ei1
         write (nu_diag,*) 'Final global energy =', ei2
         write (nu_diag,*) 'Fractional error =', abs((ei1-ei2)/ei1)
         write (nu_diag,*) 'eice2-eice1 =', ei2 - ei1
         stop
      endif
         
      if (abs(es1).gt.puny 
     &     .and. abs((es1-es2)/es1) .gt. puny) then
         write (nu_diag,*)
         write (nu_diag,*) 'Snow energy conserv error, istep =', istep1
         write (nu_diag,*) 'Initial global energy =', es1
         write (nu_diag,*) 'Final global energy =', es2
         write (nu_diag,*) 'Fractional error =', abs((es1-es2)/es1)
         write (nu_diag,*) 'esno2-esno1 =', es2 - es1
         stop
      endif

      end subroutine global_conservation

c=======================================================================

      function triarea (x, y)

!---!-------------------------------------------------------------------
!---! the area of a triangle given coordinates of its vertices
!---!-------------------------------------------------------------------

      real (kind=dbl_kind) :: triarea

      real (kind=dbl_kind), intent(in), dimension (3) ::
     &   x, y     ! x and y coordinates of vertices

      triarea = p5 * abs( (x(2)-x(1)) * (y(3)-y(1)) -
     &                    (y(2)-y(1)) * (x(3)-x(1)) )

      end function triarea

c=======================================================================

      end module ice_transport

c=======================================================================




