c $Id: $
c=======================================================================
!---! horizontal advection (via mpdata)
!---!
!---! author Elizabeth C. Hunke
!---!        Fluid Dynamics Group, Los Alamos National Laboratory
!---!
!---! Copyright, 2000.  The Regents of the University of California.
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

      implicit none

      character (len=char_len) :: 
     &   advection   ! type of advection algorithm used
                     ! upwind => 1st order mpdata scheme (donor cell)
                     ! anything else => 2nd order mpdata scheme

c=======================================================================

      contains

c=======================================================================

      subroutine transport

!---!-------------------------------------------------------------------
!---! computes the transport equations for one timestep
!---!-------------------------------------------------------------------

      use ice_flux
      use ice_timers
      use ice_state

      integer (kind=int_kind), parameter :: narr = 1 + 5*ncat ! number
                                            ! of state variable arrays
      integer (kind=int_kind) :: 
     &   i, j, layer, nc  ! standard indices
     &,  narrays          ! counter for number of state variable arrays

      real (kind=dbl_kind) :: works(imt_local,jmt_local,narr)
      real (kind=dbl_kind) :: worke(imt_local,jmt_local,ntilay)

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
     & write(6,*) "Wrong number of arrays in transport bound call"

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
        narrays = narrays + c5    ! if there were 3 arrays in this 
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

      call ice_timer_stop(3) ! advection

      !-----------------------------------------------------------------
      ! mask
      !-----------------------------------------------------------------
c$OMP PARALLEL DO PRIVATE(i,j,nc)
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

      end subroutine transport

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
c$OMP PARALLEL DO PRIVATE(i,j)
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
c$OMP PARALLEL DO PRIVATE(i,j)
       do j=1,jhi
        do i=1,ihi
         fx(i,j)=upwind(phi(i,j,n),phi(i+1,j,n),uee(i,j,n),HTE(i,j))
         fy(i,j)=upwind(phi(i,j,n),phi(i,j+1,n),vnn(i,j,n),HTN(i,j))
        enddo
       enddo
c$OMP PARALLEL DO PRIVATE(i,j)
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
c$OMP PARALLEL DO PRIVATE(i,j)
       do j=1,jhi
        do i=1,ihi
         phiavg(i,j) = p25*(phi(i,j,n) + phi(i+1,j,n)
     &               + phi(i+1,j+1,n) + phi(i,j+1,n))
        enddo
       enddo

c$OMP PARALLEL DO PRIVATE(i,j,dive,divn)
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
c$OMP PARALLEL DO PRIVATE(i,j)
       do j=1,jhi
        do i=1,ihi
         fx(i,j)=upwind(phi(i,j,n),phi(i+1,j,n),uee(i,j,n),HTE(i,j))
         fy(i,j)=upwind(phi(i,j,n),phi(i,j+1,n),vnn(i,j,n),HTN(i,j))
        enddo
       enddo

       ix=-1
c$OMP PARALLEL DO PRIVATE(i,j)
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
        write (6,*)  my_task,ix,iy,' transport unstable ',istep,k
        write (6,*)  'mpdata phi = ',phi(ix,iy,n),' n = ',n
        stop
       endif
      enddo ! narrays

      call bound_narr(narrays,phi)

      enddo ! iord
 555  continue

      end subroutine mpdata

c=======================================================================

      end module ice_transport

c=======================================================================
