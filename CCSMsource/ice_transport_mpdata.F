! $Id: ice_transport_mpdata.F,v 1.13 2006/06/08 22:16:08 eclare Exp $
!=======================================================================
!BOP
!
! !MODULE: ice_transport_mpdata - horizontal advection (via mpdata)
!
! !DESCRIPTION: 
!
! Calculates horizontal advection using mpdata (Multidimensional
! Positive Definite Advection Transport Algorithm).
!
! !REVISION HISTORY:
!
! author Elizabeth C. Hunke
!
! !INTERFACE:

      module ice_transport_mpdata
!
! !USES:
!
      use ice_model_size
      use ice_domain
      use ice_constants
      use ice_grid
      use ice_fileunits
      use ice_init, only: advection
      use ice_work, only:  work_l1, work_l2
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
      subroutine transport_mpdata
!
! !DESCRIPTION:
!
! Computes the transport equations for one timestep using mpdata. Sets
! several fields into a work array and passes it to mpdata routine.
!
! !REVISION HISTORY:
!
! author Elizabeth C. Hunke
!
! !USES:
!      
      use ice_flux
      use ice_timers
      use ice_state
      use ice_itd   ! check_state
!
! !INPUT/OUTPUT PARAMETERS:
!
!EOP
!
      integer (kind=int_kind), parameter ::
     &   narr = 1 + 5*ncat     ! number of state variable arrays

      integer (kind=int_kind) :: 
     &   i, j, k, n       ! standard indices
     &,  narrays          ! counter for number of state variable arrays

      real (kind=dbl_kind) :: works(imt_local,jmt_local,narr)
      real (kind=dbl_kind) :: worke(imt_local,jmt_local,ntilay)

      call ice_timer_start(3)  ! advection

      !-----------------------------------------------------------------
      ! Make sure state variables are not too close to zero.
      !-----------------------------------------------------------------
      call check_state

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
       enddo
      enddo

      narrays = 1
      do n=1,ncat
         do j=jlo,jhi
          do i=ilo,ihi
             works(i,j,narrays+1) = aicen(i,j,n)
             works(i,j,narrays+2) = vicen(i,j,n)
             works(i,j,narrays+3) = vsnon(i,j,n)
             works(i,j,narrays+4) = -aicen(i,j,n)*Tsfcn(i,j,n)
             works(i,j,narrays+5) = -(esnon(i,j,n)
     &                          +rhos*Lfresh*vsnon(i,j,n)) ! for Tsnow<0
          enddo
         enddo

         ! if there were 3 arrays in loop, use 3 instead and change narr 
         narrays = narrays + 5  
      enddo                     

      do k=1,ntilay
         do j=jlo,jhi
          do i=ilo,ihi
             worke(i,j,k) = -eicen(i,j,k)
          enddo
         enddo
      enddo

      if (narr /= narrays) 
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
      ! retrieve advected fields from work array
      !-----------------------------------------------------------------

      do j=1,jmt_local
       do i=1,imt_local
          aice0(i,j) = works(i,j,1)
       enddo
      enddo

      narrays = 1               ! aice0 is first array
      do n=1,ncat
         do j=1,jmt_local
          do i=1,imt_local
             aicen(i,j,n) = works(i,j,narrays+1)
             vicen(i,j,n) = works(i,j,narrays+2)
             vsnon(i,j,n) = works(i,j,narrays+3)
             if (aicen(i,j,n) > puny) then
                Tsfcn(i,j,n) = -works(i,j,narrays+4)/aicen(i,j,n)
             else
                Tsfcn(i,j,n) = Tf(i,j)
             endif
             esnon(i,j,n) = -works(i,j,narrays+5)
     &                       -rhos*Lfresh*vsnon(i,j,n) ! for Tsnow<0 
          enddo
         enddo

         narrays = narrays + 5 
      enddo

      do k=1,ntilay
         do j=1,jmt_local
          do i=1,imt_local
             eicen(i,j,k) = -worke(i,j,k)
          enddo
         enddo
      enddo

      call ice_timer_stop(3) ! advection

      !-----------------------------------------------------------------
      ! mask
      !-----------------------------------------------------------------
      do n=1,ncat
       do j=1,jmt_local
        do i=1,imt_local
         if (.not.tmask(i,j)) then
          aicen(i,j,n)    = c0
          vicen(i,j,n)    = c0
          vsnon(i,j,n)    = c0
          Tsfcn(i,j,n)    = c0
          esnon(i,j,n)    = c0
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

      end subroutine transport_mpdata

!=======================================================================
!BOP
!
! !IROUTINE: mpdata(narrays,phi) - advection according to mpdata
!
! !INTERFACE:
!
      subroutine mpdata(narrays,phi)
!
! !DESCRIPTION:
!
! Smolarkiewicz, P. K., 1984:  A fully multidimensional positive 
! definite advection transport algorithm with small implicit 
! diffusion, J. Comput. Phys., 54, 325-362.
!
! !REVISION HISTORY:
!
! author Elizabeth C. Hunke
! Spring 2003: upwind algorithm modified by William Lipscomb (LANL)
!              to improve performance.
!
! !USES:
!
      use ice_calendar
      use ice_state, only: uvel, vvel
      use ice_exit
!
! !INPUT/OUTPUT PARAMETERS:
!
      integer (kind=int_kind), intent(in) :: narrays

      real (kind=dbl_kind), intent(inout) ::
     &   phi(imt_local,jmt_local,narrays)
!
!EOP
!
      integer (kind=int_kind) :: i, j, k, ix, iy, n
      integer (kind=int_kind), parameter :: 
     &   iord=3  ! 2nd order MPDATA paramter

      real (kind=dbl_kind) ::
     &   dive(ilo:ihi,jlo:jhi)
     &,  divn(ilo:ihi,jlo:jhi)
     &,  phiavg(imt_local,jmt_local)
     &,  uee(imt_local,jmt_local,narrays)
     &,  vnn(imt_local,jmt_local,narrays)
     &,  eps
     &,  upwind, y1, y2, a, h   ! function

      upwind(y1,y2,a,h)=p5*dyn_dt*h*((a+abs(a))*y1+(a-abs(a))*y2)

      eps = eps15

      if (advection=='upwind') then

       do j=jlo,jhi
        do i=ilo,ihi
         uee(i,j,1)=p5*(uvel(i,j)+uvel(i,j-1))
         vnn(i,j,1)=p5*(vvel(i,j)+vvel(i-1,j))
        enddo
       enddo

       call bound(uee(:,:,1))
       call bound(vnn(:,:,1))

       do n = 1,narrays
        do j=1,jhi
         do i=1,ihi
            work_l1(i,j)=
     &         upwind(phi(i,j,n),phi(i+1,j,n),uee(i,j,1),HTE(i,j))
            work_l2(i,j)=
     &         upwind(phi(i,j,n),phi(i,j+1,n),vnn(i,j,1),HTN(i,j))
         enddo
        enddo
        do j=jlo,jhi
         do i=ilo,ihi
           phi(i,j,n)=phi(i,j,n)
     &      -(work_l1(i,j)-work_l1(i-1,j)+work_l2(i,j)-work_l2(i,j-1))
     &                /tarea(i,j)
         enddo
        enddo
       enddo                     ! narrays

      else   ! mpdata with antidiffusive corrections

      do n = 1,narrays
       do j=jlo,jhi
        do i=ilo,ihi
         uee(i,j,n)=p5*(uvel(i,j)+uvel(i,j-1))
         vnn(i,j,n)=p5*(vvel(i,j)+vvel(i-1,j))
        enddo
       enddo
      enddo ! narrays

      call bound_narr(narrays,uee)
      call bound_narr(narrays,vnn)

      ! upwind 
      do n = 1,narrays
       do j=1,jhi
        do i=1,ihi
           work_l1(i,j)=
     &          upwind(phi(i,j,n),phi(i+1,j,n),uee(i,j,n),HTE(i,j))
           work_l2(i,j)=
     &          upwind(phi(i,j,n),phi(i,j+1,n),vnn(i,j,n),HTN(i,j))
        enddo
       enddo
       do j=jlo,jhi
        do i=ilo,ihi
           phi(i,j,n) = phi(i,j,n)
     &       -(work_l1(i,j)-work_l1(i-1,j)+work_l2(i,j)-work_l2(i,j-1))
     &                          /tarea(i,j)
        enddo
       enddo
      enddo ! narrays

      call bound_narr(narrays,phi)

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
     1   /(phi(i+1,j,n)+phi(i,j,n)+eps) - dyn_dt*uee(i,j,n)*dive(i,j)

         vnn(i,j,n) = abs(vnn(i,j,n))*(phi(i,j+1,n)-phi(i,j,n))
     1   /(phi(i,j+1,n)+phi(i,j,n)+eps) - dyn_dt*vnn(i,j,n)*divn(i,j)
        enddo
       enddo
      enddo ! narrays

      call bound_narr(narrays,uee)
      call bound_narr(narrays,vnn)

      ! upwind with antidiffusive velocities
      do n = 1,narrays

       do j=1,jhi
        do i=1,ihi
         work_l1(i,j) = upwind(phi(i,j,n),phi(i+1,j,n),
     &                         uee(i,j,n),HTE(i,j))
         work_l2(i,j) = upwind(phi(i,j,n),phi(i,j+1,n),
     &                         vnn(i,j,n),HTN(i,j))
        enddo
       enddo

       ix=-1
       do j=jlo,jhi
        do i=ilo,ihi
         phi(i,j,n)=phi(i,j,n)
     &     -(work_l1(i,j)-work_l1(i-1,j)+work_l2(i,j)-work_l2(i,j-1))
     &                          /tarea(i,j)

         if (phi(i,j,n) < -eps12) then 
          ix=i
          iy=j
         elseif (phi(i,j,n) < c0) then
          phi(i,j,n) = c0 
         endif
        enddo
       enddo

       if (ix >= 0) then
        write (nu_diag,*)  my_task,ix,iy,' transport unstable ',istep,k
        write (nu_diag,*)  'mpdata phi = ',phi(ix,iy,n),' n = ',n
        call abort_ice ('(mpdata) transport unstable')
       endif

      enddo ! narrays

      if (k < iord) call bound_narr(narrays,phi)

      enddo ! iord

      endif ! upwind or mpdata

      end subroutine mpdata

!=======================================================================

      end module ice_transport_mpdata

!=======================================================================
