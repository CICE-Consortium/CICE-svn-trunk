c $Id: ice_itd.F,v 1.15 2004/02/25 17:35:13 eclare Exp $
!=======================================================================
!BOP
!
! !MODULE: ice_itd - initialize and redistribute ice in the ITD
!
! !DESCRIPTION:
!
! Routines to initialize the ice thickness distribution and 
! utilities to redistribute ice among categories. These routines 
! are not specific to a particular numerical implementation. \\
!
! See Bitz, C.M., and W.H. Lipscomb, 1999: 
! An energy-conserving thermodynamic model of sea ice,
! J. Geophys. Res., 104, 15,669--15,677. \\
!     
! See Bitz, C.M., M.M. Holland, A.J. Weaver, M. Eby, 2001: 
! Simulating the ice-thickness distribution in a climate model,
! J. Geophys. Res., 106, 2441--2464. \\
!
! !REVISION HISTORY:
!
! author: C. M. Bitz, UW
!         Elizabeth C. Hunke, LANL
!         William H. Lipscomb, LANL
!
! Summer 2003: Vectorized by Clifford Chen (Fujitsu) and William Lipscomb (LANL)
!
! !INTERFACE:
!
      module ice_itd
!
! !USES:
!
      use ice_kinds_mod
      use ice_model_size
      use ice_constants
      use ice_state
      use ice_fileunits
!
!EOP
!
      implicit none

      integer (kind=int_kind) ::
     &   kitd               ! type of itd conversions 
                            !   0 = delta function
                            !   1 = linear remap
     &,  ilyr1 (ncat)       ! position of the top layer in each cat
     &,  ilyrn (ncat)       ! position of the bottom layer in each cat

      real (kind=dbl_kind), parameter ::
     &   hi_min = p01       ! minimum ice thickness allowed (m)

      real (kind=dbl_kind) ::
     &   hin_max(0:ncat)    ! category limits                (m)

!-------------------------------------------------------------------
! a note regarding hi_min and hin_max(0):
! both represent a minimum ice thickness.  hin_max(0) is
! intended to be used for particular numerical implementations
! of category conversions in the ice thickness distribution.
! hi_min is a more general purpose parameter, but is specifically 
! for maintaining stability in the thermodynamics.  Currently,
! hi_min = 0.1 m
! hin_max(0) = 0.1 m for the delta function itd
! hin_max(0) = 0.0 m for linear remapping
!
! similarly, there are two values of minimum snow thickness
! (the other is defined in ice_vthermo.H since it is used only
! for thermo.) 
!
! Also note that the upper limit on the thickest category
! is only used for the linear remapping scheme
! and it is not a true upper limit on the thickness
!-------------------------------------------------------------------

!=======================================================================

      contains

!=======================================================================
!BOP
!
! !IROUTINE: init_itd - initalize area fraction and thickness boundaries for ITD
!
! !INTERFACE:
!
      subroutine init_itd
!
! !DESCRIPTION:
!
! Initialize area fraction and thickness boundaries for the itd model
!
! !REVISION HISTORY:
!
! authors: William H. Lipscomb, LANL
!          Elizabeth C. Hunke LANL
!          C. M. Bitz UW
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
     &     n    ! thickness category index

      real (kind=dbl_kind):: cc1, cc2, cc3, x1

      if (ncat.eq.1 .and. kitd.eq.1) then
         write (nu_diag,*) 'Remapping the ITD is not allowed for ncat=1'
         write (nu_diag,*) 'Use the delta function ITD option instead'
         call abort_ice ('(init_itd) Linear remapping not allowed')
      endif

      if (kitd.eq.1) then
         ! linear remapping itd category limits
         cc1 = c3/real(ncat)
         cc2 = c15*cc1      
         cc3 = c3

         hin_max(0) = c0        ! minimum ice thickness, m
      else
         ! delta function itd category limits
         cc1 = max(1.1_dbl_kind/ncat,c1*hi_min)
         cc2 = c25*cc1
         cc3 = 2.25_dbl_kind

         ! hin_max(0) should not be zero 
         ! use some caution in making it less than 0.10
         hin_max(0) = hi_min    ! minimum ice thickness, m
      endif                     ! kitd

      do n = 1, ncat
         x1 = real(n-1) / real(ncat)
         hin_max(n) = hin_max(n-1) + cc1 + cc2*(c1 + tanh(cc3*(x1-c1)))
      enddo

      if (my_task.eq.master_task) then
         write (nu_diag,*) ''
         write (nu_diag,*) 'hin_max(n-1) < Cat n < hin_max(n)'
         do n = 1, ncat
            write (nu_diag,*) hin_max(n-1),' < Cat ',n, ' < ',hin_max(n)
         enddo
      endif

      !-----------------------------------------------------------------
      ! vectors identifying first and last layer in each category
      !-----------------------------------------------------------------
      ilyr1(1) = 1                       ! if nilyr  = 4
      ilyrn(1) = nilyr                   !   ilyr1 = { 1,5,9 }
      do n = 2,ncat                      !   ilyrn = { 4,8,12} etc
         ilyr1(n) = ilyrn(n-1) + 1
         ilyrn(n) = ilyrn(n-1) + nilyr
      enddo

      end subroutine init_itd

!=======================================================================
!BOP
!
! !IROUTINE: aggregate - aggregate ice state variables
!
! !INTERFACE:
!
      subroutine aggregate
!
! !DESCRIPTION:
!
! Aggregate ice state variables over thickness categories.
!
! !REVISION HISTORY:
!
! authors: C. M. Bitz, UW
!          W. H. Lipscomb, LANL
!
! !USES:
!
      use ice_domain
      use ice_flux, only : Tf 
      use ice_grid
!
! !INPUT/OUTPUT PARAMETERS:
!
!EOP
!
      integer (kind=int_kind) :: i, j, k, n 

      integer (kind=int_kind), dimension (1:(ihi-ilo+1)*(jhi-jlo+1)) ::
     &  indxi                 ! compressed indices in i/j directions
     &, indxj

      integer (kind=int_kind) ::
     &  icells                ! number of ocean/ice cells
     &, ij                    ! combined i/j horizontal index

      !-----------------------------------------------------------------
      ! Initialize
      !-----------------------------------------------------------------
      do j=jlo,jhi
      do i=ilo,ihi
         aice0(i,j) = c1
         aice(i,j) = c0
         vice(i,j) = c0
         vsno(i,j) = c0
         eice(i,j) = c0
         esno(i,j) = c0
         Tsfc(i,j) = c0
      enddo
      enddo

      !-----------------------------------------------------------------
      ! Aggregate
      !-----------------------------------------------------------------

      icells = 0
      do j = jlo, jhi
      do i = ilo, ihi
        if (tmask(i,j)) then
          icells = icells + 1
          indxi(icells) = i
          indxj(icells) = j
        endif   ! tmask
      enddo
      enddo

      do n = 1, ncat

cdir$ ivdep      !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
         do ij = 1, icells
            i = indxi(ij)
            j = indxj(ij)
            aice(i,j) = aice(i,j) + aicen(i,j,n)
            vice(i,j) = vice(i,j) + vicen(i,j,n)
            vsno(i,j) = vsno(i,j) + vsnon(i,j,n)
            esno(i,j) = esno(i,j) + esnon(i,j,n)
            Tsfc(i,j) = Tsfc(i,j) + Tsfcn(i,j,n)*aicen(i,j,n)
         enddo                  ! ij

         do k = 1, nilyr
cdir$ ivdep      !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
            do ij = 1, icells
               i = indxi(ij)
               j = indxj(ij)
               eice(i,j) = eice(i,j) + eicen(i,j,(n-1)*nilyr+k)
            enddo
         enddo                  ! k

      enddo                     ! n

      !-----------------------------------------------------------------
      ! Temperature, open water fraction
      !-----------------------------------------------------------------
      do j=jlo,jhi
      do i=ilo,ihi
         if (aice(i,j) > c0)  then
            aice0(i,j) = max (c1 - aice(i,j), c0)
            Tsfc(i,j)  = Tsfc(i,j) / aice(i,j)
         else
            Tsfc(i,j)  = Tf(i,j)
         endif
      enddo                     ! i
      enddo                     ! j

      end subroutine aggregate

!=======================================================================
!BOP
!
! !IROUTINE: aggregate_area - aggregate ice area
!
! !INTERFACE:
!
      subroutine aggregate_area
!
! !DESCRIPTION:
!
! Aggregate ice area (but not other state variables) over thickness categories
!
! !REVISION HISTORY:
!
! authors: William H. Lipscomb, LANL
!          modified Jan 2004 by Clifford Chen, Fujitsu
!
! !USES:
!
      use ice_exit
!
! !INPUT/OUTPUT PARAMETERS:
!
!EOP
!
      integer (kind=int_kind) :: i, j, n

      logical (kind=log_kind) ::
     &   outbound          ! = .true. if aggregate ice area out of bounds

      !-----------------------------------------------------------------
      ! Initialize
      !-----------------------------------------------------------------
      do j=jlo,jhi
      do i=ilo,ihi
         aice(i,j) = c0
      enddo
      enddo

      !-----------------------------------------------------------------
      ! Aggregate
      !-----------------------------------------------------------------
      do n = 1, ncat
         do j=jlo,jhi
         do i=ilo,ihi
            aice(i,j) = aice(i,j) + aicen(i,j,n)
         enddo                  ! i
         enddo                  ! j
      enddo                     ! n

      outbound = .false.

      do j = jlo, jhi
      do i = ilo, ihi

      !-----------------------------------------------------------------
      ! Bug check
      !-----------------------------------------------------------------
         if (aice(i,j) > c1+puny .or.
     &       aice(i,j) < -puny) then
            outbound = .true.
         endif

      !-----------------------------------------------------------------
      ! open water fraction
      !-----------------------------------------------------------------
         aice0(i,j) = max (c1 - aice(i,j), c0)

      enddo                     ! i
      enddo                     ! j

      if ( outbound ) then      ! area out of bounds
        do j = jlo, jhi
        do i = ilo, ihi

         if (aice(i,j) > c1+puny .or.
     &       aice(i,j) < -puny) then
            write(nu_diag,*) ''
            write(nu_diag,*) 'aggregate ice area out of bounds'
            write(nu_diag,*) 'my_task, i, j, aice:', 
     &                        my_task, i, j, aice(i,j)
!!            write(nu_diag,*) 'aicen =', aicen(i,j,:) 
            call abort_ice('ice_itd: aggregate_area')
         endif

        enddo                   ! i
        enddo                   ! j
      endif                     ! outbound

      end subroutine aggregate_area

!=======================================================================
!BOP
!
! !IROUTINE: bound_aggregate - bound calls for aggregate ice state
!
! !INTERFACE:
!
      subroutine bound_aggregate
!
! !DESCRIPTION:
!
! Get ghost cell values for aggregate ice state variables
!
! !REVISION HISTORY:
!
! authors: William H. Lipscomb, LANL
!
! !USES:
!
! !INPUT/OUTPUT PARAMETERS:
!
!EOP
!
      use ice_grid

      call bound (aice0)
      call bound (aice)
      call bound (vice)
      call bound (vsno)
      call bound (eice)
      call bound (esno)
      call bound (Tsfc)

      end subroutine bound_aggregate

!=======================================================================
!BOP
!
! !IROUTINE: bound_state - bound calls for ice state variables
!
! !INTERFACE:
!
      subroutine bound_state
!
! !DESCRIPTION:
!
! Get ghost cell values for ice state variables in each thickness category
!
! !REVISION HISTORY:
!
! authors: William H. Lipscomb, LANL
!
! !USES:
!
! !INPUT/OUTPUT PARAMETERS:
!
!EOP
!
      use ice_grid

      call bound_narr (ncat,aicen)
      call bound_narr (ncat,vicen)
      call bound_narr (ncat,vsnon)
      call bound_narr (ntilay,eicen)
      call bound_narr (ncat,esnon)
      call bound_narr (ncat,Tsfcn)

      end subroutine bound_state

!=======================================================================
!BOP
!
! !IROUTINE: check_state - require certain fields to be monotone
!
! !INTERFACE:
!
      subroutine check_state
!
! !DESCRIPTION:
!
!  Insist that certain fields are monotone.
!  Should not be necessary if all is well, 
!  but best to keep going. Model will not conserve
!  energy and water if fields are zeroed here.
!
! !REVISION HISTORY:
!
! author: C. M. Bitz, UW
!
! !USES:
!
      use ice_flux
!
! !INPUT/OUTPUT PARAMETERS:
!
!EOP
!
      integer (kind=int_kind) ::
     &   i, j            ! horizontal indices
     &,  n               ! thickness category index
     &,  k               ! ice layer index

      integer (kind=int_kind), dimension (ilo:ihi,jlo:jhi) ::
     &   zerout    ! 0=false, 1=true

      do n = 1, ncat

         do j = jlo, jhi
         do i = ilo, ihi
            if (aicen(i,j,n) < puny .or. vicen(i,j,n) < puny) then
               zerout(i,j) = 1
            else
               zerout(i,j) = 0
            endif
         enddo
         enddo

         do k = 1, nilyr
            do j = jlo, jhi
            do i = ilo, ihi
               if (eicen(i,j,ilyr1(n)+k-1) > -puny) zerout(i,j) = 1
            enddo               ! i
            enddo               ! j
         enddo                  ! k

         do j = jlo, jhi
         do i = ilo, ihi

            if (zerout(i,j)==1) then
               aice0(i,j) = aice0(i,j) + aicen(i,j,n) 
               aicen(i,j,n) = c0
               vicen(i,j,n) = c0
               vsnon(i,j,n) = c0
               esnon(i,j,n) = c0
               Tsfcn(i,j,n) = Tf(i,j)
            elseif (vsnon(i,j,n) <= puny) then
               vsnon(i,j,n) = c0
               esnon(i,j,n) = c0
            endif

            if (vsnon(i,j,n) > puny) then
               if (-esnon(i,j,n)/vsnon(i,j,n)-Lfresh*rhos < eps04) 
     &              esnon(i,j,n) = -vsnon(i,j,n)*(Lfresh*rhos + eps04)
            endif

         enddo                  ! i
         enddo                  ! j

         do k = 1, nilyr
            do j = jlo, jhi
            do i = ilo, ihi
               if (zerout(i,j)==1) eicen(i,j,ilyr1(n)+k-1) = c0
            enddo
            enddo
         enddo                  ! k

      enddo                     ! n

      do j=jlo,jhi
      do i=ilo,ihi
         if (aice0(i,j) < puny) aice0(i,j) = c0
      enddo
      enddo

      end subroutine check_state

!=======================================================================
!BOP
!
! !IROUTINE: rebin - rebins thicknesses into defined categories
!
! !INTERFACE:
!
      subroutine rebin
!
! !DESCRIPTION:
!
! Rebins thicknesses into defined categories
!
! !REVISION HISTORY:
!
! authors: William H. Lipscomb, LANL
!          Elizabeth C. Hunke LANL 
!
! !USES:
!
      use ice_grid
!
! !INPUT/OUTPUT PARAMETERS:
!
!EOP
!
      integer (kind=int_kind) ::
     &   i,j                ! horizontal indices
     &,  n                  ! category index

      logical (kind=log_kind) ::
     &   shiftflag          ! = .true. if ice must be shifted

      integer (kind=int_kind), dimension(ilo:ihi,jlo:jhi,ncat-1) ::
     &   donor              ! donor category index

      real (kind=dbl_kind), dimension(ilo:ihi,jlo:jhi,ncat) ::
     &   hicen              ! ice thickness for each cat        (m)

      real (kind=dbl_kind), dimension(ilo:ihi,jlo:jhi,ncat-1) ::
     &   daice              ! ice area transferred
     &,  dvice              ! ice volume transferred

      !-----------------------------------------------------------------
      ! Compute ice thickness.
      !-----------------------------------------------------------------
      do n = 1, ncat
         do j = jlo, jhi
         do i = ilo, ihi
            if (aicen(i,j,n) > puny) then
               hicen(i,j,n) = vicen(i,j,n) / aicen(i,j,n)
            else
               hicen(i,j,n) = c0
            endif
         enddo                  ! i
         enddo                  ! j
      enddo                     ! n

      !-----------------------------------------------------------------
      ! make sure thickness of cat 1 is at least hin_max(0)
      !-----------------------------------------------------------------
      do j = jlo, jhi
      do i = ilo, ihi
         if (aicen(i,j,1) > puny) then
            if (hicen(i,j,1) <= hin_max(0) .and. hin_max(0) > c0 ) then
               aicen(i,j,1) = vicen(i,j,1) / hin_max(0)
               hicen(i,j,1) = hin_max(0)
            endif
         endif
      enddo                     ! i
      enddo                     ! j

      !-----------------------------------------------------------------
      ! If a category thickness is not in bounds, shift the
      ! entire area, volume, and energy to the neighboring category
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------
      ! Initialize shift arrays
      !-----------------------------------------------------------------

      do n = 1, ncat-1
         do j = jlo, jhi
         do i = ilo, ihi
            donor(i,j,n) = 0
            daice(i,j,n) = c0
            dvice(i,j,n) = c0
         enddo
         enddo
      enddo

      !-----------------------------------------------------------------
      ! Move thin categories up
      !-----------------------------------------------------------------

      do n = 1, ncat-1          ! loop over category boundaries

      !-----------------------------------------------------------------
      ! identify thicknesses that are too big
      !-----------------------------------------------------------------
         shiftflag = .false.
         do j = jlo, jhi
         do i = ilo, ihi
            if (aicen(i,j,n) > puny .and. 
     &          hicen(i,j,n) > hin_max(n)) then
               shiftflag = .true.
               donor(i,j,n) = n
               daice(i,j,n) = aicen(i,j,n)
               dvice(i,j,n) = vicen(i,j,n)
            endif
         enddo                  ! i
         enddo                  ! j

         if (shiftflag) then

      !-----------------------------------------------------------------
      ! shift ice between categories
      !-----------------------------------------------------------------
            call shift_ice (donor, daice, dvice, hicen)
             
      !-----------------------------------------------------------------
      ! reset shift parameters
      !-----------------------------------------------------------------
            do j = jlo, jhi
            do i = ilo, ihi
               donor(i,j,n) = 0
               daice(i,j,n) = c0
               dvice(i,j,n) = c0
            enddo
            enddo

         endif                  ! shiftflag

      enddo                     ! n

      !-----------------------------------------------------------------
      ! Move thick categories down
      !-----------------------------------------------------------------

      do n = ncat-1, 1, -1      ! loop over category boundaries

      !-----------------------------------------------------------------
      ! identify thicknesses that are too small
      !-----------------------------------------------------------------
         shiftflag = .false.
         do j = jlo, jhi
         do i = ilo, ihi
            if (aicen(i,j,n+1) > puny .and. 
     &          hicen(i,j,n+1) <= hin_max(n)) then
               shiftflag = .true.
               donor(i,j,n) = n+1
               daice(i,j,n) = aicen(i,j,n+1)
               dvice(i,j,n) = vicen(i,j,n+1)
            endif
         enddo                  ! i
         enddo                  ! j

         if (shiftflag) then

      !-----------------------------------------------------------------
      ! shift ice between categories
      !-----------------------------------------------------------------
            call shift_ice (donor, daice, dvice, hicen)

      !-----------------------------------------------------------------
      ! reset shift parameters
      !-----------------------------------------------------------------
            do j = jlo, jhi
            do i = ilo, ihi
               donor(i,j,n) = 0
               daice(i,j,n) = c0
               dvice(i,j,n) = c0
            enddo
            enddo

         endif                  ! shiftflag

      enddo                     ! n


      end subroutine rebin

!=======================================================================
!BOP
!
! !IROUTINE: reduce_area - reduce area when ice melts for special case ncat=1
!
! !INTERFACE:
!
      subroutine reduce_area(hice1_old, hice1)
!
! !DESCRIPTION:
!
! Reduce area when ice melts for special case of ncat=1
!
! Use CSM 1.0-like method of reducing ice area
! when melting occurs: assume only half the ice volume
! change goes to thickness decrease, the other half
! to reduction in ice fraction
!
! !REVISION HISTORY:
!
! authors: C. M. Bitz, UW 
! modified by: Elizabeth C. Hunke, LANL
!
! !USES:
!
      use ice_grid
!
! !INPUT/OUTPUT PARAMETERS:
!
      real (kind=dbl_kind), dimension(ilo:ihi,jlo:jhi),
     &     intent(in) ::
     &   hice1_old   ! old ice thickness for category 1 (m)

      real (kind=dbl_kind), dimension(ilo:ihi,jlo:jhi),
     &     intent(inout) ::
     &   hice1       ! new ice thickness for category 1 (m)
!
!EOP
!
      integer (kind=int_kind) ::
     &   i, j        ! horizontal indices

      real (kind=dbl_kind) ::
     &   hi0         ! current hi for ice fraction adjustment
     &,  dai0        ! change in aice for ice fraction adjustment
     &,  dhi         ! hice1 - hice1_old

      do j=jlo,jhi
      do i=ilo,ihi
         if (tmask(i,j)) then

      !-----------------------------------------------------------------
      ! make sure thickness of cat 1 is at least hin_max(0)
      !-----------------------------------------------------------------

            if (hice1(i,j) <= hin_max(0) .and. hin_max(0) > c0 ) then
               aicen(i,j,1) = vicen(i,j,1) / hin_max(0)
               hice1(i,j) = hin_max(0)
            endif

            if (aicen(i,j,1) > c0) then
               dhi = hice1(i,j) - hice1_old(i,j)
               if (dhi < c0) then  
                  hi0  = vicen(i,j,1) / aicen(i,j,1)
                  dai0 = vicen(i,j,1) / (hi0-p5*dhi)
     &                 - aicen(i,j,1)
                  aicen(i,j,1) = aicen(i,j,1) + dai0 
               endif
            endif

         endif                  ! tmask
      enddo                     ! i
      enddo                     ! j

      end subroutine reduce_area

!=======================================================================
!BOP
!
! !IROUTINE: shift_ice - shift ice across category boundaries 
!
! !INTERFACE:
!
      subroutine shift_ice (donor, daice, dvice, hicen)
!
! !DESCRIPTION:
!
! Shift ice across category boundaries, conserving area, volume, and
! energy.
!
! !REVISION HISTORY:
!
! authors: William H. Lipscomb, LANL
!          Elizabeth C. Hunke, LANL
!
! !USES:
!
      use ice_flux
      use ice_work, only: worka
      use ice_exit
!
! !INPUT/OUTPUT PARAMETERS:
!
      integer (kind=int_kind), dimension(ilo:ihi,jlo:jhi,ncat-1),
     &   intent(in) ::
     &   donor             ! donor category index

      real (kind=dbl_kind), dimension(ilo:ihi,jlo:jhi,ncat-1),
     &     intent(inout) ::
     &   daice             ! ice area transferred across boundary
     &,  dvice             ! ice volume transferred across boundary

      real (kind=dbl_kind), dimension(ilo:ihi,jlo:jhi,ncat),
     &     intent(inout) ::
     &   hicen             ! ice thickness for each cat        (m)
!
!EOP
!
      integer (kind=int_kind) ::
     &   i, j              ! horizontal indices
     &,  n                 ! thickness category index
     &,  n2                ! receiver category
     &,  n1                ! donor category
     &,  k                 ! ice layer index

      real (kind=dbl_kind), dimension(ilo:ihi,jlo:jhi,ncat) ::
     &   aTsfn

      real (kind=dbl_kind) ::
     &   dvsnow            ! snow volume transferred
     &,  desnow            ! snow energy transferred
     &,  deice             ! ice energy transferred
     &,  daTsf             ! aicen*Tsfcn transferred

      integer (kind=int_kind), dimension (1:(ihi-ilo+1)*(jhi-jlo+1)) ::
     &  indxi           ! compressed indices for i/j directions
     &, indxj

      integer (kind=int_kind) ::
     &  icells          ! number of cells with ice to transfer
     &, ij              ! combined i/j horizontal index

      logical (kind=log_kind) ::
     &  daice_negative       ! true if daice < -puny
     &, dvice_negative       ! true if dvice < -puny
     &, daice_greater_aicen  ! true if daice > aicen
     &, dvice_greater_vicen  ! true if dvice > vicen

      !-----------------------------------------------------------------
      ! Define a variable equal to aicen*Tsfcn
      !-----------------------------------------------------------------
      do n = 1, ncat
         do j = jlo,jhi
         do i = ilo,ihi
            aTsfn(i,j,n) = aicen(i,j,n)*Tsfcn(i,j,n)
         enddo                  ! i
         enddo                  ! j
      enddo                     ! n

      !-----------------------------------------------------------------
      ! Check for daice or dvice out of range, allowing for roundoff error
      !-----------------------------------------------------------------
      ! Note: daice < 0 or dvice < 0 usually happens when category n 
      ! has a small area, with hice(n) very close to a boundary.  Then 
      ! the coefficients of g(h) are large, and the computed daice and 
      ! dvice can be in error. If this happens, it is best to transfer 
      ! either the entire category or nothing at all, depending on which 
      ! side of the boundary hice(n) lies.
      !-----------------------------------------------------------------

      do n = 1, ncat-1

         daice_negative = .false.
         dvice_negative = .false.
         daice_greater_aicen = .false.
         dvice_greater_vicen = .false.

         do j = jlo,jhi
         do i = ilo,ihi

            if (donor(i,j,n) > 0) then 
               n1 = donor(i,j,n)

               if (daice(i,j,n) < c0) then
                  if (daice(i,j,n) > -puny) then   
                     if ((n1.eq.n .and. hicen(i,j,n1) > hin_max(n))
     &                             .or.
     &                 (n1.eq.n+1 .and. hicen(i,j,n1)<=hin_max(n))) then
                        daice(i,j,n) = aicen(i,j,n1) ! shift entire category
                        dvice(i,j,n) = vicen(i,j,n1)
                     else             
                        daice(i,j,n) = c0 ! shift no ice
                        dvice(i,j,n) = c0
                     endif
                  else
                     daice_negative = .true.
                  endif
               endif
         
               if (dvice(i,j,n) < c0) then
                  if (dvice(i,j,n) > -puny) then   
                     if ((n1.eq.n .and. hicen(i,j,n1) > hin_max(n))
     &                             .or.
     &                 (n1.eq.n+1 .and. hicen(i,j,n1)<=hin_max(n))) then
                        daice(i,j,n) = aicen(i,j,n1) ! shift entire category
                        dvice(i,j,n) = vicen(i,j,n1)
                     else             
                        daice(i,j,n) = c0 ! shift no ice
                        dvice(i,j,n) = c0
                     endif
                  else
                     dvice_negative = .true.
                  endif
               endif

            ! If daice is close to aicen, set daice = aicen.
               if (daice(i,j,n) > aicen(i,j,n1)-puny) then
                  if (daice(i,j,n) < aicen(i,j,n1)+puny) then
                     daice(i,j,n) = aicen(i,j,n1)
                     dvice(i,j,n) = vicen(i,j,n1)
                  else
                     daice_greater_aicen = .true.
                  endif
               endif    

               if (dvice(i,j,n) > vicen(i,j,n1)-puny) then
                  if (dvice(i,j,n) < vicen(i,j,n1)+puny) then   
                     daice(i,j,n) = aicen(i,j,n1)
                     dvice(i,j,n) = vicen(i,j,n1)
                  else
                     dvice_greater_vicen = .true.
                  endif
               endif
               
            endif               ! donor > 0 
         enddo                  ! i
         enddo                  ! j

      !-----------------------------------------------------------------
      ! error messages
      !-----------------------------------------------------------------

         if (daice_negative) then
            do j = jlo,jhi
            do i = ilo,ihi
               if (donor(i,j,n) > 0 .and. daice(i,j,n) <= -puny) then
                  write(nu_diag,*) my_task,':',i,j,
     &                 'ITD Neg daice =',daice(i,j,n),' boundary',n
                  call abort_ice ('ice: ITD Neg daice')
               endif
            enddo
            enddo
         endif

         if (dvice_negative) then
            do j = jlo,jhi
            do i = ilo,ihi
               if (donor(i,j,n) > 0 .and. dvice(i,j,n) <= -puny) then
                  write(nu_diag,*) my_task,':',i,j,
     &                 'ITD Neg dvice =',dvice(i,j,n),' boundary',n
                  call abort_ice ('ice: ITD Neg dvice')
               endif
            enddo
            enddo
         endif

         if (daice_greater_aicen) then
            do j = jlo,jhi
            do i = ilo,ihi
               if (donor(i,j,n) > 0) then
                  n1 = donor(i,j,n)
                  if (daice(i,j,n) >= aicen(i,j,n1)+puny) then
                     write(nu_diag,*) my_task,':',i,j,
     &                    'ITD daice > aicen, cat',n1
                     write(nu_diag,*) my_task,':',i,j,
     &                    'daice =', daice(i,j,n),
     &                    'aicen =', aicen(i,j,n1)
                     call abort_ice ('ice: ITD daice > aicen')
                  endif
               endif
            enddo
            enddo
         endif

         if (dvice_greater_vicen) then
            do j = jlo,jhi
            do i = ilo,ihi
               if (donor(i,j,n) > 0) then
                  n1 = donor(i,j,n)
                  if (dvice(i,j,n) >= vicen(i,j,n1)+puny) then
                     write(nu_diag,*) my_task,':',i,j,
     &                    'ITD dvice > vicen, cat',n1
                     write(nu_diag,*) my_task,':',i,j,
     &                    'dvice =', dvice(i,j,n),
     &                    'vicen =', vicen(i,j,n1)
                     call abort_ice ('ice: ITD dvice > vicen')
                  endif
               endif
            enddo
            enddo
         endif

      enddo                     ! boundaries 1 to ncat-1            
         
      !-----------------------------------------------------------------
      ! transfer volume and energy between categories
      !-----------------------------------------------------------------

      do n = 1, ncat-1
         icells = 0
         do j = jlo, jhi
         do i = ilo, ihi
           if (daice(i,j,n) > c0) then ! daice(n) can be < puny
             icells = icells + 1
             indxi(icells) = i
             indxj(icells) = j
           endif   ! tmask
         enddo
         enddo

cdir$ ivdep      !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
         do ij = 1, icells
            i = indxi(ij)
            j = indxj(ij)

            n1 = donor(i,j,n)
            worka(i,j) = dvice(i,j,n) / vicen(i,j,n1)
            if (n1 .eq. n) then
               n2 = n1+1
            else                ! n1 = n+1
               n2 = n
            endif
            
            aicen(i,j,n1) = aicen(i,j,n1) - daice(i,j,n)
            aicen(i,j,n2) = aicen(i,j,n2) + daice(i,j,n)
            vicen(i,j,n1) = vicen(i,j,n1) - dvice(i,j,n)
            vicen(i,j,n2) = vicen(i,j,n2) + dvice(i,j,n)
            
            dvsnow = vsnon(i,j,n1) * worka(i,j)
            vsnon(i,j,n1) = vsnon(i,j,n1) - dvsnow
            vsnon(i,j,n2) = vsnon(i,j,n2) + dvsnow
            
            daTsf = daice(i,j,n)*Tsfcn(i,j,n1)
            aTsfn(i,j,n1) = aTsfn(i,j,n1) - daTsf
            aTsfn(i,j,n2) = aTsfn(i,j,n2) + daTsf 

            desnow = esnon(i,j,n1) * worka(i,j)
            esnon(i,j,n1) = esnon(i,j,n1) - desnow
            esnon(i,j,n2) = esnon(i,j,n2) + desnow

         enddo                  ! ij

         do k = 1, nilyr
cdir$ ivdep      !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
            do ij = 1, icells
               i = indxi(ij)
               j = indxj(ij)

               n1 = donor(i,j,n)
               if (n1 .eq. n) then
                  n2 = n1+1
               else             ! n1 = n+1
                  n2 = n
               endif

               deice = eicen(i,j,ilyr1(n1)+k-1) * worka(i,j)
               eicen(i,j,ilyr1(n1)+k-1) = 
     &              eicen(i,j,ilyr1(n1)+k-1) - deice
               eicen(i,j,ilyr1(n2)+k-1) = 
     &              eicen(i,j,ilyr1(n2)+k-1) + deice
            enddo               ! ij
         enddo                  ! k

      enddo                     ! boundaries, 1 to ncat-1

      !-----------------------------------------------------------------
      ! Update ice thickness and temperature
      !-----------------------------------------------------------------

      do n = 1, ncat
         do j = jlo,jhi
         do i = ilo,ihi
            if (aicen(i,j,n) > puny) then
               hicen(i,j,n) = vicen(i,j,n) / aicen(i,j,n)
               Tsfcn(i,j,n) = aTsfn(i,j,n) / aicen(i,j,n)
            else
               hicen(i,j,n) = c0
               Tsfcn(i,j,n) = Tf(i,j)
            endif
         enddo                  ! i
         enddo                  ! j
      enddo                     ! n

      end subroutine shift_ice

!=======================================================================
!BOP
!
! !IROUTINE: column_sum - sum field over all ice categories
!
! !INTERFACE:
!
      subroutine column_sum (nsum, xin, xout)
!
! !DESCRIPTION:
!
! For each grid cell, sum field over all ice categories.
!
! !REVISION HISTORY:
!
! author: William H. Lipscomb, LANL
!
! !USES:
!
! !INPUT/OUTPUT PARAMETERS:
!
      integer (kind=int_kind), intent(in) ::
     &     nsum                              ! number of categories/layers

      real (kind=dbl_kind), intent(in) ::
     &     xin  (imt_local,jmt_local,nsum)   ! input field


      real (kind=dbl_kind), intent(out) ::
     &     xout (imt_local,jmt_local)        ! output field
!
!EOP
!
      integer (kind=int_kind) ::
     &     i, j                 ! horizontal indices
     &,    n                    ! category/layer index

      do j = 1, jmt_local
      do i = 1, imt_local
         xout(i,j) = c0
      enddo
      enddo

      do n = 1, nsum
         do j = 1, jmt_local
         do i = 1, imt_local
            xout(i,j) = xout(i,j) + xin(i,j,n)
         enddo                  ! i
         enddo                  ! j
      enddo                     ! n

      end subroutine column_sum

!=======================================================================
!BOP
!
! !IROUTINE: column_conservation_check
!
! !INTERFACE:
!
      subroutine column_conservation_check (x1, x2, max_err, fieldid)
!
! !DESCRIPTION:
!
! For each physical grid cell, check that initial and final values
! of a conserved field are equal to within a small value.
!
! !REVISION HISTORY:
!
! author: William H. Lipscomb, LANL
!
! !USES:
!
      use ice_exit
!
! !INPUT/OUTPUT PARAMETERS:
!
      real (kind=dbl_kind), intent(in) ::
     &     x1 (imt_local,jmt_local)    ! initial field

      real (kind=dbl_kind), intent(in) ::
     &     x2 (imt_local,jmt_local)    ! final field

      real (kind=dbl_kind), intent(in) :: 
     &     max_err                     ! max allowed error

      character (len=char_len), intent(in) ::
     &     fieldid                     ! field identifier
!
!EOP
!
      integer (kind=int_kind) ::
     &     i, j                        ! horizontal indices      

      logical (kind=log_kind) ::
     &   conserv_err          ! = .true. if conservation check failed

      conserv_err = .false.

      do j = jlo, jhi
      do i = ilo, ihi
         if (abs(x2(i,j)) - abs(x1(i,j)) > max_err) then
            conserv_err = .true.
         endif
      enddo
      enddo

      if ( conserv_err ) then
        do j = jlo, jhi
        do i = ilo, ihi
         if (abs(x2(i,j)) - abs(x1(i,j)) > max_err) then
            write (nu_diag,*) ''
            write (nu_diag,*) 'Conservation error: ', fieldid
            write (nu_diag,*)  my_task, ':', i, j
            write (nu_diag,*) 'Initial value =', x1(i,j)
            write (nu_diag,*) 'Final value =',   x2(i,j)
            write (nu_diag,*) 'Difference =', x2(i,j) - x1(i,j)
            call abort_ice ('ice: Conservation error')
         endif
        enddo
        enddo
      endif

      end subroutine column_conservation_check

!=======================================================================
!BOP
!
! !IROUTINE: zap_small_areas - eliminate very small ice areas
!
! !INTERFACE:
!
      subroutine zap_small_areas
!
! !DESCRIPTION:
!
! For each ice category in each grid cell, remove ice if the fractional
! area is less than puny.
! This subroutine violates mass and energy conservation by a small
!  amount.  At some point, might want to modify so that the mass and
!  energy we throw away is added back somewhere else.
!
! !REVISION HISTORY:
!
! author: William H. Lipscomb, LANL
! Nov 2003:  Modified by Julie Schramm to conserve volume and energy
!
! !USES:
	 use ice_flux, only: Tf
	 use ice_grid, only: tarean, tareas
!
! !INPUT/OUTPUT PARAMETERS:
!
!EOP
!
      integer (kind=int_kind) ::
     &     i,j       ! horizontal indices
     &,    n         ! ice category index
     &,    k         ! ice layer index

      real (kind=dbl_kind) ::
     &     dvicen(1)    ! local zapped ice volume north
     &,    dvsnon(1)    ! local zapped snow volume north
     &,    dvices(1)    ! local zapped ice volume south
     &,    dvsnos(1)    ! local zapped snow volume south
     &,    deicen(1)    ! local zapped ice energy north
     &,    deices(1)    ! local zapped ice energy south
     &,    desnon(1)    ! local zapped snow energy north
     &,    desnos(1)    ! local zapped snow energy south

      dvicen = c0
      dvices = c0
      dvsnon = c0
      dvsnos = c0
      deicen = c0
      deices = c0
      desnon = c0
      desnos = c0

      do n = 1, ncat

         do k = 1, nilyr
cdir$ ivdep      !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
            do j = jlo, jhi
            do i = ilo, ihi
               if (aicen(i,j,n) <= puny) then

                  ! zapped ice energy
                  deicen(1) = deicen(1) 
     &                      + eicen(i,j,ilyr1(n)+k-1)*tarean(i,j)
                  deices(1) = dvices(1) 
     &                      + eicen(i,j,ilyr1(n)+k-1)*tareas(i,j)
                  eicen(i,j,ilyr1(n)+k-1) = c0

               endif            ! aicen < puny
            enddo               ! i
            enddo               ! j
         enddo                  ! k

cdir$ ivdep      !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
         do j = jlo, jhi
         do i = ilo, ihi
            if (aicen(i,j,n) <= puny) then

               ! zapped ice and snow volume, snow energy
               dvicen(1) = dvicen(1) + vicen(i,j,n)*tarean(i,j)
               dvices(1) = dvices(1) + vicen(i,j,n)*tareas(i,j)
               dvsnon(1) = dvsnon(1) + vsnon(i,j,n)*tarean(i,j)
               dvsnos(1) = dvsnos(1) + vsnon(i,j,n)*tareas(i,j)
               desnon(1) = desnon(1) + esnon(i,j,n)*tarean(i,j)
               desnos(1) = desnos(1) + esnon(i,j,n)*tareas(i,j)

               aice0(i,j) = aice0(i,j) + aicen(i,j,n)
               aicen(i,j,n) = c0
               vicen(i,j,n) = c0
               vsnon(i,j,n) = c0
               esnon(i,j,n) = c0
               Tsfcn(i,j,n) = Tf(i,j)

            endif               ! aicen < puny
         enddo                  ! i
         enddo                  ! j

      enddo                     ! n

      call normalize_state (dvicen, dvices, dvsnon, dvsnos,
     &                      deicen, deices, desnon, desnos)

      end subroutine zap_small_areas

!=======================================================================
!BOP
!
! !IROUTINE: normalize_state - normalize ice state after small areas removed
!
! !INTERFACE:
!
      subroutine normalize_state (dvicen, dvices, dvsnon, dvsnos,
     &                            deicen, deices, desnon, desnos)
!
! !DESCRIPTION:
!
! Normalize the state variables so that ice and snow volume and energy 
! are conserved globally after removing small areas
!
! !REVISION HISTORY:
!
! author: C. M. Bitz, UW
!
! Nov 2003: modified by Julie Schramm (NCAR) and William Lipscomb (LANL)
!
! !USES:
!
      use ice_mpi_internal
      use ice_grid
      use ice_calendar, only: istep1, diagfreq
!
! !INPUT/OUTPUT PARAMETERS:
!
      real (kind=dbl_kind), intent(in) ::
     &     dvicen(1)    ! local nuked ice volume north
     &,    dvsnon(1)    ! local nuked snow volume north
     &,    dvices(1)    ! local nuked ice volume south
     &,    dvsnos(1)    ! local nuked snow volume south
     &,    deicen(1)    ! local nuked ice energy north
     &,    deices(1)    ! local nuked ice energy south
     &,    desnon(1)    ! local nuked snow energy north
     &,    desnos(1)    ! local nuked snow energy south
!
!EOP
!
      integer (kind=int_kind) :: i, j, n, k

      real (kind=dbl_kind) :: 
     &   dviceng, dvsnong   ! global nuked ice/snow volumes north
     &,  dvicesg, dvsnosg   ! global nuked ice/snow volumes south
     &,  deiceng, desnong   ! global nuked ice/snow energies north
     &,  deicesg, desnosg   ! global nuked ice/snow energies south
     &,  vitotn,  vstotn    ! global ice/snow volumes north 
     &,  vitots,  vstots    ! global ice/snow volumes south
     &,  eitotn,  estotn    ! global ice/snow volumes north 
     &,  eitots,  estots    ! global ice/snow volumes south
     &,  Rvicen,  Rvices    ! normalization factor for nuked ice volumes
     &,  Rvsnon,  Rvsnos    ! normalization factor for nuked sno volumes
     &,  Reicen,  Reices    ! normalization factor for nuked ice volumes
     &,  Resnon,  Resnos    ! normalization factor for nuked sno volumes

      real (kind=dbl_kind), dimension(ilo:ihi,jlo:jhi) :: 
     &   work1, work2, work3, work4

      !-----------------------------------------------------------------
      ! aggregate volume and energy
      !-----------------------------------------------------------------
      call aggregate

      !-----------------------------------------------------------------
      ! total ice and snow volume
      !-----------------------------------------------------------------
      call get_sum(0,tarean,one,vice,vitotn)
      call get_sum(0,tareas,one,vice,vitots)
      call get_sum(0,tarean,one,vsno,vstotn)
      call get_sum(0,tareas,one,vsno,vstots)

      !-----------------------------------------------------------------
      ! total ice and snow energy
      !-----------------------------------------------------------------
      call get_sum(0,tarean,one,eice,eitotn)
      call get_sum(0,tareas,one,eice,eitots)
      call get_sum(0,tarean,one,esno,estotn)
      call get_sum(0,tareas,one,esno,estots)

      !-----------------------------------------------------------------
      ! Renormalize volume and energy to account for nuked ice and snow. 
      !-----------------------------------------------------------------
      Rvicen = c1
      Rvices = c1
      Rvsnon = c1
      Rvsnos = c1

      Reicen = c1
      Reices = c1
      Resnon = c1
      Resnos = c1

      dviceng = ice_global_real_sum(1,dvicen)
      dvicesg = ice_global_real_sum(1,dvices)
      dvsnong = ice_global_real_sum(1,dvsnon)
      dvsnosg = ice_global_real_sum(1,dvsnos)

      deiceng = ice_global_real_sum(1,deicen)
      deicesg = ice_global_real_sum(1,deices)
      desnong = ice_global_real_sum(1,desnon)
      desnosg = ice_global_real_sum(1,desnos)

      if (dviceng + dvicesg + dvsnong + dvsnosg > puny) then

         if (vitotn > c0) Rvicen = Rvicen + dviceng/vitotn
         if (vitots > c0) Rvices = Rvices + dvicesg/vitots
         if (vstotn > c0) Rvsnon = Rvsnon + dvsnong/vstotn
         if (vstots > c0) Rvsnos = Rvsnos + dvsnosg/vstots

         if (eitotn > c0) Reicen = Reicen + deiceng/eitotn
         if (eitots > c0) Reices = Reices + deicesg/eitots
         if (estotn > c0) Resnon = Resnon + desnong/estotn
         if (estots > c0) Resnos = Resnos + desnosg/estots

         do j=jlo,jhi
         do i=ilo,ihi
            work1(i,j) = Rvicen*mask_n(i,j) + Rvices*mask_s(i,j)
            work2(i,j) = Rvsnon*mask_n(i,j) + Rvsnos*mask_s(i,j)
            work3(i,j) = Reicen*mask_n(i,j) + Reices*mask_s(i,j)
            work4(i,j) = Resnon*mask_n(i,j) + Resnos*mask_s(i,j)
         enddo
         enddo

         do n = 1,ncat
cdir$ ivdep      !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
            do j=jlo,jhi
            do i=ilo,ihi
               vicen(i,j,n) = vicen(i,j,n) * work1(i,j)
               vsnon(i,j,n) = vsnon(i,j,n) * work2(i,j)
               esnon(i,j,n) = esnon(i,j,n) * work4(i,j)
            enddo               ! i
            enddo               ! j
         enddo                  ! n

         do k = 1,ntilay
cdir$ ivdep      !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
            do j=jlo,jhi
            do i=ilo,ihi
               eicen(i,j,k) = eicen(i,j,k) * work3(i,j)
            enddo               ! i
            enddo               ! j
         enddo                  ! k

      endif                     ! nuked volume > puny

      if (mod(istep1,diagfreq)==0 .and. my_task==master_task) then
!!         write(nu_diag,4000) Rvicen, Rvices, Rvsnon, Rvsnos
!!         write(nu_diag,4001) Reicen, Reices, Resnon, Resnos
      endif 
 4000 format('Rnormv ', 4((1pe20.13,1x)))
 4001 format('Rnorme ', 4((1pe20.13,1x)))

      end subroutine normalize_state

!=======================================================================

      end module ice_itd

!=======================================================================
