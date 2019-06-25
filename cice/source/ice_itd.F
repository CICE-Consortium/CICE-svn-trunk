c $Id: $
c=======================================================================
!---! Routines to initialize the ice thickness distribution and 
!---! utilities to redistribute ice among categories. These routines 
!---! are not specific to a particular numerical implementation.
!---!
!---! authors C. M. Bitz
!---!         Elizabeth C. Hunke
!---!         William H. Lipscomb
c=======================================================================

      module ice_itd

      use ice_kinds_mod
      use ice_model_size
      use ice_constants
      use ice_state
      use ice_fileunits

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
     &   hs_min             ! minimum snow thickness allowed (m)
     &,  hin_max(0:ncat)    ! category limits                (m)
     &,  ain_min(0:ncat)    ! minimum fract. area allowed in each cat

      ! extra array used only for holding history information
      real (kind=dbl_kind), dimension (ilo:ihi,jlo:jhi) ::
     &   snoice             ! snow-ice formation       (m)

!---!-------------------------------------------------------------------
!---! a note regarding hi_min and hin_max(0):
!---! both represent a minimum ice thickness.  hin_max(0) is
!---! intended to be used for particular numerical implementations
!---! of category conversions in the ice thickness distribution.
!---! hi_min is a more general purpose parameter, but is specifically 
!---! for maintaining stability in the thermodynamics.  Currently,
!---! hi_min = 0.1 m
!---! hin_max(0) = 0.1 m for the delta function itd
!---! hin_max(0) = 0.0 m for linear remapping
!---!
!---! similarly, there are two values of minimum snow thickness
!---! (the other is defined in ice_vthermo.H since it is used only
!---! for thermo.) 
!---!
!---! Also note that the upper limit on the thickest category
!---! is only used for the linear remapping scheme
!---! and it is not a true upper limit on the thickness
!---!-------------------------------------------------------------------

c=======================================================================

      contains

c=======================================================================

      subroutine init_itd

!---!-------------------------------------------------------------------
!---! Initialize area fraction and thickness boundaries for the itd model
!---!
!---! authors William H. Lipscomb
!---!         Elizabeth C. Hunke and C. M. Bitz
!---!-------------------------------------------------------------------

      integer (kind=int_kind) :: nc
      real (kind=dbl_kind):: cc1, cc2, cc3, x1

      if ((ncat.eq.1) .and. (kitd.eq.1)) then
        write(nu_diag,*) 'Remapping the ITD is not allowed for ncat=1.'
        write(nu_diag,*) 'Use the delta function ITD option instead.'
        stop
      endif

      ! minimum allowable fract
      ! test carefully when changing these values
      ! round off error can cause serious problems if
      ! ain_min is too small but it should be as small
      ! as possible for accuracy
      ain_min(0)=1.0e-6_dbl_kind
      do nc=1,ncat
        ain_min(nc)=5.0e-6_dbl_kind
      enddo

      ! minimum allowable snow thickness
      hs_min = 0.00001_dbl_kind

      if (kitd.eq.1) then
        ! linear remapping itd category limits
        cc1 = c3/real(ncat)
        cc2 = c15*cc1      
        cc3 = c3

        hin_max(0) = c0  ! minimum ice thickness, m
      else
        ! delta function itd category limits
        cc1 = max(1.1_dbl_kind/ncat,c1*hi_min)
        cc2 = c25*cc1
        cc3 = 2.25_dbl_kind

        ! hin_max(0) should not be zero 
        ! use some caution in making it less than 0.10
        hin_max(0) = hi_min  ! minimum ice thickness, m
      endif

      do nc = 1, ncat
       x1 = real(nc-1) / real(ncat)
       hin_max(nc) = hin_max(nc-1) + cc1 + cc2*(c1 + tanh(cc3*(x1-c1)))
      enddo

      if (my_task.eq.master_task) then
      write(nu_diag,*) ''
      write(nu_diag,*) 'hin_max(nc-1) < Cat nc < hin_max(nc)'
      do nc = 1, ncat
         write(nu_diag,*) hin_max(nc-1),' < Cat ',nc, ' < ',hin_max(nc)
      enddo
      endif

      end subroutine init_itd

c=======================================================================

      subroutine aggregate

!---!-------------------------------------------------------------------
!---! aggregate ice state over the grid
!---!
!---! author C. M. Bitz
!---!-------------------------------------------------------------------

      use ice_domain
      use ice_flux
      use ice_grid

      integer (kind=int_kind) :: i, j, k, nc

      do j=jlo,jhi
        do i=ilo,ihi
          aice0(i,j) = c1
          aice(i,j) = c0
          vice(i,j) = c0
          vsno(i,j) = c0
          eice(i,j) = c0
          esno(i,j) = c0
          Tsfc(i,j) = Tf(i,j)
          if (tmask(i,j)) then
            Tsfc(i,j) = c0
            do nc = 1,ncat
               aice(i,j) = aice(i,j) + aicen(i,j,nc)
               vice(i,j) = vice(i,j) + vicen(i,j,nc)
               vsno(i,j) = vsno(i,j) + vsnon(i,j,nc)
               do k = 1, nilyr
                  eice(i,j) = eice(i,j) + eicen(i,j,(nc-1)*nilyr+k)
               enddo
               esno(i,j) = esno(i,j) + esnon(i,j,nc)
               Tsfc(i,j) = Tsfc(i,j) + Tsfcn(i,j,nc)*aicen(i,j,nc)
            enddo
            if (aice(i,j).gt.0.)  then
              aice0(i,j)=c1-aice(i,j)
              Tsfc(i,j) = Tsfc(i,j)/aice(i,j)
            else
              Tsfc(i,j) = Tf(i,j)
            endif
          endif ! tmask
        enddo   ! end of loop over i,j
      enddo

      call bound(aice0)
      call bound(aice)
      call bound(vice)
      call bound(vsno)
      call bound(eice)
      call bound(esno)
      call bound(Tsfc)

      end subroutine aggregate

c=======================================================================

      subroutine aggregate_pt(i,j)

!---!-------------------------------------------------------------------
!---! aggregate ice thickness distribution state 
!---!
!---! author C. M. Bitz
!---!-------------------------------------------------------------------

      use ice_flux

      integer (kind=int_kind), intent(in) :: i, j
      integer (kind=int_kind) :: k, nc

            aice(i,j) = c0
            vice(i,j) = c0
            vsno(i,j) = c0
            eice(i,j) = c0
            esno(i,j) = c0
            Tsfc(i,j) = c0
            do nc = 1,ncat
               aice(i,j) = aice(i,j) + ain(nc)
               vice(i,j) = vice(i,j) + vin(nc)
               vsno(i,j) = vsno(i,j) + vsn(nc)
               do k = 1, nilyr
                  eice(i,j) = eice(i,j) + ein(k,nc)
               enddo
               esno(i,j) = esno(i,j) + esn(nc)
               Tsfc(i,j) = Tsfc(i,j) + Tsfn(nc)*ain(nc)
            enddo
            if (aice(i,j).gt.0.)  then
              aice0(i,j)= c1-aice(i,j)
              Tsfc(i,j) = Tsfc(i,j)/aice(i,j)
            else
              Tsfc(i,j) = Tf(i,j)
            endif

      end subroutine aggregate_pt

c=======================================================================

      subroutine to_column(i,j)

!---!-------------------------------------------------------------------
!---! define column state variables for given point (i,j)
!---!
!---! authors Elizabeth C. Hunke
!---!         C. M. Bitz
!---!-------------------------------------------------------------------

      integer (kind=int_kind), intent(in) :: i,j
      integer (kind=int_kind) :: nc, layer

        ai0 = aice0(i,j)
        do nc = 1,ncat
          ain(nc)  = aicen(i,j,nc)
          vin(nc)  = vicen(i,j,nc)
          vsn(nc)  = vsnon(i,j,nc)
          Tsfn(nc) = Tsfcn(i,j,nc)
          do layer = 1,nilyr
            ein(layer,nc) = eicen(i,j,ilyr1(nc)+layer-1)
          enddo
          esn(nc) = esnon(i,j,nc)
        enddo

      end subroutine to_column

c=======================================================================

      subroutine from_column(i,j)

!---!-------------------------------------------------------------------
!---! reload state variables for given point (i,j) from column variables
!---!
!---! authors Elizabeth C. Hunke
!---!         C. M. Bitz
!---!-------------------------------------------------------------------

      integer (kind=int_kind), intent(in) :: i,j
      integer (kind=int_kind) :: nc, layer

        aice0(i,j) = ai0
        do nc = 1,ncat
          aicen(i,j,nc) = ain(nc)  
          vicen(i,j,nc) = vin(nc)  
          vsnon(i,j,nc) = vsn(nc)  
          Tsfcn(i,j,nc) = Tsfn(nc) 
          do layer = 1,nilyr
             eicen(i,j,ilyr1(nc)+layer-1) = ein(layer,nc)
          enddo
          esnon(i,j,nc) = esn(nc)
        enddo

      end subroutine from_column

c=======================================================================

      subroutine check_state(i,j)

!---!-------------------------------------------------------------------
!---! Insist that certain fields are monotone.
!---!
!---! authors C. M. Bitz
!---!         Elizabeth C. Hunke
!---!-------------------------------------------------------------------

      use ice_flux

      integer (kind=int_kind) :: i,j,layer,nc
      logical (kind=log_kind) :: zerout

      do nc=1,ncat

         zerout=.false.
         if ((aicen(i,j,nc) .lt. puny).or.(vicen(i,j,nc) .lt. puny))
     $        zerout = .true.
         do layer=ilyr1(nc),ilyrn(nc)
            if (eicen(i,j,layer) .gt. -puny) zerout=.true.
         enddo

         if (zerout) then

           esnon(i,j,nc) = c0
           do layer = ilyr1(nc),ilyrn(nc)
             eicen(i,j,layer) = c0
           enddo

           aice0(i,j) = aice0(i,j)+aicen(i,j,nc) 
           aicen(i,j,nc) = c0
           vicen(i,j,nc) = c0
           vsnon(i,j,nc) = c0
           Tsfcn(i,j,nc) = Tf(i,j)

         elseif (vsnon(i,j,nc) .lt. puny) then

           vsnon(i,j,nc) = c0
           esnon(i,j,nc) = c0

         endif

       if (vsnon(i,j,nc).ge.puny) then
         if (-esnon(i,j,nc)/vsnon(i,j,nc)-Lfresh*rhos.lt.eps04) 
     &     esnon(i,j,nc) = -vsnon(i,j,nc)*(Lfresh*rhos + eps04)
       endif

      enddo ! categories

      if (aice0(i,j).lt.puny) aice0(i,j) = c0

      end subroutine check_state

c=======================================================================

      subroutine distr_check

!---!-------------------------------------------------------------------
!---! Force ice thickness distribution to maintain two rules
!---! (1) the fractional area cannot be less than some limiting value
!---! (2) each categories thickness lies withing the 
!---!     max and min thickness range for that category 
!---!
!---! authors C. M. Bitz
!---!         Elizabeth C. Hunke
!---!-------------------------------------------------------------------

      use ice_flux
      use ice_grid

      integer (kind=int_kind) :: i, j, nc, nr, layer
      real (kind=dbl_kind) :: 
     &   hin    (ncat)       !  ice thickness for each cat (m)
     &,  hsn    (ncat)       ! snow thickness for each cat (m)


      do j=jlo,jhi
       do i=ilo,ihi
        if (tmask(i,j)) then

         call to_column(i,j)

      !-----------------------------------------------------------------
      ! Put tiny amounts of open water area out of misery by 
      ! adding to the thinnest ice provided it has ice already
      !-----------------------------------------------------------------
         if (ai0.lt.ain_min(0).and.ai0.ge.c0) then
          do nc = 1,ncat
           if (ain(nc).gt.c0) then
            ain(nc) = ain(nc)+ai0
            ai0 = c0
           endif
          enddo
         endif

      !-----------------------------------------------------------------
      ! Get rid of categories with very small areas or volumes
      !-----------------------------------------------------------------
         do nc = 1,ncat
          if ((ain(nc).gt.c0) .and.
     &       ((ain(nc).le.ain_min(nc)) .or. vin(nc).lt.puny)) then
           ai0 = ai0 + ain(nc)
           call zerocat(Tf(i,j),ain(nc),vin(nc),vsn(nc),
     $                ein(:,nc),esn(nc),Tsfn(nc),hsn(nc),hin(nc))
          elseif (vsn(nc) .lt. hs_min*ain(nc)) then
           vsn(nc) = c0
           esn(nc) = c0
          endif
         enddo

         call rebin(i,j)
         call from_column(i,j)

        endif ! tmask
       enddo  
      enddo

      end subroutine distr_check

c=======================================================================

      subroutine rebin(i,j)

!---!-------------------------------------------------------------------
!---! Rebins thicknesses into defined categories
!---!
!---! authors William H. Lipscomb
!---!         Elizabeth C. Hunke
!---!-------------------------------------------------------------------

      use ice_flux

      integer (kind=int_kind), intent(in) ::
     &   i,j

      integer (kind=int_kind) ::
     &   nc                 ! category index
     &,  shiftflag          ! = 1 if ice must be shifted
     &,  donor(ncat-1)      ! donor category index

      real (kind=dbl_kind) ::
     &   hin(ncat)          ! ice thickness for each cat        (m)
     &,  hsn(ncat)          ! ice thickness for each cat        (m)
     &,  daice(ncat-1)      ! ice area transferred
     &,  dvice(ncat-1)      ! ice volume transferred

      if (ncat.gt.1) then

      !-----------------------------------------------------------------
      ! Update thicknesses
      !-----------------------------------------------------------------
      do nc = 1, ncat
          hin(nc) = c0
          hsn(nc) = c0
          if (ain(nc).gt.eps11) then
            hin(nc) = vin(nc)/ain(nc)
            hsn(nc) = vsn(nc)/ain(nc)
          endif
      enddo
      ! make sure thickness of cat 1 is at least hin_max(0)
      if ((hin(1).le.hin_max(0)).and.(hin_max(0).gt.c0)) then
        ain(1) = vin(1)/hin_max(0)
        hin(1) = hin_max(0)
      endif

      !-----------------------------------------------------------------
      ! If a category thickness is not in bounds, shift the
      ! entire area and volume to the neighboring category
      !-----------------------------------------------------------------

      ! Move thin categories up
      shiftflag = 0
      do nc = 1, ncat-1         ! loop over category boundaries
         if (ain(nc).gt.puny .and. hin(nc).gt.hin_max(nc)) then
            shiftflag = 1
            donor(nc) = nc
            daice(nc) = ain(nc)
            dvice(nc) = vin(nc)
         else
            donor(nc) = 0
            daice(nc) = c0
            dvice(nc) = c0
         endif
      enddo
      if (shiftflag.eq.1) call shift_ice(donor,daice,dvice,hin,hsn,i,j)

      ! Move thick categories down
      shiftflag = 0
      do nc = ncat-1, 1, -1     ! loop over category boundaries
         if (ain(nc+1).gt.puny .and. hin(nc+1).le.hin_max(nc)) then
            shiftflag = 1
            donor(nc) = nc+1
            daice(nc) = ain(nc+1)
            dvice(nc) = vin(nc+1)
         else
            donor(nc) = 0
            daice(nc) = c0
            dvice(nc) = c0
         endif
      enddo
      if (shiftflag.eq.1) call shift_ice(donor,daice,dvice,hin,hsn,i,j)

      ! make sure thickness of cat 1 is at least hin_max(0)
      if ((hin(1).le.hin_max(0)).and.(hin_max(0).gt.c0)) then
        ain(1) = vin(1)/hin_max(0)
        hin(1) = hin_max(0)
      endif

      endif ! ncat > 1
          
      end subroutine rebin

c=======================================================================

      subroutine zerocat(Tf0, ai, vi,  vs, ei, es, Tsf, hs, hi)

!---!-------------------------------------------------------------------
!---! reset category variables in the case of no ice
!---!
!---! author C. M. Bitz
!---!-------------------------------------------------------------------

      real (kind=dbl_kind), intent(in) ::  Tf0      ! freezing temperature

      real (kind=dbl_kind), intent(inout) ::
     &   ai           ! fraction of ice
     &,  vi           ! volume per unit area of ice             (m)
     &,  vs           ! volume per unit area of snow            (m)
     &,  hi           ! ice thickness                           (m)
     &,  hs           ! snow thickness                          (m)
     &,  Tsf          ! ice/snow top surf. temperature          (K)
     &,  ei (nilyr)   ! energy of melting of ice per layer  (J/m^2)
     &,  es           ! energy of melting of snow           (J/m^2)

      integer (kind=int_kind) :: layer

      vi = c0
      vs = c0
      ai = c0
      hs = c0
      hi = c0
      do layer = 1,nilyr
         ei(layer) = c0
      enddo
      es = c0
      Tsf = Tf0

      end subroutine zerocat

c=======================================================================

      subroutine shift_ice (donor, daice, dvice, hin, hsn, i,j)

!---!-------------------------------------------------------------------
!---! Shift ice across category boundaries, conserving area, volume, and
!---! energy.
!---!-------------------------------------------------------------------

      integer (kind=int_kind), intent(in) ::
     &   i,j
     &,  donor(ncat-1)     ! donor category index
      real (kind=dbl_kind), intent(inout) ::
     &   daice(ncat-1)     ! ice area transferred across boundary
     &,  dvice(ncat-1)     ! ice volume transferred across boundary
     &,  hin(ncat)         !  ice thickness for each cat        (m)
     &,  hsn(ncat)         ! snow thickness for each cat        (m)

      integer (kind=int_kind) ::
     &   n                 ! thickness category index
     &,  n2                ! receiver category
     &,  n1                ! donor category
     &,  layer             ! ice layer index

      real (kind=dbl_kind) ::
     &   fvolume           ! fraction of category volume transferred
     &,  dvsnow            ! snow volume transferred
     &,  desnow            ! snow energy transferred
     &,  deice             ! ice energy transferred
     &,  aTsfn(ncat)
     &,  daTsf

      !-----------------------------------------------------------------
      ! Check for rounding errors
      !-----------------------------------------------------------------
      ! Note: daice < 0 or dvice < 0 usually happens when category n 
      ! has a small area, with hice(n) very close to a boundary.  Then 
      ! the coefficients of g(h) are large, and the computed daice and 
      ! dvice can be in error. If this happens, it is best to transfer 
      ! either the entire category or nothing at all, depending on which 
      ! side of the boundary hice(n) lies.
      !-----------------------------------------------------------------

      do n = 1, ncat-1
         if (donor(n) .ne. 0) then 
         n1 = donor(n)

         if (daice(n) .lt. c0) then
            if (abs(daice(n)) .lt. eps11) then   
               if ( (n1.eq.n .and. hin(n1).gt.hin_max(n))
     &                        .or.
     &            (n1.eq.n+1 .and. hin(n1).le.hin_max(n)) ) then
                  daice(n) = ain(n1)   ! shift entire category
                  dvice(n) = vin(n1)
               else             
                  daice(n) = c0            ! shift no ice
                  dvice(n) = c0
               endif
            else
             write(nu_diag,*) my_task,':',i,j,'ITD Neg da =',daice(n),
     &                                   ' boundary',n
             stop
            endif
         endif
         
         if (dvice(n) .lt. c0) then
            if (abs(dvice(n)) .lt. eps11) then   
               if ( (n1.eq.n .and. hin(n1).gt.hin_max(n))
     &                        .or.
     &            (n1.eq.n+1 .and. hin(n1).le.hin_max(n)) ) then
                  daice(n) = ain(n1)   ! shift entire category
                  dvice(n) = vin(n1)
               else             
                  daice(n) = c0            ! shift no ice
                  dvice(n) = c0
               endif
            else
             write(nu_diag,*) my_task,':',i,j,
     &             'ITD Neg dvice =',dvice(n),' cat',n
             stop
            endif
         endif

            ! If daice is close to ain, set daice = ain.
         if (daice(n) .gt. ain(n1)-eps11) then
            if (daice(n) .lt. ain(n1)+eps11) then
               daice(n) = ain(n1)
               dvice(n) = vin(n1)
            else
               write(nu_diag,*) my_task,':',i,j,'ITD da > ain, cat',n1
               write(nu_diag,*) my_task,':',i,j,'da =', daice(n),
     &                                     'ain =',ain(n1)
               stop
            endif
         endif    

         if (dvice(n) .gt. vin(n1)-eps11) then
            if (dvice(n) .lt. vin(n1)+eps11) then   
               daice(n) = ain(n1)
               dvice(n) = vin(n1)
            else
               write(nu_diag,*) my_task,':',i,j,'ITD dv > vin, cat',n1
               write(nu_diag,*) my_task,':',i,j,'dvice =', dvice(n),
     &                                     'vin =', vin(n1)
               stop
            endif
         endif

         endif   ! donor .ne. 0 
      enddo    ! loop over boundaries            

      !-----------------------------------------------------------------
      ! transfer volume and energy between categories
      !-----------------------------------------------------------------
      do n = 1, ncat
        aTsfn(n) = ain(n)*Tsfn(n)
      enddo

      do n = 1, ncat-1
         if (daice(n) .gt. c0) then   ! daice(n) can be < eps11
            n1 = donor(n)
            fvolume = dvice(n) / vin(n1)
            if (n1 .eq. n) then
               n2 = n1+1
            else  ! n1 = n+1
               n2 = n
            endif

            ain(n1) = ain(n1) - daice(n)
            ain(n2) = ain(n2) + daice(n)
            vin(n1) = vin(n1) - dvice(n)
            vin(n2) = vin(n2) + dvice(n)

            dvsnow = vsn(n1) * fvolume
            vsn(n1) = vsn(n1) - dvsnow
            vsn(n2) = vsn(n2) + dvsnow

            daTsf = daice(n)*Tsfn(n1)
            aTsfn(n1) = aTsfn(n1) - daTsf
            aTsfn(n2) = aTsfn(n2) + daTsf 

            desnow = esn(n1) * fvolume
            esn(n1) = esn(n1) - desnow
            esn(n2) = esn(n2) + desnow

            do layer = 1,nilyr
             deice = ein(layer,n1) * fvolume
             ein(layer,n1) = ein(layer,n1) - deice
             ein(layer,n2) = ein(layer,n2) + deice
            enddo

         endif  ! daice > 0
      enddo   ! loop over boundaries

      !-----------------------------------------------------------------
      ! Update thicknesses
      !-----------------------------------------------------------------
      do n = 1, ncat
          hin(n) = c0
          hsn(n) = c0
          if (ain(n).gt.eps11) then
            hin(n) = vin(n)/ain(n)
            hsn(n) = vsn(n)/ain(n)
            Tsfn(n) = aTsfn(n)/ain(n)
          endif
      enddo

      end subroutine shift_ice

c=======================================================================

      subroutine freeboard

!---!-----------------------------------------------------------------
!---! If there is enough snow to lower the ice/snow interface below 
!---! freeboard, convert enough snow to ice to bring the interface back 
!---! to sea-level.  Spread energy evenly among ice layers.
!---!----------------------------------------------------------------

      integer (kind=int_kind) :: i,j, nc, k
      real (kind=dbl_kind) ::
     &   hin(ncat) ! ice thickness for each cat  (m)
     &,  hsn(ncat) ! snow thickness for each cat (m)
     &,  dhi       ! change in ice thickness (m)
     &,  dhs       ! change in snow thickness (m)
     &,  de        ! energy transferred
     &,  dz        ! distance freeboard below SL (m)

      do j=jlo,jhi
       do i=ilo,ihi
         call to_column(i,j)
          do nc = 1,ncat
           hin(nc) = c0
           hsn(nc) = c0
           if (ain(nc).gt.puny) then
             hin(nc) = vin(nc)/ain(nc)
             hsn(nc) = vsn(nc)/ain(nc)
           endif

           dz = hsn(nc) - hin(nc)*(rhow-rhoi)/rhos
           if (dz .gt. puny .and. hsn(nc).gt.hs_min) then
            dhs = min(dz*rhoi/rhow, hsn(nc))
            dhi = dhs*rhos/rhoi
            snoice(i,j) = dhi    ! for history diagnostics

            ! energy transferred from snow to ice
            de = esn(nc)*dhs/hsn(nc)
            esn(nc) = esn(nc) - de
            de = de/real(nilyr)   ! energy transferred to each ice layer
            do k = 1,nilyr
             ein(k,nc) = ein(k,nc) + de
            enddo

            hin(nc) = hin(nc) + dhi
            hsn(nc) = hsn(nc) - dhs
            vin(nc) = hin(nc)*ain(nc)
            vsn(nc) = hsn(nc)*ain(nc)
           endif
          enddo
         call rebin(i,j)
         call from_column(i,j)
       enddo
      enddo

      end subroutine freeboard

c=======================================================================

      end module ice_itd

c=======================================================================
