c $Id: $
c=======================================================================
!---! Energy-conserving sea ice model
!---! Driver for thermodynamics and associated changes to the itd
!---!
!---! authors C. M. Bitz
!---!         W. H. Lipscomb
!---!         E. C. Hunke
!---!
!---! See Bitz, C.M., and W.H. Lipscomb, 1999: 
!---! An energy-conserving thermodynamic model of sea ice,
!---! J. Geophys. Res., 104, 15,669-15,677. 
c=======================================================================

      module ice_therm_driver

      use ice_kinds_mod
      use ice_model_size
      use ice_constants
      use ice_state
      use ice_flux
      use ice_diagnostics
      use ice_calendar
      use ice_therm_cice
      use ice_itd

      implicit none

      real (kind=dbl_kind), private, save ::
     &   dhice1(ilo:ihi,jlo:jhi)   ! ice thickness change, 1-category case (m)
     &,  hicen(ncat,ilo:ihi,jlo:jhi)       ! ice thickness (m)
     &,  hsnon(ncat,ilo:ihi,jlo:jhi)       ! snow thickness (m)
     &,  dhicen(ncat,ilo:ihi,jlo:jhi)      ! ice thickness change (m)
     &,  dhsnon(ncat,ilo:ihi,jlo:jhi)      ! snow thickness change (m)
     &,  qicen(nilyr,ncat,ilo:ihi,jlo:jhi) ! enthalpy per unit volume (J/m^3)
     &,  qsnon(ncat,ilo:ihi,jlo:jhi)  ! snow enthalpy per unit volume (J/m^3)

c=======================================================================

      contains

c=======================================================================

      subroutine thermo_rates

!---!-------------------------------------------------------------------
!---! compute thermodynamic growth rates and atmospheric fluxes
!---!-------------------------------------------------------------------

      use ice_domain
      use ice_timers
      use ice_grid
      use ice_ocean

      integer (kind=int_kind) :: i, j, nc

      real (kind=dbl_kind) ::
     &   Fbot            ! heat flx to ice bottom          (W/m^2)
     &,  Frsh_init       ! initial Fresh, use to get correct budgets
     &,  Fh_init         ! initial Fhnet, use to get correct budgets
     &,  Fswthru_init    ! initial Fswthru
     &,  Fsalt_init      ! initial Fsalt
     &,  hin(ncat)       ! ice thickness for each cat        (m)
     &,  hsn(ncat)       ! snow thickness for each cat       (m)

      call ice_timer_start(4)  ! column model
      do j=jlo,jhi
       do i=ilo,ihi

      !-----------------------------------------------------------------
      ! initialize 
      !-----------------------------------------------------------------
        Frsh_init = Fresh(i,j)  
        Fh_init   = Fhnet(i,j)
        Fswthru_init = Fswthru(i,j)
        Fsalt_init = Fsalt(i,j)
        call init_flux_atm(i,j)     ! atm fluxes
        call init_column_diags(i,j) ! init diagnostic var

        call ice_timer_start(5)     ! thermodynamics
        if (tmask(i,j)) then        ! ocean T cells

#ifdef oceanmixed
          call mixed_layer(i,j)     ! get ocean forcing if desired
#endif

          call to_column(i,j)       ! define column state variables
          do nc = 1,ncat
          hicen(nc,i,j) = c0
          hsnon(nc,i,j) = c0
           if (ain(nc).gt.puny) then
            hicen(nc,i,j) = vin(nc)/ain(nc)
            hsnon(nc,i,j) = vsn(nc)/ain(nc)
           endif
          enddo ! ncat

          ! initialize new ice growth and melt (thermo_vertical needs Fbot)
          call init_frzmlt(i,j,Fbot)

      !-----------------------------------------------------------------
      ! vertical thermodynamics: growth rates and fluxes for coupling
      !-----------------------------------------------------------------
          if (aice(i,j).gt.puny) then

            call thermo_cice(i,j,Fbot,
     &      dhice1(i,j),dhicen(:,i,j),dhsnon(:,i,j),
     &      Fhocn(i,j),qicen(:,:,i,j),qsnon(:,i,j))

            ! divide fields by aice
            strairxT(i,j)   = strairxT(i,j)   /aice(i,j)
            strairyT(i,j)   = strairyT(i,j)   /aice(i,j)
            Tref(i,j)      = Tref(i,j)      /aice(i,j)  

            ! history and diagnostics--first part of timestep 
            ! *_init is from last part of previous timestep
            Frsh_hist(i,j) = Fresh(i,j) - Frsh_init
            Foht_hist(i,j) = Fhnet(i,j) - Fh_init
            Fswthru_hist(i,j) = Fswthru(i,j) - Fswthru_init
            Fsalt_hist(i,j) = Fsalt(i,j) - Fsalt_init

#ifdef oceanmixed
            ! for restarting with mixed layer model
            Fhocn(i,j) = Fhocn(i,j) + Fswthru_hist(i,j) 

            ! update mixed layer with heat from ice
            sst(i,j) = sst(i,j) + (Fhnet(i,j)+Fswthru(i,j))*dt/cphm 
            sst_mixedlayer(i,j) = sst(i,j)  ! save value
#endif
          endif

          call from_column(i,j) ! update Tsfc
        endif ! tmask
        call ice_timer_stop(5)  ! thermodynamics
       enddo
      enddo
      call ice_timer_stop(4)    ! column model

      end subroutine thermo_rates

c=======================================================================

      subroutine thermo_itd

!---!-------------------------------------------------------------------
!---! changes to ice thickness distribution associated with 
!---! thermodynamic growth rates
!---! NOTE:  ocean fluxes are initialized here
!---!-------------------------------------------------------------------

      use ice_domain
      use ice_timers
      use ice_history
      use ice_grid
      use ice_itd_linear

      integer (kind=int_kind) :: i, j, nc, layer

      real (kind=dbl_kind) ::
     &   hin(ncat)   ! ice thickness for each cat        (m)
     &,  hsn(ncat)   ! snow thickness for each cat       (m)
     &,  Fnew        ! heat flx to open water        (W/m^2)

      call ice_timer_start(4)  ! column model
      do j=jlo,jhi
       do i=ilo,ihi

      !-----------------------------------------------------------------
      ! initialize 
      !-----------------------------------------------------------------
        aice_init(i,j) = aice(i,j)         ! use for correct diagnostics
        call init_flux_ocn(i,j)            ! ocean fluxes for coupler
        Fresh(i,j) = Frain(i,j)*aice(i,j)  ! rain drains thru to ocn

        if (tmask(i,j)) then
         call to_column(i,j)
         call ice_timer_start(5)   ! thermodynamics

      !-----------------------------------------------------------------
      ! update ice state
      !-----------------------------------------------------------------
         do nc = 1,ncat
           if (ain(nc).gt.puny) then
             hin(nc) = hicen(nc,i,j) + dhicen(nc,i,j)
             hsn(nc) = hsnon(nc,i,j) + dhsnon(nc,i,j)
             if (hin(nc).lt.puny) then
               call zerocat(Tf(i,j),ain(nc),vin(nc),vsn(nc),ein(:,nc),
     &                          esn(nc),Tsfn(nc),hsn(nc),hin(nc))
             else
               vin(nc) = hin(nc)*ain(nc)
               vsn(nc) = hsn(nc)*ain(nc)
               do layer = 1,nilyr
                 ein(layer,nc) = qicen(layer,nc,i,j)*vin(nc)/nilyr
               enddo
               esn(nc) = qsnon(nc,i,j)*vsn(nc)
             endif
           else
               ai0 = ai0 + ain(nc)
               call zerocat(Tf(i,j),ain(nc),vin(nc),vsn(nc),ein(:,nc),
     &                          esn(nc),Tsfn(nc),hsn(nc),hin(nc))
           endif
         enddo                    ! loop over bins
         call ice_timer_stop(5)   ! thermodynamics

!         plabel = 'pre linear itd'
!         if (istep1.ge.check_step.and.i.eq.ip.and.j.eq.jp
!     &       .and.my_task.eq.mtask) 
!     &   call print_state(plabel,i,j)

      !-----------------------------------------------------------------
      ! ice thickness redistribution
      !-----------------------------------------------------------------
         call ice_timer_start(7)  ! category conversions (advection in h)
         if (kitd.eq.1) then
           call linear_itd(i,j,dhicen(:,i,j),hicen(:,i,j),hin,hsn)
         endif
         call ice_timer_stop(7)   ! category conversions 

      !-----------------------------------------------------------------
      ! frazil ice growth 
      !-----------------------------------------------------------------
         call ice_timer_start(5)  ! thermodynamics
         call add_new_ice(i,j,Fnew,hin,hsn)
         call ice_timer_stop(5)   ! thermodynamics

      !-----------------------------------------------------------------
      ! ice thickness redistribution
      !-----------------------------------------------------------------
         call ice_timer_start(7)  ! category conversions (advection in h)
         if (ncat.eq.1) then  
           call reduce_area(i,j,dhice1(i,j),hin,hsn)
         else
           call rebin(i,j)
         endif
         call ice_timer_stop(7)   ! category conversions 

      !-----------------------------------------------------------------
      ! update cell values 
      !-----------------------------------------------------------------
         call aggregate_pt(i,j)
         ai0 = c1-aice(i,j)

#ifdef oceanmixed
         ! net ocean-ice heat flux for restarting mixed layer model
         Fhocn(i,j) = Fhocn(i,j) + Fhnet(i,j) + Fnew
#endif

      !-----------------------------------------------------------------
      ! reload variables
      !-----------------------------------------------------------------
         call from_column(i,j)
         
      !-----------------------------------------------------------------
      ! update history and diagnostic fields
      !-----------------------------------------------------------------
         ! ice concentration & volume tendency
         daidtd(i,j) = aice(i,j) ! temporarily used for initial conc.
         dvidtt(i,j) = (vice(i,j)-dvidtt(i,j))/dt
         dvidtd(i,j) = vice(i,j) ! temporarily used for initial volume
         ! fresh water and heat fluxes over current timestep
         Frsh_hist(i,j) = Frsh_hist(i,j) + Fresh(i,j)
         Foht_hist(i,j) = Foht_hist(i,j) + Fhnet(i,j)

         else
           Fhocn(i,j) = c0       ! for ocean mixed layer model
         endif                   ! tmask

      call check_state(i,j) 

       enddo
      enddo
      call ice_timer_stop(4)     ! column model

      end subroutine thermo_itd

c=======================================================================

      subroutine init_column_diags(i,j)

!---!-------------------------------------------------------------------
!---! initialize diagnostic and history variables
!---!-------------------------------------------------------------------

      use ice_history

      integer (kind=int_kind) :: i, j

        ! for internal calculations and/or coupling
        Fhocn  (i,j) = c0        ! for restarting with mixed layer model
        ! diagnostics for history
        congel (i,j) = c0
        frazil (i,j) = c0
        snoice (i,j) = c0
        meltb  (i,j) = c0
        meltt  (i,j) = c0
        sabs   (i,j) = c0
        dvidtt (i,j) = vice(i,j)
        dvidtd (i,j) = c0
        daidtd (i,j) = c0
        Frsh_hist(i,j) = c0
        Foht_hist(i,j) = c0
        Fswthru_hist(i,j) = c0

      end subroutine init_column_diags

c=======================================================================

      subroutine init_frzmlt(i,j,Fbot)

!---!-------------------------------------------------------------------
!---! initialize ocean-ice heat fluxes
!---! Assuming frzmlt is per grid box area
!---!-------------------------------------------------------------------

      integer (kind=int_kind), intent(in) :: i, j
      real (kind=dbl_kind), intent(out) ::
     &   Fbot        ! heat flx to ice bottom,         (W/m^2)

      real (kind=dbl_kind) ::
     &   ustar       ! skin friction velocity for Fbot (m/s)
     &,  deltaT      ! SST - Tf > 0

      integer (kind=int_kind) :: nc, layer

      ! Parameters for basal and lateral heat flx 
      ! 0.006 = unitless param for basal heat flx ala McPhee and Maykut
      real (kind=dbl_kind), parameter :: 
     &   cpchr =-cp_ocn*rhow*0.006_dbl_kind
#ifdef coupled
     &,  ustar_min = 5.e-3_dbl_kind
#endif
#ifndef coupled
     &,  ustar_min = 5.e-2_dbl_kind   ! for zero currents
#endif

          if (frzmlt(i,j) .ge. 0.) then     
      !-----------------------------------------------------------------
      ! freezing conditions 
      !-----------------------------------------------------------------
            Fbot = c0

          else

      !-----------------------------------------------------------------
      ! melting conditions
      !-----------------------------------------------------------------

           if (aice(i,j).gt.c0) then
            ! use boundary layer theory for Fbot
            deltaT = max((sst(i,j)-Tf(i,j)),c0) 

            ! strocnx has units N/m^2 so strocnx/rho has units m^2/s^2 
            ustar = sqrt(sqrt(strocnxT(i,j)**2+strocnyT(i,j)**2)/rhow)
            ustar = max(ustar,ustar_min)

            Fbot = cpchr*deltaT*ustar  ! < 0
            Fbot = max(Fbot,frzmlt(i,j))    ! frzmlt < Fbot < 0

!!! use all frzmlt for standalone runs
            Fbot = min(c0,frzmlt(i,j))    ! frzmlt = Fbot < 0

           else  ! aice = 0
            Fbot = c0
           endif ! aice
          endif ! sgn(frzmlt)

      end subroutine init_frzmlt

c=======================================================================

      subroutine add_new_ice(i,j,Fnew,hin,hsn)

!---!-------------------------------------------------------------------
!---! Given the volume of new ice grown in open water, compute its area
!---! and thickness and add it to the appropriate category or categories.
!---! Usually all the new ice is added to category 1.  An exception is
!---! made if there is no open water, in which case the new ice is 
!---! distributed evenly over the entire cell.
!---! New ice has hsnow = 0.
!---! Ice thicknesses may lie outside category bounds after new ice
!---! added; routine rebin must be called afterwards.
!---!
!---! authors Elizabeth C. Hunke
!---!         William H. Lipscomb
!---!-------------------------------------------------------------------

      use ice_kinds_mod

      integer (kind=int_kind), intent(in) :: i,j

      real (kind=dbl_kind), intent(out) ::
     &   Fnew             ! heat flx to open water for new ice (W/m^2)

      real (kind=dbl_kind), intent(inout) ::
     &   hin(ncat)        ! ice thickness for each cat        (m)
     &,  hsn(ncat)        ! snow thickness for each cat        (m)

      integer (kind=int_kind) ::
     &   nc               ! category indices
     &,  k                ! ice layer index

      real (kind=dbl_kind) ::
     &   vi0new           ! volume of new ice
     &,  ai0new           ! fractional area of new ice
     &,  hi0new           ! thickness of new ice
     &,  vsurp            ! volume of new ice added to cats > 1
     &,  hsurp            ! thickness of new ice added to cats > 1
     &,  qi0(nilyr)       ! frazil ice enthalpy
     &,  qi0av            ! mean value of qi0 for new ice (J kg-1)
     &,  vlyr             ! ice layer volume
     &,  area1            ! starting fractional area of existing ice

      real (kind=dbl_kind), parameter ::
     &   hfrazilmin = 0.05_dbl_kind  ! min thickness of new frazil ice (m)

      !-----------------------------------------------------------------
      ! total new ice volume, area, thickness
      !-----------------------------------------------------------------
      Fnew = max(frzmlt(i,j),c0)
      if (Fnew.gt.puny) then 
         rnilyr = real(nilyr)

         qi0av = c0
         do k = 1, nilyr
 
            ! ocn assumes new ice is fresh
            qi0(k) = -rhoi*Lfresh ! note sign convention
            ! something like this would be better
c             qi0(k) = -rhoi * (cp_ice*(Tmlt(k)-Tf(i,j))
c     &          + Lfresh*(1.-Tmlt(k)/Tf(i,j)) - cp_ocn*Tmlt(k))

            qi0av = qi0av + qi0(k)
         enddo
         qi0av = qi0av/rnilyr

         vi0new = -Fnew*dt / qi0av           ! note sign convention (qi)
      else
         vi0new = c0
      endif
      frazil(i,j) = vi0new                   ! for history diagnostics

      if (vi0new .gt. puny) then

            ! NOTE pop assumes fresh water flux due to freezing is
            ! NOT included in flux Fresh
ccc         Fresh(i,j) = Fresh(i,j) - vi0new/rhoi  ! for coupling

         ! new ice area and thickness
         ! hin_max(0) < new ice thickness < hin_max(2)
         if (ai0 .gt. puny) then
            hi0new = max(vi0new/ai0, hfrazilmin)
            if (hi0new.gt.hin_max(2) .and. ai0+puny.lt.c1) then
            ! distribute excess volume over all ice categories
               hi0new = hin_max(2)
               ai0new = ai0
               vsurp = vi0new - ai0new*hi0new
               vi0new = ai0new*hi0new
            else  ! put ice in a single category, with vsurp = 0
               ai0new = vi0new/hi0new
               vsurp = c0
            endif
             
      !-----------------------------------------------------------------
      ! combine new ice grown in open water with category 1 ice
      !-----------------------------------------------------------------

            area1 = ain(1)      ! save
            ain(1) = ain(1) + ai0new
            ai0 = ai0 - ai0new
            vin(1) = vin(1) + vi0new
      
            if (ain(1) .lt. puny) then
               print*, ''
               print*, 'istep1,my_task,i,j:', istep1,my_task,i,j
               print*, 'ain(1) =', ain(1)
               print*, 'ai0new =', ai0new
               stop
            endif

            hin(1) = vin(1) / ain(1)
            hsn(1) = hsn(1) * area1 / ain(1)
            vlyr = vi0new/rnilyr
            do k = 1, nilyr
               ein(k,1) = ein(k,1) + qi0(k)*vlyr
            enddo
            Tsfn(1) = (Tf(i,j)*ai0new + Tsfn(1)*area1) / ain(1)
            if (Tsfn(1) .gt. c0) Tsfn(1) = c0

         else                   ! ai0 < puny
            vsurp = vi0new

         endif                  ! ai0 > puny
         
      !-----------------------------------------------------------------
      ! distribute excess ice volume among ice categories by increasing
      ! ice thickness, leaving ice area unchanged
      !-----------------------------------------------------------------

         if (vsurp .gt. puny) then
            hsurp = vsurp / aice(i,j)  ! new thickness added to each cat
            do nc = 1, ncat
               if (ain(nc).gt.puny) then
                  hin(nc) = hin(nc) + hsurp  ! new cat thickness
                  vin(nc) = hin(nc)*ain(nc)  ! new cat volume
                  vlyr = hsurp/rnilyr*ain(nc) ! volume added to each layer
                  do k=1,nilyr
                     ein(k,nc) = ein(k,nc) + qi0(k)*vlyr ! new energy
                  enddo                
               endif
            enddo                
         endif                  ! vsurp > puny

      endif                     ! vi0new > puny

      end subroutine add_new_ice

c=======================================================================

      subroutine reduce_area(i,j,dhi1,hin,hsn)

!---!-------------------------------------------------------------------
!---! Reduce area when ice melts for special case of ncat=1
!---!
!---! Use CSM 1.0-like method of reducing ice area
!---! when melting occurs: assume only half the ice volume
!---! change goes to thickness decrease, the other half
!---! to reduction in ice fraction
!---!-------------------------------------------------------------------

      use ice_kinds_mod
      use ice_constants
      use ice_flux
      use ice_state

      integer (kind=int_kind), intent(in) :: i,j

      real (kind=dbl_kind), intent(in) ::
     &  dhi1         ! melt at bottom and top plus sublimation  (m)

      real (kind=dbl_kind), intent(inout) ::
     &   hin(ncat)   ! ice thickness for each cat        (m)
     &,  hsn(ncat)   ! ice thickness for each cat        (m)

      integer (kind=int_kind) :: nc

      real (kind=dbl_kind) ::
     &   hi0         ! current hi for ice fraction adjustment
     &,  dai0        ! change in aice for ice fraction adjustment

      !-----------------------------------------------------------------
      ! make sure thickness of cat 1 is at least hin_max(0)
      !-----------------------------------------------------------------
          if ((hin(1).le.hin_max(0)).and.(hin_max(0).gt.0.)) then
           ain(1) = vin(1)/hin_max(0)
           hin(1) = hin_max(0)
          endif

            if( ain(1) .gt. 0. ) then 
              if( dhi1 .lt. 0.0 ) then  
                hi0        = vin(1)/ain(1)
                dai0       = (vin(1)/(hi0-(dhi1/c2))) - ain(1)
                ain(1) = ain(1) + dai0 
              endif
            endif

      end subroutine reduce_area

c=======================================================================

      end module ice_therm_driver

c=======================================================================



