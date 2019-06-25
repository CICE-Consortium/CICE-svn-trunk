c $Id: $
c=======================================================================
!---! Linear remapping scheme for the ice thickness distribution
!---!
!---! authors William H. Lipscomb
!---!         Elizabeth C. Hunke
!---!
!---! See Lipscomb, W. H.  Remapping the thickness distribution of sea 
!---!     ice.  J. Geophys. Res., 2001.
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

      module ice_itd_linear

      use ice_model_size
      use ice_kinds_mod
      use ice_constants
      use ice_state
      use ice_itd
      use ice_calendar
      use ice_fileunits

      implicit none

c=======================================================================

      contains

c=======================================================================

      subroutine linear_itd(i,j,dhin,hin_old,hin,hsn)

!---!-------------------------------------------------------------------
!---! Ice thickness distribution scheme that shifts ice among categories.
!---!
!---! The default scheme is linear remapping, which works as follows.  See
!---! Lipscomb (2001) for more details. 
!---!
!---!   Using the thermodynamic "velocities", interpolate to find the 
!---!   velocities in thickness space at the category boundaries, and 
!---!   compute the new locations of the boundaries.  Then for each 
!---!   category, compute the thickness distribution function,  g(h), 
!---!   between hL and hR, the left and right boundaries of the category.
!---!   Assume g(h) is a linear polynomial that satisfies two conditions:
!---!   (1) The ice area implied by g(h) equals ain(n).
!---!   (2) The ice volume implied by g(h) equals ain(n)*hice(n).
!---!   Once g(h) is computed, compute the ice area and volume lying 
!---!   between the initial and new boundaries, and transfer this area 
!---!   and volume to the neighboring cell, thus restoring the initial 
!---!   boundary.
!---!-------------------------------------------------------------------

      use ice_domain

      integer (kind=int_kind), intent(in) :: 
     &  i,j               ! grid indices

      real (kind=dbl_kind), intent(in) ::
     &   dhin(ncat)       ! thickness change for remapping (m)
     &,  hin_old(ncat)    ! starting value of hin      (m)

      real (kind=dbl_kind), intent(inout) ::
     &   hin(ncat)        ! ice thickness for each cat        (m)
     &,  hsn(ncat)        ! snow thickness for each cat       (m)

      integer (kind=int_kind) :: 
     &   n, k
     &, donor(ncat-1)     ! donor category index

      real (kind=dbl_kind) ::
     &   Hbnew(0:ncat)    ! new boundary locations
     &,  slope            ! rate of change of dhice with hice
     &,  dh0              ! change in ice thickness at h = 0
     &,  da0              ! area removed at h = 0
     &,  damax            ! max allowed reduction in category 1 area
     &,  g0, g1           ! linear coefficients for g(h) 
     &,  hL               ! left end of range over which g(h) > 0
     &,  hR               ! right end of range over which g(h) > 0
     &,  etamin, etamax   ! left and right limits of integration 
     &,  work1, work2     ! temporary variables
     &,  x1               ! etamax - etamin
     &,  x2               ! (etamax^2 - etamin^2) / 2
     &,  x3               ! (etamax^3 - etamin^3) / 3
     &,  daice(ncat-1)    ! ice area transferred across boundary
     &,  dvice(ncat-1)    ! ice volume transferred across boundary

      real (kind=dbl_kind) ::
     &   vice1, vice2     ! total ice volume, initial and final
     &,  vsnow1, vsnow2   ! total snow volume, initial and final
     &,  eice1, eice2     ! ice vol * qice, initial and final
     &,  esnow1, esnow2   ! snow vol * qsnow, initial and final
     &,  totfrac          ! total fractional area; should = 1  
     &,  icevol           ! ice volume
     &,  error            ! volume/energy conservation error

      call aggregate_pt(i,j)  ! update cell values first

      if (aice(i,j).gt.eps11) then

      !-----------------------------------------------------------------
      ! quantities that the subroutine should conserve
      !-----------------------------------------------------------------
      vice1 = c0
      vsnow1 = c0
      eice1 = c0
      esnow1 = c0
      do n = 1, ncat
         vice1 = vice1 + vin(n)
         vsnow1 = vsnow1 + vsn(n)
         do k = 1, nilyr
            eice1 = eice1 + ein(k,n)
         enddo
         esnow1 = esnow1 + esn(n)
      enddo

      !-----------------------------------------------------------------
      ! new category boundaries, Hbnew
      !-----------------------------------------------------------------
       Hbnew(0) = c0
       do n = 1, ncat-1
        if (hin_old(n).gt.eps11 .and. hin_old(n+1).gt.eps11) then
               ! interpolate between adjacent category growth rates
               slope = (dhin(n+1)-dhin(n)) / 
     &                 (hin_old(n+1)-hin_old(n))
               Hbnew(n) = hin_max(n) + dhin(n) + 
     &              slope * (hin_max(n) - hin_old(n))
        elseif (hin_old(n) .gt. eps11) then   ! hin_old(n+1)=0
               Hbnew(n) = hin_max(n) + dhin(n)
        elseif (hin_old(n+1) .gt. eps11) then ! hin_old(n)=0
               Hbnew(n) = hin_max(n) + dhin(n+1)
        else
               Hbnew(n) = hin_max(n)
        endif

      !-----------------------------------------------------------------
      ! check that each boundary lies between adjacent values of hice
      !-----------------------------------------------------------------

        if (ain(n).gt.eps11 .and. hin(n).gt.Hbnew(n)) then
           write(nu_diag,*) 'istep1 = ',istep1
           write(nu_diag,*) my_task,':',i,j,'ITD hin(n) > Hbnew(n)'
           write(nu_diag,*) 'cat ',n
           write(nu_diag,*) my_task,':',i,j,'hin(n) =', hin(n)
           write(nu_diag,*) my_task,':',i,j,'Hbnew(n) =', Hbnew(n)
           stop
        elseif (ain(n+1).gt.eps11 .and. 
     &          hin(n+1).lt.Hbnew(n)) then
           write(nu_diag,*) 'istep1 = ',istep1
           write(nu_diag,*) my_task,':',i,j,'ITD hin(n+1) < Hbnew(n)'
           write(nu_diag,*) 'cat ',n
           write(nu_diag,*) my_task,':',i,j,'hin(n+1) =', hin(n+1)
           write(nu_diag,*) my_task,':',i,j,'Hbnew(n) =', Hbnew(n)
           stop
        endif

       enddo  ! categories

      !-----------------------------------------------------------------
      ! Initialize hin_max(ncat) and compute Hbnew(ncat)
      !-----------------------------------------------------------------

         if (ain(ncat) .gt. eps11) then
            Hbnew(ncat) = c3*hin(ncat) - c2*Hbnew(ncat-1)
         else
            Hbnew(ncat) = 999.9_dbl_kind  ! arbitrary big number
         endif
         if (Hbnew(ncat) .lt. hin_max(ncat-1)) 
     &        Hbnew(ncat) = hin_max(ncat-1)
         hin_max(ncat) = Hbnew(ncat)

      !-----------------------------------------------------------------
      ! area lost due to melting of thin (category 1) ice
      !-----------------------------------------------------------------
       if (ain(1) .gt. eps11) then
               dh0 = dhin(1)        

            if (dh0.lt.c0) then   ! remove area from category 1   
               dh0 = min(-dh0,hin_max(1)) ! note: use |dh0| 
               call fit_line (1, hin_old, hin_max, g0, g1, hL, hR, i,j)
               etamax = min(dh0,hR) - hL ! upper limit of integration
                                         ! note: lower limit = 0
               if (etamax .gt. c0) then
                  x1 = etamax
                  x2 = p5 * etamax*etamax
                  da0 = g1*x2 + g0*x1 ! ice area removed
                                ! constrain new thickness <= hin_old
                  damax = ain(1) * (c1-hin(1)/hin_old(1)) ! damax > 0
                  da0 = min (da0, damax)
                                ! remove area, conserving volume
                  hin(1) = hin(1) * ain(1) / (ain(1)-da0)
                  ain(1) = ain(1) - da0
               endif            ! etamax > 0
            else  ! dh0 >= 0
               Hbnew(0) = min(dh0,hin_max(1))     ! shift boundary at h = 0
            endif
       endif ! ain(1) > eps11

      !-----------------------------------------------------------------
      ! Shift ice between categories as necessary  
      !-----------------------------------------------------------------
      do n = 1, ncat-1

         if (Hbnew(n) .gt. hin_max(n)) then  ! transfer from n to n+1

       ! Make sure Hbnew(n-1) < Hbound(n), and Hbnew(n) < Hbound(n+1).
       ! In principle we could allow this, but it would make the code 
       ! more complicated.
       if (Hbnew(n-1) .gt. hin_max(n)) then
       write(nu_diag,*) 'istep1 = ',istep1
       write(nu_diag,*) my_task,':',i,j,'ITD Hbnew(n-1) > hin_max(n)'
       write(nu_diag,*) 'cat ',n
       write(nu_diag,*) my_task,':',i,j,'Hbnew(n-1) =', Hbnew(n-1)
       write(nu_diag,*) my_task,':',i,j,'hin_max(n) =', hin_max(n)
       stop
       endif
       if (Hbnew(n) .gt. hin_max(n+1)) then
       write(nu_diag,*) 'istep1 = ',istep1
       write(nu_diag,*) my_task,':',i,j,'ITD Hbnew(n) > hin_max(n+1)'
       write(nu_diag,*) 'cat ',n
       write(nu_diag,*) my_task,':',i,j,'Hbnew(n) =', Hbnew(n)
       write(nu_diag,*) my_task,':',i,j,'hin_max(n+1) =', hin_max(n+1)
       stop
       endif

            ! compute linear coefficients for g(h) in cell n
            ! actually g(eta), where eta = h - HL
            ! compute area and volume between hin_max(n) and Hbnew(n)

            call fit_line (n, hin, Hbnew, g0, g1, hL, hR, i, j)
            etamin = max (hin_max(n), hL) - hL  ! left integration limit
            etamax = min (Hbnew(n), hR) - hL   ! right integration limit
            if (etamax .gt. etamin) then
               donor(n) = n
               x1 = etamax - etamin
               work1 = etamin*etamin
               work2 = etamax*etamax
               x2 = p5 * (work2 - work1)
               work1 = work1*etamin
               work2 = work2*etamax
               x3 = p333 * (work2 - work1)
               daice(n) = g1*x2 + g0*x1               ! area transferred
               dvice(n) = g1*x3 + g0*x2 + daice(n)*hL ! volume transferred
            else                ! no ice is transferred
               donor(n) = 0
               daice(n) = c0
               dvice(n) = c0
            endif

         elseif (Hbnew(n) .lt. hin_max(n)) then ! transfer from n+1 to n

       ! check for nearest neighbor condition violation
       if (Hbnew(n) .lt. hin_max(n-1)) then
       write(nu_diag,*) 'istep1 = ',istep1
       write(nu_diag,*) my_task,':',i,j,'ITD Hbnew(n) < hin_max(n-1)'
       write(nu_diag,*) 'cat ',n
       write(nu_diag,*) my_task,':',i,j,'Hbnew(n) =', Hbnew(n)
       write(nu_diag,*) my_task,':',i,j,'hin_max(n-1) =', hin_max(n-1)
       stop
       endif
       if (Hbnew(n+1) .lt. hin_max(n)) then
       write(nu_diag,*) 'istep1 = ',istep1
       write(nu_diag,*) my_task,':',i,j,'ITD Hbnew(n+1) < hin_max(n)'
       write(nu_diag,*) 'cat ',n
       write(nu_diag,*) my_task,':',i,j,'Hbnew(n+1) =', Hbnew(n+1)
       write(nu_diag,*) my_task,':',i,j,'hin_max(n) =', hin_max(n)
       stop
       endif

            ! compute linear coefficients for g(h) in cell n+1
            ! actually g(eta), where eta = h - HL
            ! compute area and volume between hin_max(n) and Hbnew(n)

            call fit_line (n+1, hin, Hbnew, g0, g1, hL, hR, i, j)
            etamax = min (hin_max(n), hR) - hL  ! right integration limit
                                               ! note: lower limit = 0
            if (etamax .gt. 0.) then           
               donor(n) = n+1
               x1 = etamax
               x2 = p5 * etamax*etamax 
               x3 = p666 * x2 * etamax   ! = etamax^3 / 3
               daice(n) = g1*x2 + g0*x1               ! area transferred
               dvice(n) = g1*x3 + g0*x2 + daice(n)*hL ! volume transferred 
            else
               donor(n) = 0
               daice(n) = c0
               dvice(n) = c0
            endif
            
         else  ! Hbnew(n) = hin_max(n); no ice is transferred
            
            donor(n) = 0
            daice(n) = c0
            dvice(n) = c0

         endif ! Hbnew(n) > or < hin_max(n)

      enddo                     ! loop over boundaries

      ! move ice
      call shift_ice (donor, daice, dvice, hin, hsn, i,j)

      !-----------------------------------------------------------------
      ! Update thicknesses
      !-----------------------------------------------------------------
      if (hi_min.gt.c0 .and. 
     &    ain(1).gt.eps11 .and. hin(1).lt.hi_min) then
         ain(1) = ain(1) * hin(1)/hi_min
         hin(1) = hi_min
      endif

      !-----------------------------------------------------------------
      ! Reset hin_max(ncat)
      !-----------------------------------------------------------------
      if (ain(ncat) .gt. eps11) then
         hin_max(ncat) = c3*hin(ncat) - c2*hin_max(ncat-1)
      else
         hin_max(ncat) = 999.9_dbl_kind     ! arbitrary big number
      endif

      call aggregate_pt(i,j)
      ai0 = c1 - aice(i,j)

      !-----------------------------------------------------------------
      ! check that conservative quantities have not changed
      !-----------------------------------------------------------------
      vice2 = c0
      vsnow2 = c0
      eice2 = c0
      esnow2 = c0
      do n = 1, ncat
         vice2 = vice2 + vin(n)
         vsnow2 = vsnow2 + vsn(n)
         do k = 1, nilyr
            eice2 = eice2 + ein(k,n)
         enddo
         esnow2 = esnow2 + esn(n)
      enddo

      error = abs((vice2-vice1)/vice1)
      if (error .gt. eps11) then
         write(nu_diag,*) 'istep1 = ',istep1
         write(nu_diag,*) my_task,':',i,j,'Ice volume error, ITD'
         write(nu_diag,*) my_task,':',i,j,'vice1 =', vice1, 
     &                                    'vice2 =', vice2
         write(nu_diag,*) my_task,':',i,j,'Final area,thickness,volume'
         do n = 1, ncat
            write(nu_diag,*) my_task,':',i,j,n, ain(n),hin(n), 
     &                                     ain(n)*hin(n)
         enddo
      endif

      if (vsnow1 .gt. eps11) then
         error = abs((vsnow2-vsnow1)/vsnow1)
      else
         error = c0
      endif
      if (error .gt. eps11) then
         write(nu_diag,*) 'istep1 = ',istep1
         write(nu_diag,*) my_task,':',i,j,'Snow volume error, ITD'
         write(nu_diag,*) my_task,':',i,j,'vsnow1 =', vsnow1,
     &                               'vsnow2 =', vsnow2
      endif

      error = abs((eice2-eice1)/eice1)
      if (error .gt. eps11) then
         write(nu_diag,*) 'istep1 = ',istep1
         write(nu_diag,*) my_task,':',i,j,'Ice energy error, ITD'
         write(nu_diag,*) my_task,':',i,j,'eice1 =', eice1, 
     &                                    'eice2 =', eice2
      endif

      if (abs(esnow1).gt.puny) then
        error = abs((esnow2-esnow1)/esnow1)
      else
        error = abs(esnow2)
      endif
      if (error .gt. eps11) then
         write(nu_diag,*) 'istep1 = ',istep1
         write(nu_diag,*) my_task,':',i,j,'Snow energy error, ITD'
         write(nu_diag,*) my_task,':',i,j,'esnow1 =', esnow1, 
     &                               'esnow2 =', esnow2
      endif


      endif  ! aice > eps11

      end subroutine linear_itd

c=======================================================================

      subroutine fit_line (n, hin, Hb, g0, g1, hL, hR, i,j)

!---!-------------------------------------------------------------------
!---! Fit g(h) with a line, satisfying area and volume constraints.
!---! To reduce rounding errors caused by large values of g0 and g1,
!---! we compute g(eta), where eta = h - hL.
!---!-------------------------------------------------------------------

      integer (kind=int_kind), intent(in) ::
     &   i,j
     &,  n               ! thickness category for which g(h) computed 

      real (kind=dbl_kind), intent(in) ::
     &   hin(ncat)       ! ice thickness for each cat        (m)
     &,  Hb(0:ncat)      ! current category boundaries

      real (kind=dbl_kind), intent(out) ::
     &   g0, g1          ! coefficients in linear equation for g(eta)
     &,  hL              ! min value of range over which g(h) > 0
     &,  hR              ! max value of range over which g(h) > 0 

      real  (kind=dbl_kind) ::
     &   h13             ! Hb(n-1) + 1/3 * (Hb(n) - Hb(n-1))
     &,  h23             ! Hb(n-1) + 2/3 * (Hb(n) - Hb(n-1))
     &,  etaice          ! hice(n) - Hb(n-1); ice thickness in eta space
     &,  icevol          ! ice volume in eta space
     &,  dh              ! hR - hL
     &,  dhi             ! 1 / dh
     &,  dh2, dh3        ! dh^2/2, dh^3/3
     &,  work1, work2    ! temporary variables
     &,  cmp_ain         ! computed ice area; should equal ain(n)
     &,  cmp_icevol      ! computed ice volume; should equal icevol
     &,  error           ! area/volume conservation error

      if (ain(n) .gt. eps11) then

      !-----------------------------------------------------------------
      ! Initialize hL and hR
      !-----------------------------------------------------------------
      hL = Hb(n-1)
      hR = Hb(n)

      ! Change hL or hR if hin(n) falls outside central third of range
      h13 = p333 * (c2*hL + hR)     ! p333 = 1/3
      h23 = p333 * (hL + c2*hR)
      if (hin(n) .lt. h13) then
         hR = c3*hin(n) - c2*hL
      elseif (hin(n) .gt. h23) then
         hL = c3*hin(n) - c2*hR
      endif

      !-----------------------------------------------------------------
      ! Compute ice thickness and volume in eta space
      !-----------------------------------------------------------------
      etaice = hin(n) - hL
      icevol = ain(n) * etaice
     
      !-----------------------------------------------------------------
      ! Compute linear coefficients of g(eta)
      !-----------------------------------------------------------------
      dh = hR - hL
      dhi = c1 / dh
      work1 = c6 * ain(n) * dhi
      work2 = etaice * dhi
      g0 = work1 * (p666 - work2)   ! p666 = 2/3 
      g1 = c2*dhi * work1 * (work2 - p5)

      !-----------------------------------------------------------------
      ! Make sure the area and volume implied by the computed g(h) 
      ! agree with the known area and volume
      !-----------------------------------------------------------------
      dh2 = p5*dh*dh
      dh3 = p666*dh2*dh  ! = dh^3 / 3
      cmp_ain = g1*dh2 + g0*dh
      cmp_icevol  = g1*dh3 + g0*dh2

      error = abs((cmp_ain-ain(n))/ain(n))
      if (error .gt. eps11) then
       write(nu_diag,*) my_task,':',i,j,'ITD Area csv error =',error
       write(nu_diag,*) 'cat ',n
       write(nu_diag,*) my_task,':',i,j,'ain(n) =', ain(n)
       write(nu_diag,*) my_task,':',i,j,'Computed area =', cmp_ain
       stop
      endif

      error = abs((cmp_icevol-icevol)/icevol)
      if (error .gt. eps11) then
       write(nu_diag,*) my_task,':',i,j,'ITD Vol csv error =',error
       write(nu_diag,*) 'cat ',n
       write(nu_diag,*) my_task,':',i,j,'icevol =', icevol
       write(nu_diag,*) my_task,':',i,j,'Computed volume =', cmp_icevol
       stop
      endif

      else  ! ain(n) < eps11
         g0 = c0
         g1 = c0
         hL = c0
         hR = c0
      endif

      end subroutine fit_line

c=======================================================================

      end module ice_itd_linear

c=======================================================================
