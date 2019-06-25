c $Id: ice_mp.F,v 1.5 1998/04/14 14:59:51 eclare Exp $

c     program ice_mp.F
c.. main driving routine for CICE, with message passing
c
c.. author Elizabeth C. Hunke
c..        Fluid Dynamics Group, Los Alamos National Laboratory
c
c.. Copyright, 1998.  The Regents of the University of California.
c.. This software was produced under a U.S. Government contract 
c.. (W-7405-ENG-36) by Los Alamos National Laboratory, which is operated
c.. by the University of California for the U.S. Department of Energy.  
c.. The U.S. Government is licensed to use, reproduce, and distribute this 
c.. software.  Permission is granted to the public to copy and use this 
c.. software without charge, provided that this Notice and any statement
c.. of authorship are reproduced on all copies.  Neither the Government
c.. nor the University makes any warranty, express or implied, or assumes
c.. any liability or responsibility for the use of this software.
c-----------------------------------------------------------------------
      implicit none
#include "ice.h"                                 ! ice code common blocks
#include "mcl.h"                                 ! MCL specific definitions
#include "ice_mp.h"                              ! common blocks for mp
#include "netcdf.inc"                            ! include netcdf
#include "ncdf.h"                                ! ice common blocks for netcdf
c#include "timers.h"                            

      integer ier                   ! error flag for message calls
      integer i, j, k               ! generic indices
      integer rec                   ! record number for netCDF 
      real umax, umax0, hmax, hmax0

      call input_data          ! namelist variables
      call grid                ! initialize grid 
      call init_mp             ! initialize message passing 
      call init_par            ! define constants and parameters

      if (kread.eq.1) call restartfile
      if (kread.eq.0) then     ! start from scratch
        call init_var          ! define pre-initial state
        call to_coupler        ! send pre-initial state 
        call from_coupler      ! get initial state from other components
        call init_more_var     ! complete cice variable initialization
        call to_coupler        ! send initial ice state 
      endif

c.. initialize netCDF and write initial state
      call netcdf_setup(rec)   ! masks and initial state: c, h's, u
      call netcdf_fluxes_setup(rec)   ! initial state: temps, fluxes

c.. timestep loop ......................................................

      do 100 istep = 1,npt

      call calendar            ! determine the date and time
      call from_coupler        ! get updated info from flux coupler

c.. compute thermodynamics and transport every dt seconds
      call transport           ! includes thermodynamics

      call to_coupler          ! collect/send data to flux coupler

c      goto 999                ! sans dynamics
c.. elastic-viscous-plastic dynamics 
      call mass_prss           ! mass, pressure
      call viscosity           ! eta, zeta
      call stressprep          ! coefficients for stress calculations

         do k=1,ndte           ! subcycling
      call stress              ! stress tensor equation (sigma)
      call stepu               ! total surface stress, momentum equation
      call boundu(u)           ! Periodic/Neumann boundary conditions
      call boundu(v)           ! Periodic/Neumann boundary conditions
c      call div_u              ! velocity divergence
         enddo
 999  continue

c.. maximum velocity component
      umax=0.
c$DOACROSS LOCAL(i,j,umax0),REDUCTION(umax)
      do j=1,jmt
      do i=1,imt
       umax0=umax
       umax = max(umax0,abs(u(i,j)),abs(v(i,j)))
      enddo
      enddo
c.. maximum thickness
      hmax=0.
c$DOACROSS LOCAL(i,j,hmax0),REDUCTION(hmax)
      do j=1,jmt
      do i=1,imt
       hmax0=hmax
       hmax = max(hmax0,compact(i,j)*hithick(i,j))
      enddo
      enddo

c.......... write data ...........
      if (mod(istep,iplot).eq.0) then
        print*, 'max cH = ',hmax
        print*, 'max u or v = ',umax
c.. netCDF
        rec = rec + 1          ! increment record number
        call writerecs(rec)
        call writerecs_fluxes(rec)
c.. Fortran unformatted
        call dumpfile          ! dumps for restarting
      endif

 100  continue

      print*, 'dt = ',dt,'  dte = ',dte

c.. disconnect the model from the coupled system  
        ibuff(2) = 1           ! stop now flag
        call from_coupler
       
        if (ibuff_d(2) .ne. 1) then
          ibuff(2) = 1         ! stop now flag
          call to_coupler
        endif
	call mcl_finalize(ier)

      stop
      end

