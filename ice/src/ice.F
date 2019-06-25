c $Id: ice.F,v 2.0 1999/12/09 22:25:35 eclare Exp $

      program ice
c-----------------------------------------------------------------------
c.. main driving routine for CICE in stand-alone mode
c
c.. author Elizabeth C. Hunke
c..        Fluid Dynamics Group, Los Alamos National Laboratory
c
c.. Copyright, 1999.  The Regents of the University of California.
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
      include "ice.h"          ! ice code common blocks
      include "netcdf.inc"     ! include netcdf
      include "timers.h"       ! ice common blocks for timers
      integer rec              ! record number for netCDF 

      call setup_mpi           ! initial setup for message passing
c.. timing devices
      call timer_clear(-1)     ! initialize all timers
      call timer_start(0)      ! begin timing entire run

c.. initialize
      call input_data          ! namelist variables
      call grid                ! grid variables
      call init_par            ! define constants and parameters
      call init_evp            ! define evp dynamics parameters, variables
      call flux_init           ! initialize coupler fluxes

      if (kread.eq.1) then
        call restartfile       ! start from restart data
        call calendar(time_forc)
      else
        call init_var          ! start from scratch
      endif

      call ice_fraction        ! total ice area
      call sfc_albedo          ! albedo based on initial ice distribution
      call getflux             ! read data from files

      call calendar(time)      ! determine the initial date
      write_ic = .true.        ! write initial conditions
      call iceavg         
      write_ic = .false.

c.. timestep loop ......................................................
      call timer_start(1)
      do 100 istep = 1,npt

      istep1 = istep1 + 1
      time = time + dt         ! determine the time and date
      call calendar(time)      !   at the end of the timestep
      if (update_forcing) call getflux   ! read data from files

      call transport           ! advection, thermodynamics, and
                               ! category conversions

      call roughness           ! compute transfer coefficients + some fluxes
      call sfc_albedo          ! albedo based on latest ice distribution

      call evp                 ! elastic-viscous-plastic dynamics 

c.......... write data ...........
      if (mod(istep,diagfreq).eq.0) call runtime_diags  ! standard out
      call iceavg
      if (write_restart) call dumpfile              ! dumps for restarting

 100  continue                 ! end of timestep loop
      call timer_stop(1)

      call timer_stop(0)
      call timer_print(-1)     ! print all timing information

      call end_run        

      stop
      end
