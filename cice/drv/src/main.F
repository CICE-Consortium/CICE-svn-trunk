c=======================================================================
c SCCS module : main.f   r3.3
c SCCS get    : 13 May 1996 / 16:13:55
c SCCS edit   : 13 May 1996 / 16:13:52
c SCCS path   : /net/crestone.u1/nieman/csm/drv-3.n/SCCS/s.main.f
c=======================================================================
c  edited by E. C. Hunke, 8 April 1998
c
c         NCAR CSM Flux Coupler is subject to Copyright (C) 1996
c            University Corporation for Atmospheric Research
c                         All Rights Reserved
c
c             +------------------------------------------+
c             | DISTRIBUTION TERMS AND CONDITIONS NOTICE |
c             +------------------------------------------+
c
c (c) Copyright 1996 University Corporation for Atmospheric Research/
c National Center for Atmospheric Research/
c Climate and Global Dynamics Division
c
c This software, the NCAR CSM Flux Coupler was written by the 
c Oceanography Section of the Climate and Global Dynamics Division (CGD)
c of the National Center for Atmospheric Research (NCAR), which is
c operated by the University Corporation for Atmospheric Research (UCAR)
c and sponsored by the National Science Foundation (NSF).
c
c Access and use of this software shall impose the following obligations
c and understandings on the user.  The user is granted the right,
c without any fee or cost, to use, copy, modify, alter, enhance and
c distribute this software, and any derivative works thereof, and its
c supporting documentation for any purpose whatsoever, except commercial
c sales, provided that this entire notice appears in all copies of the
c software, derivative works and supporting documentation.  Further, the
c user agrees to credit UCAR/NCAR/CGD in any publications that result
c from the use of this software or in any software package that includes
c this software.  The names UCAR/NCAR/CGD, however, may not be used in
c any advertising or publicity to endorse or promote any products or
c commercial entity unless specific written permission is obtained from
c UCAR/NCAR/CGD.
c
c THE NCOM MATERIALS ARE MADE AVAILABLE WITH THE UNDERSTANDING THAT
c UCAR/NCAR/CGD IS NOT OBLIGATED TO PROVIDE (AND WILL NOT PROVIDE) THE
c USER WITH ANY SUPPORT, CONSULTING, TRAINING, OR ASSISTANCE OF ANY KIND
c WITH REGARD TO THE USE, OPERATION AND PERFORMANCE OF THIS SOFTWARE, 
c NOR TO PROVIDE THE USER WITH ANY UPDATES, REVISIONS, NEW VERSIONS, OR
c "BUG FIXES."
c
c THIS SOFTWARE IS PROVIDED BY UCAR/NCAR/CGD "AS IS" AND ANY EXPRESS OR
c IMPLIED WARRANTIES, INCLUDING BUT NOT LIMITED TO, THE IMPLIED
c WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
c DISCLAIMED.  IN NO EVENT SHALL UCAR/NCAR/CGD BE LIABLE FOR ANY
c SPECIAL, INDIRECT OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES WHATSOEVER,
c INCLUDING BUT NOT LIMITED TO CLAIMS ASSOCIATED WITH THE LOSS OF DATA
c OR PROFITS, WHICH MAY RESULT FROM AN ACTION IN CONTRACT, NEGLIGENCE OR
c OTHER TORTIOUS CLAIM THAT ARISES OUT OF OR IN CONNECTION WITH THE
c ACCESS, USE OR PERFORMANCE OF THIS SOFTWARE.
c=======================================================================
c  edited by E. C. Hunke, 8 April 1998

      PROGRAM main 

      IMPLICIT none

      !----- global variables -----
      include 'dims.h'         ! defines computational domain dimensions
      include 'index.h'        ! pointers to buffer locations
      include 'buffer.h'       ! allocates all flux & state buffers 
      include 'driver.h'       ! misc global driver data

      !----- local variables -----
      integer      n      ! inner time-stepping loop index
      character*32 tstamp ! time-stamp function (real world date & time)
      character*32 str    ! output string from tstamp function

      integer i,j
c-----------------------------------------------------------------------
c PURPOSE:
c   This is the CSM coupled model driver which coordinates the 
c   flux calculations & time integration of all component models: 
c   o atm - the atmosphere model (with an embedded lnd model)
c   o ice - the global sea ice model
c   o ocn - the global ocn model
c-----------------------------------------------------------------------

      str=tstamp()
      write(*,9002)
      write(*,9001) 'Executing CSM Flux Coupler  --  ',str
      write(*,9002)
      write(*,9999) ! copywrite notice
      write(*,9003)

      !-----------------------------------------------------------------
      ! initialize driver: 
      ! get input parms, connect to msg-passing, get restart data
      !-----------------------------------------------------------------

      call initial ()
      call msg_drvc ()
      call restart ('read')
      call srh_logic ()   ! stop, rest & hist logic: set flags

      !-----------------------------------------------------------------
      ! exchange time invariant msgs with component models
      !-----------------------------------------------------------------

      call msg_icei ()

      !-----------------------------------------------------------------
      ! set frequency for model communications
      !-----------------------------------------------------------------
      nadv_a = 72                   ! hard-wired atm steps per day
      ice_freq = nadv_a/nadv_i

      !-----------------------------------------------------------------
      ! receive IC data (states & fluxes) from component models.
      !-----------------------------------------------------------------
        call msg_icer (Si_i ,Fioi_i, Faii_i)  ! receive from ice model
        call flux_ices(Sa_i ,So_i, Fixx) ! create ocn/atm fluxes to send

      !-----------------------------------------------------------------
      ! if it is an initial run, must also exchange extra messages
      ! with other components to initialize interface fluxes computed
      ! within those components
      !-----------------------------------------------------------------
      if (rest_type(1:4) .eq. 'init') then

        !--- send state variables (and shortwave flux) to ice ---
        call mrgf_ice ( Faia_i , Faid_i , Fioo_i , Fiod_i , Fixx)
        call msg_ices (Sa_i, So_i, Fixx)

        !--- receive ice outputs and map to other grids ---
        call msg_icer (Si_i ,Fioi_i, Faii_i)

      endif

      !-----------------------------------------------------------------
      !               check time coordination of IC data
      !-----------------------------------------------------------------
         !--- in this case, all model should be coordinated in time ---
         write(*,9004)
         write(*,9001) 'Start of time coordinated integration.'
         write(*,9005)
         call timecheck(1)  ! 1 => enforce time coordination

      !-----------------------------------------------------------------
      !                 start of model integration
      !-----------------------------------------------------------------

      DO WHILE ( stop_now .eq. 0 )
       
         DO n=1,nadv_a

            !--- merge ice inputs and keep time average ---
            call mrgf_ice ( Faia_i , Faid_i , Fioo_i , Fiod_i , Fixx)

            !--- send fields to ice if it is time ---
            if (mod(n,ice_freq) .eq. 0) then
              call msg_ices (Sa_i, So_i, Fixx)
            endif

            !--- receive ice outputs if it is time and map ---
            if (mod(n,ice_freq) .eq. 0) then
              call msg_icer (Si_i ,Fioi_i, Faii_i)
            endif

            !--- drv & atm have advanced part of sync-interval ---
            sec = n*86400/nadv_a
            call srh_logic() 
         END DO

         !--- drv & component models have advanced one day ---
         nstep = nstep + 1  ! advance year,month,day,seconds
         call eday2date(nstep,cdate)
         call eday2ymd (nstep,year,month,day)
         sec = 0            
         if (init_davg.ne.0) then
           write(*,9004)
           write(*,9001) 'End of initial t-avg ocn flux calculation'
           write(*,9001) '  Start of time coordinated integration'
           write(*,9005)
           init_davg=0
         end if
         call timecheck(1)  ! verify & enforce time coordination
         call srh_logic()   ! stop, rest, hist logic: set flags

      END DO

      !-----------------------------------------------------------------
      ! send final msg to models, dis-connect drv from msg-passing 
      !-----------------------------------------------------------------

      call msg_icef ()
      call msg_drvd ()

      call shutdown () !--- driver shutdown routine ---

      str=tstamp()
      write(*,9002)
      write(*,9001) 'Terminating CSM driver  --  ',str
      write(*,9003)

      STOP
 9000 format(a)
 9001 format(2a)
 9002 format(80('='))
 9003 format(80('=')/)
 9004 format(80('-'))
 9005 format(80('-')/)
 9999 format(
     $ ' NCAR CSM Flux Coupler [4 June 96],  Copyright (C) 1996'/,
     $ ' University Corporation for Atmospheric Research, ',
     $ ' All Rights Reserved')

      END

c=======================================================================
