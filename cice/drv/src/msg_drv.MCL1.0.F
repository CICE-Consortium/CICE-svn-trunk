c=======================================================================
c SCCS module : msg_drv.MCL1.0.f   r3.2
c SCCS get    : 04 Jun 1996 / 11:26:53
c SCCS edit   : 04 Jun 1996 / 11:26:47
c SCCS path   : /net/crestone.u1/nieman/csm/drv-3.n/SCCS/s.msg_drv.MCL1.0.f
c=======================================================================
c  edited by E. C. Hunke, 8 April 1998

      SUBROUTINE msg_drvc ()

      implicit none

      !--- global ---
      include 'driver.h'       ! contains msg_group & msg_maxt
      include 'dims.h'         ! defines domain dimensions
      include 'msg.MCL1.0.h'   ! msg passing stuff

      !--- local ---
      integer      info        ! MCL lib error return code
      character*8  group_name  ! PVM group name
      integer      maxsec      ! maximum sec to wait for model check in
      integer      method      ! selects msg passing method
      integer      total_task  ! total number of tasks in system
      integer      geom        ! unused argument to mcl_alloc() 

      integer      nbytes_icec ! number of bytes in ice cntl buffer
      integer      nbytes_icei ! number of bytes in ice->drv buffer
      integer      nbytes_ices ! number of bytes in ice->drv buffer
      integer      nbytes_icer ! number of bytes in ice<-drv buffer

      parameter   (nbytes_icec = MCL_SIZEOF_INT *nibuff)
      parameter   (nbytes_icei = MCL_SIZEOF_REAL*nx_i*ny_i*3)
      parameter   (nbytes_ices = MCL_SIZEOF_REAL*nx_i*ny_i*15)
      parameter   (nbytes_icer = MCL_SIZEOF_REAL*nx_i*ny_i*17)

c-----------------------------------------------------------------------
c PURPOSE:
c   o initialize NCAR MCL 1.0 message passing lib
c-----------------------------------------------------------------------
      
      total_task = 2   ! drv + one component model makes two tasks
      group_name = msg_group
      maxsec     = 60
      method     = MCL_DATACH_PVM      
      if (msg_maxt   .gt.  maxsec     ) maxsec=msg_maxt
      if (msg_method .eq. "named-pipe") method=MCL_DATACH_PIPE
      if (msg_method .eq. "disk-drop" ) method=MCL_DATACH_DISK

c      call pvmfsettaskname ("drv",info) ! model name for XPVM runs

      !-----------------------------------------------------------------
      ! initalize connection to the component models
      !-----------------------------------------------------------------
      call mcl_initcoupler(group_name,total_task,maxsec,info)
      if (info.ge.MCL_SUCCESS) then
         write(*,9000) '(msg_drvc) success calling mcl_initcoupler'
      else
         write(*,9000) '(msg_drvc) ERROR: calling mcl_initcoupler'
         call mcl_perror(info,'(msg_drvc) after mcl_initcoupler: ')
      end if


      call mcl_con2model(model_num_ice,tid_ice,method,info)
      if (info.ge.MCL_SUCCESS) then
         write(*,9000) '(msg_drvc) success calling mcl_con2model'
      else
         write(*,9000) '(msg_drvc) ERROR: calling mcl_con2model'
         call mcl_perror(info,'(msg_drvc) after mcl_con2model: ')
      end if

      !-----------------------------------------------------------------
      ! allocate memory needed & used by MCL message passing library
      ! # bytes = (grid resolutution)*(# fields)*(size of data element)
      !-----------------------------------------------------------------

      call mcl_alloc(dd_icei,nbytes_icei,tid_ice,geom,info)
      call mcl_alloc(dd_ices,nbytes_ices,tid_ice,geom,info)
      call mcl_alloc(dd_icer,nbytes_icer,tid_ice,geom,info)

      dd_null(1)=0 
      dd_icec(1)=0 

      if (info.ge.MCL_SUCCESS) then
        write(*,9000) '(msg_drvc) success calling mcl_alloc'
      else
        write(*,9000) '(msg_drvc) ERROR: calling mcl_alloc'
        call mcl_perror(info,'(msg_drvc) after mcl_alloc:')
      end if

      !-----------------------------------------------------------------
      ! deal with an unsuccessful connection to the msg passing system
      !-----------------------------------------------------------------
      if (info.lt.MCL_SUCCESS) then
         call mcl_finalize(info)
         STOP
      end if

      RETURN
 9000 format(a)
      END

c=======================================================================
c edited by E. C. Hunke 9 April 1998

      SUBROUTINE msg_drvd ()

      implicit none

      integer  info

c-----------------------------------------------------------------------
c PURPOSE:
c   dis-connect drv to msg passing environment (MCL 1.0 msg passing lib)
c-----------------------------------------------------------------------
      
      write(*,9000) '(msg_drvd) dis-connecting from msg-passing.'
      call mcl_finalize(info)

      RETURN
 9000 format(a)
      END

c=======================================================================
