c=======================================================================
c SCCS module : driver.h   r3.2
c SCCS get    : 30 May 1996 / 16:41:25
c SCCS edit   : 30 May 1996 / 16:41:23
c SCCS path   : /net/crestone.u1/nieman/csm/drv-3.n/SCCS/s.driver.h
c=======================================================================
c Note: if integers impliment logicals, then 0 <=> false
c=======================================================================
c edited by E. C. Hunke 9 April 1998

      integer     tmpunit        ! unit # reserved for temporary usage
      parameter  (tmpunit = 90 )

      !----- restart file selection ------------------------------------
      character*32 rest_type   ! run type: initial, cont, regen, branch

      !----- stopping date selection -----------------------------------
      character*16 stop_option ! specify when to stop models
      integer      stop_date   ! date on which to stop the models
      integer      stop_nday   ! number of days to advance the models

      !----- history file selection ------------------------------------
      integer      hist_tavg   ! T => request time-avg history data

      !----- dbugging & timing info selection --------------------------
      integer      info_dbug   ! dbug level: > 0 => print extra info
      integer      info_time   ! T => print extra timing info

      !----- selection special initialization procedures ---------------
      integer     init_davg    ! T => do daily avg flux initialization

      !----- info needed by some msg passing libraries -----------------
      character*8  msg_group   ! msg passing group name
      character*16 msg_method  ! msg passing method choice
      integer      msg_maxt    ! msg passing dt for time out

      !----- rest, stop, hist, diag control flags ----------------------
      integer      stop_now    ! T => stop model now
      integer      stop_endd   ! T => stop model at end of day
      integer      rest_endd   ! T => create restart data at end of day
      integer      hist_endd   ! T => create history data at end of day

      !----- timestep & model date info --------------------------------
      integer      nstep       ! elapsed days, 0 <=> 1 Jan, year 0
      integer      nstep0      ! value of nstep on startup
      integer      nstop       ! end run when nstep = nstop
      integer      year        ! model year
      integer      month       ! model month
      integer      day         ! model day
      integer      sec         ! model elapsed seconds during model day
      integer      cdate       ! model year+month+date, coded: yyyymmdd
      integer      nadv_a      ! number of atm timesteps per day
      integer      nadv_i      ! number of ice timesteps per day
      integer      ice_freq    ! number of ice steps between ice comms

      !----- global error code -----------------------------------------
      integer      ierr        ! 0 <=> no errors or warnings

      !----- information buffer passed to & from models ----------------
      integer    nibuff        ! size of info-buffers
      parameter (nibuff=20)
      integer  ibuff_d(nibuff) ! info-buffer from drv
      integer  ibuff_i(nibuff) ! info-buffer from ice model

      !----- global data via common block ------------------------------

      COMMON / drv_global /
     &   rest_type   , 
     &   stop_option , stop_nday  , stop_date  , 
     &   hist_tavg   , 
     &   info_dbug   , info_time  , 
     &   init_davg   , 
     &   msg_group   , msg_method , msg_maxt   ,
     &   stop_now    , stop_endd  , rest_endd  ,  hist_endd ,
     &   nstep       , nstep0     , nstop      ,
     &   year        , month      , day        , sec        , 
     &   cdate       , nadv_a     , nadv_i     ,
     &   ice_freq    , ierr       ,
     &   ibuff_d     , ibuff_i 

      SAVE   / drv_global /

