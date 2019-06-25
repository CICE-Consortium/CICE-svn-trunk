c=======================================================================
c SCCS module : driver.f   r3.2
c SCCS get    : 04 Jun 1996 / 11:20:23
c SCCS edit   : 04 Jun 1996 / 11:20:20
c SCCS path   : /net/crestone.u1/nieman/csm/drv-3.n/SCCS/s.driver.f
c=======================================================================
c  edited by E. C. Hunke, 8 April 1998

      SUBROUTINE initial ()

      IMPLICIT none

      INCLUDE 'driver.h'
      INCLUDE 'dims.h'

      NAMELIST / drv_parm / 
     &   rest_type   ,   
     &   stop_option ,   
     &   stop_nday   ,   
     &   stop_date   ,
     &   hist_tavg   ,
     &   info_dbug   ,
     &   info_time   ,
     &   init_davg   ,
     &   msg_method  ,
     &   msg_group   ,
     &   msg_maxt    

      integer i, true  , false  , unset
      parameter (true=1, false=0, unset=-999999)

c-----------------------------------------------------------------------
c PURPOSE:
c   o set default input parameters
c   o read drv parameter (namelist) file
c-----------------------------------------------------------------------
      
      ierr   = 0  ! assume no errors yet
      nstop  = 0  ! initialized here for safety, will be reset

      !-----------------------------------------------------------------
      ! input parmarameter defaults
      !-----------------------------------------------------------------
      rest_type   = 'initial'
      stop_option = 'newyear'
      stop_nday   = 0
      stop_date   = 010101
      hist_tavg   = 1
      info_dbug   = 1
      info_time   = 0
      init_davg   =  unset
      msg_method  = 'unknown'
      msg_group   = 'unknown'
      msg_maxt    = 300

      !-----------------------------------------------------------------
      ! read input parms from stdin
      !-----------------------------------------------------------------
      write(*,'(a)') '(initial) reading parameters from input namelist.'
      read (*,drv_parm)

      !-----------------------------------------------------------------
      ! set defaults for initial run procedures
      !-----------------------------------------------------------------
      if ( init_davg .eq. unset) then
        if (rest_type(1:4) .eq. 'init') init_davg=true
        if (rest_type(1:4) .ne. 'init') init_davg=false
      end if

      !-----------------------------------------------------------------
      ! print values of all input parms
      !-----------------------------------------------------------------
      write(*,'(a)') '(initial) parameter values:'
      write(*,drv_parm)
      write(*,9000)

      RETURN
 9000 format(80('-'))
      END

c=======================================================================
c edited by E. C. Hunke 9 April 1998

      SUBROUTINE shutdown ()

      INCLUDE 'driver.h'

c-----------------------------------------------------------------------
c PURPOSE:
c   do all neccessary driver shut-down stuff
c-----------------------------------------------------------------------
      
      write(*,'(a,i8.8)') '(shutdown) final date: ',cdate

      RETURN
      END

c=======================================================================
c  edited by E. C. Hunke, 8 April 1998

      SUBROUTINE timecheck(iflag)

      INCLUDE 'driver.h'

      integer     iflag,nstep_x,nadv_x,date_x,year_x,month_x,day_x,sec_x
      character*3 cmonth_x

      character*3 char_month(12)
      data        char_month /'Jan','Feb','Mar','Apr','May','Jun',
     $                        'Jul','Aug','Sep','Oct','Nov','Dec'/
      save char_month

c-----------------------------------------------------------------------
c PURPOSE:
c   check for time coordination between component models
c
c OPTIONS:
c   iflag .ne. 0 <=> the drv will stop if component model dates 
c                    are uncoordinated (seconds are ignored)
c-----------------------------------------------------------------------
      
      !--- print out driver timestep & date ---
      nadv_x  = 1
      nstep_x = nstep
      date_x  = cdate
      sec_x   = sec
      call date2ymd(date_x,year_x,month_x,day_x)
      cmonth_x=char_month(month_x)
      write(*,9000) 'drv',nstep_x,nadv_x,year_x,cmonth_x,day_x,sec_x

      !--- check for time coordination (if requested) ---
      if (iflag.ne.0) then
        if (cdate .ne. ibuff_i(6)) ierr=-1
        if (sec   .ne. ibuff_i(7)) ierr=-1
      end if

      !--- print out model timesteps & dates ---
      if (    nstep.eq.nstep0 .or. nstep.eq.nstop   .or.
     &    info_dbug.gt.1      .or.  ierr.ne.0     ) then
         nadv_x  = ibuff_i(4)
         nstep_x = ibuff_i(5)
         date_x  = ibuff_i(6)
         sec_x   = ibuff_i(7)
         call date2ymd(date_x,year_x,month_x,day_x)
         cmonth_x=char_month(month_x)
         write(*,9000) 'ice',nstep_x,nadv_x,year_x,cmonth_x,day_x,sec_x
      endif

      !--- enforce time coordination ---
      if (ierr.eq.-1) then
         write(*,*)'(timecheck) ERROR: models uncoordinated in time'
         STOP 'model t'
      end if 

      RETURN
c     model=atm, step=1234567, steps/day=12, date=1984-Jan-01, 12345 sec
 9000 format('(timecheck) model=',a3,
     &       ', step='           ,i7,
     &       ', steps/day='      ,i2,
     &       ', date='           ,i4.4,'-',a3,'-',i2.2,', ',i5,' sec')
      END

c=======================================================================
c  edited by E. C. Hunke, 8 April 1998

      SUBROUTINE srh_logic()

      IMPLICIT none

      INCLUDE 'driver.h'

      integer    true,  false
      parameter (true=1,false=0)
      integer    year2,month2,day2
      integer    first_call
      logical   datevalid
      external  datevalid
      data      first_call / true /
      save      first_call

c-----------------------------------------------------------------------
c PURPOSE:
c   "srh_logic"  means  "stopping, restarting, & history data logic"
c   Here we set flags wrt stopping & restart & history data creation.
c   These flags are used by the driver, but are also sent to component 
c   models to facilitate the coordination of data sets
c
c Note:
c   o when integer values are used to implement boolean logic, 
c     zero <=> false, non-zero <=> true
c-----------------------------------------------------------------------

      !-----------------------------------------------------------------
      ! do one-time calculations
      !-----------------------------------------------------------------
      if ( first_call .eq. true ) then

        !-------------------------------------------------------
        ! compute nstop - value of nstep when integration stops
        !-------------------------------------------------------
        IF (stop_option(1:4).eq.'date') THEN
           !-- stop on given date --
           call date2eday(stop_date,nstop)
        ELSE IF (stop_option(1:5).eq.'ndays') THEN
           !-- stop after n days --
           nstop=nstep0+stop_nday
           call eday2date(nstop,stop_date)
        ELSE IF (stop_option(1:8).eq.'newmonth') THEN
           !-- stop on day 1 of next month -- 
           if (month.eq.12) call ymd2eday(year+1,      1,1,nstop)
           if (month.ne.12) call ymd2eday(year  ,month+1,1,nstop)
           call eday2date(nstop,stop_date)
        ELSE IF (stop_option(1:7).eq.'newyear') THEN
           !-- stop on January 1 of next year -- 
           call ymd2eday(year+1,1,1,nstop)
           call eday2date(nstop,stop_date)
        ELSE IF (stop_option(1:9).eq.'newdecade') THEN
           !-- stop on January 1 of next decade --
           call ymd2eday(((year/10)+1)*10,1,1,nstop)
           call eday2date(nstop,stop_date)
        ELSE 
           write(*,*) '(srh_logic) invalid stop_option = ',stop_option
           STOP 'srh_logic'
        END IF
        write(*,9001) nstop-nstep0,stop_date

        !-------------------------------------------------------
        ! warn user about the dbug of data they've requested
        !-------------------------------------------------------
        if      (info_dbug.eq.3) then
           write(*,9002) ' WARNING: info_dbug = ',info_dbug
           write(*,9002) ' * This implies hist & diag info every day'
        else if (info_dbug.gt.3) then
           write(*,9002) ' WARNING: info_dbug = ',info_dbug
           write(*,9002) ' * This implies hist & diag info every day'
           write(*,9002) ' * & every step, starting cdate = ',info_dbug
        end if

        !-------------------------------------------------------
        ! set default values for all flags
        !-------------------------------------------------------
        stop_now   = true         
        stop_endd  = false         
        rest_endd  = false         
        hist_endd  = false         
        first_call = false
      end if

      !-----------------------------------------------------------------
      ! stopping, restart & history data only happens when sec=0
      !-----------------------------------------------------------------
      if ( sec .ne. 0 ) then
        stop_now  = false         

      !-----------------------------------------------------------------
      ! it's a new day (sec=0), re-do all the calcuations
      !-----------------------------------------------------------------
      else
        stop_now  = false         
        stop_endd = false         
        rest_endd = false         
        hist_endd = false         

        !-------------------------------------------------------
        ! what is tomorrow's year, month, & day??
        !-------------------------------------------------------
        call eday2ymd (nstep+1,year2,month2,day2) 

        !-------------------------------------------------------
        ! stop data now or at the end of this day??
        !-------------------------------------------------------
        if (  nstep    .ge. nstop)  stop_now  = true         
        if ( (nstep+1) .ge. nstop)  stop_endd = true
      endif

      RETURN
   90 format(a)
 9001 format('(srh_logic) drv will advance ',i4,' days, '
     &       'stopping on ',i8.8)
 9002 format('(srh_logic) ',a,i8.8)
      END

c=======================================================================
c edited by E. C. Hunke 8 April 1998

      character*32 FUNCTION tstamp() 

      character (8) :: cTMP1, cTMP3
      character (10) :: cTMP2
      character (32) :: cDATE

c-----------------------------------------------------------------------
c  PURPOSE:
c    returns a string containing the real world time & date
c    
c  MACHINE DEPENDENCY:
c    These time & date functions are F90 standard
c-----------------------------------------------------------------------

      call date_and_time(cTMP1, cTMP2)
      cTMP3 = cTMP2(1:8)

      !--- create & write a time stamp ---
      write(cDATE,'(4a)') 'Date: ',cTMP1,'   Time: ',cTMP3
      tstamp = cDate

      RETURN
      END

c=======================================================================
c edited by E. C. Hunke 9 April 1998

      SUBROUTINE setarray(A,r,nx,ny,nz)

      IMPLICIT none

      INTEGER    nx,ny,nz
      REAL     A(nx,ny,nz),r

      INTEGER  i,j,k

c-----------------------------------------------------------------------
c PURPOSE:
c   set all elements of array A equal to r
c-----------------------------------------------------------------------

      DO k=1,nz
c$    doacross local(i,j)
      DO j=1,ny
      DO i=1,nx
        A(i,j,k)=r
      END DO
      END DO
      END DO

      RETURN
      END

c=======================================================================
