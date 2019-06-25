c $Id: ice_calendar.F,v 1.2 2001/02/23 21:06:52 schramm Exp $
c=======================================================================
!---! calendar routines for managing time
!---!
!---! authors Elizabeth C. Hunke, LANL
!---!         Tony Craig, NCAR
c=======================================================================

      module ice_calendar

      use ice_constants
      implicit none

      real (kind=dbl_kind), parameter :: 
     &   secday = 86400.0_dbl_kind ! number of seconds in a day
      integer (kind=int_kind) ::
     &   daymo(12)                 ! number of days in each month
     &,  daycal(13)                ! day number at end of month
      data daymo/31,28,31,30,31,30,31,31,30,31,30,31/
      data daycal/0,31,59,90,120,151,181,212,243,273,304,334,365/

      integer (kind=int_kind) ::
     &   istep   ! local step counter for time loop
     &,  istep0  ! step counter, number of steps taken in previous run
     &,  istep1  ! step counter, number of steps at current timestep
     &,  mday    ! day of the month
     &,  week    ! week of the year
     &,  month   ! 1 \le month \le 12
     &,  monthp  ! last month
     &,  year    ! initial year
     &,  nyr     ! year number
     &,  idate   ! date (yyyymmdd)
     &,  sec     ! elapsed seconds into date
     &,  npt     ! total number of time steps (dt) 
     &,  stop_now           ! if 1, end program execution
     &,  write_restart      ! if 1, write restart now
     &,  cpl_write_history  ! if 1, write history on command from cpl

      real (kind=dbl_kind) :: 
     &   dt                 ! thermo/transport timestep (s)
     &,  dtei               ! 1/dte, where dte is the EVP timestep (1/s)
     &,  time               ! total elapsed time (s)
     &,  time_forc          ! time of last forcing update (s)
     &,  yday               ! day of the year

      logical (kind=log_kind) ::
     &   new_year           ! new year = .true.
     &,  new_month          ! new month = .true.
     &,  new_week           ! new week = .true.
     &,  new_day            ! new day = .true.
     &,  write_ic           ! write initial condition now
     &,  write_history      ! write history now
      character (len=1) ::
     &   histfreq           !  output frequency, 'y','m','d','1'
     &,  dumpfreq           ! restart frequency, 'y','m','d'
      integer (kind=int_kind) :: diagfreq  ! diagnostic output frequency
                                           !    (10 = once per 10 dt)

c=======================================================================

      contains

c=======================================================================

      subroutine init_calendar

!---!-------------------------------------------------------------------
!---! initialize calendar variables
!---!-------------------------------------------------------------------

      istep1 = istep0   ! number of steps at current timestep
                        ! real (dumped) or imagined (use to set calendar)
      istep = 0         ! local timestep number
      time=istep0*dt    ! s
      yday=c0           ! absolute day number
      mday=0            ! day of the month
      month=0           ! month
      nyr=0             ! year
      idate=00000101    ! date
      sec=0             ! seconds into date
      stop_now = 0      ! end program execution if stop_now=1

      end subroutine init_calendar

c=======================================================================

      subroutine calendar(ttime)

      use ice_fileunits

!---!-------------------------------------------------------------------
!---! determine the date at the end of the time step
!---!-------------------------------------------------------------------

      real (kind=dbl_kind), intent(in) ::
     &   ttime                          ! time variable
      integer (kind=int_kind) ::
     &   k
     &,  nyrp,mdayp,weekp               ! previous year, day, week
      real (kind=dbl_kind) ::
     &   tday                           ! absolute day number
     &,  dayyr                          ! number of days per year

      dayyr = 365.0_dbl_kind

      nyrp=nyr
      monthp=month
      weekp=week
      mdayp=mday
      new_year=.false.
      new_month=.false.
      new_week=.false.
      new_day=.false.
      write_history=.false.
      write_restart=0

      sec = mod(ttime,secday)           ! elapsed seconds into date at 
                                        ! end of dt
      tday = (ttime-sec)/secday + c1    ! absolute day number 
      yday = mod(tday-c1,dayyr) + c1    ! day of the year 
      week = int(yday/c7)               ! week of the year
      do k = 1, 12
        if (yday .gt. float(daycal(k))) month = k ! month
      enddo
      mday = int(yday) - daycal(month)       ! day of the month
      nyr = year + int((tday-c1)/dayyr)      ! year number - 1

      idate = (nyr)*10000 + month*100 + mday ! date (yyyymmdd) 

      if (my_task.eq.master_task.and.mod(istep,diagfreq).eq.0) then
        write(nu_diag,*) ' '
        write(nu_diag,*) 'istep1:',istep1,' date: ',idate,' sec: ',sec
      endif

      if (istep.ge.npt+1) stop_now = 1
      if (nyr.ne.nyrp) new_year = .true.
      if (month.ne.monthp) new_month = .true.
      if (week.ne.weekp) new_week = .true.
      if (mday.ne.mdayp) new_day = .true.

      if (istep.gt.1) then
        if (((histfreq.eq.'y'.or.histfreq.eq.'Y').and.new_year).or. 
     &      ((histfreq.eq.'m'.or.histfreq.eq.'M').and.new_month).or.
     &      ((histfreq.eq.'w'.or.histfreq.eq.'W').and.new_week).or.
     &      ((histfreq.eq.'d'.or.histfreq.eq.'D').and.new_day).or.
     &       (histfreq.eq.'1'))
     &      write_history=.true.
        if (((dumpfreq.eq.'y'.or.dumpfreq.eq.'Y').and.new_year).or. 
     &      ((dumpfreq.eq.'m'.or.dumpfreq.eq.'M').and.new_month).or.
     &      ((dumpfreq.eq.'d'.or.dumpfreq.eq.'D').and.new_day  ))
     &      write_restart=1
      endif

      end subroutine calendar

c=======================================================================

      end module ice_calendar

c=======================================================================
