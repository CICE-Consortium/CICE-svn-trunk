c=======================================================================
c SCCS module : calendar.f   r3.1
c SCCS get    : 03 May 1996 / 17:24:45
c SCCS edit   : 03 May 1996 / 17:23:46
c SCCS path   : /net/crestone.u1/nieman/csm/drv-3.n/SCCS/s.calendar.f
c=======================================================================
c  Note: 
c   these calendar routines to conversions between
c   o the integer number of elapsed days, 
c   o the integer triple (year,month,day)
c   o the integer coded calendar date: yyyymmdd
c
c  Assumptions:
c   o there is a year 0
c   o all years have 365 days (no leap years)
c   o elapsed days = 0 <=> start of  1 Jan, year 0
c
c=======================================================================
c edited by E. C. Hunke 9 April 1998

      SUBROUTINE eday2date(eday,date)

      implicit none

      integer  eday,date

      integer  calendar(12),k,year,month,day
      data     calendar /0,31,59,90,120,151,181,212,243,273,304,334/
      save     calendar

c-----------------------------------------------------------------------
c  PURPOSE:
c    compute the calendar date: year/month/day
c  INPUT:
c    an integer number of elapsed days
c  OUTPUT:
c    coded (yyyymmdd) calendar date
c  NOTE:
c    this calendar has a year zero (but no day or month zero)
c-----------------------------------------------------------------------

      year = eday/365       ! calandar year (note: Fortran truncation)
      day  = mod(eday,365)  ! elapsed days within current year
      DO k=1,12
        IF (day .ge. calendar(k)) month=k   ! calendar month
      END DO
      day = day-calendar(month) + 1         ! calendar day

      date = year*10000 + month*100 + day  ! coded calendar date

      RETURN
      END

c=======================================================================
c edited by E. C. Hunke 9 April 1998

      SUBROUTINE eday2ymd (eday,year,month,day)

      implicit none

      integer  eday,year,month,day

      integer  calendar(12),k
      data     calendar /0,31,59,90,120,151,181,212,243,273,304,334/
      save     calendar

c-----------------------------------------------------------------------
c  PURPOSE:
c    compute the calendar date: year/month/day
c  INPUT:
c    an integer number of elapsed days
c  OUTPUT:
c    uncoded calendar date, integer year, month, & day
c  NOTE:
c    this calendar has a year zero (but no day or month zero)
c-----------------------------------------------------------------------

      year = eday/365       ! calandar year (note: Fortran truncation)
      day  = mod(eday,365)  ! elapsed days within current year
      DO k=1,12
        IF (day .ge. calendar(k)) month=k   ! calendar month
      END DO
      day = day-calendar(month) + 1         ! calendar day

      RETURN
      END

c=======================================================================
c edited by E. C. Hunke 9 April 1998

      SUBROUTINE date2ymd (date,year,month,day)

      implicit none

      integer  date,year,month,day

c-----------------------------------------------------------------------
c  PURPOSE:
c    decode the calendar date
c  INPUT:
c    calendar date in (integer) yyyymmdd format
c  OUTPUT:
c    calendar year,month,day
c-----------------------------------------------------------------------

      year =int(     date       /10000)
      month=int( mod(date,10000)/  100)
      day  =     mod(date,  100) 

      RETURN
      END
c=======================================================================
c edited by E. C. Hunke 9 April 1998

      SUBROUTINE date2eday(date,eday)

      implicit none

      integer  date,eday

      integer  calendar(12),year,month,day
      data     calendar /0,31,59,90,120,151,181,212,243,273,304,334/
      save     calendar

c-----------------------------------------------------------------------
c  PURPOSE:
c    derive elapsed days from the calendar date
c  INPUT:
c    calendar date in (integer) yyyymmdd format
c  OUTPUT:
c    elapsed days since yy-mm-dd = 00-01-01, with 0 elapsed seconds
c-----------------------------------------------------------------------

      year =int(     date       /10000)
      month=int( mod(date,10000)/  100)
      day  =     mod(date,  100) 

      eday = year*365 + calendar(month) + (day-1)

      RETURN
      END

c=======================================================================
c edited by E. C. Hunke 9 April 1998

      SUBROUTINE  ymd2eday(year,month,day,eday)

      implicit none

      integer  year,month,day,eday

      integer  calendar(12)
      data     calendar /0,31,59,90,120,151,181,212,243,273,304,334/
      save     calendar

c-----------------------------------------------------------------------
c  PURPOSE:
c    derive elapsed days from the calendar date
c  INPUT:
c    calendara year, month, & date
c  OUTPUT:
c    elapsed days since yy-mm-dd = 00-01-01, with 0 elapsed seconds
c-----------------------------------------------------------------------

      eday = year*365 + calendar(month) + (day-1)

      RETURN
      END

c=======================================================================
c edited by E. C. Hunke 9 April 1998

      logical FUNCTION datevalid(date)

      implicit none

      integer  date

      integer  year,month,day,dpm(12)
      data     dpm /31,28,31,  30,31,30,  31,31,30,  31,30,31/
      save     dpm ! J  F  M    A  M  J    J  A  S    O  N  D

c-----------------------------------------------------------------------
c  PURPOSE:
c    determine if a valid year, month & day can be decoded 
c    from the coded calendar date
c  INPUT:
c    calendar date in (integer) yyyymmdd format
c  RETURNS:
c    true or false
c-----------------------------------------------------------------------

      year =int(     date       /10000)
      month=int( mod(date,10000)/  100)
      day  =     mod(date,  100) 

      datevalid = .true.
      if (year  .lt.0) datevalid = .false.
      if (month.lt. 1) datevalid = .false.
      if (month.gt.12) datevalid = .false.
      if (day  .lt. 1) datevalid = .false.
      if (.not. datevalid ) return
      if (day .gt. dpm(month)) datevalid = .false.

      RETURN
      END

c=======================================================================
