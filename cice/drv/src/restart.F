c=======================================================================
c SCCS module : restart.f   r3.2
c SCCS get    : 13 May 1996 / 15:40:40
c SCCS edit   : 13 May 1996 / 15:40:36
c SCCS path   : /net/crestone.u1/nieman/csm/drv-3.n/SCCS/s.restart.f
c=======================================================================
c  edited by E. C. Hunke, 8 April 1998

      SUBROUTINE restart (cOption)

      implicit none

      !----- global variables -----
      include 'dims.h'              ! contains array dimensions
      include 'driver.h'            ! contains timing & data file info
      include 'buffer.h'            ! contains necessary states/fluxes

      !----- dummy arguments -----
      character*(*) cOption

      !----- local variables -----
      character*80 str         ! restart file ID string
      character*17 rest_outfile ! filename for output restart file
      character*80 loc_fn      ! local      file name
      integer      ier         ! msread/mswrite error code

c-----------------------------------------------------------------------
c PURPOSE:
c   IF (cOption='read')
c     Acquire drv restart/IC data (read a data file xor fabricate)
c     Restart data is the model output and daily avg fluxes
c     computed during the previous day of integration.
c
c  DETAILS 
c    wrt cOption='read':
c    IF (rest_type = 'initial' )
c       IC's are fabricated by this routine
c
c  MACHINE DEPENDENCY:
c       files read/written are binary unformatted and so are machine
c       dependent
c-----------------------------------------------------------------------
      
      !-----------------------------------------------------------------
      !  read a restart file
      !-----------------------------------------------------------------
      IF (cOption .eq. 'read') THEN

         !----- get IC file name -----
         IF      (rest_type(1:4) .eq. 'init') THEN
            write(*,*) '(restart) this is an initial run...'
            write(*,*) '(restart) IC data fabricated internally'
         ELSE
            write(*,*) '(restart) invalid rest_type = ',rest_type
            STOP 'restart'
         END IF

         !----- acquire IC data -----
         IF (rest_type(1:4) .eq. 'init') THEN
            write(*,*) '(restart) internally fabricating IC data.'
            call setarray( Si_i  , 0.0, nx_i,ny_i,nSi  )
            call setarray( Fioi_i, 0.0, nx_i,ny_i,nFioi)
            call setarray( Faii_i, 0.0, nx_i,ny_i,nFaii)
            str   = 'internally fabricated data'
            cdate = 000101
         END IF

         !----- reset IC data date?? -----
         call date2ymd (cdate,year,month,day)
         call date2eday(cdate,nstep)
         nstep0=nstep
         sec=0
         write(*,'(a,i8.8)') '(restart) using initial date: ',cdate

      !-----------------------------------------------------------------
      !  there aren't any other choices
      !-----------------------------------------------------------------
      ELSE
         write(*,*) '(restart) ERROR: invalid cOption = ',cOption
      END IF

      RETURN
      END

C=======================================================================
