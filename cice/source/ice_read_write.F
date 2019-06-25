c $Id: $
c=======================================================================
!---! routines for opening, reading and writing external files
!---!
!---! author Tony Craig, NCAR
c=======================================================================

      module ice_read_write

      use ice_model_size
      use ice_domain 
      use ice_mpi_internal

      implicit none

c=======================================================================

      contains

c=======================================================================

      subroutine ice_open(nu,filename,nbits)

!---!-------------------------------------------------------------------
!---! opens an unformatted file for reading
!---! nbits indicates whether the file is sequential or direct access
!---!-------------------------------------------------------------------

      integer (kind=int_kind), intent(in) :: nu,nbits
      character (*) :: filename

      if (my_task.eq.master_task) then
      if (nbits.eq.0) then
        open(nu,file=filename,form='unformatted')
      else
#ifdef OSF1
        open(nu,file=filename,recl=imt_global*jmt_global*nbits/32,
     &       form='unformatted',access='direct')
#else
        open(nu,file=filename,recl=imt_global*jmt_global*nbits/8,
     &       form='unformatted',access='direct')
#endif
      endif
      endif

      end subroutine ice_open

c=======================================================================

      subroutine ice_read(nu,nrec,work,atype,scatter,diag)

!---!-------------------------------------------------------------------
!---! reads an unformatted file
!---! work is a real array, atype indicates the format of the data
!---!-------------------------------------------------------------------

      integer (kind=int_kind), intent(in) :: nu,nrec
      integer (kind=int_kind) :: n,i,j
      character (len=4) atype     ! real,int , 
      logical (kind=log_kind) :: scatter, diag
      real (selected_real_kind(12)) :: workg (imt_global,jmt_global) ! 8 byte
      real (selected_real_kind(6)) ::  workg4(imt_global,jmt_global) ! 4 byte
      integer (selected_int_kind(6)) :: workgi(imt_global,jmt_global)! 4 byte
      real (kind=dbl_kind) ::
     &   work(ilo:ihi,jlo:jhi)
     &,  amin,amax

      if (my_task.eq.master_task) then
        if (atype.eq.'ida4') then
          read(nu,rec=nrec) workgi
          workg=real(workgi)
        elseif (atype.eq.'rda8') then
          read(nu,rec=nrec) workg
        elseif (atype.eq.'rda4') then
          read(nu,rec=nrec) workg4
          workg = workg4
        elseif (atype.eq.'ruf8') then
          read(nu) ((workg(i,j),i=1,imt_global),j=1,jmt_global)
        else
          write (6,*) ' ERROR: reading unknown atype ',atype
        endif
        if (diag) then
          amin=workg(1,1)
          amax=workg(1,1)
          do j=1,jmt_global
          do i=1,imt_global
            amin=min(amin,workg(i,j))
            amax=max(amax,workg(i,j))
          enddo
          enddo
          write (6,*) ' read_global ',nu,nrec,amin,amax
        endif
      endif

      if (scatter) call global_scatter(workg,work)

      end subroutine ice_read

c=======================================================================

      subroutine ice_write(nu,nrec,work,atype,gather)

!---!-------------------------------------------------------------------
!---! writes an unformatted file
!---! work is a real array, atype indicates the format of the data
!---!-------------------------------------------------------------------

      integer (kind=int_kind), intent(in) :: nu,nrec
      integer (kind=int_kind) :: n,i,j
      character (len=4) atype     ! real,int , 
      logical (kind=log_kind) :: gather 
      real (selected_real_kind(12)) :: workg (imt_global,jmt_global) ! 8 byte
      integer (selected_int_kind(6)) :: workgi(imt_global,jmt_global)! 4 byte
      real (kind=dbl_kind) ::
     &   work(ilo:ihi,jlo:jhi)
     &,  amin,amax

      if (gather) call global_gather(workg,work)

      if (my_task.eq.master_task) then
        if (atype.eq.'ida4') then
          workgi=nint(workg)
          write(nu,rec=nrec) workgi
        elseif (atype.eq.'rda8') then
          write(nu,rec=nrec) workg
        elseif (atype.eq.'ruf8') then
          write(nu) ((workg(i,j),i=1,imt_global),j=1,jmt_global)
        else
          write (6,*) ' ERROR: writing unknown atype ',atype
        endif
        amin=workg(1,1)
        amax=workg(1,1)
        do j=1,jmt_global
        do i=1,imt_global
          amin=min(amin,workg(i,j))
          amax=max(amax,workg(i,j))
        enddo
        enddo
        write (6,*) ' writ_global ',nu,nrec,amin,amax
      endif

      end subroutine ice_write

c=======================================================================

      end module ice_read_write

c=======================================================================
