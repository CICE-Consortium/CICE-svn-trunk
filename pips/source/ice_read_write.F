! CVS: $Id: ice_read_write.F,v 1.6 2004/02/25 17:35:14 eclare Exp $
! CVS: $Source: /home/climate/CVS-COSIM/cice/source/ice_read_write.F,v $
! CVS: $Name:  $
!=======================================================================
!BOP
!
! !MODULE: ice_read_write
!
! !DESCRIPTION:
!
! Routines for opening, reading and writing external files
!
! !REVISION HISTORY:
!
! author: Tony Craig, NCAR
!
! !INTERFACE:
!
      module ice_read_write
!
! !USES:
!
      use ice_model_size
      use ice_domain 
      use ice_mpi_internal
      use ice_fileunits
      use ice_work, only:  work_g1, work_gr
!
!EOP

      implicit none

!=======================================================================

      contains

!=======================================================================
!
!BOP
!
! !IROUTINE: ice_open - opens an unformatted file for reading
!
! !INTERFACE:
!
      subroutine ice_open(nu, filename, nbits)
!
! !DESCRIPTION:
!
! Opens an unformatted file for reading \\
! nbits indicates whether the file is sequential or direct access
!
! !REVISION HISTORY:
!
! author: Tony Craig, NCAR
!
! !USES:
!
! !INPUT/OUTPUT PARAMETERS:
!
      integer (kind=int_kind), intent(in) :: 
     &     nu            ! unit number
     &,    nbits         ! no. of bits per variable (0 for sequential access)

      character (*) :: filename
!
!EOP
!
      if (my_task.eq.master_task) then

         if (nbits.eq.0) then   ! sequential access

            open(nu,file=filename,form='unformatted')

         else                   ! direct access
! tcmod
!#ifdef OSF1
!           open(nu,file=filename,recl=imt_global*jmt_global*nbits/32,
!     &          form='unformatted',access='direct')
!#else
!           open(nu,file=filename,recl=imt_global*jmt_global*nbits/8,
!     &          form='unformatted',access='direct')
!#endif
             open(nu,file=filename,recl=imt_global*jmt_global*nbits/8,
     &            form='unformatted',access='direct')
         endif                   ! nbits = 0 

      endif                      ! my_task = master_task

      end subroutine ice_open

!=======================================================================
!BOP
!
! !IROUTINE: ice_read - reads an unformatted file
!
! !INTERFACE:
!
      subroutine ice_read(nu,  nrec,  work, atype, scatter, diag,
     &                    ignore_eof, hit_eof)
!
! !DESCRIPTION:
!
! Reads an unformatted file \\
! work is a real array, atype indicates the format of the data
!
! !REVISION HISTORY:
!
! author: Tony Craig, NCAR
!
! !USES:
!
! !INPUT/OUTPUT PARAMETERS:
!
      integer (kind=int_kind), intent(in) :: 
     &     nu                ! unit number
     &,    nrec              ! record number (0 for sequential access)

      real (kind=dbl_kind), intent(out) ::  
     &     work(ilo:ihi,jlo:jhi) ! output array (real, 8-byte)

      character (len=4) atype    ! format for input array
                                 ! (real/integer, 4-byte/8-byte)

      logical (kind=log_kind) :: 
     &     scatter               ! if true, scatter the data
     &,    diag                  ! if true, write diagnostic output

      logical (kind=log_kind), optional, intent(in)  :: ignore_eof
      logical (kind=log_kind), optional, intent(out) :: hit_eof
!
!EOP
!
      integer (kind=int_kind) :: i, j, ios

      integer(selected_int_kind(12)), dimension(:,:), allocatable :: 
     &   workgi
      integer(selected_int_kind(6)), dimension(:,:), allocatable :: 
     &   workgi4

      real (kind=dbl_kind) :: 
     &   amin, amax         ! min and max values of input array

      logical (kind=log_kind) :: ignore_eof_use

      if (my_task.eq.master_task) then

    !-------------------------------------------------------------------
    ! Read global array according to format atype
    !-------------------------------------------------------------------
        if (present(hit_eof)) hit_eof = .false.

        if (atype.eq.'ida8') then
           allocate (workgi(imt_global,jmt_global))
           read(nu,rec=nrec) workgi
           work_g1 = real(workgi)
           deallocate (workgi)
        elseif (atype.eq.'ida4') then
           allocate (workgi4(imt_global,jmt_global))
           read(nu,rec=nrec) workgi4
           work_g1 = real(workgi4)
           deallocate (workgi4)
        elseif (atype.eq.'rda8') then
           read(nu,rec=nrec) work_g1
        elseif (atype.eq.'rda4') then
           allocate (work_gr(imt_global,jmt_global))
           read(nu,rec=nrec) work_gr
           work_g1 = work_gr
           deallocate (work_gr)
        elseif (atype.eq.'ruf8') then
           if (present(ignore_eof)) then
             ignore_eof_use = ignore_eof
           else
             ignore_eof_use = .false.
           endif
           if (ignore_eof_use) then
             ! Read line from file, checking for end-of-file
             read(nu, iostat=ios) ((work_g1(i,j),i=1,imt_global),
     &                                         j=1,jmt_global)
             if (present(hit_eof)) hit_eof = ios < 0
           else
             read(nu) ((work_g1(i,j),i=1,imt_global),j=1,jmt_global)
           endif
        else
           write(nu_diag,*) ' ERROR: reading unknown atype ',atype
        endif
      endif

      if (present(hit_eof)) then
        call ice_bcast_logical(hit_eof)
        if (hit_eof) return
      endif
 
    !-------------------------------------------------------------------
    ! optional diagnostics
    !-------------------------------------------------------------------
      if (diag) then
        if (my_task.eq.master_task) then
          amin = work_g1(1,1)
          amax = work_g1(1,1)
          do j = 1,jmt_global
          do i = 1,imt_global
            amin = min(amin,work_g1(i,j))
            amax = max(amax,work_g1(i,j))
          enddo
          enddo
          write(nu_diag,*) ' read_global ',nu, nrec, amin, amax
        endif
      endif

    !-------------------------------------------------------------------
    ! Scatter data to individual processors
    !-------------------------------------------------------------------

      if (scatter) call global_scatter(work_g1,work)

      end subroutine ice_read

!=======================================================================
!BOP
!
! !IROUTINE: ice_write - writes an unformatted file
!
! !INTERFACE:
!
      subroutine ice_write(nu, nrec, work, atype, gather)
!
! !DESCRIPTION:
!
! Writes an unformatted file \\
! work is a real array, atype indicates the format of the data
!
! !REVISION HISTORY:
!
! author: Tony Craig, NCAR
!
! !USES:
!
! !INPUT/OUTPUT PARAMETERS:
!
      integer (kind=int_kind), intent(in) :: 
     &     nu                ! unit number
     &,    nrec              ! record number (0 for sequential access)

      real (kind=dbl_kind), intent(in) ::  
     &     work(ilo:ihi,jlo:jhi) ! input array (real, 8-byte)

      character (len=4) atype    ! format for output array
                                 ! (real/integer, 4-byte/8-byte)

      logical (kind=log_kind) :: gather   ! if true, gather the data
!
!EOP
!
      integer (kind=int_kind) :: i, j

      integer(selected_int_kind(12)), dimension(:,:), allocatable :: 
     &   workgi
      integer(selected_int_kind(6)), dimension(:,:), allocatable :: 
     &   workgi4

      real (kind=dbl_kind) :: 
     &   amin, amax     ! min and max values of ouput array

    !-------------------------------------------------------------------
    ! Gather data from individual processors
    !-------------------------------------------------------------------
      if (gather) call global_gather(work_g1,work)

      if (my_task.eq.master_task) then

    !-------------------------------------------------------------------
    ! Write global array according to format atype
    !-------------------------------------------------------------------
        if (atype.eq.'ida8') then
           allocate (workgi(imt_global,jmt_global))
           workgi = nint(work_g1)
           write(nu,rec=nrec) workgi
           deallocate (workgi)
        elseif (atype.eq.'ida4') then
           allocate (workgi4(imt_global,jmt_global))
           workgi4 = nint(work_g1)
           write(nu,rec=nrec) workgi4
           deallocate (workgi4)
        elseif (atype.eq.'rda8') then
           write(nu,rec=nrec) work_g1
        elseif (atype.eq.'rda4') then
           allocate (work_gr(imt_global,jmt_global))
           work_gr = work_g1
           write(nu,rec=nrec) work_gr
           deallocate (work_gr)
        elseif (atype.eq.'ruf8') then
           write(nu) ((work_g1(i,j),i=1,imt_global),j=1,jmt_global)
        else
           write(nu_diag,*) ' ERROR: writing unknown atype ',atype
        endif

    !-------------------------------------------------------------------
    ! diagnostics
    !-------------------------------------------------------------------
        amin = work_g1(1,1)
        amax = work_g1(1,1)
        do j = 1,jmt_global
        do i = 1,imt_global
           amin = min(amin,work_g1(i,j))
           amax = max(amax,work_g1(i,j))
        enddo
        enddo
        write(nu_diag,*) ' writ_global ', nu, nrec, amin, amax

      endif                     ! my_task = master_task

      end subroutine ice_write
 
!=======================================================================

      end module ice_read_write

!=======================================================================
