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
!  SVN:$Id$
!
! author: Tony Craig, NCAR
!
! 2004: Block structure added by William Lipscomb, LANL
!
! !INTERFACE:
!
      module ice_read_write
!
! !USES:
!
      use ice_kinds_mod
      use ice_constants
      use ice_communicate, only: my_task, master_task
      use ice_broadcast
      use ice_domain_size
      use ice_blocks
      use ice_fileunits
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
      if (my_task == master_task) then

         if (nbits == 0) then   ! sequential access

            open(nu,file=filename,form='unformatted')

         else                   ! direct access
            open(nu,file=filename,recl=nx_global*ny_global*nbits/8,
     &            form='unformatted',access='direct')
         endif                   ! nbits = 0

      endif                      ! my_task = master_task

      end subroutine ice_open

!=======================================================================
!BOP
!
! !IROUTINE: ice_read - read and scatter an unformatted file
!
! !INTERFACE:
!
      subroutine ice_read(nu,  nrec,  work, atype, diag,
     &                    ignore_eof, hit_eof)
!
! !DESCRIPTION:
!
! Read an unformatted file and scatter to processors\\
! work is a real array, atype indicates the format of the data
!
! !REVISION HISTORY:
!
! author: Tony Craig, NCAR
!
! !USES:
!
      use ice_domain
      use ice_gather_scatter
      use ice_work, only: work_g1, work_gr, work_gi4, work_gi8
!
! !INPUT/OUTPUT PARAMETERS:
!
      integer (kind=int_kind), intent(in) ::
     &     nu                ! unit number
     &,    nrec              ! record number (0 for sequential access)

      real (kind=dbl_kind), dimension(nx_block,ny_block,max_blocks),
     &     intent(out) ::
     &     work              ! output array (real, 8-byte)

      character (len=4) ::
     &     atype             ! format for input array
                             ! (real/integer, 4-byte/8-byte)

      logical (kind=log_kind) ::
     &     diag              ! if true, write diagnostic output

      logical (kind=log_kind), optional, intent(in)  :: ignore_eof
      logical (kind=log_kind), optional, intent(out) :: hit_eof
!
!EOP
!
      integer (kind=int_kind) :: i, j, ios

      real (kind=dbl_kind) ::
     &   amin, amax         ! min and max values of input array

      logical (kind=log_kind) :: ignore_eof_use

      if (my_task == master_task) then
         allocate(work_g1(nx_global,ny_global))
      else
         allocate(work_g1(1,1))   ! to save memory
      endif

      if (my_task == master_task) then

    !-------------------------------------------------------------------
    ! Read global array according to format atype
    !-------------------------------------------------------------------
         if (present(hit_eof)) hit_eof = .false.

         if (atype == 'ida4') then
            allocate(work_gi4(nx_global,ny_global))
            read(nu,rec=nrec) work_gi4
            work_g1 = real(work_gi4,kind=dbl_kind)
            deallocate(work_gi4)
         elseif (atype == 'ida8') then
            allocate(work_gi8(nx_global,ny_global))
            read(nu,rec=nrec) work_gi8
            work_g1 = real(work_gi8,kind=dbl_kind)
            deallocate(work_gi8)
         elseif (atype == 'rda4') then
            allocate(work_gr(nx_global,ny_global))
            read(nu,rec=nrec) work_gr
            work_g1 = work_gr
            deallocate(work_gr)
         elseif (atype == 'rda8') then
            read(nu,rec=nrec) work_g1
         elseif (atype == 'ruf8') then
            if (present(ignore_eof)) then
               ignore_eof_use = ignore_eof
            else
               ignore_eof_use = .false.
            endif
            if (ignore_eof_use) then
             ! Read line from file, checking for end-of-file
               read(nu, iostat=ios) ((work_g1(i,j),i=1,nx_global),
     &                                             j=1,ny_global)
               if (present(hit_eof)) hit_eof = ios < 0
            else
               read(nu) ((work_g1(i,j),i=1,nx_global),j=1,ny_global)
            endif
         else
            write(nu_diag,*) ' ERROR: reading unknown atype ',atype
         endif
      endif                     ! my_task = master_task

      if (present(hit_eof)) then
         call broadcast_scalar(hit_eof,master_task)
         if (hit_eof) then
            deallocate(work_g1)
            return
         endif
      endif


    !-------------------------------------------------------------------
    ! optional diagnostics
    !-------------------------------------------------------------------
      if (my_task==master_task .and. diag) then
         amin = minval(work_g1)
         amax = maxval(work_g1)
         write(nu_diag,*) ' read_global ',nu, nrec, amin, amax
      endif

    !-------------------------------------------------------------------
    ! Scatter data to individual processors.
    ! NOTE: Ghost cells are not updated.
    !-------------------------------------------------------------------

      call scatter_global(work, work_g1, master_task, distrb_info,
     &                    field_loc_noupdate, field_type_noupdate)

      deallocate(work_g1)

      end subroutine ice_read

!=======================================================================
!BOP
!
! !IROUTINE: ice_read_global - read an unformatted file
!
! !INTERFACE:
!
      subroutine ice_read_global (nu,  nrec,  work_g, atype, diag,
     &                            ignore_eof, hit_eof)
!
! !DESCRIPTION:
!
! Read an unformatted file \\
! Just like ice_read except that it returns a global array \\
! work_g is a real array, atype indicates the format of the data
!
! !REVISION HISTORY:
! Adapted by William Lipscomb, LANL, from ice_read
!
! !USES:
!
      use ice_work, only: work_gr, work_gi4, work_gi8
!
! !INPUT/OUTPUT PARAMETERS:
!
      integer (kind=int_kind), intent(in) ::
     &     nu                ! unit number
     &,    nrec              ! record number (0 for sequential access)

      real (kind=dbl_kind), dimension(nx_global,ny_global),
     &     intent(out) ::
     &     work_g            ! output array (real, 8-byte)

      character (len=4) ::
     &     atype             ! format for input array
                             ! (real/integer, 4-byte/8-byte)

      logical (kind=log_kind) ::
     &     diag              ! if true, write diagnostic output

      logical (kind=log_kind), optional, intent(in)  :: ignore_eof
      logical (kind=log_kind), optional, intent(out) :: hit_eof
!
!EOP
!
      integer (kind=int_kind) :: i, j, ios

      real (kind=dbl_kind) ::
     &   amin, amax         ! min and max values of input array

      logical (kind=log_kind) :: ignore_eof_use

      work_g(:,:) = c0

      if (my_task == master_task) then

    !-------------------------------------------------------------------
    ! Read global array according to format atype
    !-------------------------------------------------------------------
         if (present(hit_eof)) hit_eof = .false.

         if (atype == 'ida4') then
            allocate(work_gi4(nx_global,ny_global))
            read(nu,rec=nrec) work_gi4
            work_g = real(work_gi4,kind=dbl_kind)
            deallocate(work_gi4)
         elseif (atype == 'ida8') then
            allocate(work_gi8(nx_global,ny_global))
            read(nu,rec=nrec) work_gi8
            work_g = real(work_gi8,kind=dbl_kind)
            deallocate(work_gi8)
         elseif (atype == 'rda4') then
            allocate(work_gr(nx_global,ny_global))
            read(nu,rec=nrec) work_gr
            work_g = work_gr
            deallocate(work_gr)
         elseif (atype == 'rda8') then
            read(nu,rec=nrec) work_g
         elseif (atype == 'ruf8') then
            if (present(ignore_eof)) then
               ignore_eof_use = ignore_eof
            else
               ignore_eof_use = .false.
            endif
            if (ignore_eof_use) then
               ! Read line from file, checking for end-of-file
               read(nu, iostat=ios) ((work_g(i,j),i=1,nx_global),
     &                                            j=1,ny_global)
               if (present(hit_eof)) hit_eof = ios < 0
            else
               read(nu) ((work_g(i,j),i=1,nx_global),j=1,ny_global)
            endif
         else
            write(nu_diag,*) ' ERROR: reading unknown atype ',atype
         endif
      endif                     ! my_task = master_task

      if (present(hit_eof)) then
         call broadcast_scalar(hit_eof,master_task)
         if (hit_eof) return
      endif

    !-------------------------------------------------------------------
    ! optional diagnostics
    !-------------------------------------------------------------------
      if (my_task == master_task .and. diag) then
         amin = minval(work_g)
         amax = maxval(work_g)
         write(nu_diag,*) ' read_global ',nu, nrec, amin, amax
      endif

      end subroutine ice_read_global

!=======================================================================
!BOP
!
! !IROUTINE: ice_write - writes an unformatted file
!
! !INTERFACE:
!
      subroutine ice_write(nu, nrec, work, atype, diag)
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
      use ice_gather_scatter
      use ice_domain
      use ice_work, only: work_g1, work_gr, work_gi4, work_gi8
!
! !INPUT/OUTPUT PARAMETERS:
!
      integer (kind=int_kind), intent(in) ::
     &     nu                ! unit number
     &,    nrec              ! record number (0 for sequential access)

      real (kind=dbl_kind), dimension(nx_block,ny_block,max_blocks),
     &     intent(in) ::
     &     work              ! input array (real, 8-byte)

      character (len=4) ::
     &     atype             ! format for output array
                             ! (real/integer, 4-byte/8-byte)

      logical (kind=log_kind) ::
     &     diag              ! if true, write diagnostic output
!
!EOP
!
      integer (kind=int_kind) :: i, j

      real (kind=dbl_kind) ::
     &   amin, amax     ! min and max values of ouput array

    !-------------------------------------------------------------------
    ! Gather data from individual processors
    !-------------------------------------------------------------------

      if (my_task == master_task) then
         allocate(work_g1(nx_global,ny_global))
      else
         allocate(work_g1(1,1)) ! to save memory
      endif

      call gather_global(work_g1, work, master_task, distrb_info)

      if (my_task == master_task) then

    !-------------------------------------------------------------------
    ! Write global array according to format atype
    !-------------------------------------------------------------------
         if (atype == 'ida4') then
            allocate(work_gi4(nx_global,ny_global))
            work_gi4 = nint(work_g1)
            write(nu,rec=nrec) work_gi4
            deallocate(work_gi4)
         elseif (atype == 'ida8') then
            allocate(work_gi8(nx_global,ny_global))
            work_gi8 = nint(work_g1)
            write(nu,rec=nrec) work_gi8           
            deallocate(work_gi8)
         elseif (atype == 'rda4') then
            allocate(work_gr(nx_global,ny_global))
            work_gr = work_g1
            write(nu,rec=nrec) work_gr
            deallocate(work_gr)
         elseif (atype == 'rda8') then
            write(nu,rec=nrec) work_g1
         elseif (atype == 'ruf8') then
            write(nu) ((work_g1(i,j),i=1,nx_global),j=1,ny_global)
         else
            write(nu_diag,*) ' ERROR: writing unknown atype ',atype
         endif

    !-------------------------------------------------------------------
    ! diagnostics
    !-------------------------------------------------------------------
         if (diag) then
            amin = minval(work_g1)
            amax = maxval(work_g1)
            write(nu_diag,*) ' write_global ', nu, nrec, amin, amax
         endif

      endif                     ! my_task = master_task

      deallocate(work_g1)

      end subroutine ice_write

!=======================================================================

      end module ice_read_write

!=======================================================================
