!=======================================================================
!
!BOP
!
! !MODULE: ice_history - ice model history files
!
! Biogeochemistry history output
!
! The following variables are currently hard-wired as snapshots 
!   (instantaneous rather than time-averages):
!   divu, shear, sig1, sig2, trsig, mlt_onset, frz_onset, hisnap, aisnap
!
! The flags (f_<field>) can be set to '1','h','d','m','y' or 'x', where
!   n means the field will not be written.  To output the same field at
!   more than one frequency, for instance monthy and daily, set 
!   f_<field> = 'md'.
!
! !REVISION HISTORY:
!  SVN:$Id: ice_history.F90 569 2013-01-10 15:28:29Z eclare $
!
! authors Tony Craig and Bruce Briegleb, NCAR
!         Elizabeth C. Hunke and William H. Lipscomb, LANL
!         C. M. Bitz, UW
!
! 2012 Elizabeth Hunke split code from ice_history.F90
!
! !INTERFACE:
!
      module ice_history_bgc
!
! !USES:
!
      use ice_kinds_mod
      use ice_broadcast
      use ice_communicate, only: my_task, master_task
      use ice_blocks
      use ice_fileunits
      use ice_history_shared
      use ice_history_write
!
!EOP
!
      implicit none
      save
      
      !---------------------------------------------------------------
      ! flags: write to output file if true or histfreq value
      !---------------------------------------------------------------

      character (len=max_nstrm) :: &
           f_faero_atm = 'm', f_faero_ocn  = 'm', &
           f_aero      = 'm', f_aeron      = 'm'

      !---------------------------------------------------------------
      ! namelist variables
      !---------------------------------------------------------------

      namelist / icefields_bgc_nml /     &
           f_faero_atm, f_faero_ocn, &
           f_aero,      f_aeron    

      !---------------------------------------------------------------
      ! field indices
      !---------------------------------------------------------------

      ! aerosols
      integer(kind=int_kind), dimension(max_aero,max_nstrm) :: &
           n_faero_atm    , &
           n_faero_ocn    , &
           n_aerosn1      , &
           n_aerosn2      , &
           n_aeroic1      , &
           n_aeroic2

!=======================================================================

      contains

!=======================================================================
!
!BOP
!
! !IROUTINE: init_hist - initialize history files
!
! !INTERFACE:
!
      subroutine init_hist_bgc
!
! !DESCRIPTION:
!
! Initialize history files
!
! !REVISION HISTORY:
!
! authors Elizabeth C. Hunke, LANL
!
! !USES:
!
      use ice_constants
      use ice_calendar, only: nstreams
      use ice_state, only: tr_aero
      use ice_exit
!
! !INPUT/OUTPUT PARAMETERS:
!
!EOP
!
      integer (kind=int_kind) :: n, k, ns
      integer (kind=int_kind) :: nml_error ! namelist i/o error flag
      character (len=3) :: nchar

      !-----------------------------------------------------------------
      ! read namelist
      !-----------------------------------------------------------------

      call get_fileunit(nu_nml)
      if (my_task == master_task) then
         open (nu_nml, file=nml_filename, status='old',iostat=nml_error)
         if (nml_error /= 0) then
            nml_error = -1
         else
            nml_error =  1
         endif
         do while (nml_error > 0)
            read(nu_nml, nml=icefields_bgc_nml,iostat=nml_error)
            if (nml_error > 0) read(nu_nml,*)  ! for Nagware compiler
         end do
         if (nml_error == 0) close(nu_nml)
      endif
      call release_fileunit(nu_nml)

      call broadcast_scalar(nml_error, master_task)
      if (nml_error /= 0) then
         close (nu_nml)
         call abort_ice('ice: error reading icefields_bgc_nml')
      endif

      if (.not. tr_aero) then
         f_faero_atm = 'x'
         f_faero_ocn = 'x'
         f_aero      = 'x' 
         f_aeron     = 'x' ! NOTE not implemented
      endif

      call broadcast_scalar (f_faero_atm, master_task)
      call broadcast_scalar (f_faero_ocn, master_task)
      call broadcast_scalar (f_aero, master_task)
      call broadcast_scalar (f_aeron, master_task)

      ! 2D variables
      do ns = 1, nstreams

      ! Aerosols
      if (f_aero(1:1) /= 'x') then
         do n=1,n_aero
            write(nchar,'(i3.3)') n
            write(vname_in,'(a,a)') 'aerosnossl', trim(nchar)
            call define_hist_field(n_aerosn1(n,:),vname_in,"kg/kg",   &
                tstr2D, tcstr,"snow ssl aerosol mass","none", c1, c0, &
                ns, f_aero)
            write(vname_in,'(a,a)') 'aerosnoint', trim(nchar)
            call define_hist_field(n_aerosn2(n,:),vname_in,"kg/kg",   &
                tstr2D, tcstr,"snow int aerosol mass","none", c1, c0, &
                ns, f_aero)
            write(vname_in,'(a,a)') 'aeroicessl', trim(nchar)
            call define_hist_field(n_aeroic1(n,:),vname_in,"kg/kg",  &
                tstr2D, tcstr,"ice ssl aerosol mass","none", c1, c0, &
                ns, f_aero)
            write(vname_in,'(a,a)') 'aeroiceint', trim(nchar)
            call define_hist_field(n_aeroic2(n,:),vname_in,"kg/kg",  &
                tstr2D, tcstr,"ice int aerosol mass","none", c1, c0, &
                ns, f_aero)
         enddo
      endif

      if (f_faero_atm(1:1) /= 'x') then
         do n=1,n_aero
            write(nchar,'(i3.3)') n
            write(vname_in,'(a,a)') 'faero_atm', trim(nchar)
            call define_hist_field(n_faero_atm(n,:),vname_in,"kg/m^2 s", &
                tstr2D, tcstr,"aerosol deposition rate","none", c1, c0,  &
                ns, f_faero_atm)
         enddo
      endif

      if (f_faero_ocn(1:1) /= 'x') then
         do n=1,n_aero
            write(nchar,'(i3.3)') n
            write(vname_in,'(a,a)') 'faero_ocn', trim(nchar)
            call define_hist_field(n_faero_ocn(n,:),vname_in,"kg/m^2 s", &
                tstr2D, tcstr,"aerosol flux to ocean","none", c1, c0,    &
                ns, f_faero_ocn)
         enddo
      endif

      enddo ! nstreams
      
      ! 3D (category) variables must be looped separately
!      do ns = 1, nstreams
!      enddo ! ns

      end subroutine init_hist_bgc

!=======================================================================
!
!BOP
!
! !IROUTINE: accum_hist - accumulate average ice quantities or snapshots
!
! !INTERFACE:
!
      subroutine accum_hist_bgc (iblk)
!
! !DESCRIPTION:
!
! write average ice quantities or snapshots
!
! !REVISION HISTORY:
!
! author:   Elizabeth C. Hunke, LANL
!
! !USES:
!
      use ice_blocks
      use ice_domain
      use ice_state
      use ice_flux
!
! !INPUT/OUTPUT PARAMETERS:
!
      integer (kind=int_kind), intent(in) :: &
           iblk                 ! block index
!
!EOP
!
      integer (kind=int_kind) :: &
           i,j,n, &
           ilo,ihi,jlo,jhi      ! beginning and end of physical domain

      type (block) :: &
         this_block           ! block information for current block

      !---------------------------------------------------------------
      ! increment field
      !---------------------------------------------------------------

        ! Aerosols
        if (f_faero_atm(1:1) /= 'x') then
           do n=1,n_aero
              call accum_hist_field(n_faero_atm(n,:),iblk, &
                                      faero_atm(:,:,n,iblk), a2D)
           enddo
        endif
        if (f_faero_ocn(1:1) /= 'x') then
           do n=1,n_aero
              call accum_hist_field(n_faero_ocn(n,:),iblk, &
                                      faero_ocn(:,:,n,iblk), a2D)
                                    
           enddo
        endif
        if (f_aero(1:1) /= 'x') then
           do n=1,n_aero
              call accum_hist_field(n_aerosn1(n,:), iblk, &
                                 trcr(:,:,nt_aero  +4*(n-1),iblk)/rhos, a2D)
              call accum_hist_field(n_aerosn2(n,:), iblk, &
                                 trcr(:,:,nt_aero+1+4*(n-1),iblk)/rhos, a2D)
              call accum_hist_field(n_aeroic1(n,:), iblk, &
                                 trcr(:,:,nt_aero+2+4*(n-1),iblk)/rhoi, a2D)
              call accum_hist_field(n_aeroic2(n,:), iblk, &
                                 trcr(:,:,nt_aero+3+4*(n-1),iblk)/rhoi, a2D)
           enddo
        endif

      end subroutine accum_hist_bgc

!=======================================================================

      end module ice_history_bgc

!=======================================================================
