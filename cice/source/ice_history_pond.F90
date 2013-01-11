!=======================================================================
!
!BOP
!
! !MODULE: ice_history - ice model history files
!
! Driver for core history output
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
      module ice_history_pond
!
! !USES:
!
      use ice_kinds_mod
      use ice_broadcast
      use ice_communicate, only: my_task, master_task
      use ice_blocks
!      use ice_read_write
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
           f_apondn    = 'm', f_apeffn     = 'm', &
           f_hpondn    = 'm',                     &
           f_apond     = 'x', f_apond_ai   = 'x', &
           f_hpond     = 'x', f_hpond_ai   = 'x', &
           f_ipond     = 'x', f_ipond_ai   = 'x', &
           f_apeff     = 'x', f_apeff_ai   = 'x'

      !---------------------------------------------------------------
      ! namelist variables
      !---------------------------------------------------------------

      namelist / icefields_pond_nml /     &
           f_apondn,    f_apeffn   , &
           f_hpondn,                 &
           f_apond,     f_apond_ai , &  
           f_hpond,     f_hpond_ai , &  
           f_ipond,     f_ipond_ai , &  
           f_apeff,     f_apeff_ai

      !---------------------------------------------------------------
      ! field indices
      !---------------------------------------------------------------

      integer (kind=int_kind), dimension(max_nstrm) :: &
           n_apondn      , n_apeffn    , & 
           n_hpondn      , &
           n_apond       , n_apond_ai, &
           n_hpond       , n_hpond_ai, &
           n_ipond       , n_ipond_ai, &
           n_apeff       , n_apeff_ai

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
      subroutine init_hist_pond
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
      use ice_state, only: tr_pond
      use ice_exit
!
! !INPUT/OUTPUT PARAMETERS:
!
!EOP
!
      integer (kind=int_kind) :: n, k, ns
      integer (kind=int_kind) :: nml_error ! namelist i/o error flag

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
            read(nu_nml, nml=icefields_pond_nml,iostat=nml_error)
            if (nml_error > 0) read(nu_nml,*)  ! for Nagware compiler
         end do
         if (nml_error == 0) close(nu_nml)
      endif
      call release_fileunit(nu_nml)

      call broadcast_scalar(nml_error, master_task)
      if (nml_error /= 0) then
         close (nu_nml)
         call abort_ice('ice: error reading icefields_pond_nml')
      endif

      if (.not. tr_pond) then
          f_apondn    = 'x'
          f_hpondn    = 'x'
          f_apeffn    = 'x'
          f_apond     = 'x'
          f_hpond     = 'x'
          f_ipond     = 'x'
          f_apeff     = 'x'
          f_apond_ai  = 'x'
          f_hpond_ai  = 'x'
          f_ipond_ai  = 'x'
          f_apeff_ai  = 'x'
      endif

#ifndef ncdf
      f_bounds = .false.
#endif

      call broadcast_scalar (f_apondn, master_task)
      call broadcast_scalar (f_hpondn, master_task)
      call broadcast_scalar (f_apeffn, master_task)
      call broadcast_scalar (f_apond, master_task)
      call broadcast_scalar (f_hpond, master_task)
      call broadcast_scalar (f_ipond, master_task)
      call broadcast_scalar (f_apeff, master_task)
      call broadcast_scalar (f_apond_ai, master_task)
      call broadcast_scalar (f_hpond_ai, master_task)
      call broadcast_scalar (f_ipond_ai, master_task)
      call broadcast_scalar (f_apeff_ai, master_task)

      ! 2D variables
      do ns = 1, nstreams

      if (f_apond(1:1) /= 'x') &
         call define_hist_field(n_apond,"apond","1",tstr2D, tcstr, &
             "melt pond fraction of sea ice",                      &
             "none", c1, c0,                                       &
             ns, f_apond)

      if (f_apond_ai(1:1) /= 'x') &
         call define_hist_field(n_apond_ai,"apond_ai","1",tstr2D, tcstr, & 
             "melt pond fraction of grid cell",                    &
             "weighted by ice area", c1, c0,                       &
             ns, f_apond)

      if (f_hpond(1:1) /= 'x') &
         call define_hist_field(n_hpond,"hpond","m",tstr2D, tcstr, &
             "mean melt pond depth over sea ice",                  &
             "none", c1, c0,                                       &
             ns, f_hpond)

      if (f_hpond_ai(1:1) /= 'x') &
         call define_hist_field(n_hpond_ai,"hpond_ai","m",tstr2D, tcstr, & 
             "mean melt pond depth over grid cell",                &
             "weighted by ice area", c1, c0,                       &
             ns, f_hpond)

      if (f_ipond(1:1) /= 'x') &
         call define_hist_field(n_ipond,"ipond","m",tstr2D, tcstr, &
             "mean pond ice thickness over sea ice",               &
             "none", c1, c0,                                       &
             ns, f_ipond)

      if (f_ipond_ai(1:1) /= 'x') &
         call define_hist_field(n_ipond_ai,"ipond_ai","m",tstr2D, tcstr, & 
             "mean pond ice thickness over grid cell",             &
             "weighted by ice area", c1, c0,                       &
             ns, f_ipond_ai)

      if (f_apeff(1:1) /= 'x') &
         call define_hist_field(n_apeff,"apeff","1",tstr2D, tcstr, &
             "radiation-effective pond area fraction of sea ice",  &
             "none", c1, c0,  &
             ns, f_apeff)

      if (f_apeff_ai(1:1) /= 'x') &
         call define_hist_field(n_apeff_ai,"apeff_ai","1",tstr2D, tcstr, &
             "radiation-effective pond area fraction over grid cell",  &
             "weighted by ice area", c1, c0,                       &
             ns, f_apeff_ai)

      enddo ! nstreams
      
      ! 3D (category) variables must be looped separately
      do ns = 1, nstreams

        if (f_apondn(1:1) /= 'x') &
           call define_hist_field(n_apondn,"apondn","1",tstr3Dc, tcstr, &
              "melt pond fraction, category","none", c1, c0,      &            
              ns, f_apondn)

        if (f_hpondn(1:1) /= 'x') &
           call define_hist_field(n_hpondn,"hpondn","m",tstr3Dc, tcstr, &
              "melt pond depth, category","none", c1, c0,       &
              ns, f_hpondn)

        if (f_apeffn(1:1) /= 'x') &
           call define_hist_field(n_apeffn,"apeffn","1",tstr3Dc, tcstr, &
             "effective melt pond fraction, category",   &
             "none", c1, c0,                                  &
             ns, f_apeffn)

      enddo ! ns

      end subroutine init_hist_pond

!=======================================================================
!
!BOP
!
! !IROUTINE: accum_hist - accumulate average ice quantities or snapshots
!
! !INTERFACE:
!
      subroutine accum_hist_pond (iblk)
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
      use ice_shortwave, only: apeffn
      use ice_work, only: worka
!
! !INPUT/OUTPUT PARAMETERS:
!
      integer (kind=int_kind), intent(in) :: &
           iblk                 ! block index
!
!EOP
!
      integer (kind=int_kind) :: &
           i,j, &
           ilo,ihi,jlo,jhi      ! beginning and end of physical domain

      type (block) :: &
         this_block           ! block information for current block

      !---------------------------------------------------------------
      ! increment field
      !---------------------------------------------------------------

         if (tr_pond_cesm) then
         if (f_apond(1:1)/= 'x') &
             call accum_hist_field(n_apond, iblk, &
                                   trcr(:,:,nt_apnd,iblk), a2D)
         if (f_apond_ai(1:1)/= 'x') &
             call accum_hist_field(n_apond_ai, iblk, &
                                   aice(:,:,iblk) * trcr(:,:,nt_apnd,iblk), a2D)
         if (f_hpond(1:1)/= 'x') &
             call accum_hist_field(n_hpond, iblk, &
                                   trcr(:,:,nt_apnd,iblk) &
                                 * trcr(:,:,nt_hpnd,iblk), a2D)
         if (f_hpond_ai(1:1)/= 'x') &
             call accum_hist_field(n_hpond_ai, iblk, &
                                   aice(:,:,iblk) * trcr(:,:,nt_apnd,iblk) &
                                                  * trcr(:,:,nt_hpnd,iblk), a2D)

         elseif (tr_pond_lvl) then
         if (f_apond(1:1)/= 'x') &
             call accum_hist_field(n_apond, iblk, &
                            trcr(:,:,nt_alvl,iblk) * trcr(:,:,nt_apnd,iblk), a2D)
         if (f_apond_ai(1:1)/= 'x') &
             call accum_hist_field(n_apond_ai, iblk, &
                            aice(:,:,iblk) &
                          * trcr(:,:,nt_alvl,iblk) * trcr(:,:,nt_apnd,iblk), a2D)
         if (f_hpond(1:1)/= 'x') &
             call accum_hist_field(n_hpond, iblk, &
                            trcr(:,:,nt_alvl,iblk) * trcr(:,:,nt_apnd,iblk) &
                                                   * trcr(:,:,nt_hpnd,iblk), a2D)
         if (f_hpond_ai(1:1)/= 'x') &
             call accum_hist_field(n_hpond_ai, iblk, &
                            aice(:,:,iblk) &
                          * trcr(:,:,nt_alvl,iblk) * trcr(:,:,nt_apnd,iblk) &
                                                   * trcr(:,:,nt_hpnd,iblk), a2D)
         if (f_ipond(1:1)/= 'x') &
             call accum_hist_field(n_ipond, iblk, &
                            trcr(:,:,nt_alvl,iblk) * trcr(:,:,nt_apnd,iblk) &
                                                   * trcr(:,:,nt_ipnd,iblk), a2D)
         if (f_ipond_ai(1:1)/= 'x') &
             call accum_hist_field(n_ipond_ai, iblk, &
                            aice(:,:,iblk) &
                          * trcr(:,:,nt_alvl,iblk) * trcr(:,:,nt_apnd,iblk) &
                                                   * trcr(:,:,nt_ipnd,iblk), a2D)

         elseif (tr_pond_topo) then

         if (f_apond(1:1)/= 'x') &
             call accum_hist_field(n_apond, iblk, &
                                   trcr(:,:,nt_apnd,iblk), a2D)
         if (f_apond_ai(1:1)/= 'x') &
             call accum_hist_field(n_apond_ai, iblk, &
                                   aice(:,:,iblk) * trcr(:,:,nt_apnd,iblk), a2D)
         if (f_hpond(1:1)/= 'x') &
             call accum_hist_field(n_hpond, iblk, &
                                   trcr(:,:,nt_apnd,iblk) &
                                 * trcr(:,:,nt_hpnd,iblk), a2D)
         if (f_hpond_ai(1:1)/= 'x') &
             call accum_hist_field(n_hpond_ai, iblk, &
                                   aice(:,:,iblk) * trcr(:,:,nt_apnd,iblk) &
                                                  * trcr(:,:,nt_hpnd,iblk), a2D)
         if (f_ipond(1:1)/= 'x') &
             call accum_hist_field(n_ipond, iblk, &
                                   trcr(:,:,nt_apnd,iblk) &
                                 * trcr(:,:,nt_ipnd,iblk), a2D)
         if (f_ipond_ai(1:1)/= 'x') &
             call accum_hist_field(n_ipond_ai, iblk, &
                                   aice(:,:,iblk) * trcr(:,:,nt_apnd,iblk) &
                                                  * trcr(:,:,nt_ipnd,iblk), a2D)
         endif ! ponds

         this_block = get_block(blocks_ice(iblk),iblk)         
         ilo = this_block%ilo
         ihi = this_block%ihi
         jlo = this_block%jlo
         jhi = this_block%jhi

         if (f_apeff (1:1) /= 'x') then
             worka(:,:) = c0
             do j = jlo, jhi
             do i = ilo, ihi
                if (aice(i,j,iblk) > puny) worka(i,j) = apeff_ai(i,j,iblk) &
                                                      / aice(i,j,iblk)
             enddo
             enddo
             call accum_hist_field(n_apeff, iblk, worka(:,:), a2D)
         endif
         if (f_apeff_ai(1:1) /= 'x') &
             call accum_hist_field(n_apeff_ai, iblk, apeff_ai(:,:,iblk), a2D)

         ! 3D category fields
         if (f_apondn   (1:1) /= 'x') &
             call accum_hist_field(n_apondn-n2D, iblk, ncat_hist, &
                  trcrn(:,:,nt_apnd,1:ncat_hist,iblk), a3Dc)
         if (f_apeffn (1:1) /= 'x') &
             call accum_hist_field(n_apeffn-n2D,  iblk, ncat_hist, &
                  apeffn(:,:,1:ncat_hist,iblk), a3Dc)
         if (f_hpondn   (1:1) /= 'x') &
             call accum_hist_field(n_hpondn-n2D, iblk, ncat_hist, &
                    trcrn(:,:,nt_apnd,1:ncat_hist,iblk) &
                  * trcrn(:,:,nt_hpnd,1:ncat_hist,iblk), a3Dc)

      end subroutine accum_hist_pond

!=======================================================================

      end module ice_history_pond

!=======================================================================
