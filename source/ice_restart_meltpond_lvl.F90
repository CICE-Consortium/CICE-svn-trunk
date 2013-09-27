!  SVN:$Id$
!=======================================================================
!
! Level-ice melt pond scheme
!
! This meltpond parameterization was developed for use with the delta-
! Eddington radiation scheme, and only affects the radiation budget in
! the model.  That is, although the pond volume is tracked, that liquid
! water is not used elsewhere in the model for mass budgets or other
! physical processes.
!
! authors Elizabeth Hunke (LANL)
!         David Hebert (NRL Stennis)
!         Olivier Lecomte (Univ. Louvain)

      module ice_restart_meltpond_lvl

      use ice_kinds_mod
      use ice_blocks, only: nx_block, ny_block
      use ice_domain_size, only: max_blocks, ncat

      implicit none
      private
      public :: write_restart_pond_lvl, read_restart_pond_lvl

      logical (kind=log_kind), public :: & 
         restart_pond_lvl, & ! if .true., read meltponds restart file
         snowinfil           ! if .true., adjust snow depth/area in dEdd
                             !            for infiltration of melt water

      character (len=char_len), public :: &
         frzpnd              ! pond refreezing parameterization

      real (kind=dbl_kind), public :: &
         dpscale, &          ! alter e-folding time scale for flushing 
         rfracmin, &         ! minimum retained fraction of meltwater
         rfracmax, &         ! maximum retained fraction of meltwater
         pndaspect, &        ! ratio of pond depth to pond fraction
         hs1                 ! tapering parameter for snow on pond ice

      real (kind=dbl_kind), public, &
         dimension (nx_block,ny_block,ncat,max_blocks) :: &
         dhsn, &      ! depth difference for snow on sea ice and pond ice
         ffracn       ! fraction of fsurfn used to melt ipond

!=======================================================================

      contains
  
!=======================================================================
!
! Dumps all values needed for restarting
!
! authors Elizabeth C. Hunke, LANL

      subroutine write_restart_pond_lvl()

      use ice_domain_size, only: ncat
      use ice_fileunits, only: nu_dump_pond
      use ice_state, only: trcrn, nt_apnd, nt_hpnd, nt_ipnd
      use ice_restart, only: write_restart_field

      ! local variables

      logical (kind=log_kind) :: diag

      diag = .true.

      call write_restart_field(nu_dump_pond,0, trcrn(:,:,nt_apnd,:,:),'ruf8', &
                               'apnd',ncat,diag)
      call write_restart_field(nu_dump_pond,0, trcrn(:,:,nt_hpnd,:,:),'ruf8', &
                               'hpnd',ncat,diag)
      call write_restart_field(nu_dump_pond,0, trcrn(:,:,nt_ipnd,:,:),'ruf8', &
                               'ipnd',ncat,diag)
      call write_restart_field(nu_dump_pond,0,  dhsn(:,:,        :,:),'ruf8', &
                               'dhs',ncat,diag)
      call write_restart_field(nu_dump_pond,0,ffracn(:,:,        :,:),'ruf8', &
                               'ffrac',ncat,diag)

      end subroutine write_restart_pond_lvl

!=======================================================================

! Reads all values needed for a meltpond volume restart
!
! authors Elizabeth C. Hunke, LANL

      subroutine read_restart_pond_lvl()

      use ice_domain_size, only: ncat
      use ice_communicate, only: my_task, master_task
      use ice_fileunits, only: nu_diag, nu_restart_pond 
      use ice_state, only: trcrn, nt_apnd, nt_hpnd, nt_ipnd
      use ice_restart, only: read_restart_field

      ! local variables

      logical (kind=log_kind) :: &
         diag

      diag = .true.

      if (my_task == master_task) write(nu_diag,*) 'level-ice ponds'

      call read_restart_field(nu_restart_pond,0, trcrn(:,:,nt_apnd,:,:),'ruf8', &
                              'apnd',ncat,diag)
      call read_restart_field(nu_restart_pond,0, trcrn(:,:,nt_hpnd,:,:),'ruf8', &
                              'hpnd',ncat,diag)
      call read_restart_field(nu_restart_pond,0, trcrn(:,:,nt_ipnd,:,:),'ruf8', &
                              'ipnd',ncat,diag)
      call read_restart_field(nu_restart_pond,0,  dhsn(:,:,        :,:),'ruf8', &
                              'dhs',ncat,diag)
      call read_restart_field(nu_restart_pond,0,ffracn(:,:,        :,:),'ruf8', &
                              'ffrac',ncat,diag)

      end subroutine read_restart_pond_lvl

!=======================================================================

      end module ice_restart_meltpond_lvl

!=======================================================================
