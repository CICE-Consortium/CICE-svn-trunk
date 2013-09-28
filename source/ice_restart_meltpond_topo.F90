!  SVN:$Id$
!=======================================================================
!
! Topo melt pond parameterization
!
! Melt pond evolution based on the ice topography as inferred from
! the ice thickness distribution.  This code is based on (but differs
! from) that described in
!
! Flocco, D. and D. L. Feltham, 2007.  A continuum model of melt pond 
! evolution on Arctic sea ice.  J. Geophys. Res. 112, C08016, doi: 
! 10.1029/2006JC003836.
!
! Flocco, D., D. L. Feltham and A. K. Turner, 2010.  Incorporation of a
! physically based melt pond scheme into the sea ice component of a
! climate model.  J. Geophys. Res. 115, C08012, doi: 10.1029/2009JC005568.
!
! authors Daniela Flocco (UCL)
!         Adrian Turner (UCL)
! 2010 ECH added module based on original code from Daniela Flocco, UCL
! 2012 DSCHR modifications

      module ice_restart_meltpond_topo

      use ice_kinds_mod

      implicit none
      private
      public :: write_restart_pond_topo, read_restart_pond_topo
      save

      logical (kind=log_kind), public :: & 
         restart_pond_topo ! if .true., read meltponds restart file

      real (kind=dbl_kind), public :: &
         hp1               ! critical parameter for pond ice thickness

!=======================================================================

      contains
  
!=======================================================================

! Dumps all values needed for restarting
!
! authors Elizabeth C. Hunke, LANL
!         David A. Bailey, NCAR

      subroutine write_restart_pond_topo()

      use ice_domain_size, only: ncat
      use ice_fileunits, only: nu_dump_pond
      use ice_state, only: trcrn, nt_apnd, nt_hpnd, nt_ipnd
      use ice_restart, only: write_restart_field

      ! local variables

      logical (kind=log_kind) :: diag

      diag = .true.

      call write_restart_field(nu_dump_pond,0,trcrn(:,:,nt_apnd,:,:),'ruf8', &
                               'apnd',ncat,diag)
      call write_restart_field(nu_dump_pond,0,trcrn(:,:,nt_hpnd,:,:),'ruf8', &
                               'hpnd',ncat,diag)
      call write_restart_field(nu_dump_pond,0,trcrn(:,:,nt_ipnd,:,:),'ruf8', &
                               'ipnd',ncat,diag)

      end subroutine write_restart_pond_topo

!=======================================================================

! Reads all values needed for a meltpond volume restart
!
! authors Elizabeth C. Hunke, LANL
!         David A. Bailey, NCAR

      subroutine read_restart_pond_topo()

      use ice_domain_size, only: ncat
      use ice_communicate, only: my_task, master_task
      use ice_fileunits, only: nu_diag, nu_restart_pond 
      use ice_state, only: trcrn, nt_apnd, nt_hpnd, nt_ipnd
      use ice_restart, only: read_restart_field

      ! local variables

      logical (kind=log_kind) :: &
         diag

      diag = .true.

      if (my_task == master_task) write(nu_diag,*) 'min/max topo ponds'

      call read_restart_field(nu_restart_pond,0,trcrn(:,:,nt_apnd,:,:),'ruf8', &
                              'apnd',ncat,diag)
      call read_restart_field(nu_restart_pond,0,trcrn(:,:,nt_hpnd,:,:),'ruf8', &
                              'hpnd',ncat,diag)
      call read_restart_field(nu_restart_pond,0,trcrn(:,:,nt_ipnd,:,:),'ruf8', &
                              'ipnd',ncat,diag)

      end subroutine read_restart_pond_topo

!=======================================================================

      end module ice_restart_meltpond_topo

!=======================================================================
