!=======================================================================
!BOP
!
! !MODULE: ice_fileunits
!
! !DESCRIPTION:
!
! Defines unit numbers for files opened for reading or writing
!
! !REVISION HISTORY:
!  SVN:$Id$
!  SVN:$LastChangedBy$
!
!  author: Elizabeth C. Hunke, LANL
!
! !INTERFACE:
!
      module ice_fileunits
!
! !USES:
      use ice_kinds_mod
!
!EOP
!=======================================================================

      implicit none
      save

      character (len=char_len) ::
     &   diag_type            ! 'stdout' or 'file'

      integer (kind=int_kind), parameter ::
     &   nu_grid        = 11  ! grid file
     &,  nu_kmt         = 12  ! land mask file
     &,  nu_nml         = 21  ! namelist input file
     &,  nu_forcing     = 49  ! forcing data file
     &,  nu_dump        = 50  ! dump file for restarting
     &,  nu_restart     = 50  ! restart input file
     &,  nu_rst_pointer = 52  ! pointer to latest restart file

      integer (kind=int_kind) :: 
     &   nu_diag              ! diagnostics output file

      character (6), parameter ::
     &   nml_filename = 'ice_in'  ! namelist input file name

!=======================================================================

      end module ice_fileunits

!=======================================================================
