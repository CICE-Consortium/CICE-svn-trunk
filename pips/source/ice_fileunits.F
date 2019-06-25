! CVS: $Id: ice_fileunits.F,v 1.4 2003/11/18 02:00:41 lipscomb Exp $
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
!
!  author: Elizabeth C. Hunke
!          Fluid Dynamics Group, Los Alamos National Laboratory
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

      integer (kind=int_kind), parameter :: 
     &   nu_grid        = 11  ! grid file
     &,  nu_kmt         = 12  ! land mask file
     &,  nu_nml         = 21  ! namelist input file
     &,  nu_diag        = 6   ! diagnostics output file
     &,  nu_forcing     = 49  ! forcing data file
     &,  nu_dump        = 50  ! dump file for restarting
     &,  nu_restart     = 50  ! restart input file
     &,  nu_rst_pointer = 52  ! pointer to latest restart file

!=======================================================================

      end module ice_fileunits

!=======================================================================
