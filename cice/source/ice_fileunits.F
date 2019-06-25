c $Id: $
c=======================================================================
!---! Defines unit numbers for files opened for reading or writing
!---!
!---! author Elizabeth C. Hunke
!---!        Fluid Dynamics Group, Los Alamos National Laboratory
c=======================================================================

      module ice_fileunits

      use ice_kinds_mod

      implicit none
      save

      integer (kind=int_kind), parameter :: 
     &   nu_grid        = 11  ! grid file
     &,  nu_kmt         = 12  ! land mask file
     &,  nu_nml         = 21  ! namelist input file
     &,  nu_diag        = 6   ! diagnostics output file
c     &,  nu_diag        = 48  ! diagnostics output file
     &,  nu_dump        = 50  ! dump file for restarting
     &,  nu_restart     = 50  ! restart input file
     &,  nu_rst_pointer = 52  ! pointer to latest restart file

c=======================================================================

      end module ice_fileunits

c=======================================================================
