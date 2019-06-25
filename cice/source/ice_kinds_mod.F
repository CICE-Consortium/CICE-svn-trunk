c $Id: $
c=======================================================================
!---! Defines variable precision for all common data types
!---!
!---! author Elizabeth C. Hunke
!---!        Fluid Dynamics Group, Los Alamos National Laboratory
!---!
!---! Copyright, 2000.  The Regents of the University of California.
!---! This software was produced under a U.S. Government contract 
!---! (W-7405-ENG-36) by Los Alamos National Laboratory, which is 
!---! operated by the University of California for the U.S. Department 
!---! of Energy.  The U.S. Government is licensed to use, reproduce, and 
!---! distribute this software.  Permission is granted to the public to 
!---! copy and use this software without charge, provided that this 
!---! Notice and any statement of authorship are reproduced on all 
!---! copies.  Neither the Government nor the University makes any 
!---! warranty, express or implied, or assumes any liability or 
!---! responsibility for the use of this software.
!---!
!---! code originally based on kinds_mod.F in POP
c=======================================================================

      module ice_kinds_mod

      implicit none
      save

      integer, parameter :: char_len  = 80,
     &                      int_kind  = kind(1),
     &                      log_kind  = kind(.true.),
     &                      real_kind = selected_real_kind(6),
     &                      dbl_kind  = selected_real_kind(13)

      end module ice_kinds_mod

c=======================================================================

