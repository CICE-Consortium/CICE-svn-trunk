c $Id: $
c=======================================================================
!---! Defines the global domain size and number of categories
!---! and layers
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
!---! code originally based on model_size.F in POP
c=======================================================================

      module ice_model_size

      use ice_kinds_mod

      implicit none
      save

      integer (kind=int_kind), parameter :: 
     &  imt_global = 100,  ! i-axis size 
     &  jmt_global = 116,  ! j-axis size
     &  ncat       =   5,  ! number of categories
     &  nlyr       =   4,  ! number of layers in a single category
     &  ntilay     =  ncat*nlyr  ! sum of number of layers in all categories

c=======================================================================

      end module ice_model_size

c=======================================================================
