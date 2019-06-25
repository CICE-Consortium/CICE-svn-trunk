! CVS: $Id: ice_kinds_mod.F,v 1.4 2003/11/18 02:00:41 lipscomb Exp $
!=======================================================================
!BOP
!
! !MODULE: ice_kinds_mod - defines variable precision
!
! !DESCRIPTION:
!
! Defines variable precision for all common data types \\
! Code originally based on kinds_mod.F in POP
!
! !REVISION HISTORY:
!
! author: Elizabeth C. Hunke
!         Fluid Dynamics Group, Los Alamos National Laboratory
!
! !INTERFACE:
!
      module ice_kinds_mod
!
! !USES:
!
!EOP
!=======================================================================

      implicit none
      save

      integer, parameter :: char_len  = 80,
     &                      char_len_long  = 128,
     &                      int_kind  = kind(1),
     &                      log_kind  = kind(.true.),
     &                      real_kind = selected_real_kind(6),
     &                      dbl_kind  = selected_real_kind(13)
 
!=======================================================================

      end module ice_kinds_mod

!=======================================================================
