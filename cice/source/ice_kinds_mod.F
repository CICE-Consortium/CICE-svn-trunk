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
! author: Elizabeth C. Hunke and William H. Lipscomb, LANL
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
     &                      char_len_long  = 256,
     &                      int_kind  = kind(1),
     &                      log_kind  = kind(.true.),
     &                      real_kind = selected_real_kind(6),
     &                      dbl_kind  = selected_real_kind(13),
     &                      quad_kind = selected_real_kind(20)

      ! additional definitions for modules imported from POP
      integer, parameter :: i4        = selected_int_kind(6),
     &                      i8        = selected_int_kind(13),
     &                      r4        = selected_real_kind(6),
     &                      r8        = selected_real_kind(13),
     &                      r16       = selected_real_kind(20)

!=======================================================================

      end module ice_kinds_mod

!=======================================================================
