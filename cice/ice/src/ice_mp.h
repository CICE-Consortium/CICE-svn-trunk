c $Id: ice_mp.h,v 1.2 1998/04/14 14:59:52 eclare Exp $
c.. Common blocks for message passing
c-----------------------------------------------------------------------
      integer nibuff                     ! size of integer control buffer
      parameter (nibuff=20)            
      integer ibuff  (nibuff)            ! integer control buffer
      integer ibuff_d(nibuff)            ! control buffer from drv

      integer msg_id                     ! pointer to message channel 
      integer msgtype_d2i                ! message type for drv->ice
      integer msgtype_i2d                ! message type for ice->drv
      integer bufidi(MCL_SIZEOF_DESC),   ! descriptor for initial buffer
     $        bufids(MCL_SIZEOF_DESC),   ! descriptor for send    buffer 
     $        bufidr(MCL_SIZEOF_DESC),   ! descriptor for recv    buffer
     $        bufidc(MCL_SIZEOF_DESC)    ! descriptor for control buffer
      integer nsnd, nrcv                 ! number of fields sent/received
      integer nadv_i                     ! number of timesteps per day

      common /ice_mp/ ibuff, ibuff_d, msg_id, msgtype_d2i, 
     & msgtype_i2d, bufidi, bufids, bufidr, bufidc, nsnd, nrcv, 
     & nadv_i

