c=======================================================================
c SCCS module : msg.MCL1.0.h   r3.2
c SCCS get    : 03 Jun 1996 / 17:19:37
c SCCS edit   : 03 Jun 1996 / 17:19:29
c SCCS path   : /net/crestone.u1/nieman/csm/drv-3.n/SCCS/s.msg.MCL1.0.h
c=======================================================================
c edited by E. C. Hunke 9 April 1998

      !-----------------------------------------------------------------
      ! include msg passing lib header files 
      !-----------------------------------------------------------------
      include '/usr/array/PVM/include/fpvm3.h' ! pvm3 includes
      include 'mcl.h'                          ! mcl  includes

      !-----------------------------------------------------------------
      ! The driver and component models agree that the atm is "model 
      ! number" 6, the ocn is number 7, and the ice is number 8.
      ! MCL uses these agreements to establish a connection between the
      ! drv an a particular model.
      !-----------------------------------------------------------------
      integer    model_num_ice
      parameter( model_num_ice = 8)  ! model number 8 => ice model

      !-----------------------------------------------------------------
      ! MCL uses these "task identifiers" to store information about
      ! how to communicate with the various component models (tasks)
      !-----------------------------------------------------------------
      integer tid_ice  ! ice model task identifier

      !-----------------------------------------------------------------
      ! MCL uses these "data descriptors" to store pointers to
      ! dynamically allocated memory which may be used & needed by the 
      ! MCL library, depending on what msg passing method is selected.
      ! The memory is allocated by a call to mcl_alloc.  
      ! If PVM (only) is used as the msg passing method, then no memory
      ! is needed and mcl_alloc does not allocate any memory.
      !-----------------------------------------------------------------
      integer dd_icec(MCL_SIZEOF_DESC) ! data descriptor for cntl buff
      integer dd_icei(MCL_SIZEOF_DESC) ! data descriptor for init recv
      integer dd_ices(MCL_SIZEOF_DESC) ! data descriptor for ice send
      integer dd_icer(MCL_SIZEOF_DESC) ! data descriptor for ice recv
      integer dd_null(MCL_SIZEOF_DESC) ! dummy/unused data descriptor

      COMMON / msg_mcl /
     &   tid_ice ,
     &   dd_icec , dd_icei , dd_ices , dd_icer ,
     &   dd_null
      SAVE   / msg_mcl /

      !-----------------------------------------------------------------
      ! define message ID numbers 
      !-----------------------------------------------------------------
      integer    msgid_d2i_init  ! msg id for drv ->ice msg, initial msg
      integer    msgid_d2i       ! msg id for drv ->ice msg
      integer    msgid_i2d_init  ! msg id for drv<- ice msg, initial msg
      integer    msgid_i2d       ! msg id for drv<- ice msg
      parameter( msgid_d2i_init = 31 )
      parameter( msgid_d2i      = 30 )
      parameter( msgid_i2d_init = 40 )
      parameter( msgid_i2d      = 41 )

c=======================================================================
