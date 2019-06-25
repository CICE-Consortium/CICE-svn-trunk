c $Id: $
c=======================================================================
!---! Sets array sizes for local subdomain and related parallel 
!---! processing info
!---!
!---! author Elizabeth C. Hunke, LANL
!---!
!---! code originally based on domain.F in POP
c=======================================================================

      module ice_domain

      use ice_kinds_mod
      use ice_model_size

      implicit none
      save

      ! ghost cells define boundary data to avoid frequent 
      ! off-processor communication 
      integer (kind=int_kind), parameter ::  
     &  num_ghost_cells = 1

      ! imt_local, jmt_local define total size of subdomain - 
      ! the formula below automatically pads arrays if NPROC_X
      ! or NPROC_Y do not evenly divide global problem size
      integer (kind=int_kind), parameter ::
     &  imt_local = (imt_global-1)/NPROC_X + 1 + 2*num_ghost_cells,
     &  jmt_local = (jmt_global-1)/NPROC_Y + 1 + 2*num_ghost_cells

      integer (kind=int_kind), parameter :: 
     &  ilo = num_ghost_cells+1,  ! beg index of actual physical subdomain
     &  ihi = imt_local - num_ghost_cells, ! end index 
     &  jlo = num_ghost_cells+1,           ! beg index
     &  jhi = jmt_local - num_ghost_cells  ! end index 

      integer (kind=int_kind) :: 
     &  i_global(imt_local),    ! global index for local subdomain
     &  j_global(jmt_local)     ! global index for local subdomain

      integer (kind=int_kind)  ::
     &  my_task,          ! task id for local process
     &  master_task       ! task id for master process

c=======================================================================

      end module ice_domain

c=======================================================================
