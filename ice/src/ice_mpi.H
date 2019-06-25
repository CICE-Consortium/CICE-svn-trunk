c $Id: ice_mpi.H,v 2.0 1999/12/09 22:25:35 eclare Exp $
c---------------------------------------------------------------------
c  parameters and commons blocks for MPI
c---------------------------------------------------------------------
#ifdef _MPI
      integer ::
     &   nbr_west, nbr_east,        ! processor numbers for the 
     &   nbr_north, nbr_south       ! n,s,e,w neighbor processors

      integer, dimension(2,nproc_s) ::
     &   local_start                ! global i,j start for each processor

      integer, dimension(nproc_s) ::
     &   mpi_interior_int,          ! MPI data type definitions for
     &   mpi_interior_real,         ! physical domain blocks
     &   mpi_interior_int_global,   ! in case of padding
     &   mpi_interior_real_global

      common / cice_mpi /
     &   nbr_west, nbr_east,
     &   nbr_north, nbr_south,
     &   local_start,
     &   mpi_interior_int,    
     &   mpi_interior_real,   
     &   mpi_interior_int_global,
     &   mpi_interior_real_global
#endif

