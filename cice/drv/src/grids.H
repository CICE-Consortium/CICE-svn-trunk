c=======================================================================
c SCCS module : grids.h   r3.1
c SCCS get    : 03 May 1996 / 16:59:35
c SCCS edit   : 03 May 1996 / 16:59:23
c SCCS path   : /net/crestone.u1/nieman/csm/drv-3.n/SCCS/s.grids.h
c=======================================================================

      !----- model coordinate arrays and model domain masks -----
      REAL
     &   lat_i(nx_i, ny_i), lon_i(nx_i, ny_i)

      INTEGER
     &   mask_i(nx_i,ny_i)

      COMMON / grids / 
     &   lat_i, lon_i,
     &   mask_i
      
      SAVE / grids /

c-----------------------------------------------------------------------
