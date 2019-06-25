c=======================================================================
c SCCS module : buffer.h   r3.1
c SCCS get    : 03 May 1996 / 16:59:34
c SCCS edit   : 03 May 1996 / 16:59:18
c SCCS path   : /net/crestone.u1/nieman/csm/drv-3.n/SCCS/s.buffer.h
c=======================================================================
c edited by E. C. Hunke 9 April 1998

      !----- output state buffers --------------------------------------

      REAL
     &   Sa_i(nx_i,ny_i,nSa),
     &   Si_i(nx_i,ny_i,nSi),
     &   So_i(nx_i,ny_i,nSo)

      !----- input flux buffers ----------------------------------------

      REAL    
     &   Fixx(nx_i,ny_i,nFixx)

      !----- flux component buffers ------------------------------------

      !----- atm/ice flux components -----
      REAL
     &   Faia_i(nx_i,ny_i,nFaia),
     &   Faii_i(nx_i,ny_i,nFaii),  ! ice output
     &   Faid_i(nx_i,ny_i,nFaid) 

      !----- ice/ocn flux components -----
      REAL
     &   Fioi_i(nx_i,ny_i,nFioi),  ! ice output
     &   Fioo_i(nx_i,ny_i,nFioo),
     &   Fiod_i(nx_i,ny_i,nFiod)

      !----- provide global access to buffers --------------------------

      COMMON /buffer/
     &   Sa_i, 
     &   Si_i,
     &   So_i,
     &   Fixx,
     &   Faia_i,  Faii_i,  Faid_i,
     &   Fioi_i,  Fioo_i,  Fiod_i

      SAVE /buffer/
