c=======================================================================
c SCCS module : dims.h   r3.1
c SCCS get    : 03 May 1996 / 16:59:34
c SCCS edit   : 03 May 1996 / 16:59:20
c SCCS path   : /net/crestone.u1/nieman/csm/drv-3.n/SCCS/s.dims.h
c=======================================================================
c edited by E. C. Hunke 9 April 1998

      !----- computational domain dimensions ---------------------------
      INTEGER   nx_i   , ny_i
      PARAMETER(nx_i=192, ny_i=128)

      !----- output state buffers dimensions ---------------------------
      INTEGER    nSa
      INTEGER    nSi
      INTEGER    nSo

      PARAMETER (nSa= 6)
      PARAMETER (nSi= 8)
      PARAMETER (nSo= 8)

      !----- input flux buffers dimensions -----------------------------
      INTEGER    nFixx
      PARAMETER (nFixx= 8)

      !----- flux component buffers dimensions -------------------------
      INTEGER    nFaia   , nFaii   , nFaid  
      INTEGER    nFioi   , nFioo   , nFiod
 
      PARAMETER (nFaia= 3, nFaii= 4, nFaid= 2)
      PARAMETER (nFioi= 5, nFioo= 2, nFiod= 1)

