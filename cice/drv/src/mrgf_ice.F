c=======================================================================
c SCCS module : mrgf_ice.f   r3.1
c SCCS get    : 03 May 1996 / 17:25:00
c SCCS edit   : 03 May 1996 / 17:24:24
c SCCS path   : /net/crestone.u1/nieman/csm/drv-3.n/SCCS/s.mrgf_ice.f
c=======================================================================
c  edited by E. C. Hunke, 8 April 1998

      SUBROUTINE mrgf_ice(Faia,Faid,Fioo,Fiod, Fixx)

      implicit none

      !----- global -----
      include 'dims.h'          ! buffer dimensions
      include 'index.h'         ! decodes buffer contents

      !----- input -----
      real
     $  Faia(nx_i,ny_i,nFaia),  ! atm/ice flux, computed by atm
     $  Faid(nx_i,ny_i,nFaid),  ! atm/ice flux, computed by drv
     $  Fioo(nx_i,ny_i,nFioo),  ! ice/ocn flux, computed by ocn
     $  Fiod(nx_i,ny_i,nFiod)   ! ice/ocn flux, computed by drv

      !----- output -----
      real
     $  Fixx(nx_i,ny_i,nFixx)   ! all ice flux, ice input buffer

      !----- local -----
      integer
     $  i,j                     ! loop indicies

c-----------------------------------------------------------------------
c PURPOSE:
c   gather various flux component fields (all on the ice grid)
c-----------------------------------------------------------------------

c$    doacross local(i,j)
      do j=1,ny_i
      do i=1,nx_i
        Fixx(i,j,kFixx_taux ) = Faid(i,j,kFaid_taux  )
        Fixx(i,j,kFixx_tauy ) = Faid(i,j,kFaid_tauy  )
        Fixx(i,j,kFixx_dhdx ) = Fioo(i,j,kFioo_dhdx  )
        Fixx(i,j,kFixx_dhdy ) = Fioo(i,j,kFioo_dhdy  )
        Fixx(i,j,kFixx_netsw) = Faia(i,j,kFaia_netsw )
        Fixx(i,j,kFixx_lwd  ) = Faia(i,j,kFaia_lwdn  )
        Fixx(i,j,kFixx_water) = Faia(i,j,kFaia_prec  )
        Fixx(i,j,kFixx_q    ) = Fioo(i,j,kFiod_q     )
      end do
      end do

      return
      end

c=======================================================================
