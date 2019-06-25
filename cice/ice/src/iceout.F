c $Id: iceout.F,v 1.5 1998/04/14 14:59:55 eclare Exp $

c     program iceout.F
c   output files: netCDF data, Fortran unformatted dumps
c-------------------------------------------------------------------------
c.. these subroutines write/read netCDF data files ..
c-------------------------------------------------------------------------
      subroutine netcdf_setup(rec)
c.. writes land masks and initial data in netCDF format
c..       mask_file     masks
c..       compact_file  compactness
c..       thick_file    thicknesses (ice/snow and thin/thick)
c..       vel_file      velocity
      implicit none
#include "ice.h"
#include "netcdf.inc"
#include "ncdf.h"
      integer i, j
      integer rec                              ! record number

c.......... masks ............
      real*4  nc_uvm(imt,jmt)                  ! U-cell mask
      real*4  nc_hm (imt,jmt)                  ! T-cell mask
      integer mRANK                            ! rank of field
      parameter (mRANK=2)                      ! 2D
      integer mSHP(mRANK)                      ! shape of field
      integer uvmID                            ! uvm ID
      integer hmID                             ! hm ID

c.. file name
      loc_fn = mask_file 

c.. netCDF
      ! create file, enter define mode, get ID 
      print*,' Creating file ',loc_fn
      STATUS = NF_CREATE(loc_fn,NF_CLOBBER,ncID) 

      ! define dimensions
      STATUS = NF_DEF_DIM(ncID, 'x', imt, imtID)
      STATUS = NF_DEF_DIM(ncID, 'y', jmt, jmtID)

      ! define variables
      mSHP(1) = imtID
      mSHP(2) = jmtID
      STATUS = NF_DEF_VAR(ncID, 'uvm', NF_REAL, mRANK, mSHP, uvmID)
      STATUS = NF_DEF_VAR(ncID,  'hm', NF_REAL, mRANK, mSHP,  hmID)

      ! assign attributes
      STATUS = NF_PUT_ATT_TEXT(ncID, NF_GLOBAL, 'title', 10,     
     &   'Land Masks')
      STATUS = NF_PUT_ATT_TEXT(ncID, uvmID, 'long_name', 11,        
     &   'U-cell mask')
      STATUS = NF_PUT_ATT_TEXT(ncID, hmID, 'long_name', 11,        
     &   'T-cell mask')

      ! leave define mode
      STATUS = NF_ENDDEF(ncID)           

      ! store
      do j=1,jmt
       do i=1,imt
        nc_uvm(i,j) = uvm(i,j)       ! this is necessary because the
        nc_hm (i,j) =  hm(i,j)       ! ghost points aren't written
       enddo
      enddo

      STATUS = NF_PUT_VAR_REAL (ncID, uvmID, nc_uvm) 
          if (STATUS .ne. NF_NOERR) print*,'error writing uvm'
      STATUS = NF_PUT_VAR_REAL (ncID,  hmID, nc_hm) 
          if (STATUS .ne. NF_NOERR) print*,'error writing hm'

      ! close: save netCDF dataset
      STATUS = NF_CLOSE(ncID)
      if (STATUS .ne. NF_NOERR) print*,'error closing mask_file'

c......... compactness ..........
      loc_fn = compact_file 

c.. netCDF
      ! create file, enter define mode, get ID 
      print*,' Creating file ',loc_fn
      STATUS = NF_CREATE(loc_fn,NF_CLOBBER,ncID) 

      ! define dimensions
      STATUS = NF_DEF_DIM(ncID, 'x', imt, imtID)
      STATUS = NF_DEF_DIM(ncID, 'y', jmt, jmtID)
      STATUS = NF_DEF_DIM(ncID, 'time',NF_UNLIMITED, nframesID)

      ! define variables
      STATUS = NF_DEF_VAR(ncID,'time', NF_INT, 0, 0, tID_compact)

      hSHP(1) = imtID
      hSHP(2) = jmtID
      hSHP(3) = nframesID
      STATUS = NF_DEF_VAR(ncID, 'compact', NF_REAL, hRANK, hSHP, 
     &                           cmpctID)

      ! assign attributes
      STATUS = NF_PUT_ATT_TEXT(ncID, NF_GLOBAL, 'title', 15,     
     &   'Ice Compactness')
      STATUS = NF_PUT_ATT_TEXT(ncID, cmpctID, 'units', 4, 'none')    ! compact
      STATUS = NF_PUT_ATT_TEXT(ncID, cmpctID, 'long_name', 15,       
     &   'Ice Compactness')
      STATUS = NF_PUT_ATT_TEXT(ncID, tID_compact, 'units', 13,       ! time
     &   'record number') 

      ! leave define mode
      STATUS = NF_ENDDEF(ncID)           

      ! close: save netCDF dataset
      STATUS = NF_CLOSE(ncID)
      if (STATUS .ne. NF_NOERR) print*,'error closing compact_file'

c......... thickness ..........
      loc_fn = thick_file

c.. netCDF
      ! create file, enter define mode, get ID 
      print*,' Creating file ',loc_fn
      STATUS = NF_CREATE(loc_fn,NF_CLOBBER,ncID) 

      ! define dimensions
      STATUS = NF_DEF_DIM(ncID, 'x', imt, imtID)
      STATUS = NF_DEF_DIM(ncID, 'y', jmt, jmtID)
      STATUS = NF_DEF_DIM(ncID, 'time',NF_UNLIMITED, nframesID)

      ! define variables
      STATUS = NF_DEF_VAR(ncID,'time', NF_INT, 0, 0, tID_thick)
      STATUS = NF_DEF_VAR(ncID, 'hithick', NF_REAL, hRANK, hSHP, 
     &                           hithckID)
      STATUS = NF_DEF_VAR(ncID, 'hithin', NF_REAL, hRANK, hSHP, 
     &                           hithnID)
      STATUS = NF_DEF_VAR(ncID, 'hsthick', NF_REAL, hRANK, hSHP, 
     &                           hsthckID)
      STATUS = NF_DEF_VAR(ncID, 'hsthin', NF_REAL, hRANK, hSHP, 
     &                           hsthnID)

      ! assign attributes
      STATUS = NF_PUT_ATT_TEXT(ncID, NF_GLOBAL, 'title', 20,     
     &   'Ice Thickness Fields')
      STATUS = NF_PUT_ATT_TEXT(ncID, hithckID, 'units', 2, 'cm') 
      STATUS = NF_PUT_ATT_TEXT(ncID, hithnID,  'units', 2, 'cm') 
      STATUS = NF_PUT_ATT_TEXT(ncID, hsthckID, 'units', 2, 'cm') 
      STATUS = NF_PUT_ATT_TEXT(ncID, hsthnID,  'units', 2, 'cm') 
      STATUS = NF_PUT_ATT_TEXT(ncID, hithckID, 'long_name', 9,       
     &   'Thick Ice')
      STATUS = NF_PUT_ATT_TEXT(ncID, hithnID,  'long_name', 8,       
     &   'Thin Ice')
      STATUS = NF_PUT_ATT_TEXT(ncID, hsthckID, 'long_name', 17,       
     &   'Snow on Thick Ice')
      STATUS = NF_PUT_ATT_TEXT(ncID, hsthnID,  'long_name', 16,       
     &   'Snow on Thin Ice')
      STATUS = NF_PUT_ATT_TEXT(ncID, tID_thick, 'units', 13,       ! time
     &   'record number') 

      ! leave define mode
      STATUS = NF_ENDDEF(ncID)           

      ! close: save netCDF dataset
      STATUS = NF_CLOSE(ncID)
      if (STATUS .ne. NF_NOERR) print*,'error closing thick_file'

c......... velocity ..........
      loc_fn = vel_file

c.. netCDF
      ! create file, enter define mode, get ID 
      print*,' Creating file ',loc_fn
      STATUS = NF_CREATE(loc_fn,NF_CLOBBER,ncID) 

      ! define dimensions
      STATUS = NF_DEF_DIM(ncID, 'x', imt, imtID)
      STATUS = NF_DEF_DIM(ncID, 'y', jmt, jmtID)
      STATUS = NF_DEF_DIM(ncID, 'time',NF_UNLIMITED, nframesID)

      ! define variables
      STATUS = NF_DEF_VAR(ncID,'time', NF_INT, 0, 0, tID_vel)
      STATUS = NF_DEF_VAR(ncID, 'u', NF_REAL, hRANK, hSHP, 
     &                           uiceID)
      STATUS = NF_DEF_VAR(ncID, 'v', NF_REAL, hRANK, hSHP, 
     &                           viceID)

      ! assign attributes
      STATUS = NF_PUT_ATT_TEXT(ncID, NF_GLOBAL, 'title', 15,     
     &   'Velocity Fields')
      STATUS = NF_PUT_ATT_TEXT(ncID, uiceID, 'units', 4, 'cm/s') 
      STATUS = NF_PUT_ATT_TEXT(ncID, viceID, 'units', 4, 'cm/s') 
      STATUS = NF_PUT_ATT_TEXT(ncID, uiceID, 'long_name', 24,       
     &   'Ice velocity component U')
      STATUS = NF_PUT_ATT_TEXT(ncID, viceID, 'long_name', 24,       
     &   'Ice velocity component V')
      STATUS = NF_PUT_ATT_TEXT(ncID, tID_vel, 'units', 13,           ! time
     &   'record number') 

      ! leave define mode
      STATUS = NF_ENDDEF(ncID)           

      ! close: save netCDF dataset
      STATUS = NF_CLOSE(ncID)
      if (STATUS .ne. NF_NOERR) print*,'error closing vel_file'

c.. write initial data
      rec = 1
      call writerecs(rec)

      return
      end
c-------------------------------------------------------------------------
      subroutine writerecs(rec)
c.. writes compactness, thicknesses, velocity in netCDF format
      implicit none
#include "ice.h"
#include "netcdf.inc"
#include "ncdf.h"
      integer rec                              ! record number
      integer i, j                             ! indices
      integer start(hRANK)                     ! for array sections
      integer count(hRANK)

      start(1) = 1                       ! start at first index
      start(2) = 1
      start(3) = rec                     ! write current record
      count(1) = imt                     ! write entire spatial array
      count(2) = jmt
      count(3) = 1                       ! write current record

c....... compactness ...........

      do j=1,jmt
       do i=1,imt
        cmpct(i,j) = compact(i,j)        ! no ghost points
       enddo
      enddo

      ! open file
      STATUS = NF_OPEN(compact_file, NF_WRITE, ncID)

      ! store
c.. time
      STATUS = NF_PUT_VAR_INT(ncID, tID_compact, rec)
c.. compact
      STATUS = NF_PUT_VARA_REAL (ncID, cmpctID, start, count, cmpct) 

      ! close: save netCDF dataset
      STATUS = NF_CLOSE(ncID)

c....... thickness ...........

      do j=1,jmt
       do i=1,imt
        hithck(i,j) = hithick(i,j)       ! no ghost points
        hithn (i,j) = hithin (i,j)
        hsthck(i,j) = hsthick(i,j)
        hsthn (i,j) = hsthin (i,j)
       enddo
      enddo

      ! open file
      STATUS = NF_OPEN(thick_file, NF_WRITE, ncID)

      ! store
c.. time
      STATUS = NF_PUT_VAR_INT(ncID, tID_thick, rec)
c.. thicknesses
      STATUS = NF_PUT_VARA_REAL (ncID, hithckID, start, count, hithck) 
      STATUS = NF_PUT_VARA_REAL (ncID, hithnID,  start, count, hithn) 
      STATUS = NF_PUT_VARA_REAL (ncID, hsthckID, start, count, hsthck) 
      STATUS = NF_PUT_VARA_REAL (ncID, hsthnID,  start, count, hsthn) 

      ! close: save netCDF dataset
      STATUS = NF_CLOSE(ncID)

c....... velocity ...........

      do j=1,jmt
       do i=1,imt
        uice(i,j) = u(i,j)                    ! no ghost points
        vice(i,j) = v(i,j)                    ! no ghost points
       enddo
      enddo

      ! open file
      STATUS = NF_OPEN(vel_file, NF_WRITE, ncID)

      ! store
c.. time
      STATUS = NF_PUT_VAR_INT(ncID, tID_vel, rec)
c.. velocity
      STATUS = NF_PUT_VARA_REAL (ncID, uiceID, start, count, uice) 
      STATUS = NF_PUT_VARA_REAL (ncID, viceID, start, count, vice) 

      ! close: save netCDF dataset
      STATUS = NF_CLOSE(ncID)

      return
      end
c-------------------------------------------------------------------------
      subroutine netcdf_fluxes_setup(rec)
c.. writes additional data in netCDF format
c..       temps_file      temperatures
c..       fluxes_file     fluxes over/under thin ice
      implicit none
#include "ice.h"
#include "netcdf.inc"
#include "ncdf.h"
      integer i, j
      integer rec                              ! record number

c......... temperatures ..........
      loc_fn = temps_file

c.. netCDF
      ! create file, enter define mode, get ID 
      print*,' Creating file ',loc_fn
      STATUS = NF_CREATE(loc_fn,NF_CLOBBER,ncID) 

      ! define dimensions
      STATUS = NF_DEF_DIM(ncID, 'x', imt, imtID)
      STATUS = NF_DEF_DIM(ncID, 'y', jmt, jmtID)
      STATUS = NF_DEF_DIM(ncID, 'time',NF_UNLIMITED, nframesID)

      ! define variables
      hSHP(1) = imtID
      hSHP(2) = jmtID
      hSHP(3) = nframesID
      STATUS = NF_DEF_VAR(ncID, 'time', NF_INT, 0, 0, tID_temps)
      STATUS = NF_DEF_VAR(ncID, 'tthick', NF_REAL, hRANK, hSHP, 
     &                           tthickID)
      STATUS = NF_DEF_VAR(ncID, 'tthin', NF_REAL, hRANK, hSHP, 
     &                           tthinID)

      ! assign attributes
      STATUS = NF_PUT_ATT_TEXT(ncID, NF_GLOBAL, 'title', 22,     
     &   'Ice Temperature Fields')
      STATUS = NF_PUT_ATT_TEXT(ncID, tthickID, 'units', 5, 'deg C') 
      STATUS = NF_PUT_ATT_TEXT(ncID, tthinID,  'units', 5, 'deg C') 
      STATUS = NF_PUT_ATT_TEXT(ncID, tthickID, 'long_name', 30,       
     &   'Surface Temperature, Thick Ice')
      STATUS = NF_PUT_ATT_TEXT(ncID, tthinID,  'long_name', 29,       
     &   'Surface Temperature, Thin Ice')
      STATUS = NF_PUT_ATT_TEXT(ncID, tID_temps, 'units', 13,       ! time
     &   'record number') 

      ! leave define mode
      STATUS = NF_ENDDEF(ncID)           

      ! close: save netCDF dataset
      STATUS = NF_CLOSE(ncID)
      if (STATUS .ne. NF_NOERR) print*,'error closing temps_file'

c......... fluxes ..........
      loc_fn = fluxes_file

c.. netCDF
      ! create file, enter define mode, get ID 
      print*,' Creating file ',loc_fn
      STATUS = NF_CREATE(loc_fn,NF_CLOBBER,ncID) 

      ! define dimensions
      STATUS = NF_DEF_DIM(ncID, 'x', imt, imtID)
      STATUS = NF_DEF_DIM(ncID, 'y', jmt, jmtID)
      STATUS = NF_DEF_DIM(ncID, 'time',NF_UNLIMITED, nframesID)

      ! define variables
      STATUS = NF_DEF_VAR(ncID, 'time', NF_INT, 0, 0, tID_fluxes)
      STATUS = NF_DEF_VAR(ncID, 'flat', NF_REAL, hRANK, hSHP, 
     &                           flatID)
      STATUS = NF_DEF_VAR(ncID, 'fsen', NF_REAL, hRANK, hSHP, 
     &                           fsenID)
      STATUS = NF_DEF_VAR(ncID, 'flwo', NF_REAL, hRANK, hSHP, 
     &                           flwoID)
      STATUS = NF_DEF_VAR(ncID, 'fnet', NF_REAL, hRANK, hSHP, 
     &                           fnetID)

      ! assign attributes
      STATUS = NF_PUT_ATT_TEXT(ncID, NF_GLOBAL, 'title', 28,     
     &   'Heat Fluxes through Thin Ice')
      STATUS = NF_PUT_ATT_TEXT(ncID, flatID, 'units', 5, 'W/m^2') 
      STATUS = NF_PUT_ATT_TEXT(ncID, fsenID, 'units', 5, 'W/m^2') 
      STATUS = NF_PUT_ATT_TEXT(ncID, flwoID, 'units', 5, 'W/m^2') 
      STATUS = NF_PUT_ATT_TEXT(ncID, fnetID, 'units', 5, 'W/m^2') 

      STATUS = NF_PUT_ATT_TEXT(ncID, flatID, 'long_name', 16,       
     &   'Latent Heat Flux')
      STATUS = NF_PUT_ATT_TEXT(ncID, fsenID, 'long_name', 18,       
     &   'Sensible Heat Flux')
      STATUS = NF_PUT_ATT_TEXT(ncID, flwoID, 'long_name', 17,       
     &   'Outgoing Longwave')
      STATUS = NF_PUT_ATT_TEXT(ncID, fnetID, 'long_name', 24,       
     &   'Net Heat Flux Into Ocean')
      STATUS = NF_PUT_ATT_TEXT(ncID, tID_fluxes, 'units', 13,       ! time
     &   'record number') 

      ! leave define mode
      STATUS = NF_ENDDEF(ncID)           

      ! close: save netCDF dataset
      STATUS = NF_CLOSE(ncID)
      if (STATUS .ne. NF_NOERR) print*,'error closing fluxes_file'

c.. write initial data
      rec = 1
      call writerecs_fluxes(rec)

      return
      end
c-------------------------------------------------------------------------
      subroutine writerecs_fluxes(rec)
c.. writes ice temperatures, fluxes in netCDF format
      implicit none
#include "ice.h"
#include "netcdf.inc"
#include "ncdf.h"
      integer rec                              ! record number
      integer i, j                             ! indices
      integer start(hRANK)                     ! for array sections
      integer count(hRANK)

      start(1) = 1                       ! start at first index
      start(2) = 1
      start(3) = rec                     ! write current record
      count(1) = imt                     ! write entire spatial array
      count(2) = jmt
      count(3) = 1                       ! write current record

c....... temperatures ...........

      do j=1,jmt
       do i=1,imt
        tthin(i,j)  = tsfcij(i,j,0)       ! no ghost points
        tthick(i,j) = tsfcij(i,j,1)
       enddo
      enddo

      ! open file
      STATUS = NF_OPEN(temps_file, NF_WRITE, ncID)

      ! store
c.. time
      STATUS = NF_PUT_VAR_INT(ncID, tID_temps, rec)
c.. temperatures
      STATUS = NF_PUT_VARA_REAL (ncID, tthinID,  start, count, tthin) 
      STATUS = NF_PUT_VARA_REAL (ncID, tthickID, start, count, tthick) 

      ! close: save netCDF dataset
      STATUS = NF_CLOSE(ncID)

c....... fluxes over/under thin ice ...........

      do j=1,jmt
       do i=1,imt
        flat(i,j) = 4.1876e4*flatent(i,j,0)      ! W/m^2
        fsen(i,j) = 4.1876e4*fsensible(i,j,0)
        flwo(i,j) = 4.1876e4*flwout(i,j,0)
        fnet(i,j) = 4.1876e4*fhnet(i,j,0)
       enddo
      enddo

      ! open file
      STATUS = NF_OPEN(fluxes_file, NF_WRITE, ncID)

      ! store
c.. time
      STATUS = NF_PUT_VAR_INT(ncID, tID_fluxes, rec)
c.. fluxes
      STATUS = NF_PUT_VARA_REAL (ncID, flatID, start, count, flat) 
      STATUS = NF_PUT_VARA_REAL (ncID, fsenID, start, count, fsen) 
      STATUS = NF_PUT_VARA_REAL (ncID, flwoID, start, count, flwo) 
      STATUS = NF_PUT_VARA_REAL (ncID, fnetID, start, count, fnet) 

      ! close: save netCDF dataset
      STATUS = NF_CLOSE(ncID)

      return
      end
c-----------------------------------------------------------------------
c.. these subroutines write/read Fortran unformatted data files ..
c-----------------------------------------------------------------------
      subroutine dumpfile
c.. Dumps all values needed for a restart every iplot timesteps
      implicit none
#include "ice.h"
      integer i, j, k

      open (50,file=dump_file,form='unformatted')

      write (50) istep1

c   thicknesses
      write (50) ((hithick(i,j),i=1,imt),j=1,jmt)
      write (50) ((hithin(i,j),i=1,imt),j=1,jmt)
      write (50) ((hsthick(i,j),i=1,imt),j=1,jmt)
      write (50) ((hsthin(i,j),i=1,imt),j=1,jmt)
      write (50) ((compact(i,j),i=1,imt),j=1,jmt)

c   thermodynamics
      write (50) ((tsfcij(i,j,0),i=1,imt),j=1,jmt)
      write (50) ((tsfcij(i,j,1),i=1,imt),j=1,jmt)
      write (50) ((nkij(i,j,0),i=1,imt),j=1,jmt)
      write (50) ((nkij(i,j,1),i=1,imt),j=1,jmt)
      write (50) ((qstorij(i,j,0),i=1,imt),j=1,jmt)
      write (50) ((qstorij(i,j,1),i=1,imt),j=1,jmt)
        do k=0,nkmax
      write (50) ((tij(i,j,0,k),i=1,imt),j=1,jmt)
      write (50) ((tij(i,j,1,k),i=1,imt),j=1,jmt)
        enddo
      write (50) ((toceanij(i,j,0),i=1,imt),j=1,jmt)
      write (50) ((toceanij(i,j,1),i=1,imt),j=1,jmt)

c   velocity
      write (50) ((u(i,j),i=1,imt),j=1,jmt)
      write (50) ((v(i,j),i=1,imt),j=1,jmt)

c   stress
      write (50) ((sig11n(i,j),i=1,imt),j=1,jmt)
      write (50) ((sig11e(i,j),i=1,imt),j=1,jmt)
      write (50) ((sig11s(i,j),i=1,imt),j=1,jmt)
      write (50) ((sig11w(i,j),i=1,imt),j=1,jmt)
      write (50) ((sig12n(i,j),i=1,imt),j=1,jmt)
      write (50) ((sig12e(i,j),i=1,imt),j=1,jmt)
      write (50) ((sig12s(i,j),i=1,imt),j=1,jmt)
      write (50) ((sig12w(i,j),i=1,imt),j=1,jmt)
      write (50) ((sig22n(i,j),i=1,imt),j=1,jmt)
      write (50) ((sig22e(i,j),i=1,imt),j=1,jmt)
      write (50) ((sig22s(i,j),i=1,imt),j=1,jmt)
      write (50) ((sig22w(i,j),i=1,imt),j=1,jmt)

      close (50)

      return
      end
c     ------------------------------------------------------
      subroutine restartfile
c.. Restarts from a dump
      implicit none
#include "ice.h"
      integer i, j, k

      open (50,file=restrt_file,form='unformatted')

      read (50) istep0
      print*, 'dump occurred at istep=',istep0

c   thicknesses
      read (50) ((hithick(i,j),i=1,imt),j=1,jmt)
      read (50) ((hithin(i,j),i=1,imt),j=1,jmt)
      read (50) ((hsthick(i,j),i=1,imt),j=1,jmt)
      read (50) ((hsthin(i,j),i=1,imt),j=1,jmt)
      read (50) ((compact(i,j),i=1,imt),j=1,jmt)

c   thermodynamics
      read (50) ((tsfcij(i,j,0),i=1,imt),j=1,jmt)
      read (50) ((tsfcij(i,j,1),i=1,imt),j=1,jmt)
      read (50) ((nkij(i,j,0),i=1,imt),j=1,jmt)
      read (50) ((nkij(i,j,1),i=1,imt),j=1,jmt)
      read (50) ((qstorij(i,j,0),i=1,imt),j=1,jmt)
      read (50) ((qstorij(i,j,1),i=1,imt),j=1,jmt)
        do k=0,nkmax
      read (50) ((tij(i,j,0,k),i=1,imt),j=1,jmt)
      read (50) ((tij(i,j,1,k),i=1,imt),j=1,jmt)
        enddo
      read (50) ((toceanij(i,j,0),i=1,imt),j=1,jmt)
      read (50) ((toceanij(i,j,1),i=1,imt),j=1,jmt)

c   velocity
      read (50) ((u(i,j),i=1,imt),j=1,jmt)
      read (50) ((v(i,j),i=1,imt),j=1,jmt)

c   stress
      read (50) ((sig11n(i,j),i=1,imt),j=1,jmt)
      read (50) ((sig11e(i,j),i=1,imt),j=1,jmt)
      read (50) ((sig11s(i,j),i=1,imt),j=1,jmt)
      read (50) ((sig11w(i,j),i=1,imt),j=1,jmt)
      read (50) ((sig12n(i,j),i=1,imt),j=1,jmt)
      read (50) ((sig12e(i,j),i=1,imt),j=1,jmt)
      read (50) ((sig12s(i,j),i=1,imt),j=1,jmt)
      read (50) ((sig12w(i,j),i=1,imt),j=1,jmt)
      read (50) ((sig22n(i,j),i=1,imt),j=1,jmt)
      read (50) ((sig22e(i,j),i=1,imt),j=1,jmt)
      read (50) ((sig22s(i,j),i=1,imt),j=1,jmt)
      read (50) ((sig22w(i,j),i=1,imt),j=1,jmt)

      close (50)

      return
      end


