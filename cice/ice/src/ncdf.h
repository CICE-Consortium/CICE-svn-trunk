c $Id: ncdf.h,v 1.4 1998/04/14 14:59:56 eclare Exp $
c.. declarations and parameters for netCDF 
c-----------------------------------------------------------------------
      integer ncID                             ! file ID
      character (80) :: loc_fn                 ! netCDF file name

      integer STATUS                           ! error flag

      integer nframes                          ! time parameter
      integer jmtID, imtID, nframesID          ! dimension IDs
      integer tID_compact                      ! time ID in compact.nc
      integer tID_thick                        ! time ID in thick.nc
      integer tID_vel                          ! time ID in vel.nc
      integer tID_temps                        ! time ID in temps.nc
      integer tID_fluxes                       ! time ID in fluxes.nc

c.. compactness, thickness and velocity fields
      integer hRANK                            ! rank of field
      parameter (hRANK=3)                      ! 2D + time
      integer hSHP(hRANK)                      ! shape of field

      real*4 cmpct(imt,jmt)                    ! compact field
      integer cmpctID                          ! compact ID
      real*4 hithck(imt,jmt)                   ! hithick field
      integer hithckID                         ! hithick ID
      real*4 hithn(imt,jmt)                    ! hithin field
      integer hithnID                          ! hithin ID
      real*4 hsthck(imt,jmt)                   ! hsthick field
      integer hsthckID                         ! hsthick ID
      real*4 hsthn(imt,jmt)                    ! hsthin field 
      integer hsthnID                          ! hsthin ID
      real*4 uice(imt,jmt)                     ! u field
      integer uiceID                           ! u ID
      real*4 vice(imt,jmt)                     ! v field
      integer viceID                           ! v ID

c.. temperatures and heat fluxes
      real*4 tthin(imt,jmt)                    ! thin ice sfc temp field
      integer tthinID                          ! thin ice sfc temp ID
      real*4 tthick(imt,jmt)                   ! thick ice sfc temp field
      integer tthickID                         ! thick ice sfc temp ID
      real*4 flat(imt,jmt)                     ! latent heat flux field
      integer flatID                           ! latent heat flux ID
      real*4 fsen(imt,jmt)                     ! sensible heat flux field
      integer fsenID                           ! sensible heat flux ID
      real*4 flwo(imt,jmt)                     ! outgoing longwave field
      integer flwoID                           ! outgoing longwave ID
      real*4 fnet(imt,jmt)                     ! net flux to ocean field
      integer fnetID                           ! net flux to ocean ID

c..   variable IDs and dimensions
      common /ncdf/ imtID, jmtID, nframesID, 
     & tID_compact, tID_thick, tID_vel, tID_temps, tID_fluxes,
     & cmpctID, hithckID, hithnID, hsthckID, hsthnID, 
     & uiceID, viceID, tthinID, tthickID,
     & flatID, fsenID, flwoID, fnetID

