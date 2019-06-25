c $Id: history.H,v 2.1 1999/12/10 22:26:52 eclare Exp $

c vars for history file namelist

      integer    avgsiz  ! number of fields to be averaged
      parameter (avgsiz = 23)

      real      avgct    ! average sample counter
      real      aa       ! field accumulations and averages
      logical   iout     ! logical if field should be written

      common/comavg/ avgct,
     &               aa(ilo:ihi,jlo:jhi,avgsiz),
     &               iout(avgsiz)

      real cona(avgsiz)
      real conb(avgsiz)
      character*16 vname(avgsiz)     ! variable names
      character*16 vunit(avgsiz)     ! variable units
      character*40 vdesc(avgsiz)    ! variable descriptions

      data cona / 23*1.0 /    ! note avgsiz
      data conb / 23*0.0 /    ! note avgsiz

      data vname / 
     &     'cmpct', 'hi1',   'hi2',   'hs1',  'hs2',
     &     'uice',  'vice',  'tsfc1', 'tsfc2',
     &     'flat',  'fsen',  'flwo',  'fnet',  
     &     'fslt',  'frsh',  'sig1',  'sig2',
     &     'precp', 'fswdn', 'flwdn', 'sst', 'sss', 'c2' /

      data vdesc / 
     &     'Total Ice Fraction',
     &     'Ice Thickness, cat 1',
     &     'Ice Thickness, cat 2',
     &     'Snow Depth, cat 1',
     &     'Snow Depth, cat 2',
     &     'Ice Velocity U',
     &     'Ice Velocity V',
     &     'Surface Temp, cat 1',
     &     'Surface Temp, cat 2',
     &     'Latent Heat Flux',
     &     'Sensible Heat Flux',
     &     'Outward Longwave',
     &     'Ocean Heat Used',
     &     'Salt Flux',
     &     'Fresh Water Flux',
     &     'Normalized Principal Stress 1',
     &     'Normalized Principal Stress 2',
     &     'Snowfall Rate',
     &     'Downward Shortwave',
     &     'Downward Longwave',
     &     'Sea Surface Temperature',
     &     'Sea Surface Salinity',
     &     'Ice Compactness, cat 2' /

      data vunit / 
     &     ' ',       'm',      'm',        'm',
     &     'm',       'm/s',    'm/s',      'deg C',
     &     'deg C',    'W/m2',   'W/m2',     'W/m2',
     &     'W/m2',    'kg/m^2 s',   'kg/m^2 s',  
     &     ' ',       ' ',
     &     'kg/m^2 s','W/m^2', 'W/m^2', 'deg C', 'ppt', ' '  /

      logical
     &     f_cmpct, f_hi1,  f_hi2,   f_hs1,   f_hs2, 
     &     f_uice,  f_vice, f_tsfc1, f_tsfc2, 
     &     f_flat,  f_fsen, f_flwo,  f_fnet,  
     &     f_salt,  f_fresh,
     &     f_sig1,  f_sig2,
     &     f_snow,  f_fsw,  f_flw,   f_sst,   f_sss, f_c2

      namelist / icefields /
     &     f_cmpct, f_hi1,  f_hi2,   f_hs1,   f_hs2, 
     &     f_uice,  f_vice, f_tsfc1, f_tsfc2, 
     &     f_flat,  f_fsen, f_flwo,  f_fnet,  
     &     f_salt,  f_fresh,
     &     f_sig1,  f_sig2,
     &     f_snow,  f_fsw,  f_flw,   f_sst,   f_sss, f_c2

      equivalence (iout(1),f_cmpct)
      equivalence (iout(2),f_hi1)
      equivalence (iout(3),f_hi2)
      equivalence (iout(4),f_hs1)
      equivalence (iout(5),f_hs2)
      equivalence (iout(6),f_uice)
      equivalence (iout(7),f_vice)
      equivalence (iout(8),f_tsfc1)
      equivalence (iout(9),f_tsfc2)
      equivalence (iout(10),f_flat)
      equivalence (iout(11),f_fsen)
      equivalence (iout(12),f_flwo)
      equivalence (iout(13),f_fnet)
      equivalence (iout(14),f_salt)
      equivalence (iout(15),f_fresh)
      equivalence (iout(16),f_sig1)
      equivalence (iout(17),f_sig2)
      equivalence (iout(18),f_snow)
      equivalence (iout(19),f_fsw)
      equivalence (iout(20),f_flw)
      equivalence (iout(21),f_sst)
      equivalence (iout(22),f_sss)
      equivalence (iout(23),f_c2)












