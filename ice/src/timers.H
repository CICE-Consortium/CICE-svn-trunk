c $Id: timers.H,v 2.0 1999/12/09 22:25:36 eclare Exp $
c-----------------------------------------------------------------------
      integer nb_chronos
      parameter ( nb_chronos = 8 )
      double precision tchrono(0:nb_chronos), uchrono(0:nb_chronos)
      character*8 timern(0:nb_chronos)
      data timern /'Total','TimeLoop','Dynamics','Thermo','Advectn',
     &  'Cat Conv','Stress','Stepu','Bound' /
      common /time_chrono/ tchrono,uchrono
