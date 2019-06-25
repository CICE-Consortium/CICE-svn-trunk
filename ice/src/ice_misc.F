c $Id: ice_misc.F,v 2.1 2000/01/14 23:07:07 eclare Exp $

c--------------------------------------------------------------
c.. miscellaneous subroutines
c--------------------------------------------------------------
      subroutine roughness
c.. computes transfer coefficients for latent and sensible heat
c.. fluxes, wind stress, and the atmospheric reference temperature
      implicit none
      include "ice.h"          ! ice code common blocks
      include "thermw.h"       ! thermodynamics common blocks

      !--- constants --------------------------------------
      real,parameter :: umin   =  1.0            ! minimum wind speed (m/s)
      real,parameter :: gravit =  9.80616        ! Acceleration of gravity
      real,parameter :: xkar   =  0.4            ! Von Karman constant
      real,parameter :: zref   = 10.0            ! reference height (m)
      real,parameter :: zvir   =  0.606          ! rh2o/rair - 1.0
      real,parameter :: ztref  =  2.0            ! atm reference height T (m)
      real,parameter :: cpwv   =  1.81e3         ! Specific heat of water vapor
      real,parameter :: cpvir  = cpwv/cpair - 1. ! Defined as cpwv/cpair - 1.

      !--- local variables --------------------------------
      integer :: i, j, n, k
      real    :: Tsfk   ! surface ice temperature (K)
      real    :: vmag   ! surface wind magnitude, limited (m/s)
      real    :: zvirT  ! zvir * potential temperature (for convenience)(K)
      real    :: thvbot ! virtual temperature      (K)
      real    :: ssq    ! sea surface humidity     (kg/kg)
      real    :: dssqdt ! derivative of ssq wrt Ts (kg/kg/K)
      real    :: delt   ! potential T difference   (K)
      real    :: delq   ! humidity difference      (kg/kg)
      real    :: alz    ! ln(zbot/z10)
      real    :: tau    ! stress at zbot
      real    :: xkargz ! xkar * gravit * zlvl
      real    :: rdn    ! sqrt of neutral exchange coefficient (momentum)
      real    :: rhn    ! sqrt of neutral exchange coefficient (heat)
      real    :: ren    ! sqrt of neutral exchange coefficient (water)
      real    :: hol    ! H (at zlvl) over L
      real    :: stable ! stability factor
      real    :: xsq    ! temporary variable
      real    :: xqq    ! temporary variable
      real    :: psimh  ! stability function at zbot (momentum)
      real    :: psixh  ! stability function at zbot (heat and water)
      real    :: rd     ! sqrt of exchange coefficient (momentum)
      real    :: rh     ! sqrt of exchange coefficient (heat)
      real    :: re     ! sqrt of exchange coefficient (water)            
      real    :: bn     ! exchange coef funct for interpolation
      real    :: bh     ! exchange coef funct for interpolation
      real    :: fac    ! interpolation factor
      real    :: ln0    ! log factor for interpolation
      real    :: ln3    ! log factor for interpolation
      real    :: cp     ! specific heat of moist air
      real    :: ustr   ! ustar
      real    :: tstr   ! tstar
      real    :: qstr   ! qstar

      real ustarn(ilo:ihi,jlo:jhi,ncat)  ! ustar for each category
      real tstarn(ilo:ihi,jlo:jhi,ncat)  ! tstar
      real qstarn(ilo:ihi,jlo:jhi,ncat)  ! qstar

      real ustar(ilo:ihi,jlo:jhi)  ! ustar, merged over categories
      real tstar(ilo:ihi,jlo:jhi)  ! tstar
      real qstar(ilo:ihi,jlo:jhi)  ! qstar

      real Trefn(ilo:ihi,jlo:jhi,ncat)  ! reference temp, ea category

      !--- local functions --------------------------------
      real    :: Tk      ! dummy arg ~ temperature (K)
      real    :: qsat    ! the saturation humididty of air (kg/m^3)
      real    :: dqsatdt ! derivivative of qsat wrt surface temperature
      real    :: xd      ! dummy argument  
      real    :: psimhu  ! unstable part of psimh
      real    :: psixhu  ! unstable part of psimx

      qsat   (Tk) = 640380. / exp(5107.4/Tk)
      dqsatdt(Tk) = (5107.4 / Tk**2) * 640380. / exp(5107.4/Tk)
      psimhu (xd) = log((1.+xd*(2.+xd))*(1.+xd*xd)/8.) 
     &            - 2.*atan(xd) + 1.571
      psixhu (xd) = 2.*log((1. + xd*xd)/2.)

c$OMP PARALLEL DO PRIVATE(i,j,n,k,Tsfk,vmag,zvirT,thvbot,ssq,
c$OMP&   dssqdt,delt,delq,alz,xkargz,cp,rdn,rhn,ren,ustr,tstr,
c$OMP&   qstr,hol,stable,xsq,xqq,psimh,psixh,rd,rh,re,bn,bh,
c$OMP&   ln0,ln3,fac)
      do n=1,ncat
       do j=jlo,jhi
        do i=ilo,ihi

         if (hm(i,j).gt.tiny) then

      !--- define some needed variables ---
      Tsfk   = Tsfc(i,j,n) + Tf0kel         ! ice sfc temp (K)
      vmag   = max(umin, wind(i,j))
      zvirT  = zvir * potT(i,j)
      thvbot = potT(i,j) + zvirT * Qa(i,j)  ! virtual pot temp (K)
      ssq    = qsat(Tsfk) / rhoa(i,j)       ! sea surf hum (kg/kg)
      dssqdt = dqsatdt(Tsfk) / rhoa(i,j)    ! deriv of ssq wrt Ts 
      delt   = potT(i,j) - Tsfk             ! pot temp diff (K)
      delq   = Qa(i,j) - ssq                ! spec hum dif (kg/kg)
      alz    = log(zlvl(i,j)/zref) 
      xkargz = xkar * gravit * zlvl(i,j)
      cp     = cpair*(1. + cpvir*ssq) 

      !------------------------------------------------------------
      !     first estimate of Z/L and ustar, tstar and qstar
      !------------------------------------------------------------

      !--- neutral coefficients, z/L = 0.0 ---
      rdn = xkar/log(zref/iceruf(n))
      rhn = rdn
      ren = rdn

      !--- ustar,tstar,qstar ----
      ustr = rdn * vmag
      tstr = rhn * delt  
      qstr = ren * delq  

      !------------------------------------------------------------
      !     iterate to converge on Z/L, ustar, tstar and qstar
      !------------------------------------------------------------

      do k=1,5

      !--- compute stability & evaluate all stability functions ---
      hol    = xkargz/ustr**2 * (tstr + zvirT*qstr)/thvbot
      hol    = sign( min(abs(hol),10.), hol )
      stable = 0.5 + sign(0.5 , hol)
      xsq    = max(sqrt(abs(1. - 16.*hol)) , 1.)
      xqq    = sqrt(xsq)
      psimh  = -5.*hol*stable + (1.-stable)*psimhu(xqq)
      psixh  = -5.*hol*stable + (1.-stable)*psixhu(xqq)

      !--- shift all coeffs to measurement height and stability ---
      rd = rdn / (1.+rdn/xkar*(alz-psimh))
      rh = rhn / (1.+rhn/xkar*(alz-psixh))
      re = ren / (1.+ren/xkar*(alz-psixh))

      !--- update ustar, tstar, qstar using updated, shifted coeffs --
      ustr = rd * vmag 
      tstr = rh * delt 
      qstr = re * delq 

      enddo    ! end iteration

      !--- save arrays
      ustarn(i,j,n) = ustr
      tstarn(i,j,n) = tstr
      qstarn(i,j,n) = qstr

      !--- prepare coefficients for turbulent flux calculation
      shcoef(i,j,n) = rhoa(i,j)*wind(i,j)*cp           *rd*rh
      lhcoef(i,j,n) = rhoa(i,j)*wind(i,j)*(Lvap+Lfresh)*rd*re

      !--- prepare factor for reference temperature interpolation.
      !--- Assume that 
      !--- cn = rdn*rdn, cm=rd*rd and ch=rh*rd, and therefore 
      !--- 1/sqrt(cn(i,j))=1/rdn and sqrt(cm(i,j))/ch(i,j)=1/rh 
      bn = xkar/rdn
      bh = xkar/rh

      !--- Interpolation factor for stable and unstable cases
      ln0 = log(1.0 + (ztref/zlvl(i,j))*(exp(bn) - 1.0))
      ln3 = log(1.0 + (ztref/zlvl(i,j))*(exp(bn - bh) - 1.0))
      fac = ( ln0 - ztref/zlvl(i,j)*(bn-bh) )/bh*stable 
     &    + (ln0-ln3)/bh*(1.-stable)
      fac = min(max(fac,0.),1.)
      Trefn(i,j,n) = Tsfk + (Tair(i,j) - Tsfk)*fac

         else  ! land points

      ustarn(i,j,n) = 0.
      tstarn(i,j,n) = 0.
      qstarn(i,j,n) = 0.
       Trefn(i,j,n) = 0.
 
         endif
        enddo
       enddo
      enddo

      !------------------------------------------------------------
      !             merge coefficients across categories
      !------------------------------------------------------------

      call merge_state(ustarn,ustar)  ! does not include category 0,
c      call merge_state(tstarn,tstar)  ! open water
c      call merge_state(qstarn,qstar)
      call merge_state( Trefn,Tref )

cc$OMP PARALLEL DO PRIVATE(i,j,n)     ! includes open water cat 0
c      do j=jlo,jhi
c       do i=ilo,ihi
c        ustar(i,j) = 0.
c        tstar(i,j) = 0.
c        qstar(i,j) = 0.
c          fac(i,j) = 0.
c         do n = 0,ncat
c          ustar(i,j) = ustar(i,j) + compact(i,j,n)*ustarn(i,j,n)
c          tstar(i,j) = tstar(i,j) + compact(i,j,n)*tstarn(i,j,n)
c          qstar(i,j) = qstar(i,j) + compact(i,j,n)*qstarn(i,j,n)
c           Tref(i,j) =  Tref(i,j) + compact(i,j,n)* Trefn(i,j,n)
c         enddo
c       enddo
c      enddo

      !------------------------------------------------------------
      !         compute the momentum flux
      !------------------------------------------------------------

c$OMP PARALLEL DO PRIVATE(i,j,tau,vmag)
      do j=jlo,jhi
       do i=ilo,ihi
         if (hm(i,j).gt.tiny) then
      tau = rhoa(i,j) * ustar(i,j)**2
      vmag = max(umin,wind(i,j))
      strairx(i,j) = tau * uatm(i,j) / vmag  ! kg/m s^2
      strairy(i,j) = tau * vatm(i,j) / vmag
         else  ! land points
      strairx(i,j) = 0.
      strairy(i,j) = 0.
         endif
       enddo
      enddo

      call t2ugrid(strairx)
      call t2ugrid(strairy)
      call t2ugrid(wind) 

      return
      end
c--------------------------------------------------------------
      subroutine calendar(ttime)
c.. determine the date at the end of the time step
      implicit none
      include "ice.h"          ! ice code common blocks

      real ttime                                ! time variable
      integer k
      integer nyrp,mdayp,weekp                  ! previous year, day, week
      real tday                                 ! absolute day number
      real dayyr                                ! number of days per year

      dayyr = secyr/86400. 

      nyrp=nyr
      monthp=month
      weekp=week
      mdayp=mday
      new_year=.false.
      new_month=.false.
      new_week=.false.
      new_day=.false.
      write_history=.false.
      write_restart=.false.
      update_forcing=.false.

      sec = mod(ttime,86400.)           ! elapsed seconds into date at 
                                        ! end of dt
      tday = (ttime-sec)/86400. + 1.    ! absolute day number 
      yday = mod(tday-1.,dayyr) + 1.    ! day of the year 
      week = int(yday/7.)               ! week of the year
      do k = 1, 12
        if (yday .gt. float(daycal(k))) month = k ! month
      enddo
      mday = int(yday) - daycal(month)       ! day of the month
      nyr = year + int((tday-1.)/dayyr)      ! year number - 1

      idate = (nyr)*10000 + month*100 + mday ! date (yyyymmdd) 

      if (my_task.eq.master_task.and.mod(istep,diagfreq).eq.0) then
        write(6,*) ' '
        write(6,*) 'istep1:',istep1,' date: ',idate,' sec: ',sec
      endif

      if (nyr.ne.nyrp) new_year = .true.
      if (month.ne.monthp) new_month = .true.
      if (week.ne.weekp) new_week = .true.
      if (mday.ne.mdayp) new_day = .true.

      if (((outfreq.eq.'y'.or.outfreq.eq.'Y').and.new_year).or. 
     &      ((outfreq.eq.'m'.or.outfreq.eq.'M').and.new_month).or.
     &      ((outfreq.eq.'w'.or.outfreq.eq.'W').and.new_week).or.
     &      ((outfreq.eq.'d'.or.outfreq.eq.'D').and.new_day).or.
     &       (outfreq.eq.'1'))
     &      write_history=.true.
      if (((dumpfreq.eq.'y'.or.dumpfreq.eq.'Y').and.new_year).or. 
     &      ((dumpfreq.eq.'m'.or.dumpfreq.eq.'M').and.new_month).or.
     &      ((dumpfreq.eq.'d'.or.dumpfreq.eq.'D').and.new_day  ))
c     &        .and.mod(mday,5).eq.0))
     &      write_restart=.true.
      if (new_day) update_forcing=.true.

      return
      end
c-----------------------------------------------------------------------
      subroutine ice_global_real_minmax(nc,work,string)
      implicit none
      include "ice.h"
      integer nc
      real work(nc),amin,amax
      character*8 string
      amin = ice_global_real_minval(nc,work)
      amax = ice_global_real_maxval(nc,work)
      if (my_task.eq.master_task) 
     &    write (6,*) ' min/max ',string,amin,amax
      return
      end
c--------------------------------------------------------------
      function ice_global_real_minval(nc,work)
c.. computes minimum over the global grid
      implicit none
      include "ice.h"          ! ice code common blocks
      include "ice_mpi.h"      ! common blocks for mpi
#ifdef _MPI
      include "mpif.h"         ! MPI library definitions
#endif
      integer nc,n
      real work(nc),local_val

      local_val = minval(work)
#ifdef _MPI
      call MPI_ALLREDUCE(local_val, ice_global_real_minval, 1,
     &                   MPI_REAL8, MPI_MIN, MPI_COMM_ICE, ierr)
#else
      ice_global_real_minval=local_val
#endif
      return
      end
c--------------------------------------------------------------
      function ice_global_real_maxval(nc,work)
c.. computes maximum over the global grid
      implicit none
      include "ice.h"          ! ice code common blocks
      include "ice_mpi.h"      ! common blocks for mpi
#ifdef _MPI
      include "mpif.h"         ! MPI library definitions
#endif
      integer nc,n
      real work(nc),local_val

      local_val = maxval(work)
#ifdef _MPI
      call MPI_ALLREDUCE(local_val, ice_global_real_maxval, 1,
     &                   MPI_REAL8, MPI_MAX, MPI_COMM_ICE, ierr)
#else
      ice_global_real_maxval=local_val
#endif
      return
      end
c--------------------------------------------------------------
      function ice_global_real_sum(nc,work)
c.. sums given array over the global grid
      implicit none
      include "ice.h"          ! ice code common blocks
      include "ice_mpi.h"      ! common blocks for mpi
#ifdef _MPI
      include "mpif.h"         ! MPI library definitions
#endif
      integer nc,n,i1,i2,j1,j2
      real work(nc),local_val

      local_val = sum(work)
#ifdef _MPI
      call MPI_ALLREDUCE(local_val, ice_global_real_sum, 1,
     &                   MPI_REAL8, MPI_SUM, MPI_COMM_ICE, ierr)
#else
      ice_global_real_sum=local_val
#endif
      return
      end
c--------------------------------------------------------------
      subroutine ice_bcast_rscalar(val)
      implicit none
      include 'ice.h'
      include "ice_mpi.h"
#ifdef _MPI
      include "mpif.h"              ! MPI specific definitions
#endif

      real val
#ifdef _MPI
      call MPI_BCAST(val,1,MPI_REAL8,master_task,MPI_COMM_ICE,ierr)
#endif
      return
      end
c--------------------------------------------------------------
      subroutine ice_bcast_iscalar(ival)
      implicit none
      include 'ice.h'
      include "ice_mpi.h"
#ifdef _MPI
      include "mpif.h"              ! MPI specific definitions
#endif

      integer ival
#ifdef _MPI
      call MPI_BCAST(ival,1,MPI_INTEGER,master_task,MPI_COMM_ICE,ierr)
#endif
      return
      end
c--------------------------------------------------------------
      subroutine ice_open(nu,filename,nbits)
c.. opens an unformatted file for reading
c.. nbits indicates whether the file is sequential or direct access
      implicit none
      include "ice.h"          ! ice code common blocks
      include "ice_mpi.h"      ! common blocks for mpi
      integer nu,nbits
      character (*) :: filename
      if (my_task.eq.master_task) then
      if (nbits.eq.0) then
        open(nu,file=filename,form='unformatted')
      else
        open(nu,file=filename,recl=imt_global*jmt_global*nbits/8,
     &       form='unformatted',access='direct')
      endif
      endif
      return
      end
c--------------------------------------------------------------
      subroutine ice_read(nu,nrec,work,type)
c.. reads an unformatted file
c.. work is a real array, type indicates the format of the data
      implicit none
      include "ice.h"          ! ice code common blocks
      include "ice_mpi.h"      ! common blocks for mpi
#ifdef _MPI
      include "mpif.h"         ! MPI library definitions
#endif
      integer nu,nrec,n,i,j
      character*4 type     ! real,int , 
      real*8    workg (imt_global,jmt_global)
      real*4    workg4(imt_global,jmt_global)
      integer*4 workgi(imt_global,jmt_global)
      real work(ilo:ihi,jlo:jhi)
      real amin,amax

      if (my_task.eq.master_task) then
        if (type.eq.'ida4') then
          read(nu,rec=nrec) workgi
          workg=real(workgi)
        elseif (type.eq.'rda8') then
          read(nu,rec=nrec) workg
        elseif (type.eq.'rda4') then
          read(nu,rec=nrec) workg4
          workg = workg4
        elseif (type.eq.'ruf8') then
          read(nu) ((workg(i,j),i=1,imt_global),j=1,jmt_global)
        else
          write (6,*) ' ERROR: reading unknown type ',type
        endif
        amin=workg(1,1)
        amax=workg(1,1)
        do j=1,jmt_global
        do i=1,imt_global
          amin=min(amin,workg(i,j))
          amax=max(amax,workg(i,j))
        enddo
        enddo
        write (6,*) ' read_global ',nu,nrec,amin,amax
      endif

      call global_scatter(workg,work)

      return
      end
c--------------------------------------------------------------
      subroutine ice_write(nu,nrec,work,type)
c.. writes an unformatted file
c.. work is a real array, type indicates the format of the data
      implicit none
      include "ice.h"          ! ice code common blocks
      include "ice_mpi.h"      ! common blocks for mpi
#ifdef _MPI
      include "mpif.h"         ! MPI library definitions
#endif
      integer nu,nrec,n,i,j
      character*4 type     ! real,int , 
      real*8    workg (imt_global,jmt_global)
      integer*4 workgi(imt_global,jmt_global)
      real work(ilo:ihi,jlo:jhi)
      real amin,amax

      call global_gather(workg,work)

      if (my_task.eq.master_task) then
        if (type.eq.'ida4') then
          workgi=nint(workg)
          write(nu,rec=nrec) workgi
        elseif (type.eq.'rda8') then
          write(nu,rec=nrec) workg
        elseif (type.eq.'ruf8') then
          write(nu) ((workg(i,j),i=1,imt_global),j=1,jmt_global)
        else
          write (6,*) ' ERROR: writing unknown type ',type
        endif
      endif

      return
      end
c--------------------------------------------------------------
      subroutine end_run
      implicit none
      include "ice.h"          ! ice code common blocks
      include "ice_mpi.h"      ! common blocks for mpi
#ifdef _MPI
      include "mpif.h"         ! MPI library definitions

      call MPI_FINALIZE(ierr)
#endif

      return
      end
c-----------------------------------------------------------------------
      subroutine div_u
c.. calculates divergence of U
c.. this routine is provided for your convenience--it is not necessary
      implicit none
      include "ice.h"          ! ice code common blocks
      integer i, j
      real factor, UTE, UTW, VTN, VTS

c$OMP PARALLEL DO PRIVATE(i,j,factor,UTE,UTW,VTN,VTS)
      do j=jlo,jhi
       do i=ilo,ihi

      factor = 0.5*dxtr(i,j)*dytr(i,j)     ! 1/(2 dx dy)
      UTE =  factor*HTE(i,j)
      UTW = -factor*HTE(i-1,j)
      VTN =  factor*HTN(i,j)
      VTS = -factor*HTN(i,j-1)
      
      divu(i,j) = UTE*(u(i,j)+u(i,j-1)) + UTW*(u(i-1,j) + u(i-1,j-1))
     1         + VTN*(v(i,j)+v(i-1,j)) + VTS*(v(i,j-1) + v(i-1,j-1))

       enddo
      enddo

      return
      end
c-----------------------------------------------------------------------








