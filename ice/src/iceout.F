c $Id: iceout.F,v 2.1 1999/12/10 22:26:53 eclare Exp $

c-------------------------------------------------------------------------
c   output files: netCDF data, Fortran unformatted dumps
c-------------------------------------------------------------------------
c.. these subroutines write netCDF data files ..
c-------------------------------------------------------------------------
      subroutine iceavg
c.. Average ice quantities based upon 'iaction':       
c.. by Bruce P. Briegleb, 31 August 1999 
c.. modified by C.M. Bitz, T. Craig and E. Hunke...
      implicit none
      include "ice.h"
      include "history.h"
      integer k        ! generic counters

      avgct = avgct + 1.0

      if (iout( 1)) aa(:,:, 1) = aa(:,:, 1) + ifrc(ilo:ihi,jlo:jhi)
      if (iout( 2)) aa(:,:, 2) = aa(:,:, 2) + hice(ilo:ihi,jlo:jhi,1)
      if (iout( 3)) aa(:,:, 3) = aa(:,:, 3) + hice(ilo:ihi,jlo:jhi,ncat)
      if (iout( 4)) aa(:,:, 4) = aa(:,:, 4) + hsnow(ilo:ihi,jlo:jhi,1)
      if (iout( 5)) aa(:,:, 5) = aa(:,:, 5) + 
     &   hsnow(ilo:ihi,jlo:jhi,ncat)
      if (iout( 6)) aa(:,:, 6) = aa(:,:, 6) + u(ilo:ihi,jlo:jhi)
      if (iout( 7)) aa(:,:, 7) = aa(:,:, 7) + v(ilo:ihi,jlo:jhi)
      if (iout( 8)) aa(:,:, 8) = aa(:,:, 8) + Tsfc(ilo:ihi,jlo:jhi,1)
      if (iout( 9)) aa(:,:, 9) = aa(:,:, 9) + Tsfc(ilo:ihi,jlo:jhi,ncat)
      if (iout(10)) aa(:,:,10) = aa(:,:,10) + 
     &   flatent(ilo:ihi,jlo:jhi,1) + flatent(ilo:ihi,jlo:jhi,2)
      if (iout(11)) aa(:,:,11) = aa(:,:,11) + 
     &   fsensible(ilo:ihi,jlo:jhi,1) + fsensible(ilo:ihi,jlo:jhi,2)
      if (iout(12)) aa(:,:,12) = aa(:,:,12) + 
     &   flwout(ilo:ihi,jlo:jhi,1) + flwout(ilo:ihi,jlo:jhi,2)
      if (iout(13)) aa(:,:,13) = aa(:,:,13) + 
     &   fhnet(ilo:ihi,jlo:jhi,1) + fhnet(ilo:ihi,jlo:jhi,2)
      if (iout(14)) aa(:,:,14) = aa(:,:,14) +
     &   fsalt(ilo:ihi,jlo:jhi,1) + fsalt(ilo:ihi,jlo:jhi,2)
      if (iout(15)) aa(:,:,15) = aa(:,:,15) +
     &   fresh(ilo:ihi,jlo:jhi,1) + fresh(ilo:ihi,jlo:jhi,2)
      call principal_stress
      if (iout(16)) aa(:,:,16) = aa(:,:,16) + sig1(ilo:ihi,jlo:jhi)
      if (iout(17)) aa(:,:,17) = aa(:,:,17) + sig2(ilo:ihi,jlo:jhi)
      if (iout(18)) aa(:,:,18) = aa(:,:,18) + snow(ilo:ihi,jlo:jhi)
      if (iout(19)) aa(:,:,19) = aa(:,:,19) + fsw(ilo:ihi,jlo:jhi)
      if (iout(20)) aa(:,:,20) = aa(:,:,20) + flw(ilo:ihi,jlo:jhi)
      if (iout(21)) aa(:,:,21) = aa(:,:,21) + sst(ilo:ihi,jlo:jhi)
      if (iout(22)) aa(:,:,22) = aa(:,:,22) + sss(ilo:ihi,jlo:jhi)
      if (iout(23)) aa(:,:,23) = aa(:,:,23) +
     &   compact(ilo:ihi,jlo:jhi,ncat)

      if (write_history .or. write_ic) then
        do k=1,avgsiz
          if (iout(k)) aa(:,:,k) = (cona(k)*aa(:,:,k)/avgct + conb(k))
        enddo
        call icecdf
        aa = 0.0
        avgct = 0.
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine icecdf
      implicit none
      include "ice.h"
      include "history.h"
      include "netcdf.inc"

      integer :: i,j,n,ncid,status,imtid,jmtid,timid,varid
      integer :: length,lenstr
      integer, dimension(3) :: dimid,start,count
      real*4 :: ltime
      real*4,dimension(imt_global) :: x
      real*4,dimension(jmt_global) :: y
      character(len=80) :: ncfile,title
      real*4 :: spval
      parameter (spval=1.0e36)

      real,dimension(ilo:ihi,jlo:jhi) :: work
      real,dimension(imt_global,jmt_global) :: arg
      real*4,dimension(imt_global,jmt_global) :: ar

      if (my_task.eq.master_task) then

        do i=1,imt_global
          x(i)=i
        enddo
        do j=1,jmt_global
          y(j)=j
        enddo

        ltime=time/86400

        if (write_ic) then
          write(ncfile,'(a,i8.8,''.incond.nc'')') 
     &    history_file(1:lenstr(history_file)),idate
        else
          write(ncfile,'(a,i8.8,''.'',i6.6,''.nc'')') 
     &    history_file(1:lenstr(history_file)),idate,sec
        endif

        status = nf_create(ncfile, nf_clobber, ncid)
        write(6,*) ' '
        write(6,*) 'Writing ',ncfile(1:lenstr(ncfile))
        write(6,*) ' '

c define dimensions

        status = nf_def_dim(ncid,'x',imt_global,imtid)
        status = nf_def_dim(ncid,'y',jmt_global,jmtid)
        status = nf_def_dim(ncid,'time',NF_UNLIMITED,timid)

c define coordinate variables

        dimid(1) = imtid
        status = nf_def_var(ncid,'iindex',nf_float,1,dimid,varid)
        status = nf_put_att_text(ncid,varid,'long_name',
     $                           7,'i index')
        status = nf_put_att_text(ncid,varid,'units',
     $                           5,'index')

        dimid(1) = jmtid
        status = nf_def_var(ncid,'jindex',nf_float,1,dimid,varid)
        status = nf_put_att_text(ncid,varid,'long_name',
     $                           7,'j index')
        status = nf_put_att_text(ncid,varid,'units',
     $                           5,'index')

        dimid(1) = timid
        status = nf_def_var(ncid,'time',nf_float,1,dimid,varid)
        status = nf_put_att_text(ncid,varid,'long_name',
     $                           10,'model time')
        status = nf_put_att_text(ncid,varid,'units',
     $                  30,'days since 0000-00-00 00:00:00')  ! for now
        status = nf_put_att_text(ncid,varid,'calendar',
     $                           6,'noleap')
        if (outfreq .eq. '1') then
        status = nf_put_att_text(ncid,varid,'time_rep',
     $                          13,'instantaneous')
        else
        status = nf_put_att_text(ncid,varid,'time_rep',
     $                           8,'averaged')
        endif

c define variables and variable attributes

        dimid(1) = imtid
        dimid(2) = jmtid
        dimid(3) = timid
        do n=1,avgsiz
          if (iout(n)) then
            status  = nf_def_var(ncid, vname(n), nf_float, 
     $                         3, dimid, varid)
            length = lenstr(vunit(n))
            status = nf_put_att_text(ncid,varid,
     $             'units',length,vunit(n))
            length = lenstr(vdesc(n))
            status = nf_put_att_text(ncid,varid,
     $             'long_name',length,vdesc(n))
            status = nf_put_att_text(ncid,varid,'coordinates',
     $                             8,'i j time')
            status = nf_put_att_real(ncid,varid,'missing_value',
     $                             nf_float,1,spval)
          endif
        enddo

c global attributes
c..  the user should change these to something useful ..
        title  = 'sea ice model output'
        length = lenstr(title)
        status = nf_put_att_text(ncid,nf_global,'title',length,title)
        title = 'sea ice model: CICE v2'
        length = lenstr(title)
        status = nf_put_att_text(ncid,nf_global,'source',length,title)
        title = 'log file ice.log.'                      ! for now
        length = lenstr(title)
        status = nf_put_att_text(ncid,nf_global,'history',length,title)
        title = 'of Elizabeth Hunke'
        length = lenstr(title)
        status = 
     &       nf_put_att_text(ncid,nf_global,'conventions',length,title)

c end define mode

        status = nf_enddef(ncid)

c write out coordinate variables

        status = nf_inq_varid(ncid,'iindex',varid)
        status = nf_put_var_real(ncid,varid,x)

        status = nf_inq_varid(ncid,'jindex',varid)
        status = nf_put_var_real(ncid,varid,y)

        status = nf_inq_varid(ncid,'time',varid)
        status = nf_put_vara_real(ncid,varid,1,1,ltime)

      endif          ! master_task

c write out variable data 

      start=1
      count(1)=imt_global
      count(2)=jmt_global
      count(3)=1
      do n=1,avgsiz
        if (iout(n)) then
          work=aa(:,:,n)
          call global_gather(arg,work)
          ar=arg
          if (my_task.eq.master_task) then
            status = nf_inq_varid(ncid,vname(n),varid)
            status = nf_put_vara_real(ncid,varid,start,count,ar)
          endif
        endif
      enddo

c close output dataset

      if (my_task.eq.master_task) then
        status = nf_close(ncid)
      endif

      return
      end
c-----------------------------------------------------------------------
      integer function lenstr(label)

c compute length of string by finding first non-blank
c character from the right.

      implicit none
      character*(*) label
      integer length  ! length of character string 
      integer len     ! intrinsic length function
      integer n       ! loop index

      length = len(label)
      do n=length,1,-1
        if( label(n:n) .ne. ' ' ) goto 10 
      enddo
 10   lenstr = n

      return 
      end

c-----------------------------------------------------------------------
c.. these subroutines write/read Fortran unformatted data files ..
c-----------------------------------------------------------------------
      subroutine dumpfile
c.. Dumps all values needed for a restart 
      implicit none
      include "ice.h"          ! ice code common blocks
      integer i, j, k, n, lenstr
      character(len=80) :: filename

      real work(ilo:ihi,jlo:jhi)

      write(filename,'(a,i8.8,''.'',i6.6)') 
     &  dump_file(1:lenstr(dump_file)),idate,sec
      call ice_open(50,filename,0)
      if (my_task.eq.master_task) then
        write (50) istep1,time,time_forc
        write(6,*) 'Writing ',filename(1:lenstr(filename))
        write(6,*) 'Restart written ',istep1,time,time_forc
      endif

c   new ice
      work(ilo:ihi,jlo:jhi)=hice(ilo:ihi,jlo:jhi,0)
      call ice_write(50,0,work,'ruf8')
      work(ilo:ihi,jlo:jhi)=compact(ilo:ihi,jlo:jhi,0)
      call ice_write(50,0,work,'ruf8')

      do n=1,ncat

c   thicknesses
      work(ilo:ihi,jlo:jhi)=hice(ilo:ihi,jlo:jhi,n)
      call ice_write(50,0,work,'ruf8')
      work(ilo:ihi,jlo:jhi)=hsnow(ilo:ihi,jlo:jhi,n)
      call ice_write(50,0,work,'ruf8')
      work(ilo:ihi,jlo:jhi)=compact(ilo:ihi,jlo:jhi,n)
      call ice_write(50,0,work,'ruf8')

c   thermodynamics
      work(ilo:ihi,jlo:jhi)=Tsfc(ilo:ihi,jlo:jhi,n)
      call ice_write(50,0,work,'ruf8')
      do k=0,nlyr
        work(ilo:ihi,jlo:jhi)=qice(ilo:ihi,jlo:jhi,n,k)
        call ice_write(50,0,work,'ruf8')
      enddo
      work(ilo:ihi,jlo:jhi)=shcoef(ilo:ihi,jlo:jhi,n)
      call ice_write(50,0,work,'ruf8')
      work(ilo:ihi,jlo:jhi)=lhcoef(ilo:ihi,jlo:jhi,n)
      call ice_write(50,0,work,'ruf8')

      enddo    ! categories

c   wind stress
      work(ilo:ihi,jlo:jhi)=strairx(ilo:ihi,jlo:jhi)
      call ice_write(50,0,work,'ruf8')
      work(ilo:ihi,jlo:jhi)=strairy(ilo:ihi,jlo:jhi)
      call ice_write(50,0,work,'ruf8')

c   velocity
      work(ilo:ihi,jlo:jhi)=u(ilo:ihi,jlo:jhi)
      call ice_write(50,0,work,'ruf8')
      work(ilo:ihi,jlo:jhi)=v(ilo:ihi,jlo:jhi)
      call ice_write(50,0,work,'ruf8')

c   stress
      work(ilo:ihi,jlo:jhi)=sig11ne(ilo:ihi,jlo:jhi)
      call ice_write(50,0,work,'ruf8')
      work(ilo:ihi,jlo:jhi)=sig11se(ilo:ihi,jlo:jhi)
      call ice_write(50,0,work,'ruf8')
      work(ilo:ihi,jlo:jhi)=sig11sw(ilo:ihi,jlo:jhi)
      call ice_write(50,0,work,'ruf8')
      work(ilo:ihi,jlo:jhi)=sig11nw(ilo:ihi,jlo:jhi)
      call ice_write(50,0,work,'ruf8')
      work(ilo:ihi,jlo:jhi)=sig12ne(ilo:ihi,jlo:jhi)
      call ice_write(50,0,work,'ruf8')
      work(ilo:ihi,jlo:jhi)=sig12se(ilo:ihi,jlo:jhi)
      call ice_write(50,0,work,'ruf8')
      work(ilo:ihi,jlo:jhi)=sig12sw(ilo:ihi,jlo:jhi)
      call ice_write(50,0,work,'ruf8')
      work(ilo:ihi,jlo:jhi)=sig12nw(ilo:ihi,jlo:jhi)
      call ice_write(50,0,work,'ruf8')
      work(ilo:ihi,jlo:jhi)=sig22ne(ilo:ihi,jlo:jhi)
      call ice_write(50,0,work,'ruf8')
      work(ilo:ihi,jlo:jhi)=sig22se(ilo:ihi,jlo:jhi)
      call ice_write(50,0,work,'ruf8')
      work(ilo:ihi,jlo:jhi)=sig22sw(ilo:ihi,jlo:jhi)
      call ice_write(50,0,work,'ruf8')
      work(ilo:ihi,jlo:jhi)=sig22nw(ilo:ihi,jlo:jhi)
      call ice_write(50,0,work,'ruf8')

      if (my_task.eq.master_task) close(50)
      return
      end
c     ------------------------------------------------------
      subroutine restartfile
c.. Restarts from a dump
      implicit none
      include "ice.h"          ! ice code common blocks
      integer i, j, k, n, lenstr

      real work(ilo:ihi,jlo:jhi)

      call ice_open(50,restrt_file,0)
      if (my_task.eq.master_task) then
        read (50) istep0,time,time_forc
        write(6,*) 'read ',restrt_file(1:lenstr(restrt_file))
        write(6,*) 'dump read at istep=',istep0,time,time_forc
      endif
      call ice_bcast_iscalar(istep0)
      istep1=istep0
      call ice_bcast_rscalar(time)
      call ice_bcast_rscalar(time_forc)

c   new ice
      call ice_read(50,0,work,'ruf8')
      hice(ilo:ihi,jlo:jhi,0)=work(ilo:ihi,jlo:jhi)
      call ice_read(50,0,work,'ruf8')
      compact(ilo:ihi,jlo:jhi,0)=work(ilo:ihi,jlo:jhi)

      do n=1,ncat

c   thicknesses
      call ice_read(50,0,work,'ruf8')
      hice(ilo:ihi,jlo:jhi,n)=work(ilo:ihi,jlo:jhi)
      call ice_read(50,0,work,'ruf8')
      hsnow(ilo:ihi,jlo:jhi,n)=work(ilo:ihi,jlo:jhi)
      call ice_read(50,0,work,'ruf8')
      compact(ilo:ihi,jlo:jhi,n)=work(ilo:ihi,jlo:jhi)

c   thermodynamics
      call ice_read(50,0,work,'ruf8')
      Tsfc(ilo:ihi,jlo:jhi,n)=work(ilo:ihi,jlo:jhi)
      do k=0,nlyr
        call ice_read(50,0,work,'ruf8')
        qice(ilo:ihi,jlo:jhi,n,k)=work(ilo:ihi,jlo:jhi)
      enddo
      call ice_read(50,0,work,'ruf8')
      shcoef(ilo:ihi,jlo:jhi,n)=work(ilo:ihi,jlo:jhi)
      call ice_read(50,0,work,'ruf8')
      lhcoef(ilo:ihi,jlo:jhi,n)=work(ilo:ihi,jlo:jhi)

      enddo    ! categories

c   wind stress
      call ice_read(50,0,work,'ruf8')
      strairx(ilo:ihi,jlo:jhi)=work(ilo:ihi,jlo:jhi)
      call ice_read(50,0,work,'ruf8')
      strairy(ilo:ihi,jlo:jhi)=work(ilo:ihi,jlo:jhi)

c   velocity
      call ice_read(50,0,work,'ruf8')
      u(ilo:ihi,jlo:jhi)=work(ilo:ihi,jlo:jhi)
      call ice_read(50,0,work,'ruf8')
      v(ilo:ihi,jlo:jhi)=work(ilo:ihi,jlo:jhi)

c   stress
      call ice_read(50,0,work,'ruf8')
      sig11ne(ilo:ihi,jlo:jhi)=work(ilo:ihi,jlo:jhi)
      call ice_read(50,0,work,'ruf8')
      sig11se(ilo:ihi,jlo:jhi)=work(ilo:ihi,jlo:jhi)
      call ice_read(50,0,work,'ruf8')
      sig11sw(ilo:ihi,jlo:jhi)=work(ilo:ihi,jlo:jhi)
      call ice_read(50,0,work,'ruf8')
      sig11nw(ilo:ihi,jlo:jhi)=work(ilo:ihi,jlo:jhi)
      call ice_read(50,0,work,'ruf8')
      sig12ne(ilo:ihi,jlo:jhi)=work(ilo:ihi,jlo:jhi)
      call ice_read(50,0,work,'ruf8')
      sig12se(ilo:ihi,jlo:jhi)=work(ilo:ihi,jlo:jhi)
      call ice_read(50,0,work,'ruf8')
      sig12sw(ilo:ihi,jlo:jhi)=work(ilo:ihi,jlo:jhi)
      call ice_read(50,0,work,'ruf8')
      sig12nw(ilo:ihi,jlo:jhi)=work(ilo:ihi,jlo:jhi)
      call ice_read(50,0,work,'ruf8')
      sig22ne(ilo:ihi,jlo:jhi)=work(ilo:ihi,jlo:jhi)
      call ice_read(50,0,work,'ruf8')
      sig22se(ilo:ihi,jlo:jhi)=work(ilo:ihi,jlo:jhi)
      call ice_read(50,0,work,'ruf8')
      sig22sw(ilo:ihi,jlo:jhi)=work(ilo:ihi,jlo:jhi)
      call ice_read(50,0,work,'ruf8')
      sig22nw(ilo:ihi,jlo:jhi)=work(ilo:ihi,jlo:jhi)

      if (my_task.eq.master_task) close(50)

c.. update boundary conditions
      do n = 0,ncat
        call bound(hice(:,:,n))
        call bound(hsnow(:,:,n))
        call bound(compact(:,:,n))
        call bound(qice(:,:,n,1))
        call bound(qice(:,:,n,2))
        call bound(Tsfc(:,:,n))
      end do
      call bound(u)      
      call bound(v)      

      return
      end
c-----------------------------------------------------------------------
      subroutine runtime_diags
c.. Writes diagnostic info (max, min, global sums, etc) to standard out
      implicit none
      include "ice.h"          ! ice code common blocks
      real umaxn, hmaxn, shmaxn, arean
      real umaxs, hmaxs, shmaxs, areas
      real sfreshn,ssaltn,sheatn
      real sfreshs,ssalts,sheats
      real work(imt_local,jmt_local)
      real work1(ilo:ihi,jlo:jhi),work2(ilo:ihi,jlo:jhi)
      integer n

c.. maximum velocity
      work  = max(abs(mask_n*u),abs(mask_n*v))
      umaxn = ice_global_real_maxval(imt_local*jmt_local,WORK)
      work  = max(abs(mask_s*u),abs(mask_s*v))
      umaxs = ice_global_real_maxval(imt_local*jmt_local,WORK)
c.. maximum effective thickness
      work  = mask_n(:,:)*compact(:,:,ncat)*hice(:,:,ncat)
      hmaxn = ice_global_real_maxval(imt_local*jmt_local,WORK)
      work  = mask_s(:,:)*compact(:,:,ncat)*hice(:,:,ncat)
      hmaxs = ice_global_real_maxval(imt_local*jmt_local,WORK)
c.. total ice area
      work2 = tarea(ilo:ihi,jlo:jhi) * ifrc(ilo:ihi,jlo:jhi)
      work1 = mask_n(ilo:ihi,jlo:jhi)*work2
      arean = ice_global_real_sum((jhi-jlo+1)*(ihi-ilo+1),WORK1)
     & / 1.e6  ! km^2
      work2 = mask_s(ilo:ihi,jlo:jhi)*work2
      areas = ice_global_real_sum((jhi-jlo+1)*(ihi-ilo+1),WORK2)
     & / 1.e6  ! km^2
c.. total ice volume
      work2 = 0.
      do n = 1,ncat
      work2 = work2 + tarea(ilo:ihi,jlo:jhi)
     &      * compact(ilo:ihi,jlo:jhi,n)*hice(ilo:ihi,jlo:jhi,n)
      enddo
      work1 = mask_n(ilo:ihi,jlo:jhi)*work2
      shmaxn= ice_global_real_sum((jhi-jlo+1)*(ihi-ilo+1),WORK1)
      work2 = mask_s(ilo:ihi,jlo:jhi)*work2
      shmaxs= ice_global_real_sum((jhi-jlo+1)*(ihi-ilo+1),WORK2)
c.. various fluxes
      work2 = mask_n(ilo:ihi,jlo:jhi)*fresh(ilo:ihi,jlo:jhi,ncat)
      sfreshn = ice_global_real_sum((jhi-jlo+1)*(ihi-ilo+1),WORK2)
      work2 = mask_n(ilo:ihi,jlo:jhi)*fsalt(ilo:ihi,jlo:jhi,ncat)
      ssaltn = ice_global_real_sum((jhi-jlo+1)*(ihi-ilo+1),WORK2)
      work2 = mask_n(ilo:ihi,jlo:jhi)*fhnet(ilo:ihi,jlo:jhi,ncat)
      sheatn = ice_global_real_sum((jhi-jlo+1)*(ihi-ilo+1),WORK2)
      work2 = mask_s(ilo:ihi,jlo:jhi)*fresh(ilo:ihi,jlo:jhi,ncat)
      sfreshs = ice_global_real_sum((jhi-jlo+1)*(ihi-ilo+1),WORK2)
      work2 = mask_s(ilo:ihi,jlo:jhi)*fsalt(ilo:ihi,jlo:jhi,ncat)
      ssalts = ice_global_real_sum((jhi-jlo+1)*(ihi-ilo+1),WORK2)
      work2 = mask_s(ilo:ihi,jlo:jhi)*fhnet(ilo:ihi,jlo:jhi,ncat)
      sheats = ice_global_real_sum((jhi-jlo+1)*(ihi-ilo+1),WORK2)

      if (my_task.eq.master_task) then
        write (6,899) 'Arctic','Antarctic'
        write (6,900) 'max cH (m)           = ',hmaxn,hmaxs
        write (6,901) 'total area (km^2)    = ',arean,areas
        write (6,901) 'total volume (m^3)   = ',shmaxn,shmaxs
        write (6,900) 'max u, v (m/s)       = ',umaxn,umaxs
        write (6,901) 'fresh h2o (kg/m^2 s) = ',sfreshn,sfreshs
        write (6,901) 'salt flux (kg/m^2 s) = ',ssaltn,ssalts
        write (6,901) 'heat used (W/m^2)    = ',sheatn,sheats
      endif

 899  format (27x,a24,2x,a24)
 900  format (a25,2x,f24.17,2x,f24.17)
 901  format (a25,2x,1pe24.17,2x,1pe24.17)
      return
      end
c----------------------------------------------------------------------
      subroutine grid_netcdf
c.. writes grid and land masks in netCDF format
      implicit none
      include "ice.h"          ! ice code common blocks
      include "netcdf.inc"     ! include netcdf
      integer ncID                             ! file ID
      character (80) :: loc_fn                 ! netCDF file name
      integer jmtID, imtID                     ! dimension IDs
      integer STATUS                           ! error flag
      integer i, j
      integer mRANK                            ! rank of field
      parameter (mRANK=2)                      ! 2D
      integer mSHP(mRANK)                      ! shape of field
      real work(ilo:ihi,jlo:jhi)
      real workg(imt_global,jmt_global)

c.......... grid latitude/longitude ............
      real*4  nc_lat(imt_global,jmt_global)    ! latitude
      real*4  nc_lon(imt_global,jmt_global)    ! longitude
      integer latID                            ! latitude ID
      integer lonID                            ! longitude ID

c.......... masks ............
      real*4  nc_uvm(imt_global,jmt_global)    ! U-cell mask
      real*4  nc_hm (imt_global,jmt_global)    ! T-cell mask
      integer uvmID                            ! uvm ID
      integer hmID                             ! hm ID

      work(ilo:ihi,jlo:jhi) = ULAT(ilo:ihi,jlo:jhi)
      call global_gather(workg,work)
      nc_lat=workg
      work(ilo:ihi,jlo:jhi) = ULONG(ilo:ihi,jlo:jhi)
      call global_gather(workg,work)
      nc_lon=workg
      work(ilo:ihi,jlo:jhi) = uvm(ilo:ihi,jlo:jhi)
      call global_gather(workg,work)
      nc_uvm=workg
      work(ilo:ihi,jlo:jhi) = hm(ilo:ihi,jlo:jhi)
      call global_gather(workg,work)
      nc_hm=workg

      if (my_task.eq.master_task) then
c.. file name
      loc_fn = 'icegrid.nc'
      write (6,*) ' Creating file ',loc_fn
      STATUS = NF_CREATE(loc_fn,NF_CLOBBER,ncID) 
      STATUS = NF_DEF_DIM(ncID, 'x', imt_global, imtID)
      STATUS = NF_DEF_DIM(ncID, 'y', jmt_global, jmtID)
      mSHP(1) = imtID
      mSHP(2) = jmtID
      STATUS = NF_DEF_VAR(ncID, 'lat', NF_REAL, mRANK, mSHP, latID)
      STATUS = NF_DEF_VAR(ncID, 'lon', NF_REAL, mRANK, mSHP, lonID)
      STATUS = NF_DEF_VAR(ncID, 'uvm', NF_REAL, mRANK, mSHP, uvmID)
      STATUS = NF_DEF_VAR(ncID,  'hm', NF_REAL, mRANK, mSHP,  hmID)
      STATUS = NF_PUT_ATT_TEXT(ncID, NF_GLOBAL, 'title', 22,     
     &   'Lat/Lon and Land Masks')
      STATUS = NF_PUT_ATT_TEXT(ncID, latID, 'long_name', 8,        
     &   'Latitude')
      STATUS = NF_PUT_ATT_TEXT(ncID, lonID, 'long_name', 9,        
     &   'Longitude')
      STATUS = NF_PUT_ATT_TEXT(ncID, uvmID, 'long_name', 11,        
     &   'U-cell mask')
      STATUS = NF_PUT_ATT_TEXT(ncID, hmID, 'long_name', 11,        
     &   'T-cell mask')
      STATUS = NF_ENDDEF(ncID)           
      STATUS = NF_PUT_VAR_REAL (ncID, latID, nc_lat) 
          if (STATUS .ne. NF_NOERR) write (6,*) 'error writing lat'
      STATUS = NF_PUT_VAR_REAL (ncID, lonID, nc_lon) 
          if (STATUS .ne. NF_NOERR) write (6,*) 'error writing lon'
      STATUS = NF_PUT_VAR_REAL (ncID, uvmID, nc_uvm) 
          if (STATUS .ne. NF_NOERR) write (6,*) 'error writing uvm'
      STATUS = NF_PUT_VAR_REAL (ncID,  hmID, nc_hm) 
          if (STATUS .ne. NF_NOERR) write (6,*) 'error writing hm'
      STATUS = NF_CLOSE(ncID)
      if (STATUS .ne. NF_NOERR) write (6,*) 'error closing mask_file'

      endif

      return
      end
c-------------------------------------------------------------------------
