!=======================================================================
!BOP
!
! !MODULE: ice_dyn_eap - elastic-anisotropic-plastic sea ice 
!                        dynamics model
!
! !DESCRIPTION:
!
! Elastic-anisotropic sea ice dynamics model
! Computes ice velocity and deformation
!
! See:
!
! Wilchinsky, A.V. and D.L. Feltham (2006). Modelling the rheology of 
! sea ice as a collection of diamond-shaped floes. 
! Journal of Non-Newtonian Fluid Mechanics, 138(1), 22-32.
!
! Tsamados, M., D.L. Feltham, and A.V. Wilchinsky (2012). Impact on new
! anisotropic rheology on simulations of Arctic sea ice. JGR, in press.
!
! !REVISION HISTORY:
!  SVN:$Id: ice_dyn_eap.F90 $
!
! authors: Michel Tsamados, CPOM 
!          David Schroeder, CPOM
!
! !INTERFACE:
!
      module ice_dyn_eap
!
! !USES:
!
      use ice_kinds_mod
      use ice_blocks, only: nx_block, ny_block
      use ice_domain_size, only: max_blocks
!
!EOP
!
      implicit none
      private
      public :: eap, init_eap, write_restart_eap
      save

      ! Look-up table needed for calculating structure tensor
      integer (int_kind), parameter :: & 
        nx_yield            =  21, &
        ny_yield            =  21, &
        na_yield            =  11

      real (kind=dbl_kind), dimension (nx_yield,ny_yield,na_yield) :: & 
        s11r, s12r, s22r, s11s, s12s, s22s           

      real (kind=dbl_kind), dimension (nx_block,ny_block,max_blocks) :: &
         a11_1, a11_2, a11_3, a11_4,                  & ! components of 
         a12_1, a12_2, a12_3, a12_4                     ! structure tensor

      ! history
      real (kind=dbl_kind), dimension(nx_block,ny_block,max_blocks), public :: &
         e11      , & ! components of strain rate tensor (1/s)
         e12      , & 
         e22      , &
         yieldstress11, & ! components of yield stress tensor (kg/s^2)
         yieldstress12, &
         yieldstress22, &
         s11      , & ! components of stress tensor (kg/s^2)
         s12      , &
         s22      , &
         a11      , & ! components of structure tensor ()
         a12

!=======================================================================

      contains

!=======================================================================
!BOP
!
! !IROUTINE: eap - elastic-anisotropic-plastic dynamics driver
!
! !INTERFACE:
!
      subroutine eap (dt)
!
! !DESCRIPTION:
!
! Elastic-anisotropic-plastic dynamics driver
! based on subroutine evp
!
#ifdef CICE_IN_NEMO
! via NEMO (unless calc_strair is true).  These values are supplied  
! rotated on u grid and multiplied by aice.  strairxT = 0 in this  
! case so operations in evp_prep1 are pointless but carried out to  
! minimise code changes.
#endif
!
! !REVISION HISTORY:
!
! same as module
!
! !USES:
!
      use ice_atmo, only: Cdn_ocn
      use ice_boundary, only: ice_halo, ice_HaloMask, ice_HaloUpdate, &
          ice_HaloDestroy
      use ice_blocks, only: block, get_block
      use ice_constants, only: field_loc_center, field_loc_NEcorner, &
          field_type_scalar, field_type_vector, c0, p5
      use ice_domain, only: nblocks, blocks_ice, halo_info, maskhalo_dyn
      use ice_dyn_shared, only: fcor_blk, ndte, dtei, a_min, m_min, &
          cosw, sinw, denom1, uvel_init, vvel_init, arlx1i, &
          evp_prep1, evp_prep2, stepu, evp_finish
      use ice_flux, only: rdg_conv, rdg_shear, prs_sig, strairxT, strairyT, &
          strairx, strairy, uocn, vocn, ss_tltx, ss_tlty, iceumask, fm, &
          strtltx, strtlty, strocnx, strocny, strintx, strinty, &
          strocnxT, strocnyT, &
          stressp_1, stressp_2, stressp_3, stressp_4, &
          stressm_1, stressm_2, stressm_3, stressm_4, &
          stress12_1, stress12_2, stress12_3, stress12_4
      use ice_grid, only: tmask, umask, dxt, dyt, dxhy, dyhx, cxp, cyp, cxm, cym, &
          tarear, uarear, tinyarea, to_ugrid, t2ugrid_vector, u2tgrid_vector
      use ice_mechred, only: ice_strength
      use ice_state, only: aice, vice, vsno, uvel, vvel, divu, shear, &
          aice_init, aice0, aicen, vicen, strength
      use ice_timers, only: timer_dynamics, timer_bound, &
          ice_timer_start, ice_timer_stop
!
! !INPUT/OUTPUT PARAMETERS:
!
      real (kind=dbl_kind), intent(in) :: &
         dt      ! time step
!
! local variables
!
      integer (kind=int_kind) :: & 
         ksub           , & ! subcycle step
         iblk           , & ! block index
         ilo,ihi,jlo,jhi, & ! beginning and end of physical domain
         i, j

      integer (kind=int_kind), dimension(max_blocks) :: & 
         icellt   , & ! no. of cells where icetmask = 1
         icellu       ! no. of cells where iceumask = 1

      integer (kind=int_kind), dimension (nx_block*ny_block, max_blocks) :: &
         indxti   , & ! compressed index in i-direction
         indxtj   , & ! compressed index in j-direction
         indxui   , & ! compressed index in i-direction
         indxuj       ! compressed index in j-direction

      real (kind=dbl_kind), dimension (nx_block,ny_block,max_blocks) :: &
         tmass    , & ! total mass of ice and snow (kg/m^2)
         waterx   , & ! for ocean stress calculation, x (m/s)
         watery   , & ! for ocean stress calculation, y (m/s)
         forcex   , & ! work array: combined atm stress and ocn tilt, x
         forcey   , & ! work array: combined atm stress and ocn tilt, y
         aiu      , & ! ice fraction on u-grid
         umass    , & ! total mass of ice and snow (u grid)
         umassdti     ! mass of U-cell/dte (kg/m^2 s)

      real (kind=dbl_kind), allocatable :: fld2(:,:,:,:)

      real (kind=dbl_kind), dimension(nx_block,ny_block,8):: &
         str          ! stress combinations for momentum equation

      integer (kind=int_kind), dimension (nx_block,ny_block,max_blocks) :: &
         icetmask   ! ice extent mask (T-cell)

      integer (kind=int_kind), dimension (nx_block,ny_block,max_blocks) :: &
         halomask     ! mask for masked halo creation

      type (ice_halo) :: &
         halo_info_mask !  ghost cell update info for masked halo

      type (block) :: &
         this_block           ! block information for current block
      
      call ice_timer_start(timer_dynamics) ! dynamics

      !-----------------------------------------------------------------
      ! Initialize
      !-----------------------------------------------------------------

      allocate(fld2(nx_block,ny_block,2,max_blocks))

       ! This call is needed only if dt changes during runtime.
!      call set_evp_parameters (dt)

      !$OMP PARALLEL DO PRIVATE(iblk,i,j,ilo,ihi,jlo,jhi,this_block)
      do iblk = 1, nblocks
         do j = 1, ny_block 
         do i = 1, nx_block 
            rdg_conv (i,j,iblk) = c0 
            rdg_shear(i,j,iblk) = c0 
            divu (i,j,iblk) = c0 
            shear(i,j,iblk) = c0 
            prs_sig(i,j,iblk) = c0 
            e11(i,j,iblk) = c0
            e12(i,j,iblk) = c0
            e22(i,j,iblk) = c0
            s11(i,j,iblk) = c0
            s12(i,j,iblk) = c0
            s22(i,j,iblk) = c0
            yieldstress11(i,j,iblk) = c0
            yieldstress12(i,j,iblk) = c0
            yieldstress22(i,j,iblk) = c0
         enddo
         enddo

      !-----------------------------------------------------------------
      ! preparation for dynamics
      !-----------------------------------------------------------------

         this_block = get_block(blocks_ice(iblk),iblk)         
         ilo = this_block%ilo
         ihi = this_block%ihi
         jlo = this_block%jlo
         jhi = this_block%jhi

         call evp_prep1 (nx_block,           ny_block,           & 
                         ilo, ihi,           jlo, jhi,           &
                         aice    (:,:,iblk), vice    (:,:,iblk), & 
                         vsno    (:,:,iblk), tmask   (:,:,iblk), & 
                         strairxT(:,:,iblk), strairyT(:,:,iblk), & 
                         strairx (:,:,iblk), strairy (:,:,iblk), & 
                         tmass   (:,:,iblk), icetmask(:,:,iblk))

      enddo                     ! iblk
      !$OMP END PARALLEL DO

      call ice_timer_start(timer_bound)
      call ice_HaloUpdate (icetmask,          halo_info, &
                           field_loc_center,  field_type_scalar)
      call ice_timer_stop(timer_bound)

      !-----------------------------------------------------------------
      ! convert fields from T to U grid
      !-----------------------------------------------------------------

      call to_ugrid(tmass,umass)
      call to_ugrid(aice_init, aiu)

#ifdef CICE_IN_NEMO
      !----------------------------------------------------------------
      ! Set wind stress to values supplied via NEMO
      ! This wind stress is rotated on u grid and multiplied by aice
      !----------------------------------------------------------------
      if (.not. calc_strair) then       
         strairx(:,:,:) = strax(:,:,:)
         strairy(:,:,:) = stray(:,:,:)
      else
#endif
         call t2ugrid_vector(strairx)
         call t2ugrid_vector(strairy)
#ifdef CICE_IN_NEMO
      endif
#endif

      !$OMP PARALLEL DO PRIVATE(iblk,ilo,ihi,jlo,jhi,this_block)
      do iblk = 1, nblocks

      !-----------------------------------------------------------------
      ! more preparation for dynamics
      !-----------------------------------------------------------------

         this_block = get_block(blocks_ice(iblk),iblk)         
         ilo = this_block%ilo
         ihi = this_block%ihi
         jlo = this_block%jlo
         jhi = this_block%jhi

         call evp_prep2 (nx_block,             ny_block,             & 
                         ilo, ihi,             jlo, jhi,             &
                         icellt(iblk),         icellu(iblk),         & 
                         indxti      (:,iblk), indxtj      (:,iblk), & 
                         indxui      (:,iblk), indxuj      (:,iblk), & 
                         aiu       (:,:,iblk), umass     (:,:,iblk), & 
                         umassdti  (:,:,iblk), fcor_blk  (:,:,iblk), & 
                         umask     (:,:,iblk),                       & 
                         uocn      (:,:,iblk), vocn      (:,:,iblk), & 
                         strairx   (:,:,iblk), strairy   (:,:,iblk), & 
                         ss_tltx   (:,:,iblk), ss_tlty   (:,:,iblk), &  
                         icetmask  (:,:,iblk), iceumask  (:,:,iblk), & 
                         fm        (:,:,iblk), dt,                   & 
                         strtltx   (:,:,iblk), strtlty   (:,:,iblk), & 
                         strocnx   (:,:,iblk), strocny   (:,:,iblk), & 
                         strintx   (:,:,iblk), strinty   (:,:,iblk), & 
                         waterx    (:,:,iblk), watery    (:,:,iblk), & 
                         forcex    (:,:,iblk), forcey    (:,:,iblk), & 
                         stressp_1 (:,:,iblk), stressp_2 (:,:,iblk), & 
                         stressp_3 (:,:,iblk), stressp_4 (:,:,iblk), & 
                         stressm_1 (:,:,iblk), stressm_2 (:,:,iblk), & 
                         stressm_3 (:,:,iblk), stressm_4 (:,:,iblk), & 
                         stress12_1(:,:,iblk), stress12_2(:,:,iblk), & 
                         stress12_3(:,:,iblk), stress12_4(:,:,iblk), & 
                         uvel_init (:,:,iblk), vvel_init (:,:,iblk), &
                         uvel      (:,:,iblk), vvel      (:,:,iblk))

      !-----------------------------------------------------------------
      ! Initialize structure tensor
      !-----------------------------------------------------------------

         do j = 1, ny_block
         do i = 1, nx_block
            if (icetmask(i,j,iblk)==0) then
            ! structure tensor
               a11_1(i,j,iblk) = p5
               a11_2(i,j,iblk) = p5
               a11_3(i,j,iblk) = p5
               a11_4(i,j,iblk) = p5
               a12_1(i,j,iblk) = c0
               a12_2(i,j,iblk) = c0
               a12_3(i,j,iblk) = c0
               a12_4(i,j,iblk) = c0
            endif                  ! icetmask
         enddo                     ! i
         enddo                     ! j

      !-----------------------------------------------------------------
      ! ice strength
      ! New strength used in Ukita Moritz rheology         
      !-----------------------------------------------------------------

         call ice_strength (nx_block, ny_block,   & 
                            ilo, ihi, jlo, jhi,   &
                            icellt(iblk),         & 
                            indxti      (:,iblk), & 
                            indxtj      (:,iblk), & 
                            aice    (:,:,  iblk), & 
                            vice    (:,:,  iblk), & 
                            aice0   (:,:,  iblk), & 
                            aicen   (:,:,:,iblk), &  
                            vicen   (:,:,:,iblk), & 
                            strength(:,:,  iblk))

         ! load velocity into array for boundary updates
         fld2(:,:,1,iblk) = uvel(:,:,iblk)
         fld2(:,:,2,iblk) = vvel(:,:,iblk)

      enddo  ! iblk
      !$OMP END PARALLEL DO

      call ice_timer_start(timer_bound)
      call ice_HaloUpdate (strength,           halo_info, &
                           field_loc_center,   field_type_scalar)
      ! velocities may have changed in evp_prep2
      call ice_HaloUpdate (fld2,               halo_info, &
                           field_loc_NEcorner, field_type_vector)

      ! unload
      !$OMP PARALLEL DO PRIVATE(iblk)
      do iblk = 1,nblocks
         uvel(:,:,iblk) = fld2(:,:,1,iblk)
         vvel(:,:,iblk) = fld2(:,:,2,iblk)
      enddo
      !$OMP END PARALLEL DO

      if (maskhalo_dyn) &
         call ice_HaloMask(halo_info_mask, halo_info, icetmask)
      call ice_timer_stop(timer_bound)

      do ksub = 1,ndte        ! subcycling

      !-----------------------------------------------------------------
      ! stress tensor equation, total surface stress
      !-----------------------------------------------------------------

         !$OMP PARALLEL DO PRIVATE(iblk)
         do iblk = 1, nblocks

            call stress_eap  (nx_block,             ny_block,             &
                              ksub,                 ndte,                 &
                              icellt(iblk),                               &
                              indxti      (:,iblk), indxtj      (:,iblk), &
                              arlx1i,               denom1,         &
                              uvel      (:,:,iblk), vvel      (:,:,iblk), &
                              dxt       (:,:,iblk), dyt       (:,:,iblk), &
                              dxhy      (:,:,iblk), dyhx      (:,:,iblk), &
                              cxp       (:,:,iblk), cyp       (:,:,iblk), &
                              cxm       (:,:,iblk), cym       (:,:,iblk), &
                              tarear    (:,:,iblk), tinyarea  (:,:,iblk), &
                              strength  (:,:,iblk),                       &
                              a11       (:,:,iblk), a12  (:,:,iblk),      &
                              a11_1     (:,:,iblk), a11_2   (:,:,iblk),   &
                              a11_3     (:,:,iblk), a11_4   (:,:,iblk),   &
                              a12_1     (:,:,iblk), a12_2   (:,:,iblk),   &
                              a12_3     (:,:,iblk), a12_4   (:,:,iblk),   &
                              stressp_1 (:,:,iblk), stressp_2 (:,:,iblk), &
                              stressp_3 (:,:,iblk), stressp_4 (:,:,iblk), &
                              stressm_1 (:,:,iblk), stressm_2 (:,:,iblk), &
                              stressm_3 (:,:,iblk), stressm_4 (:,:,iblk), &
                              stress12_1(:,:,iblk), stress12_2(:,:,iblk), &
                              stress12_3(:,:,iblk), stress12_4(:,:,iblk), &
                              shear     (:,:,iblk), divu      (:,:,iblk), &
                              e11       (:,:,iblk), e12       (:,:,iblk), &
                              e22       (:,:,iblk),                       &
                              s11       (:,:,iblk), s12       (:,:,iblk), &
                              s22       (:,:,iblk),                       &
                              yieldstress11 (:,:,iblk),                   &
                              yieldstress12 (:,:,iblk),                           &
                              yieldstress22 (:,:,iblk),                   &
                              prs_sig   (:,:,iblk),                       &
                              rdg_conv  (:,:,iblk), rdg_shear (:,:,iblk), &
                              str       (:,:,:))

      !-----------------------------------------------------------------
      ! momentum equation
      !-----------------------------------------------------------------

            call stepu (nx_block,            ny_block,           & 
                        icellu       (iblk), Cdn_ocn (:,:,iblk), & 
                        indxui     (:,iblk), indxuj    (:,iblk), & 
                        aiu      (:,:,iblk), str     (:,:,:),    & 
                        uocn     (:,:,iblk), vocn    (:,:,iblk), &     
                        waterx   (:,:,iblk), watery  (:,:,iblk), & 
                        forcex   (:,:,iblk), forcey  (:,:,iblk), & 
                        umassdti (:,:,iblk), fm      (:,:,iblk), & 
                        uarear   (:,:,iblk),                     & 
                        strocnx  (:,:,iblk), strocny (:,:,iblk), & 
                        strintx  (:,:,iblk), strinty (:,:,iblk), & 
                        uvel_init(:,:,iblk), vvel_init(:,:,iblk),&
                        uvel     (:,:,iblk), vvel    (:,:,iblk))

            ! load velocity into array for boundary updates
            fld2(:,:,1,iblk) = uvel(:,:,iblk)
            fld2(:,:,2,iblk) = vvel(:,:,iblk)

      !-----------------------------------------------------------------
      ! evolution of structure tensor A
      !-----------------------------------------------------------------

            call stepa (nx_block,          ny_block,                &
                        dtei,              icellt     (iblk),       &
                        indxti   (:,iblk), indxtj    (:,iblk),      &
                        a11    (:,:,iblk), a12  (:,:,iblk),         &
                        a11_1  (:,:,iblk), a11_2   (:,:,iblk),      &
                        a11_3  (:,:,iblk), a11_4   (:,:,iblk),      &
                        a12_1  (:,:,iblk), a12_2   (:,:,iblk),      &
                        a12_3  (:,:,iblk), a12_4   (:,:,iblk),      &
                        stressp_1(:,:,iblk), stressp_2(:,:,iblk),   &
                        stressp_3(:,:,iblk), stressp_4(:,:,iblk),   &
                        stressm_1(:,:,iblk), stressm_2(:,:,iblk),   &
                        stressm_3(:,:,iblk), stressm_4(:,:,iblk),   &
                        stress12_1(:,:,iblk), stress12_2(:,:,iblk), &
                        stress12_3(:,:,iblk), stress12_4(:,:,iblk))

         enddo
         !$OMP END PARALLEL DO

         call ice_timer_start(timer_bound)
         if (maskhalo_dyn) then
            call ice_HaloUpdate (fld2,               halo_info_mask, &
                                 field_loc_NEcorner, field_type_vector)
         else
            call ice_HaloUpdate (fld2,               halo_info, &
                                 field_loc_NEcorner, field_type_vector)
         endif

         ! unload
         !$OMP PARALLEL DO PRIVATE(iblk)
         do iblk = 1,nblocks
            uvel(:,:,iblk) = fld2(:,:,1,iblk)
            vvel(:,:,iblk) = fld2(:,:,2,iblk)
         enddo
         !$OMP END PARALLEL DO
         call ice_timer_stop(timer_bound)

      enddo                     ! subcycling

      deallocate(fld2)
      if (maskhalo_dyn) call ice_HaloDestroy(halo_info_mask)

      !-----------------------------------------------------------------
      ! ice-ocean stress
      !-----------------------------------------------------------------

      !$OMP PARALLEL DO PRIVATE(iblk)
      do iblk = 1, nblocks

         call evp_finish                               & 
              (nx_block,           ny_block,           & 
               icellu      (iblk), Cdn_ocn (:,:,iblk), & 
               indxui    (:,iblk), indxuj    (:,iblk), & 
               uvel    (:,:,iblk), vvel    (:,:,iblk), & 
               uocn    (:,:,iblk), vocn    (:,:,iblk), & 
               aiu     (:,:,iblk), fm      (:,:,iblk), &
               strintx (:,:,iblk), strinty (:,:,iblk), &
               strairx (:,:,iblk), strairy (:,:,iblk), & 
               strocnx (:,:,iblk), strocny (:,:,iblk), & 
               strocnxT(:,:,iblk), strocnyT(:,:,iblk))

      enddo
      !$OMP END PARALLEL DO

      call u2tgrid_vector(strocnxT)    ! shift
      call u2tgrid_vector(strocnyT)

      call ice_timer_stop(timer_dynamics)    ! dynamics

      end subroutine eap

!=======================================================================
!BOP
!
! !IROUTINE: init_eap - initialize parameters needed for eap dynamics
!
! !INTERFACE:
!
      subroutine init_eap (dt)
!
! !DESCRIPTION:
!
! Initialize parameters and variables needed for the eap dynamics
! (based on init_evp)
!
! !REVISION HISTORY:
!
! same as module
!
! !USES:
!
      use ice_blocks, only: nx_block, ny_block
      use ice_communicate, only: my_task, master_task
      use ice_constants, only: c0, p5
      use ice_domain, only: nblocks
      use ice_dyn_shared, only: init_evp
      use ice_exit, only: abort_ice
      use ice_restart, only: runtype
      use ice_fileunits, only: nu_diag, nu_eap, eap_filename, &
          get_fileunit, release_fileunit
!
! !INPUT/OUTPUT PARAMETERS:
!
      real (kind=dbl_kind), intent(in) :: &
         dt      ! time step
!
! local variables
!
      integer (kind=int_kind) :: &
         i, j, k, &
         iblk, &         ! block index
         eap_error       ! eap input file i/o error flag

      call init_evp (dt)

      !$OMP PARALLEL DO PRIVATE(iblk,i,j)
      do iblk = 1, nblocks
      do j = 1, ny_block
      do i = 1, nx_block
         e11(i,j,iblk) = c0
         e12(i,j,iblk) = c0
         e22(i,j,iblk) = c0
         s11(i,j,iblk) = c0
         s12(i,j,iblk) = c0
         s22(i,j,iblk) = c0
         yieldstress11(i,j,iblk) = c0
         yieldstress12(i,j,iblk) = c0
         yieldstress22(i,j,iblk) = c0
         a11_1 (i,j,iblk) = p5
         a11_2 (i,j,iblk) = p5
         a11_3 (i,j,iblk) = p5
         a11_4 (i,j,iblk) = p5
         a12_1 (i,j,iblk) = c0
         a12_2 (i,j,iblk) = c0
         a12_3 (i,j,iblk) = c0
         a12_4 (i,j,iblk) = c0
      enddo                     ! i
      enddo                     ! j
      enddo                     ! iblk
      !$OMP END PARALLEL DO

      !-----------------------------------------------------------------
      ! read stresses for eap dynamics (see Appendix A1)
      !-----------------------------------------------------------------

          call get_fileunit(nu_eap)
          open(nu_eap,file=eap_filename,status='old',iostat=eap_error)
          if (my_task == master_task) write(nu_diag,*) 'Reading eap stresses'
          do k = 1, na_yield
          do i = 1, nx_yield
          do j = 1, ny_yield
             if (eap_error == 0) read(nu_eap,*,iostat=eap_error) &
                      s11r(i,j,k), s12r(i,j,k), s22r(i,j,k), &
                      s11s(i,j,k), s12s(i,j,k), s22s(i,j,k)
          end do
          end do
          end do
          if (my_task == master_task) then
             if (eap_error == 0) close(nu_eap)
             if (eap_error /= 0) &
                          call abort_ice('ice: error reading eap stresses')
          endif
          call release_fileunit(nu_eap)

      if (runtype == 'continue') call read_restart_eap

      end subroutine init_eap

!=======================================================================
!BOP
!
! !IROUTINE: stress_eap - computes strain rates and 
!                         internal stress components
!
! !INTERFACE:
!
      subroutine stress_eap  (nx_block,   ny_block,       &
                              ksub,       ndte,           &
                              icellt,                     &
                              indxti,     indxtj,         &
                              arlx1i,     denom1,         &
                              uvel,       vvel,           &
                              dxt,        dyt,            &
                              dxhy,       dyhx,           &
                              cxp,        cyp,            &
                              cxm,        cym,            &
                              tarear,     tinyarea,       &
                              strength,                   &
                              a11, a12,                   &
                              a11_1, a11_2, a11_3, a11_4, &
                              a12_1, a12_2, a12_3, a12_4, &
                              stressp_1,  stressp_2,      &
                              stressp_3,  stressp_4,      &
                              stressm_1,  stressm_2,      &
                              stressm_3,  stressm_4,      &
                              stress12_1, stress12_2,     &
                              stress12_3, stress12_4,     &
                              shear,      divu,           &
                              e11,        e12,            &
                              e22,                        &
                              s11,        s12,            &
                              s22,                        &
                              yieldstress11,              &
                              yieldstress12,              &
                              yieldstress22,              &
                              prs_sig,                    &
                              rdg_conv,   rdg_shear,      &
                              str)

!
! !DESCRIPTION:
!
! Computes the rates of strain and internal stress components for
! each of the four corners on each T-grid cell.
! Computes stress terms for the momentum equation
! (based on subroutine stress)
!
! !REVISION HISTORY:
!
! same as module
!
! !USES
!
      use ice_constants, only: c0, p027, p055, p111, p166, &
          p2, p222, p25, p333, p5, puny
!
! !INPUT/OUTPUT PARAMETERS:
!
      integer (kind=int_kind), intent(in) :: & 
         nx_block, ny_block, & ! block dimensions
         ksub              , & ! subcycling step
         ndte              , & ! number of subcycles
         icellt                ! no. of cells where icetmask = 1

      integer (kind=int_kind), dimension (nx_block*ny_block), & 
         intent(in) :: &
         indxti   , & ! compressed index in i-direction
         indxtj       ! compressed index in j-direction

      real (kind=dbl_kind), intent(in) :: &
         arlx1i   , & ! dte/2T (original) or 1/alpha1 (revised)
         denom1       ! constant for stress equation

      real (kind=dbl_kind), dimension (nx_block,ny_block), intent(in) :: &
         strength , & ! ice strength (N/m)
         uvel     , & ! x-component of velocity (m/s)
         vvel     , & ! y-component of velocity (m/s)
         dxt      , & ! width of T-cell through the middle (m)
         dyt      , & ! height of T-cell through the middle (m)
         dxhy     , & ! 0.5*(HTE - HTE)
         dyhx     , & ! 0.5*(HTN - HTN)
         cyp      , & ! 1.5*HTE - 0.5*HTE
         cxp      , & ! 1.5*HTN - 0.5*HTN
         cym      , & ! 0.5*HTE - 1.5*HTE
         cxm      , & ! 0.5*HTN - 1.5*HTN
         tarear   , & ! 1/tarea
         tinyarea     ! puny*tarea

      real (kind=dbl_kind), dimension (nx_block,ny_block), & 
         intent(inout) :: &
         stressp_1, stressp_2, stressp_3, stressp_4, & ! sigma11+sigma22
         stressm_1, stressm_2, stressm_3, stressm_4, & ! sigma11-sigma22
         stress12_1,stress12_2,stress12_3,stress12_4   ! sigma12

      real (kind=dbl_kind), dimension (nx_block,ny_block), &
         intent(inout) :: &
         a11, a12, a11_1, a11_2, a11_3, a11_4, & ! structure tensor
         a12_1, a12_2, a12_3, a12_4              ! structure tensor

      real (kind=dbl_kind), dimension (nx_block,ny_block), & 
         intent(inout) :: &
         prs_sig  , & ! replacement pressure, for stress calc
         shear    , & ! strain rate II component (1/s)
         divu     , & ! strain rate I component, velocity divergence (1/s)
         e11      , & ! components of strain rate tensor (1/s)
         e12      , & ! 
         e22      , & ! 
         s11      , & ! components of stress tensor (kg/s^2)
         s12      , & ! 
         s22      , & ! 
         yieldstress11, & ! components of yield stress tensor (kg/s^2)
         yieldstress12, &
         yieldstress22, &
         rdg_conv , & ! convergence term for ridging (1/s)
         rdg_shear    ! shear term for ridging (1/s)

      real (kind=dbl_kind), dimension(nx_block,ny_block,8), & 
         intent(out) :: &
         str          ! stress combinations
!
! local variables
!
      integer (kind=int_kind) :: &
         i, j, ij

      real (kind=dbl_kind), dimension (nx_block,ny_block):: &
         stressptmp_1, stressptmp_2, stressptmp_3, stressptmp_4, & ! sigma11+sigma22
         stressmtmp_1, stressmtmp_2, stressmtmp_3, stressmtmp_4, & ! sigma11-sigma22
         stress12tmp_1,stress12tmp_2,stress12tmp_3,stress12tmp_4   ! sigma12

      real (kind=dbl_kind) :: &
        divune, divunw, divuse, divusw            , & ! divergence
        tensionne, tensionnw, tensionse, tensionsw, & ! tension
        shearne, shearnw, shearse, shearsw        , & ! shearing
        ssigpn, ssigps, ssigpe, ssigpw            , &
        ssigmn, ssigms, ssigme, ssigmw            , &
        ssig12n, ssig12s, ssig12e, ssig12w        , &
        ssigp1, ssigp2, ssigm1, ssigm2, ssig121, ssig122, &
        csigpne, csigpnw, csigpse, csigpsw        , &
        csigmne, csigmnw, csigmse, csigmsw        , &
        csig12ne, csig12nw, csig12se, csig12sw    , &
        str12ew, str12we, str12ns, str12sn        , &
        strp_tmp, strm_tmp

      real (kind=dbl_kind) :: &
        alpharne, alpharnw, alpharsw, alpharse,     &
        alphasne, alphasnw, alphassw, alphasse

      !-----------------------------------------------------------------
      ! Initialize
      !-----------------------------------------------------------------

      str(:,:,:) = c0

!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
      do ij = 1, icellt
         i = indxti(ij)
         j = indxtj(ij)

      !-----------------------------------------------------------------
      ! strain rates
      ! NOTE these are actually strain rates * area  (m^2/s)
      !-----------------------------------------------------------------
         ! divergence  =  e_11 + e_22
         divune    = cyp(i,j)*uvel(i  ,j  ) - dyt(i,j)*uvel(i-1,j  ) &
                   + cxp(i,j)*vvel(i  ,j  ) - dxt(i,j)*vvel(i  ,j-1)
         divunw    = cym(i,j)*uvel(i-1,j  ) + dyt(i,j)*uvel(i  ,j  ) &
                   + cxp(i,j)*vvel(i-1,j  ) - dxt(i,j)*vvel(i-1,j-1)
         divusw    = cym(i,j)*uvel(i-1,j-1) + dyt(i,j)*uvel(i  ,j-1) &
                   + cxm(i,j)*vvel(i-1,j-1) + dxt(i,j)*vvel(i-1,j  )
         divuse    = cyp(i,j)*uvel(i  ,j-1) - dyt(i,j)*uvel(i-1,j-1) &
                   + cxm(i,j)*vvel(i  ,j-1) + dxt(i,j)*vvel(i  ,j  )

         ! tension strain rate  =  e_11 - e_22
         tensionne = -cym(i,j)*uvel(i  ,j  ) - dyt(i,j)*uvel(i-1,j  ) &
                   +  cxm(i,j)*vvel(i  ,j  ) + dxt(i,j)*vvel(i  ,j-1)
         tensionnw = -cyp(i,j)*uvel(i-1,j  ) + dyt(i,j)*uvel(i  ,j  ) &
                   +  cxm(i,j)*vvel(i-1,j  ) + dxt(i,j)*vvel(i-1,j-1)
         tensionsw = -cyp(i,j)*uvel(i-1,j-1) + dyt(i,j)*uvel(i  ,j-1) &
                   +  cxp(i,j)*vvel(i-1,j-1) - dxt(i,j)*vvel(i-1,j  )
         tensionse = -cym(i,j)*uvel(i  ,j-1) - dyt(i,j)*uvel(i-1,j-1) &
                   +  cxp(i,j)*vvel(i  ,j-1) - dxt(i,j)*vvel(i  ,j  )

         ! shearing strain rate  =  e_12
         shearne = -cym(i,j)*vvel(i  ,j  ) - dyt(i,j)*vvel(i-1,j  ) &
                 -  cxm(i,j)*uvel(i  ,j  ) - dxt(i,j)*uvel(i  ,j-1)
         shearnw = -cyp(i,j)*vvel(i-1,j  ) + dyt(i,j)*vvel(i  ,j  ) &
                 -  cxm(i,j)*uvel(i-1,j  ) - dxt(i,j)*uvel(i-1,j-1)
         shearsw = -cyp(i,j)*vvel(i-1,j-1) + dyt(i,j)*vvel(i  ,j-1) &
                 -  cxp(i,j)*uvel(i-1,j-1) + dxt(i,j)*uvel(i-1,j  )
         shearse = -cym(i,j)*vvel(i  ,j-1) - dyt(i,j)*vvel(i-1,j-1) &
                 -  cxp(i,j)*uvel(i  ,j-1) + dxt(i,j)*uvel(i  ,j  )

      !-----------------------------------------------------------------
      ! Stress updated depending on strain rate and structure tensor
      !-----------------------------------------------------------------

         ! ne
         call update_stress_rdg (ksub, ndte, divune, tensionne, &
                                 shearne, a11_1(i,j), a12_1(i,j), &
                                 stressptmp_1(i,j), stressmtmp_1(i,j), &
                                 stress12tmp_1(i,j), strength(i,j), &
                                 alpharne, alphasne)
         ! nw
         call update_stress_rdg (ksub, ndte, divunw, tensionnw, &
                                 shearnw, a11_2(i,j), a12_2(i,j), &
                                 stressptmp_2(i,j), stressmtmp_2(i,j), &
                                 stress12tmp_2(i,j), strength(i,j), &
                                 alpharnw, alphasnw)
         ! sw
         call update_stress_rdg (ksub, ndte, divusw, tensionsw, &
                                 shearsw, a11_3(i,j), a12_3(i,j), &
                                 stressptmp_3(i,j), stressmtmp_3(i,j), &
                                 stress12tmp_3(i,j), strength(i,j), &
                                 alpharsw, alphassw)
         ! se
         call update_stress_rdg (ksub, ndte, divuse, tensionse, &
                                 shearse, a11_4(i,j), a12_4(i,j), &
                                 stressptmp_4(i,j), stressmtmp_4(i,j), &
                                 stress12tmp_4(i,j), strength(i,j), &
                                 alpharse, alphasse)

      !-----------------------------------------------------------------
      ! on last subcycle, save quantities for mechanical redistribution
      !-----------------------------------------------------------------
         if (ksub == ndte) then

            ! diagnostic only
            ! shear = sqrt(tension**2 + shearing**2)
            shear(i,j) = p25*tarear(i,j)*sqrt( &
                 (tensionne + tensionnw + tensionse + tensionsw)**2 &
                +  (shearne +   shearnw +   shearse +   shearsw)**2)

            divu(i,j) = p25*(divune + divunw + divuse + divusw) * tarear(i,j)
            rdg_conv(i,j)  = -min(p25*(alpharne + alpharnw &
                                     + alpharsw + alpharse),c0) * tarear(i,j)
            !rdg_shear=0 for computing closing_net in ridge_prep
            !rdg_shear(i,j) = p25*(alphasne + alphasnw &
            !                    + alphassw + alphasse) * tarear(i,j)
         endif

         e11(i,j) = p5*p25*(divune + divunw + divuse + divusw + &
                    tensionne + tensionnw + tensionse + tensionsw) * tarear(i,j)

         e12(i,j) = p5*p25*(shearne + shearnw + shearse + shearsw) * tarear(i,j)

         e22(i,j) = p5*p25*(divune + divunw + divuse + divusw - &
                    tensionne - tensionnw - tensionse - tensionsw) * tarear(i,j)

         prs_sig(i,j) = strength(i,j) 

      !-----------------------------------------------------------------
      ! elastic relaxation, see Eq. A12-A14
      !-----------------------------------------------------------------

         stressp_1(i,j) = (stressp_1(i,j) + stressptmp_1(i,j)*arlx1i) &
                          * denom1
         stressp_2(i,j) = (stressp_2(i,j) + stressptmp_2(i,j)*arlx1i) &
                          * denom1
         stressp_3(i,j) = (stressp_3(i,j) + stressptmp_3(i,j)*arlx1i) &
                          * denom1
         stressp_4(i,j) = (stressp_4(i,j) + stressptmp_4(i,j)*arlx1i) &
                          * denom1

         stressm_1(i,j) = (stressm_1(i,j) + stressmtmp_1(i,j)*arlx1i) &
                          * denom1
         stressm_2(i,j) = (stressm_2(i,j) + stressmtmp_2(i,j)*arlx1i) &
                          * denom1
         stressm_3(i,j) = (stressm_3(i,j) + stressmtmp_3(i,j)*arlx1i) &
                          * denom1
         stressm_4(i,j) = (stressm_4(i,j) + stressmtmp_4(i,j)*arlx1i) &
                          * denom1

         stress12_1(i,j) = (stress12_1(i,j) + stress12tmp_1(i,j)*arlx1i) &
                          * denom1
         stress12_2(i,j) = (stress12_2(i,j) + stress12tmp_2(i,j)*arlx1i) &
                          * denom1
         stress12_3(i,j) = (stress12_3(i,j) + stress12tmp_3(i,j)*arlx1i) &
                          * denom1
         stress12_4(i,j) = (stress12_4(i,j) + stress12tmp_4(i,j)*arlx1i) &
                          * denom1

          s11(i,j) = p5 * p25 * (stressp_1(i,j) + stressp_2(i,j) &
                               + stressp_3(i,j) + stressp_4(i,j) &
                               + stressm_1(i,j) + stressm_2(i,j) &
                               + stressm_3(i,j) + stressm_4(i,j))
          s22(i,j) = p5 * p25 * (stressp_1(i,j) + stressp_2(i,j) &
                               + stressp_3(i,j) + stressp_4(i,j) &
                               - stressm_1(i,j) - stressm_2(i,j) &
                               - stressm_3(i,j) - stressm_4(i,j))
          s12(i,j) = p25 *      (stress12_1(i,j) + stress12_2(i,j) &
                               + stress12_3(i,j) + stress12_4(i,j))

          yieldstress11(i,j) = p5 * p25 * (stressptmp_1(i,j) + stressptmp_2(i,j) &
                                         + stressptmp_3(i,j) + stressptmp_4(i,j) &
                                         + stressmtmp_1(i,j) + stressmtmp_2(i,j) &
                                         + stressmtmp_3(i,j) + stressmtmp_4(i,j))
          yieldstress22(i,j) = p5 * p25 * (stressptmp_1(i,j) + stressptmp_2(i,j) &
                                         + stressptmp_3(i,j) + stressptmp_4(i,j) &
                                         - stressmtmp_1(i,j) - stressmtmp_2(i,j) &
                                         - stressmtmp_3(i,j) - stressmtmp_4(i,j))
          yieldstress12(i,j) = p25 *      (stress12tmp_1(i,j) + stress12tmp_2(i,j) &
                                         + stress12tmp_3(i,j) + stress12tmp_4(i,j))

      !-----------------------------------------------------------------
      ! Eliminate underflows.
      ! The following code is commented out because it is relatively
      ! expensive and most compilers include a flag that accomplishes
      ! the same thing more efficiently.  This code is cheaper than
      ! handling underflows if the compiler lacks a flag; uncomment
      ! it in that case.  The compiler flag is often described with the
      ! phrase "flush to zero".
      !-----------------------------------------------------------------

!      stressp_1(i,j) = sign(max(abs(stressp_1(i,j)),puny),stressp_1(i,j))
!      stressp_2(i,j) = sign(max(abs(stressp_2(i,j)),puny),stressp_2(i,j))
!      stressp_3(i,j) = sign(max(abs(stressp_3(i,j)),puny),stressp_3(i,j))
!      stressp_4(i,j) = sign(max(abs(stressp_4(i,j)),puny),stressp_4(i,j))

!      stressm_1(i,j) = sign(max(abs(stressm_1(i,j)),puny),stressm_1(i,j))
!      stressm_2(i,j) = sign(max(abs(stressm_2(i,j)),puny),stressm_2(i,j))
!      stressm_3(i,j) = sign(max(abs(stressm_3(i,j)),puny),stressm_3(i,j))
!      stressm_4(i,j) = sign(max(abs(stressm_4(i,j)),puny),stressm_4(i,j))

!      stress12_1(i,j) = sign(max(abs(stress12_1(i,j)),puny),stress12_1(i,j))
!      stress12_2(i,j) = sign(max(abs(stress12_2(i,j)),puny),stress12_2(i,j))
!      stress12_3(i,j) = sign(max(abs(stress12_3(i,j)),puny),stress12_3(i,j))
!      stress12_4(i,j) = sign(max(abs(stress12_4(i,j)),puny),stress12_4(i,j))

      !-----------------------------------------------------------------
      ! combinations of the stresses for the momentum equation ! kg/s^2
      !-----------------------------------------------------------------

         ssigpn  = stressp_1(i,j) + stressp_2(i,j)
         ssigps  = stressp_3(i,j) + stressp_4(i,j)
         ssigpe  = stressp_1(i,j) + stressp_4(i,j)
         ssigpw  = stressp_2(i,j) + stressp_3(i,j)
         ssigp1  =(stressp_1(i,j) + stressp_3(i,j))*p055
         ssigp2  =(stressp_2(i,j) + stressp_4(i,j))*p055

         ssigmn  = stressm_1(i,j) + stressm_2(i,j)
         ssigms  = stressm_3(i,j) + stressm_4(i,j)
         ssigme  = stressm_1(i,j) + stressm_4(i,j)
         ssigmw  = stressm_2(i,j) + stressm_3(i,j)
         ssigm1  =(stressm_1(i,j) + stressm_3(i,j))*p055
         ssigm2  =(stressm_2(i,j) + stressm_4(i,j))*p055

         ssig12n = stress12_1(i,j) + stress12_2(i,j)
         ssig12s = stress12_3(i,j) + stress12_4(i,j)
         ssig12e = stress12_1(i,j) + stress12_4(i,j)
         ssig12w = stress12_2(i,j) + stress12_3(i,j)
         ssig121 =(stress12_1(i,j) + stress12_3(i,j))*p111
         ssig122 =(stress12_2(i,j) + stress12_4(i,j))*p111

         csigpne = p111*stressp_1(i,j) + ssigp2 + p027*stressp_3(i,j)
         csigpnw = p111*stressp_2(i,j) + ssigp1 + p027*stressp_4(i,j)
         csigpsw = p111*stressp_3(i,j) + ssigp2 + p027*stressp_1(i,j)
         csigpse = p111*stressp_4(i,j) + ssigp1 + p027*stressp_2(i,j)
         
         csigmne = p111*stressm_1(i,j) + ssigm2 + p027*stressm_3(i,j)
         csigmnw = p111*stressm_2(i,j) + ssigm1 + p027*stressm_4(i,j)
         csigmsw = p111*stressm_3(i,j) + ssigm2 + p027*stressm_1(i,j)
         csigmse = p111*stressm_4(i,j) + ssigm1 + p027*stressm_2(i,j)
         
         csig12ne = p222*stress12_1(i,j) + ssig122 &
                  + p055*stress12_3(i,j)
         csig12nw = p222*stress12_2(i,j) + ssig121 &
                  + p055*stress12_4(i,j)
         csig12sw = p222*stress12_3(i,j) + ssig122 &
                  + p055*stress12_1(i,j)
         csig12se = p222*stress12_4(i,j) + ssig121 &
                  + p055*stress12_2(i,j)

         str12ew = p5*dxt(i,j)*(p333*ssig12e + p166*ssig12w)
         str12we = p5*dxt(i,j)*(p333*ssig12w + p166*ssig12e)
         str12ns = p5*dyt(i,j)*(p333*ssig12n + p166*ssig12s)
         str12sn = p5*dyt(i,j)*(p333*ssig12s + p166*ssig12n)

      !-----------------------------------------------------------------
      ! for dF/dx (u momentum)
      !-----------------------------------------------------------------
         strp_tmp  = p25*dyt(i,j)*(p333*ssigpn  + p166*ssigps)
         strm_tmp  = p25*dyt(i,j)*(p333*ssigmn  + p166*ssigms)

         ! northeast (i,j)
         str(i,j,1) = -strp_tmp - strm_tmp - str12ew &
              + dxhy(i,j)*(-csigpne + csigmne) + dyhx(i,j)*csig12ne

         ! northwest (i+1,j)
         str(i,j,2) = strp_tmp + strm_tmp - str12we &
              + dxhy(i,j)*(-csigpnw + csigmnw) + dyhx(i,j)*csig12nw

         strp_tmp  = p25*dyt(i,j)*(p333*ssigps  + p166*ssigpn)
         strm_tmp  = p25*dyt(i,j)*(p333*ssigms  + p166*ssigmn)

         ! southeast (i,j+1)
         str(i,j,3) = -strp_tmp - strm_tmp + str12ew &
              + dxhy(i,j)*(-csigpse + csigmse) + dyhx(i,j)*csig12se

         ! southwest (i+1,j+1)
         str(i,j,4) = strp_tmp + strm_tmp + str12we &
              + dxhy(i,j)*(-csigpsw + csigmsw) + dyhx(i,j)*csig12sw

      !-----------------------------------------------------------------
      ! for dF/dy (v momentum)
      !-----------------------------------------------------------------
         strp_tmp  = p25*dxt(i,j)*(p333*ssigpe  + p166*ssigpw)
         strm_tmp  = p25*dxt(i,j)*(p333*ssigme  + p166*ssigmw)

         ! northeast (i,j)
         str(i,j,5) = -strp_tmp + strm_tmp - str12ns &
              - dyhx(i,j)*(csigpne + csigmne) + dxhy(i,j)*csig12ne

         ! southeast (i,j+1)
         str(i,j,6) = strp_tmp - strm_tmp - str12sn &
              - dyhx(i,j)*(csigpse + csigmse) + dxhy(i,j)*csig12se

         strp_tmp  = p25*dxt(i,j)*(p333*ssigpw  + p166*ssigpe)
         strm_tmp  = p25*dxt(i,j)*(p333*ssigmw  + p166*ssigme)

         ! northwest (i+1,j)
         str(i,j,7) = -strp_tmp + strm_tmp + str12ns &
              - dyhx(i,j)*(csigpnw + csigmnw) + dxhy(i,j)*csig12nw

         ! southwest (i+1,j+1)
         str(i,j,8) = strp_tmp - strm_tmp + str12sn &
              - dyhx(i,j)*(csigpsw + csigmsw) + dxhy(i,j)*csig12sw

      enddo                     ! ij

      end subroutine stress_eap

!=======================================================================
!BOP
!
! !IROUTINE: update_stress_rdg
!
! !INTERFACE:
!
      subroutine update_stress_rdg (ksub, ndte, divu, tension, &
                                   shear, a11, a12, &
                                   stressp,  stressm, &
                                   stress12, strength, &
                                   alphar, alphas)
!
! !DESCRIPTION:
!
! Updates the stress depending on values of strain rate and structure
! tensor and for ksub=ndte it computes closing and sliding rate
!
! !REVISION HISTORY:
!
! same as module

      use ice_constants, only: c0, p025, p05, p1, p5, c1, c2, c12, puny, &
          pi, pih, pi2, piq
!
! !INPUT/OUTPUT PARAMETERS:
!
      integer (kind=int_kind), intent(in) :: &
         ksub, &
         ndte

      real (kind=dbl_kind), intent(in) :: &
         a11,     a12,                   &
         divu,    tension, shear,        &
         strength

      real (kind=dbl_kind), intent(out) :: &
         stressp, stressm, stress12, &
         alphar, alphas     
!
!     local variables
!
      real (kind=dbl_kind), dimension(2,2) :: &
         Q, Qd, atemp,                   &
         dtemp,               &
         sig, sigprime

      real (kind=dbl_kind) :: &
         stemp11r, stemp12r, stemp22r,   &
         stemp11s, stemp12s, stemp22s

      real (kind=dbl_kind) :: &
         rotstemp11r, rotstemp12r, rotstemp22r,   &
         rotstemp11s, rotstemp12s, rotstemp22s

      real (kind=dbl_kind) :: &
         gamma, alpha, x, y, xmin, ymin, dx, dy, da, invdx, invdy, invda, invsin

      real (kind=dbl_kind) :: &
         invleng, dtemp1, dtemp2, atempprime

      real (kind=dbl_kind) :: &
         kx ,ky, ka, fx, fy, fa

      real (kind=dbl_kind) :: &
	 invstressconviso

      real (kind=dbl_kind), parameter :: &
         kfriction = 0.45_dbl_kind

! Factor to maintain the same stress as in EVP (see Section 3)
! Can be set to 1 otherwise
      invstressconviso = c1/(c1+kfriction*kfriction)

! compute eigenvalues, eigenvectors and angles for structure tensor, strain rates

! 1) structure tensor

         atemp(1,1) = a11
         atemp(1,2) = a12
!         atemp(2,1) = a12
         atemp(2,2) = 1-a11

! gamma: angle between general coordiantes and principal axis 
         gamma = p5*atan2((c2*atemp(1,2)),(atemp(1,1) - atemp(2,2)))

! rotational tensor from general coordinates into principal axis
         Q(1,1) = cos(gamma)
         Q(1,2) = sin(gamma)
!         Q(2,1) = -Q(1,2)
         Q(2,2) = Q(1,1)

! rotation Q*atemp*Q^T
         atempprime = Q(1,1)*(Q(1,1)*atemp(1,1) + 2*Q(1,2)*atemp(1,2)) &
                      + Q(1,2)*Q(1,2)*atemp(2,2) 

! make first principal value the largest
         if (atempprime < p5) atempprime = c1 - atempprime

! 2) strain rate

         dtemp(1,1) = p5*(divu + tension)
         dtemp(1,2) = shear*p5
!         dtemp(2,1) = shear*p5
         dtemp(2,2) = p5*(divu - tension)

! alpha: angle between general coordiantes and principal axis
         alpha = p5*atan2((c2*dtemp(1,2)),(dtemp(1,1) - dtemp(2,2)))

! y: angle between major principal axis of strain rate and structure tensor
! to make sure y between 0 and pi/2
	 if (alpha > gamma) alpha = alpha - pi
	 if (alpha < gamma-pi) alpha = alpha + pi

         y = gamma - alpha

! rotational tensor (anticlockwise) from general coordinates into principal axis
         Qd(1,1) = cos(alpha)
         Qd(1,2) = sin(alpha)
!         Qd(2,1) = -Qd(1,2)
         Qd(2,2) = Qd(1,1)
         
         dtemp1 = Qd(1,1)*(Qd(1,1)*dtemp(1,1) + 2*Qd(1,2)*dtemp(1,2)) &
                           + Qd(1,2)*Qd(1,2)*dtemp(2,2)
         dtemp2 = Qd(1,2)*(Qd(1,2)*dtemp(1,1) - 2*Qd(1,1)*dtemp(1,2)) &
                           + Qd(1,1)*Qd(1,1)*dtemp(2,2)

! In cos and sin values
         if ((ABS(dtemp1) > puny).or.(ABS(dtemp2) > puny)) then
           invleng = c1/sqrt(dtemp1*dtemp1 + dtemp2*dtemp2)
           dtemp1 = dtemp1*invleng
           dtemp2 = dtemp2*invleng
           x = atan2(dtemp2,dtemp1)
         else
           x = c0
         endif

! to ensure the angle lies between pi/4 and 9 pi/4 
         if (x < p5*pih) x = x + pi2

! Now calculate updated stress tensor
         xmin = piq+pi
         ymin = c0
         dx   = pi2*p025
         dy   = pih*p1
         da   = p05
         invdx = c1/dx
         invdy = c1/dy
         invda = c1/da
         invsin = c1/sin(2*pi/c12)

         kx = int((x-xmin)*invdx) + 1
         ky = int((y-ymin)*invdy) + 1
         ka = int((atempprime-p5)*invda) + 1
         fx = (x - xmin) - (kx - 1)*dx
         fy = (y - ymin) - (ky - 1)*dy
         fa = (atempprime - p5) - (ka - 1)*da

! Determine sigma_r(A1,Zeta,y) by interpolation in kx,ky,ka space i
! (see Section A1) 

         stemp11r = s11r(kx,ky,ka)*(1-fx*invdx)*(1-fy*invdy)*(1-fa*invda) &     
                    + s11r(kx+1,ky,ka)*fx*invdx*(1-fy*invdy)*(1-fa*invda) &
                    + s11r(kx,ky+1,ka)*(1-fx*invdx)*fy*invdy*(1-fa*invda) &
                    + s11r(kx,ky,ka+1)*(1-fx*invdx)*(1-fy*invdy)*fa*invda &
                    + s11r(kx+1,ky,ka+1)*fx*invdx*(1-fy*invdy)*fa*invda &
                    + s11r(kx,ky+1,ka+1)*(1-fx*invdx)*fy*invdy*fa*invda &
                    + s11r(kx+1,ky+1,ka)*fx*invdx*fy*invdy*(1-fa*invda) &
                    + s11r(kx+1,ky+1,ka+1)*fx*invdx*fy*invdy*fa*invda

         stemp12r = s12r(kx,ky,ka)*(1-fx*invdx)*(1-fy*invdy)*(1-fa*invda) &
                    + s12r(kx+1,ky,ka)*fx*invdx*(1-fy*invdy)*(1-fa*invda) &
                    + s12r(kx,ky+1,ka)*(1-fx*invdx)*fy*invdy*(1-fa*invda) &
                    + s12r(kx,ky,ka+1)*(1-fx*invdx)*(1-fy*invdy)*fa*invda &
                    + s12r(kx+1,ky,ka+1)*fx*invdx*(1-fy*invdy)*fa*invda &
                    + s12r(kx,ky+1,ka+1)*(1-fx*invdx)*fy*invdy*fa*invda &
                    + s12r(kx+1,ky+1,ka)*fx*invdx*fy*invdy*(1-fa*invda) &
                    + s12r(kx+1,ky+1,ka+1)*fx*invdx*fy*invdy*fa*invda

         stemp22r = s22r(kx,ky,ka)*(1-fx*invdx)*(1-fy*invdy)*(1-fa*invda) &
                    + s22r(kx+1,ky,ka)*fx*invdx*(1-fy*invdy)*(1-fa*invda) &
                    + s22r(kx,ky+1,ka)*(1-fx*invdx)*fy*invdy*(1-fa*invda) &
                    + s22r(kx,ky,ka+1)*(1-fx*invdx)*(1-fy*invdy)*fa*invda &
                    + s22r(kx+1,ky,ka+1)*fx*invdx*(1-fy*invdy)*fa*invda &
                    + s22r(kx,ky+1,ka+1)*(1-fx*invdx)*fy*invdy*fa*invda &
                    + s22r(kx+1,ky+1,ka)*fx*invdx*fy*invdy*(1-fa*invda) &
                    + s22r(kx+1,ky+1,ka+1)*fx*invdx*fy*invdy*fa*invda

! Determine sigma_s(A1,Zeta,y) by interpolation in kx,ky,ka space (see Eq. A1)

         stemp11s = s11s(kx,ky,ka)*(1-fx*invdx)*(1-fy*invdy)*(1-fa*invda) &
                    + s11s(kx+1,ky,ka)*fx*invdx*(1-fy*invdy)*(1-fa*invda) &
                    + s11s(kx,ky+1,ka)*(1-fx*invdx)*fy*invdy*(1-fa*invda) &
                    + s11s(kx,ky,ka+1)*(1-fx*invdx)*(1-fy*invdy)*fa*invda &
                    + s11s(kx+1,ky,ka+1)*fx*invdx*(1-fy*invdy)*fa*invda &
                    + s11s(kx,ky+1,ka+1)*(1-fx*invdx)*fy*invdy*fa*invda &
                    + s11s(kx+1,ky+1,ka)*fx*invdx*fy*invdy*(1-fa*invda) &
                    + s11s(kx+1,ky+1,ka+1)*fx*invdx*fy*invdy*fa*invda

         stemp12s = s12s(kx,ky,ka)*(1-fx*invdx)*(1-fy*invdy)*(1-fa*invda) &
                    + s12s(kx+1,ky,ka)*fx*invdx*(1-fy*invdy)*(1-fa*invda) &
                    + s12s(kx,ky+1,ka)*(1-fx*invdx)*fy*invdy*(1-fa*invda) &
                    + s12s(kx,ky,ka+1)*(1-fx*invdx)*(1-fy*invdy)*fa*invda &
                    + s12s(kx+1,ky,ka+1)*fx*invdx*(1-fy*invdy)*fa*invda &
                    + s12s(kx,ky+1,ka+1)*(1-fx*invdx)*fy*invdy*fa*invda &
                    + s12s(kx+1,ky+1,ka)*fx*invdx*fy*invdy*(1-fa*invda) &
                    + s12s(kx+1,ky+1,ka+1)*fx*invdx*fy*invdy*fa*invda

         stemp22s = s22s(kx,ky,ka)*(1-fx*invdx)*(1-fy*invdy)*(1-fa*invda) &
                    + s22s(kx+1,ky,ka)*fx*invdx*(1-fy*invdy)*(1-fa*invda) &
                    + s22s(kx,ky+1,ka)*(1-fx*invdx)*fy*invdy*(1-fa*invda) &
                    + s22s(kx,ky,ka+1)*(1-fx*invdx)*(1-fy*invdy)*fa*invda &
                    + s22s(kx+1,ky,ka+1)*fx*invdx*(1-fy*invdy)*fa*invda &
                    + s22s(kx,ky+1,ka+1)*(1-fx*invdx)*fy*invdy*fa*invda &
                    + s22s(kx+1,ky+1,ka)*fx*invdx*fy*invdy*(1-fa*invda) &
                    + s22s(kx+1,ky+1,ka+1)*fx*invdx*fy*invdy*fa*invda

! Calculate mean ice stress over a collection of floes (Equation 3)

         stressp  = strength*(stemp11r + kfriction*stemp11s + stemp22r &
                    + kfriction*stemp22s) *invsin *invstressconviso
         stress12 = strength*(stemp12r + kfriction*stemp12s) &
                    *invsin *invstressconviso
         stressm  = strength*(stemp11r + kfriction*stemp11s - stemp22r &
                    - kfriction*stemp22s)*invsin *invstressconviso

! Back - rotation of the stress from principal aces into general coordinates

! Update stress
         sig(1,1) = p5*(stressp + stressm)
         sig(1,2) = stress12
!         sig(2,1) = stress12
         sig(2,2) = p5*(stressp - stressm)

         sigprime(1,1) = Q(1,1)*(Q(1,1)*sig(1,1) - 2*Q(1,2)*sig(1,2)) &
                         + Q(1,2)*Q(1,2)*sig(2,2)
         sigprime(1,2) = Q(1,1)*(Q(1,1)*sig(1,2) + Q(1,2)*sig(1,1) &
                         - Q(1,2)*sig(2,2)) - Q(1,2)*Q(1,2)*sig(1,2)
         sigprime(2,2) = Q(1,2)*(Q(1,2)*sig(1,1) + 2*Q(1,1)*sig(1,2)) &
                         + Q(1,1)*Q(1,1)*sig(2,2)

         stressp  = sigprime(1,1) + sigprime(2,2)
         stress12 = sigprime(1,2)
         stressm  = sigprime(1,1) - sigprime(2,2)

! Compute ridging and sliding functions in general coordinates (Equation 11)
         if (ksub == ndte) then

           rotstemp11r = Q(1,1)*(Q(1,1)*stemp11r-2*Q(1,2)*stemp12r) &
                         + Q(1,2)*Q(1,2)*stemp22r
           rotstemp12r = Q(1,1)*(Q(1,1)*stemp12r+Q(1,2)*stemp11r-Q(1,2)*stemp22r) &
                         - Q(1,2)*Q(1,2)*stemp12r
           rotstemp22r = Q(1,2)*(Q(1,2)*stemp11r+2*Q(1,1)*stemp12r) & 
                         + Q(1,1)*Q(1,1)*stemp22r

           rotstemp11s = Q(1,1)*(Q(1,1)*stemp11s-2*Q(1,2)*stemp12s) &
                         + Q(1,2)*Q(1,2)*stemp22s
           rotstemp12s = Q(1,1)*(Q(1,1)*stemp12s+Q(1,2)*stemp11s-Q(1,2)*stemp22s) &
                         - Q(1,2)*Q(1,2)*stemp12s
           rotstemp22s = Q(1,2)*(Q(1,2)*stemp11s+2*Q(1,1)*stemp12s) & 
                         + Q(1,1)*Q(1,1)*stemp22s

           alphar =  rotstemp11r*dtemp(1,1) + 2*rotstemp12r*dtemp(1,2) &
                     + rotstemp22r*dtemp(2,2)
           alphas =  rotstemp11s*dtemp(1,1) + 2*rotstemp12s*dtemp(1,2) &
                     + rotstemp22s*dtemp(2,2)

         endif

      end subroutine update_stress_rdg

!=======================================================================
!BOP
!
! !IROUTINE: stepa - computes structure tensor
!
! !INTERFACE:
!
      subroutine stepa  (nx_block,   ny_block,       &
                         dtei,       icellt,         &
                         indxti,     indxtj,         &
                         a11, a12,                   &
                         a11_1, a11_2, a11_3, a11_4, &
                         a12_1, a12_2, a12_3, a12_4, &
                         stressp_1,  stressp_2,      &
                         stressp_3,  stressp_4,      &
                         stressm_1,  stressm_2,      &
                         stressm_3,  stressm_4,      &
                         stress12_1, stress12_2,     &
                         stress12_3, stress12_4)
! !DESCRIPTION:
!
! Solves evolution equation for structure tensor (A19, A20) 
!
! !REVISION HISTORY:
!
! same as module
!
! !USES
!
      use ice_constants, only: p001, p2, p25, p5
!
! !INPUT/OUTPUT PARAMETERS:
!
      integer (kind=int_kind), intent(in) :: &
         nx_block, ny_block, & ! block dimensions
         icellt                ! no. of cells where icetmask = 1

      real (kind=dbl_kind), intent(in) :: &
         dtei        ! 1/dte, where dte is subcycling timestep (1/s)

      integer (kind=int_kind), dimension (nx_block*ny_block), &
         intent(in) :: &
         indxti   , & ! compressed index in i-direction
         indxtj       ! compressed index in j-direction

      real (kind=dbl_kind), dimension (nx_block,ny_block), &
         intent(in) :: &
         ! ice stress tensor (kg/s^2) in each corner of T cell
         stressp_1, stressp_2, stressp_3, stressp_4, & ! sigma11+sigma22
         stressm_1, stressm_2, stressm_3, stressm_4, & ! sigma11-sigma22
         stress12_1, stress12_2, stress12_3, stress12_4    ! sigma12

      real (kind=dbl_kind), dimension (nx_block,ny_block), &
         intent(inout) :: &
         ! structure tensor () in each corner of T cell
         a11, a12, a11_1, a11_2, a11_3, a11_4, & ! components of 
         a12_1, a12_2, a12_3, a12_4              ! structure tensor ()
!
! local variables
!
      integer (kind=int_kind) :: &
         i, j, ij

      real (kind=dbl_kind) :: &
        mresult11, mresult12

      real (kind=dbl_kind), parameter :: &
        kth  = p2*p001             

      do ij = 1, icellt
         i = indxti(ij)
         j = indxtj(ij)

! ne 
         call calc_ffrac(1, stressp_1(i,j), stressm_1(i,j),  &
                           stress12_1(i,j),                  &
                           a11_1(i,j),                       &
                           mresult11)

         call calc_ffrac(2, stressp_1(i,j), stressm_1(i,j),  &
                           stress12_1(i,j),                  &
                           a12_1(i,j),                       &
                           mresult12)

         a11_1(i,j) = (a11_1(i,j)*dtei + kth*p5 - mresult11) / &
                      (dtei + kth) ! implicit
         a12_1(i,j) = (a12_1(i,j)*dtei - mresult12)/(dtei + kth) ! implicit

 
! nw
         call calc_ffrac(1, stressp_2(i,j), stressm_2(i,j),  &
                           stress12_2(i,j),                  &
                           a11_2(i,j),                       &
                           mresult11)

         call calc_ffrac(2, stressp_2(i,j), stressm_2(i,j),  &
                           stress12_2(i,j),                  &
                           a12_2(i,j),                       &
                           mresult12)

         a11_2(i,j) = (a11_2(i,j)*dtei + kth*p5 - mresult11) / &
                      (dtei + kth) ! implicit
         a12_2(i,j) = (a12_2(i,j)*dtei - mresult12)/(dtei + kth) ! implicit

! sw
         call calc_ffrac(1, stressp_3(i,j), stressm_3(i,j),  &
                           stress12_3(i,j),                  &
                           a11_3(i,j),                       &
                           mresult11)

         call calc_ffrac(2, stressp_3(i,j), stressm_3(i,j),  &
                           stress12_3(i,j),                  &
                           a12_3(i,j),                       &
                           mresult12)

         a11_3(i,j) = (a11_3(i,j)*dtei + kth*p5 - mresult11) / &
                      (dtei + kth) ! implicit
         a12_3(i,j) = (a12_3(i,j)*dtei - mresult12)/(dtei + kth) ! implicit
                                
! se
         call calc_ffrac(1, stressp_4(i,j), stressm_4(i,j),  &
                           stress12_4(i,j),                  &
                           a11_4(i,j),                       &
                           mresult11)

         call calc_ffrac(2, stressp_4(i,j), stressm_4(i,j),  &
                           stress12_4(i,j),                  &
                           a12_4(i,j),                       &
                           mresult12)

         a11_4(i,j) = (a11_4(i,j)*dtei + kth*p5 - mresult11) / &
                      (dtei + kth) ! implicit
         a12_4(i,j) = (a12_4(i,j)*dtei - mresult12)/(dtei + kth) ! implicit

! average structure tensor

         a11(i,j) = p25*(a11_1(i,j) + a11_2(i,j) + a11_3(i,j) + a11_4(i,j))
         a12(i,j) = p25*(a12_1(i,j) + a12_2(i,j) + a12_3(i,j) + a12_4(i,j))
               
      enddo                     ! ij
      
      end subroutine stepa

!=======================================================================
!BOP
!
! !IROUTINE: calc_ffrac 
!
! !INTERFACE:
!
      subroutine calc_ffrac (blockno, stressp, stressm, &
                             stress12,                  &
                             a1x,                       &
                             mresult)
!
! !DESCRIPTION:
!
! computes term in evolution equation for structure tensor which determines
! the ice floe re-orientation due to fracture
! Eq. 7: Ffrac = -kf(A-S) or = 0 depending on sigma_1 and sigma_2
!
! !REVISION HISTORY:
!
! same as module
!
! !USES:
!
      use ice_constants, only: c0, p001, p1, p5, c2, c3
!
! !INPUT/OUTPUT PARAMETERS:
!
      integer(kind=int_kind), intent(in) :: &
         blockno

      real (kind=dbl_kind), intent(in) :: &
         stressp, stressm, stress12, a1x

      real (kind=dbl_kind), intent(out) :: &
         mresult
!
! local variables
!
      real (kind=dbl_kind), dimension(2,2) :: &
         Q, sigma

      real (kind=dbl_kind) :: &
         gamma, sigma_1, sigma_2, m

      real (kind=dbl_kind), parameter :: &
         kfrac = p001, &
         threshold = c3*p1

       sigma(1,1) = p5*(stressp+stressm) 
       sigma(1,2) = stress12 
       sigma(2,1) = sigma(1,2) 
       sigma(2,2) = p5*(stressp-stressm) 

       gamma = p5*atan2((c2*sigma(1,2)),(sigma(1,1)-sigma(2,2)))

! rotational tensor to get into sigma principal axis

       Q(1,1) = cos(gamma)
       Q(1,2) = sin(gamma)
       Q(2,1) = -Q(1,2)
       Q(2,2) = Q(1,1)

       sigma_1 = Q(1,1)*(Q(1,1)*sigma(1,1)+2*Q(1,2)*sigma(1,2)) + &
            Q(1,2)*Q(1,2)*sigma(2,2) ! S(1,1)
       sigma_2 = Q(1,2)*(Q(1,2)*sigma(1,1)-2*Q(1,1)*sigma(1,2)) + &
            Q(1,1)*Q(1,1)*sigma(2,2) ! S(2,2)

! Pure divergence
       if ((sigma_1 >= c0).and.(sigma_2 >= c0))  then
         mresult = c0

! Unconfined compression: cracking of blocks not along the axial splitting direction
! which leads to the loss of their shape, so we again model it through diffusion
       elseif ((sigma_1 >= c0).and.(sigma_2 < c0))  then
         if (blockno == 1) mresult = kfrac * (a1x - Q(1,2)*Q(1,2))
         if (blockno == 2) mresult = kfrac * (a1x + Q(1,1)*Q(1,2))

! Shear faulting
       elseif (sigma_2 == c0) then
         mresult =c0
       elseif ((sigma_1 <= c0).and.(sigma_1/sigma_2 <= threshold)) then
         if (blockno == 1) mresult = kfrac * (a1x - Q(1,2)*Q(1,2))
         if (blockno == 2) mresult = kfrac * (a1x + Q(1,1)*Q(1,2))

! Horizontal spalling
       else 
         mresult = c0
       endif

      end subroutine calc_ffrac

!=======================================================================
!---! these subroutines write/read Fortran unformatted data files ..
!=======================================================================
!
!BOP
!
! !IROUTINE: dumpfile - dumps all fields required for restart
!
! !INTERFACE:
!
      subroutine write_restart_eap (filename_spec)
!
! !DESCRIPTION:
!
! Dumps all values needed for a restart
!
! !REVISION HISTORY:
!
! same as module
!
! !USES:
!
      use ice_calendar, only: sec, month, mday, nyr, istep1, &
                              time, time_forc, idate, year_init
      use ice_communicate, only: my_task, master_task
      use ice_fileunits, only: nu_diag, nu_dump_eap
      use ice_read_write, only: ice_open, ice_write
      use ice_restart, only: lenstr, restart_dir, restart_file
!
! !INPUT/OUTPUT PARAMETERS:
!
      character(len=char_len_long), intent(in), optional :: filename_spec

!EOP
!
      integer (kind=int_kind) :: &
          i, j, k, n, it, iblk, & ! counting indices
          iyear, imonth, iday     ! year, month, day

      character(len=char_len_long) :: filename

      logical (kind=log_kind) :: diag

      ! construct path/file
      if (present(filename_spec)) then
         filename = trim(filename_spec)
      else
         iyear = nyr + year_init - 1
         imonth = month
         iday = mday
         
         write(filename,'(a,a,a,i4.4,a,i2.2,a,i2.2,a,i5.5)') &
              restart_dir(1:lenstr(restart_dir)), &
              restart_file(1:lenstr(restart_file)),'.eap.', &
              iyear,'-',month,'-',mday,'-',sec
      end if
         
      ! begin writing restart data
      call ice_open(nu_dump_eap,filename,0)

      if (my_task == master_task) then
        write(nu_dump_eap) istep1,time,time_forc
        write(nu_diag,*) 'Writing ',filename(1:lenstr(filename))
      endif

      diag = .true.

      !-----------------------------------------------------------------
      ! structure tensor 
      !-----------------------------------------------------------------
      call ice_write(nu_dump_eap,0,a11_1,'ruf8',diag)
      call ice_write(nu_dump_eap,0,a11_3,'ruf8',diag)
      call ice_write(nu_dump_eap,0,a11_2,'ruf8',diag)
      call ice_write(nu_dump_eap,0,a11_4,'ruf8',diag)

      call ice_write(nu_dump_eap,0,a12_1,'ruf8',diag)
      call ice_write(nu_dump_eap,0,a12_3,'ruf8',diag)
      call ice_write(nu_dump_eap,0,a12_2,'ruf8',diag)
      call ice_write(nu_dump_eap,0,a12_4,'ruf8',diag)

      if (my_task == master_task) close(nu_dump_eap)

      end subroutine write_restart_eap

!=======================================================================
!BOP
!
! !IROUTINE: read_restart_eap - reads all fields required for restart
!
! !INTERFACE:
!
      subroutine read_restart_eap(filename_spec)
!
! !DESCRIPTION:
!
! Reads all values needed for elastic anisotropic plastic dynamics restart
!
! !REVISION HISTORY:
!
! same as module
!
! !USES:
! 
      use ice_blocks, only: nghost
      use ice_calendar, only: istep1, time, time_forc
      use ice_communicate, only: my_task, master_task
      use ice_constants, only: c0
      use ice_domain, only: distrb_info, nblocks
      use ice_domain_size, only: nx_global, ny_global
      use ice_exit, only: abort_ice
      use ice_gather_scatter, only: scatter_global_stress
      use ice_read_write, only: ice_open, ice_read_global
      use ice_restart, only: lenstr, restart_file, &
                             pointer_file, runtype
      use ice_fileunits, only: nu_diag, nu_rst_pointer, nu_restart_eap
!
! !INPUT/OUTPUT PARAMETERS:
!
      character(len=char_len_long), intent(in), optional :: filename_spec

!EOP
!
      integer (kind=int_kind) :: &
          i, j, k, n, it, iblk ! counting indices

      character(len=char_len_long) :: &
         filename, filename0, string1, string2

      logical (kind=log_kind) :: &
         diag

      real (kind=dbl_kind), dimension(:,:), allocatable :: &
         work_g1, work_g2

      if (my_task == master_task) then
         open(nu_rst_pointer,file=pointer_file)
         read(nu_rst_pointer,'(a)') filename0
         filename = trim(filename0)
         close(nu_rst_pointer)

         ! reconstruct path/file
         n = index(filename0,trim(restart_file))
         if (n == 0) call abort_ice('eap restart: filename discrepancy')
         string1 = trim(filename0(1:n-1))
         string2 = trim(filename0(n+lenstr(restart_file):lenstr(filename0)))
         write(filename,'(a,a,a,a)') &
            string1(1:lenstr(string1)), &
            restart_file(1:lenstr(restart_file)),'.eap', &
            string2(1:lenstr(string2))
      endif ! master_task

      if (runtype == 'continue') then 

      call ice_open(nu_restart_eap,filename,0)

      if (my_task == master_task) then
        read(nu_restart_eap) istep1,time,time_forc
        write(nu_diag,*) 'Reading ',filename(1:lenstr(filename))
      endif

      diag = .true.

      !-----------------------------------------------------------------
      ! Structure tensor must be read and scattered in pairs in order
      ! to properly match corner values across a tripole grid cut.
      !-----------------------------------------------------------------
         if (my_task == master_task) write(nu_diag,*) &
           'structure tensor'
    
         allocate (work_g1(nx_global,ny_global), &
                work_g2(nx_global,ny_global))

         call ice_read_global(nu_restart_eap,0,work_g1,'ruf8',diag) ! a11_1
         call ice_read_global(nu_restart_eap,0,work_g2,'ruf8',diag) ! a11_3
         call scatter_global_stress(a11_1, work_g1, work_g2, &
                                 master_task, distrb_info)
         call scatter_global_stress(a11_3, work_g2, work_g1, &
                                 master_task, distrb_info)

         call ice_read_global(nu_restart_eap,0,work_g1,'ruf8',diag) ! a11_2
         call ice_read_global(nu_restart_eap,0,work_g2,'ruf8',diag) ! a11_4
         call scatter_global_stress(a11_2, work_g1, work_g2, &
                                 master_task, distrb_info)
         call scatter_global_stress(a11_4, work_g2, work_g1, &
                                 master_task, distrb_info)

         call ice_read_global(nu_restart_eap,0,work_g1,'ruf8',diag) ! a12_1
         call ice_read_global(nu_restart_eap,0,work_g2,'ruf8',diag) ! a12_3
         call scatter_global_stress(a12_1, work_g1, work_g2, &
                                 master_task, distrb_info)
         call scatter_global_stress(a12_3, work_g2, work_g1, &
                                 master_task, distrb_info)

         call ice_read_global(nu_restart_eap,0,work_g1,'ruf8',diag) ! a12_2
         call ice_read_global(nu_restart_eap,0,work_g2,'ruf8',diag) ! a12_4
         call scatter_global_stress(a12_2, work_g1, work_g2, &
                                 master_task, distrb_info)
         call scatter_global_stress(a12_4, work_g2, work_g1, &
                                 master_task, distrb_info)

         deallocate (work_g1, work_g2)

      if (my_task == master_task) close(nu_restart_eap)
      endif

      !-----------------------------------------------------------------
      ! Ensure unused values in west and south ghost cells are 0
      !-----------------------------------------------------------------

         !$OMP PARALLEL DO PRIVATE(iblk,i,j)
         do iblk = 1, nblocks
            do j = 1, nghost
            do i = 1, nx_block
               a11_1 (i,j,iblk) = c0
               a11_2 (i,j,iblk) = c0
               a11_3 (i,j,iblk) = c0
               a11_4 (i,j,iblk) = c0
               a12_1 (i,j,iblk) = c0
               a12_2 (i,j,iblk) = c0
               a12_3 (i,j,iblk) = c0
               a12_4 (i,j,iblk) = c0
            enddo
            enddo
            do j = 1, ny_block
            do i = 1, nghost
               a11_1 (i,j,iblk) = c0
               a11_2 (i,j,iblk) = c0
               a11_3 (i,j,iblk) = c0
               a11_4 (i,j,iblk) = c0
               a12_1 (i,j,iblk) = c0
               a12_2 (i,j,iblk) = c0
               a12_3 (i,j,iblk) = c0
               a12_4 (i,j,iblk) = c0
            enddo
            enddo
         enddo
         !$OMP END PARALLEL DO

      end subroutine read_restart_eap

!=======================================================================

      end module ice_dyn_eap

!=======================================================================
