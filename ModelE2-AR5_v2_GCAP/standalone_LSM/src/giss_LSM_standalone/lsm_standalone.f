!@sum Routine to run the standalone version of the 
!@sum NASA GISS Land Surface Model (LSM) coupled to the
!@sum Ent Dynamic Global Terrestrial Ecosystem Model (DGTEM)
!@sum Options for meteorological forcings curently include:
!@sum     1) GSWP2 reananalysis data
!@sum     2) FLUXNET data for a single cell at the ecosystem scale
!%sum     3) GISS GCM output
!@auth I. Aleinov, M.J. Puma, Y. Kim
!-----------------------------------------------------------------------
#define USE_GSWP_FORCINGS
!#define RESOL_2x2_5
!#define ECOSYSTEM_SCALE  ! Single cell w/ FLUXNET met; Keep USE_GSWP_FORCINGS


#ifndef ECOSYSTEM_SCALE
#define PRINT_DIAGS ! Option to print offline global monthly diagnostics
#endif

!-----------------------------------------------------------------------
      module lsm
      use ent_mod
      use ghy_h, only : ngm, imt, nlsn, LS_NFRAC
      use sle001, only : advnc
      use filemanager, only : openunit, closeunit

      implicit none
!     Set boundary for global run    
#ifdef RESOL_2x2_5
      integer, parameter :: im=144, jm=90
#else
      integer, parameter :: im=72, jm=46
#endif
!     Set boundary for parallelization      
      integer,parameter :: j_0h=1, j_1h=jm
      integer,parameter :: i_0h=1, i_1h=im
!     Parameters for ecosystem-scale runs
      integer :: i_site, j_site
      real*8 :: latd
      character*120 :: LAI_file='fluxnet_LAI'
      character*120 :: vheight='fluxnet_vheight'
!     Logical flags for the Ent DGTEM
      logical :: force_VEG
      logical :: do_soilinit, do_soilresp
      logical :: do_phenology_activegrowth, do_structuralgrowth
      logical :: do_frost_hardiness, do_init_geo
      logical :: do_patchdynamics, do_spinup

      type t_lsm_state
        real*8, pointer :: w_ij(:,:,:,:)      ! soil water [m]
        real*8, pointer :: ht_ij(:,:,:,:)     ! soil heat [J]
        integer, pointer :: nsn_ij(:, :, :)   ! number of snow layers [-]     
        real*8, pointer :: dzsn_ij(:, :, :, :)! snow layer thicknesses [m]
        real*8, pointer :: wsn_ij(:, :, :, :) ! snow layer water-equiv.depth[m]
        real*8, pointer :: hsn_ij(:, :, :, :) ! snow layer heat contents [J]
        real*8, pointer :: fr_snow_ij(:, :, :)! fraction of land with snowcover
        real*8, pointer :: Qf_ij(:,:)! Foliage surf. vapor mixing ratio [kg/kg]
        type (entcelltype_public), pointer :: entcells(:,:)
      end type t_lsm_state

      type t_lsm_bc
        real*8, pointer ::
     &       fearth(:,:),
     &       top_index_ij(:,:),                 
     &       top_dev_ij(:,:),                   
     &       dz_ij(:,:,:),                   
     &       q_ij(:,:,:,:),              
     &       qk_ij(:,:,:,:),             
     &       sl_ij(:,:)        
      end type t_lsm_bc

      contains

!-----------------------------------------------------------------------

      subroutine read_input_parameters(itime1,jyear,jday,tyr_sec,dtsec
     &                     ,spin_start_yr, spin_end_yr, num_times_spin)
      use filemanager
      use parser
      use param

!@var Timestep size of forcings data [seconds]
      real*8  :: dtsec
!@var Simulation end time [number of timesteps of size dtsec]
      integer :: itime1
!@var Counter for time of year [seconds]
      integer :: tyr_sec
!@var Gregorian day of year (julian day of year for modelE met forcings)
      integer :: jday
!@var Gregorian year (julian year for modelE met forcings)
      integer :: jyear

      integer :: spin_start_yr, spin_end_yr, num_times_spin
      character*120 :: ifile="ent_input"
      integer iu_IFILE
!     Binary flags for the Ent DGTEM
      integer :: force_VEG_bi
      integer :: do_soilinit_bi,do_soilresp_bi
      integer :: do_phenology_activegrowth_bi,do_structuralgrowth_bi
      integer :: do_frost_hardiness_bi,do_init_geo_bi
      integer :: do_patchdynamics_bi,do_spinup_bi

      namelist /PARAMETERS/ i_site,j_site,latd,
     &     itime1,jyear,jday,tyr_sec,dtsec,
     &     spin_start_yr, spin_end_yr, num_times_spin,
     &     force_VEG_bi,
     &     do_soilinit_bi, do_soilresp_bi, 
     &     do_phenology_activegrowth_bi,do_structuralgrowth_bi,
     &     do_frost_hardiness_bi,do_init_geo_bi,
     &     do_patchdynamics_bi,do_spinup_bi

      ! read input file with run settings
      print *,"reading input parameters from ", ifile
      call openunit(trim(ifile),iu_IFILE,.false.,.true.)
      call parse_params(iu_IFILE)

      call sync_param( "itime1", itime1 )
      call sync_param( "jyear", jyear )
      call sync_param( "jday", jday )
      call sync_param( "tyr_sec", tyr_sec )
      call sync_param( "latd", latd )
      call sync_param( "i_site", i_site )
      call sync_param( "j_site", j_site )
      call sync_param( "dtsec", dtsec )
      call sync_param( "spin_start_yr", spin_start_yr )
      call sync_param( "spin_end_yr", spin_end_yr )
      call sync_param( "num_times_spin", num_times_spin )
      call sync_param( "force_VEG_bi", force_VEG_bi)
      call sync_param( "do_soilinit_bi",do_soilinit_bi)
      call sync_param( "do_soilresp_bi",do_soilresp_bi)
      call sync_param( "do_phenology_activegrowth_bi",
     &                 do_phenology_activegrowth_bi)
      call sync_param( "do_structuralgrowth_bi",do_structuralgrowth_bi)
      call sync_param( "do_frost_hardiness_bi",do_frost_hardiness_bi)
      call sync_param( "do_init_geo_bi",do_init_geo_bi)
      call sync_param( "do_patchdynamics_bi",do_patchdynamics_bi)
      call sync_param( "do_spinup_bi",do_spinup_bi)
      print *,"finished reading input parameters from ", ifile
      
!     Convert binary flags to logicals
      if (force_VEG_bi==1) then
         force_VEG=.true.
      else
         force_VEG=.false.
      endif

      if (do_soilinit_bi==1) then
         do_soilinit=.true.
      else
         do_soilinit=.false.
      endif

      if (do_soilresp_bi==1) then
         do_soilresp=.true.
      else
         do_soilresp=.false.
      endif

      if (do_phenology_activegrowth_bi==1) then
         do_phenology_activegrowth=.true.
      else
         do_phenology_activegrowth=.false.
      endif

      if (do_structuralgrowth_bi==1) then
         do_structuralgrowth=.true.
      else
         do_structuralgrowth=.false.
      endif

      if (do_frost_hardiness_bi==1) then
         do_frost_hardiness=.true.
      else
         do_frost_hardiness=.false.
      endif

      if (do_init_geo_bi==1) then
         do_init_geo=.true.
      else
         do_init_geo=.false.
      endif

      if (do_patchdynamics_bi==1) then
         do_patchdynamics=.true.
      else
         do_patchdynamics=.false.
      endif

      if (do_spinup_bi==1) then
         do_spinup=.true.
      else
         do_spinup=.false.
      endif

      call closeunit(iu_IFILE)

!YKIM
!possible combinations of options
!do_phenology_activegrowth=true & do_struturalgrowth=true
!do_phenology_activegrowth=true & do_struturalgrowth=false
!do_phenology_activegrowth=false & do_struturalgrowth=false
!impossible combinations of options
!do_phenology_activegrowth=false & do_struturalgrowth=true
      if (.not.do_phenology_activegrowth .and.
     &   do_structuralgrowth) then
         print*,"impossible combinations of input parameters"
         stop
      endif

      return

      end subroutine read_input_parameters

!-----------------------------------------------------------------------

      subroutine lsm_init_state( s, fearth, jday, jyear )
      type(t_lsm_state) :: s
      real*8 :: fearth(:,:)
      integer jday, jyear
      integer i,j
      integer sumj(jm)
      real*8, dimension(N_COVERTYPES,1:im,1:jm) :: laidata  !cohort
      real*8, dimension(N_COVERTYPES,1:im,1:jm) :: veg_height
#ifdef ECOSYSTEM_SCALE
      real*8::
     &     w_ij_temp(0:ngm,ls_nfrac,im,jm),
     &     ht_ij_temp(0:ngm,ls_nfrac,im,jm),
     &     nsn_ij_temp(     2,im,jm),
     &     dzsn_ij_temp(nlsn,2,im,jm),
     &     wsn_ij_temp(nlsn,2,im,jm),
     &     hsn_ij_temp(nlsn,2,im,jm),
     &     fr_snow_ij_temp(2,im,jm),
     &     Qf_ij_temp(im,jm) 

      allocate(
     &     s%w_ij(0:ngm,ls_nfrac,1:1,1:1),
     &     s%ht_ij(0:ngm,ls_nfrac,1:1,1:1),
     &     s%nsn_ij(     2,1:1,1:1),
     &     s%dzsn_ij(nlsn,2,1:1,1:1),
     &     s%wsn_ij(nlsn,2,1:1,1:1),
     &     s%hsn_ij(nlsn,2,1:1,1:1),
     &     s%fr_snow_ij(2,1:1,1:1),
     &     s%Qf_ij(1:1,1:1) )

#else
      allocate(
     &     s%w_ij(0:ngm,ls_nfrac,im,j_0h:j_1h),
     &     s%ht_ij(0:ngm,ls_nfrac,im,j_0h:j_1h),
     &     s%nsn_ij(     2,im,j_0h:j_1h),
     &     s%dzsn_ij(nlsn,2,im,j_0h:j_1h),
     &     s%wsn_ij(nlsn,2,im,j_0h:j_1h),
     &     s%hsn_ij(nlsn,2,im,j_0h:j_1h),
     &     s%fr_snow_ij(2,im,j_0h:j_1h),
     &     s%Qf_ij(im,j_0h:j_1h) )

#endif

#ifdef ECOSYSTEM_SCALE
      read(952)
     &       w_ij_temp(0:ngm,1:ls_nfrac,1:im,1:jm),
     &       ht_ij_temp(0:ngm,1:ls_nfrac,1:im,1:jm)

      ! something is wrong with snow - remove it for now
      nsn_ij_temp = 0
      dzsn_ij_temp = 0.d0
      wsn_ij_temp = 0.d0
      hsn_ij_temp = 0.d0
      fr_snow_ij_temp = 0.d0
      Qf_ij_temp = 0.1d0
      laidata = 0.d0
      veg_height = 0.d0

      s%w_ij(0:ngm,1:ls_nfrac,1,1)= w_ij_temp(0:ngm,1:ls_nfrac,
     &       i_site,j_site)
      s%ht_ij(0:ngm,1:ls_nfrac,1,1)=ht_ij_temp(0:ngm,1:ls_nfrac,
     &       i_site,j_site)
      s%nsn_ij(1:2, 1,1)= nsn_ij_temp(    1: 2,i_site,j_site)
      s%dzsn_ij(1:nlsn, 1:2, 1,1)=dzsn_ij_temp(1:nlsn,1:2,i_site,j_site)
      s%wsn_ij(1:nlsn, 1:2, 1,1)=wsn_ij_temp(1:nlsn,1:2,i_site,j_site)
      s%hsn_ij(1:nlsn, 1:2, 1,1)=hsn_ij_temp(1:nlsn,1:2,i_site,j_site)
      s%fr_snow_ij(1:2, 1,1)=fr_snow_ij_temp(1:2,i_site,j_site)
      s%Qf_ij(1,1)= Qf_ij_temp(i_site,j_site) 


      allocate( s%entcells(1,1) )
      call ent_cell_nullify( s%entcells )
      call ent_cell_construct( s%entcells(1,1) )


      call set_veg_data_single( s%entcells,
     &     im, jm, i_site, j_site, jday, jyear, laidata)


#else
      read(952)
     &       s%w_ij(0:ngm,1:ls_nfrac,:,:),         
     &       s%ht_ij(0:ngm,1:ls_nfrac,:,:),        
     &       s%nsn_ij    (1:2, :,:),              
     &       s%dzsn_ij   (1:nlsn, 1:2, :,:),      
     &       s%wsn_ij    (1:nlsn, 1:2, :,:),      
     &       s%hsn_ij    (1:nlsn, 1:2, :,:),      
     &       s%fr_snow_ij(1:2, :,:)
!     &      ,s%Qf_ij(:,:)

      ! something is wrong with snow - remove it for now
      s%nsn_ij = 0
      s%dzsn_ij = 0.d0
      s%wsn_ij = 0.d0
      s%hsn_ij = 0.d0
      s%fr_snow_ij = 0.d0
      s%Qf_ij(:,:)=0.1d0

      allocate( s%entcells(im,j_0h:j_1h) )
      call ent_cell_nullify( s%entcells )

      sumj(:) = 0
      do j=1,jm
        do i=1,im
          if ( fearth(i,j) > 0.d0 )
     &         call ent_cell_construct( s%entcells(i,j) )
          if ( fearth(i,j) > 0.d0 ) sumj(j) = sumj(j) + 1
        enddo
      enddo

      print *,"sumj= ", sumj

      call set_vegetation_data( s%entcells,
     &     IM, JM, 1, IM, 1, JM, jday, jyear )
#endif

      end subroutine lsm_init_state

!-----------------------------------------------------------------------

      subroutine lsm_init_bc( bc )
      type (t_lsm_bc) :: bc
      integer :: iu_TOPO, iu_TOP_INDEX, iu_SOIL
      character*80 :: title
      real*4, allocatable :: temp_local(:,:,:)
      real*4, allocatable :: temp4(:,:)

#ifdef ECOSYSTEM_SCALE
        real*8 :: dz_ij_temp(im,jm,1:ngm),
     &       q_ij_temp(im,jm,1:imt,1:ngm),
     &       qk_ij_temp(im,jm,1:imt,1:ngm),
     &       fearth_temp(im,jm),
     &       top_index_ij_temp(im,jm), 
     &       top_dev_ij_temp(im,jm), 
     &       sl_ij_temp(im,jm)

        allocate(
     &     bc%fearth      (1:1,1:1),
     &     bc%top_index_ij(1:1,1:1),
     &     bc%top_dev_ij  (1:1,1:1),
     &     bc%dz_ij       (1:1,1:1, ngm),
     &     bc%q_ij        (1:1,1:1, imt, ngm),
     &     bc%qk_ij       (1:1,1:1, imt, ngm),
     &     bc%sl_ij       (1:1,1:1) )
     
        read(953)
     &       fearth_temp(:,:),
     &       top_index_ij_temp(:,:),                 
     &       top_dev_ij_temp(:,:),                   
     &       dz_ij_temp(:,:,1:ngm),                   
     &       q_ij_temp(:,:,1:imt,1:ngm),              
     &       qk_ij_temp(:,:,1:imt,1:ngm),             
     &       sl_ij_temp(:,:) 
     
      
        bc%fearth(1,1)           = fearth_temp(i_site,j_site)
        bc%top_index_ij(1,1)     = top_index_ij_temp(i_site,j_site)
        bc%top_dev_ij(1,1)       = top_dev_ij_temp(i_site,j_site)
        bc%dz_ij(1,1,1:ngm)      = dz_ij_temp(i_site,j_site,1:ngm)
        bc%q_ij(1,1,1:imt,1:ngm) = q_ij_temp(i_site,j_site,1:imt,1:ngm)
        bc%qk_ij(1,1,1:imt,1:ngm)= qk_ij_temp(i_site,j_site,1:imt,1:ngm)
        bc%sl_ij(1,1)            = sl_ij_temp(i_site,j_site)

#else
        allocate(
     &     bc%fearth      (im, j_0h:j_1h),
     &     bc%top_index_ij(im, j_0h:j_1h),
     &     bc%top_dev_ij  (im, j_0h:j_1h),
     &     bc%dz_ij       (im, j_0h:j_1h, ngm),
     &     bc%q_ij        (im, j_0h:j_1h, imt, ngm),
     &     bc%qk_ij       (im, j_0h:j_1h, imt, ngm),
     &     bc%sl_ij       (im, j_0h:j_1h) )

cddd        read(953)
cddd     &       bc%fearth(:,:),
cddd     &       bc%top_index_ij(:,:),                 
cddd     &       bc%top_dev_ij(:,:),                   
cddd     &       bc%dz_ij(:,:,1:ngm),                   
cddd     &       bc%q_ij(:,:,1:imt,1:ngm),              
cddd     &       bc%qk_ij(:,:,1:imt,1:ngm),             
cddd     &       bc%sl_ij(:,:) 

        allocate( temp4(im, jm) )
        call openunit("TOPO",iu_TOPO,.true.,.true.)
        read(iu_TOPO)
        read(iu_TOPO)
        read(iu_TOPO) title, temp4
        bc%fearth(:,:) = temp4(:,:)
        call closeunit(iu_TOPO)
        call openunit("TOP_INDEX",iu_TOP_INDEX,.true.,.true.)
        read(iu_TOP_INDEX) title, temp4
        bc%top_index_ij(:,:) = temp4(:,:)
        read(iu_TOP_INDEX) title, temp4
        bc%top_dev_ij(:,:) = temp4(:,:)
        call closeunit(iu_TOP_INDEX)
        deallocate( temp4 )

        call openunit("SOIL",iu_SOIL,.true.,.true.)
        allocate(temp_local(1:im,1:jm,11*ngm+1))
        read(iu_SOIL) temp_local
        bc%dz_ij(:,:,:)   = temp_local(:,:,1:ngm)
        bc%q_ij(:,:,:,:) =
     &       reshape( temp_local(:,:,1+ngm:) ,
     &       (/im,jm,imt,ngm/) )
        bc%qk_ij(:,:,:,:) =
     &       reshape( temp_local(:,:,1+ngm+ngm*imt:) ,
     &       (/im,jm,imt,ngm/) )
        bc%sl_ij(:,:)  =
     &       temp_local(:,:,1+ngm+ngm*imt+ngm*imt)
        deallocate(temp_local)
        call closeunit (iu_SOIL)

#endif

      end subroutine lsm_init_bc

!-----------------------------------------------------------------------

      subroutine get_forcings(
     &         force_Ca              ,
     &         force_cos_zen_angle   ,
     &         force_vis_rad         ,
     &         force_direct_vis_rad  ,
     &         force_prec_ms         ,
     &         force_eprec_w         ,
     &         force_sprec_ms        ,
     &         force_seprec_w        ,
     &         force_srheat          ,
     &         force_trheat          ,
     &         force_ts              ,
     &         force_qs              ,
     &         force_ps              ,
     &         force_rhosrf          ,
     &         force_cdh             ,
     &         force_qm1             ,
     &         force_ws              ,
     &         force_pbl_args_ws0    ,
     &         force_pbl_args_tprime ,
     &         force_pbl_args_qprime )
      real*8
     &         force_Ca (1:im,1:jm)             ,
     &         force_cos_zen_angle (1:im,1:jm)  ,
     &         force_vis_rad (1:im,1:jm)        ,
     &         force_direct_vis_rad (1:im,1:jm) ,
     &         force_prec_ms (1:im,1:jm)        ,
     &         force_eprec_w (1:im,1:jm)        ,
     &         force_sprec_ms (1:im,1:jm)       ,
     &         force_seprec_w (1:im,1:jm)       ,
     &         force_srheat (1:im,1:jm)         ,
     &         force_trheat (1:im,1:jm)         ,
     &         force_ts (1:im,1:jm)             ,
     &         force_qs (1:im,1:jm)             ,
     &         force_ps (1:im,1:jm)             ,
     &         force_rhosrf (1:im,1:jm)         ,
     &         force_cdh (1:im,1:jm)            ,
     &         force_qm1 (1:im,1:jm)            ,
     &         force_ws (1:im,1:jm)             ,
     &         force_pbl_args_ws0 (1:im,1:jm)   ,
     &         force_pbl_args_tprime (1:im,1:jm),
     &         force_pbl_args_qprime (1:im,1:jm)

      print *,"reading forcings"

      read(951)
     &         force_Ca (1:im,1:jm)             ,
     &         force_cos_zen_angle (1:im,1:jm)  ,
     &         force_vis_rad (1:im,1:jm)        ,
     &         force_direct_vis_rad (1:im,1:jm) ,
     &         force_prec_ms (1:im,1:jm)        ,
     &         force_eprec_w (1:im,1:jm)        ,
     &         force_sprec_ms (1:im,1:jm)       ,
     &         force_seprec_w (1:im,1:jm)       ,
     &         force_srheat (1:im,1:jm)         ,
     &         force_trheat (1:im,1:jm)         ,
     &         force_ts (1:im,1:jm)             ,
     &         force_qs (1:im,1:jm)             ,
     &         force_ps (1:im,1:jm)             ,
     &         force_rhosrf (1:im,1:jm)         ,
     &         force_cdh (1:im,1:jm)            ,
     &         force_qm1 (1:im,1:jm)            ,
     &         force_ws (1:im,1:jm)             ,
     &         force_pbl_args_ws0 (1:im,1:jm)   ,
     &         force_pbl_args_tprime (1:im,1:jm),
     &         force_pbl_args_qprime (1:im,1:jm)

      print *,"forcings ok"

      end subroutine get_forcings

!-----------------------------------------------------------------------

      subroutine lsm_initialize( dt_in , n_month , tbcs_ij, tcanopy_ij
#ifdef PRINT_DIAGS
     &                  , aevap_mnth,aevapw_mnth,aevapd_mnth,aevapb_mnth
     &                  , aintercep_mnth,aruns_mnth,arunu_mnth
     &                  , agpp_mnth,arauto_mnth,asoilresp_mnth
     &                  , aevapvg_mnth,aevapvs_mnth,aevapbs_mnth
     &                  , asrht_mnth,atrht_mnth,aalbedo_mnth
     &                  , aclab_mnth,asoilCpoolsum_mnth,aepp_mnth
     &                  , atrg_mnth,ashg_mnth,alhg_mnth,aepc_mnth
     &                  , n_count_mnth
     &                 , w_b1_mnth,w_b2_mnth,w_b3_mnth
     &                 , w_b4_mnth,w_b5_mnth,w_b6_mnth
     &                 , w_v1_mnth,w_v2_mnth,w_v3_mnth
     &                 , w_v4_mnth,w_v5_mnth,w_v6_mnth
     &                 , ht_b1_mnth,ht_b2_mnth,ht_b3_mnth
     &                 , ht_b4_mnth,ht_b5_mnth,ht_b6_mnth
     &                 , ht_v1_mnth,ht_v2_mnth,ht_v3_mnth
     &                 , ht_v4_mnth,ht_v5_mnth,ht_v6_mnth
     &                 , dzsn_b1,dzsn_b2,dzsn_b3,dzsn_v1,dzsn_v2,dzsn_v3
     &                 , wsn_b1 ,wsn_b2 ,wsn_b3 ,wsn_v1 ,wsn_v2 ,wsn_v3
     &                 , hsn_b1 ,hsn_b2 ,hsn_b3 ,hsn_v1 ,hsn_v2 ,hsn_v3
     &                 , nsn_b,nsn_v,frsnow_b,frsnow_v,Qf_mnth
     &                 , abetad_mnth
     &  , aClivepool_leaf_m,aClivepool_froot_m,aClivepool_wood_m
     &  , aCdeadpool_surfmet_m,aCdeadpool_surfstr_m,aCdeadpool_soilmet_m
     &  , aCdeadpool_soilstr_m,aCdeadpool_cwd_m,aCdeadpool_surfmic_m
     &  , aCdeadpool_soilmic_m,aCdeadpool_slow_m,aCdeadpool_passive_m
     &  , alai_m, canopyH2O_m, canopyheat_m
#endif
     &                           )

      use sle001, only : hl0, dt

#ifdef ECOSYSTEM_SCALE
      real*8, intent(out) :: tbcs_ij, tcanopy_ij
#else
      real*8,dimension(im,jm), intent(out) :: tbcs_ij, tcanopy_ij
#endif
      real*8 :: dt_in
      integer :: n_month
#ifdef PRINT_DIAGS
      real*8, dimension(im,jm,12),intent(out) ::
     &                    aevap_mnth,aevapw_mnth,aevapd_mnth,aevapb_mnth
     &                  , aintercep_mnth,aruns_mnth,arunu_mnth
     &                  , agpp_mnth,arauto_mnth,asoilresp_mnth
     &                  , aevapvg_mnth,aevapvs_mnth,aevapbs_mnth
     &                  , asrht_mnth,atrht_mnth,aalbedo_mnth
     &                  , aclab_mnth,asoilCpoolsum_mnth,aepp_mnth
     &                  , atrg_mnth,ashg_mnth,alhg_mnth,aepc_mnth
     &                  , w_b1_mnth,w_b2_mnth,w_b3_mnth
     &                  , w_b4_mnth,w_b5_mnth,w_b6_mnth
     &                  , w_v1_mnth,w_v2_mnth,w_v3_mnth
     &                  , w_v4_mnth,w_v5_mnth,w_v6_mnth
     &                  , ht_b1_mnth,ht_b2_mnth,ht_b3_mnth
     &                  , ht_b4_mnth,ht_b5_mnth,ht_b6_mnth
     &                  , ht_v1_mnth,ht_v2_mnth,ht_v3_mnth
     &                  , ht_v4_mnth,ht_v5_mnth,ht_v6_mnth
     &                 , dzsn_b1,dzsn_b2,dzsn_b3,dzsn_v1,dzsn_v2,dzsn_v3
     &                 , wsn_b1 ,wsn_b2 ,wsn_b3 ,wsn_v1 ,wsn_v2 ,wsn_v3
     &                 , hsn_b1 ,hsn_b2 ,hsn_b3 ,hsn_v1 ,hsn_v2 ,hsn_v3
     &                 , nsn_b,nsn_v,frsnow_b,frsnow_v,Qf_mnth
     &                 , abetad_mnth
     &  , aClivepool_leaf_m,aClivepool_froot_m,aClivepool_wood_m
     &  , aCdeadpool_surfmet_m,aCdeadpool_surfstr_m,aCdeadpool_soilmet_m
     &  , aCdeadpool_soilstr_m,aCdeadpool_cwd_m,aCdeadpool_surfmic_m
     &  , aCdeadpool_soilmic_m,aCdeadpool_slow_m,aCdeadpool_passive_m
     &  , alai_m, canopyH2O_m, canopyheat_m
      real*8,dimension(12),intent(out) :: n_count_mnth
#endif

      call ent_init_config(
     &     do_soilresp=do_soilresp
     &     ,do_phenology_activegrowth=do_phenology_activegrowth
     &     ,do_structuralgrowth=do_structuralgrowth
     &     ,do_frost_hardiness=do_frost_hardiness
     &     ,do_init_geo=do_init_geo
     &     ,do_patchdynamics=do_patchdynamics)

      n_month = 0          ! intialize month number
      dt = dt_in           ! intialize timestep size [sec]
#ifdef ECOSYSTEM_SCALE
      tbcs_ij = -1d30 ! initialize ground surface temperature[C]
      tcanopy_ij = -1d30 !initialize canopy temperature[C]
#else
      tbcs_ij(:,:) = -1d30 ! initialize ground surface temperature[C]
      tcanopy_ij(:,:) = -1d30 !initialize canopy temperature[C]
#endif
      call hl0

#ifdef PRINT_DIAGS
      call zero_diags(   aevap_mnth,aevapw_mnth,aevapd_mnth,aevapb_mnth
     &                 , aintercep_mnth,aruns_mnth,arunu_mnth
     &                 , agpp_mnth,arauto_mnth,asoilresp_mnth
     &                 , aevapvg_mnth,aevapvs_mnth,aevapbs_mnth
     &                 , asrht_mnth,atrht_mnth,aalbedo_mnth
     &                 , aclab_mnth,asoilCpoolsum_mnth,aepp_mnth
     &                 , atrg_mnth,ashg_mnth,alhg_mnth,aepc_mnth
     &                 , n_count_mnth
     &                 , w_b1_mnth,w_b2_mnth,w_b3_mnth
     &                 , w_b4_mnth,w_b5_mnth,w_b6_mnth
     &                 , w_v1_mnth,w_v2_mnth,w_v3_mnth
     &                 , w_v4_mnth,w_v5_mnth,w_v6_mnth
     &                 , ht_b1_mnth,ht_b2_mnth,ht_b3_mnth
     &                 , ht_b4_mnth,ht_b5_mnth,ht_b6_mnth
     &                 , ht_v1_mnth,ht_v2_mnth,ht_v3_mnth
     &                 , ht_v4_mnth,ht_v5_mnth,ht_v6_mnth
     &                 , dzsn_b1,dzsn_b2,dzsn_b3,dzsn_v1,dzsn_v2,dzsn_v3
     &                 , wsn_b1 ,wsn_b2 ,wsn_b3 ,wsn_v1 ,wsn_v2 ,wsn_v3
     &                 , hsn_b1 ,hsn_b2 ,hsn_b3 ,hsn_v1 ,hsn_v2 ,hsn_v3
     &                 , nsn_b,nsn_v,frsnow_b,frsnow_v,Qf_mnth
     &                 , abetad_mnth
     &  , aClivepool_leaf_m,aClivepool_froot_m,aClivepool_wood_m
     &  , aCdeadpool_surfmet_m,aCdeadpool_surfstr_m,aCdeadpool_soilmet_m
     &  , aCdeadpool_soilstr_m,aCdeadpool_cwd_m,aCdeadpool_surfmic_m
     &  , aCdeadpool_soilmic_m,aCdeadpool_slow_m,aCdeadpool_passive_m
     &  , alai_m,canopyH2O_m,canopyheat_m)
#endif

      end subroutine lsm_initialize

!-----------------------------------------------------------------------

      subroutine lsm_run(s,bc,jday,jyear,tyr_sec,dtsec,j0,j1,k_mnth
     &                 , iu_LAI, iu_vht
#ifdef PRINT_DIAGS
     &                 , aevap_mnth,aevapw_mnth,aevapd_mnth,aevapb_mnth
     &                 , aintercep_mnth,aruns_mnth,arunu_mnth
     &                 , agpp_mnth,arauto_mnth,asoilresp_mnth
     &                 , aevapvg_mnth,aevapvs_mnth,aevapbs_mnth
     &                 , asrht_mnth,atrht_mnth,aalbedo_mnth
     &                 , aclab_mnth,asoilCpoolsum_mnth,aepp_mnth
     &                 , atrg_mnth,ashg_mnth,alhg_mnth,aepc_mnth
     &                 , n_count_mnth
     &                 , w_b1_mnth,w_b2_mnth,w_b3_mnth
     &                 , w_b4_mnth,w_b5_mnth,w_b6_mnth
     &                 , w_v1_mnth,w_v2_mnth,w_v3_mnth
     &                 , w_v4_mnth,w_v5_mnth,w_v6_mnth
     &                 , ht_b1_mnth,ht_b2_mnth,ht_b3_mnth
     &                 , ht_b4_mnth,ht_b5_mnth,ht_b6_mnth
     &                 , ht_v1_mnth,ht_v2_mnth,ht_v3_mnth
     &                 , ht_v4_mnth,ht_v5_mnth,ht_v6_mnth
     &                 , dzsn_b1,dzsn_b2,dzsn_b3,dzsn_v1,dzsn_v2,dzsn_v3
     &                 , wsn_b1 ,wsn_b2 ,wsn_b3 ,wsn_v1 ,wsn_v2 ,wsn_v3
     &                 , hsn_b1 ,hsn_b2 ,hsn_b3 ,hsn_v1 ,hsn_v2 ,hsn_v3
     &                 , nsn_b,nsn_v,frsnow_b,frsnow_v,Qf_mnth
     &                 , abetad_mnth
     &  , aClivepool_leaf_m,aClivepool_froot_m,aClivepool_wood_m
     &  , aCdeadpool_surfmet_m,aCdeadpool_surfstr_m,aCdeadpool_soilmet_m
     &  , aCdeadpool_soilstr_m,aCdeadpool_cwd_m,aCdeadpool_surfmic_m
     &  , aCdeadpool_soilmic_m,aCdeadpool_slow_m,aCdeadpool_passive_m
     &  , alai_m,canopyH2O_m,canopyheat_m
#endif
#ifdef USE_GSWP_FORCINGS
     &   , iu_vector,tbcs_ij,tcanopy_ij
     &   , end_of_input_flag
#endif
     &   , end_of_day_flag
     &     )

#ifdef USE_GSWP_FORCINGS
      use drv_gswp_force, only : get_gswp_forcings
#endif
      use sle001, only : tp,tbcs,aevap,aevapw,aevapd,aevapb,aintercep
     &     ,aruns,arunu,agpp,arauto,asoilresp
     &     ,aevapvg,aevapvs,aevapbs,asrht,atrht,aalbedo
     &     ,aclab,asoilCpoolsum,aepp,atrg,ashg,alhg,aepc
     &     ,abetad,alai
     &     ,aClivepool_leaf,aClivepool_froot,aClivepool_wood
     &     ,aCdeadpool_surfmet,aCdeadpool_surfstr,aCdeadpool_soilmet
     &     ,aCdeadpool_soilstr,aCdeadpool_cwd,aCdeadpool_surfmic
     &     ,aCdeadpool_soilmic,aCdeadpool_slow,aCdeadpool_passive
     &     ,alai, acnc, acna

      use domain_decomp, only : mype
      type(t_lsm_state) :: s
      type (t_lsm_bc) :: bc

      integer, intent(in) :: jday, jyear, tyr_sec, j0, j1
      integer, intent(in) :: k_mnth !current month number
      real*8, intent(in) :: dtsec
      integer :: iu_LAI, iu_vht
      real*8 fb,fv
      integer i,j
      integer, save :: jday_old = -32768

#ifdef ECOSYSTEM_SCALE 
      real*8 ::
     &     force_Ca (1:1,1:1)             ,
     &     force_cos_zen_angle (1:1,1:1)  ,
     &     force_vis_rad (1:1,1:1)        ,
     &     force_direct_vis_rad (1:1,1:1) ,
     &     force_prec_ms (1:1,1:1)        ,
     &     force_eprec_w (1:1,1:1)        ,
     &     force_sprec_ms (1:1,1:1)       ,
     &     force_seprec_w (1:1,1:1)       ,
     &     force_srheat (1:1,1:1)         ,
     &     force_trheat (1:1,1:1)         ,
     &     force_ts (1:1,1:1)             ,
     &     force_qs (1:1,1:1)             ,
     &     force_ps (1:1,1:1)             ,
     &     force_rhosrf (1:1,1:1)         ,
     &     force_cdh (1:1,1:1)            ,
     &     force_qm1 (1:1,1:1)            ,
     &     force_ws (1:1,1:1)             ,
     &     force_pbl_args_ws0 (1:1,1:1)   ,
     &     force_pbl_args_tprime (1:1,1:1),
     &     force_pbl_args_qprime (1:1,1:1)
#else
      real*8 ::
     &     force_Ca (1:im,1:jm)             ,
     &     force_cos_zen_angle (1:im,1:jm)  ,
     &     force_vis_rad (1:im,1:jm)        ,
     &     force_direct_vis_rad (1:im,1:jm) ,
     &     force_prec_ms (1:im,1:jm)        ,
     &     force_eprec_w (1:im,1:jm)        ,
     &     force_sprec_ms (1:im,1:jm)       ,
     &     force_seprec_w (1:im,1:jm)       ,
     &     force_srheat (1:im,1:jm)         ,
     &     force_trheat (1:im,1:jm)         ,
     &     force_ts (1:im,1:jm)             ,
     &     force_qs (1:im,1:jm)             ,
     &     force_ps (1:im,1:jm)             ,
     &     force_rhosrf (1:im,1:jm)         ,
     &     force_cdh (1:im,1:jm)            ,
     &     force_qm1 (1:im,1:jm)            ,
     &     force_ws (1:im,1:jm)             ,
     &     force_pbl_args_ws0 (1:im,1:jm)   ,
     &     force_pbl_args_tprime (1:im,1:jm),
     &     force_pbl_args_qprime (1:im,1:jm)
#endif

#ifdef PRINT_DIAGS
      real*8, dimension(im,jm,12),intent(out) ::
     &                    aevap_mnth,aevapw_mnth,aevapd_mnth,aevapb_mnth
     &                  , aintercep_mnth,aruns_mnth,arunu_mnth
     &                  , agpp_mnth,arauto_mnth,asoilresp_mnth
     &                  , aevapvg_mnth,aevapvs_mnth,aevapbs_mnth
     &                  , asrht_mnth,atrht_mnth,aalbedo_mnth
     &                  , aclab_mnth,asoilCpoolsum_mnth,aepp_mnth
     &                  , atrg_mnth,ashg_mnth,alhg_mnth,aepc_mnth
     &                  ,w_b1_mnth,w_b2_mnth,w_b3_mnth
     &                  ,w_b4_mnth,w_b5_mnth,w_b6_mnth
     &                  ,w_v1_mnth,w_v2_mnth,w_v3_mnth
     &                  ,w_v4_mnth,w_v5_mnth,w_v6_mnth
     &                  ,ht_b1_mnth,ht_b2_mnth,ht_b3_mnth
     &                  ,ht_b4_mnth,ht_b5_mnth,ht_b6_mnth
     &                  ,ht_v1_mnth,ht_v2_mnth,ht_v3_mnth
     &                  ,ht_v4_mnth,ht_v5_mnth,ht_v6_mnth
     &                  ,dzsn_b1,dzsn_b2,dzsn_b3,dzsn_v1,dzsn_v2,dzsn_v3
     &                  ,wsn_b1 ,wsn_b2 ,wsn_b3 ,wsn_v1 ,wsn_v2 ,wsn_v3
     &                  ,hsn_b1 ,hsn_b2 ,hsn_b3 ,hsn_v1 ,hsn_v2 ,hsn_v3
     &                  ,nsn_b,nsn_v,frsnow_b,frsnow_v,Qf_mnth
     &                  ,abetad_mnth
     &  , aClivepool_leaf_m,aClivepool_froot_m,aClivepool_wood_m
     &  , aCdeadpool_surfmet_m,aCdeadpool_surfstr_m,aCdeadpool_soilmet_m
     &  , aCdeadpool_soilstr_m,aCdeadpool_cwd_m,aCdeadpool_surfmic_m
     &  , aCdeadpool_soilmic_m,aCdeadpool_slow_m,aCdeadpool_passive_m
     &  , alai_m,canopyH2O_m,canopyheat_m

      real*8, dimension(12) :: n_count_mnth
#endif

#ifdef USE_GSWP_FORCINGS
      integer, parameter :: varnums = 9      ! # of GSWP2 variables
      logical :: end_of_input_flag
#ifdef ECOSYSTEM_SCALE
      integer :: iu_vector
      real*8, intent(inout) :: tbcs_ij,tcanopy_ij
      real*8, dimension(N_COVERTYPES,1,1) :: laidata  !cohort
      real*8, dimension(N_COVERTYPES,1,1) :: veg_height
#else
      integer,dimension(varnums) :: iu_vector
      real*8,dimension(im,jm), intent(inout) :: tbcs_ij,tcanopy_ij
#endif


#else
      integer,dimension(varnums) :: iu_vector
      real*8,dimension(im,jm), intent(inout) :: tbcs_ij,tcanopy_ij
#endif

      logical end_of_day_flag

      !call sysusage(mype+4,1)

#ifdef USE_GSWP_FORCINGS
      call get_gswp_forcings(
     &     i_site, j_site, latd,
     &     jday,
     &     jyear,
     &     tyr_sec,
     &     dtsec,
     &     iu_vector,
     &     tbcs_ij,
     &     tcanopy_ij,
     &     end_of_input_flag,
#else
          call get_forcings(
#endif
#ifdef ECOSYSTEM_SCALE
     &         force_Ca (1:1,1:1)             ,
     &         force_cos_zen_angle (1:1,1:1)  ,
     &         force_vis_rad (1:1,1:1)        ,
     &         force_direct_vis_rad (1:1,1:1) ,
     &         force_prec_ms (1:1,1:1)        ,
     &         force_eprec_w (1:1,1:1)        ,
     &         force_sprec_ms (1:1,1:1)       ,
     &         force_seprec_w (1:1,1:1)       ,
     &         force_srheat (1:1,1:1)         ,
     &         force_trheat (1:1,1:1)         ,
     &         force_ts (1:1,1:1)             ,
     &         force_qs (1:1,1:1)             ,
     &         force_ps (1:1,1:1)             ,
     &         force_rhosrf (1:1,1:1)         ,
     &         force_cdh (1:1,1:1)            ,
     &         force_qm1 (1:1,1:1)            ,
     &         force_ws (1:1,1:1)             ,
     &         force_pbl_args_ws0 (1:1,1:1)   ,
     &         force_pbl_args_tprime (1:1,1:1),
     &         force_pbl_args_qprime (1:1,1:1)! xxx
     &         )
#else
     &         force_Ca (1:im,1:jm)             ,
     &         force_cos_zen_angle (1:im,1:jm)  ,
     &         force_vis_rad (1:im,1:jm)        ,
     &         force_direct_vis_rad (1:im,1:jm) ,
     &         force_prec_ms (1:im,1:jm)        ,
     &         force_eprec_w (1:im,1:jm)        ,
     &         force_sprec_ms (1:im,1:jm)       ,
     &         force_seprec_w (1:im,1:jm)       ,
     &         force_srheat (1:im,1:jm)         ,
     &         force_trheat (1:im,1:jm)         ,
     &         force_ts (1:im,1:jm)             ,
     &         force_qs (1:im,1:jm)             ,
     &         force_ps (1:im,1:jm)             ,
     &         force_rhosrf (1:im,1:jm)         ,
     &         force_cdh (1:im,1:jm)            ,
     &         force_qm1 (1:im,1:jm)            ,
     &         force_ws (1:im,1:jm)             ,
     &         force_pbl_args_ws0 (1:im,1:jm)   ,
     &         force_pbl_args_tprime (1:im,1:jm),
     &         force_pbl_args_qprime (1:im,1:jm)! xxx
     &         )
#endif
      !call sysusage(mype+4,2)
      !call sysusage(mype+12,1)

#ifdef ECOSYSTEM_SCALE
      if (force_VEG) 
     & call update_FLUXNET_LAI(iu_LAI,iu_vht,s%entcells(:,:)
     &       , im, jm, i_site, j_site, jday, jyear, laidata, veg_height)
#endif
      ! update vegetation only once per day
      if ( jday /= jday_old) then
#ifdef ECOSYSTEM_SCALE
        call update_veg_data_single( s%entcells(:,:)
     &       , im, jm, i_site, j_site, jday, jyear, laidata, veg_height)

#else ! GLOBAL_SCALE
        call update_vegetation_data( s%entcells(:,j0:j1),
     &       im, jm, 1, im, j0, j1, jday, jyear )

#endif
        jday_old = jday
      endif

      !call sysusage(mype+12,2)
      !call sysusage(mype+8,1)

#ifdef ECOSYSTEM_SCALE
      ! really fb, fv are not needed for Ent, but just in case...
      call ent_get_exports( s%entcells(1,1),
     &         fraction_of_vegetated_soil=fv
     &         )
      fb = 1.d0 - fv

      call advnc(
!-------------- Ent specific
     &         s%entcells(1,1), force_Ca(1,1),
     &         force_cos_zen_angle(1,1), force_vis_rad(1,1),
     &         force_direct_vis_rad(1,1),
     &         s%Qf_ij(1,1),
!-------------- old vegetation scheme (not implemented at the moment)
!     &         vegcell,
!-------------- prognostic vars
     &         s%w_ij(0:ngm,1:LS_NFRAC,1,1),         
     &         s%ht_ij(0:ngm,1:LS_NFRAC,1,1),        
     &         s%nsn_ij    (1:2, 1, 1),              
     &         s%dzsn_ij   (1:nlsn, 1:2, 1, 1),      
     &         s%wsn_ij    (1:nlsn, 1:2, 1, 1),      
     &         s%hsn_ij    (1:nlsn, 1:2, 1, 1),      
     &         s%fr_snow_ij(1:2, 1, 1),          
!-------------- BC's    
     &         bc%top_index_ij(1, 1),                 
     &         bc%top_dev_ij(1, 1),                   
     &         bc%dz_ij(1,1,1:ngm),                   
     &         bc%q_ij(1,1,1:imt,1:ngm),              
     &         bc%qk_ij(1,1,1:imt,1:ngm),             
     &         bc%sl_ij(1,1),                         
     &         fb,                                 
     &         fv,                        
!-------------- forcings         
     &         force_prec_ms (1,1)        ,
     &         force_eprec_w (1,1)        ,
     &         force_sprec_ms (1,1)       ,
     &         force_seprec_w (1,1)       ,
     &         0.d0,
     &         0.d0,
     &         force_srheat (1,1)         ,
     &         force_trheat (1,1)         ,
     &         force_ts (1,1)             ,
     &         force_qs (1,1)             ,
     &         force_ps (1,1)             ,
     &         force_rhosrf (1,1)         ,
     &         force_cdh (1,1)            ,
     &         force_qm1 (1,1)            ,
     &         force_ws (1,1)             ,
     &         force_pbl_args_ws0 (1,1)   ,
     &         force_pbl_args_tprime (1,1),
     &         force_pbl_args_qprime (1,1),
     &         end_of_day_flag)

      tbcs_ij=tbcs
      tcanopy_ij=tp(0,2)

      write(9995,'(150(1pe16.8))')
     &    force_prec_ms (1,1)        ,
     &    force_eprec_w (1,1)        ,
     &    force_sprec_ms (1,1)       ,
     &    force_seprec_w (1,1)       ,
     &    force_srheat (1,1)         ,
     &    force_trheat (1,1)         ,
     &    force_ts (1,1)             ,
     &    force_qs (1,1)             ,
     &    force_ps (1,1)             ,
     &    force_rhosrf (1,1)         ,
     &    force_cdh (1,1)            ,
     &    force_qm1 (1,1)            ,
     &    force_ws (1,1)             ,!13
     &    aevap,aevapw,aevapd,aevapb,aintercep,aruns,arunu, !20
     &    agpp,arauto,asoilresp,abetad, !24
     &    aevapvg,aevapvs,aevapbs,asrht,atrht,aalbedo, !30
     &    aclab,asoilCpoolsum,aepp,atrg,ashg,alhg,aepc,!37
     &    s%w_ij(1,1,1,1),s%w_ij(2,1,1,1),s%w_ij(3,1,1,1),!40
     &    s%w_ij(4,1,1,1),s%w_ij(5,1,1,1),s%w_ij(6,1,1,1),!43
     &    s%w_ij(1,2,1,1),s%w_ij(2,2,1,1),s%w_ij(3,2,1,1),!46
     &    s%w_ij(4,2,1,1),s%w_ij(5,2,1,1),s%w_ij(6,2,1,1),!49
     &    s%ht_ij(1,1,1,1),s%ht_ij(2,1,1,1),s%ht_ij(3,1,1,1),!52
     &    s%ht_ij(4,1,1,1),s%ht_ij(5,1,1,1),s%ht_ij(6,1,1,1),!55
     &    s%ht_ij(1,2,1,1),s%ht_ij(2,2,1,1),s%ht_ij(3,2,1,1),!58
     &    s%ht_ij(4,2,1,1),s%ht_ij(5,2,1,1),s%ht_ij(6,2,1,1),!61
     &    s%dzsn_ij(1,1,1,1),s%dzsn_ij(2,1,1,1),s%dzsn_ij(3,1,1,1),!64
     &    s%dzsn_ij(1,2,1,1),s%dzsn_ij(2,2,1,1),s%dzsn_ij(3,2,1,1),!67
     &    s%wsn_ij(1,1,1,1),s%wsn_ij(2,1,1,1),s%wsn_ij(3,1,1,1),!70
     &    s%wsn_ij(1,2,1,1),s%wsn_ij(2,2,1,1),s%wsn_ij(3,2,1,1),!73
     &    s%hsn_ij(1,1,1,1),s%hsn_ij(2,1,1,1),s%hsn_ij(3,1,1,1),!76
     &    s%hsn_ij(1,2,1,1),s%hsn_ij(2,2,1,1),s%hsn_ij(3,2,1,1),!79
     &    real(s%nsn_ij(1,1,1)),real(s%nsn_ij(2,1,1)),!81
     &    s%fr_snow_ij(1,1,1), s%fr_snow_ij(2,1,1),!83
     &    s%Qf_ij(1,1),!84
     &    aClivepool_leaf,aClivepool_froot,aClivepool_wood,!87
     &    aCdeadpool_surfmet,aCdeadpool_surfstr,aCdeadpool_soilmet,!90
     &    aCdeadpool_soilstr,aCdeadpool_cwd,aCdeadpool_surfmic,!93
     &    aCdeadpool_soilmic,aCdeadpool_slow,aCdeadpool_passive,!96
     &    alai,s%w_ij(0,2,1,1),s%ht_ij(0,2,1,1), !99
     &    acnc,acna, !101
     &    force_vis_rad(1,1),force_direct_vis_rad(1,1),tp(0,2), !104  
     &    tp(1,1),tp(2,1),tp(3,1),tp(4,1),tp(5,1),tp(6,1), !110
     &    tp(1,2),tp(2,2),tp(3,2),tp(4,2),tp(5,2),tp(6,2)  !116
#else


#ifdef PRINT_DIAGS
          n_count_mnth(k_mnth) = n_count_mnth(k_mnth) + 1.d0
#endif

      loop_j: do j = j0,j1  ! do j = 1,jm
        loop_i: do i = 1,im
          if ( bc%fearth(i,j) <= 0.d0 ) cycle loop_i
          if (i==35 .and. j==84) cycle loop_i ! not a land cell for GSWP2 data
!          write(933,*) "lsm_run: pricessing i,j ", i, j

          ! really fb, fv are not needed for Ent, but just in case...
          call ent_get_exports( s%entcells(i,j),
     &         fraction_of_vegetated_soil=fv
     &         )
          fb = 1.d0 - fv
!          print *, 'i,j,s%w_ij(1,1,i,j),s%w_ij(1,2,i,j)
!     &              ,s%ht_ij(1,1,i,j),s%ht_ij(1,2,i,j)
!     &              ,force_trheat(i,j)'
!     &                ,i,j
!     &                ,s%w_ij(1,1,i,j),s%w_ij(1,2,i,j)
!     &                ,s%ht_ij(1,1,i,j),s%ht_ij(1,2,i,j)
!     &                ,force_trheat(i,j)

          call advnc(
!-------------- Ent specific
     &         s%entcells(i,j), force_Ca(i,j),
     &         force_cos_zen_angle(i,j), force_vis_rad(i,j),
     &         force_direct_vis_rad(i,j),
     &         s%Qf_ij(i,j),
!-------------- old vegetation scheme (not implemented at the moment)
!     &         vegcell,
!-------------- prognostic vars
     &         s%w_ij(0:ngm,1:LS_NFRAC,i,j),         
     &         s%ht_ij(0:ngm,1:LS_NFRAC,i,j),        
     &         s%nsn_ij    (1:2, i, j),              
     &         s%dzsn_ij   (1:nlsn, 1:2, i, j),      
     &         s%wsn_ij    (1:nlsn, 1:2, i, j),      
     &         s%hsn_ij    (1:nlsn, 1:2, i, j),      
     &         s%fr_snow_ij(1:2, i, j),          
!-------------- BC's    
     &         bc%top_index_ij(i, j),                 
     &         bc%top_dev_ij(i, j),                   
     &         bc%dz_ij(i,j,1:ngm),                   
     &         bc%q_ij(i,j,1:imt,1:ngm),              
     &         bc%qk_ij(i,j,1:imt,1:ngm),             
     &         bc%sl_ij(i,j),                         
     &         fb,                                 
     &         fv,                        
!-------------- forcings         
     &         force_prec_ms (i,j)        ,
     &         force_eprec_w (i,j)        ,
     &         force_sprec_ms (i,j)       ,
     &         force_seprec_w (i,j)       ,
     &         0.d0,
     &         0.d0,
     &         force_srheat (i,j)         ,
     &         force_trheat (i,j)         ,
     &         force_ts (i,j)             ,
     &         force_qs (i,j)             ,
     &         force_ps (i,j)             ,
     &         force_rhosrf (i,j)         ,
     &         force_cdh (i,j)            ,
     &         force_qm1 (i,j)            ,
     &         force_ws (i,j)             ,
     &         force_pbl_args_ws0 (i,j)   ,
     &         force_pbl_args_tprime (i,j),
     &         force_pbl_args_qprime (i,j),
     &         end_of_day_flag )

!         Assign surface temperature for next timestep's Ch calculation  
          tbcs_ij(i,j)=tbcs
          tcanopy_ij(i,j)=tp(0,2)

!         Assign accumulators to output diagnostocs
#ifdef PRINT_DIAGS
          aevap_mnth(i,j,k_mnth)  = aevap_mnth(i,j,k_mnth) + aevap
          aevapw_mnth(i,j,k_mnth) = aevapw_mnth(i,j,k_mnth) + aevapw
          aevapd_mnth(i,j,k_mnth) = aevapd_mnth(i,j,k_mnth) + aevapd
          aevapb_mnth(i,j,k_mnth) = aevapb_mnth(i,j,k_mnth) + aevapb
          aintercep_mnth(i,j,k_mnth) = aintercep_mnth(i,j,k_mnth) + 
     &                                 aintercep
          aruns_mnth(i,j,k_mnth)  = aruns_mnth(i,j,k_mnth) + aruns
          arunu_mnth(i,j,k_mnth)  = arunu_mnth(i,j,k_mnth) + arunu
          agpp_mnth(i,j,k_mnth)   = agpp_mnth(i,j,k_mnth) + agpp
          arauto_mnth(i,j,k_mnth) = arauto_mnth(i,j,k_mnth) + arauto
          asoilresp_mnth(i,j,k_mnth)= asoilresp_mnth(i,j,k_mnth)+ 
     &                                asoilresp
          aevapvg_mnth(i,j,k_mnth)  = aevapvg_mnth(i,j,k_mnth) + aevapvg
          aevapvs_mnth(i,j,k_mnth)  = aevapvs_mnth(i,j,k_mnth) + aevapvs
          aevapbs_mnth(i,j,k_mnth)  = aevapbs_mnth(i,j,k_mnth) + aevapbs

          asrht_mnth(i,j,k_mnth) = asrht_mnth(i,j,k_mnth) + asrht
          atrht_mnth(i,j,k_mnth) = atrht_mnth(i,j,k_mnth) + atrht
          aalbedo_mnth(i,j,k_mnth) = aalbedo_mnth(i,j,k_mnth) + aalbedo

          aclab_mnth(i,j,k_mnth) = aclab_mnth(i,j,k_mnth) + aclab
          alai_m(i,j,k_mnth) = alai_m(i,j,k_mnth) + alai

          asoilCpoolsum_mnth(i,j,k_mnth) = 
     &          asoilCpoolsum_mnth(i,j,k_mnth) + asoilCpoolsum
          aepp_mnth(i,j,k_mnth) = aepp_mnth(i,j,k_mnth) + aepp
          atrg_mnth(i,j,k_mnth) = atrg_mnth(i,j,k_mnth) + atrg
          ashg_mnth(i,j,k_mnth) = ashg_mnth(i,j,k_mnth) + ashg
          alhg_mnth(i,j,k_mnth) = alhg_mnth(i,j,k_mnth) + alhg
          aepc_mnth(i,j,k_mnth) = aepc_mnth(i,j,k_mnth) + aepc
          abetad_mnth(i,j,k_mnth) = abetad_mnth(i,j,k_mnth) + abetad

          aClivepool_leaf_m(i,j,k_mnth) = 
     &        aClivepool_leaf_m(i,j,k_mnth)  + aClivepool_leaf
          aClivepool_froot_m(i,j,k_mnth) = 
     &        aClivepool_froot_m(i,j,k_mnth) + aClivepool_froot
          aClivepool_wood_m(i,j,k_mnth) =
     &        aClivepool_wood_m(i,j,k_mnth)  + aClivepool_wood 

          aCdeadpool_surfmet_m(i,j,k_mnth) =
     &        aCdeadpool_surfmet_m(i,j,k_mnth) + aCdeadpool_surfmet 
          aCdeadpool_surfstr_m(i,j,k_mnth) =
     &        aCdeadpool_surfstr_m(i,j,k_mnth) + aCdeadpool_surfstr 
          aCdeadpool_soilmet_m(i,j,k_mnth) =
     &        aCdeadpool_soilmet_m(i,j,k_mnth) + aCdeadpool_soilmet 
          aCdeadpool_soilstr_m(i,j,k_mnth) =
     &        aCdeadpool_soilstr_m(i,j,k_mnth) + aCdeadpool_soilstr 
          aCdeadpool_cwd_m(i,j,k_mnth) =
     &        aCdeadpool_cwd_m(i,j,k_mnth) + aCdeadpool_cwd 
          aCdeadpool_surfmic_m(i,j,k_mnth) =
     &        aCdeadpool_surfmic_m(i,j,k_mnth) + aCdeadpool_surfmic 
          aCdeadpool_soilmic_m(i,j,k_mnth) =
     &        aCdeadpool_soilmic_m(i,j,k_mnth) + aCdeadpool_soilmic 
          aCdeadpool_slow_m(i,j,k_mnth) =
     &        aCdeadpool_slow_m(i,j,k_mnth) + aCdeadpool_slow 
          aCdeadpool_passive_m(i,j,k_mnth) =
     &        aCdeadpool_passive_m(i,j,k_mnth) + aCdeadpool_passive 

          ! Accumulate prognostic variables
          w_b1_mnth(i,j,k_mnth)=w_b1_mnth(i,j,k_mnth)+s%w_ij(1,1,i,j)!Bare soil
          w_b2_mnth(i,j,k_mnth)=w_b2_mnth(i,j,k_mnth)+s%w_ij(2,1,i,j)!Bare soil
          w_b3_mnth(i,j,k_mnth)=w_b3_mnth(i,j,k_mnth)+s%w_ij(3,1,i,j)!Bare soil
          w_b4_mnth(i,j,k_mnth)=w_b4_mnth(i,j,k_mnth)+s%w_ij(4,1,i,j)!Bare soil
          w_b5_mnth(i,j,k_mnth)=w_b5_mnth(i,j,k_mnth)+s%w_ij(5,1,i,j)!Bare soil
          w_b6_mnth(i,j,k_mnth)=w_b6_mnth(i,j,k_mnth)+s%w_ij(6,1,i,j)!Bare soil

          w_v1_mnth(i,j,k_mnth)=w_v1_mnth(i,j,k_mnth)+s%w_ij(1,2,i,j)!Veg soil
          w_v2_mnth(i,j,k_mnth)=w_v2_mnth(i,j,k_mnth)+s%w_ij(2,2,i,j)!Veg soil
          w_v3_mnth(i,j,k_mnth)=w_v3_mnth(i,j,k_mnth)+s%w_ij(3,2,i,j)!Veg soil
          w_v4_mnth(i,j,k_mnth)=w_v4_mnth(i,j,k_mnth)+s%w_ij(4,2,i,j)!Veg soil
          w_v5_mnth(i,j,k_mnth)=w_v5_mnth(i,j,k_mnth)+s%w_ij(5,2,i,j)!Veg soil
          w_v6_mnth(i,j,k_mnth)=w_v6_mnth(i,j,k_mnth)+s%w_ij(6,2,i,j)!Veg soil

          ht_b1_mnth(i,j,k_mnth)=ht_b1_mnth(i,j,k_mnth)+s%ht_ij(1,1,i,j)!Bare
          ht_b2_mnth(i,j,k_mnth)=ht_b2_mnth(i,j,k_mnth)+s%ht_ij(2,1,i,j)!Bare
          ht_b3_mnth(i,j,k_mnth)=ht_b3_mnth(i,j,k_mnth)+s%ht_ij(3,1,i,j)!Bare
          ht_b4_mnth(i,j,k_mnth)=ht_b4_mnth(i,j,k_mnth)+s%ht_ij(4,1,i,j)!Bare
          ht_b5_mnth(i,j,k_mnth)=ht_b5_mnth(i,j,k_mnth)+s%ht_ij(5,1,i,j)!Bare
          ht_b6_mnth(i,j,k_mnth)=ht_b6_mnth(i,j,k_mnth)+s%ht_ij(6,1,i,j)!Bare

          ht_v1_mnth(i,j,k_mnth)=ht_v1_mnth(i,j,k_mnth)+s%ht_ij(1,2,i,j)!Veg
          ht_v2_mnth(i,j,k_mnth)=ht_v2_mnth(i,j,k_mnth)+s%ht_ij(2,2,i,j)!Veg
          ht_v3_mnth(i,j,k_mnth)=ht_v3_mnth(i,j,k_mnth)+s%ht_ij(3,2,i,j)!Veg
          ht_v4_mnth(i,j,k_mnth)=ht_v4_mnth(i,j,k_mnth)+s%ht_ij(4,2,i,j)!Veg
          ht_v5_mnth(i,j,k_mnth)=ht_v5_mnth(i,j,k_mnth)+s%ht_ij(5,2,i,j)!Veg
          ht_v6_mnth(i,j,k_mnth)=ht_v6_mnth(i,j,k_mnth)+s%ht_ij(6,2,i,j)!Veg

          dzsn_b1(i,j,k_mnth) = dzsn_b1(i,j,k_mnth)+s%dzsn_ij(1,1,i,j)
          dzsn_b2(i,j,k_mnth) = dzsn_b2(i,j,k_mnth)+s%dzsn_ij(2,1,i,j)
          dzsn_b3(i,j,k_mnth) = dzsn_b3(i,j,k_mnth)+s%dzsn_ij(3,1,i,j)
          dzsn_v1(i,j,k_mnth) = dzsn_v1(i,j,k_mnth)+s%dzsn_ij(1,2,i,j)
          dzsn_v2(i,j,k_mnth) = dzsn_v2(i,j,k_mnth)+s%dzsn_ij(2,2,i,j)
          dzsn_v3(i,j,k_mnth) = dzsn_v3(i,j,k_mnth)+s%dzsn_ij(3,2,i,j)

          wsn_b1(i,j,k_mnth) = wsn_b1(i,j,k_mnth)+s%wsn_ij(1,1,i,j)
          wsn_b2(i,j,k_mnth) = wsn_b2(i,j,k_mnth)+s%wsn_ij(2,1,i,j)
          wsn_b3(i,j,k_mnth) = wsn_b3(i,j,k_mnth)+s%wsn_ij(3,1,i,j)
          wsn_v1(i,j,k_mnth) = wsn_v1(i,j,k_mnth)+s%wsn_ij(1,2,i,j)
          wsn_v2(i,j,k_mnth) = wsn_v2(i,j,k_mnth)+s%wsn_ij(2,2,i,j)
          wsn_v3(i,j,k_mnth) = wsn_v3(i,j,k_mnth)+s%wsn_ij(3,2,i,j)

          hsn_b1(i,j,k_mnth) = hsn_b1(i,j,k_mnth)+s%hsn_ij(1,1,i,j)
          hsn_b2(i,j,k_mnth) = hsn_b2(i,j,k_mnth)+s%hsn_ij(2,1,i,j)
          hsn_b3(i,j,k_mnth) = hsn_b3(i,j,k_mnth)+s%hsn_ij(3,1,i,j)
          hsn_v1(i,j,k_mnth) = hsn_v1(i,j,k_mnth)+s%hsn_ij(1,2,i,j)
          hsn_v2(i,j,k_mnth) = hsn_v2(i,j,k_mnth)+s%hsn_ij(2,2,i,j)
          hsn_v3(i,j,k_mnth) = hsn_v3(i,j,k_mnth)+s%hsn_ij(3,2,i,j)

          nsn_b(i,j,k_mnth) = nsn_b(i,j,k_mnth)+ s%nsn_ij(1,i,j)
          nsn_v(i,j,k_mnth) = nsn_v(i,j,k_mnth)+ s%nsn_ij(2,i,j)

          frsnow_b(i,j,k_mnth)=frsnow_b(i,j,k_mnth)+s%fr_snow_ij(1,i,j)
          frsnow_v(i,j,k_mnth)=frsnow_v(i,j,k_mnth)+s%fr_snow_ij(2,i,j)

          Qf_mnth(i,j,k_mnth)=Qf_mnth(i,j,k_mnth)+ s%Qf_ij(i,j)
          canopyH2O_m(i,j,k_mnth)=canopyH2O_m(i,j,k_mnth) + 
     &                                   s%w_ij(0,2,i,j)
          canopyheat_m(i,j,k_mnth)=canopyheat_m(i,j,k_mnth) +
     &                                   s%ht_ij(0,2,i,j)
#endif
!          write(936,*) i,j, tbcs, s%w_ij(0:ngm,1:LS_NFRAC,i,j),
!     &     s%ht_ij(0:ngm,1:LS_NFRAC,i,j)

        enddo loop_i
      enddo loop_j


#endif

      !call sysusage(mype+8,2)
      end subroutine lsm_run

!-----------------------------------------------------------------------

      subroutine set_vegetation_data( entcells,
     &     im, jm, i0, i1, j0, j1, jday, year )
!@sum read standard GISS vegetation BC's and pass them to Ent for
!@+   initialization of Ent cells. Halo cells ignored, i.e.
!@+   entcells should be a slice without halo
      use ent_prescribed_drv, only : init_canopy_physical,prescr_vegdata
      use ent_prescr_veg, only : init_params !old name:prescr_calcconst
      type(entcelltype_public), intent(out) :: entcells(I0:I1,J0:J1)
      integer, intent(in) :: im, jm, i0, i1, j0, j1, jday, year
      !Local variables
      real*8, dimension(N_COVERTYPES,I0:I1,J0:J1) :: vegdata !cohort
      real*8, dimension(N_BANDS,N_COVERTYPES,I0:I1,J0:J1) :: albedodata !patch, NOTE:snow
      real*8, dimension(N_COVERTYPES,I0:I1,J0:J1) :: laidata  !cohort
      real*8, dimension(N_COVERTYPES,I0:I1,J0:J1) :: hdata    !cohort
      real*8, dimension(N_COVERTYPES) :: nmdata    !cohort
      real*8, dimension(N_COVERTYPES,N_DEPTH) :: rootprofdata !Root fraction of veg type.
      real*8, dimension(N_COVERTYPES,I0:I1,J0:J1) :: popdata !Dummy population density:  0-bare soil, 1-vegetated
      real*8, dimension(N_COVERTYPES,I0:I1,J0:J1) :: dbhdata !Diameter at breast height for woody veg.(cm)
      real*8, dimension(N_COVERTYPES,I0:I1,J0:J1) :: craddata !Crown radius (m)
      real*8, dimension(N_COVERTYPES,N_BPOOLS,I0:I1,J0:J1) :: cpooldata !Carbon pools in individuals
      integer, dimension(N_COVERTYPES) :: soildata ! soil types 1-bright 2-dark
      real*8, dimension(N_SOIL_TEXTURES,I0:I1,J0:J1) :: soil_texture
      real*8, dimension(I0:I1,J0:J1) :: Ci_ini,CNC_ini,Tcan_ini,Qf_ini
      real*8, dimension(N_PFT,PTRACE,NPOOLS-NLIVE,N_CASA_LAYERS,
     &     I0:I1,J0:J1):: Tpool_ini

      !-----Local---------
      integer i,j
      real*8 heat_capacity

      call init_params() !old name: prescr_calcconst()

      call prescr_vegdata(jday, year, 
     &     IM,JM,I0,I1,J0,J1,vegdata,albedodata,laidata,hdata,nmdata,
     &     popdata,dbhdata,craddata,cpooldata,rootprofdata,
     &     soildata,soil_texture,Tpool_ini,
     &     do_soilinit,do_phenology_activegrowth,do_init_geo,.true.)
      call init_canopy_physical(I0, I1, J0, J1,
     &     Ci_ini, CNC_ini, Tcan_ini, Qf_ini)

      !Translate gridded data to Entdata structure
      !GISS data:  a patch per vegetation cover fraction, one cohort per patch
      call ent_cell_set(entcells, vegdata, popdata, laidata,
     &     hdata, dbhdata, craddata, cpooldata, nmdata, rootprofdata, 
     &     soildata, albedodata, soil_texture,
     &     Ci_ini, CNC_ini, Tcan_ini, Qf_ini, Tpool_ini, .true.)
 
      ! just in case, do nothing, just set heat capacities
      call ent_prescribe_vegupdate(entcells)

      end subroutine set_vegetation_data

!-----------------------------------------------------------------------

      subroutine update_vegetation_data( entcells,
     &     im, jm, i0, i1, j0, j1, jday, year )
!@sum read standard GISS vegetation BC's and pass them to Ent for
!@+   initialization of Ent cells. Halo cells ignored, i.e.
!@+   entcells should be a slice without halo
      type(entcelltype_public), intent(out) :: entcells(I0:I1,J0:J1)
      integer, intent(in) :: im, jm, i0, i1, j0, j1, jday, year
      !Local variables
      real*8, dimension(N_BANDS,N_COVERTYPES,I0:I1,J0:J1) :: albedodata !patch, NOTE:snow
      real*8, dimension(N_COVERTYPES,I0:I1,J0:J1) :: laidata  !cohort
      !-----Local---------
      integer hemi(I0:I1,J0:J1)
      integer i,j
      real*8 , dimension(N_COVERTYPES):: LAI_buf

!Update Ent exactly like in ent_prog

            !* Set hemisphere flags.
      if ( J0<=JM/2 )   hemi(:,J0:min(JM/2,J1))   = -1    ! S.
      if ( J1>=JM/2+1 ) hemi(:,max(JM/2+1,J0):J1) =  1    ! N.
      if (force_VEG) then
          call ent_prescribe_vegupdate(entcells,hemi,jday,year,
     &         do_giss_phenology=
     &         (.not.do_phenology_activegrowth).and.(.not.force_VEG),
     &         do_giss_albedo=.true.,
     &         do_giss_lai=
     &         (.not.do_phenology_activegrowth).and.(.not.force_VEG),
     &         update_crops=.false.)
       else !Don't include laidata and hdata in parameter list.
          call ent_prescribe_vegupdate(entcells,hemi,jday,year,
     &         do_giss_phenology=
     &         (.not.do_phenology_activegrowth).and.(.not.force_VEG),
     &         do_giss_albedo=.true.,
     &         do_giss_lai=
     &         (.not.do_phenology_activegrowth).and.(.not.force_VEG),
     &         update_crops=.false.)
       endif

      return

      end subroutine update_vegetation_data

!-----------------------------------------------------------------------
#ifdef ECOSYSTEM_SCALE
      subroutine set_veg_data_single( entcells,
     &     im, jm, isite, jsite, jday, year,laidata) 
!@sum read standard GISS vegetation BC's and pass them to Ent for
!@+   initialization of Ent cells. 
      use ent_prescribed_drv, only : init_canopy_physical,prescr_vegdata
      use ent_prescr_veg, only : prescr_calcconst
      type(entcelltype_public), intent(out) :: entcells(1:1,1:1)
      integer, intent(in) :: im, jm, isite, jsite, jday, year
      !Local variables
      real*8, dimension(N_COVERTYPES,1:im,1:jm) :: vegdata !cohort
      real*8, dimension(N_BANDS,N_COVERTYPES,1:im,1:jm) :: albedodata !patch, NOTE:snow
      real*8, dimension(N_COVERTYPES,1:im,1:jm) :: laidata  !cohort
      real*8, dimension(N_COVERTYPES,1:im,1:jm) :: hdata    !cohort
      real*8, dimension(N_COVERTYPES) :: nmdata    !cohort
      real*8, dimension(N_COVERTYPES,N_DEPTH) :: rootprofdata !Root fraction of veg type.
      real*8, dimension(N_COVERTYPES,1:im,1:jm) :: popdata !Dummy population density:  0-bare soil, 1-vegetated
      real*8, dimension(N_COVERTYPES,1:im,1:jm) :: dbhdata !Diameter at breast height for woody veg.(cm)
      real*8, dimension(N_COVERTYPES,1:im,1:jm) :: craddata !Crown radius (m)
      real*8, dimension(N_COVERTYPES,N_BPOOLS,1:im,1:jm) :: 
     &                 cpooldata !Carbon pools in individuals
      integer, dimension(N_COVERTYPES) :: soildata ! soil types 1-bright 2-dark
      real*8, dimension(N_SOIL_TEXTURES,1:im,1:jm) ::
     &           soil_texture
      real*8, dimension(1:im,1:jm) :: Ci_ini,CNC_ini,Tcan_ini,Qf_ini
      real*8, dimension(N_PFT,PTRACE,NPOOLS-NLIVE,N_CASA_LAYERS,
     &     1:im,1:jm):: Tpool_ini

      real*8, dimension(N_COVERTYPES,1,1) :: vegdata_single !cohort
      real*8, dimension(N_BANDS,N_COVERTYPES,1,1) :: albedodata_single !patch, NOTE:snow
      real*8, dimension(N_COVERTYPES,1,1) :: laidata_single  !cohort
      real*8, dimension(N_COVERTYPES,1,1) :: hdata_single    
      real*8, dimension(N_COVERTYPES,1,1) :: popdata_single  
      real*8, dimension(N_COVERTYPES,1,1) :: dbhdata_single  
      real*8, dimension(N_COVERTYPES,1,1) :: craddata_single 
      real*8, dimension(N_COVERTYPES,N_BPOOLS,1,1) :: cpooldata_single !Carbon pools in individuals
      real*8, dimension(N_SOIL_TEXTURES,1,1) :: soil_texture_single
      real*8, dimension(1,1) :: Ci_ini_single,CNC_ini_single,
     &                          Tcan_ini_single,Qf_ini_single
      real*8, dimension(N_PFT,PTRACE,NPOOLS-NLIVE,N_CASA_LAYERS,
     &     1,1):: Tpool_ini_single
    
      !-----Local---------
      integer i,j
      real*8 heat_capacity

      call prescr_calcconst()

      call prescr_vegdata(jday, year, 
     &     IM,JM,1,im,1,jm,vegdata,albedodata,laidata,hdata,nmdata,
     &     popdata,dbhdata,craddata,cpooldata,rootprofdata,
     &     soildata,soil_texture,Tpool_ini,
     &     do_soilinit,do_phenology_activegrowth,do_init_geo,.true.)

      call init_canopy_physical(1, im, 1, jm,
     &     Ci_ini, CNC_ini, Tcan_ini, Qf_ini)

      vegdata_single(:,1,1)=vegdata(:, isite, jsite)
      albedodata_single(:,:,1,1)=albedodata(:,:,isite, jsite)
      laidata_single(:,1,1)=laidata(:,isite, jsite)
      hdata_single(:,1,1)=hdata(:,isite, jsite)
      popdata_single(:,1,1)=popdata(:,isite, jsite)
      dbhdata_single(:,1,1)=dbhdata(:,isite, jsite)
      craddata_single(:,1,1)=craddata(:,isite, jsite)
      cpooldata_single(:,:,1,1)=cpooldata(:,:,isite, jsite)
      soil_texture_single(:,1,1)=soil_texture(:,i_site, j_site)
      Tpool_ini_single(:,:,:,:,1,1)=Tpool_ini(:,:,:,:, isite, jsite)
      Ci_ini_single(1,1)=Ci_ini(isite, jsite)
      CNC_ini_single(1,1)=CNC_ini(isite, jsite)
      Tcan_ini_single(1,1)=Tcan_ini(isite, jsite)
      Qf_ini_single(1,1)=Qf_ini(isite, jsite)

      !Translate gridded data to Entdata structure
      !GISS data:  a patch per vegetation cover fraction, one cohort per patch

      call ent_cell_set(entcells,vegdata_single,popdata_single,
     &     laidata_single,hdata_single,dbhdata_single,craddata_single,
     &     cpooldata_single,nmdata,rootprofdata, 
     &     soildata, albedodata_single, soil_texture_single,
     &     Ci_ini_single, CNC_ini_single, Tcan_ini_single, 
     &     Qf_ini_single, Tpool_ini_single, .true.)
 
      ! just in case, do nothing, just set heat capacities
      call ent_prescribe_vegupdate(entcells)

      end subroutine set_veg_data_single

!-----------------------------------------------------------------------

      subroutine update_veg_data_single( entcells
     &     , im, jm, isite, jsite, jday, year, laidata, veg_height) 

!@sum read standard GISS vegetation BC's and pass them to Ent for
!@+   initialization of Ent cells.
      use ent_prescribed_drv, only:
     &     prescr_get_laidata,prescr_veg_albedodata
      use ent_pfts
      !use ent_prescr_veg, only: prescr_get_laidata,prescr_veg_albedodata
      type(entcelltype_public), intent(out) :: entcells(1:1,1:1)
      integer, intent(in) :: im, jm, isite, jsite, jday, year
      !Local variables
      real*8, dimension(N_BANDS,N_COVERTYPES,1:im,1:jm) :: albedodata !patch, NOTE:snow
      real*8, dimension(N_COVERTYPES,1,1) :: laidata  !cohort
      real*8, dimension(N_COVERTYPES,1,1) :: veg_height
      !-----Local---------
      integer hemi(1,1)
      integer i,j
      real*8, dimension(N_BANDS,N_COVERTYPES,1,1) :: albedodata_single !patch

      if ( jsite<=JM/2 )   hemi(:,1) = -1    ! S.
      if ( jsite>=JM/2+1 ) hemi(:,1) =  1    ! N.

      if (force_VEG) then
      call ent_prescribe_vegupdate(entcells,hemi,jday,year,
     &         do_giss_phenology=
     &         (.not.do_phenology_activegrowth).and.(.not.force_VEG),
     &         do_giss_albedo=.true.,
     &         do_giss_lai=
     &         (.not.do_phenology_activegrowth).and.(.not.force_VEG),
     &         update_crops=.false.
     &        ,laidata=laidata(COVEROFFSET+1:COVEROFFSET+N_PFT,:,:)
     &        ,hdata=veg_height(COVEROFFSET+1:COVEROFFSET+N_PFT,:,:))
      else
      call ent_prescribe_vegupdate(entcells,hemi,jday,year,
     &         do_giss_phenology=
     &         (.not.do_phenology_activegrowth).and.(.not.force_VEG),
     &         do_giss_albedo=.true.,
     &         do_giss_lai=
     &         (.not.do_phenology_activegrowth).and.(.not.force_VEG),
     &         update_crops=.false.)
      endif

      return
      end subroutine update_veg_data_single

!-----------------------------------------------------------------------

      subroutine update_FLUXNET_LAI( iu_LAI,iu_vht,entcells
     &     , im, jm, isite, jsite, jday, year,laidata,veg_height )
!@sum read FLUXNET LAI data and pass to ent
      use ent_pfts
      type(entcelltype_public), intent(out) :: entcells(1:1,1:1)
      integer, intent(in) :: im, jm, isite, jsite, jday, year
      real*8, dimension(N_COVERTYPES,1,1) :: laidata  !cohort
      real*8, dimension(N_COVERTYPES,1,1):: veg_height
      !-----Local---------
      integer hemi(1,1)
      integer i,j
      integer :: iu_LAI, iu_vht
      real*8 , dimension(N_COVERTYPES):: LAI_buf, vht_buf

!!!! HACK : trying to update Ent exactly like in ent_prog
      !* Set hemisphere flags.
      if ( jsite<=JM/2 )   hemi(:,1) = -1    ! S.
      if ( jsite>=JM/2+1 ) hemi(:,1) =  1    ! N.

!     Read Leaf Area Index data
      read(iu_LAI,*) LAI_buf

!     Read vegetation height data      
!      rewind(iu_vht) !KIM- height is the timeseries now!!!
      read(iu_vht,*) vht_buf(:)

      do i=1,N_COVERTYPES
         laidata(i,1,1) = LAI_buf(i)
         veg_height(i,1,1) = vht_buf(i)
      end do

      call ent_prescribe_vegupdate(entcells,hemi,jday,year,
     &         do_giss_phenology=
     &         (.not.do_phenology_activegrowth).and.(.not.force_VEG),
     &         do_giss_albedo=.true.,
     &         do_giss_lai=
     &         (.not.do_phenology_activegrowth).and.(.not.force_VEG),
     &         update_crops=.false.
     &        ,laidata=laidata(COVEROFFSET+1:COVEROFFSET+N_PFT,:,:)
     &        ,hdata=veg_height(COVEROFFSET+1:COVEROFFSET+N_PFT,:,:))

      !print *, 'LAI=',laidata(6,1,1)

      end subroutine update_FLUXNET_LAI
#endif
!-----------------------------------------------------------------------

      subroutine update_water_heat( s,bc,iu_water)
!@sum read climatological soil moist values daily and compute heat in layers
!!!      use ghy_com, only : q_ij, dz_ij, w_ij, ht_ij
      use sle001, only : thm

      type(t_lsm_state) :: s
      type (t_lsm_bc) :: bc

      integer,dimension(12) :: iu_water
      real*8,parameter::rhow = 1000.d0
      real*8,parameter::lhm = 334590.d0 !J/kg latent heat of melting at 0 deg C
      real*8,parameter::shw_kg = 4185.d0 !J/kg/deg C heat cap. of water @ 20 C
      real*8,parameter::shi_kg = 2060.d0 !J/kg/deg C heat cap. of pure ice @0C
      real*8, parameter :: shc_soil_texture(imt)
     &     = (/2d6,2d6,2d6,2.5d6,2.4d6/)!spec. heat cap. of soil text.(J/K/M^3)

      real*8 ::shwv,shiv,shc_layer
      real*8,dimension(ngm) :: ht_cap,w_stor
      real*8,dimension(ngm,LS_NFRAC-1,im,jm) :: w_ij_temp, ht_ij_temp
#ifdef ECOSYSTEM_SCALE
      real*8,dimension(ngm,LS_NFRAC-1) ::w_old
#else
      real*8,dimension(ngm,LS_NFRAC-1,im,jm) ::w_old
#endif
      real*8 ::tp,fice
      integer :: i,j,k,m,mm
      ! volumetric quantities
      shwv=shw_kg*rhow
      shiv=shi_kg*rhow

      do k=1,ngm
         do i = 1,im
            ! input soil moisture values in mm; convert to m below
            read(iu_water(k),*) (w_ij_temp(k,1,i,j),j=1,jm)
            read(iu_water(k+6),*) (w_ij_temp(k,2,i,j),j=1,jm)
         enddo
      enddo

      do k=1,ngm
         do m=1,(LS_NFRAC-1)
#ifdef ECOSYSTEM_SCALE
            w_old(k,m)=s%w_ij(k,m,1,1)
            s%w_ij(k,m,1,1)= w_ij_temp(k,m,i_site,j_site)/1000.d0
#else
            w_old(k,m,:,:)=s%w_ij(k,m,:,:)
            s%w_ij(k,m,:,:)= w_ij_temp(k,m,:,:)/1000.d0
#endif
         enddo
      enddo
      ! set soil moisture under lakes to bare soil moisture
#ifdef ECOSYSTEM_SCALE
      s%w_ij(:,3,1,1)= s%w_ij(:,1,1,1)
#else
      s%w_ij(:,3,:,:)= s%w_ij(:,1,:,:)
#endif

#ifdef ECOSYSTEM_SCALE
      do k=1,ngm
        ! compute max water storage and heat capacity
         w_stor(k) = 0.d0
         do mm=1,imt-1
             w_stor(k) = w_stor(k) + 
     &              bc%q_ij(1,1,mm,k)*thm(0,mm)*bc%dz_ij(1,1,k)
         enddo
         shc_layer = 0.d0
         do mm=1,imt
            shc_layer = shc_layer + 
     &              bc%q_ij(1,1,mm,k)*shc_soil_texture(mm)
         enddo
         ht_cap(k) = (bc%dz_ij(1,1,k)-w_stor(k)) * shc_layer

         ! loop over land fractions
         do m=1,(LS_NFRAC-1)
            call heat_to_temperature( tp, fice,
     &            s%ht_ij(k,m,1,1), w_old(k,m), ht_cap(ngm) )

            s%ht_ij(k,m,1,1)= s%ht_ij(k,m,1,1)
     &            + (fice*shiv+(1.d0-fice)*shwv)*
     &              (s%w_ij(k,m,1,1)- w_old(k,m)) * tp
         enddo
      enddo

#else
      do i=1,im
         do j=1,jm
            do k=1,ngm
              ! compute max water storage and heat capacity
               w_stor(k) = 0.d0
               do mm=1,imt-1
                   w_stor(k) = w_stor(k) + 
     &               bc%q_ij(i,j,mm,k)*thm(0,mm)*bc%dz_ij(i,j,k)
               enddo
               shc_layer = 0.d0
               do mm=1,imt
                  shc_layer = shc_layer + 
     &                 bc%q_ij(i,j,mm,k)*shc_soil_texture(mm)
               enddo
               ht_cap(k) = (bc%dz_ij(i,j,k)-w_stor(k)) * shc_layer

               do m=1,(LS_NFRAC-1)
                  call heat_to_temperature( tp, fice,
     &               s%ht_ij(k,m,i,j), w_old(k,m,i,j), ht_cap(ngm) )

                  s%ht_ij(k,m,i,j)= s%ht_ij(k,m,i,j)
     &                 + (fice*shiv+(1.d0-fice)*shwv)*
     &                 (s%w_ij(k,m,i,j)- w_old(k,m,i,j)) * tp

               enddo

            enddo
         enddo
      enddo
#endif

      end subroutine update_water_heat

!-----------------------------------------------------------------------

      subroutine heat_to_temperature(tp, fice, ht, w, ht_cap)
      real*8, intent(out) :: tp, fice
      real*8, intent(in) :: ht, w, ht_cap

      real*8,parameter::rhow = 1000.d0
      real*8,parameter::lhm = 334590.d0 !J/kg latent heat of melting at 0 deg C
      real*8,parameter::shw_kg = 4185.d0 !J/kg C heat capacity of water at 20 C
      real*8,parameter::shi_kg = 2060.d0 !J/kg heat capacity of pure ice at 0C
      ! volumetric quantities
      real*8, parameter :: lhmv=lhm*rhow, shwv=shw_kg*rhow,
     &     shiv=shi_kg*rhow

      fice = 0.d0
      if( lhmv*w+ht < 0.d0 ) then ! all frozen
        tp = ( ht + w*lhmv )/( ht_cap + w*shiv )
        fice = 1.d0
      else if( ht > 0.d0 ) then ! all melted
        tp = ht /( ht_cap + w*shwv )
      else  ! part frozen
        tp = 0.d0
        if( w > 1d-12 ) fice = -ht /( lhmv*w )
      endif

      end subroutine heat_to_temperature

!-----------------------------------------------------------------------
#ifdef PRINT_DIAGS  /*  Monthly diagnostic accumulators */

      subroutine print_diags(     
     &                   aevap_mnth,aevapw_mnth,aevapd_mnth,aevapb_mnth
     &                 , aintercep_mnth,aruns_mnth,arunu_mnth
     &                 , agpp_mnth,arauto_mnth,asoilresp_mnth
     &                 , aevapvg_mnth,aevapvs_mnth,aevapbs_mnth
     &                 , asrht_mnth,atrht_mnth,aalbedo_mnth
     &                 , aclab_mnth,asoilCpoolsum_mnth,aepp_mnth
     &                 , atrg_mnth,ashg_mnth,alhg_mnth,aepc_mnth
     &                 , n_count_mnth
     &                 , w_b1_mnth,w_b2_mnth,w_b3_mnth
     &                 , w_b4_mnth,w_b5_mnth,w_b6_mnth
     &                 , w_v1_mnth,w_v2_mnth,w_v3_mnth
     &                 , w_v4_mnth,w_v5_mnth,w_v6_mnth
     &                 , ht_b1_mnth,ht_b2_mnth,ht_b3_mnth
     &                 , ht_b4_mnth,ht_b5_mnth,ht_b6_mnth
     &                 , ht_v1_mnth,ht_v2_mnth,ht_v3_mnth
     &                 , ht_v4_mnth,ht_v5_mnth,ht_v6_mnth
     &                 , dzsn_b1,dzsn_b2,dzsn_b3,dzsn_v1,dzsn_v2,dzsn_v3
     &                 , wsn_b1 ,wsn_b2 ,wsn_b3 ,wsn_v1 ,wsn_v2 ,wsn_v3
     &                 , hsn_b1 ,hsn_b2 ,hsn_b3 ,hsn_v1 ,hsn_v2 ,hsn_v3
     &                 , nsn_b,nsn_v,frsnow_b,frsnow_v,Qf_mnth
     &                 , abetad_mnth
     &  , aClivepool_leaf_m,aClivepool_froot_m,aClivepool_wood_m
     &  , aCdeadpool_surfmet_m,aCdeadpool_surfstr_m,aCdeadpool_soilmet_m
     &  , aCdeadpool_soilstr_m,aCdeadpool_cwd_m,aCdeadpool_surfmic_m
     &  , aCdeadpool_soilmic_m,aCdeadpool_slow_m,aCdeadpool_passive_m
     &  , alai_m,canopyH2O_m,canopyheat_m
     &    )

      use domain_decomp, only : array_gather, mype

      character*80 :: title       ! title of diagnostic to print to file
      real*8, dimension(im,jm,12) ::
     &                    aevap_mnth,aevapw_mnth,aevapd_mnth,aevapb_mnth
     &                  , aintercep_mnth,aruns_mnth,arunu_mnth
     &                  , agpp_mnth,arauto_mnth,asoilresp_mnth
     &                  , aevapvg_mnth,aevapvs_mnth,aevapbs_mnth
     &                  , asrht_mnth,atrht_mnth,aalbedo_mnth
     &                  , aclab_mnth,asoilCpoolsum_mnth,aepp_mnth
     &                  , atrg_mnth,ashg_mnth,alhg_mnth,aepc_mnth
     &                  , w_b1_mnth,w_b2_mnth,w_b3_mnth
     &                  , w_b4_mnth,w_b5_mnth,w_b6_mnth
     &                  , w_v1_mnth,w_v2_mnth,w_v3_mnth
     &                  , w_v4_mnth,w_v5_mnth,w_v6_mnth
     &                  , ht_b1_mnth,ht_b2_mnth,ht_b3_mnth
     &                  , ht_b4_mnth,ht_b5_mnth,ht_b6_mnth
     &                  , ht_v1_mnth,ht_v2_mnth,ht_v3_mnth
     &                  , ht_v4_mnth,ht_v5_mnth,ht_v6_mnth
     &                 , dzsn_b1,dzsn_b2,dzsn_b3,dzsn_v1,dzsn_v2,dzsn_v3
     &                 , wsn_b1 ,wsn_b2 ,wsn_b3 ,wsn_v1 ,wsn_v2 ,wsn_v3
     &                 , hsn_b1 ,hsn_b2 ,hsn_b3 ,hsn_v1 ,hsn_v2 ,hsn_v3
     &                 , nsn_b,nsn_v,frsnow_b,frsnow_v,Qf_mnth
     &                 , abetad_mnth
     &  , aClivepool_leaf_m,aClivepool_froot_m,aClivepool_wood_m
     &  , aCdeadpool_surfmet_m,aCdeadpool_surfstr_m,aCdeadpool_soilmet_m
     &  , aCdeadpool_soilstr_m,aCdeadpool_cwd_m,aCdeadpool_surfmic_m
     &  , aCdeadpool_soilmic_m,aCdeadpool_slow_m,aCdeadpool_passive_m
     &  , alai_m,canopyH2O_m,canopyheat_m

      real*8,dimension(12),intent(in) :: n_count_mnth
      integer :: i_mnth      ! month counter


#ifdef USE_ESMF
          !print *, n_count_mnth
          do i_mnth=1,12
             title = 'GISS LSM: Total evaporation [kg/m2/month]'
             call array_gather( aevap_mnth(:,:,i_mnth) )
             if ( mype == 0 ) write(980) title,
     &                        real(aevap_mnth(:,:,i_mnth),kind=4)
             title = 'GISS LSM: Canopy evaporation (evapw)[kg/m2/month]'
             call array_gather( aevapw_mnth(:,:,i_mnth) )
             if ( mype == 0 ) write(981) title,
     &                        real(aevapw_mnth(:,:,i_mnth),kind=4)
             title = 'GISS LSM: Transpiration (evapd)[kg/m2/month]'
             call array_gather( aevapd_mnth(:,:,i_mnth) )
             if ( mype == 0 ) write(982) title,
     &                         real(aevapd_mnth(:,:,i_mnth),kind=4)
             title = 'GISS LSM: Bare soil evap. (evapb)[kg/m2/month]'
             call array_gather( aevapb_mnth(:,:,i_mnth) )
             if ( mype == 0 ) write(983) title,
     &                         real(aevapb_mnth(:,:,i_mnth),kind=4)
             title = 'GISS LSM: Canopy interception'
             call array_gather( aintercep_mnth(:,:,i_mnth) )
             if ( mype == 0 ) write(984) title,
     &                        real(aintercep_mnth(:,:,i_mnth),kind=4)
             title = 'GISS LSM: Surface runoff [kg/m2/mnth]'
             call array_gather( aruns_mnth(:,:,i_mnth) )
             if ( mype == 0 ) write(985) title,
     &                         real(aruns_mnth(:,:,i_mnth),kind=4)
             title = 'GISS LSM: Subsurface flow out of soil[kg/m2/mnth]'
             call array_gather( arunu_mnth(:,:,i_mnth) )
             if ( mype == 0 ) write(986) title,
     &                         real(arunu_mnth(:,:,i_mnth),kind=4)
             title = 'GISS LSM: Gross primary productivity[kgC/m2/mnth]'
             call array_gather( agpp_mnth(:,:,i_mnth) )
             if ( mype == 0 ) write(987) title,
     &                         real(agpp_mnth(:,:,i_mnth),kind=4)
             title = 'GISS LSM: Autotrophic respiration[kgC/m2/mnth]'
             call array_gather( arauto_mnth(:,:,i_mnth) )
             if ( mype == 0 ) write(988) title,
     &                         real(arauto_mnth(:,:,i_mnth),kind=4)
             title = 'GISS LSM: Soil respiration[kgC/m2/mnth]'
             call array_gather( asoilresp_mnth(:,:,i_mnth) )
             if ( mype == 0 ) write(989) title,
     &                         real(asoilresp_mnth(:,:,i_mnth),kind=4)
             title = 'GISS LSM: Evap. from vegetated soil[kg/m2/month]'
             call array_gather( aevapvg_mnth(:,:,i_mnth) )
             if ( mype == 0 ) write(990) title,
     &                         real(aevapvg_mnth(:,:,i_mnth),kind=4)
             title = 'GISS LSM: Evap. from snow on veg[kg/m2/month]'
             call array_gather( aevapvs_mnth(:,:,i_mnth) )
             if ( mype == 0 ) write(991) title,
     &                         real(aevapvs_mnth(:,:,i_mnth),kind=4)
             title = 'GISS LSM: Evap. from snow on bare soil[kg/m2/mn]'
             call array_gather( aevapbs_mnth(:,:,i_mnth) )
             if ( mype == 0 ) write(992) title,
     &                         real(aevapbs_mnth(:,:,i_mnth),kind=4)
             title = 'GISS LSM: Net shortwave radiation [W/m2]'
             call array_gather( asrht_mnth(:,:,i_mnth) )
             asrht_mnth(:,:,i_mnth) = 
     &                asrht_mnth(:,:,i_mnth)/n_count_mnth(i_mnth)
             if ( mype == 0 ) write(993) title,
     &                         real(asrht_mnth(:,:,i_mnth),kind=4)
             title = 'GISS LSM: Net longwave radiation [W/m2]'
             call array_gather( atrht_mnth(:,:,i_mnth) )
              atrht_mnth(:,:,i_mnth) = 
     &                atrht_mnth(:,:,i_mnth)/n_count_mnth(i_mnth)
             if ( mype == 0 ) write(994) title,
     &                         real(atrht_mnth(:,:,i_mnth),kind=4)
             title = 'GISS LSM: Grid-cell mean albedo'
              aalbedo_mnth(:,:,i_mnth) = 
     &                aalbedo_mnth(:,:,i_mnth)/n_count_mnth(i_mnth)
             call array_gather( aalbedo_mnth(:,:,i_mnth) )
             if ( mype == 0 ) write(995) title,
     &                          real(aalbedo_mnth(:,:,i_mnth),kind=4)
             title = 'GISS LSM: Plant labile carbon [kgC/m2]'
             aclab_mnth(:,:,i_mnth) = 
     &                aclab_mnth(:,:,i_mnth)/n_count_mnth(i_mnth)
             call array_gather( aclab_mnth(:,:,i_mnth) )
             if ( mype == 0 ) write(996) title,
     &                         real(aclab_mnth(:,:,i_mnth),kind=4)
             title = 'GISS LSM: Sum of soil carbon pools[kgC/m2]'
             asoilCpoolsum_mnth(:,:,i_mnth) = 
     &             asoilCpoolsum_mnth(:,:,i_mnth)/n_count_mnth(i_mnth)
             call array_gather( asoilCpoolsum_mnth(:,:,i_mnth) )
             if ( mype == 0 ) write(997) title,
     &                      real(asoilCpoolsum_mnth(:,:,i_mnth),kind=4)
             title = 'GISS LSM: Penman potential evap. [kg/m2/month]'
             call array_gather( aepp_mnth(:,:,i_mnth) )
             if ( mype == 0 ) write(998) title,
     &                      real(aepp_mnth(:,:,i_mnth),kind=4)
             title = 'GISS LSM: Thermal heat from ground [J/m2/mnth]'
             call array_gather( atrg_mnth(:,:,i_mnth) )
             if ( mype == 0 ) write(999) title,
     &                      real(atrg_mnth(:,:,i_mnth),kind=4)
             title = 'GISS LSM: Sensible heat from ground [J/m2/mnth]'
             call array_gather( ashg_mnth(:,:,i_mnth) )
             if ( mype == 0 ) write(1000) title,
     &                      real(ashg_mnth(:,:,i_mnth),kind=4)
             title = 'GISS LSM: Latent heat flux [W/m2]'
             call array_gather( alhg_mnth(:,:,i_mnth) )
             alhg_mnth(:,:,i_mnth) = 
     &                alhg_mnth(:,:,i_mnth)/n_count_mnth(i_mnth)
             if ( mype == 0 ) write(1001) title,
     &                      real(alhg_mnth(:,:,i_mnth),kind=4)
             title = 'GISS LSM: Pot. evap. from canopy [kg/m2/mnth]'
             call array_gather( aepc_mnth(:,:,i_mnth) )
             if ( mype == 0 ) write(1002) title,
     &                      real(aepc_mnth(:,:,i_mnth),kind=4)

             title = 'GISS LSM: Soil water; Bare layer 1 [m]'
             call array_gather( w_b1_mnth(:,:,i_mnth) )
             w_b1_mnth(:,:,i_mnth) = 
     &                w_b1_mnth(:,:,i_mnth)/n_count_mnth(i_mnth)
             if ( mype == 0 ) write(1003) title,
     &                      real(w_b1_mnth(:,:,i_mnth),kind=4)

             title = 'GISS LSM: Soil water; Bare layer 2 [m]'
             call array_gather( w_b2_mnth(:,:,i_mnth) )
             w_b2_mnth(:,:,i_mnth) = 
     &                w_b2_mnth(:,:,i_mnth)/n_count_mnth(i_mnth)
             if ( mype == 0 ) write(1004) title,
     &                      real(w_b2_mnth(:,:,i_mnth),kind=4)

             title = 'GISS LSM: Soil water; Bare layer 3 [m]'
             call array_gather( w_b3_mnth(:,:,i_mnth) )
             w_b3_mnth(:,:,i_mnth) = 
     &                w_b3_mnth(:,:,i_mnth)/n_count_mnth(i_mnth)
             if ( mype == 0 ) write(1005) title,
     &                      real(w_b3_mnth(:,:,i_mnth),kind=4)

             title = 'GISS LSM: Soil water; Bare layer 4 [m]'
             call array_gather( w_b4_mnth(:,:,i_mnth) )
             w_b4_mnth(:,:,i_mnth) = 
     &                w_b4_mnth(:,:,i_mnth)/n_count_mnth(i_mnth)
             if ( mype == 0 ) write(1006) title,
     &                      real(w_b4_mnth(:,:,i_mnth),kind=4)

             title = 'GISS LSM: Soil water; Bare layer 5 [m]'
             call array_gather( w_b5_mnth(:,:,i_mnth) )
             w_b5_mnth(:,:,i_mnth) = 
     &                w_b5_mnth(:,:,i_mnth)/n_count_mnth(i_mnth)
             if ( mype == 0 ) write(1007) title,
     &                      real(w_b5_mnth(:,:,i_mnth),kind=4)

             title = 'GISS LSM: Soil water; Bare layer 6 [m]'
             call array_gather( w_b6_mnth(:,:,i_mnth) )
             w_b6_mnth(:,:,i_mnth) = 
     &                w_b6_mnth(:,:,i_mnth)/n_count_mnth(i_mnth)
             if ( mype == 0 ) write(1008) title,
     &                      real(w_b6_mnth(:,:,i_mnth),kind=4)

             title = 'GISS LSM: Soil water; Veg layer 1 [m]'
             call array_gather( w_v1_mnth(:,:,i_mnth) )
             w_v1_mnth(:,:,i_mnth) = 
     &                w_v1_mnth(:,:,i_mnth)/n_count_mnth(i_mnth)
             if ( mype == 0 ) write(1009) title,
     &                      real(w_v1_mnth(:,:,i_mnth),kind=4)

             title = 'GISS LSM: Soil water; Veg layer 2 [m]'
             call array_gather( w_v2_mnth(:,:,i_mnth) )
             w_v2_mnth(:,:,i_mnth) = 
     &                w_v2_mnth(:,:,i_mnth)/n_count_mnth(i_mnth)
             if ( mype == 0 ) write(1010) title,
     &                      real(w_v2_mnth(:,:,i_mnth),kind=4)

             title = 'GISS LSM: Soil water; Veg layer 3 [m]'
             call array_gather( w_v3_mnth(:,:,i_mnth) )
             w_v3_mnth(:,:,i_mnth) = 
     &                w_v3_mnth(:,:,i_mnth)/n_count_mnth(i_mnth)
             if ( mype == 0 ) write(1011) title,
     &                      real(w_v3_mnth(:,:,i_mnth),kind=4)

             title = 'GISS LSM: Soil water; Veg layer 4 [m]'
             call array_gather( w_v4_mnth(:,:,i_mnth) )
             w_v4_mnth(:,:,i_mnth) = 
     &                w_v4_mnth(:,:,i_mnth)/n_count_mnth(i_mnth)
             if ( mype == 0 ) write(1012) title,
     &                      real(w_v4_mnth(:,:,i_mnth),kind=4)

             title = 'GISS LSM: Soil water; Veg layer 5 [m]'
             call array_gather( w_v5_mnth(:,:,i_mnth) )
             w_v5_mnth(:,:,i_mnth) = 
     &                w_v5_mnth(:,:,i_mnth)/n_count_mnth(i_mnth)
             if ( mype == 0 ) write(1013) title,
     &                      real(w_v5_mnth(:,:,i_mnth),kind=4)

             title = 'GISS LSM: Soil water; Veg layer 6 [m]'
             call array_gather( w_v6_mnth(:,:,i_mnth) )
             w_v6_mnth(:,:,i_mnth) = 
     &                w_v6_mnth(:,:,i_mnth)/n_count_mnth(i_mnth)
             if ( mype == 0 ) write(1014) title,
     &                      real(w_v6_mnth(:,:,i_mnth),kind=4)

             title = 'GISS LSM: Soil heat; Bare layer 1 [J/m2]'
             call array_gather( ht_b1_mnth(:,:,i_mnth) )
             ht_b1_mnth(:,:,i_mnth) = 
     &                ht_b1_mnth(:,:,i_mnth)/n_count_mnth(i_mnth)
             if ( mype == 0 ) write(1015) title,
     &                      real(ht_b1_mnth(:,:,i_mnth),kind=4)

             title = 'GISS LSM: Soil heat; Bare layer 2 [J/m2]'
             call array_gather( ht_b2_mnth(:,:,i_mnth) )
             ht_b2_mnth(:,:,i_mnth) = 
     &                ht_b2_mnth(:,:,i_mnth)/n_count_mnth(i_mnth)
             if ( mype == 0 ) write(1016) title,
     &                      real(ht_b2_mnth(:,:,i_mnth),kind=4)

             title = 'GISS LSM: Soil heat; Bare layer 3 [J/m2]'
             call array_gather( ht_b3_mnth(:,:,i_mnth) )
             ht_b3_mnth(:,:,i_mnth) = 
     &                ht_b3_mnth(:,:,i_mnth)/n_count_mnth(i_mnth)
             if ( mype == 0 ) write(1017) title,
     &                      real(ht_b3_mnth(:,:,i_mnth),kind=4)

             title = 'GISS LSM: Soil heat; Bare layer 4 [J/m2]'
             call array_gather( ht_b4_mnth(:,:,i_mnth) )
             ht_b4_mnth(:,:,i_mnth) = 
     &                ht_b4_mnth(:,:,i_mnth)/n_count_mnth(i_mnth)
             if ( mype == 0 ) write(1018) title,
     &                      real(ht_b4_mnth(:,:,i_mnth),kind=4)

             title = 'GISS LSM: Soil heat; Bare layer 5 [J/m2]'
             call array_gather( ht_b5_mnth(:,:,i_mnth) )
             ht_b5_mnth(:,:,i_mnth) = 
     &                ht_b5_mnth(:,:,i_mnth)/n_count_mnth(i_mnth)
             if ( mype == 0 ) write(1019) title,
     &                      real(ht_b5_mnth(:,:,i_mnth),kind=4)

             title = 'GISS LSM: Soil heat; Bare layer 6 [J/m2]'
             call array_gather( ht_b6_mnth(:,:,i_mnth) )
             ht_b6_mnth(:,:,i_mnth) = 
     &                ht_b6_mnth(:,:,i_mnth)/n_count_mnth(i_mnth)
             if ( mype == 0 ) write(1020) title,
     &                      real(ht_b6_mnth(:,:,i_mnth),kind=4)

             title = 'GISS LSM: Soil heat; Veg layer 1 [J/m2]'
             call array_gather( ht_v1_mnth(:,:,i_mnth) )
             ht_v1_mnth(:,:,i_mnth) = 
     &                ht_v1_mnth(:,:,i_mnth)/n_count_mnth(i_mnth)
             if ( mype == 0 ) write(1021) title,
     &                      real(ht_v1_mnth(:,:,i_mnth),kind=4)

             title = 'GISS LSM: Soil heat; Veg layer 2 [J/m2]'
             call array_gather( ht_v2_mnth(:,:,i_mnth) )
             ht_v2_mnth(:,:,i_mnth) = 
     &                ht_v2_mnth(:,:,i_mnth)/n_count_mnth(i_mnth)
             if ( mype == 0 ) write(1022) title,
     &                      real(ht_v2_mnth(:,:,i_mnth),kind=4)

             title = 'GISS LSM: Soil heat; Veg layer 3 [J/m2]'
             call array_gather( ht_v3_mnth(:,:,i_mnth) )
             ht_v3_mnth(:,:,i_mnth) = 
     &                ht_v3_mnth(:,:,i_mnth)/n_count_mnth(i_mnth)
             if ( mype == 0 ) write(1023) title,
     &                      real(ht_v3_mnth(:,:,i_mnth),kind=4)

             title = 'GISS LSM: Soil heat; Veg layer 4 [J/m2]'
             call array_gather( ht_v4_mnth(:,:,i_mnth) )
             ht_v4_mnth(:,:,i_mnth) = 
     &                ht_v4_mnth(:,:,i_mnth)/n_count_mnth(i_mnth)
             if ( mype == 0 ) write(1024) title,
     &                      real(ht_v4_mnth(:,:,i_mnth),kind=4)

             title = 'GISS LSM: Soil heat; Veg layer 5 [J/m2]'
             call array_gather( ht_v5_mnth(:,:,i_mnth) )
             ht_v5_mnth(:,:,i_mnth) = 
     &                ht_v5_mnth(:,:,i_mnth)/n_count_mnth(i_mnth)
             if ( mype == 0 ) write(1025) title,
     &                      real(ht_v5_mnth(:,:,i_mnth),kind=4)

             title = 'GISS LSM: Soil heat; Veg layer 6 [J/m2]'
             call array_gather( ht_v6_mnth(:,:,i_mnth) )
             ht_v6_mnth(:,:,i_mnth) = 
     &                ht_v6_mnth(:,:,i_mnth)/n_count_mnth(i_mnth)
             if ( mype == 0 ) write(1026) title,
     &                      real(ht_v6_mnth(:,:,i_mnth),kind=4)

             title = 'GISS LSM: Snow layer thickness; bare layer 1 [m]'
             call array_gather( dzsn_b1(:,:,i_mnth) )
             dzsn_b1(:,:,i_mnth) = 
     &                dzsn_b1(:,:,i_mnth)/n_count_mnth(i_mnth)
             if ( mype == 0 ) write(1027) title,
     &                      real(dzsn_b1(:,:,i_mnth),kind=4)

             title = 'GISS LSM:  Snow layer thickness; bare layer 2 [m]'
             call array_gather( dzsn_b2(:,:,i_mnth) )
             dzsn_b2(:,:,i_mnth) = 
     &                dzsn_b2(:,:,i_mnth)/n_count_mnth(i_mnth)
             if ( mype == 0 ) write(1028) title,
     &                      real(dzsn_b2(:,:,i_mnth),kind=4)

             title = 'GISS LSM: Snow layer thickness; bare layer 3 [m]'
             call array_gather( dzsn_b3(:,:,i_mnth) )
             dzsn_b3(:,:,i_mnth) = 
     &                dzsn_b3(:,:,i_mnth)/n_count_mnth(i_mnth)
             if ( mype == 0 ) write(1029) title,
     &                      real(dzsn_b3(:,:,i_mnth),kind=4)

             title = 'GISS LSM: Snow layer thickness; bare layer 1 [m]'
             call array_gather( dzsn_v1(:,:,i_mnth) )
             dzsn_v1(:,:,i_mnth) = 
     &                dzsn_v1(:,:,i_mnth)/n_count_mnth(i_mnth)
             if ( mype == 0 ) write(1030) title,
     &                      real(dzsn_v1(:,:,i_mnth),kind=4)

             title = 'GISS LSM: Snow layer thickness; bare layer 2 [m]'
             call array_gather( dzsn_v2(:,:,i_mnth) )
             dzsn_v2(:,:,i_mnth) = 
     &                dzsn_v2(:,:,i_mnth)/n_count_mnth(i_mnth)
             if ( mype == 0 ) write(1031) title,
     &                      real(dzsn_v2(:,:,i_mnth),kind=4)

             title = 'GISS LSM: Snow layer thickness; bare layer 3 [m]'
             call array_gather( dzsn_v3(:,:,i_mnth) )
             dzsn_v3(:,:,i_mnth) = 
     &                dzsn_v3(:,:,i_mnth)/n_count_mnth(i_mnth)
             if ( mype == 0 ) write(1032) title,
     &                      real(dzsn_v3(:,:,i_mnth),kind=4)

             title = 'GISS LSM: Snow water equival; bare layer 1 [m]'
             call array_gather( wsn_b1(:,:,i_mnth) )
             wsn_b1(:,:,i_mnth) = 
     &                wsn_b1(:,:,i_mnth)/n_count_mnth(i_mnth)
             if ( mype == 0 ) write(1033) title,
     &                      real(wsn_b1(:,:,i_mnth),kind=4)

             title = 'GISS LSM: Snow water equival; bare layer 2 [m]'
             call array_gather( wsn_b2(:,:,i_mnth) )
             wsn_b2(:,:,i_mnth) = 
     &                wsn_b2(:,:,i_mnth)/n_count_mnth(i_mnth)
             if ( mype == 0 ) write(1034) title,
     &                      real(wsn_b2(:,:,i_mnth),kind=4)

             title = 'GISS LSM: Snow water equival; bare layer 3 [m]'
             call array_gather( wsn_b3(:,:,i_mnth) )
             wsn_b3(:,:,i_mnth) = 
     &                wsn_b3(:,:,i_mnth)/n_count_mnth(i_mnth)
             if ( mype == 0 ) write(1035) title,
     &                      real(wsn_b3(:,:,i_mnth),kind=4)

             title = 'GISS LSM: Snow water equival; veg layer 1 [m]'
             call array_gather( wsn_v1(:,:,i_mnth) )
             wsn_v1(:,:,i_mnth) = 
     &                wsn_v1(:,:,i_mnth)/n_count_mnth(i_mnth)
             if ( mype == 0 ) write(1036) title,
     &                      real(wsn_v1(:,:,i_mnth),kind=4)

             title = 'GISS LSM: Snow water equival; veg layer 2 [m]'
             call array_gather( wsn_v2(:,:,i_mnth) )
             wsn_v2(:,:,i_mnth) = 
     &                wsn_v2(:,:,i_mnth)/n_count_mnth(i_mnth)
             if ( mype == 0 ) write(1037) title,
     &                      real(wsn_v2(:,:,i_mnth),kind=4)

             title = 'GISS LSM: Snow water equival; veg layer 3 [m]'
             call array_gather( wsn_v3(:,:,i_mnth) )
             wsn_v3(:,:,i_mnth) = 
     &                wsn_v3(:,:,i_mnth)/n_count_mnth(i_mnth)
             if ( mype == 0 ) write(1038) title,
     &                      real(wsn_v3(:,:,i_mnth),kind=4)

             title = 'GISS LSM: Snow layer heat; bare layer 1 [J/m2]'
             call array_gather( hsn_b1(:,:,i_mnth) )
             hsn_b1(:,:,i_mnth) = 
     &                hsn_b1(:,:,i_mnth)/n_count_mnth(i_mnth)
             if ( mype == 0 ) write(1039) title,
     &                      real(hsn_b1(:,:,i_mnth),kind=4)

             title = 'GISS LSM: Snow layer heat; bare layer 2 [J/m2]'
             call array_gather( hsn_b2(:,:,i_mnth) )
             hsn_b2(:,:,i_mnth) = 
     &                hsn_b2(:,:,i_mnth)/n_count_mnth(i_mnth)
             if ( mype == 0 ) write(1040) title,
     &                      real(hsn_b2(:,:,i_mnth),kind=4)

             title = 'GISS LSM:  Snow layer heat; bare layer 3 [J/m2]'
             call array_gather( hsn_b3(:,:,i_mnth) )
             hsn_b3(:,:,i_mnth) = 
     &                hsn_b3(:,:,i_mnth)/n_count_mnth(i_mnth)
             if ( mype == 0 ) write(1041) title,
     &                      real(hsn_b3(:,:,i_mnth),kind=4)

             title = 'GISS LSM:  Snow layer heat; veg layer 1 [J/m2]'
             call array_gather( hsn_v1(:,:,i_mnth) )
             hsn_v1(:,:,i_mnth) = 
     &                hsn_v1(:,:,i_mnth)/n_count_mnth(i_mnth)
             if ( mype == 0 ) write(1042) title,
     &                      real(hsn_v1(:,:,i_mnth),kind=4)

             title = 'GISS LSM:  Snow layer heat; veg layer 2 [J/m2]'
             call array_gather( hsn_v2(:,:,i_mnth) )
             hsn_v2(:,:,i_mnth) = 
     &                hsn_v2(:,:,i_mnth)/n_count_mnth(i_mnth)
             if ( mype == 0 ) write(1043) title,
     &                      real(hsn_v2(:,:,i_mnth),kind=4)

             title = 'GISS LSM:  Snow layer heat; veg layer 3 [J/m2]'
             call array_gather( hsn_v3(:,:,i_mnth) )
             hsn_v3(:,:,i_mnth) = 
     &                hsn_v3(:,:,i_mnth)/n_count_mnth(i_mnth)
             if ( mype == 0 ) write(1044) title,
     &                      real(hsn_v3(:,:,i_mnth),kind=4)

             title = 'GISS LSM: # of snow layers on bare land'
             call array_gather( nsn_b(:,:,i_mnth) )
             nsn_b(:,:,i_mnth) = 
     &                nsn_b(:,:,i_mnth)/n_count_mnth(i_mnth)
             if ( mype == 0 ) write(1045) title,
     &                      real(nsn_b(:,:,i_mnth),kind=4)

             title = 'GISS LSM: # of snow layers on vegetated land'
             call array_gather( nsn_v(:,:,i_mnth) )
             nsn_v(:,:,i_mnth) = 
     &                nsn_v(:,:,i_mnth)/n_count_mnth(i_mnth)
             if ( mype == 0 ) write(1046) title,
     &                      real(nsn_v(:,:,i_mnth),kind=4)

             title = 'GISS LSM: Snow-covered bare fraction [-]'
             call array_gather( frsnow_b(:,:,i_mnth) )
             frsnow_b(:,:,i_mnth) = 
     &                frsnow_b(:,:,i_mnth)/n_count_mnth(i_mnth)
             if ( mype == 0 ) write(1047) title,
     &                      real(frsnow_b(:,:,i_mnth),kind=4)

             title = 'GISS LSM:  Snow-covered vegetated fraction [-]'
             call array_gather( frsnow_v(:,:,i_mnth) )
             frsnow_v(:,:,i_mnth) = 
     &                frsnow_v(:,:,i_mnth)/n_count_mnth(i_mnth)
             if ( mype == 0 ) write(1048) title,
     &                      real(frsnow_v(:,:,i_mnth),kind=4)

             title = 'GISS LSM:Foliage surf. vapor mixing ratio [kg/kg]'
             call array_gather( Qf_mnth(:,:,i_mnth) )
             Qf_mnth(:,:,i_mnth) = 
     &                Qf_mnth(:,:,i_mnth)/n_count_mnth(i_mnth)
             if ( mype == 0 ) write(1049) title,
     &                      real(Qf_mnth(:,:,i_mnth),kind=4)

             title = 'GISS LSM:root zone betad [-]'
             call array_gather( abetad_mnth(:,:,i_mnth) )
             abetad_mnth(:,:,i_mnth) = 
     &                abetad_mnth(:,:,i_mnth)/n_count_mnth(i_mnth)
             if ( mype == 0 ) write(1050) title,
     &                      real(abetad_mnth(:,:,i_mnth),kind=4)

             title = 'GISS LSM:Live leaf carbon pool [kg/m2]'
             call array_gather( aClivepool_leaf_m(:,:,i_mnth) )
             aClivepool_leaf_m(:,:,i_mnth) = 
     &              aClivepool_leaf_m(:,:,i_mnth)/n_count_mnth(i_mnth)
             if ( mype == 0 ) write(1051) title,
     &                      real(aClivepool_leaf_m(:,:,i_mnth),kind=4)

             title = 'GISS LSM:Live froot carbon pool [kg/m2]'
             call array_gather( aClivepool_froot_m(:,:,i_mnth) )
             aClivepool_froot_m(:,:,i_mnth) = 
     &              aClivepool_froot_m(:,:,i_mnth)/n_count_mnth(i_mnth)
             if ( mype == 0 ) write(1052) title,
     &                    real(aClivepool_froot_m(:,:,i_mnth),kind=4)

             title = 'GISS LSM:Live wood carbon pool [kg/m2]'
             call array_gather( aClivepool_wood_m(:,:,i_mnth) )
             aClivepool_wood_m(:,:,i_mnth) = 
     &              aClivepool_wood_m(:,:,i_mnth)/n_count_mnth(i_mnth)
             if ( mype == 0 ) write(1053) title,
     &                    real(aClivepool_wood_m(:,:,i_mnth),kind=4)

             title = 'GISS LSM:Dead surface metabolic C pool [kg/m2]'
             call array_gather( aCdeadpool_surfmet_m(:,:,i_mnth) )
             aCdeadpool_surfmet_m(:,:,i_mnth) = 
     &            aCdeadpool_surfmet_m(:,:,i_mnth)/n_count_mnth(i_mnth)
             if ( mype == 0 ) write(1054) title,
     &                    real(aCdeadpool_surfmet_m(:,:,i_mnth),kind=4)

             title = 'GISS LSM:Dead surface structural C pool [kg/m2]'
             call array_gather(aCdeadpool_surfstr_m(:,:,i_mnth) )
             aCdeadpool_surfstr_m(:,:,i_mnth) = 
     &           aCdeadpool_surfstr_m(:,:,i_mnth)/n_count_mnth(i_mnth)
             if ( mype == 0 ) write(1055) title,
     &                    real(aCdeadpool_surfstr_m(:,:,i_mnth),kind=4)

             title = 'GISS LSM:Dead soil metabolic C pool [kg/m2]'
             call array_gather(aCdeadpool_soilmet_m(:,:,i_mnth) )
             aCdeadpool_soilmet_m(:,:,i_mnth) = 
     &           aCdeadpool_soilmet_m(:,:,i_mnth)/n_count_mnth(i_mnth)
             if ( mype == 0 ) write(1056) title,
     &                    real(aCdeadpool_soilmet_m(:,:,i_mnth),kind=4)

             title = 'GISS LSM:Dead soil structural C pool [kg/m2]'
             call array_gather(aCdeadpool_soilstr_m(:,:,i_mnth) )
             aCdeadpool_soilstr_m(:,:,i_mnth) = 
     &           aCdeadpool_soilstr_m(:,:,i_mnth)/n_count_mnth(i_mnth)
             if ( mype == 0 ) write(1057) title,
     &                    real(aCdeadpool_soilstr_m(:,:,i_mnth),kind=4)

             title = 'GISS LSM:Coarse woody debris C pool [kg/m2]'
             call array_gather( aCdeadpool_cwd_m(:,:,i_mnth) )
             aCdeadpool_cwd_m(:,:,i_mnth) = 
     &            aCdeadpool_cwd_m(:,:,i_mnth)/n_count_mnth(i_mnth)
             if ( mype == 0 ) write(1058) title,
     &                    real(aCdeadpool_cwd_m(:,:,i_mnth),kind=4)

             title = 'GISS LSM:Dead surface microbial C pool [kg/m2]'
             call array_gather(aCdeadpool_surfmic_m(:,:,i_mnth) )
             aCdeadpool_surfmic_m(:,:,i_mnth) = 
     &            aCdeadpool_surfmic_m(:,:,i_mnth)/n_count_mnth(i_mnth)
             if ( mype == 0 ) write(1059) title,
     &                    real(aCdeadpool_surfmic_m(:,:,i_mnth),kind=4)

             title = 'GISS LSM:Dead soil microbial C pool [kg/m2]'
             call array_gather( aCdeadpool_soilmic_m(:,:,i_mnth) )
             aCdeadpool_soilmic_m(:,:,i_mnth) = 
     &            aCdeadpool_soilmic_m(:,:,i_mnth)/n_count_mnth(i_mnth)
             if ( mype == 0 ) write(1060) title,
     &                    real(aCdeadpool_soilmic_m(:,:,i_mnth),kind=4)

             title = 'GISS LSM:Dead slow (~10 yr) C pool [kg/m2]'
             call array_gather(aCdeadpool_slow_m(:,:,i_mnth) )
             aCdeadpool_slow_m(:,:,i_mnth) = 
     &            aCdeadpool_slow_m(:,:,i_mnth)/n_count_mnth(i_mnth)
             if ( mype == 0 ) write(1061) title,
     &                    real(aCdeadpool_slow_m(:,:,i_mnth),kind=4)

             title = 'GISS LSM:Dead passive (~100 yr) C pool [kg/m2]'
             call array_gather(aCdeadpool_passive_m(:,:,i_mnth) )
             aCdeadpool_passive_m(:,:,i_mnth) = 
     &            aCdeadpool_passive_m(:,:,i_mnth)/n_count_mnth(i_mnth)
             if ( mype == 0 ) write(1062) title,
     &                    real(aCdeadpool_passive_m(:,:,i_mnth),kind=4)

             title = 'GISS LSM: Leaf Area Index [-]'
             call array_gather(alai_m(:,:,i_mnth) )
             alai_m(:,:,i_mnth)=alai_m(:,:,i_mnth)/n_count_mnth(i_mnth)
             if ( mype == 0 ) write(1063) title,
     &                    real(alai_m(:,:,i_mnth),kind=4)

             title = 'GISS LSM: Water stored on the canopy [m]'
             call array_gather(canopyH2O_m(:,:,i_mnth) )
             canopyH2O_m(:,:,i_mnth)=canopyH2O_m(:,:,i_mnth)
     &                               /n_count_mnth(i_mnth)
             if ( mype == 0 ) write(1064) title,
     &                    real(canopyH2O_m(:,:,i_mnth),kind=4)

             title = 'GISS LSM: Heat of the canopy [J/m2]'
             call array_gather(canopyheat_m(:,:,i_mnth) )
             canopyheat_m(:,:,i_mnth)=canopyheat_m(:,:,i_mnth)
     &                                /n_count_mnth(i_mnth)
             if ( mype == 0 ) write(1065) title,
     &                    real(canopyheat_m(:,:,i_mnth),kind=4)


          enddo
#else
!     Annual monthly diagnostic accumulators
          do i_mnth=1,12
              title = 'GISS LSM: Total evaporation'
              write(980) title,real(aevap_mnth(:,:,i_mnth),kind=4)
              title = 'GISS LSM: Canopy evaporation (evapw)'
              write(981) title,real(aevapw_mnth(:,:,i_mnth),kind=4)
              title = 'GISS LSM: Transpiration (evapd)'
              write(982) title,real(aevapd_mnth(:,:,i_mnth),kind=4)
              title = 'GISS LSM: Bare soil evaporation (evapb)'
              write(983) title,real(aevapb_mnth(:,:,i_mnth),kind=4)
              title = 'GISS LSM: Canopy interception'
              write(984) title,real(aintercep_mnth(:,:,i_mnth),kind=4)
              title = 'GISS LSM: Surface runoff (runs)'
              write(985) title,real(aruns_mnth(:,:,i_mnth),kind=4)
              title = 'GISS LSM: Subsurface flow out of soil column'
              write(986) title,real(arunu_mnth(:,:,i_mnth),kind=4)
              title = 'GISS LSM: Gross primary productivity'
              write(987) title,real(agpp_mnth(:,:,i_mnth),kind=4)
              title = 'GISS LSM: Autotrophic respiration'
              write(988) title,real(arauto_mnth(:,:,i_mnth),kind=4)
              title = 'GISS LSM: Soil respiration'
              write(989) title,real(asoilresp_mnth(:,:,i_mnth),kind=4)
              title = 'GISS LSM: Evap. from vegetated soil'
              write(990) title,real(aevapvg_mnth(:,:,i_mnth),kind=4)
              title = 'GISS LSM: Evap. from snow on vegetation'
              write(991) title,real(aevapvs_mnth(:,:,i_mnth),kind=4)
              title = 'GISS LSM: Evap. from snow on bare soil'
              write(992) title,real(aevapbs_mnth(:,:,i_mnth),kind=4)
              title = 'GISS LSM: Net shortwave radiation [W/m2]'
              asrht_mnth(:,:,i_mnth) = 
     &                asrht_mnth(:,:,i_mnth)/n_count_mnth(i_mnth)
              write(993) title,real(asrht_mnth(:,:,i_mnth),kind=4)
              title = 'GISS LSM: Net longwave radiation [W/m2]'
              atrht_mnth(:,:,i_mnth) = 
     &                atrht_mnth(:,:,i_mnth)/n_count_mnth(i_mnth)
              write(994) title,real(atrht_mnth(:,:,i_mnth),kind=4)
              title = 'GISS LSM: Grid-cell mean albedo'
              aalbedo_mnth(:,:,i_mnth) = 
     &                aalbedo_mnth(:,:,i_mnth)/n_count_mnth(i_mnth)
              write(995) title,real(aalbedo_mnth(:,:,i_mnth),kind=4)
              title = 'GISS LSM: Plant labile carbon [kgC/m2]'
              aclab_mnth(:,:,i_mnth) = 
     &                aclab_mnth(:,:,i_mnth)/n_count_mnth(i_mnth)
              write(996) title,real(aclab_mnth(:,:,i_mnth),kind=4)
              title = 'GISS LSM: Sum of soil carbon pools[kgC/m2]'
              asoilCpoolsum_mnth(:,:,i_mnth) = 
     &          asoilCpoolsum_mnth(:,:,i_mnth)/n_count_mnth(i_mnth)
              write(997) title,real(asoilCpoolsum_mnth(:,:,i_mnth)
     &                              ,kind=4)
              title = 'GISS LSM: Penman potential evap. [kg/m2/month]'
              write(998) title,real(aepp_mnth(:,:,i_mnth),kind=4)
              title = 'GISS LSM: Thermal heat from ground [J/m2/mnth]'
              write(999) title,real(atrg_mnth(:,:,i_mnth),kind=4)
              title = 'GISS LSM: Sensible heat from ground [J/m2/mnth]'
              write(1000) title,real(ashg_mnth(:,:,i_mnth),kind=4)
              title = 'GISS LSM: Latent heat flux [W/m2]'
              alhg_mnth(:,:,i_mnth) = 
     &                alhg_mnth(:,:,i_mnth)/n_count_mnth(i_mnth)
              write(1001) title,real(alhg_mnth(:,:,i_mnth),kind=4)
              title = 'GISS LSM: Pot. evap. from canopy [kg/m2/mnth]'
              write(1002) title,real(aepc_mnth(:,:,i_mnth),kind=4)
             title = 'GISS LSM: Soil water; Bare layer 1 [m]'
             w_b1_mnth(:,:,i_mnth) = 
     &                w_b1_mnth(:,:,i_mnth)/n_count_mnth(i_mnth)
             write(1003) title,real(w_b1_mnth(:,:,i_mnth),kind=4)
             title = 'GISS LSM: Soil water; Bare layer 2 [m]'
             w_b2_mnth(:,:,i_mnth) = 
     &                w_b2_mnth(:,:,i_mnth)/n_count_mnth(i_mnth)
             write(1004) title,real(w_b2_mnth(:,:,i_mnth),kind=4)
             title = 'GISS LSM: Soil water; Bare layer 3 [m]'
             w_b3_mnth(:,:,i_mnth) = 
     &                w_b3_mnth(:,:,i_mnth)/n_count_mnth(i_mnth)
             write(1005) title,real(w_b3_mnth(:,:,i_mnth),kind=4)
             title = 'GISS LSM: Soil water; Bare layer 4 [m]'
             w_b4_mnth(:,:,i_mnth) = 
     &                w_b4_mnth(:,:,i_mnth)/n_count_mnth(i_mnth)
             write(1006) title,real(w_b4_mnth(:,:,i_mnth),kind=4)
             title = 'GISS LSM: Soil water; Bare layer 5 [m]'
             w_b5_mnth(:,:,i_mnth) = 
     &                w_b5_mnth(:,:,i_mnth)/n_count_mnth(i_mnth)
             write(1007) title,real(w_b5_mnth(:,:,i_mnth),kind=4)
             title = 'GISS LSM: Soil water; Bare layer 6 [m]'
             w_b6_mnth(:,:,i_mnth) = 
     &                w_b6_mnth(:,:,i_mnth)/n_count_mnth(i_mnth)
             write(1008) title,real(w_b6_mnth(:,:,i_mnth),kind=4)
             title = 'GISS LSM: Soil water; Veg layer 1 [m]'
             w_v1_mnth(:,:,i_mnth) = 
     &                w_v1_mnth(:,:,i_mnth)/n_count_mnth(i_mnth)
             write(1009) title,real(w_v1_mnth(:,:,i_mnth),kind=4)
             title = 'GISS LSM: Soil water; Veg layer 2 [m]'
             w_v2_mnth(:,:,i_mnth) = 
     &                w_v2_mnth(:,:,i_mnth)/n_count_mnth(i_mnth)
             write(1010) title,real(w_v2_mnth(:,:,i_mnth),kind=4)
             title = 'GISS LSM: Soil water; Veg layer 3 [m]'
             w_v3_mnth(:,:,i_mnth) = 
     &                w_v3_mnth(:,:,i_mnth)/n_count_mnth(i_mnth)
             write(1011) title,real(w_v3_mnth(:,:,i_mnth),kind=4)
             title = 'GISS LSM: Soil water; Veg layer 4 [m]'
             w_v4_mnth(:,:,i_mnth) = 
     &                w_v4_mnth(:,:,i_mnth)/n_count_mnth(i_mnth)
             write(1012) title,real(w_v4_mnth(:,:,i_mnth),kind=4)
             title = 'GISS LSM: Soil water; Veg layer 5 [m]'
             w_v5_mnth(:,:,i_mnth) = 
     &                w_v5_mnth(:,:,i_mnth)/n_count_mnth(i_mnth)
             write(1013) title,real(w_v5_mnth(:,:,i_mnth),kind=4)
             title = 'GISS LSM: Soil water; Veg layer 6 [m]'
             w_v6_mnth(:,:,i_mnth) = 
     &                w_v6_mnth(:,:,i_mnth)/n_count_mnth(i_mnth)
             write(1014) title,real(w_v6_mnth(:,:,i_mnth),kind=4)
             title = 'GISS LSM: Soil heat; Bare layer 1 [J/m2]'
             ht_b1_mnth(:,:,i_mnth) = 
     &                ht_b1_mnth(:,:,i_mnth)/n_count_mnth(i_mnth)
             write(1015) title,real(ht_b1_mnth(:,:,i_mnth),kind=4)
             title = 'GISS LSM: Soil heat; Bare layer 2 [J/m2]'
             ht_b2_mnth(:,:,i_mnth) = 
     &                ht_b2_mnth(:,:,i_mnth)/n_count_mnth(i_mnth)
             write(1016) title,real(ht_b2_mnth(:,:,i_mnth),kind=4)
             title = 'GISS LSM: Soil heat; Bare layer 3 [J/m2]'
             ht_b3_mnth(:,:,i_mnth) = 
     &                ht_b3_mnth(:,:,i_mnth)/n_count_mnth(i_mnth)
             write(1017) title,real(ht_b3_mnth(:,:,i_mnth),kind=4)
             title = 'GISS LSM: Soil heat; Bare layer 4 [J/m2]'
             ht_b4_mnth(:,:,i_mnth) = 
     &                ht_b4_mnth(:,:,i_mnth)/n_count_mnth(i_mnth)
             write(1018) title,real(ht_b4_mnth(:,:,i_mnth),kind=4)
             title = 'GISS LSM: Soil heat; Bare layer 5 [J/m2]'
             ht_b5_mnth(:,:,i_mnth) = 
     &                ht_b5_mnth(:,:,i_mnth)/n_count_mnth(i_mnth)
             write(1019) title,real(ht_b5_mnth(:,:,i_mnth),kind=4)
             title = 'GISS LSM: Soil heat; Bare layer 6 [J/m2]'
             ht_b6_mnth(:,:,i_mnth) = 
     &                ht_b6_mnth(:,:,i_mnth)/n_count_mnth(i_mnth)
             write(1020) title,real(ht_b6_mnth(:,:,i_mnth),kind=4)
             title = 'GISS LSM: Soil heat; Veg layer 1 [J/m2]'
             ht_v1_mnth(:,:,i_mnth) = 
     &                ht_v1_mnth(:,:,i_mnth)/n_count_mnth(i_mnth)
             write(1021) title,real(ht_v1_mnth(:,:,i_mnth),kind=4)
             title = 'GISS LSM: Soil heat; Veg layer 2 [J/m2]'
             ht_v2_mnth(:,:,i_mnth) = 
     &                ht_v2_mnth(:,:,i_mnth)/n_count_mnth(i_mnth)
             write(1022) title,real(ht_v2_mnth(:,:,i_mnth),kind=4)
             title = 'GISS LSM: Soil heat; Veg layer 3 [J/m2]'
             ht_v3_mnth(:,:,i_mnth) = 
     &                ht_v3_mnth(:,:,i_mnth)/n_count_mnth(i_mnth)
             write(1023) title,real(ht_v3_mnth(:,:,i_mnth),kind=4)
             title = 'GISS LSM: Soil heat; Veg layer 4 [J/m2]'
             ht_v4_mnth(:,:,i_mnth) = 
     &                ht_v4_mnth(:,:,i_mnth)/n_count_mnth(i_mnth)
             write(1024) title,real(ht_v4_mnth(:,:,i_mnth),kind=4)
             title = 'GISS LSM: Soil heat; Veg layer 5 [J/m2]'
             ht_v5_mnth(:,:,i_mnth) = 
     &                ht_v5_mnth(:,:,i_mnth)/n_count_mnth(i_mnth)
             write(1025) title,real(ht_v5_mnth(:,:,i_mnth),kind=4)
             title = 'GISS LSM: Soil heat; Veg layer 6 [J/m2]'
             ht_v6_mnth(:,:,i_mnth) = 
     &                ht_v6_mnth(:,:,i_mnth)/n_count_mnth(i_mnth)
             write(1026) title,real(ht_v6_mnth(:,:,i_mnth),kind=4)
             title = 'GISS LSM: Snow layer thickness; bare layer 1 [m]'
             dzsn_b1(:,:,i_mnth) = 
     &                dzsn_b1(:,:,i_mnth)/n_count_mnth(i_mnth)
             write(1027) title,real(dzsn_b1(:,:,i_mnth),kind=4)
             title = 'GISS LSM:  Snow layer thickness; bare layer 2 [m]'
             dzsn_b2(:,:,i_mnth) = 
     &                dzsn_b2(:,:,i_mnth)/n_count_mnth(i_mnth)
             write(1028) title,real(dzsn_b2(:,:,i_mnth),kind=4)
             title = 'GISS LSM: Snow layer thickness; bare layer 3 [m]'
             dzsn_b3(:,:,i_mnth) = 
     &                dzsn_b3(:,:,i_mnth)/n_count_mnth(i_mnth)
             write(1029) title,real(dzsn_b3(:,:,i_mnth),kind=4)
             title = 'GISS LSM: Snow layer thickness; bare layer 1 [m]'
             dzsn_v1(:,:,i_mnth) = 
     &                dzsn_v1(:,:,i_mnth)/n_count_mnth(i_mnth)
             write(1030) title,real(dzsn_v1(:,:,i_mnth),kind=4)
             title = 'GISS LSM: Snow layer thickness; bare layer 2 [m]'
             dzsn_v2(:,:,i_mnth) = 
     &                dzsn_v2(:,:,i_mnth)/n_count_mnth(i_mnth)
             write(1031) title,real(dzsn_v2(:,:,i_mnth),kind=4)
             title = 'GISS LSM: Snow layer thickness; bare layer 3 [m]'
             dzsn_v3(:,:,i_mnth) = 
     &                dzsn_v3(:,:,i_mnth)/n_count_mnth(i_mnth)
             write(1032) title,real(dzsn_v3(:,:,i_mnth),kind=4)
             title = 'GISS LSM: Snow water equival; bare layer 1 [m]'
             wsn_b1(:,:,i_mnth) = 
     &                wsn_b1(:,:,i_mnth)/n_count_mnth(i_mnth)
             write(1033) title,real(wsn_b1(:,:,i_mnth),kind=4)
             title = 'GISS LSM: Snow water equival; bare layer 2 [m]'
             wsn_b2(:,:,i_mnth) = 
     &                wsn_b2(:,:,i_mnth)/n_count_mnth(i_mnth)
             write(1034) title,real(wsn_b2(:,:,i_mnth),kind=4)
             title = 'GISS LSM: Snow water equival; bare layer 3 [m]'
             wsn_b3(:,:,i_mnth) = 
     &                wsn_b3(:,:,i_mnth)/n_count_mnth(i_mnth)
             write(1035) title,real(wsn_b3(:,:,i_mnth),kind=4)
             title = 'GISS LSM: Snow water equival; veg layer 1 [m]'
             wsn_v1(:,:,i_mnth) = 
     &                wsn_v1(:,:,i_mnth)/n_count_mnth(i_mnth)
             write(1036) title,real(wsn_v1(:,:,i_mnth),kind=4)
             title = 'GISS LSM: Snow water equival; veg layer 2 [m]'
             wsn_v2(:,:,i_mnth) = 
     &                wsn_v2(:,:,i_mnth)/n_count_mnth(i_mnth)
             write(1037) title,real(wsn_v2(:,:,i_mnth),kind=4)
             title = 'GISS LSM: Snow water equival; veg layer 3 [m]'
             wsn_v3(:,:,i_mnth) = 
     &                wsn_v3(:,:,i_mnth)/n_count_mnth(i_mnth)
             write(1038) title,real(wsn_v3(:,:,i_mnth),kind=4)
             title = 'GISS LSM: Snow layer heat; bare layer 1 [J/m2]'
             hsn_b1(:,:,i_mnth) = 
     &                hsn_b1(:,:,i_mnth)/n_count_mnth(i_mnth)
             write(1039) title,real(hsn_b1(:,:,i_mnth),kind=4)
             title = 'GISS LSM: Snow layer heat; bare layer 2 [J/m2]'
             hsn_b2(:,:,i_mnth) = 
     &                hsn_b2(:,:,i_mnth)/n_count_mnth(i_mnth)
             write(1040) title,real(hsn_b2(:,:,i_mnth),kind=4)
             title = 'GISS LSM:  Snow layer heat; bare layer 3 [J/m2]'
             hsn_b3(:,:,i_mnth) = 
     &                hsn_b3(:,:,i_mnth)/n_count_mnth(i_mnth)
             write(1041) title,real(hsn_b3(:,:,i_mnth),kind=4)
             title = 'GISS LSM:  Snow layer heat; veg layer 1 [J/m2]'
             hsn_v1(:,:,i_mnth) = 
     &                hsn_v1(:,:,i_mnth)/n_count_mnth(i_mnth)
             write(1042) title,real(hsn_v1(:,:,i_mnth),kind=4)
             title = 'GISS LSM:  Snow layer heat; veg layer 2 [J/m2]'
             hsn_v2(:,:,i_mnth) = 
     &                hsn_v2(:,:,i_mnth)/n_count_mnth(i_mnth)
             write(1043) title,real(hsn_v2(:,:,i_mnth),kind=4)
             title = 'GISS LSM:  Snow layer heat; veg layer 3 [J/m2]'
             hsn_v3(:,:,i_mnth) = 
     &                hsn_v3(:,:,i_mnth)/n_count_mnth(i_mnth)
             write(1044) title,real(hsn_v3(:,:,i_mnth),kind=4)
             title = 'GISS LSM: # of snow layers on bare land'
             nsn_b(:,:,i_mnth) = 
     &                nsn_b(:,:,i_mnth)/n_count_mnth(i_mnth)
             write(1045) title,real(nsn_b(:,:,i_mnth),kind=4)
             title = 'GISS LSM: # of snow layers on vegetated land'
             nsn_v(:,:,i_mnth) = 
     &                nsn_v(:,:,i_mnth)/n_count_mnth(i_mnth)
             write(1046) title,real(nsn_v(:,:,i_mnth),kind=4)
             title = 'GISS LSM: Snow-covered bare fraction [-]'
             frsnow_b(:,:,i_mnth) = 
     &                frsnow_b(:,:,i_mnth)/n_count_mnth(i_mnth)
             write(1047) title,real(frsnow_b(:,:,i_mnth),kind=4)
             title = 'GISS LSM:  Snow-covered vegetated fraction [-]'
             frsnow_v(:,:,i_mnth) = 
     &                frsnow_v(:,:,i_mnth)/n_count_mnth(i_mnth)
             write(1048) title,real(frsnow_v(:,:,i_mnth),kind=4)
             title = 'GISS LSM:Foliage surf. vapor mixing ratio [kg/kg]'
             Qf_mnth(:,:,i_mnth) = 
     &                Qf_mnth(:,:,i_mnth)/n_count_mnth(i_mnth)
             write(1049) title,real(Qf_mnth(:,:,i_mnth),kind=4)
             title = 'GISS LSM:root zone betad [-]'
             abetad_mnth(:,:,i_mnth) = 
     &                abetad_mnth(:,:,i_mnth)/n_count_mnth(i_mnth)
             write(1050) title,real(abetad_mnth(:,:,i_mnth),kind=4)
             title = 'GISS LSM:Live leaf carbon pool [kg/m2]'
             aClivepool_leaf_m(:,:,i_mnth) = 
     &              aClivepool_leaf_m(:,:,i_mnth)/n_count_mnth(i_mnth)
             write(1051) title,
     &              real(aClivepool_leaf_m(:,:,i_mnth),kind=4)
             title = 'GISS LSM:Live froot carbon pool [kg/m2]'
             aClivepool_froot_m(:,:,i_mnth) = 
     &              aClivepool_froot_m(:,:,i_mnth)/n_count_mnth(i_mnth)
             write(1052) title,
     &                    real(aClivepool_froot_m(:,:,i_mnth),kind=4)
             title = 'GISS LSM:Live wood carbon pool [kg/m2]'
             aClivepool_wood_m(:,:,i_mnth) = 
     &              aClivepool_wood_m(:,:,i_mnth)/n_count_mnth(i_mnth)
             write(1053) title,
     &                    real(aClivepool_wood_m(:,:,i_mnth),kind=4)
             title = 'GISS LSM:Dead surface metabolic C pool [kg/m2]'
             aCdeadpool_surfmet_m(:,:,i_mnth) = 
     &            aCdeadpool_surfmet_m(:,:,i_mnth)/n_count_mnth(i_mnth)
             write(1054) title,
     &                    real(aCdeadpool_surfmet_m(:,:,i_mnth),kind=4)
             title = 'GISS LSM:Dead surface structural C pool [kg/m2]'
             aCdeadpool_surfstr_m(:,:,i_mnth) = 
     &           aCdeadpool_surfstr_m(:,:,i_mnth)/n_count_mnth(i_mnth)
             write(1055) title,
     &                    real(aCdeadpool_surfstr_m(:,:,i_mnth),kind=4)
             title = 'GISS LSM:Dead soil metabolic C pool [kg/m2]'
             aCdeadpool_soilmet_m(:,:,i_mnth) = 
     &           aCdeadpool_soilmet_m(:,:,i_mnth)/n_count_mnth(i_mnth)
             write(1056) title,
     &                    real(aCdeadpool_soilmet_m(:,:,i_mnth),kind=4)
             title = 'GISS LSM:Dead soil structural C pool [kg/m2]'
             aCdeadpool_soilstr_m(:,:,i_mnth) = 
     &           aCdeadpool_soilstr_m(:,:,i_mnth)/n_count_mnth(i_mnth)
             write(1057) title,
     &                    real(aCdeadpool_soilstr_m(:,:,i_mnth),kind=4)
             title = 'GISS LSM:Coarse woody debris C pool [kg/m2]'
             aCdeadpool_cwd_m(:,:,i_mnth) = 
     &            aCdeadpool_cwd_m(:,:,i_mnth)/n_count_mnth(i_mnth)
             write(1058) title,
     &                    real(aCdeadpool_cwd_m(:,:,i_mnth),kind=4)
             title = 'GISS LSM:Dead surface microbial C pool [kg/m2]'
             aCdeadpool_surfmic_m(:,:,i_mnth) = 
     &            aCdeadpool_surfmic_m(:,:,i_mnth)/n_count_mnth(i_mnth)
             write(1059) title,
     &                    real(aCdeadpool_surfmic_m(:,:,i_mnth),kind=4)
             title = 'GISS LSM:Dead soil microbial C pool [kg/m2]'
             aCdeadpool_soilmic_m(:,:,i_mnth) = 
     &            aCdeadpool_soilmic_m(:,:,i_mnth)/n_count_mnth(i_mnth)
             write(1060) title,
     &                    real(aCdeadpool_soilmic_m(:,:,i_mnth),kind=4)
             title = 'GISS LSM:Dead slow (~10 yr) C pool [kg/m2]'
             aCdeadpool_slow_m(:,:,i_mnth) = 
     &            aCdeadpool_slow_m(:,:,i_mnth)/n_count_mnth(i_mnth)
             write(1061) title,
     &                    real(aCdeadpool_slow_m(:,:,i_mnth),kind=4)
             title = 'GISS LSM:Dead passive (~100 yr) C pool [kg/m2]'
             aCdeadpool_passive_m(:,:,i_mnth) = 
     &            aCdeadpool_passive_m(:,:,i_mnth)/n_count_mnth(i_mnth)
             write(1062) title,
     &                    real(aCdeadpool_passive_m(:,:,i_mnth),kind=4)
             title = 'GISS LSM:Leaf Area Index [-]'
             alai_m(:,:,i_mnth)=alai_m(:,:,i_mnth)/n_count_mnth(i_mnth)
             write(1063) title,real(alai_m(:,:,i_mnth),kind=4)

             title = 'GISS LSM: Water stored on the canopy [m]'
             canopyH2O_m(:,:,i_mnth)=canopyH2O_m(:,:,i_mnth)
     &                               /n_count_mnth(i_mnth)
             write(1064) title,real(canopyH2O_m(:,:,i_mnth),kind=4)

             title = 'GISS LSM: Heat of the canopy [J/m2]'
             canopyheat_m(:,:,i_mnth)=canopyheat_m(:,:,i_mnth)
     &                                /n_count_mnth(i_mnth)
             write(1065) title,real(canopyheat_m(:,:,i_mnth),kind=4)

          enddo
#endif
      end subroutine print_diags

#endif 

!-----------------------------------------------------------------------
#ifdef PRINT_DIAGS  /*  Monthly diagnostic accumulators */

      subroutine zero_diags(
     &                   aevap_mnth,aevapw_mnth,aevapd_mnth,aevapb_mnth
     &                 , aintercep_mnth,aruns_mnth,arunu_mnth
     &                 , agpp_mnth,arauto_mnth,asoilresp_mnth
     &                 , aevapvg_mnth,aevapvs_mnth,aevapbs_mnth
     &                 , asrht_mnth,atrht_mnth,aalbedo_mnth
     &                 , aclab_mnth,asoilCpoolsum_mnth,aepp_mnth
     &                 , atrg_mnth,ashg_mnth,alhg_mnth,aepc_mnth
     &                 , n_count_mnth
     &                 , w_b1_mnth,w_b2_mnth,w_b3_mnth
     &                 , w_b4_mnth,w_b5_mnth,w_b6_mnth
     &                 , w_v1_mnth,w_v2_mnth,w_v3_mnth
     &                 , w_v4_mnth,w_v5_mnth,w_v6_mnth
     &                 , ht_b1_mnth,ht_b2_mnth,ht_b3_mnth
     &                 , ht_b4_mnth,ht_b5_mnth,ht_b6_mnth
     &                 , ht_v1_mnth,ht_v2_mnth,ht_v3_mnth
     &                 , ht_v4_mnth,ht_v5_mnth,ht_v6_mnth
     &                 , dzsn_b1,dzsn_b2,dzsn_b3,dzsn_v1,dzsn_v2,dzsn_v3
     &                 , wsn_b1 ,wsn_b2 ,wsn_b3 ,wsn_v1 ,wsn_v2 ,wsn_v3
     &                 , hsn_b1 ,hsn_b2 ,hsn_b3 ,hsn_v1 ,hsn_v2 ,hsn_v3
     &                 , nsn_b,nsn_v,frsnow_b,frsnow_v,Qf_mnth
     &                 , abetad_mnth
     &  , aClivepool_leaf_m,aClivepool_froot_m,aClivepool_wood_m
     &  , aCdeadpool_surfmet_m,aCdeadpool_surfstr_m,aCdeadpool_soilmet_m
     &  , aCdeadpool_soilstr_m,aCdeadpool_cwd_m,aCdeadpool_surfmic_m
     &  , aCdeadpool_soilmic_m,aCdeadpool_slow_m,aCdeadpool_passive_m
     &  , alai_m,canopyH2O_m,canopyheat_m
     &    )     

      real*8, dimension(im,jm,12) ::
     &                    aevap_mnth,aevapw_mnth,aevapd_mnth,aevapb_mnth
     &                  , aintercep_mnth,aruns_mnth,arunu_mnth
     &                  , agpp_mnth,arauto_mnth,asoilresp_mnth
     &                  , aevapvg_mnth,aevapvs_mnth,aevapbs_mnth
     &                  , asrht_mnth,atrht_mnth,aalbedo_mnth
     &                  , aclab_mnth,asoilCpoolsum_mnth,aepp_mnth
     &                  , atrg_mnth,ashg_mnth,alhg_mnth,aepc_mnth
     &                  , w_b1_mnth,w_b2_mnth,w_b3_mnth
     &                  , w_b4_mnth,w_b5_mnth,w_b6_mnth
     &                  , w_v1_mnth,w_v2_mnth,w_v3_mnth
     &                  , w_v4_mnth,w_v5_mnth,w_v6_mnth
     &                  , ht_b1_mnth,ht_b2_mnth,ht_b3_mnth
     &                  , ht_b4_mnth,ht_b5_mnth,ht_b6_mnth
     &                  , ht_v1_mnth,ht_v2_mnth,ht_v3_mnth
     &                  , ht_v4_mnth,ht_v5_mnth,ht_v6_mnth
     &                 , dzsn_b1,dzsn_b2,dzsn_b3,dzsn_v1,dzsn_v2,dzsn_v3
     &                 , wsn_b1 ,wsn_b2 ,wsn_b3 ,wsn_v1 ,wsn_v2 ,wsn_v3
     &                 , hsn_b1 ,hsn_b2 ,hsn_b3 ,hsn_v1 ,hsn_v2 ,hsn_v3
     &                 , nsn_b,nsn_v,frsnow_b,frsnow_v,Qf_mnth
     &                 , abetad_mnth
     &  , aClivepool_leaf_m,aClivepool_froot_m,aClivepool_wood_m
     &  , aCdeadpool_surfmet_m,aCdeadpool_surfstr_m,aCdeadpool_soilmet_m
     &  , aCdeadpool_soilstr_m,aCdeadpool_cwd_m,aCdeadpool_surfmic_m
     &  , aCdeadpool_soilmic_m,aCdeadpool_slow_m,aCdeadpool_passive_m
     &  , alai_m,canopyH2O_m,canopyheat_m
      real*8,dimension(12),intent(inout) :: n_count_mnth

!         Reset annual monthly diagnostic/prognostic accumulators
          aevap_mnth(:,:,:) = 0.d0
          aevapw_mnth(:,:,:) = 0.d0
          aevapd_mnth(:,:,:) = 0.d0
          aevapb_mnth(:,:,:) = 0.d0
          aintercep_mnth(:,:,:) = 0.d0
          aruns_mnth(:,:,:) = 0.d0
          arunu_mnth(:,:,:) = 0.d0
          agpp_mnth(:,:,:) = 0.d0
          arauto_mnth(:,:,:) = 0.d0
          asoilresp_mnth(:,:,:) = 0.d0
          aevapvg_mnth(:,:,:) = 0.d0
          aevapvs_mnth(:,:,:) = 0.d0
          aevapbs_mnth(:,:,:) = 0.d0
          asrht_mnth(:,:,:) = 0.d0
          atrht_mnth(:,:,:) = 0.d0
          aalbedo_mnth(:,:,:) = 0.d0
          aclab_mnth(:,:,:) = 0.d0
          asoilCpoolsum_mnth(:,:,:) = 0.d0
          aepp_mnth(:,:,:) = 0.d0
          atrg_mnth(:,:,:) = 0.d0
          ashg_mnth(:,:,:) = 0.d0
          alhg_mnth(:,:,:) = 0.d0
          aepc_mnth(:,:,:) = 0.d0
          w_b1_mnth(:,:,:)=0.d0
          w_b2_mnth(:,:,:)=0.d0
          w_b3_mnth(:,:,:)=0.d0
          w_b4_mnth(:,:,:)=0.d0
          w_b5_mnth(:,:,:)=0.d0
          w_b6_mnth(:,:,:)=0.d0
          w_v1_mnth(:,:,:)=0.d0
          w_v2_mnth(:,:,:)=0.d0
          w_v3_mnth(:,:,:)=0.d0
          w_v4_mnth(:,:,:)=0.d0
          w_v5_mnth(:,:,:)=0.d0
          w_v6_mnth(:,:,:)=0.d0
          ht_b1_mnth(:,:,:)=0.d0
          ht_b2_mnth(:,:,:)=0.d0
          ht_b3_mnth(:,:,:)=0.d0
          ht_b4_mnth(:,:,:)=0.d0
          ht_b5_mnth(:,:,:)=0.d0
          ht_b6_mnth(:,:,:)=0.d0
          ht_v1_mnth(:,:,:)=0.d0
          ht_v2_mnth(:,:,:)=0.d0
          ht_v3_mnth(:,:,:)=0.d0
          ht_v4_mnth(:,:,:)=0.d0
          ht_v5_mnth(:,:,:)=0.d0
          ht_v6_mnth(:,:,:)=0.d0
          dzsn_b1(:,:,:)=0.d0
          dzsn_b2(:,:,:)=0.d0
          dzsn_b3(:,:,:)=0.d0
          dzsn_v1(:,:,:)=0.d0
          dzsn_v2(:,:,:)=0.d0
          dzsn_v3(:,:,:)=0.d0
          wsn_b1(:,:,:)=0.d0
          wsn_b2(:,:,:)=0.d0
          wsn_b3(:,:,:)=0.d0
          wsn_v1(:,:,:)=0.d0
          wsn_v2(:,:,:)=0.d0
          wsn_v3(:,:,:)=0.d0
          hsn_b1(:,:,:)=0.d0
          hsn_b2(:,:,:)=0.d0
          hsn_b3(:,:,:)=0.d0
          hsn_v1(:,:,:)=0.d0
          hsn_v2(:,:,:)=0.d0
          hsn_v3(:,:,:)=0.d0
          nsn_b(:,:,:)=0.d0
          nsn_v (:,:,:)=0.d0
          frsnow_b(:,:,:)=0.d0
          frsnow_v(:,:,:)=0.d0
          Qf_mnth(:,:,:)=0.d0
          abetad_mnth(:,:,:)=0.d0
          aClivepool_leaf_m(:,:,:)=0.d0
          aClivepool_froot_m(:,:,:)=0.d0
          aClivepool_wood_m(:,:,:)=0.d0
          aCdeadpool_surfmet_m(:,:,:)=0.d0
          aCdeadpool_surfstr_m(:,:,:)=0.d0
          aCdeadpool_soilmet_m(:,:,:)=0.d0
          aCdeadpool_soilstr_m(:,:,:)=0.d0
          aCdeadpool_cwd_m(:,:,:)=0.d0
          aCdeadpool_surfmic_m(:,:,:)=0.d0
          aCdeadpool_soilmic_m(:,:,:)=0.d0
          aCdeadpool_slow_m(:,:,:)=0.d0
          aCdeadpool_passive_m(:,:,:)=0.d0
          alai_m(:,:,:)=0.d0
          canopyH2O_m(:,:,:)=0.d0
          canopyheat_m(:,:,:)=0.d0
          n_count_mnth(:) = 0.d0 !Reset counter

      end subroutine zero_diags

#endif
!-----------------------------------------------------------------------

      end module lsm

!-----------------------------------------------------------------------

cddd      MODULE TRIDIAG_MOD
cddd!@sum TRIDIAG_MOD contains subroutine TRIDIAG
cddd      PRIVATE
cddd      PUBLIC TRIDIAG
cddd
cddd      contains
cddd
cddd      SUBROUTINE TRIDIAG(A,B,C,R,U,N)
cddd!@sum  TRIDIAG  solves a tridiagonal matrix equation (A,B,C)U=R
cddd!@auth Numerical Recipes
cddd!@ver  1.0
cddd      IMPLICIT NONE
cdddc      INTEGER, PARAMETER :: NMAX = 8000  !@var NMAX workspace
cddd      INTEGER, INTENT(IN):: N         !@var N    dimension of arrays
cddd      REAL*8, INTENT(IN) :: A(N)   !@var A    coefficients of u_i-1
cddd      REAL*8, INTENT(IN) :: B(N)   !@var B    coefficients of u_i
cddd      REAL*8, INTENT(IN) :: C(N)   !@var C    coefficients of u_i+1
cddd      REAL*8, INTENT(IN) :: R(N)   !@var R    RHS vector
cddd      REAL*8, INTENT(OUT):: U(N)   !@var U    solution vector
cddd      REAL*8 :: BET                !@var BET  work variable
cddd      REAL*8 :: GAM(Size(A))       !@var GAM  work array
cddd      INTEGER :: J                 !@var J    loop variable
cddd
cdddc      IF ( N > NMAX )
cdddc     &     call stop_model("TRIDIAG: N > NMAX, increase NMAX",255)
cddd      BET=B(1)
cddd      IF (BET.eq.0) call stop_model("TRIDIAG: DENOMINATOR = ZERO",255)
cddd      U(1)=R(1)/BET
cddd      DO J=2,N
cddd        GAM(J)=C(J-1)/BET
cddd        BET=B(J)-A(J)*GAM(J)
cddd        IF (BET.eq.0) call stop_model("TRIDIAG: DENOMINATOR = ZERO",255)
cddd        U(J)=(R(J)-A(J)*U(J-1))/BET
cddd      END DO
cddd      DO J=N-1,1,-1
cddd        U(J)=U(J)-GAM(J+1)*U(J+1)
cddd      END DO
cddd      RETURN
cddd      END SUBROUTINE TRIDIAG
cddd
cddd      end MODULE TRIDIAG_MOD

!-----------------------------------------------------------------------

      subroutine LSM_standalone
      use parser
      use param
      use filemanager, only : openunit,closeunit
      use lsm
      use ent_mod
#ifdef USE_GSWP_FORCINGS
      use drv_gswp_force, only : init_forcings, get_month
#endif
      use domain_decomp, only : app_init, mype
      implicit none
#ifdef USE_ESMF
#include "mpi_defs.h"
#include "mpif.h"
#endif

      type (t_lsm_bc) :: bc
      type(t_lsm_state) :: lsm_state
!@var Starting simulation timestep number [-]
      integer, parameter :: itime0=1
!@var Simulation timestep number of size dtsec [-]
      integer :: itime
!@var Timestep number of size 3 hours [-]
      integer :: itime_3hr
!@var Timestep size of forcings data [seconds]
      real*8  :: dtsec
!@var Simulation end time [number of timesteps of size dtsec]
      integer :: itime1
!@var Counter for time of year [seconds]
      integer :: tyr_sec
!@var Gregorian day of year (julian day of year for modelE met forcings)
      integer :: jday
!@var Gregorian year (julian year for modelE met forcings)
      integer :: jyear
!&var Current month number & month counter
      integer :: n_month, i_mnth
!@var Counter for spinup years
      integer :: i_spin 
      integer :: spin_start_yr, spin_end_yr, num_times_spin

      integer, save :: jday_old = -32768
      real*8  :: time
      integer :: iu_LAI, iu_vht
      
#ifdef PRINT_DIAGS  /*  Monthly diagnostic accumulators */
      character*80 :: title       ! title of diagnostic to print to file
      real*8, dimension(im,jm,12) ::
     &                    aevap_mnth,aevapw_mnth,aevapd_mnth,aevapb_mnth
     &                  , aintercep_mnth,aruns_mnth,arunu_mnth
     &                  , agpp_mnth,arauto_mnth,asoilresp_mnth
     &                  , aevapvg_mnth,aevapvs_mnth,aevapbs_mnth
     &                  , asrht_mnth,atrht_mnth,aalbedo_mnth
     &                  , aclab_mnth,asoilCpoolsum_mnth,aepp_mnth
     &                  , atrg_mnth,ashg_mnth,alhg_mnth,aepc_mnth
     &                  , w_b1_mnth,w_b2_mnth,w_b3_mnth
     &                  , w_b4_mnth,w_b5_mnth,w_b6_mnth
     &                  , w_v1_mnth,w_v2_mnth,w_v3_mnth
     &                  , w_v4_mnth,w_v5_mnth,w_v6_mnth
     &                  , ht_b1_mnth,ht_b2_mnth,ht_b3_mnth
     &                  , ht_b4_mnth,ht_b5_mnth,ht_b6_mnth
     &                  , ht_v1_mnth,ht_v2_mnth,ht_v3_mnth
     &                  , ht_v4_mnth,ht_v5_mnth,ht_v6_mnth
     &                 , dzsn_b1,dzsn_b2,dzsn_b3,dzsn_v1,dzsn_v2,dzsn_v3
     &                 , wsn_b1 ,wsn_b2 ,wsn_b3 ,wsn_v1 ,wsn_v2 ,wsn_v3
     &                 , hsn_b1 ,hsn_b2 ,hsn_b3 ,hsn_v1 ,hsn_v2 ,hsn_v3
     &                 , nsn_b,nsn_v,frsnow_b,frsnow_v,Qf_mnth
     &                 , abetad_mnth
     &  , aClivepool_leaf_m,aClivepool_froot_m,aClivepool_wood_m
     &  , aCdeadpool_surfmet_m,aCdeadpool_surfstr_m,aCdeadpool_soilmet_m
     &  , aCdeadpool_soilstr_m,aCdeadpool_cwd_m,aCdeadpool_surfmic_m
     &  , aCdeadpool_soilmic_m,aCdeadpool_slow_m,aCdeadpool_passive_m
     &  , alai_m,canopyH2O_m,canopyheat_m
      real*8,dimension(12) :: n_count_mnth
#endif

#ifdef USE_GSWP_FORCINGS
!@var Number of GSWP2 variables
      integer, parameter :: varnums = 9
      logical :: end_of_input_flag
      real*8 :: num_input_steps
#ifdef ECOSYSTEM_SCALE
      real*8, save :: tbcs_ij, tcanopy_ij
      integer,save :: iu_vector
#else
      real*8,dimension(im,jm), save :: tbcs_ij, tcanopy_ij
      integer,dimension(varnums),save :: iu_vector
#endif

#endif
      !integer fd, offset, idim, jdim
      integer :: j0, j1, ierr
      logical :: end_of_day_flag

!     Default values of run settings
#ifdef USE_GSWP_FORCINGS     
      jyear = 1983
      jday = 1
      tyr_sec = (jday-1)*24*3600
      i_spin = 1
      spin_start_yr = 1983
      spin_end_yr = 1985
      num_times_spin = 2
#else
!     Default of GCM setup
      jday = 335
      jyear = 1949
      jday = 1
      tyr_sec = (jday-1)*24*3600
#endif

!     Read input file with run settings
      call read_input_parameters(itime1,jyear,jday,tyr_sec,dtsec
     &       ,spin_start_yr, spin_end_yr, num_times_spin)

      if (force_VEG) then
        call openunit(trim(LAI_file),iu_LAI,.false.,.true.)
        call openunit(trim(vheight),iu_vht,.false.,.true.)
      else
        iu_LAI = 0
        iu_vht = 0
      end if

      call app_init(jm,j0,j1)
      !print *,"domain decomposition: j0,j1 = ", j0,j1

!     Initialize GHY and ENT here ...
      call lsm_initialize( dtsec,n_month,tbcs_ij,tcanopy_ij
#ifdef PRINT_DIAGS
     &                 , aevap_mnth,aevapw_mnth,aevapd_mnth,aevapb_mnth
     &                 , aintercep_mnth,aruns_mnth,arunu_mnth
     &                 , agpp_mnth,arauto_mnth,asoilresp_mnth
     &                 , aevapvg_mnth,aevapvs_mnth,aevapbs_mnth
     &                 , asrht_mnth,atrht_mnth,aalbedo_mnth
     &                 , aclab_mnth,asoilCpoolsum_mnth,aepp_mnth
     &                 , atrg_mnth,ashg_mnth,alhg_mnth,aepc_mnth
     &                 , n_count_mnth
     &                 , w_b1_mnth,w_b2_mnth,w_b3_mnth
     &                 , w_b4_mnth,w_b5_mnth,w_b6_mnth
     &                 , w_v1_mnth,w_v2_mnth,w_v3_mnth
     &                 , w_v4_mnth,w_v5_mnth,w_v6_mnth
     &                 , ht_b1_mnth,ht_b2_mnth,ht_b3_mnth
     &                 , ht_b4_mnth,ht_b5_mnth,ht_b6_mnth
     &                 , ht_v1_mnth,ht_v2_mnth,ht_v3_mnth
     &                 , ht_v4_mnth,ht_v5_mnth,ht_v6_mnth
     &                 , dzsn_b1,dzsn_b2,dzsn_b3,dzsn_v1,dzsn_v2,dzsn_v3
     &                 , wsn_b1 ,wsn_b2 ,wsn_b3 ,wsn_v1 ,wsn_v2 ,wsn_v3
     &                 , hsn_b1 ,hsn_b2 ,hsn_b3 ,hsn_v1 ,hsn_v2 ,hsn_v3
     &                 , nsn_b,nsn_v,frsnow_b,frsnow_v,Qf_mnth
     &                 , abetad_mnth
     &  , aClivepool_leaf_m,aClivepool_froot_m,aClivepool_wood_m
     &  , aCdeadpool_surfmet_m,aCdeadpool_surfstr_m,aCdeadpool_soilmet_m
     &  , aCdeadpool_soilstr_m,aCdeadpool_cwd_m,aCdeadpool_surfmic_m
     &  , aCdeadpool_soilmic_m,aCdeadpool_slow_m,aCdeadpool_passive_m
     &  , alai_m,canopyH2O_m,canopyheat_m
#endif
     &                     )

      call lsm_init_bc( bc )
      call lsm_init_state( lsm_state, bc%fearth, jday, jyear )

#ifdef USE_GSWP_FORCINGS
      ! initialize GSWP forcings
      call init_forcings(i_site,j_site,latd
     &              ,jday,jyear,(jday-1)*24*3600,dtsec,iu_vector
     &              ,tbcs_ij,tcanopy_ij)

      end_of_input_flag = .false.
#endif

      !call sysusage(mype,0)
      !call sysusage(mype,1)

      !call sysusage(mype+4,0)
      !call sysusage(mype+8,0)
      !call sysusage(mype+12,0)
      num_input_steps = dble(itime1)*dtsec/10800.d0
      time = 0.d0
      
      ! Main Time Loop
      loop_time_main: do itime=itime0,itime1

        !call MPI_Barrier(MPI_COMM_WORLD, ierr)

        itime_3hr = tyr_sec/10800 + 1 ! 10800 is # of sec in 3 hr
        call get_month(jyear,itime_3hr,n_month) ! month # from start of year

        if ( jday .ne. jday_old ) then ! new day
          end_of_day_flag = .true.
          jday_old = jday
        else
          end_of_day_flag = .false.
        endif

!       Determine if last timestep of GSWP2 input data
#ifndef ECOSYSTEM_SCALE
        if(dble(itime)*dtsec/10800.d0 > (num_input_steps-1.d0) )then
           end_of_input_flag = .true.
        endif
#endif

!       Run the LSM one timestep of size dtsec
        call lsm_run(lsm_state,bc,jday,jyear,tyr_sec,dtsec
     &                 , j0,j1,n_month
     &                 , iu_LAI,iu_vht
#ifdef PRINT_DIAGS
     &                 , aevap_mnth,aevapw_mnth,aevapd_mnth,aevapb_mnth
     &                 , aintercep_mnth,aruns_mnth,arunu_mnth
     &                 , agpp_mnth,arauto_mnth,asoilresp_mnth
     &                 , aevapvg_mnth,aevapvs_mnth,aevapbs_mnth
     &                 , asrht_mnth,atrht_mnth,aalbedo_mnth
     &                 , aclab_mnth,asoilCpoolsum_mnth,aepp_mnth
     &                 , atrg_mnth,ashg_mnth,alhg_mnth,aepc_mnth
     &                 , n_count_mnth
     &                 , w_b1_mnth,w_b2_mnth,w_b3_mnth
     &                 , w_b4_mnth,w_b5_mnth,w_b6_mnth
     &                 , w_v1_mnth,w_v2_mnth,w_v3_mnth
     &                 , w_v4_mnth,w_v5_mnth,w_v6_mnth
     &                 , ht_b1_mnth,ht_b2_mnth,ht_b3_mnth
     &                 , ht_b4_mnth,ht_b5_mnth,ht_b6_mnth
     &                 , ht_v1_mnth,ht_v2_mnth,ht_v3_mnth
     &                 , ht_v4_mnth,ht_v5_mnth,ht_v6_mnth
     &                 , dzsn_b1,dzsn_b2,dzsn_b3,dzsn_v1,dzsn_v2,dzsn_v3
     &                 , wsn_b1 ,wsn_b2 ,wsn_b3 ,wsn_v1 ,wsn_v2 ,wsn_v3
     &                 , hsn_b1 ,hsn_b2 ,hsn_b3 ,hsn_v1 ,hsn_v2 ,hsn_v3
     &                 , nsn_b,nsn_v,frsnow_b,frsnow_v,Qf_mnth
     &                 , abetad_mnth
     &  , aClivepool_leaf_m,aClivepool_froot_m,aClivepool_wood_m
     &  , aCdeadpool_surfmet_m,aCdeadpool_surfstr_m,aCdeadpool_soilmet_m
     &  , aCdeadpool_soilstr_m,aCdeadpool_cwd_m,aCdeadpool_surfmic_m
     &  , aCdeadpool_soilmic_m,aCdeadpool_slow_m,aCdeadpool_passive_m
     &  , alai_m,canopyH2O_m,canopyheat_m
#endif
#ifdef USE_GSWP_FORCINGS
     &       ,iu_vector,tbcs_ij,tcanopy_ij
     &       ,end_of_input_flag
#endif
     &       ,end_of_day_flag
     &       )

!       Update time-related variables
        tyr_sec = tyr_sec + nint(dtsec)
        jday = tyr_sec/(24*3600) + 1
        time = time + dtsec

!       End of year check: reset time variables and print diags
        if (      (jday > 365 .and. mod(jyear,4).ne.0) ! non-leap
     &       .or.  jday > 366 ) then                   ! leap
          jday = 1
          tyr_sec = 0
          if (do_spinup) then     !Spin-up check
             if (i_spin<num_times_spin)then !Reset-year check
                jyear = jyear + 1
                if (jyear>spin_end_yr)then
                   i_spin = i_spin + 1
                   jyear = spin_start_yr
                endif
#ifdef ECOSYSTEM_SCALE
                rewind(iu_vector) !Rewind meteorological input files
                rewind(iu_LAI)    !Rewind LAI input file
                rewind(iu_vht)    !Rwind vheight, which is time series now!
#endif
             else                 !No reset, just increase year 
                jyear = jyear + 1
             endif

          else                    !No spinup, increas yr & rewind yearly files
             jyear = jyear + 1
#ifdef ECOSYSTEM_SCALE
             rewind(iu_LAI)       !Rewind LAI input file
             rewind(iu_vht)
#endif
          endif

#ifdef PRINT_DIAGS
          call print_diags(
     &                   aevap_mnth,aevapw_mnth,aevapd_mnth,aevapb_mnth
     &                 , aintercep_mnth,aruns_mnth,arunu_mnth
     &                 , agpp_mnth,arauto_mnth,asoilresp_mnth
     &                 , aevapvg_mnth,aevapvs_mnth,aevapbs_mnth
     &                 , asrht_mnth,atrht_mnth,aalbedo_mnth
     &                 , aclab_mnth,asoilCpoolsum_mnth,aepp_mnth
     &                 , atrg_mnth,ashg_mnth,alhg_mnth,aepc_mnth
     &                 , n_count_mnth
     &                 , w_b1_mnth,w_b2_mnth,w_b3_mnth
     &                 , w_b4_mnth,w_b5_mnth,w_b6_mnth
     &                 , w_v1_mnth,w_v2_mnth,w_v3_mnth
     &                 , w_v4_mnth,w_v5_mnth,w_v6_mnth
     &                 , ht_b1_mnth,ht_b2_mnth,ht_b3_mnth
     &                 , ht_b4_mnth,ht_b5_mnth,ht_b6_mnth
     &                 , ht_v1_mnth,ht_v2_mnth,ht_v3_mnth
     &                 , ht_v4_mnth,ht_v5_mnth,ht_v6_mnth
     &                 , dzsn_b1,dzsn_b2,dzsn_b3,dzsn_v1,dzsn_v2,dzsn_v3
     &                 , wsn_b1 ,wsn_b2 ,wsn_b3 ,wsn_v1 ,wsn_v2 ,wsn_v3
     &                 , hsn_b1 ,hsn_b2 ,hsn_b3 ,hsn_v1 ,hsn_v2 ,hsn_v3
     &                 , nsn_b,nsn_v,frsnow_b,frsnow_v,Qf_mnth
     &                 , abetad_mnth
     &  , aClivepool_leaf_m,aClivepool_froot_m,aClivepool_wood_m
     &  , aCdeadpool_surfmet_m,aCdeadpool_surfstr_m,aCdeadpool_soilmet_m
     &  , aCdeadpool_soilstr_m,aCdeadpool_cwd_m,aCdeadpool_surfmic_m
     &  , aCdeadpool_soilmic_m,aCdeadpool_slow_m,aCdeadpool_passive_m
     &  , alai_m,canopyH2O_m,canopyheat_m)

          call zero_diags(
     &                   aevap_mnth,aevapw_mnth,aevapd_mnth,aevapb_mnth
     &                 , aintercep_mnth,aruns_mnth,arunu_mnth
     &                 , agpp_mnth,arauto_mnth,asoilresp_mnth
     &                 , aevapvg_mnth,aevapvs_mnth,aevapbs_mnth
     &                 , asrht_mnth,atrht_mnth,aalbedo_mnth
     &                 , aclab_mnth,asoilCpoolsum_mnth,aepp_mnth
     &                 , atrg_mnth,ashg_mnth,alhg_mnth,aepc_mnth
     &                 , n_count_mnth
     &                 , w_b1_mnth,w_b2_mnth,w_b3_mnth
     &                 , w_b4_mnth,w_b5_mnth,w_b6_mnth
     &                 , w_v1_mnth,w_v2_mnth,w_v3_mnth
     &                 , w_v4_mnth,w_v5_mnth,w_v6_mnth
     &                 , ht_b1_mnth,ht_b2_mnth,ht_b3_mnth
     &                 , ht_b4_mnth,ht_b5_mnth,ht_b6_mnth
     &                 , ht_v1_mnth,ht_v2_mnth,ht_v3_mnth
     &                 , ht_v4_mnth,ht_v5_mnth,ht_v6_mnth
     &                 , dzsn_b1,dzsn_b2,dzsn_b3,dzsn_v1,dzsn_v2,dzsn_v3
     &                 , wsn_b1 ,wsn_b2 ,wsn_b3 ,wsn_v1 ,wsn_v2 ,wsn_v3
     &                 , hsn_b1 ,hsn_b2 ,hsn_b3 ,hsn_v1 ,hsn_v2 ,hsn_v3
     &                 , nsn_b,nsn_v,frsnow_b,frsnow_v,Qf_mnth
     &                 , abetad_mnth
     &  , aClivepool_leaf_m,aClivepool_froot_m,aClivepool_wood_m
     &  , aCdeadpool_surfmet_m,aCdeadpool_surfstr_m,aCdeadpool_soilmet_m
     &  , aCdeadpool_soilstr_m,aCdeadpool_cwd_m,aCdeadpool_surfmic_m
     &  , aCdeadpool_soilmic_m,aCdeadpool_slow_m,aCdeadpool_passive_m
     &  , alai_m,canopyH2O_m,canopyheat_m)

#endif
        endif ! End of year check
      enddo loop_time_main

      !call sysusage(mype,2)
      !call sysusage(mype,3)

      !call sysusage(mype+4,3)
      !call sysusage(mype+8,3)
      !call sysusage(mype+12,3)

#ifdef MMAP_OUTPUT
      offset = 0
      idim = im 
      jdim = jm 
      j0 = 1
      j1 = jm

      if (force_VEG) then
        call closeunit(iu_LAI)
        call closeunit(iu_vht)
      end if

      call openunit("diag1",fd,.true.,.false.)
      write(fd) diags_ij(:,:,1)
      call closeunit(fd)

      call mmap_open( fd, "diag1$", "w")
      if ( fd < 0 ) then
        print *,"can't open"
        call stop_model("Can't open",255)
      endif
      call mmap_write(diags_ij(1,1,1), fd, offset,
     &     idim, jdim, j0, j1)
      call mmap_close( fd )
#endif

#ifdef USE_ESMF
      call MPI_Finalize(ierr)
#endif

      end subroutine LSM_standalone
