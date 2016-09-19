!sum Routine to run the standalone NASA GISS Land Surface Model (LSM) 
!@+   coupled to the Ent Dynamic Global Terrestrial Ecosystem Model (DGTEM).
!@+   Options for meteorological forcings curently include:
!@+     1) Global forcing.
!@+     2) Single grid cell forcing read from global files.
!@+          Specify run boundaries in rundeck, I0,I1,J0,J1.
!@+     3) Single ecosystem forcing read from single file forcing.
!@+          Compiler switch LSM_DRV=SINGLE.
!@auth I. Aleinov, M.J. Puma, Y. Kim, N.Y. Kiang
!-----------------------------------------------------------------------
!@+   Subroutine to read forcings from different datasets are in different
!@+   fortran module with the same module name and same main interface
!@+   subroutine names, but different file names.  Vegetation forcings are
!@+   handled in a separate module, drv_veg.
!@+   OPTS_giss_LSM_standalone compiler switch LSM_DRV can be:
!@+      default, not set:   old lsm_standalone.f & drv_gswp_forcings.f
!@+      GLOBAL:             drv_met is drv_met_gswp_nk.f
!@+      SINGLE:             drv_met is drv_met_single_nk.f
!-----------------------------------------------------------------------
!@+   Given uniqueness of forcing datasets, some specific compiler directives
!@+   may be necessary.  Define here.
#define USE_GSWP_FORCINGS
#define PRINT_DIAGS
!-----------------------------------------------------------------------
      module lsm

      use CONSTANT, only : undef
      use FILEMANAGER
      use ent_mod
      use ghy_h, only : ngm, imt, nlsn, LS_NFRAC
      use sle001, only : advnc
      use lsm_phys_util, only : UNDEFINT

      implicit none

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

#ifdef PRINT_DIAGS  /*  Monthly diagnostic accumulators - Global to module */
      character*80 :: title       ! title of diagnostic to print to file
      real*8, allocatable, dimension(:,:,:) :: !dimension(I0:I1,J0:J1,12)
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

      contains

!-----------------------------------------------------------------------
#ifdef READ_WITH_SYNCPARAM
      subroutine read_input_parameters(
!     Called by process_input_parameters to read input parameter file
!     and process binary flags to logical.  See parameter descriptions there.
     &     IM,JM,I0,I1,J0,J1
     &     ,i0v, i1v, j0v, j1v         ,i0vh, i1vh, j0vh, j1vh
     &     ,i0vc, i1vc, j0vc, j1vc     ,i0cc, i1cc, j0cc, j1cc
     &     ,i0h, i1h, j0h, j1h         ,i0s, i1s, j0s, j1s
     &     ,i0sc, i1sc, j0sc, j1sc     ,i0st, i1st, j0st, j1st
     &     ,i0d, i1d, j0d, j1d
     &     ,year1,year2, tyr_sec1,tyr_sec2 
     &     ,dtsec, Mdtsec, Mtstart, Vdtsec, Vrows, Vtstart
     &     ,spin_start_yr, spin_end_yr, num_times_spin
     &     ,force_VEG
     &     ,do_soilinit, do_soilresp
     &     ,do_phenology_activegrowth, do_structuralgrowth
     &     ,do_frost_hardiness, do_patchdynamics, do_init_geo, do_spinup
     &     ,BASEFOLD_OUT)

      use filemanager
      use parser
      use param

      integer ::  IM,JM,I0,I1,J0,J1
     &     ,i0v, i1v, j0v, j1v         ,i0vh, i1vh, j0vh, j1vh
     &     ,i0vc, i1vc, j0vc, j1vc     ,i0cc, i1cc, j0cc, j1cc
     &     ,i0h, i1h, j0h, j1h         ,i0s, i1s, j0s, j1s
     &     ,i0sc, i1sc, j0sc, j1sc     ,i0st, i1st, j0st, j1st
     &     ,i0d, i1d, j0d, j1d
     &     ,year1,year2, tyr_sec1,tyr_sec2 
      real*8 :: dtsec
      integer :: Mdtsec, Mtstart, Vdtsec, Vrows, Vtstart
     &     ,spin_start_yr, spin_end_yr, num_times_spin
!     real*8 :: latd
      logical :: force_VEG
     &     ,do_soilinit, do_soilresp
     &     ,do_phenology_activegrowth, do_structuralgrowth
     &     ,do_frost_hardiness, do_patchdynamics, do_init_geo, do_spinup
      character*80 :: BASEFOLD_OUT

!------ Local ------
      character*120 :: ifile="ent_input"
      integer :: iu_IFILE
      integer :: force_VEG_bi
     &     ,do_soilinit_bi, do_soilresp_bi
     &     ,do_phenology_activegrowth_bi, do_structuralgrowth_bi
     &     ,do_frost_hardiness_bi, do_patchdynamics_bi, do_init_geo_bi
     &     ,do_spinup_bi

      namelist /PARAMETERS/ IM, JM, I0, I1, J0, J1
     &     ,i0v, i1v, j0v, j1v
     &     ,i0v, i1v, j0v, j1v         ,i0vh, i1vh, j0vh, j1vh
     &     ,i0vc, i1vc, j0vc, j1vc     ,i0cc, i1cc, j0cc, j1cc
     &     ,i0h, i1h, j0h, j1h         ,i0s, i1s, j0s, j1s
     &     ,i0sc, i1sc, j0sc, j1sc     ,i0st, i1st, j0st, j1st
     &     ,i0d, i1d, j0d, j1d
     &     ,year1,year2, tyr_sec1,tyr_sec2 
     &     ,dtsec, Mdtsec, Mtstart,Vdtsec, Vrows, Vtstart
     &     ,spin_start_yr, spin_end_yr, num_times_spin
     &     ,force_VEG_bi
     &     ,do_soilinit_bi,do_soilresp_bi
     &     ,do_phenology_activegrowth_bi,do_structuralgrowth_bi
     &     ,do_frost_hardiness_bi,do_patchdynamics_bi,do_init_geo_bi
     &     ,do_spinup_bi
     &     ,BASEFOLD_OUT


! Read input file with run settings
      print *,"reading input parameters from ", ifile
      call openunit(trim(ifile),iu_IFILE,.false.,.true.)
      call parse_params(iu_IFILE)
      call closeunit(iu_IFILE)
      print *, "parsed input parameters"

!Convert binary integer flags to logical
      call sync_param( "IM", IM )
      print *, "IM", IM
      call sync_param( "JM", JM )
      call sync_param( "I0", I0 )
      call sync_param( "I1", I1 )
      call sync_param( "J0", J0 )
      call sync_param( "J1", J1 )
      call sync_param( "i0v", i0v )
      call sync_param( "i1v", i1v )
      call sync_param( "j0v", j0v )
      call sync_param( "j1v", j1v ) 
      call sync_param( "i0vh", i0vh )
      call sync_param( "i1vh", i1vh )
      call sync_param( "j0vh", j0vh )
      call sync_param( "j1vh", j1vh )
      call sync_param( "i0vc", i0vc )
      call sync_param( "i1vc", i1vc )
      call sync_param( "j0vc", j0vc )
      call sync_param( "j1vc", j1vc )
      call sync_param( "i0cc", i0cc )
      call sync_param( "i1cc", i1cc )
      call sync_param( "j0cc", j0cc )
      call sync_param( "j1cc", j1cc )
      call sync_param( "i0h", i0h )
      call sync_param( "i1h", i1h )
      call sync_param( "j0h", j0h )
      call sync_param( "j1h", j1h )
      call sync_param( "i0s", i0s )
      call sync_param( "i1s", i1s )
      call sync_param( "j0s", j0s )
      call sync_param( "j1s", j1s )
      call sync_param( "i0sc", i0sc )
      call sync_param( "i1sc", i1sc )
      call sync_param( "j0sc", j0sc )
      call sync_param( "j1sc", j1sc )
      call sync_param( "i0st", i0st )
      call sync_param( "i1st", i1st )
      call sync_param( "j0st", j0st )
      call sync_param( "j1st", j1st )
      call sync_param( "i0d", i0d )
      call sync_param( "i1d", i1d )
      call sync_param( "j0d", j0d )
      call sync_param( "j1d", j1d )
      call sync_param( "year1", year1 )
      call sync_param( "year2", year2 )
      call sync_param( "tyr_sec1", tyr_sec1 )
      call sync_param( "tyr_sec2", tyr_sec2 )
      call sync_param( "dtsec", dtsec )
      call sync_param( "Mdtsec", Mdtsec )
      call sync_param( "Mtstart", Mtstart )
      call sync_param( "Vdtsec", Vdtsec )
      call sync_param( "Vrows", Vrows )
      call sync_param( "Vtstart", Vtstart )
      call sync_param( "spin_start_yr", spin_start_yr )
      call sync_param( "spin_end_yr", spin_end_yr )
      call sync_param( "num_times_spin", num_times_spin )
      call sync_param( "force_VEG", force_VEG_bi )
      call sync_param( "do_soilinit", do_soilinit_bi )
      call sync_param( "do_soilresp", do_soilresp_bi )
      call sync_param( "do_phenology_activegrowth", 
     &     do_phenology_activegrowth_bi )
      call sync_param( "do_structuralgrowth", do_structuralgrowth_bi )
      call sync_param( "do_frost_hardiness", do_frost_hardiness_bi )
      call sync_param( "do_patchdynamics", do_patchdynamics_bi )
      call sync_param( "do_init_geo", do_init_geo_bi )
      call sync_param( "do_spinup", do_spinup_bi )

      call sync_param( "BASEFOLD_OUT", basefold_out )

!     Convert binary flags to logicals.
!     If using GISS GCM sync_param, need to convert binary flags to logicals.
      if (force_VEG_bi==1) force_VEG=.true.
      if (do_soilinit_bi==1) do_soilinit=.true. 
      if (do_soilresp_bi==1) do_soilresp=.true. 
      if (do_phenology_activegrowth_bi==1)  
     &     do_phenology_activegrowth=.true.  
      if (do_structuralgrowth_bi==1) do_structuralgrowth=.true.
      if (do_frost_hardiness_bi==1) do_frost_hardiness=.true.
      if (do_patchdynamics_bi==1) do_patchdynamics=.true.
      if (do_init_geo_bi==1) do_init_geo=.true.
      if (do_spinup_bi==1) do_spinup=.true.

      print *,"finished reading input parameters from ", ifile

      end subroutine read_input_parameters
#endif

!-----------------------------------------------------------------------

      subroutine process_input_parameters(
     &     IM,JM,I0,I1,J0,J1
     &     ,i0v, i1v, j0v, j1v         ,i0vh, i1vh, j0vh, j1vh
     &     ,i0vc, i1vc, j0vc, j1vc     ,i0cc, i1cc, j0cc, j1cc
     &     ,i0h, i1h, j0h, j1h         ,i0s, i1s, j0s, j1s
     &     ,i0sc, i1sc, j0sc, j1sc     ,i0st, i1st, j0st, j1st
     &     ,i0d, i1d, j0d, j1d
     &     ,year1,year2, tyr_sec1,tyr_sec2 
     &     ,dtsec, Mdtsec, Mtstart, Vdtsec, Vrows, Vtstart
     &     ,spin_start_yr, spin_end_yr, num_times_spin
!     &    !, latd
     &     ,force_VEG
     &     ,do_soilinit, do_soilresp
     &     ,do_phenology_activegrowth, do_structuralgrowth
     &     ,do_frost_hardiness, do_patchdynamics, do_init_geo, do_spinup
     &     ,BASEFOLD_OUT)
 
! Read run parameters, calculate default values for parameters not input,
! and check for errors in inputs.
! NOTE:  No error check is done for FILE_dtsec or FILE_rows, because these depend
!        on the format of driver files.  Some driver files include a whole year, 
!        while other datasets may split the drivers into monthly files.

      implicit none
!@var Grid dimensions
      integer :: IM,JM,I0,I1,J0,J1
!@var Grid indices for files for various boundary condition or drivers.
      integer :: i0v, i1v,j0v, j1v     !Vegetation LAI driver file boundaries. 
      integer :: i0vh, i1vh,j0vh, j1vh !Vegetation height driver file boundaries.
      integer :: i0vc, i1vc,j0vc, j1vc !Vegetation cover file boundaries.
      integer :: i0cc, i1cc,j0cc, j1cc !Crop cover file boundaries.
      integer :: i0sc, i1sc,j0sc, j1sc !Soil carbon file boundaries.
      integer :: i0st, i1st,j0st, j1st !Soil texture file boundaries.
      integer :: i0h, i1h,j0h, j1h     !Hydrology LSM boundary condition files boundaries.
      integer :: i0s, i1s,j0s, j1s     !Hydrology state files boundaries.
      integer :: i0d, i1d,j0d, j1d     !Meteorological driver files boundaries.

!@var Gregorian year start and finish (julian year for modelE met forcings)
      integer :: year1,year2
!@var Run start time in year1, end time in year2 (sec from beginning of year)
!     Seconds make it consistent with restart files.
      integer :: tyr_sec1,tyr_sec2
!@var Timestep size of model run [seconds]
      real*8  :: dtsec
!@var Timestep size of met forcings data [seconds], GSWP2 is 3 hr in seconds
      integer  :: Mdtsec
!@var Number of rows in a met data forcing file, at time steps of Mdtsec.
!      integer :: Mrows
!@var Year A.D. of start of met data forcings.  THIS ASSUMES ALL FILES START Jan 1.
!+    IF DATA DO NOT START ON JAN 1, THEN DUMMY ROWS MUST BE ADDED.
      integer :: Mtstart
!@var Timestep size of veg forcings data [seconds]
      integer :: Vdtsec
!@var Number of rows in a veg data forcing file, at time steps of VEGFILE_dtsec.
      integer :: Vrows
!@var Year A.D. of start of met data forcings.  THIS ASSUMES ALL FILES START Jan 1.
!+    IF DATA DO NOT START ON JAN 1, THEN DUMMY ROWS MUST BE ADDED.
      integer :: Vtstart
!!@var Spin-up bounds.
      integer :: spin_start_yr, spin_end_yr, num_times_spin
!     Flag for LSM drivers.
!      real*8 :: latd            !Latitude, used only for single grid runs.
                                !Should replace with I0,I1,J0,J1
      logical :: force_VEG
      logical :: do_soilinit,do_soilresp
      logical :: do_phenology_activegrowth,do_structuralgrowth
      logical :: do_frost_hardiness,do_patchdynamics,do_init_geo
      logical :: do_spinup
      character*80 :: BASEFOLD_OUT
      !---Local-----------
      character*80 :: ifile="ent_input"
      integer :: iu_IFILE
!     Binary flags for the Ent DGTEM
      integer :: force_VEG_bi
      integer :: do_soilinit_bi,do_soilresp_bi
      integer :: do_phenology_activegrowth_bi,do_structuralgrowth_bi
      integer :: do_frost_hardiness_bi,do_patchdynamics_bi
      integer :: do_init_geo_bi
      integer :: do_spinup_bi

!#ifndef READ_WITH_SYNCPARAM
      namelist /PARAMETERS/ IM, JM, I0, I1, J0, J1
     &     ,i0v, i1v, j0v, j1v
     &     ,i0v, i1v, j0v, j1v         ,i0vh, i1vh, j0vh, j1vh
     &     ,i0vc, i1vc, j0vc, j1vc     ,i0cc, i1cc, j0cc, j1cc
     &     ,i0h, i1h, j0h, j1h         ,i0s, i1s, j0s, j1s
     &     ,i0sc, i1sc, j0sc, j1sc     ,i0st, i1st, j0st, j1st
     &     ,i0d, i1d, j0d, j1d
     &     ,year1,year2, tyr_sec1,tyr_sec2 
     &     ,dtsec, Mdtsec, Mtstart,Vdtsec, Vrows, Vtstart
     &     ,spin_start_yr, spin_end_yr, num_times_spin
     &     ,force_VEG_bi
     &     ,do_soilinit_bi,do_soilresp_bi
     &     ,do_phenology_activegrowth_bi,do_structuralgrowth_bi
     &     ,do_frost_hardiness_bi,do_patchdynamics_bi,do_init_geo_bi
     &     ,do_spinup_bi
     &     ,BASEFOLD_OUT
     &     ,force_VEG
     &     ,do_soilinit, do_soilresp
     &     ,do_phenology_activegrowth,do_structuralgrowth
     &     ,do_frost_hardiness,do_patchdynamics,do_init_geo,do_spinup
     &     ,BASEFOLD_OUT

      call zero_bounds(I0, I1, J0, J1
     &     ,i0v, i1v, j0v, j1v
     &     ,i0vh, i1vh, j0vh, j1vh
     &     ,i0vc, i1vc, j0vc, j1vc
     &     ,i0cc, i1cc, j0cc, j1cc
     &     ,i0sc, i1sc, j0sc, j1sc
     &     ,i0st, i1st, j0st, j1st
     &     ,i0h, i1h, j0h, j1h
     &     ,i0s, i1s, j0s, j1s
     &     ,i0d, i1d, j0d, j1d)

      !Initialize optional input parameters.- run times
      year1 = UNDEFINT
      year2 = UNDEFINT
      tyr_sec1 = UNDEFINT
      tyr_sec2 = UNDEFINT
      dtsec = 1800.             !Default half-hour (seconds)
      Mdtsec = UNDEFINT
      Mtstart = UNDEFINT
      Vdtsec = UNDEFINT
      Vrows = UNDEFINT
      Vtstart = UNDEFINT
      spin_start_yr = UNDEFINT
      spin_end_yr = UNDEFINT
      num_times_spin = 0

      !Initialize array bounds
      IM=72    !Grid resolution
      JM=46

      ! Default flags
      force_VEG=.false.            
      do_soilinit=.false.
      do_soilresp=.false.
      do_phenology_activegrowth=.false.
      do_structuralgrowth=.false.
      do_frost_hardiness=.true.  !This should be default true!
      do_patchdynamics=.false.
      do_init_geo=.false.
      do_spinup=.false.

!     Default latd (single grid runs only)
!     latd = undef

!#ifdef READ_WITH_SYNCPARAM
!       call read_input_parameters(
!     &        IM,JM,I0,I1,J0,J1
!     &        ,i0v, i1v, j0v, j1v         ,i0vh, i1vh, j0vh, j1vh
!     &        ,i0vc, i1vc, j0vc, j1vc     ,i0cc, i1cc, j0cc, j1cc
!     &        ,i0h, i1h, j0h, j1h         ,i0s, i1s, j0s, j1s
!     &        ,i0sc, i1sc, j0sc, j1sc     ,i0st, i1st, j0st, j1st
!     &        ,i0d, i1d, j0d, j1d
!     &        ,year1,year2, tyr_sec1,tyr_sec2 
!     &        ,dtsec, Mdtsec, Mtstart, Vdtsec, Vrows, Vtstart
!     &        ,spin_start_yr, spin_end_yr, num_times_spin
!     &        ,force_VEG
!     &        ,do_soilinit, do_soilresp
!     &        ,do_phenology_activegrowth, do_structuralgrowth
!     &        ,do_frost_hardiness, do_patchdynamics, do_spinup
!     &        ,BASEFOLD_OUT)
!#else
!       print *,"reading input parameters from ", ifile
!       call openunit(trim(ifile),iu_IFILE,.false.,.true.)
!       read(iu_IFILE, NML=PARAMETERS, ERR=10)
!       call closeunit(iu_IFILE)
!       print *,"finished reading input parameters from ", ifile
!#endif

      ! Read input file with run settings
      print *,"reading input parameters from ", ifile
      call openunit(trim(ifile),iu_IFILE,.false.,.true.)
      read(iu_IFILE, NML=PARAMETERS, ERR=10)
      call closeunit(iu_IFILE)
      print *,"finished reading input parameters from ", ifile

      call check_time_bounds(dtsec,year1,year2,tyr_sec1,tyr_sec2
     &      ,spin_start_yr, spin_end_yr, num_times_spin)

#ifdef LSM_DRV_SINGLE      
      call set_single_bounds(force_VEG, do_soilinit
     &     ,IM,JM,I0,I1,J0,J1
     &     ,i0v, i1v, j0v, j1v         ,i0vh, i1vh, j0vh, j1vh
     &     ,i0vc, i1vc, j0vc, j1vc     ,i0cc, i1cc, j0cc, j1cc
     &     ,i0sc, i1sc, j0sc, j1sc     ,i0st, i1st, j0st, j1st
     &     ,i0h, i1h, j0h, j1h         ,i0d, i1d, j0d, j1d)
         !Override all grid bounds with defaults for single-grid,
         !single-file driver bounds.
#endif

      call check_grid_bounds(force_VEG, do_soilinit
     &     ,IM,JM,I0,I1,J0,J1
     &     ,i0v, i1v, j0v, j1v         ,i0vh, i1vh, j0vh, j1vh
     &     ,i0vc, i1vc, j0vc, j1vc     ,i0cc, i1cc, j0cc, j1cc
     &     ,i0sc, i1sc, j0sc, j1sc     ,i0st, i1st, j0st, j1st
     &     ,i0h, i1h, j0h, j1h         ,i0s, i1s, j0s, j1s
     &     ,i0d, i1d, j0d, j1d)

!YKIM
!possible combinations of options
!do_phenology_activegrowth=true & do_structuralgrowth=true
!do_phenology_activegrowth=true & do_structuralgrowth=false
!do_phenology_activegrowth=false & do_structuralgrowth=false
!impossible combinations of options
!do_phenology_activegrowth=false & do_structuralgrowth=true
      if (.not.do_phenology_activegrowth .and.
     &   do_structuralgrowth) then
         print*,"Impossible combinations of input flags."
         call stop_model("Impossible combinations of input flags.",255)
      endif

      print *,"Input parameters:"
      print *, "Boundary parameters: "
     &     ,IM,JM,I0,I1,J0,J1, "\n"
     &     ,i0v, i1v, j0v, j1v      ,i0vh, i1vh, j0vh, j1vh, "\n"
     &     ,i0vc, i1vc, j0vc, j1vc  ,i0cc, i1cc, j0cc, j1cc, "\n"
     &     ,i0sc, i1sc, j0sc, j1sc  ,i0st, i1st, j0st, j1st, "\n"
     &     ,i0h, i1h, j0h, j1h      ,i0s, i1s, j0s, j1s, "\n"
     &     ,i0d, i1d, j0d, j1d
      print *,"Time parameters: "
     &     ,year1,year2, tyr_sec1,tyr_sec2 
     &     ,dtsec, Mdtsec, Mtstart, Vdtsec, Vrows, Vtstart
     &     ,spin_start_yr, spin_end_yr, num_times_spin
!     &    !, latd
      print *,"Ent flags: "
     &     ,force_VEG
     &     ,do_soilinit, do_soilresp
     &     ,do_phenology_activegrowth, do_structuralgrowth
     &     ,do_frost_hardiness, do_patchdynamics, do_init_geo, do_spinup
      print *,"BASEFOLD_OUT:", BASEFOLD_OUT

      return
 10   continue
      print *,"error reading namelist file:", ifile
      stop 255

      end subroutine process_input_parameters

!-----------------------------------------------------------------------
      subroutine zero_bounds( I0, I1, J0, J1
     &     ,i0v, i1v, j0v, j1v
     &     ,i0vh, i1vh, j0vh, j1vh
     &     ,i0vc, i1vc, j0vc, j1vc
     &     ,i0cc, i1cc, j0cc, j1cc
     &     ,i0sc, i1sc, j0sc, j1sc
     &     ,i0st, i1st, j0st, j1st
     &     ,i0h, i1h, j0h, j1h
     &     ,i0s, i1s, j0s, j1s
     &     ,i0d, i1d, j0d, j1d)
      !@var Grid dimensions
      integer :: IM,JM,I0,I1,J0,J1
      !@var Grid indices for single cell
      integer :: i0v, i1v,j0v, j1v     !Vegetation LAI driver file boundaries. 
      integer :: i0vh, i1vh,j0vh, j1vh !Vegetation height driver file boundaries.
      integer :: i0vc, i1vc,j0vc, j1vc !Vegetation cover file boundaries.
      integer :: i0cc, i1cc,j0cc, j1cc !Crop cover file boundaries.
      integer :: i0sc, i1sc,j0sc, j1sc !Soil carbon file boundaries.
      integer :: i0st, i1st,j0st, j1st !Soil texture file boundaries.
      integer :: i0h, i1h,j0h, j1h     !Hydrology LSM boundary condition files boundaries.
      integer :: i0s, i1s,j0s, j1s     !Hydrology state files boundaries.
      integer :: i0d, i1d,j0d, j1d     !Meteorological driver files boundaries.
      !-----Local------
      integer, parameter :: U = UNDEFINT

      I0=U     !Run boundaries
      I1=U
      J0=U
      J1=U

      i0v=U    !LAI driver file if force_VEG=.true.
      i1v=U
      j0v=U
      j1v=U

      i0vh=U   !Veg HEIGHT driver file if force_VEG=.true.
      i1vh=U
      j0vh=U
      j1vh=U

      i0vc=U   !Vegetation cover file if not MIXED_VEG
      i1vc=U
      j0vc=U
      j1vc=U

      i0cc=U   !Crop cover file
      i1cc=U
      j0cc=U
      j1cc=U

      i0sc=U   !Soil carbon file
      i1sc=U
      j0sc=U
      j1sc=U

      i0st=U   !Soil texture file
      i1st=U
      j0st=U
      j1st=U

      i0h=U    !Hydrology surface topo, boundary conditions files
      i1h=U
      j0h=U
      j1h=U

      i0s=U    !Hydrology state files
      i1s=U
      j0s=U
      j1s=U

      i0d=U    !Meteorological drivers file(s)
      i1d=U
      j0d=U
      j1d=U
      end subroutine zero_bounds

!-----------------------------------------------------------------------
      subroutine check_time_bounds(dtsec,jyear1,jyear2,tyr_sec1,tyr_sec2
     &      ,spin_start_yr, spin_end_yr, num_times_spin)
      use drv_met, only : MAX_METDRV_YEAR, MAX_METDRV_YEARSEC
      use lsm_phys_util, only : LEAPYR
      real*8, intent(in) :: dtsec
      integer,intent(inout) :: jyear1,jyear2,tyr_sec1,tyr_sec2
     &      ,spin_start_yr, spin_end_yr, num_times_spin
      !----------
      logical :: flag = .false.
      integer :: days,dtint

      dtint = dtsec

      !Calculate DEFAULT values for parameters not input.
      if (jyear2.eq.UNDEFINT) jyear2=jyear1 !If year2 not specified, default same as year1
      if (tyr_sec1.lt.0) tyr_sec1 = 0 !Default beginning of year.
      if (LEAPYR(jyear2)) then 
          days = 366
      else
          days = 365
      endif
      if (tyr_sec2.lt.0)  then  !Default to last time step of year
         tyr_sec2 =days*24*60*60 - dtint
      else if (tyr_sec2.eq.days*24*60*60) then
         tyr_sec2 =days*24*60*60 - dtint
         print *,"Correcting tyr_sec2 to last time step of year",
     &     tyr_sec2
      endif
      if ((mod(tyr_sec1,dtint).ne.0).or.(mod(tyr_sec2,dtint).ne.0)) then
         print *,"Error: time bounds not multiple of dtsec"
         call stop_model("Time bounds not multiple of dtsec",255)
      endif

      !Check times make sense.
      if (jyear2.lt.jyear1) flag=.true.
      if ((jyear1.eq.jyear2).and.(tyr_sec1.gt.tyr_sec2)) flag=.true.
      if (flag) then
         print *,"Inconsistent time bound inputs."
         call stop_model("Inconsistent time bound inputs.",255)
      endif
      if (num_times_spin.gt.0) then
         if ((spin_start_yr.eq.UNDEFINT).or.(spin_end_yr.eq.UNDEFINT)) 
     &        flag=.true.
         if (spin_end_yr.lt.spin_start_yr) flag=.true.
      endif
      if (flag) then
         print *,"Inconsistent spin bound inputs."
         call stop_model("Inconsistent spin bound inputs.",255)
      endif
      if ((jyear2.gt.MAX_METDRV_YEAR).and.
     &     (tyr_sec2.gt.MAX_METDRV_YEARSEC)) then
         print *,"End time is greater than available input met data."
         call stop_model("End time exceeds input data time.",255)
      endif
      end subroutine check_time_bounds

!-----------------------------------------------------------------------
      subroutine set_single_bounds(
      !Override all grid bounds with defaults for single-grid,
      !single-file driver bounds.
     &     force_VEG, do_soilinit
     &     ,IM,JM,I0,I1,J0,J1
     &     ,i0v, i1v, j0v, j1v         ,i0vh, i1vh, j0vh, j1vh
     &     ,i0vc, i1vc, j0vc, j1vc     ,i0cc, i1cc, j0cc, j1cc
     &     ,i0sc, i1sc, j0sc, j1sc     ,i0st, i1st, j0st, j1st
     &     ,i0h, i1h, j0h, j1h         !,i0s, i1s, j0s, j1s
     &     ,i0d, i1d, j0d, j1d)
      logical,intent(in) :: force_VEG, do_soilinit
      integer,intent(inout) :: IM,JM,I0,I1,J0,J1
     &     ,i0v, i1v, j0v, j1v         ,i0vh, i1vh, j0vh, j1vh
     &     ,i0vc, i1vc, j0vc, j1vc     ,i0cc, i1cc, j0cc, j1cc
     &     ,i0sc, i1sc, j0sc, j1sc     ,i0st, i1st, j0st, j1st
     &     ,i0h, i1h, j0h, j1h         !,i0s, i1s, j0s, j1s
     &     ,i0d, i1d, j0d, j1d

      if (force_VEG) then
         i0v=I0                 !LAI
         i1v=I1
         j0v=J0
         j1v=J1

         i0vh=I0                !Height
         i1vh=I1
         j0vh=J0
         j1vh=J1
      else
         i0v=UNDEFINT
         i1v=UNDEFINT
         j0v=UNDEFINT
         j1v=UNDEFINT

         i0vh=UNDEFINT
         i1vh=UNDEFINT
         j0vh=UNDEFINT
         j1vh=UNDEFINT
      endif

      i0vc=I0                   !Veg cover
      i1vc=I1
      j0vc=J0
      j1vc=J1

      i0cc=I0                   !Crop cover
      i1cc=I1
      j0cc=J0
      j1cc=J1

      if (do_soilinit) then
         i0sc=I0                !Soil carbon
         i1sc=I1
         j0sc=J0
         j1sc=J1
      else
         i0sc=UNDEFINT
         i1sc=UNDEFINT
         j0sc=UNDEFINT
         j1sc=UNDEFINT
      endif

      !Soil texture may be from global file. Don't override rundeck.
      !i0st=I0                   !Soil texture
      !i1st=I1
      !j0st=J0
      !j1st=J1

      !Hydrology state may be from global file.
      !i0h=I0                    !Hydrology state
      !i1h=I1
      !j0h=J0
      !j1h=J1

      !Land surface boundary conditions may be from global file.
      !i0s=I0                    !Land surface boundary conditions
      !i1s=I1
      !j0s=J0
      !j1s=J1

      i0d=I0                    !Meteorological drivers
      i1d=I1
      j0d=J0
      j1d=J1
 
      end subroutine set_single_bounds

!-----------------------------------------------------------------------
      subroutine check_grid_bounds(force_VEG, do_soilinit
     &     ,IM,JM,I0,I1,J0,J1
     &     ,i0v, i1v, j0v, j1v         ,i0vh, i1vh, j0vh, j1vh
     &     ,i0vc, i1vc, j0vc, j1vc     ,i0cc, i1cc, j0cc, j1cc
     &     ,i0sc, i1sc, j0sc, j1sc     ,i0st, i1st, j0st, j1st
     &     ,i0h, i1h, j0h, j1h         ,i0s, i1s, j0s, j1s
     &     ,i0d, i1d, j0d, j1d)
      logical,intent(in) :: force_VEG, do_soilinit
      integer,intent(inout) ::  IM,JM,I0,I1,J0,J1
     &     ,i0v, i1v, j0v, j1v         ,i0vh, i1vh, j0vh, j1vh
     &     ,i0vc, i1vc, j0vc, j1vc     ,i0cc, i1cc, j0cc, j1cc
     &     ,i0sc, i1sc, j0sc, j1sc     ,i0st, i1st, j0st, j1st
     &     ,i0h, i1h, j0h, j1h         ,i0s, i1s, j0s, j1s
     &     ,i0d, i1d, j0d, j1d

      !----Local--------------------
      logical :: flag = .false.
      integer, parameter :: U = UNDEFINT

      !Calculate DEFAULTS for bounds that were not input
      !Note:  This assumes the user input all of a set of i0,i1,j0,j1, or none.
      !       This does not check for errors if, e.g. I1 was input but I0 was not.
      !       If any one is missing, then all are considered missing.
      if ((I0.eq.U).or.(I1.eq.U).or.(J0.eq.U).or.(J1.eq.U)) then
         I0=1               !Run boundaries:  Assume global for LSM_standalone
         I1=IM
         J0=1
         J1=JM
      endif

      if (((i0v.eq.U).or.(i1v.eq.U).or.(j0v.eq.U).or.(j1v.eq.U))
     &     .and.(force_VEG)) then
         i0v=I0                 !LAI driver file if force_VEG=.true.
         i1v=I1                 !Assume same as run bounds.
         j0v=J0
         j1v=J1
      endif

      if (((i0vh.eq.U).or.(i1vh.eq.U).or.(j0vh.eq.U).or.(j1vh.eq.U))
     &     .and.(force_VEG)) then
         i0vh=i0v               !Veg HEIGHT driver file if force_VEG=.true.
         i1vh=i1v               !Assume same as LAI driver file.
         j0vh=j0v
         j1vh=j1v
      endif

      if ((i0v.eq.U).or.(i1v.eq.U).or.(j0v.eq.U).or.(j1v.eq.U)) then
         i0vc=i0v               !Vegetation cover file if not MIXED_VEG
         i1vc=i1v               !Assume same as LAI and HEIGHT driver files.
         j0vc=j0v
         j1vc=j1v
      endif

      if ((i0cc.eq.U).or.(i1cc.eq.U).or.(j0cc.eq.U).or.(j1cc.eq.U)) then
         i0cc=i0vc               !Crop cover file.  Assume same as veg cover file.
         i1cc=i1vc
         j0cc=j0vc
         j1cc=j1vc
      endif

      if ((i0h.eq.U).or.(i1h.eq.U).or.(j0h.eq.U).or.(j1h.eq.U)) then
         i0h=1                  !Hydrology surface topo, boundary conditions files
         i1h=IM                 !Assume from global dataset
         j0h=1
         j1h=JM
      endif

      if (((i0sc.eq.U).or.(i1sc.eq.U).or.(j0sc.eq.U).or.(j1sc.eq.U))
     &     .and.(do_soilinit)) then
         i0sc=I0                 !Soil carbon file.
         i1sc=I1                 !Assume same as run bounds.
         j0sc=J0
         j1sc=J1
      endif

      if ((i0st.eq.U).or.(i1st.eq.U).or.(j0st.eq.U).or.(j1st.eq.U)) then
         i0st=1                  !Soil texture.
         i1st=IM                 !Assume from global dataset. Better if site-based.
         j0st=1
         j1st=JM
      endif


      if ((i0s.eq.U).or.(i1s.eq.U).or.(j0s.eq.U).or.(j1s.eq.U)) then
         i0s=1                  !Hydrology state files
         i1s=IM                 !Assume from global dataset.
         j0s=1
         j1s=JM
      endif

      if ((i0d.eq.U).or.(i1d.eq.U).or.(j0d.eq.U).or.(j1d.eq.U)) then
         i0d=I0                 !Meteorological drivers file(s)
         i1d=I1                 !Assume same as run bounds.
         j0d=J0
         j1d=J1
      endif

      !Check run boundaries and input files within grid resolution
      if ( (I1>IM).or.(J1>JM) ) flag=.true.
      if (force_VEG) then
         if ( (i1v>IM).or.(j1v>JM).or.(i1vh>IM).or.(j1vh>JM)
     &        .or.(i1vc>IM).or.(j1vc>JM).or.(i1cc>IM).or.(j1cc>JM) ) 
     &        flag=.true.
      endif
      if ( (i1sc>IM).or.(j1sc>JM) )   flag=.true. 
      if ( (i1st>IM).or.(j1st>JM) )   flag=.true. 
      if ( (i1h>IM).or.(j1h>JM).or.(i1d>IM).or.(j1d>JM) )   flag=.true. 
      if ( (i1s>IM).or.(j1s>JM) )   flag=.true. 

      if (flag) then
        print *,"Array boundaries outside grid resolution"

        call stop_model("Array boundaries outside grid resolution",255)
      endif

      !Check run boundaries within input file boundaries.
      if (force_VEG) then
         if ( (I0<i0v).or.(J0<j0v).or.(I0<i0vh).or.(J0<j0vh)
     &        .or.(I0<i0vc).or.(J0<j0vc).or.(I0<i0cc).or.(J0<j0cc) ) 
     &        flag=.true.
      endif
      if ( (I0<i0sc).or.(J0<j0sc) )   flag=.true. 
      if ( (I0<i0st).or.(J0<j0st) )   flag=.true. 
      if ( (I0<i0h).or.(J0<j0h).or.(I0<i0d).or.(j0<j0d) )   flag=.true. 
      if ( (I0<i0s).or.(J0<j0s) )   flag=.true. 

      if (force_VEG) then
         if ( (I1>i1v).or.(J1>j1v).or.(I1>i1vh).or.(J1>j1vh)
     &        .or.(I1>i1vc).or.(J1>j1vc).or.(I1>i1cc).or.(J1>j1cc) ) 
     &        flag=.true.
      endif
      if ( (I1>i1sc).or.(J1>j1sc) )   flag=.true. 
      if ( (I1>i1st).or.(J1>j1st) )   flag=.true. 
      if ( (I1>i1h).or.(J1>j1h).or.(I1>i1d).or.(j1>j1d) )   flag=.true. 
      if ( (I1>i1s).or.(J1>j1s) )   flag=.true. 

      if (flag) then
        print *,"Run boundaries outside input file boundaries"
        call stop_model("Run bounds outside input file bounds",255)
      endif

      end subroutine check_grid_bounds
!-----------------------------------------------------------------------
      subroutine lsm_allocate_state(IM, JM, I0, I1, J0, J1,fearth,s)
      use ghy_h, only : ngm, nlsn, LS_NFRAC
      integer :: IM, JM !Grid resolution
      integer :: I0, I1, J0, J1 !Run boundaries
      real*8 :: fearth(I0:I1,J0:J1)
      type(t_lsm_state) :: s
      !----Local-------
      integer :: i,j

      !Allocate run-size arrays for Ent.
      allocate( s%entcells(I0:I1,J0:J1) )
      call ent_cell_nullify( s%entcells )

      do j=J0,J1
        do i=I0,I1
          if ( fearth(i,j) > 0.d0 )
     &         call ent_cell_construct( s%entcells(i,j) )
        enddo
      enddo

      !Allocate run-size arrays for hydrology
      allocate(
     &     s%w_ij(0:ngm,ls_nfrac,I0:I1,J0:J1),
     &     s%ht_ij(0:ngm,ls_nfrac,I0:I1,J0:J1),
     &     s%nsn_ij(     2,I0:I1,J0:J1),
     &     s%dzsn_ij(nlsn,2,I0:I1,J0:J1),
     &     s%wsn_ij(nlsn,2,I0:I1,J0:J1),
     &     s%hsn_ij(nlsn,2,I0:I1,J0:J1),
     &     s%fr_snow_ij(2,I0:I1,J0:J1),
     &     s%Qf_ij(I0:I1,J0:J1) )

      s%w_ij(:,:,:,:) = undef
      s%ht_ij(:,:,:,:) = undef
      s%nsn_ij(:,:,:) = UNDEFINT
      s%dzsn_ij(:,:,:,:) = undef
      s%wsn_ij(:,:,:,:) = undef
      s%hsn_ij(:,:,:,:) = undef
      s%fr_snow_ij(:,:,:) = undef
      s%Qf_ij(:,:) = undef

      end subroutine lsm_allocate_state
!-----------------------------------------------------------------------

      subroutine lsm_init_state_hydro( IM, JM, I0, I1, J0, J1
     &     ,i0s, i1s, j0s, j1s
     &     ,s)
      !lsm_init_state   Initialize surface hydrology state.
      use ghy_h, only : ngm, nlsn, LS_NFRAC
      integer, intent(in) :: IM, JM !Grid resolution
      integer, intent(in) :: I0, I1, J0, J1 !Run boundaries
      integer, intent(in) :: i0s, i1s, j0s, j1s !Hydrology state input file boundaries
      type(t_lsm_state) :: s
      !----
      integer i,j
!      integer sumj(jm)
!      real*8, dimension(N_COVERTYPES,I0:I1,J0:J1) :: laidata  !cohort
!      real*8, dimension(N_COVERTYPES,I0:I1,J0:J1) :: veg_height
      real*8:: !Input files (generally same as GCM)
     &     w_ij_read(0:ngm,ls_nfrac,i0s:i1s,j0s:j1s),
     &     ht_ij_read(0:ngm,ls_nfrac,i0s:i1s,j0s:j1s),
     &     nsn_ij_read(     2,i0s:i1s,j0s:j1s),
     &     dzsn_ij_read(nlsn,2,i0s:i1s,j0s:j1s),
     &     wsn_ij_read(nlsn,2,i0s:i1s,j0s:j1s),
     &     hsn_ij_read(nlsn,2,i0s:i1s,j0s:j1s),
     &     fr_snow_ij_read(2,i0s:i1s,j0s:j1s),
     &     Qf_ij_read(i0s:i1s,j0s:j1s) 


       !tbcs_lsm(I0:I1,J0:J1) = -1d30 ! initialize ground surface temperature[C]
      !tcanopy_lsm(I0:I1,J0:J1) = -1d30 !initialize canopy temperature[C]

!      print *,"sumj= ", sumj !##What's this for?

      !Read single to global size hydrology arrays.
      read(952)  !## Replace this with iu_state ##
     &       w_ij_read(0:ngm,1:ls_nfrac,:,:),         
     &       ht_ij_read(0:ngm,1:ls_nfrac,:,:),        
     &       nsn_ij_read    (1:2, :,:),              
     &       dzsn_ij_read   (1:nlsn, 1:2, :,:),      
     &       wsn_ij_read    (1:nlsn, 1:2, :,:),      
     &       hsn_ij_read    (1:nlsn, 1:2, :,:),      
     &       fr_snow_ij_read(1:2, :,:)
!     &      Qf_ij_read(:,:)  !##CHECK

      !Assign cells: extract run region from input files.
      s%w_ij(0:ngm,1:ls_nfrac,:,:)= 
     &     w_ij_read(0:ngm,1:ls_nfrac,I0:I1,J0:J1)
      s%ht_ij(0:ngm,1:ls_nfrac,:,:)=
     &     ht_ij_read(0:ngm,1:ls_nfrac,I0:I1,J0:J1)
      s%nsn_ij(1:2,:,: )= nsn_ij_read(1: 2,I0:I1,J0:J1)
      s%dzsn_ij(1:nlsn, 1:2,:,: )=
     &     dzsn_ij_read(1:nlsn,1:2,I0:I1,J0:J1)
      s%wsn_ij(1:nlsn, 1:2, :,:)=
     &     wsn_ij_read(1:nlsn,1:2,I0:I1,J0:J1)
      s%hsn_ij(1:nlsn, 1:2, :,:)=
     &     hsn_ij_read(1:nlsn,1:2,I0:I1,J0:J1)
      s%fr_snow_ij(1:2,:,: )=fr_snow_ij_read(1:2,I0:I1,J0:J1)
      s%Qf_ij(:,:)= Qf_ij_read(I0:I1,J0:J1) 

      ! something is wrong with snow - remove it for now ##
      s%nsn_ij    (1:2, :,:) = 0
      s%dzsn_ij   (1:nlsn, 1:2, :,:)= 0.d0
      s%wsn_ij    (1:nlsn, 1:2, :,:)= 0.d0
      s%hsn_ij    (1:nlsn, 1:2, :,:)= 0.d0
      s%fr_snow_ij(1:2, :,:) = 0.d0
      s%Qf_ij(:,:)=0.1d0

      end subroutine lsm_init_state_hydro
!======================================================================!

      subroutine lsm_set_bc_gcmdump(
     &     IM,JM,I0,I1,J0,J1,i0f,j0f,i1f,j1f, bc )
      !Old scheme for reading from fort.953 GCM dump
      integer,intent(in) :: IM,JM,I0,I1,J0,J1  
      integer,intent(in) :: i0f,j0f,i1f,j1f !Input file boundaries
      type (t_lsm_bc) :: bc
      !------
      integer :: iu_TOPO, iu_TOP_INDEX, iu_SOIL
      character*80 :: title
      real*4, allocatable :: temp_local(:,:,:)
      real*4, allocatable :: temp4(:,:)

      !Land arrays for input files.
      real*8 :: fearth_read (i0f:i1f,j0f:j1f),
     &     top_index_ij_read(i0f:i1f,j0f:j1f), 
     &     top_dev_ij_read  (i0f:i1f,j0f:j1f), 
     &     dz_ij_read       (i0f:i1f,j0f:j1f,1:ngm),
     &     q_ij_read        (i0f:i1f,j0f:j1f,1:imt,1:ngm),
     &     qk_ij_read       (i0f:i1f,j0f:j1f,1:imt,1:ngm),
     &     sl_ij_read       (i0f:i1f,j0f:j1f)


      !Run-size arrays (single to global size)
      allocate(
     &     bc%fearth      (I0:I1, J0:J1),
     &     bc%top_index_ij(I0:I1, J0:J1),
     &     bc%top_dev_ij  (I0:I1, J0:J1),
     &     bc%dz_ij       (I0:I1, J0:J1, ngm),
     &     bc%q_ij        (I0:I1, J0:J1, imt, ngm),
     &     bc%qk_ij       (I0:I1, J0:J1, imt, ngm),
     &     bc%sl_ij       (I0:I1, J0:J1) )


      !Read global land arrays
      !## Replace with iu_bc ##
      read(953)  
     &     fearth_read(:,:),
     &     top_index_ij_read(:,:),                 
     &     top_dev_ij_read(:,:),                   
     &     dz_ij_read(:,:,1:ngm),                   
     &     q_ij_read(:,:,1:imt,1:ngm),              
     &     qk_ij_read(:,:,1:imt,1:ngm),             
     &     sl_ij_read(:,:) 

      !Assign cells: extract run boundaries from input files.
      bc%fearth(:,:) = fearth_read(I0:I1,J0:J1)
      bc%top_index_ij(:,:) = top_index_ij_read(I0:I1,J0:J1)
      bc%top_dev_ij(:,:) = top_dev_ij_read(I0:I1,J0:J1)              
      bc%dz_ij(:,:,1:ngm) = dz_ij_read(I0:I1,J0:J1,1:ngm)
      bc%q_ij(:,:,1:imt,1:ngm) = 
     &     q_ij_read(I0:I1,J0:J1,1:imt,1:ngm)
      bc%qk_ij(:,:,1:imt,1:ngm) = 
     &     qk_ij_read(I0:I1,J0:J1,1:imt,1:ngm)
      bc%sl_ij(:,:) = sl_ij_read(I0:I1,J0:J1)

      end subroutine lsm_set_bc_gcmdump
!-----------------------------------------------------------------------

      subroutine lsm_set_bc(IM,JM,I0,I1,J0,J1,i0f,i1f,j0f,j1f, bc )
      !lsm_set_bc  Set up geological parameters:  topo, top_index, soil type.
      integer,intent(in) :: IM,JM,I0,I1,J0,J1 !Run boundaries.
      integer,intent(in) :: i0f,i1f,j0f,j1f !Input file bounds >= run bounds.
      type (t_lsm_bc) :: bc
      !------
      integer :: iu_TOPO, iu_TOP_INDEX, iu_SOIL
      integer :: IR, JR
      character*80 :: title
      real*4, allocatable :: temp_local(:,:,:)
      real*4, allocatable :: temp4(:,:)

      !Run-size arrays (single to global)
      allocate(
     &     bc%fearth      (I0:I1, J0:J1),
     &     bc%top_index_ij(I0:I1, J0:J1),
     &     bc%top_dev_ij  (I0:I1, J0:J1),
     &     bc%dz_ij       (I0:I1, J0:J1, ngm),
     &     bc%q_ij        (I0:I1, J0:J1, imt, ngm),
     &     bc%qk_ij       (I0:I1, J0:J1, imt, ngm),
     &     bc%sl_ij       (I0:I1, J0:J1) )

      !Read files.  Assign cells: extract run boundaries from input files.

      allocate( temp4(i0f:i1f,j0f:j1f) )
      call openunit("TOPO",iu_TOPO,.true.,.true.)
      read(iu_TOPO)
      read(iu_TOPO)
      read(iu_TOPO) title, temp4
      bc%fearth(I0:I1,J0:J1) = temp4(I0:I1, J0:J1)
      write(1220,*) temp4(I0:I1,J0:J1)
      write(1221,*) bc%fearth(I0:I1,J0:J1)
      call closeunit(iu_TOPO)

      call openunit("TOP_INDEX",iu_TOP_INDEX,.true.,.true.)
      read(iu_TOP_INDEX) title, temp4
      bc%top_index_ij(:,:) = temp4(I0:I1, J0:J1)
      read(iu_TOP_INDEX) title, temp4
      bc%top_dev_ij(:,:) = temp4(I0:I1, J0:J1)
      call closeunit(iu_TOP_INDEX)
      deallocate( temp4 )

      IR = I1 - I0 + 1
      JR = J1 - J0 + 1
      call openunit("SOIL",iu_SOIL,.true.,.true.)
      allocate(temp_local(i0f:i1f,j0f:j1f,11*ngm+1))
      read(iu_SOIL) temp_local
      bc%dz_ij(:,:,:)   = temp_local(I0:I1, J0:J1,1:ngm)
      bc%q_ij(:,:,:,:) =
     &     reshape( temp_local(I0:I1,J0:J1,1+ngm:) ,
     &     (/IR,JR,imt,ngm/) )
      bc%qk_ij(:,:,:,:) =
     &       reshape( temp_local(I0:I1, J0:J1,1+ngm+ngm*imt:) ,
     &     (/IR,JR,imt,ngm/) )
      bc%sl_ij(:,:)  =
     &     temp_local(I0:I1, J0:J1,1+ngm+ngm*imt+ngm*imt)
      deallocate(temp_local)
      call closeunit (iu_SOIL)

      end subroutine lsm_set_bc

!-----------------------------------------------------------------------
! Routine for testing with GCM met outputs      
      subroutine get_forcings_GCMdump(
     &     I0,I1,J0,J1,
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
      integer,intent(in) :: I0,I1,J0,J1
      real*8
     &         force_Ca (I0:I1,J0:J1)             ,
     &         force_cos_zen_angle (I0:I1,J0:J1)  ,
     &         force_vis_rad (I0:I1,J0:J1)        ,
     &         force_direct_vis_rad (I0:I1,J0:J1) ,
     &         force_prec_ms (I0:I1,J0:J1)        ,
     &         force_eprec_w (I0:I1,J0:J1)        ,
     &         force_sprec_ms (I0:I1,J0:J1)       ,
     &         force_seprec_w (I0:I1,J0:J1)       ,
     &         force_srheat (I0:I1,J0:J1)         ,
     &         force_trheat (I0:I1,J0:J1)         ,
     &         force_ts (I0:I1,J0:J1)             ,
     &         force_qs (I0:I1,J0:J1)             ,
     &         force_ps (I0:I1,J0:J1)             ,
     &         force_rhosrf (I0:I1,J0:J1)         ,
     &         force_cdh (I0:I1,J0:J1)            ,
     &         force_qm1 (I0:I1,J0:J1)            ,
     &         force_ws (I0:I1,J0:J1)             ,
     &         force_pbl_args_ws0 (I0:I1,J0:J1)   ,
     &         force_pbl_args_tprime (I0:I1,J0:J1),
     &         force_pbl_args_qprime (I0:I1,J0:J1)

      print *,"reading forcings"

      read(951)
     &         force_Ca (I0:I1,J0:J1)             ,
     &         force_cos_zen_angle (I0:I1,J0:J1)  ,
     &         force_vis_rad (I0:I1,J0:J1)        ,
     &         force_direct_vis_rad (I0:I1,J0:J1) ,
     &         force_prec_ms (I0:I1,J0:J1)        ,
     &         force_eprec_w (I0:I1,J0:J1)        ,
     &         force_sprec_ms (I0:I1,J0:J1)       ,
     &         force_seprec_w (I0:I1,J0:J1)       ,
     &         force_srheat (I0:I1,J0:J1)         ,
     &         force_trheat (I0:I1,J0:J1)         ,
     &         force_ts (I0:I1,J0:J1)             ,
     &         force_qs (I0:I1,J0:J1)             ,
     &         force_ps (I0:I1,J0:J1)             ,
     &         force_rhosrf (I0:I1,J0:J1)         ,
     &         force_cdh (I0:I1,J0:J1)            ,
     &         force_qm1 (I0:I1,J0:J1)            ,
     &         force_ws (I0:I1,J0:J1)             ,
     &         force_pbl_args_ws0 (I0:I1,J0:J1)   ,
     &         force_pbl_args_tprime (I0:I1,J0:J1),
     &         force_pbl_args_qprime (I0:I1,J0:J1)

      print *,"GCM test forcings ok"

      end subroutine get_forcings_GCMdump

!-----------------------------------------------------------------------

      subroutine lsm_set_config( dt_in
     &      ,force_VEG
     &     ,do_soilinit, do_soilresp
     &     ,do_phenology_activegrowth, do_structuralgrowth
     &     ,do_frost_hardiness, do_patchdynamics,do_init_geo,do_spinup)
!Set up coupled model configurations.

      use sle001, only : hl0, dt

      real*8 :: dt_in
!      integer :: n_month
!     Logical flags for the Ent DGTEM
      logical :: force_VEG
      logical :: do_soilinit, do_soilresp
      logical :: do_phenology_activegrowth, do_structuralgrowth
      logical :: do_frost_hardiness, do_patchdynamics, do_init_geo
      logical :: do_spinup

      dt = dt_in           ! intialize timestep size [sec]

      call ent_init_config(
     &     do_soilresp=do_soilresp
     &     ,do_phenology_activegrowth=do_phenology_activegrowth
     &     ,do_structuralgrowth=do_structuralgrowth
     &     ,do_frost_hardiness=do_frost_hardiness
     &     ,do_patchdynamics=do_patchdynamics
     &     ,do_init_geo=do_init_geo)

      call hl0

#ifdef PRINT_DIAGS
      print *, "Zeroing diags"
      call zero_diags()!Diag vars are global to module.
#endif

      end subroutine lsm_set_config

!-----------------------------------------------------------------------

      subroutine lsm_run(s,bc,jday,jyear,tyr_sec,dtsec,
     &     I0,I1,J0,J1,              !k_mnth
     &     tbcs_lsm,tcanopy_lsm    ,   
     &     force_Ca              ,
     &     force_cos_zen_angle   ,
     &     force_vis_rad         ,
     &     force_direct_vis_rad  ,
     &     force_prec_ms         ,
     &     force_eprec_w         ,
     &     force_sprec_ms        ,
     &     force_seprec_w        ,
     &     force_srheat          ,
     &     force_trheat          ,
     &     force_ts              ,
     &     force_qs              ,
     &     force_ps              ,
     &     force_rhosrf          ,
     &     force_cdh             ,
     &     force_qm1             ,
     &     force_ws              ,
     &     force_pbl_args_ws0    ,
     &     force_pbl_args_tprime ,
     &     force_pbl_args_qprime ,
     &     end_of_day_flag)

      use lsm_phys_util, only: month
      use drv_met, only : get_forcings_met !get_gswp_forcings
      use domain_decomp, only : mype
      use sle001, only : tp,tbcs

      implicit none
      type(t_lsm_state) :: s
      type (t_lsm_bc) :: bc

      integer, intent(in) :: jday, jyear, tyr_sec
      real*8, intent(in) :: dtsec !ModelE time step is in real*8
      integer, intent(in) :: I0,I1,J0,J1
      real*8,dimension(I0:I1,J0:J1), intent(inout) :: 
     &     tbcs_lsm,tcanopy_lsm
      real*8,dimension(I0:I1,J0:J1), intent(in) :: 
     &     force_Ca             ,
     &     force_cos_zen_angle  ,
     &     force_vis_rad        ,
     &     force_direct_vis_rad ,
     &     force_prec_ms        ,
     &     force_eprec_w        ,
     &     force_sprec_ms       ,
     &     force_seprec_w       ,
     &     force_srheat         ,
     &     force_trheat         ,
     &     force_ts             ,
     &     force_qs             ,
     &     force_ps             ,
     &     force_rhosrf         ,
     &     force_cdh            ,
     &     force_qm1            ,
     &     force_ws             ,
     &     force_pbl_args_ws0   ,
     &     force_pbl_args_tprime,
     &     force_pbl_args_qprime
      logical, intent(in) :: end_of_day_flag
      !----Local Variables------------------------
#ifdef PRINT_DIAGS      
      integer :: k_mnth !current month number
!      real*8, dimension(12) :: n_count_mnth
#endif      
      real*8 fb,fv
      integer i,j

      !call sysusage(mype+4,1)
      !call sysusage(mype+4,2)
      !call sysusage(mype+12,1)
      !call sysusage(mype+12,2)
      !call sysusage(mype+8,1)

#ifdef PRINT_DIAGS
      !call get_month(jyear,itime_3hr,k_mnth) ! month # from start of year
      k_mnth = month(jyear,tyr_sec)
      !Update global counter
      n_count_mnth(k_mnth) = n_count_mnth(k_mnth) + 1.d0
#endif


      ! Loop over run boundaries.
      loop_j: do j = J0,J1      ! do j = 1,jm
        loop_i: do i = I0,I1
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
          tbcs_lsm(i,j)=tbcs
          tcanopy_lsm(i,j)=tp(0,2)

#ifdef PRINT_DIAGS
!     Assign accumulators to output diagnostocs
      call accumulate_diags(i,j,k_mnth,s)
#endif

!          write(936,*) i,j, tbcs, s%w_ij(0:ngm,1:LS_NFRAC,i,j),
!     &     s%ht_ij(0:ngm,1:LS_NFRAC,i,j)

        enddo loop_i
      enddo loop_j

      !call sysusage(mype+8,2)
      end subroutine lsm_run

!-----------------------------------------------------------------------

      subroutine lsm_init_state_vegsoil(
     &     force_VEG,do_soilinit, do_phenology_activegrowth, do_init_geo
     &     ,jday, jyear
     &     ,IM,JM,I0,I1,J0,J1
     &     ,i0v, i1v, j0v, j1v         ,i0vh, i1vh, j0vh, j1vh
     &     ,i0vc, i1vc, j0vc, j1vc     ,i0cc, i1cc, j0cc, j1cc
     &     ,i0sc, i1sc,j0sc, j1sc      ,i0st, i1st,j0st, j1st
     &     ,fearth, entcells)
!@sum See Ent_standalone ent_prog.f/subroutine ent_init_vegstruct.
!@+   Initialization of Ent cells:  vegetation structure, soil type and texture.
!@+   Must be called for all run types.  If force_VEG, then first call to
!@+   read veg forcings will re-set LAI, height, and Clabile.
!@+   Halo cells ignored, i.e. entcells should be a slice without halo
      use ent_prescribed_drv, only : init_canopy_physical,
     &     prescr_vegdata, ent_init_params
       !use ent_prescr_veg, only : prescr_calc_shc
      !use ent_prescr_veg, only : prescr_calcconst
      implicit none
      logical, intent(in) :: force_VEG, do_soilinit
     &     ,do_phenology_activegrowth, do_init_geo
      integer, intent(in) :: jday, jyear
      integer, intent(in) :: IM,JM,I0,I1,J0,J1
     &     ,i0v, i1v, j0v, j1v         ,i0vh, i1vh, j0vh, j1vh
     &     ,i0vc, i1vc, j0vc, j1vc     ,i0cc, i1cc, j0cc, j1cc
     &     ,i0sc, i1sc,j0sc, j1sc      ,i0st, i1st,j0st, j1st 
      real*8 :: fearth(I0:I1,J0:J1)
      type(entcelltype_public), intent(inout) :: entcells(I0:I1,J0:J1)
      !--------Local variables-------
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
      integer i,j,k
      real*8 heat_capacity
      character*15 num
      
      !For diagnostics of initialization - NK
      character*80 TITLE
      real*4 :: LAYER(I0:I1,J0:J1) 
      real*4 :: LAYER2(I0:I1,J0:J1)
      real*4 :: LAYER3(I0:I1,J0:J1)
      real*4 :: LSOIL(N_DEPTH,I0:I1,J0:J1)
      real*8 :: elai, ebeta(N_DEPTH)

      !call prescr_calcconst()
      call ent_init_params()

!#ifndef MIXED_CANOPY
      print *,'Calling prescr_vegdata',do_init_geo
      call prescr_vegdata(jday, jyear, 
     &     IM,JM,I0,I1,J0,J1,
! LATER: Need to introduce file bounds for generic file reading. ##-NK
!     &     I0v,I1v,J0v,J1v,      I0vh,I1vh,J0vh,J1vh,
!     &     I0vc,I1vc,J0vc,J1vc,  I0cc,I1cc,J0cc,J1cc,
!     &     I0sc, I1sc,J0sc, J1sc,I0st, I1st,J0st, J1st, 
     &     vegdata,albedodata,laidata,hdata,nmdata,
     &     popdata,dbhdata,craddata,cpooldata,rootprofdata,
     &     soildata,soil_texture,Tpool_ini,
     &     do_soilinit,do_phenology_activegrowth, do_init_geo
     &     ,do_read_from_files=.true.)

      print *, 'Writing input diagnostics.'
      open(2001,file='gvegdata.ijk',form="unformatted",status="unknown")
      open(2002,file='ghdata.ijk',form="unformatted",status="unknown")
      open(2003,file='glaidata.ijk',form="unformatted",status="unknown")
      open(2004,file='gpopdata.ijk',form="unformatted",status="unknown")
      open(2005,file='gdbhdata.ijk',form="unformatted",status="unknown")

      LAYER(:,:) = 0.0
      LAYER2(:,:) = 0.0
      LAYER3(:,:) = 0.0
      do k=1,N_COVERTYPES
         write(num,*) k
         TITLE = "COVER "//trim(num)
         LAYER(:,:) = vegdata(k,:,:)
         LAYER3(:,:) = LAYER3(:,:) + vegdata(k,:,:)
         write(2001) TITLE, LAYER
         TITLE = "HEIGHT "//trim(num)
         LAYER(:,:) = hdata(k,:,:)
         write(2002) TITLE, LAYER
         TITLE = "LAI "//trim(num)
         LAYER(:,:) = laidata(k,:,:)
         write(2003) TITLE, LAYER
         LAYER2 = LAYER2 + laidata(k,:,:)*vegdata(k,:,:)
         TITLE = "POPDENS "//trim(num)
         LAYER(:,:) = popdata(k,:,:)
         write(2004) TITLE, LAYER
         TITLE = "DBH "//trim(num)
         LAYER(:,:) = dbhdata(k,:,:)
         write(2005) TITLE, LAYER
      enddo
      do i=I0,I1
         do j=j0,j1
            if (LAYER3(i,j).gt.0.) then
               LAYER2(i,j) = LAYER2(i,j) / LAYER3(i,j)
            else
               if (LAYER2(i,j).gt.0) then
                  print *, "WARNING: LAI where all bare.",i,j
               else
                  LAYER2(i,j) = 0.0
               endif
            endif
         enddo
      enddo
      TITLE = "LAI total (incl. bare+sand)"
      write(2003) TITLE, LAYER2

      close(2001)
      close(2002)
      close(2004)
      close(2005)

      call init_canopy_physical(I0, I1, J0, J1,
     &     Ci_ini, CNC_ini, Tcan_ini, Qf_ini)

      !Translate gridded data to Entdata structure
      !GISS data:  a patch per vegetation cover fraction, one cohort per patch
      call ent_cell_set(entcells, vegdata, popdata, laidata,
     &     hdata, dbhdata, craddata, cpooldata, nmdata, rootprofdata, 
     &     soildata, albedodata, soil_texture,
     &     Ci_ini, CNC_ini, Tcan_ini, Qf_ini, Tpool_ini
     &     , reinitialize=.true.)

      print *,"Writing some diags."
      ! Loop over run boundaries.
      loop_j: do j = J0,J1      ! do j = 1,jm
        loop_i: do i = I0,I1
          elai=0.d0
          ebeta(:)=0.d0
          LAYER(i,j) = elai
          LSOIL(:,i,j) = ebeta(:)
          if ( fearth(i,j).le.0.d0 ) cycle loop_i
          if (i==35 .and. j==84) cycle loop_i ! not a land cell for GSWP2 data
          !print *,"ent_get_exports",i,j
          call ent_get_exports( entcells(i,j),
     &         beta_soil_layers=ebeta, leaf_area_index=elai
     &         )
          LAYER(i,j) = elai
          LSOIL(:,i,j) = ebeta(:)
        enddo loop_i
      enddo loop_j
      TITLE = "LAI ent_cell_set:"
      write(2003) TITLE, LAYER
      do k=1,N_DEPTH
         write(num,*) k
         TITLE = "beta soil "//trim(num)
         write(2003) TITLE, LSOIL(k,:,:)
      enddo
      !call flush(6)
      close(2003)



      ! just in case, do nothing, just set heat capacities
      ! Fixed - now heat capacites shc automatically updated by summarize_entcell -NK
!      call ent_prescribe_vegupdate(entcells)

!#else
! MIXED_CANOPY
!         !* Patches with mixed-community structure.*!
!         call ent_struct_construct(entcells, IM,JM,I0,I1,J0,J1)
!         write(*,*) 'Constructed ent struct...'
!         call ent_struct_initphys(cells,IM,JM,I0,I1,J0,J1,
!     &        do_soilinit, .true.)'
!#endif
 
      end subroutine lsm_init_state_vegsoil

!-----------------------------------------------------------------------
      subroutine prescribe_soilcarbon(do_soilinit
     &     ,IM,JM,I0,I1,J0,J1,entcells)
      ! Set soil carbon pools from file.
      use ent_prescribed_drv, only:  prescr_get_soilpools,
     &     prescr_get_soil_C_total
      logical, intent(in) :: do_soilinit
      integer, intent(in) :: IM,JM,I0,I1,J0,J1
      type(entcelltype_public), intent(inout) :: entcells(I0:I1,J0:J1)
      !-----
      real*8 :: soil_C_total(N_CASA_LAYERS,I0:I1,J0:J1)
      real*8, dimension(N_PFT,PTRACE,NPOOLS-NLIVE,N_CASA_LAYERS,
     &     I0:I1,J0:J1):: Tpools !g/m2
      integer :: i,j

      print *,'Resetting soil carbon pools'
      if (.not.do_soilinit) then
         Tpools(:,:,:,:,:,:) = 0.d0
      else
        call prescr_get_soil_C_total(IM,JM,I0,I1,J0,J1,soil_C_total)
        call prescr_get_soilpools(I0,I1,J0,J1,soil_C_total,Tpools)
      endif

      print *,'Calling ent_cell_set_soilcarbon'
      call ent_cell_set_soilcarbon( entcells,Tpools)

      print *,'Reset soil carbon pools'
      end subroutine prescribe_soilcarbon
!-----------------------------------------------------------------------
      subroutine lsm_init_forcings_veg(force_VEG
     &     ,do_phenology_activegrowth
     &     ,tyr,updt,vdt, vtstart
     &     ,VROWMAX_IN, N_PFT, hemi, jday, jyear
     &     ,I0,I1,J0,J1,I0f,I1f,J0f,J1f,entcells)
!@sum ent_init_forcings_veg.  Initialize vegetation structure from
!+    forcings files.

      use drv_veg, only :  allocate_iu_veg,open_forcings_veg
     &     ,init_interp_forcings_veg, close_forcings_veg
     &     ,rewind_files_veg
      use ent_mod,only : ent_prescribe_vegupdate
      implicit none
      logical,intent(in) :: force_VEG
     &     ,do_phenology_activegrowth
      integer,intent(in) :: tyr,updt,vdt, vtstart
      integer,intent(in) :: VROWMAX_IN, N_PFT
      integer,intent(in) :: hemi(I0:I1,J0:J1)
      integer,intent(in) :: jday, jyear
      integer,intent(in) :: I0,I1,J0,J1,I0f,I1f,J0f,J1f
      type(entcelltype_public),intent(inout) :: entcells(:,:)
      !-----
      real*8, dimension(:,:,:),pointer :: LAI  !cohort
      real*8, dimension(:,:,:),pointer :: height    !cohort
      
      if (force_VEG) then
         allocate(LAI(N_PFT,I0:I1,J0:J1), height(N_PFT,I0:I1,J0:J1))
         call allocate_iu_veg(N_PFT)
         call open_forcings_veg(VROWMAX_IN) !May want to skip header rows.
         call init_interp_forcings_veg(tyr,updt,vdt, vtstart
     &        ,VROWMAX_IN, N_PFT
     &        ,I0,I1,J0,J1,I0f,I1f,J0f,J1f
     &        ,LAI,height)
         call rewind_files_veg(vdt,vtstart)
         call ent_prescribe_vegupdate(entcells,hemi,jday,jyear,
     &        do_giss_phenology=
     &        (.not.do_phenology_activegrowth).and.(.not.force_VEG),
     &        do_giss_albedo=.true.,
     &        do_giss_lai=
     &        (.not.do_phenology_activegrowth).and.(.not.force_VEG),
     &        update_crops=.false.,
     &        laidata=LAI,      ! pft LAI
     &        hdata=height, init=.true. ) ! pft height
         deallocate(LAI, height)
      endif
      end subroutine lsm_init_forcings_veg
!-----------------------------------------------------------------------

      subroutine deallocate_LSM(bc, s)
      implicit none
      type (t_lsm_bc) :: bc
      type(t_lsm_state) :: s

      deallocate( bc%fearth 
     &     ,bc%top_index_ij
     &     ,bc%top_dev_ij 
     &     ,bc%dz_ij 
     &     ,bc%q_ij 
     &     ,bc%qk_ij
     &     ,bc%sl_ij )

      end subroutine deallocate_LSM
!-----------------------------------------------------------------------
!Original in Ent_standalone, without reading veg data.
!      subroutine update_vegetation_data( entcells,
!     &     im,jm,i0,i1,j0,j1, hemi,jday, year, laidata, hdata, init )
!!@sum read standard GISS vegetation BC's and pass them to Ent for
!!@+   initialization of Ent cells. Halo cells ignored, i.e.
!!@+   entcells should be a slice without halo
!      use ent_prescribed_drv, only:
!     &     prescr_get_laidata
!      type(entcelltype_public), intent(out) :: entcells(I0:I1,J0:J1)
!      integer, intent(in) :: im, jm, i0, i1, j0, j1, jday, year
!      
!      logical, intent(in) :: init
!      !-----Local---------
!      real*8, dimension(N_BANDS,N_COVERTYPES,I0:I1,J0:J1) :: albedodata !patch, NOTE:snow
!      real*8, dimension(N_COVERTYPES,I0:I1,J0:J1) :: laidata  !cohort
!      integer hemi(I0:I1,J0:J1)
!      integer i,j
!      real*8 , dimension(N_COVERTYPES):: LAI_buf
!
!!Update Ent exactly like in ent_prog
!      !## Maybe move, so don't need to calculate hemi every time step. ##
!            !* Set hemisphere flags.
!
!      if (force_VEG) then
!         call ent_prescribe_vegupdate(entcells,hemi,jday,year,
!     &        do_giss_phenology=
!     &        (.not.do_phenology_activegrowth).and.(.not.force_VEG),
!     &        do_giss_albedo=.true.,
!     &        do_giss_lai=
!     &        (.not.do_phenology_activegrowth).and.(.not.force_VEG),
!     &        update_crops=.false., 
!     &        laidata=laidata, hdata=hdata,init=init)
!       else !If not force_VEG, don't include laidata and hdata in parameter list.
!          call ent_prescribe_vegupdate(entcells,hemi,jday,year,
!     &         do_giss_phenology=
!     &         (.not.do_phenology_activegrowth).and.(.not.force_VEG),
!     &         do_giss_albedo=.true.,
!     &         do_giss_lai=
!     &         (.not.do_phenology_activegrowth).and.(.not.force_VEG),
!     &         update_crops=.false., init=init)
!       endif
!
!       end subroutine update_vegetation_data
!-----------------------------------------------------------------------
      
      subroutine update_veg(
      !Update vegetation, from forcings files or from prescribed calculations.
      !Call get_interp_forcings_veg, 
      !then follow exact same calls as Ent_standalone ent_prog.
     i     hemi,jday,jyear,tyr,
     i     IM,JM,I0,I1,J0,J1,I0v,I1v,J0v,J1v,
     i     I0vh,I1vh,J0vh,J1vh,I0vc,I1vc,J0vc,J1vc,I0cc,I1cc,J0cc,J1cc,
     i     force_VEG,do_phenology_activegrowth,
     i     update_crops,init,
     i     updt,vdt,
     o     entcells)
      use drv_veg
      integer, intent(in) :: hemi(I0:I1,J0:J1)
      integer, intent(in) :: jday,jyear,tyr
      integer, intent(in) :: IM,JM,I0,I1,J0,J1,I0v,I1v,J0v,J1v
      integer, intent(in) :: I0vh,I1vh,J0vh,J1vh,I0vc,I1vc,J0vc,J1vc
     i     ,I0cc,I1cc,J0cc,J1cc
      logical, intent(in) :: force_VEG
     i     ,do_phenology_activegrowth, update_crops,init
      integer, intent(in) :: updt, vdt !seconds, update and row interval for veg
      type(entcelltype_public), intent(inout) :: entcells(:,:)
      !----Local-----
      real*8, dimension(:,:,:),pointer :: LAI, height

!     if (update_VEG_interval) 
!     --Assume outside check performed for appropriate update interval.
         if (force_VEG) then
            call get_interp_forcings_veg(tyr,updt,vdt
     &           ,I0,I1,J0,J1,I0v,I1v,J0v,J1v, LAI,height)

            call ent_prescribe_vegupdate(entcells,hemi,jday,jyear,
     &           do_giss_phenology=
     &           (.not.do_phenology_activegrowth).and.(.not.force_VEG),
     &           do_giss_albedo=.true.,
     &           do_giss_lai=
     &           (.not.do_phenology_activegrowth).and.(.not.force_VEG),
     &           update_crops=.false.,
     &           laidata=LAI,
     &           hdata=height,
     &           init=init ) ! pft height
            !print *,'ent_prescribe_vegupdate with LAI and height'
         else   !If not force_VEG, don't include laidata and hdata in parameter list.
            call ent_prescribe_vegupdate(entcells,hemi,jday,jyear,
     &           do_giss_phenology=
     &           (.not.do_phenology_activegrowth).and.(.not.force_VEG),
     &           do_giss_albedo=.true.,
     &           do_giss_lai=
     &           (.not.do_phenology_activegrowth).and.(.not.force_VEG),
     &           update_crops=.false.,init=init)
            !print *,'ent_prescribe_vegupdate with calculated lai, ht'
         endif
      !endif                     !update_VEG_interval


!#ifdef LSM_DRV_SINGLE
!      if (force_VEG) 
!     & call update_FLUXNET_LAI(iu_LAI,iu_vht,s%entcells(:,:)
!     &       , im, jm, i1f, j1f, jday, jyear, laidata, veg_height)
!    #endif
!      ! update vegetation only once per day
!      if ( jday /= jday_old) then
!#ifdef LSM_DRV_SINGLE
!        call update_veg_data_single( s%entcells(:,:)
!     &       , im, jm, i1f, j1f, jday, jyear, laidata, veg_height)
!#else ! GLOBAL_SCALE
!        call update_vegetation_data( s%entcells(:,j0:j1),
!     &       im, jm, 1, im, j0, j1, jday, jyear,init )
!#endif
      end subroutine update_veg
!-----------------------------------------------------------------------

!      subroutine update_FLUXNET_LAI( iu_LAI,iu_vht,entcells
!     &     , im, jm, I0,I1,J0,J1, jday, year,laidata,veg_height )
!!@sum read FLUXNET LAI data and pass to ent
!      use ent_pfts
!      integer :: iu_LAI, iu_vht
!      type(entcelltype_public), intent(inout) :: entcells(I0:I1,J0:J1)
!      integer, intent(in) :: im, jm,I0,I1,J0,J1
!      integer, intent(in) :: jday, year
!      real*8, dimension(N_COVERTYPES,1,1) :: laidata  !cohort
!      real*8, dimension(N_COVERTYPES,1,1):: veg_height
!      !-----Local---------
!      integer hemi(1,1)
!      integer i,j
!      real*8 , dimension(N_COVERTYPES):: LAI_buf, vht_buf
!
!!!! HACK : trying to update Ent exactly like in ent_prog
!      !* Set hemisphere flags.
!      if ( jsite<=JM/2 )   hemi(:,J0:min(JM/2,J1)) = -1    ! SOUTH
!      if ( jsite>=JM/2+1 ) hemi(:,max(JM/2+1,J0)) =  1    ! NORTH
!
!     Read Leaf Area Index data
!      read(iu_LAI,*) LAI_buf
!
!     Read vegetation height data      
!!      rewind(iu_vht) !KIM- height is the timeseries now!!!
!      read(iu_vht,*) vht_buf(:)
!
!      do i=1,N_COVERTYPES
!         laidata(i,1,1) = LAI_buf(i)
!         veg_height(i,1,1) = vht_buf(i)
!      end do
!
!      call ent_prescribe_vegupdate(entcells,hemi,jday,year,
!     &         do_giss_phenology=
!     &         (.not.do_phenology_activegrowth).and.(.not.force_VEG),
!     &         do_giss_albedo=.true.,
!     &         do_giss_lai=
!     &         (.not.do_phenology_activegrowth).and.(.not.force_VEG),
!     &         update_crops=.false.
!     &        ,laidata=laidata(COVEROFFSET+1:COVEROFFSET+N_PFT,:,:)
!     &        ,hdata=veg_height(COVEROFFSET+1:COVEROFFSET+N_PFT,:,:))
!
!      !print *, 'LAI=',laidata(6,1,1)
!
!      end subroutine update_FILE_LAI
!?#endif
!-----------------------------------------------------------------------

!=======================================================================
#ifdef PRINT_DIAGS
      !/*  Monthly diagnostic accumulators */
      subroutine allocate_diags(IM,JM,I0,I1,J0,J1)
      !Allocate global diagnostics arrays.
      integer,intent(in) :: IM,JM,I0,I1,J0,J1
      !---- Local ---------
      integer :: JA, JB
 
      JA=J0
      JB=J1
#ifdef USE_ESMF
      JA=1
      JB=JM
#endif

      allocate( aevap_mnth(I0:I1,JA:JB,12) )
      allocate( aevapw_mnth(I0:I1,JA:JB,12) )
      allocate( aevapd_mnth(I0:I1,JA:JB,12) )
      allocate( aevapb_mnth(I0:I1,JA:JB,12) )
      allocate( aintercep_mnth(I0:I1,JA:JB,12) )
      allocate( aruns_mnth(I0:I1,JA:JB,12) )
      allocate( arunu_mnth(I0:I1,JA:JB,12) )
      allocate( agpp_mnth(I0:I1,JA:JB,12) )
      allocate( arauto_mnth(I0:I1,JA:JB,12) )
      allocate( asoilresp_mnth(I0:I1,JA:JB,12) )
      allocate( aevapvg_mnth(I0:I1,JA:JB,12) )
      allocate( aevapvs_mnth(I0:I1,JA:JB,12) )
      allocate( aevapbs_mnth(I0:I1,JA:JB,12) )
      allocate( asrht_mnth(I0:I1,JA:JB,12) )
      allocate( atrht_mnth(I0:I1,JA:JB,12) )
      allocate( aalbedo_mnth(I0:I1,JA:JB,12) )
      allocate( aclab_mnth(I0:I1,JA:JB,12) )
      allocate( asoilCpoolsum_mnth(I0:I1,JA:JB,12) )
      allocate( aepp_mnth(I0:I1,JA:JB,12) )
      allocate( atrg_mnth(I0:I1,JA:JB,12) )
      allocate( ashg_mnth(I0:I1,JA:JB,12) )
      allocate( alhg_mnth(I0:I1,JA:JB,12) )
      allocate( aepc_mnth(I0:I1,JA:JB,12) )
      allocate( w_b1_mnth(I0:I1,JA:JB,12) )
      allocate( w_b2_mnth(I0:I1,JA:JB,12) )
      allocate( w_b3_mnth(I0:I1,JA:JB,12) )
      allocate( w_b4_mnth(I0:I1,JA:JB,12) )
      allocate( w_b5_mnth(I0:I1,JA:JB,12) )
      allocate( w_b6_mnth(I0:I1,JA:JB,12) )
      allocate( w_v1_mnth(I0:I1,JA:JB,12) )
      allocate( w_v2_mnth(I0:I1,JA:JB,12) )
      allocate( w_v3_mnth(I0:I1,JA:JB,12) )
      allocate( w_v4_mnth(I0:I1,JA:JB,12) )
      allocate( w_v5_mnth(I0:I1,JA:JB,12) )
      allocate( w_v6_mnth(I0:I1,JA:JB,12) )
      allocate( ht_b1_mnth(I0:I1,JA:JB,12) )
      allocate( ht_b2_mnth(I0:I1,JA:JB,12) )
      allocate( ht_b3_mnth(I0:I1,JA:JB,12) )
      allocate( ht_b4_mnth(I0:I1,JA:JB,12) )
      allocate( ht_b5_mnth(I0:I1,JA:JB,12) )
      allocate( ht_b6_mnth(I0:I1,JA:JB,12) )
      allocate( ht_v1_mnth(I0:I1,JA:JB,12) )
      allocate( ht_v2_mnth(I0:I1,JA:JB,12) )
      allocate( ht_v3_mnth(I0:I1,JA:JB,12) )
      allocate( ht_v4_mnth(I0:I1,JA:JB,12) )
      allocate( ht_v5_mnth(I0:I1,JA:JB,12) )
      allocate( ht_v6_mnth(I0:I1,JA:JB,12) )
      allocate( dzsn_b1(I0:I1,JA:JB,12) )
      allocate( dzsn_b2(I0:I1,JA:JB,12) )
      allocate( dzsn_b3(I0:I1,JA:JB,12) )
      allocate( dzsn_v1(I0:I1,JA:JB,12) )
      allocate( dzsn_v2(I0:I1,JA:JB,12) )
      allocate( dzsn_v3(I0:I1,JA:JB,12) )
      allocate( wsn_b1 (I0:I1,JA:JB,12) )
      allocate( wsn_b2 (I0:I1,JA:JB,12) )
      allocate( wsn_b3 (I0:I1,JA:JB,12) )
      allocate( wsn_v1 (I0:I1,JA:JB,12) )
      allocate( wsn_v2 (I0:I1,JA:JB,12) )
      allocate( wsn_v3(I0:I1,JA:JB,12) )
      allocate( hsn_b1 (I0:I1,JA:JB,12) )
      allocate( hsn_b2 (I0:I1,JA:JB,12) )
      allocate( hsn_b3 (I0:I1,JA:JB,12) )
      allocate( hsn_v1 (I0:I1,JA:JB,12) )
      allocate( hsn_v2 (I0:I1,JA:JB,12) )
      allocate( hsn_v3(I0:I1,JA:JB,12) )
      allocate( nsn_b(I0:I1,JA:JB,12) )
      allocate( nsn_v(I0:I1,JA:JB,12) )
      allocate( frsnow_b(I0:I1,JA:JB,12) )
      allocate( frsnow_v(I0:I1,JA:JB,12) )
      allocate( Qf_mnth(I0:I1,JA:JB,12) )
      allocate( abetad_mnth(I0:I1,JA:JB,12) )
      allocate( aClivepool_leaf_m(I0:I1,JA:JB,12) )
      allocate( aClivepool_froot_m(I0:I1,JA:JB,12) )
      allocate( aClivepool_wood_m(I0:I1,JA:JB,12) )
      allocate( aCdeadpool_surfmet_m(I0:I1,JA:JB,12) )
      allocate( aCdeadpool_surfstr_m(I0:I1,JA:JB,12) )
      allocate( aCdeadpool_soilmet_m(I0:I1,JA:JB,12) )
      allocate( aCdeadpool_soilstr_m(I0:I1,JA:JB,12) )
      allocate( aCdeadpool_cwd_m(I0:I1,JA:JB,12) )
      allocate( aCdeadpool_surfmic_m(I0:I1,JA:JB,12) )
      allocate( aCdeadpool_soilmic_m(I0:I1,JA:JB,12) )
      allocate( aCdeadpool_slow_m(I0:I1,JA:JB,12) )
      allocate( aCdeadpool_passive_m(I0:I1,JA:JB,12) )
      allocate( alai_m(I0:I1,JA:JB,12) )
      allocate( canopyH2O_m(I0:I1,JA:JB,12) )
      allocate( canopyheat_m(I0:I1,JA:JB,12) )

      end subroutine allocate_diags
!-----------------------------------------------------------------------
      subroutine deallocate_diags()
      !Allocate global diagnostics arrays.

      deallocate( aevap_mnth )
      deallocate( aevapw_mnth )
      deallocate( aevapd_mnth )
      deallocate( aevapb_mnth )
      deallocate( aintercep_mnth )
      deallocate( aruns_mnth )
      deallocate( arunu_mnth )
      deallocate( agpp_mnth )
      deallocate( arauto_mnth )
      deallocate( asoilresp_mnth )
      deallocate( aevapvg_mnth )
      deallocate( aevapvs_mnth )
      deallocate( aevapbs_mnth )
      deallocate( asrht_mnth )
      deallocate( atrht_mnth )
      deallocate( aalbedo_mnth )
      deallocate( aclab_mnth )
      deallocate( asoilCpoolsum_mnth )
      deallocate( aepp_mnth )
      deallocate( atrg_mnth )
      deallocate( ashg_mnth )
      deallocate( alhg_mnth )
      deallocate( aepc_mnth )
      deallocate( w_b1_mnth )
      deallocate( w_b2_mnth )
      deallocate( w_b3_mnth )
      deallocate( w_b4_mnth )
      deallocate( w_b5_mnth )
      deallocate( w_b6_mnth )
      deallocate( w_v1_mnth )
      deallocate( w_v2_mnth )
      deallocate( w_v3_mnth )
      deallocate( w_v4_mnth )
      deallocate( w_v5_mnth )
      deallocate( w_v6_mnth )
      deallocate( ht_b1_mnth )
      deallocate( ht_b2_mnth )
      deallocate( ht_b3_mnth )
      deallocate( ht_b4_mnth )
      deallocate( ht_b5_mnth )
      deallocate( ht_b6_mnth )
      deallocate( ht_v1_mnth )
      deallocate( ht_v2_mnth )
      deallocate( ht_v3_mnth )
      deallocate( ht_v4_mnth )
      deallocate( ht_v5_mnth )
      deallocate( ht_v6_mnth )
      deallocate( dzsn_b1 )
      deallocate( dzsn_b2 )
      deallocate( dzsn_b3 )
      deallocate( dzsn_v1 )
      deallocate( dzsn_v2 )
      deallocate( dzsn_v3 )
      deallocate( wsn_b1  )
      deallocate( wsn_b2  )
      deallocate( wsn_b3  )
      deallocate( wsn_v1  )
      deallocate( wsn_v2  )
      deallocate( wsn_v3 )
      deallocate( hsn_b1  )
      deallocate( hsn_b2  )
      deallocate( hsn_b3  )
      deallocate( hsn_v1  )
      deallocate( hsn_v2  )
      deallocate( hsn_v3 )
      deallocate( nsn_b )
      deallocate( nsn_v )
      deallocate( frsnow_b )
      deallocate( frsnow_v )
      deallocate( Qf_mnth )
      deallocate( abetad_mnth )
      deallocate( aClivepool_leaf_m )
      deallocate( aClivepool_froot_m )
      deallocate( aClivepool_wood_m )
      deallocate( aCdeadpool_surfmet_m )
      deallocate( aCdeadpool_surfstr_m )
      deallocate( aCdeadpool_soilmet_m )
      deallocate( aCdeadpool_soilstr_m )
      deallocate( aCdeadpool_cwd_m )
      deallocate( aCdeadpool_surfmic_m )
      deallocate( aCdeadpool_soilmic_m )
      deallocate( aCdeadpool_slow_m )
      deallocate( aCdeadpool_passive_m )
      deallocate( alai_m )
      deallocate( canopyH2O_m )
      deallocate( canopyheat_m )

      end subroutine deallocate_diags
!-----------------------------------------------------------------------

      subroutine print_diags(logcomment,year)
      use domain_decomp, only : array_gather, mype
      use lsm_phys_util, only : mon
      implicit none
      !Global real*8,dimension(12),intent(in) :: n_count_mnth
      character*80,intent(in) :: logcomment
      integer,intent(in) :: year
      !---- Local ------
      character*80 :: title       ! title of diagnostic to print to file
      integer :: i_mnth      ! month counter
      character*4 :: yearstr
      character*9 :: yrmn !year month string, YYYY MMM

      if (mype==0) then
         print *, logcomment
      endif

      write(yearstr,'(i4)') year

#ifdef USE_ESMF
          do i_mnth=1,12
             !print *, 'n_count_mnth(i_mnth)',i_mnth,n_count_mnth(i_mnth)
            if (n_count_mnth(i_mnth)>0) then

             yrmn = yearstr//' '//mon(i_mnth)//' '

             title = yrmn//' Total evaporation [kg/m2/month]'
             call array_gather( aevap_mnth(:,:,i_mnth) )
             if ( mype == 0 ) write(980) title,
     &                        real(aevap_mnth(:,:,i_mnth),kind=4)
             title = yrmn//' Canopy evaporation (evapw)[kg/m2/month]'
             call array_gather( aevapw_mnth(:,:,i_mnth) )
             if ( mype == 0 ) write(981) title,
     &                        real(aevapw_mnth(:,:,i_mnth),kind=4)
             title = yrmn//' Transpiration (evapd)[kg/m2/month]'
             call array_gather( aevapd_mnth(:,:,i_mnth) )
             if ( mype == 0 ) write(982) title,
     &                         real(aevapd_mnth(:,:,i_mnth),kind=4)
             title = yrmn//' Bare soil evap. (evapb)[kg/m2/month]'
             call array_gather( aevapb_mnth(:,:,i_mnth) )
             if ( mype == 0 ) write(983) title,
     &                         real(aevapb_mnth(:,:,i_mnth),kind=4)
             title = yrmn//' Canopy interception'
             call array_gather( aintercep_mnth(:,:,i_mnth) )
             if ( mype == 0 ) write(984) title,
     &                        real(aintercep_mnth(:,:,i_mnth),kind=4)
             title = yrmn//' Surface runoff [kg/m2/mnth]'
             call array_gather( aruns_mnth(:,:,i_mnth) )
             if ( mype == 0 ) write(985) title,
     &                         real(aruns_mnth(:,:,i_mnth),kind=4)
             title = yrmn//' Subsurface flow out of soil[kg/m2/mnth]'
             call array_gather( arunu_mnth(:,:,i_mnth) )
             if ( mype == 0 ) write(986) title,
     &                         real(arunu_mnth(:,:,i_mnth),kind=4)
             title = yrmn//' Gross primary productivity[kgC/m2/mnth]'
             call array_gather( agpp_mnth(:,:,i_mnth) )
             if ( mype == 0 ) write(987) title,
     &                         real(agpp_mnth(:,:,i_mnth),kind=4)
             title = yrmn//' Autotrophic respiration[kgC/m2/mnth]'
             call array_gather( arauto_mnth(:,:,i_mnth) )
             if ( mype == 0 ) write(988) title,
     &                         real(arauto_mnth(:,:,i_mnth),kind=4)
             title = yrmn//' Soil respiration[kgC/m2/mnth]'
             call array_gather( asoilresp_mnth(:,:,i_mnth) )
             if ( mype == 0 ) write(989) title,
     &                         real(asoilresp_mnth(:,:,i_mnth),kind=4)
             title = yrmn//' Evap. from vegetated soil[kg/m2/month]'
             call array_gather( aevapvg_mnth(:,:,i_mnth) )
             if ( mype == 0 ) write(990) title,
     &                         real(aevapvg_mnth(:,:,i_mnth),kind=4)
             title = yrmn//' Evap. from snow on veg[kg/m2/month]'
             call array_gather( aevapvs_mnth(:,:,i_mnth) )
             if ( mype == 0 ) write(991) title,
     &                         real(aevapvs_mnth(:,:,i_mnth),kind=4)
             title = yrmn//' Evap. from snow on bare soil[kg/m2/mn]'
             call array_gather( aevapbs_mnth(:,:,i_mnth) )
             if ( mype == 0 ) write(992) title,
     &                         real(aevapbs_mnth(:,:,i_mnth),kind=4)
             title = yrmn//' Net shortwave radiation [W/m2]'
             call array_gather( asrht_mnth(:,:,i_mnth) )
             asrht_mnth(:,:,i_mnth) = 
     &                asrht_mnth(:,:,i_mnth)/n_count_mnth(i_mnth)
             if ( mype == 0 ) write(993) title,
     &                         real(asrht_mnth(:,:,i_mnth),kind=4)
             title = yrmn//' Net longwave radiation [W/m2]'
             call array_gather( atrht_mnth(:,:,i_mnth) )
              atrht_mnth(:,:,i_mnth) = 
     &                atrht_mnth(:,:,i_mnth)/n_count_mnth(i_mnth)
             if ( mype == 0 ) write(994) title,
     &                         real(atrht_mnth(:,:,i_mnth),kind=4)
             title = yrmn//' Grid-cell mean albedo'
              aalbedo_mnth(:,:,i_mnth) = 
     &                aalbedo_mnth(:,:,i_mnth)/n_count_mnth(i_mnth)
             call array_gather( aalbedo_mnth(:,:,i_mnth) )
             if ( mype == 0 ) write(995) title,
     &                          real(aalbedo_mnth(:,:,i_mnth),kind=4)
             title = yrmn//' Plant labile carbon [kgC/m2]'
             aclab_mnth(:,:,i_mnth) = 
     &                aclab_mnth(:,:,i_mnth)/n_count_mnth(i_mnth)
             call array_gather( aclab_mnth(:,:,i_mnth) )
             if ( mype == 0 ) write(996) title,
     &                         real(aclab_mnth(:,:,i_mnth),kind=4)
             title = yrmn//' Sum of soil carbon pools[kgC/m2]'
             asoilCpoolsum_mnth(:,:,i_mnth) = 
     &             asoilCpoolsum_mnth(:,:,i_mnth)/n_count_mnth(i_mnth)
             call array_gather( asoilCpoolsum_mnth(:,:,i_mnth) )
             if ( mype == 0 ) write(997) title,
     &                      real(asoilCpoolsum_mnth(:,:,i_mnth),kind=4)
             title = yrmn//' Penman potential evap. [kg/m2/month]'
             call array_gather( aepp_mnth(:,:,i_mnth) )
             if ( mype == 0 ) write(998) title,
     &                      real(aepp_mnth(:,:,i_mnth),kind=4)
             title = yrmn//' Thermal heat from ground [J/m2/mnth]'
             call array_gather( atrg_mnth(:,:,i_mnth) )
             if ( mype == 0 ) write(999) title,
     &                      real(atrg_mnth(:,:,i_mnth),kind=4)
             title = yrmn//' Sensible heat from ground [J/m2/mnth]'
             call array_gather( ashg_mnth(:,:,i_mnth) )
             if ( mype == 0 ) write(1000) title,
     &                      real(ashg_mnth(:,:,i_mnth),kind=4)
             title = yrmn//' Latent heat flux [W/m2]'
             call array_gather( alhg_mnth(:,:,i_mnth) )
             alhg_mnth(:,:,i_mnth) = 
     &                alhg_mnth(:,:,i_mnth)/n_count_mnth(i_mnth)
             if ( mype == 0 ) write(1001) title,
     &                      real(alhg_mnth(:,:,i_mnth),kind=4)
             title = yrmn//' Pot. evap. from canopy [kg/m2/mnth]'
             call array_gather( aepc_mnth(:,:,i_mnth) )
             if ( mype == 0 ) write(1002) title,
     &                      real(aepc_mnth(:,:,i_mnth),kind=4)

             title = yrmn//' Soil water; Bare layer 1 [m]'
             call array_gather( w_b1_mnth(:,:,i_mnth) )
             w_b1_mnth(:,:,i_mnth) = 
     &                w_b1_mnth(:,:,i_mnth)/n_count_mnth(i_mnth)
             if ( mype == 0 ) write(1003) title,
     &                      real(w_b1_mnth(:,:,i_mnth),kind=4)

             title = yrmn//' Soil water; Bare layer 2 [m]'
             call array_gather( w_b2_mnth(:,:,i_mnth) )
             w_b2_mnth(:,:,i_mnth) = 
     &                w_b2_mnth(:,:,i_mnth)/n_count_mnth(i_mnth)
             if ( mype == 0 ) write(1004) title,
     &                      real(w_b2_mnth(:,:,i_mnth),kind=4)

             title = yrmn//' Soil water; Bare layer 3 [m]'
             call array_gather( w_b3_mnth(:,:,i_mnth) )
             w_b3_mnth(:,:,i_mnth) = 
     &                w_b3_mnth(:,:,i_mnth)/n_count_mnth(i_mnth)
             if ( mype == 0 ) write(1005) title,
     &                      real(w_b3_mnth(:,:,i_mnth),kind=4)

             title = yrmn//' Soil water; Bare layer 4 [m]'
             call array_gather( w_b4_mnth(:,:,i_mnth) )
             w_b4_mnth(:,:,i_mnth) = 
     &                w_b4_mnth(:,:,i_mnth)/n_count_mnth(i_mnth)
             if ( mype == 0 ) write(1006) title,
     &                      real(w_b4_mnth(:,:,i_mnth),kind=4)

             title = yrmn//' Soil water; Bare layer 5 [m]'
             call array_gather( w_b5_mnth(:,:,i_mnth) )
             w_b5_mnth(:,:,i_mnth) = 
     &                w_b5_mnth(:,:,i_mnth)/n_count_mnth(i_mnth)
             if ( mype == 0 ) write(1007) title,
     &                      real(w_b5_mnth(:,:,i_mnth),kind=4)

             title = yrmn//' Soil water; Bare layer 6 [m]'
             call array_gather( w_b6_mnth(:,:,i_mnth) )
             w_b6_mnth(:,:,i_mnth) = 
     &                w_b6_mnth(:,:,i_mnth)/n_count_mnth(i_mnth)
             if ( mype == 0 ) write(1008) title,
     &                      real(w_b6_mnth(:,:,i_mnth),kind=4)

             title = yrmn//' Soil water; Veg layer 1 [m]'
             call array_gather( w_v1_mnth(:,:,i_mnth) )
             w_v1_mnth(:,:,i_mnth) = 
     &                w_v1_mnth(:,:,i_mnth)/n_count_mnth(i_mnth)
             if ( mype == 0 ) write(1009) title,
     &                      real(w_v1_mnth(:,:,i_mnth),kind=4)

             title = yrmn//' Soil water; Veg layer 2 [m]'
             call array_gather( w_v2_mnth(:,:,i_mnth) )
             w_v2_mnth(:,:,i_mnth) = 
     &                w_v2_mnth(:,:,i_mnth)/n_count_mnth(i_mnth)
             if ( mype == 0 ) write(1010) title,
     &                      real(w_v2_mnth(:,:,i_mnth),kind=4)

             title = yrmn//' Soil water; Veg layer 3 [m]'
             call array_gather( w_v3_mnth(:,:,i_mnth) )
             w_v3_mnth(:,:,i_mnth) = 
     &                w_v3_mnth(:,:,i_mnth)/n_count_mnth(i_mnth)
             if ( mype == 0 ) write(1011) title,
     &                      real(w_v3_mnth(:,:,i_mnth),kind=4)

             title = yrmn//' Soil water; Veg layer 4 [m]'
             call array_gather( w_v4_mnth(:,:,i_mnth) )
             w_v4_mnth(:,:,i_mnth) = 
     &                w_v4_mnth(:,:,i_mnth)/n_count_mnth(i_mnth)
             if ( mype == 0 ) write(1012) title,
     &                      real(w_v4_mnth(:,:,i_mnth),kind=4)

             title = yrmn//' Soil water; Veg layer 5 [m]'
             call array_gather( w_v5_mnth(:,:,i_mnth) )
             w_v5_mnth(:,:,i_mnth) = 
     &                w_v5_mnth(:,:,i_mnth)/n_count_mnth(i_mnth)
             if ( mype == 0 ) write(1013) title,
     &                      real(w_v5_mnth(:,:,i_mnth),kind=4)

             title = yrmn//' Soil water; Veg layer 6 [m]'
             call array_gather( w_v6_mnth(:,:,i_mnth) )
             w_v6_mnth(:,:,i_mnth) = 
     &                w_v6_mnth(:,:,i_mnth)/n_count_mnth(i_mnth)
             if ( mype == 0 ) write(1014) title,
     &                      real(w_v6_mnth(:,:,i_mnth),kind=4)

             title = yrmn//' Soil heat; Bare layer 1 [J/m2]'
             call array_gather( ht_b1_mnth(:,:,i_mnth) )
             ht_b1_mnth(:,:,i_mnth) = 
     &                ht_b1_mnth(:,:,i_mnth)/n_count_mnth(i_mnth)
             if ( mype == 0 ) write(1015) title,
     &                      real(ht_b1_mnth(:,:,i_mnth),kind=4)

             title = yrmn//' Soil heat; Bare layer 2 [J/m2]'
             call array_gather( ht_b2_mnth(:,:,i_mnth) )
             ht_b2_mnth(:,:,i_mnth) = 
     &                ht_b2_mnth(:,:,i_mnth)/n_count_mnth(i_mnth)
             if ( mype == 0 ) write(1016) title,
     &                      real(ht_b2_mnth(:,:,i_mnth),kind=4)

             title = yrmn//' Soil heat; Bare layer 3 [J/m2]'
             call array_gather( ht_b3_mnth(:,:,i_mnth) )
             ht_b3_mnth(:,:,i_mnth) = 
     &                ht_b3_mnth(:,:,i_mnth)/n_count_mnth(i_mnth)
             if ( mype == 0 ) write(1017) title,
     &                      real(ht_b3_mnth(:,:,i_mnth),kind=4)

             title = yrmn//' Soil heat; Bare layer 4 [J/m2]'
             call array_gather( ht_b4_mnth(:,:,i_mnth) )
             ht_b4_mnth(:,:,i_mnth) = 
     &                ht_b4_mnth(:,:,i_mnth)/n_count_mnth(i_mnth)
             if ( mype == 0 ) write(1018) title,
     &                      real(ht_b4_mnth(:,:,i_mnth),kind=4)

             title = yrmn//' Soil heat; Bare layer 5 [J/m2]'
             call array_gather( ht_b5_mnth(:,:,i_mnth) )
             ht_b5_mnth(:,:,i_mnth) = 
     &                ht_b5_mnth(:,:,i_mnth)/n_count_mnth(i_mnth)
             if ( mype == 0 ) write(1019) title,
     &                      real(ht_b5_mnth(:,:,i_mnth),kind=4)

             title = yrmn//' Soil heat; Bare layer 6 [J/m2]'
             call array_gather( ht_b6_mnth(:,:,i_mnth) )
             ht_b6_mnth(:,:,i_mnth) = 
     &                ht_b6_mnth(:,:,i_mnth)/n_count_mnth(i_mnth)
             if ( mype == 0 ) write(1020) title,
     &                      real(ht_b6_mnth(:,:,i_mnth),kind=4)

             title = yrmn//' Soil heat; Veg layer 1 [J/m2]'
             call array_gather( ht_v1_mnth(:,:,i_mnth) )
             ht_v1_mnth(:,:,i_mnth) = 
     &                ht_v1_mnth(:,:,i_mnth)/n_count_mnth(i_mnth)
             if ( mype == 0 ) write(1021) title,
     &                      real(ht_v1_mnth(:,:,i_mnth),kind=4)

             title = yrmn//' Soil heat; Veg layer 2 [J/m2]'
             call array_gather( ht_v2_mnth(:,:,i_mnth) )
             ht_v2_mnth(:,:,i_mnth) = 
     &                ht_v2_mnth(:,:,i_mnth)/n_count_mnth(i_mnth)
             if ( mype == 0 ) write(1022) title,
     &                      real(ht_v2_mnth(:,:,i_mnth),kind=4)

             title = yrmn//' Soil heat; Veg layer 3 [J/m2]'
             call array_gather( ht_v3_mnth(:,:,i_mnth) )
             ht_v3_mnth(:,:,i_mnth) = 
     &                ht_v3_mnth(:,:,i_mnth)/n_count_mnth(i_mnth)
             if ( mype == 0 ) write(1023) title,
     &                      real(ht_v3_mnth(:,:,i_mnth),kind=4)

             title = yrmn//' Soil heat; Veg layer 4 [J/m2]'
             call array_gather( ht_v4_mnth(:,:,i_mnth) )
             ht_v4_mnth(:,:,i_mnth) = 
     &                ht_v4_mnth(:,:,i_mnth)/n_count_mnth(i_mnth)
             if ( mype == 0 ) write(1024) title,
     &                      real(ht_v4_mnth(:,:,i_mnth),kind=4)

             title = yrmn//' Soil heat; Veg layer 5 [J/m2]'
             call array_gather( ht_v5_mnth(:,:,i_mnth) )
             ht_v5_mnth(:,:,i_mnth) = 
     &                ht_v5_mnth(:,:,i_mnth)/n_count_mnth(i_mnth)
             if ( mype == 0 ) write(1025) title,
     &                      real(ht_v5_mnth(:,:,i_mnth),kind=4)

             title = yrmn//' Soil heat; Veg layer 6 [J/m2]'
             call array_gather( ht_v6_mnth(:,:,i_mnth) )
             ht_v6_mnth(:,:,i_mnth) = 
     &                ht_v6_mnth(:,:,i_mnth)/n_count_mnth(i_mnth)
             if ( mype == 0 ) write(1026) title,
     &                      real(ht_v6_mnth(:,:,i_mnth),kind=4)

             title = yrmn//' Snow layer thickness; bare layer 1 [m]'
             call array_gather( dzsn_b1(:,:,i_mnth) )
             dzsn_b1(:,:,i_mnth) = 
     &                dzsn_b1(:,:,i_mnth)/n_count_mnth(i_mnth)
             if ( mype == 0 ) write(1027) title,
     &                      real(dzsn_b1(:,:,i_mnth),kind=4)

             title = yrmn//'  Snow layer thickness; bare layer 2 [m]'
             call array_gather( dzsn_b2(:,:,i_mnth) )
             dzsn_b2(:,:,i_mnth) = 
     &                dzsn_b2(:,:,i_mnth)/n_count_mnth(i_mnth)
             if ( mype == 0 ) write(1028) title,
     &                      real(dzsn_b2(:,:,i_mnth),kind=4)

             title = yrmn//' Snow layer thickness; bare layer 3 [m]'
             call array_gather( dzsn_b3(:,:,i_mnth) )
             dzsn_b3(:,:,i_mnth) = 
     &                dzsn_b3(:,:,i_mnth)/n_count_mnth(i_mnth)
             if ( mype == 0 ) write(1029) title,
     &                      real(dzsn_b3(:,:,i_mnth),kind=4)

             title = yrmn//' Snow layer thickness; bare layer 1 [m]'
             call array_gather( dzsn_v1(:,:,i_mnth) )
             dzsn_v1(:,:,i_mnth) = 
     &                dzsn_v1(:,:,i_mnth)/n_count_mnth(i_mnth)
             if ( mype == 0 ) write(1030) title,
     &                      real(dzsn_v1(:,:,i_mnth),kind=4)

             title = yrmn//' Snow layer thickness; bare layer 2 [m]'
             call array_gather( dzsn_v2(:,:,i_mnth) )
             dzsn_v2(:,:,i_mnth) = 
     &                dzsn_v2(:,:,i_mnth)/n_count_mnth(i_mnth)
             if ( mype == 0 ) write(1031) title,
     &                      real(dzsn_v2(:,:,i_mnth),kind=4)

             title = yrmn//' Snow layer thickness; bare layer 3 [m]'
             call array_gather( dzsn_v3(:,:,i_mnth) )
             dzsn_v3(:,:,i_mnth) = 
     &                dzsn_v3(:,:,i_mnth)/n_count_mnth(i_mnth)
             if ( mype == 0 ) write(1032) title,
     &                      real(dzsn_v3(:,:,i_mnth),kind=4)

             title = yrmn//' Snow water equival; bare layer 1 [m]'
             call array_gather( wsn_b1(:,:,i_mnth) )
             wsn_b1(:,:,i_mnth) = 
     &                wsn_b1(:,:,i_mnth)/n_count_mnth(i_mnth)
             if ( mype == 0 ) write(1033) title,
     &                      real(wsn_b1(:,:,i_mnth),kind=4)

             title = yrmn//' Snow water equival; bare layer 2 [m]'
             call array_gather( wsn_b2(:,:,i_mnth) )
             wsn_b2(:,:,i_mnth) = 
     &                wsn_b2(:,:,i_mnth)/n_count_mnth(i_mnth)
             if ( mype == 0 ) write(1034) title,
     &                      real(wsn_b2(:,:,i_mnth),kind=4)

             title = yrmn//' Snow water equival; bare layer 3 [m]'
             call array_gather( wsn_b3(:,:,i_mnth) )
             wsn_b3(:,:,i_mnth) = 
     &                wsn_b3(:,:,i_mnth)/n_count_mnth(i_mnth)
             if ( mype == 0 ) write(1035) title,
     &                      real(wsn_b3(:,:,i_mnth),kind=4)

             title = yrmn//' Snow water equival; veg layer 1 [m]'
             call array_gather( wsn_v1(:,:,i_mnth) )
             wsn_v1(:,:,i_mnth) = 
     &                wsn_v1(:,:,i_mnth)/n_count_mnth(i_mnth)
             if ( mype == 0 ) write(1036) title,
     &                      real(wsn_v1(:,:,i_mnth),kind=4)

             title = yrmn//' Snow water equival; veg layer 2 [m]'
             call array_gather( wsn_v2(:,:,i_mnth) )
             wsn_v2(:,:,i_mnth) = 
     &                wsn_v2(:,:,i_mnth)/n_count_mnth(i_mnth)
             if ( mype == 0 ) write(1037) title,
     &                      real(wsn_v2(:,:,i_mnth),kind=4)

             title = yrmn//' Snow water equival; veg layer 3 [m]'
             call array_gather( wsn_v3(:,:,i_mnth) )
             wsn_v3(:,:,i_mnth) = 
     &                wsn_v3(:,:,i_mnth)/n_count_mnth(i_mnth)
             if ( mype == 0 ) write(1038) title,
     &                      real(wsn_v3(:,:,i_mnth),kind=4)

             title = yrmn//' Snow layer heat; bare layer 1 [J/m2]'
             call array_gather( hsn_b1(:,:,i_mnth) )
             hsn_b1(:,:,i_mnth) = 
     &                hsn_b1(:,:,i_mnth)/n_count_mnth(i_mnth)
             if ( mype == 0 ) write(1039) title,
     &                      real(hsn_b1(:,:,i_mnth),kind=4)

             title = yrmn//' Snow layer heat; bare layer 2 [J/m2]'
             call array_gather( hsn_b2(:,:,i_mnth) )
             hsn_b2(:,:,i_mnth) = 
     &                hsn_b2(:,:,i_mnth)/n_count_mnth(i_mnth)
             if ( mype == 0 ) write(1040) title,
     &                      real(hsn_b2(:,:,i_mnth),kind=4)

             title = yrmn//'  Snow layer heat; bare layer 3 [J/m2]'
             call array_gather( hsn_b3(:,:,i_mnth) )
             hsn_b3(:,:,i_mnth) = 
     &                hsn_b3(:,:,i_mnth)/n_count_mnth(i_mnth)
             if ( mype == 0 ) write(1041) title,
     &                      real(hsn_b3(:,:,i_mnth),kind=4)

             title = yrmn//'  Snow layer heat; veg layer 1 [J/m2]'
             call array_gather( hsn_v1(:,:,i_mnth) )
             hsn_v1(:,:,i_mnth) = 
     &                hsn_v1(:,:,i_mnth)/n_count_mnth(i_mnth)
             if ( mype == 0 ) write(1042) title,
     &                      real(hsn_v1(:,:,i_mnth),kind=4)

             title = yrmn//'  Snow layer heat; veg layer 2 [J/m2]'
             call array_gather( hsn_v2(:,:,i_mnth) )
             hsn_v2(:,:,i_mnth) = 
     &                hsn_v2(:,:,i_mnth)/n_count_mnth(i_mnth)
             if ( mype == 0 ) write(1043) title,
     &                      real(hsn_v2(:,:,i_mnth),kind=4)

             title = yrmn//'  Snow layer heat; veg layer 3 [J/m2]'
             call array_gather( hsn_v3(:,:,i_mnth) )
             hsn_v3(:,:,i_mnth) = 
     &                hsn_v3(:,:,i_mnth)/n_count_mnth(i_mnth)
             if ( mype == 0 ) write(1044) title,
     &                      real(hsn_v3(:,:,i_mnth),kind=4)

             title = yrmn//' # of snow layers on bare land'
             call array_gather( nsn_b(:,:,i_mnth) )
             nsn_b(:,:,i_mnth) = 
     &                nsn_b(:,:,i_mnth)/n_count_mnth(i_mnth)
             if ( mype == 0 ) write(1045) title,
     &                      real(nsn_b(:,:,i_mnth),kind=4)

             title = yrmn//' # of snow layers on vegetated land'
             call array_gather( nsn_v(:,:,i_mnth) )
             nsn_v(:,:,i_mnth) = 
     &                nsn_v(:,:,i_mnth)/n_count_mnth(i_mnth)
             if ( mype == 0 ) write(1046) title,
     &                      real(nsn_v(:,:,i_mnth),kind=4)

             title = yrmn//' Snow-covered bare fraction [-]'
             call array_gather( frsnow_b(:,:,i_mnth) )
             frsnow_b(:,:,i_mnth) = 
     &                frsnow_b(:,:,i_mnth)/n_count_mnth(i_mnth)
             if ( mype == 0 ) write(1047) title,
     &                      real(frsnow_b(:,:,i_mnth),kind=4)

             title = yrmn//'  Snow-covered vegetated fraction [-]'
             call array_gather( frsnow_v(:,:,i_mnth) )
             frsnow_v(:,:,i_mnth) = 
     &                frsnow_v(:,:,i_mnth)/n_count_mnth(i_mnth)
             if ( mype == 0 ) write(1048) title,
     &                      real(frsnow_v(:,:,i_mnth),kind=4)

             title = yrmn//'Foliage surf. vapor mixing ratio [kg/kg]'
             call array_gather( Qf_mnth(:,:,i_mnth) )
             Qf_mnth(:,:,i_mnth) = 
     &                Qf_mnth(:,:,i_mnth)/n_count_mnth(i_mnth)
             if ( mype == 0 ) write(1049) title,
     &                      real(Qf_mnth(:,:,i_mnth),kind=4)

             title = yrmn//'root zone betad [-]'
             call array_gather( abetad_mnth(:,:,i_mnth) )
             abetad_mnth(:,:,i_mnth) = 
     &                abetad_mnth(:,:,i_mnth)/n_count_mnth(i_mnth)
             if ( mype == 0 ) write(1050) title,
     &                      real(abetad_mnth(:,:,i_mnth),kind=4)

             title = yrmn//'Live leaf carbon pool [kg/m2]'
             call array_gather( aClivepool_leaf_m(:,:,i_mnth) )
             aClivepool_leaf_m(:,:,i_mnth) = 
     &              aClivepool_leaf_m(:,:,i_mnth)/n_count_mnth(i_mnth)
             if ( mype == 0 ) write(1051) title,
     &                      real(aClivepool_leaf_m(:,:,i_mnth),kind=4)

             title = yrmn//'Live froot carbon pool [kg/m2]'
             call array_gather( aClivepool_froot_m(:,:,i_mnth) )
             aClivepool_froot_m(:,:,i_mnth) = 
     &              aClivepool_froot_m(:,:,i_mnth)/n_count_mnth(i_mnth)
             if ( mype == 0 ) write(1052) title,
     &                    real(aClivepool_froot_m(:,:,i_mnth),kind=4)

             title = yrmn//'Live wood carbon pool [kg/m2]'
             call array_gather( aClivepool_wood_m(:,:,i_mnth) )
             aClivepool_wood_m(:,:,i_mnth) = 
     &              aClivepool_wood_m(:,:,i_mnth)/n_count_mnth(i_mnth)
             if ( mype == 0 ) write(1053) title,
     &                    real(aClivepool_wood_m(:,:,i_mnth),kind=4)

             title = yrmn//'Dead surface metabolic C pool [kg/m2]'
             call array_gather( aCdeadpool_surfmet_m(:,:,i_mnth) )
             aCdeadpool_surfmet_m(:,:,i_mnth) = 
     &            aCdeadpool_surfmet_m(:,:,i_mnth)/n_count_mnth(i_mnth)
             if ( mype == 0 ) write(1054) title,
     &                    real(aCdeadpool_surfmet_m(:,:,i_mnth),kind=4)

             title = yrmn//'Dead surface structural C pool [kg/m2]'
             call array_gather(aCdeadpool_surfstr_m(:,:,i_mnth) )
             aCdeadpool_surfstr_m(:,:,i_mnth) = 
     &           aCdeadpool_surfstr_m(:,:,i_mnth)/n_count_mnth(i_mnth)
             if ( mype == 0 ) write(1055) title,
     &                    real(aCdeadpool_surfstr_m(:,:,i_mnth),kind=4)

             title = yrmn//'Dead soil metabolic C pool [kg/m2]'
             call array_gather(aCdeadpool_soilmet_m(:,:,i_mnth) )
             aCdeadpool_soilmet_m(:,:,i_mnth) = 
     &           aCdeadpool_soilmet_m(:,:,i_mnth)/n_count_mnth(i_mnth)
             if ( mype == 0 ) write(1056) title,
     &                    real(aCdeadpool_soilmet_m(:,:,i_mnth),kind=4)

             title = yrmn//'Dead soil structural C pool [kg/m2]'
             call array_gather(aCdeadpool_soilstr_m(:,:,i_mnth) )
             aCdeadpool_soilstr_m(:,:,i_mnth) = 
     &           aCdeadpool_soilstr_m(:,:,i_mnth)/n_count_mnth(i_mnth)
             if ( mype == 0 ) write(1057) title,
     &                    real(aCdeadpool_soilstr_m(:,:,i_mnth),kind=4)

             title = yrmn//'Coarse woody debris C pool [kg/m2]'
             call array_gather( aCdeadpool_cwd_m(:,:,i_mnth) )
             aCdeadpool_cwd_m(:,:,i_mnth) = 
     &            aCdeadpool_cwd_m(:,:,i_mnth)/n_count_mnth(i_mnth)
             if ( mype == 0 ) write(1058) title,
     &                    real(aCdeadpool_cwd_m(:,:,i_mnth),kind=4)

             title = yrmn//'Dead surface microbial C pool [kg/m2]'
             call array_gather(aCdeadpool_surfmic_m(:,:,i_mnth) )
             aCdeadpool_surfmic_m(:,:,i_mnth) = 
     &            aCdeadpool_surfmic_m(:,:,i_mnth)/n_count_mnth(i_mnth)
             if ( mype == 0 ) write(1059) title,
     &                    real(aCdeadpool_surfmic_m(:,:,i_mnth),kind=4)

             title = yrmn//'Dead soil microbial C pool [kg/m2]'
             call array_gather( aCdeadpool_soilmic_m(:,:,i_mnth) )
             aCdeadpool_soilmic_m(:,:,i_mnth) = 
     &            aCdeadpool_soilmic_m(:,:,i_mnth)/n_count_mnth(i_mnth)
             if ( mype == 0 ) write(1060) title,
     &                    real(aCdeadpool_soilmic_m(:,:,i_mnth),kind=4)

             title = yrmn//'Dead slow (~10 yr) C pool [kg/m2]'
             call array_gather(aCdeadpool_slow_m(:,:,i_mnth) )
             aCdeadpool_slow_m(:,:,i_mnth) = 
     &            aCdeadpool_slow_m(:,:,i_mnth)/n_count_mnth(i_mnth)
             if ( mype == 0 ) write(1061) title,
     &                    real(aCdeadpool_slow_m(:,:,i_mnth),kind=4)

             title = yrmn//'Dead passive (~100 yr) C pool [kg/m2]'
             call array_gather(aCdeadpool_passive_m(:,:,i_mnth) )
             aCdeadpool_passive_m(:,:,i_mnth) = 
     &            aCdeadpool_passive_m(:,:,i_mnth)/n_count_mnth(i_mnth)
             if ( mype == 0 ) write(1062) title,
     &                    real(aCdeadpool_passive_m(:,:,i_mnth),kind=4)

             title = yrmn//' Leaf Area Index [-]'
             call array_gather(alai_m(:,:,i_mnth) )
             alai_m(:,:,i_mnth)=alai_m(:,:,i_mnth)/n_count_mnth(i_mnth)
             if ( mype == 0 ) write(1063) title,
     &                    real(alai_m(:,:,i_mnth),kind=4)

             title = yrmn//' Water stored on the canopy [m]'
             call array_gather(canopyH2O_m(:,:,i_mnth) )
             canopyH2O_m(:,:,i_mnth)=canopyH2O_m(:,:,i_mnth)
     &                               /n_count_mnth(i_mnth)
             if ( mype == 0 ) write(1064) title,
     &                    real(canopyH2O_m(:,:,i_mnth),kind=4)

             title = yrmn//' Heat of the canopy [J/m2]'
             call array_gather(canopyheat_m(:,:,i_mnth) )
             canopyheat_m(:,:,i_mnth)=canopyheat_m(:,:,i_mnth)
     &                                /n_count_mnth(i_mnth)
             if ( mype == 0 ) write(1065) title,
     &                    real(canopyheat_m(:,:,i_mnth),kind=4)
           endif !if n_count_mnth>0
          enddo
#else
!     Annual monthly diagnostic accumulators
          do i_mnth=1,12
              yrmn = yearstr//' '//mon(i_mnth)//' '
              title = yrmn//' Total evaporation'
              write(980) title,real(aevap_mnth(:,:,i_mnth),kind=4)
              title = yrmn//' Canopy evaporation (evapw)'
              write(981) title,real(aevapw_mnth(:,:,i_mnth),kind=4)
              title = yrmn//' Transpiration (evapd)'
              write(982) title,real(aevapd_mnth(:,:,i_mnth),kind=4)
              title = yrmn//' Bare soil evaporation (evapb)'
              write(983) title,real(aevapb_mnth(:,:,i_mnth),kind=4)
              title = yrmn//' Canopy interception'
              write(984) title,real(aintercep_mnth(:,:,i_mnth),kind=4)
              title = yrmn//' Surface runoff (runs)'
              write(985) title,real(aruns_mnth(:,:,i_mnth),kind=4)
              title = yrmn//' Subsurface flow out of soil column'
              write(986) title,real(arunu_mnth(:,:,i_mnth),kind=4)
              title = yrmn//' Gross primary productivity '//
     &             '(1E-09 kgC/m2/s)'
              write(987) title,real(agpp_mnth(:,:,i_mnth),kind=4)
              title = yrmn//' Autotrophic respiration '//
     &             '(1E-09 kgC/m2/s)'
              write(988) title,real(arauto_mnth(:,:,i_mnth),kind=4)
              title = yrmn//' Soil respiration '//
     &             '(1E-09 kgC/m2/s)'
              write(989) title,real(asoilresp_mnth(:,:,i_mnth),kind=4)
              title = yrmn//' Evap. from vegetated soil'
              write(990) title,real(aevapvg_mnth(:,:,i_mnth),kind=4)
              title = yrmn//' Evap. from snow on vegetation'
              write(991) title,real(aevapvs_mnth(:,:,i_mnth),kind=4)
              title = yrmn//' Evap. from snow on bare soil'
              write(992) title,real(aevapbs_mnth(:,:,i_mnth),kind=4)
              title = yrmn//' Net shortwave radiation [W/m2]'
              asrht_mnth(:,:,i_mnth) = 
     &                asrht_mnth(:,:,i_mnth)/n_count_mnth(i_mnth)
              write(993) title,real(asrht_mnth(:,:,i_mnth),kind=4)
              title = yrmn//' Net longwave radiation [W/m2]'
              atrht_mnth(:,:,i_mnth) = 
     &                atrht_mnth(:,:,i_mnth)/n_count_mnth(i_mnth)
              write(994) title,real(atrht_mnth(:,:,i_mnth),kind=4)
              title = yrmn//' Grid-cell mean albedo'
              aalbedo_mnth(:,:,i_mnth) = 
     &                aalbedo_mnth(:,:,i_mnth)/n_count_mnth(i_mnth)
              write(995) title,real(aalbedo_mnth(:,:,i_mnth),kind=4)
              title = yrmn//' Plant labile carbon [kgC/m2]'
              aclab_mnth(:,:,i_mnth) = 
     &                aclab_mnth(:,:,i_mnth)/n_count_mnth(i_mnth)
              write(996) title,real(aclab_mnth(:,:,i_mnth),kind=4)
              title = yrmn//' Sum of soil carbon pools[kgC/m2]'
              asoilCpoolsum_mnth(:,:,i_mnth) = 
     &          asoilCpoolsum_mnth(:,:,i_mnth)/n_count_mnth(i_mnth)
              write(997) title,real(asoilCpoolsum_mnth(:,:,i_mnth)
     &                              ,kind=4)
              title = yrmn//' Penman potential evap. [kg/m2/month]'
              write(998) title,real(aepp_mnth(:,:,i_mnth),kind=4)
              title = yrmn//' Thermal heat from ground [J/m2/mnth]'
              write(999) title,real(atrg_mnth(:,:,i_mnth),kind=4)
              title = yrmn//' Sensible heat from ground [J/m2/mnth]'
              write(1000) title,real(ashg_mnth(:,:,i_mnth),kind=4)
              title = yrmn//' Latent heat flux [W/m2]'
              alhg_mnth(:,:,i_mnth) = 
     &                alhg_mnth(:,:,i_mnth)/n_count_mnth(i_mnth)
              write(1001) title,real(alhg_mnth(:,:,i_mnth),kind=4)
              title = yrmn//' Pot. evap. from canopy [kg/m2/mnth]'
              write(1002) title,real(aepc_mnth(:,:,i_mnth),kind=4)
             title = yrmn//' Soil water; Bare layer 1 [m]'
             w_b1_mnth(:,:,i_mnth) = 
     &                w_b1_mnth(:,:,i_mnth)/n_count_mnth(i_mnth)
             write(1003) title,real(w_b1_mnth(:,:,i_mnth),kind=4)
             title = yrmn//' Soil water; Bare layer 2 [m]'
             w_b2_mnth(:,:,i_mnth) = 
     &                w_b2_mnth(:,:,i_mnth)/n_count_mnth(i_mnth)
             write(1004) title,real(w_b2_mnth(:,:,i_mnth),kind=4)
             title = yrmn//' Soil water; Bare layer 3 [m]'
             w_b3_mnth(:,:,i_mnth) = 
     &                w_b3_mnth(:,:,i_mnth)/n_count_mnth(i_mnth)
             write(1005) title,real(w_b3_mnth(:,:,i_mnth),kind=4)
             title = yrmn//' Soil water; Bare layer 4 [m]'
             w_b4_mnth(:,:,i_mnth) = 
     &                w_b4_mnth(:,:,i_mnth)/n_count_mnth(i_mnth)
             write(1006) title,real(w_b4_mnth(:,:,i_mnth),kind=4)
             title = yrmn//' Soil water; Bare layer 5 [m]'
             w_b5_mnth(:,:,i_mnth) = 
     &                w_b5_mnth(:,:,i_mnth)/n_count_mnth(i_mnth)
             write(1007) title,real(w_b5_mnth(:,:,i_mnth),kind=4)
             title = yrmn//' Soil water; Bare layer 6 [m]'
             w_b6_mnth(:,:,i_mnth) = 
     &                w_b6_mnth(:,:,i_mnth)/n_count_mnth(i_mnth)
             write(1008) title,real(w_b6_mnth(:,:,i_mnth),kind=4)
             title = yrmn//' Soil water; Veg layer 1 [m]'
             w_v1_mnth(:,:,i_mnth) = 
     &                w_v1_mnth(:,:,i_mnth)/n_count_mnth(i_mnth)
             write(1009) title,real(w_v1_mnth(:,:,i_mnth),kind=4)
             title = yrmn//' Soil water; Veg layer 2 [m]'
             w_v2_mnth(:,:,i_mnth) = 
     &                w_v2_mnth(:,:,i_mnth)/n_count_mnth(i_mnth)
             write(1010) title,real(w_v2_mnth(:,:,i_mnth),kind=4)
             title = yrmn//' Soil water; Veg layer 3 [m]'
             w_v3_mnth(:,:,i_mnth) = 
     &                w_v3_mnth(:,:,i_mnth)/n_count_mnth(i_mnth)
             write(1011) title,real(w_v3_mnth(:,:,i_mnth),kind=4)
             title = yrmn//' Soil water; Veg layer 4 [m]'
             w_v4_mnth(:,:,i_mnth) = 
     &                w_v4_mnth(:,:,i_mnth)/n_count_mnth(i_mnth)
             write(1012) title,real(w_v4_mnth(:,:,i_mnth),kind=4)
             title = yrmn//' Soil water; Veg layer 5 [m]'
             w_v5_mnth(:,:,i_mnth) = 
     &                w_v5_mnth(:,:,i_mnth)/n_count_mnth(i_mnth)
             write(1013) title,real(w_v5_mnth(:,:,i_mnth),kind=4)
             title = yrmn//' Soil water; Veg layer 6 [m]'
             w_v6_mnth(:,:,i_mnth) = 
     &                w_v6_mnth(:,:,i_mnth)/n_count_mnth(i_mnth)
             write(1014) title,real(w_v6_mnth(:,:,i_mnth),kind=4)
             title = yrmn//' Soil heat; Bare layer 1 [J/m2]'
             ht_b1_mnth(:,:,i_mnth) = 
     &                ht_b1_mnth(:,:,i_mnth)/n_count_mnth(i_mnth)
             write(1015) title,real(ht_b1_mnth(:,:,i_mnth),kind=4)
             title = yrmn//' Soil heat; Bare layer 2 [J/m2]'
             ht_b2_mnth(:,:,i_mnth) = 
     &                ht_b2_mnth(:,:,i_mnth)/n_count_mnth(i_mnth)
             write(1016) title,real(ht_b2_mnth(:,:,i_mnth),kind=4)
             title = yrmn//' Soil heat; Bare layer 3 [J/m2]'
             ht_b3_mnth(:,:,i_mnth) = 
     &                ht_b3_mnth(:,:,i_mnth)/n_count_mnth(i_mnth)
             write(1017) title,real(ht_b3_mnth(:,:,i_mnth),kind=4)
             title = yrmn//' Soil heat; Bare layer 4 [J/m2]'
             ht_b4_mnth(:,:,i_mnth) = 
     &                ht_b4_mnth(:,:,i_mnth)/n_count_mnth(i_mnth)
             write(1018) title,real(ht_b4_mnth(:,:,i_mnth),kind=4)
             title = yrmn//' Soil heat; Bare layer 5 [J/m2]'
             ht_b5_mnth(:,:,i_mnth) = 
     &                ht_b5_mnth(:,:,i_mnth)/n_count_mnth(i_mnth)
             write(1019) title,real(ht_b5_mnth(:,:,i_mnth),kind=4)
             title = yrmn//' Soil heat; Bare layer 6 [J/m2]'
             ht_b6_mnth(:,:,i_mnth) = 
     &                ht_b6_mnth(:,:,i_mnth)/n_count_mnth(i_mnth)
             write(1020) title,real(ht_b6_mnth(:,:,i_mnth),kind=4)
             title = yrmn//' Soil heat; Veg layer 1 [J/m2]'
             ht_v1_mnth(:,:,i_mnth) = 
     &                ht_v1_mnth(:,:,i_mnth)/n_count_mnth(i_mnth)
             write(1021) title,real(ht_v1_mnth(:,:,i_mnth),kind=4)
             title = yrmn//' Soil heat; Veg layer 2 [J/m2]'
             ht_v2_mnth(:,:,i_mnth) = 
     &                ht_v2_mnth(:,:,i_mnth)/n_count_mnth(i_mnth)
             write(1022) title,real(ht_v2_mnth(:,:,i_mnth),kind=4)
             title = yrmn//' Soil heat; Veg layer 3 [J/m2]'
             ht_v3_mnth(:,:,i_mnth) = 
     &                ht_v3_mnth(:,:,i_mnth)/n_count_mnth(i_mnth)
             write(1023) title,real(ht_v3_mnth(:,:,i_mnth),kind=4)
             title = yrmn//' Soil heat; Veg layer 4 [J/m2]'
             ht_v4_mnth(:,:,i_mnth) = 
     &                ht_v4_mnth(:,:,i_mnth)/n_count_mnth(i_mnth)
             write(1024) title,real(ht_v4_mnth(:,:,i_mnth),kind=4)
             title = yrmn//' Soil heat; Veg layer 5 [J/m2]'
             ht_v5_mnth(:,:,i_mnth) = 
     &                ht_v5_mnth(:,:,i_mnth)/n_count_mnth(i_mnth)
             write(1025) title,real(ht_v5_mnth(:,:,i_mnth),kind=4)
             title = yrmn//' Soil heat; Veg layer 6 [J/m2]'
             ht_v6_mnth(:,:,i_mnth) = 
     &                ht_v6_mnth(:,:,i_mnth)/n_count_mnth(i_mnth)
             write(1026) title,real(ht_v6_mnth(:,:,i_mnth),kind=4)
             title = yrmn//' Snow layer thickness; bare layer 1 [m]'
             dzsn_b1(:,:,i_mnth) = 
     &                dzsn_b1(:,:,i_mnth)/n_count_mnth(i_mnth)
             write(1027) title,real(dzsn_b1(:,:,i_mnth),kind=4)
             title = yrmn//'  Snow layer thickness; bare layer 2 [m]'
             dzsn_b2(:,:,i_mnth) = 
     &                dzsn_b2(:,:,i_mnth)/n_count_mnth(i_mnth)
             write(1028) title,real(dzsn_b2(:,:,i_mnth),kind=4)
             title = yrmn//' Snow layer thickness; bare layer 3 [m]'
             dzsn_b3(:,:,i_mnth) = 
     &                dzsn_b3(:,:,i_mnth)/n_count_mnth(i_mnth)
             write(1029) title,real(dzsn_b3(:,:,i_mnth),kind=4)
             title = yrmn//' Snow layer thickness; bare layer 1 [m]'
             dzsn_v1(:,:,i_mnth) = 
     &                dzsn_v1(:,:,i_mnth)/n_count_mnth(i_mnth)
             write(1030) title,real(dzsn_v1(:,:,i_mnth),kind=4)
             title = yrmn//' Snow layer thickness; bare layer 2 [m]'
             dzsn_v2(:,:,i_mnth) = 
     &                dzsn_v2(:,:,i_mnth)/n_count_mnth(i_mnth)
             write(1031) title,real(dzsn_v2(:,:,i_mnth),kind=4)
             title = yrmn//' Snow layer thickness; bare layer 3 [m]'
             dzsn_v3(:,:,i_mnth) = 
     &                dzsn_v3(:,:,i_mnth)/n_count_mnth(i_mnth)
             write(1032) title,real(dzsn_v3(:,:,i_mnth),kind=4)
             title = yrmn//' Snow water equival; bare layer 1 [m]'
             wsn_b1(:,:,i_mnth) = 
     &                wsn_b1(:,:,i_mnth)/n_count_mnth(i_mnth)
             write(1033) title,real(wsn_b1(:,:,i_mnth),kind=4)
             title = yrmn//' Snow water equival; bare layer 2 [m]'
             wsn_b2(:,:,i_mnth) = 
     &                wsn_b2(:,:,i_mnth)/n_count_mnth(i_mnth)
             write(1034) title,real(wsn_b2(:,:,i_mnth),kind=4)
             title = yrmn//' Snow water equival; bare layer 3 [m]'
             wsn_b3(:,:,i_mnth) = 
     &                wsn_b3(:,:,i_mnth)/n_count_mnth(i_mnth)
             write(1035) title,real(wsn_b3(:,:,i_mnth),kind=4)
             title = yrmn//' Snow water equival; veg layer 1 [m]'
             wsn_v1(:,:,i_mnth) = 
     &                wsn_v1(:,:,i_mnth)/n_count_mnth(i_mnth)
             write(1036) title,real(wsn_v1(:,:,i_mnth),kind=4)
             title = yrmn//' Snow water equival; veg layer 2 [m]'
             wsn_v2(:,:,i_mnth) = 
     &                wsn_v2(:,:,i_mnth)/n_count_mnth(i_mnth)
             write(1037) title,real(wsn_v2(:,:,i_mnth),kind=4)
             title = yrmn//' Snow water equival; veg layer 3 [m]'
             wsn_v3(:,:,i_mnth) = 
     &                wsn_v3(:,:,i_mnth)/n_count_mnth(i_mnth)
             write(1038) title,real(wsn_v3(:,:,i_mnth),kind=4)
             title = yrmn//' Snow layer heat; bare layer 1 [J/m2]'
             hsn_b1(:,:,i_mnth) = 
     &                hsn_b1(:,:,i_mnth)/n_count_mnth(i_mnth)
             write(1039) title,real(hsn_b1(:,:,i_mnth),kind=4)
             title = yrmn//' Snow layer heat; bare layer 2 [J/m2]'
             hsn_b2(:,:,i_mnth) = 
     &                hsn_b2(:,:,i_mnth)/n_count_mnth(i_mnth)
             write(1040) title,real(hsn_b2(:,:,i_mnth),kind=4)
             title = yrmn//'  Snow layer heat; bare layer 3 [J/m2]'
             hsn_b3(:,:,i_mnth) = 
     &                hsn_b3(:,:,i_mnth)/n_count_mnth(i_mnth)
             write(1041) title,real(hsn_b3(:,:,i_mnth),kind=4)
             title = yrmn//'  Snow layer heat; veg layer 1 [J/m2]'
             hsn_v1(:,:,i_mnth) = 
     &                hsn_v1(:,:,i_mnth)/n_count_mnth(i_mnth)
             write(1042) title,real(hsn_v1(:,:,i_mnth),kind=4)
             title = yrmn//'  Snow layer heat; veg layer 2 [J/m2]'
             hsn_v2(:,:,i_mnth) = 
     &                hsn_v2(:,:,i_mnth)/n_count_mnth(i_mnth)
             write(1043) title,real(hsn_v2(:,:,i_mnth),kind=4)
             title = yrmn//'  Snow layer heat; veg layer 3 [J/m2]'
             hsn_v3(:,:,i_mnth) = 
     &                hsn_v3(:,:,i_mnth)/n_count_mnth(i_mnth)
             write(1044) title,real(hsn_v3(:,:,i_mnth),kind=4)
             title = yrmn//' # of snow layers on bare land'
             nsn_b(:,:,i_mnth) = 
     &                nsn_b(:,:,i_mnth)/n_count_mnth(i_mnth)
             write(1045) title,real(nsn_b(:,:,i_mnth),kind=4)
             title = yrmn//' # of snow layers on vegetated land'
             nsn_v(:,:,i_mnth) = 
     &                nsn_v(:,:,i_mnth)/n_count_mnth(i_mnth)
             write(1046) title,real(nsn_v(:,:,i_mnth),kind=4)
             title = yrmn//' Snow-covered bare fraction [-]'
             frsnow_b(:,:,i_mnth) = 
     &                frsnow_b(:,:,i_mnth)/n_count_mnth(i_mnth)
             write(1047) title,real(frsnow_b(:,:,i_mnth),kind=4)
             title = yrmn//'  Snow-covered vegetated fraction [-]'
             frsnow_v(:,:,i_mnth) = 
     &                frsnow_v(:,:,i_mnth)/n_count_mnth(i_mnth)
             write(1048) title,real(frsnow_v(:,:,i_mnth),kind=4)
             title = yrmn//'Foliage surf. vapor mixing ratio [kg/kg]'
             Qf_mnth(:,:,i_mnth) = 
     &                Qf_mnth(:,:,i_mnth)/n_count_mnth(i_mnth)
             write(1049) title,real(Qf_mnth(:,:,i_mnth),kind=4)
             title = yrmn//' root zone betad [-]'
             abetad_mnth(:,:,i_mnth) = 
     &                abetad_mnth(:,:,i_mnth)/n_count_mnth(i_mnth)
             write(1050) title,real(abetad_mnth(:,:,i_mnth),kind=4)
             title = yrmn//' Live leaf carbon pool [kgC/m2]'
             aClivepool_leaf_m(:,:,i_mnth) = 
     &              aClivepool_leaf_m(:,:,i_mnth)/n_count_mnth(i_mnth)
             write(1051) title,
     &              real(aClivepool_leaf_m(:,:,i_mnth),kind=4)
             title = yrmn//' Live froot carbon pool [kgC/m2]'
             aClivepool_froot_m(:,:,i_mnth) = 
     &              aClivepool_froot_m(:,:,i_mnth)/n_count_mnth(i_mnth)
             write(1052) title,
     &                    real(aClivepool_froot_m(:,:,i_mnth),kind=4)
             title = yrmn//' Live wood carbon pool [kgC/m2]'
             aClivepool_wood_m(:,:,i_mnth) = 
     &              aClivepool_wood_m(:,:,i_mnth)/n_count_mnth(i_mnth)
             write(1053) title,
     &                    real(aClivepool_wood_m(:,:,i_mnth),kind=4)
             title = yrmn//' Dead surface metabolic C pool [kgC/m2]'
             aCdeadpool_surfmet_m(:,:,i_mnth) = 
     &            aCdeadpool_surfmet_m(:,:,i_mnth)/n_count_mnth(i_mnth)
             write(1054) title,
     &                    real(aCdeadpool_surfmet_m(:,:,i_mnth),kind=4)
             title = yrmn//' Dead surface structural C pool [kgC/m2]'
             aCdeadpool_surfstr_m(:,:,i_mnth) = 
     &           aCdeadpool_surfstr_m(:,:,i_mnth)/n_count_mnth(i_mnth)
             write(1055) title,
     &                    real(aCdeadpool_surfstr_m(:,:,i_mnth),kind=4)
             title = yrmn//' Dead soil metabolic C pool [kgC/m2]'
             aCdeadpool_soilmet_m(:,:,i_mnth) = 
     &           aCdeadpool_soilmet_m(:,:,i_mnth)/n_count_mnth(i_mnth)
             write(1056) title,
     &                    real(aCdeadpool_soilmet_m(:,:,i_mnth),kind=4)
             title = yrmn//' Dead soil structural C pool [kgC/m2]'
             aCdeadpool_soilstr_m(:,:,i_mnth) = 
     &           aCdeadpool_soilstr_m(:,:,i_mnth)/n_count_mnth(i_mnth)
             write(1057) title,
     &                    real(aCdeadpool_soilstr_m(:,:,i_mnth),kind=4)
             title = yrmn//' Coarse woody debris C pool [kgC/m2]'
             aCdeadpool_cwd_m(:,:,i_mnth) = 
     &            aCdeadpool_cwd_m(:,:,i_mnth)/n_count_mnth(i_mnth)
             write(1058) title,
     &                    real(aCdeadpool_cwd_m(:,:,i_mnth),kind=4)
             title = yrmn//' Dead surface microbial C pool [kgC/m2]'
             aCdeadpool_surfmic_m(:,:,i_mnth) = 
     &            aCdeadpool_surfmic_m(:,:,i_mnth)/n_count_mnth(i_mnth)
             write(1059) title,
     &                    real(aCdeadpool_surfmic_m(:,:,i_mnth),kind=4)
             title = yrmn//' Dead soil microbial C pool [kgC/m2]'
             aCdeadpool_soilmic_m(:,:,i_mnth) = 
     &            aCdeadpool_soilmic_m(:,:,i_mnth)/n_count_mnth(i_mnth)
             write(1060) title,
     &                    real(aCdeadpool_soilmic_m(:,:,i_mnth),kind=4)
             title = yrmn//' Dead slow (~10 yr) C pool [kgC/m2]'
             aCdeadpool_slow_m(:,:,i_mnth) = 
     &            aCdeadpool_slow_m(:,:,i_mnth)/n_count_mnth(i_mnth)
             write(1061) title,
     &                    real(aCdeadpool_slow_m(:,:,i_mnth),kind=4)
             title = yrmn//' Dead passive (~100 yr) C pool [kgC/m2]'
             aCdeadpool_passive_m(:,:,i_mnth) = 
     &            aCdeadpool_passive_m(:,:,i_mnth)/n_count_mnth(i_mnth)
             write(1062) title,
     &                    real(aCdeadpool_passive_m(:,:,i_mnth),kind=4)
             title = yrmn//' Leaf Area Index [-]'
             alai_m(:,:,i_mnth)=alai_m(:,:,i_mnth)/n_count_mnth(i_mnth)
             write(1063) title,real(alai_m(:,:,i_mnth),kind=4)

             title = yrmn//' Water stored on the canopy [m]'
             canopyH2O_m(:,:,i_mnth)=canopyH2O_m(:,:,i_mnth)
     &                               /n_count_mnth(i_mnth)
             write(1064) title,real(canopyH2O_m(:,:,i_mnth),kind=4)

             title = yrmn//' Heat of the canopy [J/m2]'
             canopyheat_m(:,:,i_mnth)=canopyheat_m(:,:,i_mnth)
     &                                /n_count_mnth(i_mnth)
             write(1065) title,real(canopyheat_m(:,:,i_mnth),kind=4)

             !* Other sums of diagnostics for convenience.
             title = yrmn//'  Total plant carbon [kgC/m2]'
             write(1066) title, 
     &            real(aClivepool_leaf_m(:,:,i_mnth),kind=4) +
     &            real(aClivepool_froot_m(:,:,i_mnth),kind=4)+
     &            real(aClivepool_wood_m(:,:,i_mnth),kind=4)
             title = yrmn//'  Total soil carbon [kgC/m2]'
             write(1067) title, 
     &            real(aCdeadpool_surfmet_m(:,:,i_mnth),kind=4) +
     &            real(aCdeadpool_surfstr_m(:,:,i_mnth),kind=4) +
     &            real(aCdeadpool_soilmet_m(:,:,i_mnth),kind=4) +
     &            real(aCdeadpool_soilstr_m(:,:,i_mnth),kind=4) +
     &            real(aCdeadpool_cwd_m(:,:,i_mnth),kind=4) +
     &            real(aCdeadpool_surfmic_m(:,:,i_mnth),kind=4) +
     &            real(aCdeadpool_soilmic_m(:,:,i_mnth),kind=4) +
     &            real(aCdeadpool_slow_m(:,:,i_mnth),kind=4) +
     &            real(aCdeadpool_passive_m(:,:,i_mnth),kind=4)
          enddo
#endif
      if (mype==0) then
          print *, 'Finished printing diags.'
      endif
      end subroutine print_diags

!-----------------------------------------------------------------------
     /*  Zero Monthly diagnostic accumulators */
      subroutine zero_diags()     

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
!-----------------------------------------------------------------------

      subroutine accumulate_diags(i,j,k_mnth,s )
      use sle001, only : tp,tbcs,aevap,aevapw,aevapd,aevapb,aintercep
     &     ,aruns,arunu,agpp,arauto,asoilresp
     &     ,aevapvg,aevapvs,aevapbs,asrht,atrht,aalbedo
     &     ,aclab,asoilCpoolsum,aepp,atrg,ashg,alhg,aepc
     &     ,abetad !,alai - listed second time below
     &     ,aClivepool_leaf,aClivepool_froot,aClivepool_wood
     &     ,aCdeadpool_surfmet,aCdeadpool_surfstr,aCdeadpool_soilmet
     &     ,aCdeadpool_soilstr,aCdeadpool_cwd,aCdeadpool_surfmic
     &     ,aCdeadpool_soilmic,aCdeadpool_slow,aCdeadpool_passive
     &     ,alai, acnc, acna

      integer, intent(in) :: i,j,k_mnth
      type(t_lsm_state) :: s

       !print *, "Accumulating diags: i,j,k_mnth",i,j,k_mnth

          aevap_mnth(i,j,k_mnth)  = aevap_mnth(i,j,k_mnth) + aevap
          aevapw_mnth(i,j,k_mnth) = aevapw_mnth(i,j,k_mnth) + aevapw
          aevapd_mnth(i,j,k_mnth) = aevapd_mnth(i,j,k_mnth) + aevapd
          aevapb_mnth(i,j,k_mnth) = aevapb_mnth(i,j,k_mnth) + aevapb
          aintercep_mnth(i,j,k_mnth) = aintercep_mnth(i,j,k_mnth) + 
     &                                 aintercep
          aruns_mnth(i,j,k_mnth)  = aruns_mnth(i,j,k_mnth) + aruns
          arunu_mnth(i,j,k_mnth)  = arunu_mnth(i,j,k_mnth) + arunu
          agpp_mnth(i,j,k_mnth)   = agpp_mnth(i,j,k_mnth) + agpp*
     &         (1.d-3)/86400 *1.d9  !g/m2/day to 1e-9 kg/m2/s
          arauto_mnth(i,j,k_mnth) = arauto_mnth(i,j,k_mnth) + arauto*
     &         (1.d-3)/86400 *1.d9  !g/m2/day to 1e-9 kg/m2/s
          asoilresp_mnth(i,j,k_mnth)= asoilresp_mnth(i,j,k_mnth)+ 
     &                                asoilresp*
     &         (1.d-3)/86400 *1.d9  !g/m2/day to 1e-9 kg/m2/s
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

          !* Convert from gC to kgC
          aClivepool_leaf_m(i,j,k_mnth) = 
     &        aClivepool_leaf_m(i,j,k_mnth)  + aClivepool_leaf*1.d-3
          aClivepool_froot_m(i,j,k_mnth) = 
     &        aClivepool_froot_m(i,j,k_mnth) + aClivepool_froot*1.d-3
          aClivepool_wood_m(i,j,k_mnth) =
     &        aClivepool_wood_m(i,j,k_mnth)  + aClivepool_wood*1.d-3

          aCdeadpool_surfmet_m(i,j,k_mnth) =
     &        aCdeadpool_surfmet_m(i,j,k_mnth) +aCdeadpool_surfmet*1.d-3
          aCdeadpool_surfstr_m(i,j,k_mnth) =
     &        aCdeadpool_surfstr_m(i,j,k_mnth) +aCdeadpool_surfstr*1.d-3
          aCdeadpool_soilmet_m(i,j,k_mnth) =
     &        aCdeadpool_soilmet_m(i,j,k_mnth) +aCdeadpool_soilmet*1.d-3
          aCdeadpool_soilstr_m(i,j,k_mnth) =
     &        aCdeadpool_soilstr_m(i,j,k_mnth) +aCdeadpool_soilstr*1.d-3
          aCdeadpool_cwd_m(i,j,k_mnth) =
     &        aCdeadpool_cwd_m(i,j,k_mnth) +aCdeadpool_cwd*1.d-3
          aCdeadpool_surfmic_m(i,j,k_mnth) =
     &        aCdeadpool_surfmic_m(i,j,k_mnth) +aCdeadpool_surfmic*1.d-3
          aCdeadpool_soilmic_m(i,j,k_mnth) =
     &        aCdeadpool_soilmic_m(i,j,k_mnth) +aCdeadpool_soilmic*1.d-3
          aCdeadpool_slow_m(i,j,k_mnth) =
     &        aCdeadpool_slow_m(i,j,k_mnth) +aCdeadpool_slow*1.d-3
          aCdeadpool_passive_m(i,j,k_mnth) =
     &        aCdeadpool_passive_m(i,j,k_mnth) +aCdeadpool_passive*1.d-3

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
          end subroutine accumulate_diags

!-----------------------------------------------------------------------
#ifdef LSM_DRV_SINGLE      
      subroutine print_diags_single(i,j)
      implicit none
      integer, intent(in) :: i, j
      write(9995,'(150(1pe16.8))')
     &    force_prec_ms (i,j)        ,
     &    force_eprec_w (i,j)        ,
     &    force_sprec_ms (i,j)       ,
     &    force_seprec_w (i,j)       ,
     &    force_srheat (i,j)         ,
     &    force_trheat (i,j)         ,
     &    force_ts (i,j)             ,
     &    force_qs (i,j)             ,
     &    force_ps (i,j)             ,
     &    force_rhosrf (i,j)         ,
     &    force_cdh (i,j)            ,
     &    force_qm1 (i,j)            ,
     &    force_ws (i,j)             ,!13
     &    aevap,aevapw,aevapd,aevapb,aintercep,aruns,arunu, !20
     &    agpp,arauto,asoilresp,abetad, !24
     &    aevapvg,aevapvs,aevapbs,asrht,atrht,aalbedo, !30
     &    aclab,asoilCpoolsum,aepp,atrg,ashg,alhg,aepc,!37
     &    s%w_ij(1,1,i,j),s%w_ij(2,1,i,j),s%w_ij(3,1,i,j),!40
     &    s%w_ij(4,1,i,j),s%w_ij(5,1,i,j),s%w_ij(6,1,i,j),!43
     &    s%w_ij(1,2,i,j),s%w_ij(2,2,i,j),s%w_ij(3,2,i,j),!46
     &    s%w_ij(4,2,i,j),s%w_ij(5,2,i,j),s%w_ij(6,2,i,j),!49
     &    s%ht_ij(1,1,i,j),s%ht_ij(2,1,i,j),s%ht_ij(3,1,i,j),!52
     &    s%ht_ij(4,1,i,j),s%ht_ij(5,1,i,j),s%ht_ij(6,1,i,j),!55
     &    s%ht_ij(1,2,i,j),s%ht_ij(2,2,i,j),s%ht_ij(3,2,i,j),!58
     &    s%ht_ij(4,2,i,j),s%ht_ij(5,2,i,j),s%ht_ij(6,2,i,j),!61
     &    s%dzsn_ij(1,1,i,j),s%dzsn_ij(2,1,i,j),s%dzsn_ij(3,1,i,j),!64
     &    s%dzsn_ij(1,2,i,j),s%dzsn_ij(2,2,i,j),s%dzsn_ij(3,2,i,j),!67
     &    s%wsn_ij(1,1,i,j),s%wsn_ij(2,1,i,j),s%wsn_ij(3,1,i,j),!70
     &    s%wsn_ij(1,2,i,j),s%wsn_ij(2,2,i,j),s%wsn_ij(3,2,i,j),!73
     &    s%hsn_ij(1,1,i,j),s%hsn_ij(2,1,i,j),s%hsn_ij(3,1,i,j),!76
     &    s%hsn_ij(1,2,i,j),s%hsn_ij(2,2,i,j),s%hsn_ij(3,2,i,j),!79
     &    real(s%nsn_ij(1,i,j)),real(s%nsn_ij(2,i,j)),!81
     &    s%fr_snow_ij(1,i,j), s%fr_snow_ij(2,i,j),!83
     &    s%Qf_ij(i,j),!84
     &    aClivepool_leaf,aClivepool_froot,aClivepool_wood,!87
     &    aCdeadpool_surfmet,aCdeadpool_surfstr,aCdeadpool_soilmet,!90
     &    aCdeadpool_soilstr,aCdeadpool_cwd,aCdeadpool_surfmic,!93
     &    aCdeadpool_soilmic,aCdeadpool_slow,aCdeadpool_passive,!96
     &    alai,s%w_ij(0,2,i,j),s%ht_ij(0,2,i,j), !99
     &    acnc,acna, !101
     &    force_vis_rad(i,j),force_direct_vis_rad(i,j),tp(0,2), !104  
     &    tp(1,1),tp(2,1),tp(3,1),tp(4,1),tp(5,1),tp(6,1), !110 
     &    tp(1,2),tp(2,2),tp(3,2),tp(4,2),tp(5,2),tp(6,2)  !116
      end subroutine print_diags_single
#endif
!LSM_DRV_SINGLE

#endif
!PRINT_DIAGS

!-----------------------------------------------------------------------
      function calc_steps(dtsec,tyr_sec2,year1,year2,
     &     spin_start_yr,spin_end_yr,num_times_spin)
      !Calculate number of simulation steps or other steps of timestep dtsec.
      use lsm_phys_util, only : LEAPYR
      implicit none
      integer :: calc_steps
      real*8 :: dtsec
      integer :: tyr_sec2,year1,year2,
     &     spin_start_yr,spin_end_yr,num_times_spin
      integer :: i,j,steps,days

      days = 0
      do i=spin_start_yr,spin_end_yr
           if (LEAPYR(i)) then
              days = days + 366
           else
              days = days + 365
           endif
      enddo
      days = days * num_times_spin
      do i=year1,year2-1
           if (LEAPYR(i)) then
              days = days + 366
           else
              days = days + 365
           endif
      enddo
      calc_steps = floor((days*24.*60.*60. + tyr_sec2)/dtsec) + 1

      end function calc_steps
!-----------------------------------------------------------------------

      end module lsm

!=======================================================================

      subroutine LSM_standalone
      use parser
      use param
      use filemanager, only : openunit,closeunit
      use lsm
      use ent_mod
      use drv_veg, only : rewind_files_veg
      use drv_met, only : init_forcings_met, get_forcings_met !get_gswp_forcings
     &     ,rewind_files_met, deallocate_drv_met
      use drv_veg, only : deallocate_drv_veg
      use lsm_phys_util, only: month
      use domain_decomp, only : app_init, mype
      implicit none
#ifdef USE_ESMF
#include "mpi_defs.h"
#include "mpif.h"
#endif

      !LSM main state.
      type (t_lsm_bc) :: bc
      type(t_lsm_state) :: lsm_state

      !Input files.
!      character*120 :: LAI_file='fluxnet_LAI'
!      character*120 :: vheight_file='fluxnet_vheight'


!--- Run parameters:  array bounds, time, spin-up specs.------------------------

      integer :: IM, JM         !Grid resolution
      integer :: I0,I1,J0,J1    !Run boundaries

!Input file array boundaries.  If MIXED_VEG, array bounds for veg not needed.
      integer :: i0v, i1v,j0v, j1v     !Vegetation LAI driver file boundaries. 
      integer :: i0vh, i1vh,j0vh, j1vh !Vegetation height driver file boundaries.
      integer :: i0vc, i1vc,j0vc, j1vc !Vegetation cover file boundaries.
      integer :: i0cc, i1cc,j0cc, j1cc !Crop cover file boundaries.
      integer :: i0h, i1h,j0h, j1h     !Hydrology LSM topo files boundaries.
      integer :: i0s, i1s,j0s, j1s     !Hydrology state files boundaries.
      integer :: i0sc, i1sc,j0sc, j1sc !Soil carbon state file boundaries.
      integer :: i0st, i1st,j0st, j1st !Soil texture file boundaries.
      integer :: i0d, i1d,j0d, j1d     !Meteorological driver files boundaries.

!@var Gregorian year start and finish (julian year for modelE met forcings)
      integer :: year1,year2
!@var Start time in year1, end time in year2 (seconds from beginning of year)
!     Seconds make it consistent with restart files.
      integer :: tyr_sec1,tyr_sec2
!@var Timestep size of model [seconds] ModelE has dt in real*8.
      real*8  :: dtsec
!@var Timestep size of met forcings data [seconds], GSWP2 is 3 hr in seconds
      integer  :: METFILE_dtsec
!@var Number of rows in a met data forcing file, at time steps of METFILE_dtsec.
!      integer :: METFILE_rows
!@var Year of start of met data forcings.  THIS ASSUMES ALL FILES START Jan 1.
!+    IF DATA DO NOT START ON JAN 1, THEN DUMMY ROWS MUST BE ADDED.
      integer :: METFILE_tstart
!@var Timestep size of veg forcings data [seconds]
      integer :: VEGFILE_dtsec
!@var Number of rows in a veg data forcing file, at time steps of VEGFILE_dtsec.
      integer :: VEGFILE_rows
!@var Year of start of met data forcings.  THIS ASSUMES ALL FILES START Jan 1.
!+    IF DATA DO NOT START ON JAN 1, THEN DUMMY ROWS MUST BE ADDED.
      integer :: VEGFILE_tstart
!@var Counter for spinup years
      integer :: i_spin 
      integer :: spin_start_yr, spin_end_yr, num_times_spin

!     Logical flags and input paramters for the Ent DGTEM
      logical :: force_VEG
      logical :: do_soilinit, do_soilresp
      logical :: do_phenology_activegrowth, do_structuralgrowth
      logical :: do_frost_hardiness, do_patchdynamics, do_init_geo
      logical :: do_spinup

!     Base directory for met drivers
      character*80 :: BASEDIR_METDRV

      !------ End, run parameters --------------------------------------------

      !Forcing variables and LSM prognostic variables to drive Ent.
      real*8,dimension(:,:),allocatable ::  !dimension(I0:I1,J0:J1)
     &     tbcs_lsm,tcanopy_lsm,
     &     force_Ca             ,
     &     force_cos_zen_angle  ,
     &     force_vis_rad        ,
     &     force_direct_vis_rad ,
     &     force_prec_ms        ,
     &     force_eprec_w        ,
     &     force_sprec_ms       ,
     &     force_seprec_w       ,
     &     force_srheat         ,
     &     force_trheat         ,
     &     force_ts             ,
     &     force_qs             ,
     &     force_ps             ,
     &     force_rhosrf         ,
     &     force_cdh            ,
     &     force_qm1            ,
     &     force_ws             ,
     &     force_pbl_args_ws0   ,
     &     force_pbl_args_tprime,
     &     force_pbl_args_qprime

      !-------Local---------------------

      integer, parameter :: itime0=1 !Starting simulation timestep number [-]
      integer, parameter :: DAYSEC = 86400 !Seconds in day
      real*8, parameter :: DAYSECREAL = 86400. !Seconds in day

!      real*8  :: time
!      Timestep number of size 3 hours [-]
      integer :: itime_3hr
!      logical :: end_of_input_flag
      real*8 :: num_input_steps
      integer :: max_input_steps
      integer :: itime, max_itime
      integer :: ierr
      logical :: end_of_day_flag
      integer :: i_mnth !Month counter
      integer, allocatable :: hemi(:,:)

      integer :: tyr_sec        !Time from beginning of year [seconds]
      integer, save :: jday_old = -32768  !Save jday last state
      integer :: jday           !Day of year in run
      integer :: jyear          !Year A.D. in run
      character*80 :: logcomment !String for passing to subroutine for printing.

!     Read input file with run settings
      call process_input_parameters(
     &     IM,JM,I0,I1,J0,J1
     &     ,i0v, i1v, j0v, j1v         ,i0vh, i1vh, j0vh, j1vh
     &     ,i0vc, i1vc, j0vc, j1vc     ,i0cc, i1cc, j0cc, j1cc
     &     ,i0h, i1h, j0h, j1h         ,i0s, i1s, j0s, j1s
     &     ,i0sc, i1sc,j0sc, j1sc      ,i0st, i1st, j0st, j1st
     &     ,i0d, i1d, j0d, j1d
     &     ,year1,year2, tyr_sec1,tyr_sec2 
     &     ,dtsec,METFILE_dtsec, METFILE_tstart
     &     ,VEGFILE_dtsec, VEGFILE_rows, VEGFILE_tstart
     &     ,spin_start_yr, spin_end_yr, num_times_spin
     &     ,force_VEG
     &     ,do_soilinit, do_soilresp
     &     ,do_phenology_activegrowth, do_structuralgrowth
     &     ,do_frost_hardiness, do_patchdynamics, do_init_geo,do_spinup
     &     ,BASEDIR_METDRV )

      tyr_sec = tyr_sec1  ! + year1 * YEARSEC
      jday = floor(tyr_sec1/DAYSECREAL) + 1
      jday_old = jday !initialize
      if (do_spinup) then
         jyear=spin_start_yr
      else
         jyear = year1
      endif


      call app_init(jm,j0,j1)
      print *,mype,"jm =",jm,"domain decomposition: j0,j1 = ", j0,j1

      allocate(
     &     tbcs_lsm(I0:I1,J0:J1), tcanopy_lsm(I0:I1,J0:J1),
     &     force_Ca             (I0:I1,J0:J1),
     &     force_cos_zen_angle  (I0:I1,J0:J1),
     &     force_vis_rad        (I0:I1,J0:J1),
     &     force_direct_vis_rad (I0:I1,J0:J1),
     &     force_prec_ms        (I0:I1,J0:J1),
     &     force_eprec_w        (I0:I1,J0:J1),
     &     force_sprec_ms       (I0:I1,J0:J1),
     &     force_seprec_w       (I0:I1,J0:J1),
     &     force_srheat         (I0:I1,J0:J1),
     &     force_trheat         (I0:I1,J0:J1),
     &     force_ts             (I0:I1,J0:J1),
     &     force_qs             (I0:I1,J0:J1),
     &     force_ps             (I0:I1,J0:J1),
     &     force_rhosrf         (I0:I1,J0:J1),
     &     force_cdh            (I0:I1,J0:J1),
     &     force_qm1            (I0:I1,J0:J1),
     &     force_ws             (I0:I1,J0:J1),
     &     force_pbl_args_ws0   (I0:I1,J0:J1),
     &     force_pbl_args_tprime(I0:I1,J0:J1),
     &     force_pbl_args_qprime(I0:I1,J0:J1) )

#ifdef PRINT_DIAGS
      call allocate_diags(IM,JM,I0,I1,J0,J1)
#endif

!     ---- Set up run configuration.  Initialize GHY and ENT.---


      !* Set run configuration.
      call lsm_set_config(dtsec 
     &     ,force_VEG, do_soilinit, do_soilresp
     &     ,do_phenology_activegrowth, do_structuralgrowth
     &     ,do_frost_hardiness, do_patchdynamics,do_init_geo, do_spinup)
      print *, mype,"Set LSM run config."

      !* Set geological and soil structure
      call lsm_set_bc(IM,JM,I0,I1,J0,J1,i0h,i1h,j0h,j1h, bc ) 
      print *, mype,"Set LSM topo."

      !* Allocate lsm state
      call lsm_allocate_state(IM,JM,I0,I1,J0,J1,bc%fearth,lsm_state )
      !* Initialize hydrology and veg/soil
      call lsm_init_state_hydro( IM, JM, I0, I1, J0, J1
     &     ,i0s, i1s, j0s, j1s, lsm_state )
      print *, mype,"Initialized LSM hydro state."

      call lsm_init_state_vegsoil(force_VEG, do_soilinit
     &     ,do_phenology_activegrowth, do_init_geo
     &     , jday, jyear
     &     ,IM,JM,I0,I1,J0,J1
     &     ,i0v, i1v, j0v, j1v         ,i0vh, i1vh, j0vh, j1vh
     &     ,i0vc, i1vc, j0vc, j1vc     ,i0cc, i1cc, j0cc, j1cc
     &     ,i0sc, i1sc,j0sc, j1sc      ,i0st, i1st,j0st, j1st
     &     ,bc%fearth, lsm_state%entcells)
      print *,mype,"Initialized LSM veg and soil state."

      !* Initialize surface meteorology from driver file. 
      !* Read first driver value then rewind to beginning.
      call init_forcings_met(IM,JM,I0,I1,J0,J1,i0d,i1d,j0d,j1d,
     &     jday,jyear,tyr_sec,dtsec, BASEDIR_METDRV,
     &     tbcs_lsm,tcanopy_lsm,
     &     force_Ca             ,     force_cos_zen_angle  ,
     &     force_vis_rad        ,     force_direct_vis_rad , 
     &     force_prec_ms        ,     force_eprec_w        ,
     &     force_sprec_ms       ,     force_seprec_w       ,
     &     force_srheat         ,     force_trheat         ,
     &     force_ts             ,     force_qs             ,
     &     force_ps             ,     force_rhosrf         ,
     &     force_cdh            ,     force_qm1            ,
     &     force_ws             ,     force_pbl_args_ws0   ,
     &     force_pbl_args_tprime,     force_pbl_args_qprime)
      print *,mype,"Initialized met forcings."
!      end_of_input_flag = .false.

      !* Initialize vegetation pools, replace LAI, height, and Clabile.
      allocate( hemi(I0:I1,J0:J1))
      if ( J0<=JM/2 )   hemi(:,J0:min(JM/2,J1))   = -1    ! S.
      if ( J1>=JM/2+1 ) hemi(:,max(JM/2+1,J0):J1) =  1    ! N.
      if (force_VEG) then
          call lsm_init_forcings_veg(force_VEG
     &     ,do_phenology_activegrowth
     &     ,tyr_sec,DAYSEC,VEGFILE_dtsec, VEGFILE_tstart
     &     ,VEGFILE_ROWS, N_PFT, hemi, jday, jyear
     &     ,I0,I1,J0,J1,I0v,I1v,J0v,J1v,lsm_state%entcells)
          print *, mype,"Initialized forced veg."
      endif

      !call sysusage(mype,0)
      !call sysusage(mype,1)
      !call sysusage(mype+4,0)
      !call sysusage(mype+8,0)
      !call sysusage(mype+12,0)

      !-------- Main Time Loop --------------------
!     num_input_steps = dble(itime1)*dtsec/INPUT_dtsec
      max_input_steps = floor( (365*DAYSECREAL*(
     &     + (spin_end_yr-spin_start_yr+1)*num_times_spin
     &     + (year2-year1+1))
     &     + tyr_sec2 - METFILE_dtsec)
     &     / METFILE_dtsec)
      max_itime = calc_steps(dtsec,tyr_sec2,year1,year2,
     &     spin_start_yr,spin_end_yr,num_times_spin)

!      max_itime = floor((365*DAYSECREAL*(
!     &     + (spin_end_yr-spin_start_yr+1)*num_times_spin
!     &     + (year2-year1+1))
!     &     + tyr_sec2 - dtsec)
!     &     / dtsec)
      i_spin = 1
      tyr_sec = tyr_sec1

      print *, mype,"Starting run... max_input_steps, max_itime:",
     &   max_input_steps,max_itime
      if (mype==0) then
         print *, "jyear",jyear," month",month(jyear,tyr_sec)
     &        ," jday",jday, tyr_sec !write to log file
      endif
#ifndef USE_ESMF         
         write(0,*), "jyear",jyear," month",month(jyear,tyr_sec) !write to screen
     &          ," jday",jday, tyr_sec !write to log file
#endif

      !flush(6)
      loop_time_main: do itime=1,max_itime
        !call MPI_Barrier(MPI_COMM_WORLD, ierr)

      if ( jday .ne. jday_old ) then ! new day
         end_of_day_flag = .true.
         jday_old = jday
         if (mype==0) then
            print *, "jyear",jyear," month",month(jyear,tyr_sec)
     &          ," jday",jday, tyr_sec !write to log file
         endif
#ifndef USE_ESMF         
         write(0,*), "jyear",jyear," month",month(jyear,tyr_sec) !write to screen
     &          ," jday",jday, tyr_sec !write to log file
#endif

      else
         end_of_day_flag = .false.
      endif

!     * Update veg once a day.  
      !(Could move this to after lsm_run call, but may screw up diags).
!      if (.not.do_phenology_activegrowth) then !prescribed or file updates
      !Call this even for do_phenology_activegrowth=true to temporarily
      !update with Matthews albedoes by pft. - NK
         call update_veg(
     i        hemi,jday,jyear,tyr_sec,
     i        IM,JM,I0,I1,J0,J1,I0v,I1v,J0v,J1v,I0vh,I1vh,J0vh,J1vh,
     i        I0vc,I1vc,J0vc,J1vc,I0cc,I1cc,J0cc,J1cc,
     i        force_VEG,do_phenology_activegrowth,
     i        .false.,.false.,  !update_crop flag, init flag
     i        DAYSEC, VEGFILE_dtsec,
     o        lsm_state%entcells)
!      endif

      jday_old = jday

!     * Get met forcings
!      call get_gswp_forcings(
      call get_forcings_met(
     &     IM,JM,I0,I1,J0,J1, i0d,i1d,j0d,j1d, ! latd,
     &     jday, jyear, tyr_sec, dtsec, 
!     &     end_of_input_flag, 
     &     tbcs_lsm              ,         tcanopy_lsm           ,   
     &     force_Ca              ,         force_cos_zen_angle   ,
     &     force_vis_rad         ,         force_direct_vis_rad  ,
     &     force_prec_ms         ,         force_eprec_w         ,
     &     force_sprec_ms        ,         force_seprec_w        ,
     &     force_srheat          ,         force_trheat          ,
     &     force_ts              ,         force_qs              ,
     &     force_ps              ,         force_rhosrf          ,
     &     force_cdh             ,         force_qm1             ,
     &     force_ws              ,         force_pbl_args_ws0    ,
     &     force_pbl_args_tprime ,         force_pbl_args_qprime )

!     * Run the LSM one timestep of size dtsec
      call lsm_run(lsm_state,bc,jday,jyear,tyr_sec,dtsec,
     &     I0,I1,J0,J1,
     &     tbcs_lsm              ,         tcanopy_lsm           ,   
     &     force_Ca              ,         force_cos_zen_angle   ,
     &     force_vis_rad         ,         force_direct_vis_rad  ,
     &     force_prec_ms         ,         force_eprec_w         ,
     &     force_sprec_ms        ,         force_seprec_w        ,
     &     force_srheat          ,         force_trheat          ,
     &     force_ts              ,         force_qs              ,
     &     force_ps              ,         force_rhosrf          ,
     &     force_cdh             ,         force_qm1             ,
     &     force_ws              ,         force_pbl_args_ws0    ,
     &     force_pbl_args_tprime ,         force_pbl_args_qprime ,
     &     end_of_day_flag)
      
!     * Update time-related variables
      tyr_sec = tyr_sec + nint(dtsec)
      !jday = ceiling(tyr_sec/(24*3600.))
      jday = floor(tyr_sec/DAYSECREAL) + 1
      
!     * End of year check: reset time variables and print diags
        if (      (jday > 365 .and. mod(jyear,4).ne.0) ! non-leap
     &       .or.  jday > 366 ) then                   ! leap
#ifdef PRINT_DIAGS
          !Print diags at end of year
          !print *, 'Printing diags at end of year'
          logcomment = 'Printing diags at end of year'
          call print_diags(logcomment, jyear )
          call zero_diags()
#endif
          !* Update time counters for new year.
          jyear = jyear + 1
          jday = 1
          tyr_sec = 0
          if (do_spinup) then     !Spin-up check
             if (jyear>spin_end_yr)then
                   i_spin = i_spin + 1
                   if (i_spin.le.num_times_spin)then !Reset-year check
                      jyear = spin_start_yr
                      call rewind_files_met
                      if (force_VEG) then
                         call rewind_files_veg(VEGFILE_dtsec,
     &                                         VEGFILE_tstart)
                      endif
                   else
                	print *,'Finished spin-up.'
                	!Reset soil carbon after spinup if prescribed init.
                        if (do_soilinit) then
                           call prescribe_soilcarbon(do_soilinit
     &                          ,IM,JM,I0,I1,J0,J1,lsm_state%entcells)
                        endif
                        ! Set year to year1 of run (spin-up years independent of run years)
                        ! Generally this should be the year after the last spin-up year.
                	jyear = year1
                	do_spinup=.false.
                  endif
             endif
          endif


        endif ! End of year check

        !print *,"jyear, jday, tyr_sec",jyear,jday, tyr_sec
      enddo loop_time_main !----------------------------------------------------

#ifdef PRINT_DIAGS
          !Print diags if end of year not reached.
          if ((tyr_sec.gt.0.).and.(tyr_sec.lt.(365.*DAYSECREAL))) then
             !print *, 'Printing diags for partial year run'
             logcomment = 'Printing diags for partial year run'
             call print_diags(logcomment,jyear)
             !call zero_diags()
          endif
#endif

      !call sysusage(mype,2)
      !call sysusage(mype,3)
      !call sysusage(mype+4,3)
      !call sysusage(mype+8,3)
      !call sysusage(mype+12,3)

#ifdef USE_ESMF
      call MPI_Finalize(ierr)
#endif
      deallocate(
     &     tbcs_lsm,tcanopy_lsm,
     &     force_Ca             ,
     &     force_cos_zen_angle  ,
     &     force_vis_rad        ,
     &     force_direct_vis_rad ,
     &     force_prec_ms        ,
     &     force_eprec_w        ,
     &     force_sprec_ms       ,
     &     force_seprec_w       ,
     &     force_srheat         ,
     &     force_trheat         ,
     &     force_ts             ,
     &     force_qs             ,
     &     force_ps             ,
     &     force_rhosrf         ,
     &     force_cdh            ,
     &     force_qm1            ,
     &     force_ws             ,
     &     force_pbl_args_ws0   ,
     &     force_pbl_args_tprime,
     &     force_pbl_args_qprime )

      !call deallocate_drv_met
      !call deallocate_drv_veg(force_VEG)
      !call deallocate_LSM(bc, lsm_state)

       print *,"Completed run:  jyear",jyear,"jday",jday
     &     ,"tyr_sec",tyr_sec,mype
      end subroutine LSM_standalone
