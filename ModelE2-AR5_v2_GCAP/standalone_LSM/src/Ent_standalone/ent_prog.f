      module ent_prog_mod
!@sum Utilities for off-line run of Ent model.
#define ENT_STANDALONE_DIAG TRUE
#define PRINT_DRIVERS

      use ent_mod
      use FILEMANAGER
      implicit none

      contains

!************************************************************************
      subroutine ent_read_state( cells )
!@sum read ent state from the file
      type(entcelltype_public), intent(out) :: cells(:,:)
      !---
      integer, parameter :: MAX_BUFFER=100000 ! need realistic estimate
      real*8 buffer(MAX_BUFFER)
      integer iu_entstate
      integer ic, jc, i, j

      ic = size(cells,1)
      jc = size(cells,2)

      call openunit('ent_state',iu_entstate,.true.,.true.)
      do j=1,jc
        do i=1,ic
          read(iu_entstate) buffer
          ! check length of buffer : if( buffer(1) > MAX_BUFFER ) ??
          call ent_cell_unpack(buffer, cells(i,j))
        enddo
      enddo
      call closeunit(iu_entstate)

      end subroutine ent_read_state

!************************************************************************
      subroutine ent_write_state( cells )
!@sum write ent state to the file
      type(entcelltype_public), intent(in) :: cells(:,:)
      !---
      real*8, pointer :: buffer(:)
      integer iu_entstate
      integer ic, jc, i, j

      ic = size(cells,1)
      jc = size(cells,2)

      call openunit('ent_state_new',iu_entstate,.true.,.false.)
      do j=1,jc
        do i=1,ic
          call ent_cell_pack(buffer, cells(i,j))
          write(iu_entstate) buffer
          deallocate(buffer)
        enddo
      enddo
      call closeunit(iu_entstate)

      end subroutine ent_write_state

!************************************************************************
#ifdef MIXED_CANOPY
      subroutine ent_struct_construct (cells, IM,JM,I0,I1,J0,J1 )
!@sum ent_struct_construct  Read in csv file of Ent vegetation data structure.
!+    Loops at the entcell level and call ent_mod's ent_struct_setup to do
!+    patches and cohorts.
      !Same as do_ent_struct in Entcover_standalone ent_struct.f, except here
      ! cells are allocated before routine is called.
      implicit none
      integer :: IM,JM,I0,I1,J0,J1
      type(entcelltype_public), intent(inout) :: cells(I0:I1,J0:J1)
      !---------
      integer :: iu_entstruct
      integer :: iu_entstruct_ascii
      character*80, parameter :: ENTSTRUCT="ENTSTRUCT"
      integer :: N_CASAlayers    !N_CASA_layers
      integer :: i,j
      character :: next
      integer counter

!      call ent_cell_construct( cells ) !Gets called before this routine.

      call openunit(trim(ENTSTRUCT),iu_entstruct,.false.,.true.)
      read(iu_entstruct,*) N_CASAlayers
      write(*,*) 'N_CASAlayers',N_CASAlayers

      
      counter = 0
      do
        counter = counter + 1
        if ( counter >= 100000) 
     &       call stop_model("do_ent_struct: infinite loop?",255)        
        read(iu_entstruct,*) next
        if (next.eq.'$') then   
           write(*,*) 'Error in logic: end of entcell encountered.'
           exit
        else if (next.eq.'#') then !End of data
           write (*,*) 'Done.'
           exit
        else if (next.eq.'*') then !skip comment
        else if (next.eq.'e') then !new entcell
          write(*,*) 'e'
          read(iu_entstruct,*) i,j
          write(*,*) 'i,j',i,j
          call check_ij_bounds(i,j,I0,I1,J0,J1)
          call ent_struct_setup(cells(i,j),iu_entstruct)
        endif
      enddo
      call closeunit(iu_entstruct)

      write(*,*) 'Writing cellsbuf to ent_state_new...'
      call ent_write_state(cells)

      write(*,*) 'Writing  ascii file to ent_struct_ascii...'
      call openunit('ent_struct_ascii',
     &     iu_entstruct_ascii,.false.,.false.)
      do i=I0,I1
         do j=J0,J1
            call ent_cell_print(iu_entstruct_ascii, cells(i,j))
         enddo
      enddo
      call closeunit(iu_entstruct_ascii)
      end subroutine ent_struct_construct

!************************************************************************
      subroutine ent_struct_initphys(cells,IM,JM,I0,I1,J0,J1,
     &     do_soilinit,
!     &     do_phenology_activegrowth,
     &     do_read_from_files)
!@sum ent_struct_initphys.  Set soil physics and initial biophysics values.
      use ent_prescribed_drv, only : ent_struct_get_phys
      use ent_mod, only : ent_struct_initphys_cells

      type(entcelltype_public), intent(inout) :: cells(I0:I1,J0:J1)
      integer,intent(in) :: IM,JM,I0,I1,J0,J1
      logical,intent(in) :: do_soilinit
!      logical,intent(in) :: do_phenology_activegrowth
      logical,intent(in) :: do_read_from_files
      !---Local-----
      real*8, dimension(I0:I1,J0:J1) :: Ci_ini,CNC_ini,Tcan_ini,Qf_ini
      real*8, dimension(N_CASA_LAYERS,I0:I1,J0:J1) :: soil_C_total
!     integer, dimension(N_COVERTYPES) :: soil_color !soil types 1-bright 2-dark
      real*8, dimension(N_SOIL_TEXTURES,I0:I1,J0:J1) :: soil_texture
      real*8, dimension(N_PFT,PTRACE,NPOOLS-NLIVE,N_CASA_LAYERS,
     &                  I0:I1,J0:J1):: Tpooldata  !g/m2
      integer :: i,j
      
      call ent_struct_get_phys(IM,JM,I0, I1, J0, J1,
     &     Ci_ini, CNC_ini, Tcan_ini, Qf_ini,
     &     soil_texture,soil_C_total,Tpooldata,
     &     do_soilinit,do_read_from_files)

      do i=I0,I1
         do j=J0,J1
            call ent_struct_initphys_cells( cells(i,j),
     &           soil_texture(:,i,j), 
     &           Ci_ini(i,j),CNC_ini(i,j),Tcan_ini(i,j),Qf_ini(i,j),
     &           Tpooldata(:,:,:,:,i,j), .false.)
         end do
      end do

      end subroutine ent_struct_initphys

#endif !MIXED_CANOPY
!************************************************************************
#ifdef MIXED_CANOPY
      subroutine check_ij_bounds(i,j,I0,I1,J0,J1)
      use filemanager
      implicit none
      integer :: i,j,I0,I1,J0,J1

      if ((i.lt.I0).or.(i.gt.I1)) then
         call stop_model("ent_make_struct i out of bounds",255)
      elseif ((j.lt.J0).or.(j.gt.J1)) then
         call stop_model("ent_make_struct j out of bounds",255)
      endif

      end subroutine check_ij_bounds
#endif !MIXED_CANOPY
!************************************************************************

      subroutine ent_init_vegstruct( cells,
     &     IM, JM, I0, I1, J0, J1, jday, year,
     &     do_soilinit,do_phenology_activegrowth)!,mixed_VEG)
      use ent_prescribed_drv, only : init_canopy_physical,prescr_vegdata
     &     ,ent_init_params
      !use ent_prescr_veg, only : prescr_calcconst
#ifdef ENT_STANDALONE_DIAG
      use ent_prescr_veg, only : print_ent_pfts
      use ent_pfts, only : alamin, alamax
#endif
      use ent_const

      implicit none
      type(entcelltype_public), intent(inout) :: cells(I0:I1,J0:J1)
      integer, intent(in) :: IM, JM, I0, I1, J0, J1, jday, year
      logical, intent(in) :: do_soilinit, do_phenology_activegrowth
!      logical, intent(in) :: mixed_VEG
      !---Local variables-----
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
     &                  I0:I1,J0:J1):: Tpooldata  !in g/m2 
      !------
      integer :: iu,i,j
      
      !hemi(J0:JM/2)   = -1    ! S.
      !hemi(JM/2+1:J1) =  1    ! N.

      !Read land surface parameters or use defaults
      !prescr data sets:
      !call prescr_calcconst()
      call ent_init_params()

#ifdef ENT_STANDALONE_DIAG
      call print_ent_pfts()
      write(*,*) "alamin: ", alamin
      write(*,*) "alamax: ", alamax
#endif
 
#ifdef PFT_MODEL_ENT
      write(*,*) "PFT_MODEL_ENT defined"
#endif


!      if (.not.mixed_VEG) then
#ifndef MIXED_CANOPY

         !* Mosaicked subgrid fractions. *!
         !* NOTE:  This initializes with default LAI.
         call prescr_vegdata(jday, year, 
     &     IM,JM,I0,I1,J0,J1,vegdata,albedodata,laidata,hdata,nmdata,
     &     popdata,dbhdata,craddata,cpooldata,rootprofdata,
     &     soildata,soil_texture,Tpooldata,
     &     do_soilinit,do_phenology_activegrowth,do_init_geo=.false.,
     &     do_read_from_files=.true.)
         !print *,"popdata in ent_GISS_init: ",popdata
         print *,"N_PFT=",N_PFT
         print *, "Tpooldata in ent_GISS_init: ",Tpooldata
         print *,"soil_texture:",soil_texture
         print *,"vegdata:",vegdata
         print *,"laidata:",laidata
         
         !Translate gridded data to Entdata structure
         !GISS data:  a patch per veg cover fraction, one cohort per patch
         call init_canopy_physical(I0, I1, J0, J1,
     &        Ci_ini, CNC_ini, Tcan_ini, Qf_ini)

!         print *,'In ent_init_vegstruct before ent_cell_set'
!         call openunit('entcells_print_init',
!     &        iu,.false.,.false.)
!         write(iu,*)'In ent_init_vegstruct before ent_cell_set'
!         call loop_ent_cell_print(iu,I0,I1,J0,J1,cells)

         call ent_cell_set(cells, vegdata, popdata, laidata, 
     &        hdata, dbhdata, craddata, cpooldata, nmdata, rootprofdata, 
     &        soildata, albedodata, soil_texture,
     &        Ci_ini, CNC_ini, Tcan_ini, Qf_ini, Tpooldata,
     &        reinitialize=.true.)  
!         print *,'In ent_init_vegstruct after ent_cell_set'
!         write(iu,*)'In ent_init_vegstruct after ent_cell_set'
!         call loop_ent_cell_print(iu,I0,I1,J0,J1,cells)
!         call closeunit(iu)
         
!      else !mixed_VEG
#else 
!MIXED_CANOPY
         !* Patches with mixed-community structure.*!
         
         call ent_struct_construct(cells, IM,JM,I0,I1,J0,J1)
         write(*,*) 'Constructed ent struct...'
         call ent_struct_initphys(cells,IM,JM,I0,I1,J0,J1,
     &        do_soilinit,
!     &        do_phenology_activegrowth,
     &        .true.)
!      endif

#endif 
!MIXED_CANOPY

      end subroutine ent_init_vegstruct

!************************************************************************
      subroutine loop_ent_cell_print(iu,I0,I1,J0,J1,cells)
      integer :: iu, I0,I1,J0,J1,i,j
      type(entcelltype_public), intent(in) :: cells(I0:I1,J0:J1)
      
        do i=I0,I1
           do j=J0,J1
              call ent_cell_print(iu, cells(i,j))
           enddo
        enddo

      end subroutine loop_ent_cell_print


      subroutine loop_ent_cell_print_diag(iu,I0,I1,J0,J1,cells)
      integer :: iu, I0,I1,J0,J1,i,j
      type(entcelltype_public), intent(in) :: cells(I0:I1,J0:J1)
      
        do i=I0,I1
           do j=J0,J1
              call ent_cell_print_diag(iu, cells(i,j))
           enddo
        enddo

      end subroutine loop_ent_cell_print_diag

!************************************************************************

      end module ent_prog_mod
!************************************************************************




!************************************************************************
      !program ent_prog
      subroutine Ent_standalone
      ! this driver just passes the bounds of the region, so that all
      ! arrays can be allocated on a stack
      implicit none
      integer IM, JM, I0, I1, J0, J1, i0f, i1f, j0f, j1f
      integer jday, year, jday2, year2
      real*8 dt !seconds
      logical :: force_VEG
      logical :: do_soilinit,do_soilresp
      logical :: do_phenology_activegrowth, do_structuralgrowth
      logical :: do_frost_hardiness
      logical :: do_patchdynamics, do_init_geo, do_spinup
      integer :: skip !#HACK to skip records at beginning of forcing file
!      logical :: mixed_VEG

      print *, "Started."

      !* Default configuration
      dt = 1800.d0              !second. Default value may be over-ridden by ent_input
      force_VEG = .false.
      do_soilinit = .true.
      do_soilresp = .true.
      do_phenology_activegrowth = .false.
      do_structuralgrowth = .false.
      do_frost_hardiness = .true.
      do_patchdynamics = .false.
      do_init_geo = .false.
      do_spinup = .false.
      skip = 0
!      mixed_VEG = .false.

! Default values of IM, JM, I0, I1, J0, J1, jday, year, Ci_ini, CNC_ini, 
! Tcan_ini, Qf_ini should be set either here or inside of 
! "read_input_parameters"


      !* Set world bounds (should correspond to format of input files)
      IM = 72; JM = 46

      !* Set bounds of simulated region (i0=i1; j0=j1 for single cell)
      !* (all arrays will be allocated to this size)
      !* Full grid
      !i0=1; i1=IM
      !j0=1; j1=JM

      !* E.g. if one grid cell
      !i0 = 39; i1 = 39
      !j0 = 36; j1 = 36

      !* Single grid cell corresponding to Ponca lat/long
!      i0 = 17; i1 = 17
!      j0 = 33; j1 = 33

      !* Default, entire grid.
      i0 = 1; i1 = 72
      j0 = 1; j1 = 46

      !* dims to default forcings file
      i0f = 1; i1f = 72
      j0f = 1; j1f = 46

      !* Default date start for GISS GCM 10-day forcings
      jday = 152                !June 1
      year = 1980
      jday2 = jday + 10
      year2 = -1 !Initialize

      write(*,*) 'Got here before read_input_parameters'
      call read_input_parameters(IM, JM, I0, I1, J0, J1
     &     ,jday, year, jday2, year2,dt
     &     ,i0f, i1f, j0f, j1f
     &     ,force_VEG
     &     ,do_soilinit, do_soilresp, do_phenology_activegrowth
     &     ,do_structuralgrowth, do_frost_hardiness
     &     ,do_patchdynamics, do_init_geo, skip, do_spinup)!, mixed_VEG)
      if (year2.eq.-1) year2=year  !If year2 was not specified, default same as year.

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

      print *,"ent_input: "
     &     , jday, year, jday2, year2,dt
     &     , i0f, i1f, j0f, j1f
     &     , force_VEG
     &     , do_soilinit, do_soilresp, do_phenology_activegrowth
     &     , do_structuralgrowth, do_frost_hardiness
     &     , do_patchdynamics, do_init_geo,do_spinup !, mixed_VEG


      print *,"starting program"
      call run_offline(IM, JM, I0, I1, J0, J1, jday, year,jday2,year2,dt
     &     ,i0f, i1f, j0f, j1f
     &     ,force_VEG
     &     ,do_soilinit, do_soilresp, do_phenology_activegrowth
     &     ,do_structuralgrowth, do_frost_hardiness 
     &     ,do_patchdynamics, do_init_geo
     &     ,skip,do_spinup)!,mixed_VEG)


      print *,"Ent run completed."
      !end program ent_prog
      end subroutine Ent_standalone


!************************************************************************
      subroutine read_input_parameters(IM, JM, I0, I1, J0, J1,
     &     jday, year, jday2, year2, dt,i0f, i1f, j0f, j1f, 
     &     force_VEG, do_soilinit,
     &     do_soilresp, do_phenology_activegrowth,
     &     do_structuralgrowth,do_frost_hardiness,
     &     do_patchdynamics,do_init_geo,skip,do_spinup)!,mixed_VEG)

      use filemanager
      implicit none
      integer IM, JM, I0, I1, J0, J1, i0f, i1f, j0f, j1f
      integer jday, year,jday2, year2
      real*8 dt
      logical force_VEG
      logical do_soilinit, do_soilresp, do_frost_hardiness
      logical do_phenology_activegrowth, do_structuralgrowth
      logical do_patchdynamics, do_init_geo, do_spinup
!      logical mixed_VEG
      integer skip !#HACK to skip records at beginning of forcing file
      !----
      integer iu_ent_input
      character*80, parameter :: file_ent_input="ent_input"
      namelist /input_parameters/ IM, JM, I0, I1, J0, J1
     &     ,jday, year, jday2, year2,dt
     &     ,i0f, i1f, j0f, j1f, force_VEG, do_soilinit
     &     ,do_soilresp, do_phenology_activegrowth
     &     ,do_structuralgrowth, do_frost_hardiness 
     &     ,do_patchdynamics,do_init_geo, skip, do_spinup
!, mixed_VEG

      print *,"reading input parameters from ", file_ent_input
      call openunit(trim(file_ent_input),iu_ent_input,.false.,.true.)
      read(iu_ent_input, NML=input_parameters, ERR=10)
      call closeunit(iu_ent_input)

      return
 10   continue
      print *,"error reading namelist file:", file_ent_input
      stop 255
      end subroutine read_input_parameters


!************************************************************************
      subroutine run_offline(IM, JM, I0, I1, J0, J1
     &     ,jday, year, jday2, year2,dt
     &     ,i0f, i1f, j0f, j1f
     &     ,force_VEG, do_soilinit
     &     ,do_soilresp, do_phenology_activegrowth
     &     ,do_structuralgrowth, do_frost_hardiness
     &     ,do_patchdynamics,do_init_geo ,skip,do_spinup)!,mixed_VEG)
      !****************************************************************
      !* Example program to run Ent coupled to a GCM.
      !* This version assumes an explicit scheme for calculation of
      !* canopy conductance, photosynthesis, and temperature.
      !* - For parallelization, the GCM provides the grid bounds to Ent.
      !* - Ent initializes vegetation structure parameters within these bounds.
      !* - GCM provides initial surface meteorological state variables.
      !* - Ent initializes GCANOPY, Ci, Qf, zeroes GPP
      !* - Simulation loop:
      !*   a. GCM provides meteorological drivers and canopy temperature
      !*   b. Ent updates.
      !*   c. Program gets from Ent: GCANOPY, Ci, Qf, GPP, TRANS_SW, albedo
      !*   d. GCM updates Qf, Tcanopy, given GCANOPY, TRANS_SW, albedo
      !*      and runs tracers on C and N.
      !****************************************************************

      !GCM_coupler:  need to write these routines specific to GCM.
!      use ent_GCM_coupler, only: 
!     &     GCM__get_grid, GCM_get_time, GCM_getdrv_cell, GCM_EWB

      use ent_mod
      !use ent_prescrveg
      use ent_prog_mod
      use ent_forcings
      use filemanager
      use ent_const 
      !use ent_const, only : JEQUATOR
      
      implicit none
      integer, intent(in) :: IM, JM, I0, I1, J0, J1, i0f, i1f, j0f, j1f
      integer, intent(in) :: jday, year,jday2,year2
      logical, intent(in) :: force_VEG, do_soilinit
      logical, intent(in) :: do_soilresp,do_frost_hardiness
      logical, intent(in) :: do_phenology_activegrowth
      logical, intent(in) :: do_structuralgrowth
      logical, intent(in) :: do_patchdynamics
      logical, intent(in) :: do_init_geo

      integer, intent(in) :: skip !#HACK to skip records at beginning of forcing file
      logical, intent(in) :: do_spinup
!      logical, intent(in) :: mixed_VEG
      real*8, intent(in) :: dt
      !---Local----
      integer :: jdaycount  !Only needed for prognostic phenology
      logical :: do_rewind

      ! Run parameters - NOT NEEDED
!      integer, parameter :: BIOPHYSICS_ONLY = 1
!      integer, parameter :: GISSVEG_CASA = 2
!      integer, parameter :: RUN_TYPE = BIOPHYSICS_ONLY
!      integer, parameter :: RUN_TYPE = GISSVEG_CASA

      ! time settings
!      real*8, parameter :: dt = 1800.d0 !
!      real*8, parameter :: max_time = dt * 2!run for x time steps
!      real*8, parameter :: max_time = dt * 48.*10. !run for n days
!      real*8, parameter :: max_time = dt * 48.*365. !run for 1 year
      real*8 :: max_time !In seconds
      real*8,parameter :: max_days = 365 !In days
      real*8, parameter :: save_interval=86400.d0*10d0 !save every 10 days

      ! model-dependent aarray sizes (maybe should be USEs from ent_...)
      !integer,parameter :: NSOILLAYERS = 6
      !integer,parameter :: NALBBANDS = 6

      !Coupling variables
      !Cell-level summary values - CALCULATED BY GCM/EWB OR OFF-LINE FILE
      real*8, dimension(I0:I1,J0:J1) ::
     &     TairC                !Air temperature (Celsius) !KIM - for phenology
     &     ,TcanopyC            !Canopy temperature (Celsius)
     &     ,Qf                  !Foliage surface specif humidity (kg vapor/ kg air)
     &     ,P_mbar              !Atmospheric pressure (mb)
     &     ,Ca                  !@Atmos CO2 conc at surface height (mol/m3).
     &     ,Ch                  !Ground to surface heat transfer coefficient 
     &     ,U                   !Surface layer wind speed (m s-1)
      !Radiation may later be broken down into hyperspectral increments.
      ! in an array
     &     ,IPARdif             !Incident diffuse PAR (vis) (W m-2)
     &     ,IPARdir             !Incident direct PAR (vis) (W m-2)
     &     ,CosZen              !cos of solar zenith angle
!     &     ,Soiltemp            !soil temp avg top 30 cm (C) -PK 6/27/06
!     &     ,Soilmoist           !soil vol moist avg top 30 cm
!      real*8, dimension(N_CASA_LAYERS,I0:I1,J0:J1) ::
      real*8, dimension(N_DEPTH,I0:I1,J0:J1) ::  !Changed to GCM layering -NK
     &      Soiltemp           !soil temp 
     &     ,Soilmoist          !soil volum moist 
      real*8, dimension(N_DEPTH,I0:I1,J0:J1) ::
     &     Soilmp              !Soil matric potential
     &     ,fice                !Fraction of soil water that is ice.
      real*8, dimension(:,:,:), pointer ::
     &     LAI                  ! prescribed LAI (for all pft's in the cell)
     &     ,height              ! prescribed height (for all pft's in the cell)
!      real*8, dimension(:,:,:), pointer ::
!     &     heightm              ! prescribed height (for all pft's in the cell)
      

      !Coupling and diagnostic variables specific to Ent (output)
      real*8, dimension(I0:I1,J0:J1) ::
     &     GCANOPY              !Canopy conductance of water vapor (m s-1). 
     &     ,Ci                  !Internal foliage CO2 concentration (mol/m3)
     &     ,GPP                 !Gross primary productivity (kg[C]/m2/s).
     &     ,TRANS_SW            !Transmittance of shortwave through canopy to soil 
     &     ,z0, CO2flux         !Variables from Ent to GCM
     &     ,C_labile, R_auto    !Variables from Ent to GCM
      !un-comment next one if want to export (see get_ent_exports below) -PK
!     &     ,Soil_resp           !soil respiration, patch level (added for ent_get_exports) -PK
      real*8, dimension(N_DEPTH,I0:I1,J0:J1) ::
     &     betadl
      real*8, dimension(N_BANDS,I0:I1,J0:J1) ::
     &     albedo               !Variables from Ent to GCM
      !un-comment next one if want to export (see get_ent_exports below) -PK
!      real*8, dimension(PTRACE,NPOOLS,N_CASA_LAYERS,I0:I1,J0:J1) ::
!     &     Tpool  !now explicitly depth-structured -PK 7/07
      !---------------------------------------------------------------------

      type(entcelltype_public) :: cells(I0:I1,J0:J1)
      !integer iu_forcings
      integer iu_results
      real*8 time, time_since_last_save
      integer hemi(I0:I1,J0:J1) !hemisphere flags = 1 for N., =-1 for S.
      logical :: update_day    !For prescribed phenology, litter
      real*8 :: fw(I0:I1,J0:J1) ! fraction of wet canopy

      !---Other local vars
      integer :: i,j,iu_entcells !For printing entcells - NK

      print *,"started run_offline with:", IM, JM, I0, I1, J0, J1

      fw(:,:) = 0.d0   ! all canopy is dry

      !* Now allocate optional arrays
      if ( force_VEG ) then
        allocate( LAI(N_PFT,I0:I1,J0:J1) )
        allocate( height(N_PFT,I0:I1,J0:J1) )
        !allocate( heightm(N_PFT,I0:I1,J0:J1) )
      else
        nullify( LAI )
        nullify( height )
        !nullify(heightm)
      endif

      !* Now Ent should be initialized before any other calls to ent_*
      call ent_init_config(
     &     do_soilresp=do_soilresp
     &     ,do_phenology_activegrowth=do_phenology_activegrowth
     &     ,do_structuralgrowth=do_structuralgrowth
     &     ,do_frost_hardiness=do_frost_hardiness
     &     ,do_patchdynamics=do_patchdynamics)
!     &     ,mixed_VEG=mixed_VEG)

      !* Set hemisphere flags.
      if ( J0<=JM/2 )   hemi(:,J0:min(JM/2,J1))   = -1    ! S.
      if ( J1>=JM/2+1 ) hemi(:,max(JM/2+1,J0):J1) =  1    ! N.
      print *,"set hemi ok"
      
      !* Initialize ent cells (makes head entcell).
      call ent_cell_construct( cells )
      print *,"ent_cell_construct passed ok" 

      !if binary intialization file exists
      ! call ent_read_state( cells )
      !else initialize from scratch
      !jday = 152  !June 1
      !year = 1980

      !* Initialize vegetation cover.  
      call ent_init_vegstruct( cells, IM, JM, I0, I1, J0, J1, 
     &        jday, year, !jday and year only need for Matthews veg.
     &        do_soilinit,do_phenology_activegrowth)!,mixed_VEG)
      print *,"ent_init_vegstruct passed ok"
    
      !* Open files with forcings (one record per time step)
      !!call openunit('ent_forcings',iu_forcings,.true.,.true.)
!      write(995,*) 'Opening forcings files'
      call open_forcings_file( i0f, i1f, j0f, j1f, 
     &     force_VEG,N_CASA_LAYERS,skip)!,mixed_VEG)
 
      ! open file for ent results (to save on each time step)
      call openunit('ent_results',iu_results,.true.,.false.)

#ifdef ENT_STANDALONE_DIAG
      call openunit('entcells_print',
     &       iu_entcells,.false.,.false.)
      call loop_ent_cell_print_diag( iu_entcells, I0,I1,J0,J1,cells ) 
#ifdef PRINT_DRIVERS
      if (N_PFT.eq.8) then
         write(980,*) "I0 I1 J0 J1 N_DEPTH TairC  TcanopyC Qf P_mbar
     &Ca Ch U IPARdif IPARdir CosZen 
     &Smp1 Smp2 Smp3 Smp4 Smp5 Smp6
     &Sm1 Sm2 Sm3 Sm4 Sm5 Sm6
     &St1 St2 St3 St4 St5 St6
     &fice1 fice2 fice3 fice4 fice5 fice6 
     &LAI1 LAI2 LAI3 LAI4 LAI5 LAI6 LAI7 LAI8
     & h1 h2 h3 h4 h5 h6 h7 h8" 
      else
         write(980,*) "I0 I1 J0 J1 N_DEPTH TairC  TcanopyC Qf P_mbar
     &Ca Ch U IPARdif IPARdir CosZen
     &Smp1 Smp2 Smp3 Smp4 Smp5 Smp6
     &Sm1 Sm2 Sm3 Sm4 Sm5 Sm6
     &St1 St2 St3 St4 St5 St6
     &fice1 fice2 fice3 fice4 fice5 fice6 
     &LAI1 LAI2 LAI3 LAI4 LAI5 LAI6 LAI7 LAI8
     &LAI9  LAI10  LAI11  LAI12  LAI13  LAI14 LAI15 LAI16
     &h1 h2 h3 h4 h5 h6 h7 h8 h9 h10 h11 h12 h13 h14 h15 h16"
      endif

!!!! NEVER ever use explicit number for i/o unit (except 99) !!!
!!! or at least set it to a very big number to avoid conflicts with filemanager
      write(981,*) "GCANOPY Ci Qf GPP TRANS_SW z0 CO2flux 
     & betadl1 betadl2 betadl3 betadl4 betadl5 betadl6
     & albedo1 albedo2 albedo3 albedo4 albedo5 albedo6"
#endif
      !*********** DIAGNOSTICS FOR PLOTTING in ent.f *******************!
      if (N_CASA_LAYERS == 1) then  !labels depend on number of layers -PK 7/07
         write(995,*)"patchnum IPARdir IPARdif coszen pft n lai heightm
     & leaf froot wood surfmet surfstr soilmet                      
     & soilstr cwd surfmic soilmic slow passive
     & C_fol C_w C_froot C_root C_lab C_repro TRANS_SW Ci GPP Rauto
     & Soilresp NPP CO2flux GCANOPY IPP senescefrac Sacclim 
     & c_total c_growth litter betad"
      else if (N_CASA_LAYERS == 2) then
         write(995,*)"patchnum IPARdir IPARdif coszen pft n lai heightm
     & leaf1 froot1 wood1 surfmet1 surfstr1 soilmet1                    !suffixes 1,2 for 1st, 2nd soil bgc layers -PK  
     & soilstr1 cwd1 surfmic1 soilmic1 slow1 passive1
     & leaf2 froot2 wood2 surfmet2 surfstr2 soilmet2  
     & soilstr2 cwd2 surfmic2 soilmic2 slow2 passive2 
     & C_fol C_w C_froot C_root C_lab C_repro TRANS_SW
     & Ci GPP Rauto Soilresp NPP CO2flux GCANOPY IPP senescefrac"
      end if

!      write(996,*) "area IPARdir IPARdif LAI fv
!     & leaf froot wood surfmet surfstr soilmet                      
!     & soilstr cwd surfmic soilmic slow passive
!     & C_fol  C_w  C_froot  C_root  C_lab 
!     & TRANS_SW Ci  GPP R_auto Soil_resp NPP CO2flux GCANOPY"

!      write(998, *) "TcanopyC cop%GPP cop%R_root cop%R_auto cop%C_fol"

!#endif
!      write(982,*)"CiPa gasc tk vegpar:sigma sqrtexpr kdf rhor kbl alai 
!     &     nm vh Ntot vegalbedo ppar:pcp Kc Ko n1 n2 m1 m2 msat N0 Oi k"
#endif

      
      !* Start loop over time.
      !max_time = 86400*(jday2-jday+1)*(year2-year+1)
      !YK- in this way,we can start the simulation in the middle of year
      !and end in the middle of the other year, 
      !although still assuming 365 days for a year
      !max_time = 86400*((365-jday+1)+jday2+365*(year2-year-1)) !Old
      !This max_time gives total time between jday1 in year1 and jday2 in year2
      !max_time = 86400*(365*(year2-year1-1) - (jday1-1) - (365-jday2))
      !This max_time gives total time from beginning of year1 to jday2 in year2
      max_time = 86400*(365*(year2-year) + jday2)
      print *,"max_time: ",max_time
      !time = 0.d0
      time = 86400*(jday-1) !Allows starting run in middle of year at jday1.
      time_since_last_save = 0.d0
      jdaycount = jday
      update_day = .false. !Initialize
      do_rewind=.false.   !Initialize

      !* Initialize vegetation structure for forced veg structure.
        if (force_VEG) then
          call read_forcings(I0, I1, J0, J1, 
     &          N_DEPTH, N_CASA_LAYERS,N_PFT,
     &          TairC, TcanopyC, Qf, P_mbar, Ca, Ch, U,
     &          IPARdif, IPARdir, CosZen,
     &          Soiltemp, Soilmoist, Soilmp, fice, LAI, height, 
     &          force_VEG, do_rewind.and.do_spinup) !,mixd_VEG)
           print *, 'call read_forcings ok - init'

          !* Really should call separate ent_prescribe_init with prescribed laimax
#ifdef ENT_STANDALONE_DIAG
      call loop_ent_cell_print_diag( iu_entcells, I0,I1,J0,J1,cells ) 
#endif
          call ent_prescribe_vegupdate(cells,hemi,jdaycount,year,
     &         do_giss_phenology=
     &         (.not.do_phenology_activegrowth).and.(.not.force_VEG),
     &         do_giss_albedo=.true.,
     &         do_giss_lai=
     &         (.not.do_phenology_activegrowth).and.(.not.force_VEG),
     &         update_crops=.false.,
!     &         mixed_veg=mixed_VEG,
     &         laidata=LAI,     ! pft LAI
     &         hdata=height, init=.true. ) ! pft height
          
           print *, 'call ent_prescribe_vegupdate passed okay - init'
#ifdef ENT_STANDALONE_DIAG
      call loop_ent_cell_print_diag( iu_entcells, I0,I1,J0,J1,cells ) 
#endif
           print *, 'LAI',LAI
           print *, 'HEIGHT',height
           print *,"ent_init_vegstruct passed ok"

          call rewind_files(.true., N_CASA_LAYERS)
        end if

      do while( time < max_time )

        print *,"---------------------------------------------------"
        print *,"started time step, time=", time, "dt=", dt
!        write(99,*)"---------------------------------------------------"
!        write(99,*)"started time step, time=", time, "dt=", dt

#ifdef DEBUG        
      call assign_dummyvals(I0,I1,J0,J1,N_DEPTH,
     &  TairC,TcanopyC,Qf,P_mbar,Ca,Ch,U,
     &  IPARdif,IPARdir,CosZen, Soilmoist, Soilmp,fice)
#endif        

        print *, 'calling read_forcings' 
        call read_forcings(I0, I1, J0, J1, 
     &     N_DEPTH, N_CASA_LAYERS,N_PFT,
     &     TairC, TcanopyC, Qf, P_mbar, Ca, Ch, U,
     &     IPARdif, IPARdir, CosZen,
     &     Soiltemp, Soilmoist, Soilmp, fice, LAI, height, 
     &     force_VEG, do_rewind.and.do_spinup) !,mixed_VEG)
        print *, 'call read_forcings ok' 
        print *, 'force_VEG =',force_VEG

#ifdef DEBUG        
        print *, 'forcings:', I0, I1, J0, J1, N_DEPTH,
     &          TairC,
     &          TcanopyC, Qf, P_mbar, Ca, Ch, U,
     &          IPARdif, IPARdir, CosZen,
     &          Soilmp(:,I0:I1,J0:J1), Soilmoist, Soiltemp, 
     &          fice(:,I0:I1,J0:J1),LAI(:,I0:I1,J0:J1),
     &          height(:,I0:I1,J0:J1)
#endif

#ifdef ENT_STANDALONE_DIAG
#ifdef PRINT_DRIVERS
        if (force_VEG) then
           write(980,'(5(i5),100(1pe16.8))') I0, I1, J0, J1, N_DEPTH,
     &          TairC,
     &          TcanopyC, Qf, P_mbar, Ca, Ch, U,
     &          IPARdif, IPARdir, CosZen,
     &          Soilmp(:,I0:I1,J0:J1), Soilmoist, Soiltemp, 
     &          fice(:,I0:I1,J0:J1),LAI(:,I0:I1,J0:J1),
     &          height(:,I0:I1,J0:J1)
        else
           write(980,'(5(i5),100(1pe16.8))') I0, I1, J0, J1, N_DEPTH,
     &          TairC,
     &          TcanopyC, Qf, P_mbar, Ca, Ch, U,
     &          IPARdif, IPARdir, CosZen,
     &          Soilmp(:,I0:I1,J0:J1), Soilmoist, Soiltemp, 
     &          fice(:,I0:I1,J0:J1)
        end if
#endif
#endif

        !* Set forcings

        !#### IF STARTING FROM RESTART, NEED TO INIALIZE ##
        !####     Ci, CNC, Tcan, Qf                      ##

        print *, 'calling ent_set_forcings'  
#ifdef DEBUG        
        call openunit('entcells_print',
     &       iu_entcells,.false.,.false.)
        call loop_ent_cell_print(iu_entcells,I0,I1,J0,J1,cells)
        !call closeunit(iu_entcells)
#endif

        call ent_set_forcings( cells,
     &       air_temperature=TairC,
     &       canopy_temperature=TcanopyC,
     &       canopy_air_humidity=Qf,
     &       surf_pressure=P_mbar,
     &       surf_CO2=Ca,
     &       heat_transfer_coef=Ch,
     &       wind_speed=U,
     &       total_visible_rad=IPARdir+IPARdif,
     &       direct_visible_rad=IPARdir,
     &       cos_solar_zenith_angle=CosZen,
     &       canopy_wet_fraction=fw,
     &       soil_temp=Soiltemp,  
     &       soil_moist=Soilmoist,
     &       soil_matric_pot=Soilmp,
     &       soil_ice_fraction=fice
     &       )
        print *, 'call ent_set_forcings ok'  

      !* NEW STREAMLINED CONTROL *!
      if (update_day) then
!          print *, 'Got here, updating veg.'
        if (force_VEG) then
          call ent_prescribe_vegupdate(cells,hemi,jdaycount,year,
     &         do_giss_phenology=
     &         (.not.do_phenology_activegrowth).and.(.not.force_VEG),
     &         do_giss_albedo=.true.,
     &         do_giss_lai=
     &         (.not.do_phenology_activegrowth).and.(.not.force_VEG),
     &         update_crops=.false.,
     &         laidata=LAI,     ! pft LAI
     &         hdata=height, init=.false. )   ! pft height
          print *,'ent_prescribe_vegupdate with LAI and height'
        else !Don't include laidata and hdata in parameter list.
          call ent_prescribe_vegupdate(cells,hemi,jdaycount,year,
     &         do_giss_phenology=
     &         (.not.do_phenology_activegrowth).and.(.not.force_VEG),
     &         do_giss_albedo=.true.,
     &         do_giss_lai=
     &         (.not.do_phenology_activegrowth).and.(.not.force_VEG),
     &         update_crops=.false.,init=.false.)
          print *,'ent_prescribe_vegupdate with calculated lai, ht'
        endif
      endif                     !update_day

      call ent_run(cells,dt,update_day) !Change name to ent_processes, calls ent_integrate
!      print *, 'call ent_run ok' 
 
 
      !***************************!

      !************************************************************************
!      if (.false.) then
!#ifdef PFT_MODEL_ENT
!        call ent_run (cells, dt, time)
!#endif
!        if (RUN_TYPE.eq.BIOPHYSICS_ONLY) then
!          if (update_day) then
!            if ( force_VEG ) then !if reading in external LAI  -PK 8/16/07 
!              call ent_prescribe_vegupdate(cells,hemi,jday,year,
!     &             do_giss_phenology=.true.,
!     &             do_giss_lai=(.not.force_VEG),
!     &             update_crops=.false.,
!     &             laidata=LAI, ! pft LAI
!     &             hdata=height ) ! pft height
!                                !jday = mod(jday,365) + 1  !Use this if jday is not read in with forcings.
!            else
!              call ent_prescribe_vegupdate(cells,hemi,jday,year,
!     &             do_giss_phenology=.true.,
!     &             do_giss_lai=(.not.force_VEG),
!     &             update_crops=.false.)
!                   !jday = mod(jday,365) + 1  !Use this if jday is not read in with forcings.
!            end if              !force_VEG
!          end if                !update_day
!          call ent_fast_processes( cells, dt )
!        else if (RUN_TYPE.eq.GISSVEG_CASA) THEN
!          if (update_day) then
!            if ( force_VEG ) then !if reading in external LAI  -PK 8/16/07 
!             call ent_prescribe_vegupdate(cells,hemi,jday,year,
!     &           do_giss_phenology=.true.,
!     &           do_giss_lai=(.not.force_VEG),
!     &           update_crops=.false.,
!     &           laidata=LAI,         ! pft LAI
!     &           hdata=height )       ! pft height
!                 !jday = mod(jday,365) + 1  !Use this if jday is not read in with forcings.
!            else
!                 call ent_prescribe_vegupdate(cells,hemi,jday,year,
!     &           do_giss_phenology=.true.,
!     &           do_giss_lai=(.not.force_VEG),
!     &           update_crops=.false.)
!                 !jday = mod(jday,365) + 1  !Use this if jday is not read in with forcings.
!            end if  !force_VEG
!           end if  !update_day
!          !call ent_run( cells, dt)
!          call ent_fast_processes( cells, dt )
!          
!          print *, 'call ent_run ok' !test -PK 7/11/06
!          print *, 'time=',time
!        end if
!
!      end if !.false 
      !************************************************************************

#ifdef DEBUG
        print *,"Got here before ent_get_exports."
#endif
        !* Extract results from Ent structure for GCM or diagnostics.
        call ent_get_exports( cells,
     &       canopy_conductance=GCANOPY,
     &       beta_soil_layers=betadl,
     &       shortwave_transmit=TRANS_SW,
     &       leafinternal_CO2=Ci,
     &       foliage_humidity=Qf,
     &       canopy_gpp=GPP,
     &       roughness_length=z0,
     &       flux_CO2=CO2flux,
     &       C_labile=C_labile,
     &       R_auto=R_auto,
     &       albedo=albedo
!     &       ,soilresp=Soil_resp
!     &       ,soilcpools=Tpool
     &     )

        !* Save results to a file.
        write(iu_results) GCANOPY,Ci,Qf,GPP,TRANS_SW,z0,CO2flux,
     &       betadl,albedo


#ifdef ENT_STANDALONE_DIAG
#ifdef PRINT_DRIVERS
        write(981,'(100e16.6)') GCANOPY,Ci,Qf,GPP,TRANS_SW,z0,CO2flux,
     &       betadl,albedo

!        print *,"Got here after iu_results."
!        print *,"canopy_conductance=",GCANOPY(I0,J0)
!        print *,"beta_soil_layers=",betadl(:,I0,J0)
!        print *,"shortwave_transmit=",TRANS_SW(I0,J0)
!        print *,"leafinternal_CO2=",Ci(I0,J0)
!        print *,"foliage_humidity=",Qf(I0,J0)
!        print *,"canopy_gpp=",GPP(I0,J0)
!        print *,"roughness_length=",z0(I0,J0)
!        print *,"flux_CO2=",CO2flux(I0,J0)
!        print *,"albedo=",albedo(:,I0,J0)
#endif
#endif


        !* Write ent state to a restart file.
        time_since_last_save = time_since_last_save + dt
        if ( time_since_last_save > save_interval ) then
          call ent_write_state( cells )
          time_since_last_save = 0.d0
        endif

        time = time + dt  !Moved time update to before check of update_day-NK
        update_day=(mod(time,86400.d0) .eq. 0.d0).and.(time.ne.0.d0)

        if (update_day) then
          jdaycount = jdaycount + 1
          if (jdaycount>365) then !numdays=365 for 1-year data, etc.
!YK - to accomodate the case in which the starting jday is not 1.
!          if (jdaycount-jday+1 > 365) then
            jdaycount = 1       !reset
            do_rewind =.true.   !so if do_spinup is true, forcing file will be reset.
          else
            do_rewind = .false.
          endif
          !jdaycount = mod(jdaycount,365) !mod(365) goes from 0 to 364 for days 365 and 1-364.  Use this if jday is not read in with forcings.
        endif

        
      enddo

      call closeunit(iu_results)
      !!call closeunit(iu_forcings)
      call close_forcings_file( force_VEG, N_CASA_LAYERS)!, mixed_VEG )

      end subroutine run_offline

!***************************************************************************
      

!***************************************************************************

      subroutine assign_dummyvals(I0,I1,J0,J1,N_DEPTH,
     &     TairC,
     &     TcanopyC,Qf,P_mbar,Ca,Ch,U,
     &     IPARdif,IPARdir,CosZen, Soilmoist, Soilmp,fice)
!! Dummy forcings (in the absence of real ones)
      integer, intent(in) :: I0,I1,J0,J1,N_DEPTH
      real*8,intent(out) :: TairC,TcanopyC,Qf,P_mbar,Ca,Ch,U
      real*8,intent(out) :: IPARdif,IPARdir,CosZen
      real*8, dimension(:,:,:) ::Soilmoist, Soilmp,fice

      TairC = 20.
      TcanopyC = 20.            ! TPJUN1950_datafile - canopy temperature
      Qf = 0.                   !Needs to be passed in from land surface.
      P_mbar = 900.             ! PSJUN1950_datafile - surface pressure
      Ca = 350.0e-6*(P_mbar*100)/(8.3145*298.15) !mol m-3
      Ch = 1.d-2                ! CHJUN1950_datafile - heat transf coeff
      U = 2.                    ! VSJUN1950_datafile -  V wind speed
                                ! USJUN1950_datafile - U wind speed
      IPARdif = 200.            ! VISJUN1950_datafile - total vis
      IPARdir = 50.             ! DVISJUN1950_datafile - direct vis
      CosZen = .5
      Soilmoist(:,:,:) = .01
      Soilmp(:,:,:) = -.1 ! MP1JUN1950_datafile - soil matric potential ???
                                ! should be 6 layers !
      fice(N_DEPTH,I0:I1,J0:J1) = 0. ! SIC1JUN1950_datafile - soil ice ????
                                ! should be 6 layers !
      end subroutine assign_dummyvals

!************************************************************************

  

