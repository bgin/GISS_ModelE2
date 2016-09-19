      module ent_struct_mod
!@sum Utilities for off-line run of Ent model.

      use ent_mod
      use FILEMANAGER
      use ent_types
      use ent_const
      use cohorts
      use patches
      use entcells
!      use ent_mod, only : ent_cell_pack

      implicit none

      contains

!************************************************************************
      !* Same subroutine as in ent_prog.f *!
      subroutine read_input_parameters(IM, JM, I0, I1, J0, J1,
     &     jday, year, jday2, year2, dt,i0f, i1f, j0f, j1f, 
     &     force_VEG, do_soilinit,
     &     do_soilresp, do_phenology, do_frost_hardiness,
     &     do_patchdynamics,skip,do_spinup)

      use filemanager
      implicit none
      integer IM, JM, I0, I1, J0, J1, i0f, i1f, j0f, j1f
      integer jday, year,jday2, year2
      real*8 dt
      logical force_VEG
      logical do_soilinit, do_soilresp, do_phenology, do_frost_hardiness
      logical do_patchdynamics, do_spinup
      integer skip !#HACK to skip records at beginning of forcing file
      !----
      integer iu_ent_input
      character*80, parameter :: file_ent_input="ent_input"
      namelist /input_parameters/ IM, JM, I0, I1, J0, J1
     &     ,jday, year, jday2, year2,dt
     &     ,i0f, i1f, j0f, j1f, force_VEG, do_soilinit
     &     ,do_soilresp, do_phenology, do_frost_hardiness 
     &     ,do_patchdynamics, skip, do_spinup

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

      subroutine ent_write_state( cells )
!@sum write ent state to the file
      type(entcelltype_public), intent(in) :: cells(:,:)
      !---
      real*8, pointer :: buffer(:)
      integer iu_entstate
      integer ic, jc, i, j

      ic = size(cells,1)
      jc = size(cells,2)

      write(*,*) 'ic, jc', ic, jc

      call openunit('ent_struct',iu_entstate,.true.,.false.)
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
      subroutine do_ent_struct_old (IM,JM,I0,I1,J0,J1 )
      implicit none
      integer :: IM,JM,I0,I1,J0,J1
      !---------
      type(entcelltype_public) cells(I0:I1,J0:J1)
      integer :: iu_entstruct
      integer :: iu_entstruct_ascii
      character*80, parameter :: ENTSTRUCT="ENTSTRUCT"
      integer :: N_CASAlayers    !N_CASA_layers
      integer :: i,j
      character :: next
      integer counter

      call ent_cell_construct( cells )

      !call ent_struct_setup (IM,JM,I0,I1,J0,J1,cells )
      ! This calls ent_struct_readcsv (IM,JM,I0,I1,J0,J1,cells )

!      open(iu_entstruct,file="ent_struct_init.csv",status='unknown')
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

      end subroutine do_ent_struct_old
!************************************************************************
      subroutine do_ent_struct (IM,JM,I0,I1,J0,J1 )
      implicit none
      integer :: IM,JM,I0,I1,J0,J1
      !---------
      type(entcelltype_public) cells(I0:I1,J0:J1)
      integer :: iu_entstruct
      integer :: iu_entstruct_ascii
      character*80, parameter :: ENTSTRUCT="ENTSTRUCT"
      integer :: N_CASAlayers    !N_CASA_layers
      integer :: i,j
      character :: next
      integer counter

      call ent_cell_construct( cells )

!      open(iu_entstruct,file="ent_struct_init.csv",status='unknown')
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

      end subroutine do_ent_struct
!!************************************************************************

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
!************************************************************************

      end module ent_struct_mod


!************************************************************************
      !program ent_struct
      subroutine Ent_struct
!@sum ent_struct
!@+   Reads in a text file description of Ent vegetation data structure
!@+   with entcells, patches, and mixed cohorts, and generates an
!@+   entcell array data structure.  Writes the data structure to a restart.
      use ent_struct_mod
      use ent_mod, only: ent_struct_setup

      implicit none
      integer :: IM, JM, I0,I1,J0,J1,i0f, i1f, j0f, j1f
      integer :: jday, year,jday2, year2
      real*8  :: dt             !seconds
      logical :: force_VEG
      logical :: do_soilinit, do_soilresp, do_phenology
      logical :: do_frost_hardiness, do_patchdynamics, do_spinup
      integer :: skip !#HACK to skip records at beginning of forcing file

      !* Default configuration
      dt = 1800.d0    !seconds. Default value may be over-ridden by ent_input
      force_VEG = .false.
      do_soilinit = .true.
      do_soilresp = .true.
      do_phenology = .false.
      do_frost_hardiness = .true.
      do_patchdynamics = .false.
      do_spinup = .false.
      skip = 0

      write(*,*) 'Got here before read_input_parameters in ent_struct.f'
       call read_input_parameters(IM, JM, I0, I1, J0, J1
     &     ,jday, year, jday2, year2,dt
     &     ,i0f, i1f, j0f, j1f
     &     ,force_VEG
     &     ,do_soilinit, do_soilresp, do_phenology, do_frost_hardiness
     &     ,do_patchdynamics, skip, do_spinup)
      if (year2.eq.-1) year2=year  !If year2 was not specified, default same as year.

      print *,"ent_input: "
     &     , jday, year, jday2, year2,dt
     &     , i0f, i1f, j0f, j1f
     &     , force_VEG
     &     , do_soilinit, do_soilresp, do_phenology, do_frost_hardiness
     &     , do_patchdynamics, do_spinup


      call do_ent_struct(IM,JM,I0,I1,J0,J1)

      write(*,*) 'Done reading and saving ent data structure.'
      !end program ent_struct
      end subroutine ent_struct


