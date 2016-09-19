! for the moment, this file is included in other files
!#include "rundeck_opts.h"

      module sstmod

!@sum  Module sstmod contains the arrays/subroutines needed to prescribe
!@+    ocean surface temperature from input files.
!@auth Original Development Team
!@auth M. Kelley restructuring and netcdf-based input options

      use timestream_mod, only : timestream
      implicit none
      save

!@var sst sea surface temperature (C)
      real*8, dimension(:,:), allocatable :: sst

!@var SSTstream interface for reading and time-interpolating SST files
!@+   See general usage notes in timestream_mod.
!@+   Note regarding comparison of results to runs that use traditional I/O:
!@+   if SST datafiles contain monthly means and the piecewise parabolic
!@+   method is used for monthly->daily interpolation, the presence of OSST_eom
!@+   in the rundeck will prompt read_stream to read end-of-month values from
!@+   OSST_eom rather than computing them on the fly.  OSST and OSST_eom may
!@+   refer to the same file or directory.  The on-the-fly result will differ
!@+   due to roundoff effects.
      type(timestream) :: SSTstream

!@var tocean_4io an array for restart file compatibility with ML ocean
!@+   (see OCNML2.f for definition of its tocean array)
      real*8, dimension(:,:,:), allocatable :: tocean_4io

      contains

      subroutine alloc_sstmod
!@sum alloc_sstmod allocates arrays in module sstmod
      use domain_decomp_atm, only : grid,get
      implicit none
      integer :: i_0h,i_1h,j_0h,j_1h,ier
      call get(grid,j_strt_halo=j_0h,j_stop_halo=j_1h)
      i_0h = grid%i_strt_halo
      i_1h = grid%i_stop_halo
      allocate(sst(i_0h:i_1h,j_0h:j_1h))
      allocate(tocean_4io(3,i_0h:i_1h,j_0h:j_1h))
      end subroutine alloc_sstmod

      subroutine init_sstmod
!@sum init_sstmod initializes the SSTstream object
      use domain_decomp_atm, only : grid,get
      use model_com, only : focean,jyear,jday
      use timestream_mod, only : init_stream
      implicit none
      call init_stream(grid,SSTstream,'OSST','sst',-100d0,100d0,'ppm',
     &       jyear,jday,msk=focean)
      end subroutine init_sstmod

      subroutine set_gtemp_sst
!@sum set_gtemp_sst copies sst into the ocean position in the gtemp array
      use constant, only : tf
      use domain_decomp_atm, only : grid,get
      use fluxes, only : gtemp,gtempr
      use geom, only : imaxj
      use model_com, only : focean
#ifdef SCM
      USE MODEL_COM, only : I_TARG,J_TARG
      USE SCMCOM, only : iu_scm_prt,SCM_SURFACE_FLAG,ATSKIN
#endif
      implicit none
      integer :: i,j,j_0,j_1, i_0,i_1
      call get(grid,i_strt=i_0,i_stop=i_1,j_strt=j_0,j_stop=j_1)
      do j=j_0,j_1
      do i=i_0,imaxj(j)
        if (focean(i,j).gt.0) then
          gtemp(1,1,i,j)=sst(i,j)
          gtemp(2,1,i,j)=tocean_4io(2,i,j) ! to preserve identical diagnostics
          gtempr(1,i,j) =sst(i,j)+tf
#ifdef SCM
c         keep ocean temp fixed for SCM case where surface
c         temp is supplied
          if (I.eq.I_TARG.and.J.eq.J_TARG) then
            if (SCM_SURFACE_FLAG.ge.1) then
              GTEMP(1,1,I,J) = ATSKIN
              GTEMPR(1,I,J) = ATSKIN + TF
            endif
          endif
#endif
        endif
      enddo
      enddo
      end subroutine set_gtemp_sst

      subroutine def_rsf_sstmod(fid)
!@sum  def_rsf_sstmod defines sstmod array structure in restart files
!@auth M. Kelley
!@ver  beta
      use domain_decomp_atm, only : grid
      use pario, only : defvar
      implicit none
      integer fid   !@var fid file id
      call defvar(grid,fid,tocean_4io,'tocean(d3,dist_im,dist_jm)')
      return
      end subroutine def_rsf_sstmod

      subroutine new_io_sstmod(fid,iaction)
!@sum  new_io_sstmod read/write sstmod arrays from/to restart files
!@auth M. Kelley
!@ver  beta new_ prefix avoids name clash with the default version
      use model_com, only : ioread,iowrite
      use domain_decomp_atm, only : grid
      use pario, only : write_dist_data,read_dist_data
      implicit none
      integer fid   !@var fid unit number of read/write
      integer iaction !@var iaction flag for reading or writing to file
      select case (iaction)
      case (iowrite)            ! output to restart file
        tocean_4io(1,:,:) = sst(:,:)
        call write_dist_data(grid, fid, 'tocean', tocean_4io, jdim=3)
      case (ioread)            ! input from restart file
        call read_dist_data(grid, fid, 'tocean', tocean_4io, jdim=3)
        sst(:,:) = tocean_4io(1,:,:)
      end select
      return
      end subroutine new_io_sstmod

      end module sstmod

      subroutine read_sst(end_of_day)
!@sum read_sst invokes procedures to read sea surface temperature from
!@+   input files and perform time interpolation
      use domain_decomp_atm, only : get,grid
      use model_com, only : im,jm,focean,itime,itimei,jyear,jday
      use geom, only : imaxj
      use seaice, only : tfrez
      use fluxes, only : sss
      use timestream_mod, only : read_stream
      use sstmod, only : SSTstream,SST
      implicit none

      logical, intent(in) :: end_of_day
      real*8 :: tfo
      integer i,j

      integer :: j_0,j_1, i_0,i_1
      logical :: have_north_pole, have_south_pole

      call get(grid,i_strt=i_0,i_stop=i_1,j_strt=j_0,j_stop=j_1,
     &         have_south_pole=have_south_pole,
     &         have_north_pole=have_north_pole)

      if(.not.(end_of_day.or.itime.eq.itimei)) return

C**** read and time-interpolate
      call read_stream(grid,SSTstream,jyear,jday,SST)

c**** bounds on sst
      do j=j_0,j_1
      do i=i_0,imaxj(j)
        if (focean(i,j).gt.0) then
          tfo=tfrez(sss(i,j))
          if (sst(i,j).lt.tfo) sst(i,j)=tfo
        else
          sst(i,j) = 0.
        endif
      end do
      end do

c**** replicate values at pole
      if(have_north_pole) then
        if (focean(1,jm).gt.0) then
          do i=2,im
            sst(i,jm)=sst(1,jm)
          end do
        end if
      end if
      if(have_south_pole) then
        if (focean(1,1).gt.0) then
          do i=2,im
            sst(i,1)=sst(1,1)
          end do
        end if
      end if

      return
      end subroutine read_sst
