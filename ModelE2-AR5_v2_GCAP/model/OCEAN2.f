! for the moment, this file is included in other files
!#include "rundeck_opts.h"
!#ifdef TRACERS_ATM_ONLY
!#undef TRACERS_ON
!#undef TRACERS_WATER
!#endif

      subroutine init_ocean(iniocean,istart)
!@sum init_OCEAN initializes ocean variables
!@auth Original Development Team
!@ver  1.0
      use domain_decomp_atm, only : grid, get
      use model_com, only : kocean,focean,ioread
#ifdef TRACERS_WATER
      use tracer_com, only : trw0
      use fluxes, only : gtracer
#endif
      use fluxes, only : sss,uosurf,vosurf,uisurf,visurf,ogeoza
      use seaice, only : qsfix, osurf_tilt
      use ocnml, only : init_ocnml,set_gtemp_ocnml
      use sstmod, only : init_sstmod,set_gtemp_sst
      use pario, only : par_open, par_close
      implicit none
      logical, intent(in) :: iniocean  ! true if starting from ic.
      integer, intent(in) :: istart
!@var sss0 default sea surface salinity (psu)
      real*8, parameter :: sss0=34.7d0
      integer :: fid
      integer :: i,j
      integer :: i_0,i_1, j_0,j_1

      call get(grid,I_STRT=I_0,I_STOP=I_1,j_strt=j_0,j_stop=j_1)

      if (istart.le.0) then
        if(kocean.ge.1) call init_ODEEP(.false.)
        return
      end if

C**** Cold start
      if (istart.le.2) then
        fid = par_open(grid,'GIC','read')
        call new_io_ocean (fid,ioread)
        call par_close(grid,fid)
      end if

c**** set fluxed arrays for oceans
      do j=j_0,j_1
      do i=i_0,i_1
        if (focean(i,j).gt.0) then
          sss(i,j) = sss0
#ifdef TRACERS_WATER
          gtracer(:,1,i,j)=trw0(:)
#endif
        else
          sss(i,j) = 0.
        end if
c**** for the time being assume zero surface velocities for drag calc
        uosurf(i,j)=0. ; vosurf(i,j)=0.
        uisurf(i,j)=0. ; visurf(i,j)=0.
c**** also zero out surface height variations
        ogeoza(i,j)=0.
      end do
      end do
c**** keep salinity in sea ice constant for fixed-sst and qflux models
      qsfix = .true.
c**** make sure to use geostrophy for ocean tilt term in ice dynamics
c**** (if required). since ocean currents are zero, this implies no sea
c**** surface tilt term.
      osurf_tilt = 0


C**** 
      if (kocean.eq.0) then
        call set_gtemp_sst
        call init_sstmod
      else
        call set_gtemp_ocnml
        ! read ML depths and OHT
        call init_ocnml(iniOCEAN,istart)
      endif

      return
      end subroutine init_ocean

      subroutine alloc_ocean(grid)
!@sum alloc_ocean calls allocation routines for either the
!@+     (1) prescribed ocean module (kocean=0)
!@+     (2) mixed-layer ocean module (kocean=1)
      use domain_decomp_atm, only : dist_grid
      use sstmod, only  : alloc_sstmod
      use ocnml, only : alloc_ocnml
      use param, only : get_param
      use model_com, only : kocean
      implicit none
      type (dist_grid), intent(in) :: grid
      call get_param('kocean',kocean)
      if(kocean.eq.0) then
        call alloc_sstmod
      else
        call alloc_ocnml
        call alloc_odeep(grid)
      endif
      end subroutine alloc_ocean

      subroutine daily_ocean(end_of_day)
!@sum daily_ocean calls daily update routines for either the
!@+     (1) prescribed ocean module (kocean=0)
!@+     (2) mixed-layer ocean module (kocean=1)
      use model_com, only : kocean
      use sstmod, only : set_gtemp_sst
      use ocnml, only : daily_ocnml,set_gtemp_ocnml
      implicit none
      logical, intent(in) :: end_of_day

      if (kocean.ge.1) then
        ! update prescribed ml depth, perform associated adjustments
        call daily_ocnml(end_of_day)
        call set_gtemp_ocnml
      else
        ! update prescribed sst
        call read_sst(end_of_day)
        call set_gtemp_sst
      end if

      return
      end subroutine daily_ocean

      subroutine oceans
!@sum ocean calls routines to apply surface fluxes to either the
!@+     (1) prescribed ocean (kocean=0, no-op)
!@+     (2) mixed-layer ocean (kocean=1)
!@auth original development team
!@ver  1.0
      use model_com, only : focean,kocean
      use geom, only : imaxj,axyp
      use diag_com, only : oa
      use fluxes, only : dmsi,dhsi,dssi
#ifdef TRACERS_WATER
     *     ,dtrsi
#endif
      use fluxes, only : eflowo,egmelt
      use domain_decomp_atm, only : grid,get
      use ocnml, only : run_ocnml,set_gtemp_ocnml
      implicit none

      integer i,j, j_0,j_1, i_0,i_1

      call get(grid,i_strt=i_0,i_stop=i_1,j_strt=j_0,j_stop=j_1)

      if(kocean.ge.1) then
        ! surface fluxes affect predicted ocean temperature
        call run_ocnml
        call set_gtemp_ocnml
      else
        ! add rvr e to surf. energy budget, set ice formation rate = 0
        do j=j_0,j_1
        do i=i_0,imaxj(j)
          if (focean(i,j).gt.0) then
            oa(i,j,4)=oa(i,j,4)+
     &         (eflowo(i,j)+egmelt(i,j))/(focean(i,j)*axyp(i,j))
            dmsi(1,i,j)=0.
            dmsi(2,i,j)=0.
            dhsi(1,i,j)=0.
            dhsi(2,i,j)=0.
            dssi(1,i,j)=0.
            dssi(2,I,J)=0.
#ifdef TRACERS_WATER
            dtrsi(:,1,I,J)=0.
            dtrsi(:,2,I,J)=0.
#endif
          end if
        end do
        end do
      endif
      return
      end subroutine oceans

      subroutine precip_oc
!@sum  precip_oc driver for applying precipitation fluxes to mixed-layer ocean.
!@+    This routine could be folded into oceans, but exists separately for
!@+    historical/diagnostic reasons.
!@auth original development team
!@ver  1.0
      use diag_com, only : oa
      use fluxes, only : eprec
      use model_com, only : kocean,focean
      use ocnml, only : precip_ocnml,set_gtemp_ocnml
      implicit none
      where(focean.gt.0.) oa(:,:,4) = oa(:,:,4)+eprec(:,:)
      if(kocean.ge.1) then
        call precip_ocnml
        call set_gtemp_ocnml
      endif
      return
      end subroutine precip_oc

      subroutine def_rsf_ocean(fid)
!@sum  def_rsf_ocean defines ocean array structure in restart files
!@auth M. Kelley
!@ver  beta
      use model_com, only : kocean
      use sstmod, only : def_rsf_sstmod
      !use ocnml, only : def_rsf_ocnml
      use domain_decomp_atm, only : grid
      implicit none
      integer fid   !@var fid file id
      if(kocean.eq.0) then
        call def_rsf_sstmod(fid)
      else
        call def_rsf_ocnml(fid)
      endif
      return
      end subroutine def_rsf_ocean

      subroutine new_io_ocean(fid,iaction)
!@sum  new_io_ocean read/write ocean arrays from/to restart files
!@auth M. Kelley
!@ver  beta new_ prefix avoids name clash with the default version
      use model_com, only : kocean
      use sstmod, only : new_io_sstmod
      !use ocnml, only : new_io_ocnml
      implicit none
      integer fid   !@var fid unit number of read/write
      integer iaction !@var iaction flag for reading or writing to file
      if(kocean.eq.0) then
        call new_io_sstmod(fid,iaction)
      else
        call new_io_ocnml(fid,iaction)
      endif
      return
      end subroutine new_io_ocean
