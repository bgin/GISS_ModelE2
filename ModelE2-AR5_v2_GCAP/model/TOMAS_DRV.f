#include "rundeck_opts.h"

!@sum  TOMAS_DRV: TwO-Moment Aerosol Sectional (TOMAS) microphysics driver 
!@+     aerosol microphysics (nucleation,coagulation, condensation) and 
!@+     SO4 formation from (clouds) aqueous chemistry. 
!@auth Peter Adams/Jeff Pierce/Yunha Lee (implemented into ModelE by Yunha Lee)
!@ver  1.0
!@calls various subroutines under TOMAS_DRV.f and TOMAS_microphysics.f 

      MODULE TOMAS_AEROSOL

      USE MODEL_COM, only : im,jm,lm     ! dimensions
      USE TRACER_COM, only : ntm,nbins
      IMPLICIT NONE 

C-----INCLUDE FILES--------------------------------------------------
!@param ibins : number of size bins used in TOMAS (equal to nbins defined in TRACER_COM) 
!@param icomp : number of size-resolved chemical species (SO4,SS,ECOB,ECIL,OCOB,OCIL,DUST,NH4,AER-WATER
!@param idiag : number of aerosol diagnostic species (NH4 and AER-WATER) 
      integer,parameter :: ibins=nbins
      integer,parameter :: icomp=9
      integer,parameter :: idiag=2
!@param ncomp : number of size-resolved chemical species that actually used in GCM 
!@+    NH4 is not size-resolved out of TOMAS algorithm, which makes ncomp differs from icomp
      integer,parameter ::  ncomp=8
!@param srtso4,srtna, etc : index number for icomp 
      integer srtso4, srtna, srth2o, srtecob, srtecil, srtocob,
     &	      srtocil, srtnh4, srtdust
      parameter (srtso4=1,
     &           srtna =2,
     &           srtecob=3,
     &           srtecil=4,
     &           srtocob=5,
     &           srtocil=6,
     &		 srtdust=7,
     &           srtnh4=8,
     &           srth2o=9)

!@var Nk and Mk contain the number and mass size distributions of the
!@+   aerosol.  Units are #/grid cell or kg/grid cell, respectively.
!@+   Nkd and Mkd store values of Nk and Mk for diagnostic purposes.
!@var Gc are gas phase concentrations (kg/grid cell) of species
!@+   corresponding to the aerosol species (e.g. H2SO4 for sulfate).

      real*8 Nk(ibins), Mk(ibins,icomp), Gc(icomp-1)
      real*8 Nkd(ibins), Mkd(ibins,icomp), Gcd(icomp-1) 

C The following variables describe the grid cell in which the
C microphysics is operating.
!@var boxvol : volume of grid cell (cm3)
!@var boxmass : air mass of grid cell (kg)
!@var temp : temperature (K) of grid cell
!@var pres : air pressure (Pa) of grid cell 
!@var rh : relative humidity (0-1)
      real*8 boxvol    
      real*8 boxmass    
      real*8 temp      
      real*8 pres       
      real*8 rh        

C Physical properties of aerosol components
!@param molwt : molecular weight of chemical species (icomp) 
      real molwt(icomp)
      data molwt/96., 58.45, 200., 200., 200., 200., 100.,18.,18./

!@param bin_nuc/tern_nuc/ion_nuc/actv_nuc : Flag for which nucleation parameterizations to use (1=on)
      integer bin_nuc, tern_nuc, ion_nuc, actv_nuc   
      parameter(bin_nuc=1, tern_nuc=0, ion_nuc=0, actv_nuc=0) 

!@var soa_amp : mass growth amplification factor (determined by the 
!@+             amount of soa that needs to be condensed
      real*8 soa_amp

!@param tau_soa : 1st order timescale in which SOA condenses (0.5 days)
      real*8 tau_soa
      parameter(tau_soa=0.5d0)

!@var SOArate is the rate of SOA condensing (kg/s)
!      real*8, ALLOCATABLE, DIMENSION(:,:)  :: SOA_chem !SOA formation [kg]
      real*8  SOArate
!@var surf_area : aerosol surface area [micon^2 cm^-3]
      real*8 surf_area 
!@var ionrate : ion pair formation rate [ion pairs cm^-3 s^-1]
      real*8 ionrate 
!@var binact10/binact02 : lookup table of activated size bin at supersatuaration of 1.0% and 0.2%
      integer, dimension(101,101,101) :: binact10,binact02
!@var fraction10/fraction02 : lookup table of chemical composition fraction in the activated size bin at supersatuaration of 1.0% and 0.2%
      real*8,dimension(101,101,101) :: fraction10,fraction02
!@param ptype : number of aerosol microphysics process 
      integer, parameter :: ptype=7
#ifdef TOMAS_HETCHEM
     *    +1
#endif
!@var AEROD : saving aerosol microphyiscs dianogstics
      real*8, ALLOCATABLE,dimension(:,:,:,:,:) :: AEROD
!@var AQSO4oxid_mc/AQSO4oxid_ls : aqueous h2so4 formation for convective and large-scale clouds 
      real*8, ALLOCATABLE,DIMENSION(:,:,:) :: AQSO4oxid_mc,AQSO4oxid_ls 
!@var h2so4_chem : h2so4 formation rate from so2+oh [kg of H2SO4/sec]
      real*8, ALLOCATABLE,DIMENSION(:,:,:)  ::  h2so4_chem
#ifdef TOMAS_HETCHEM
!@var h2so4_chem : h2so4 formation rate from so4 hetchem [kg of H2SO4/sec]
      real*8, ALLOCATABLE,DIMENSION(:,:,:,:)  ::  h2so4_hetchem
#endif
!@var N_subgridcg : aerosol number emission rate changed by subgrid coagulation 
      real*8, ALLOCATABLE,DIMENSION(:,:,:,:,:) :: N_subgridcg 
!@var M_subgridcg : aerosol mass emission rate changed by subgrid coagulation
      real*8, ALLOCATABLE,DIMENSION(:,:,:,:,:,:)  :: M_subgridcg 
!@var trm_emis : TRM before emission and used in subgridcoagualtion process 
      real*8, ALLOCATABLE,DIMENSION(:,:,:,:)  :: trm_emis
!@var CCN_TOMAS [CM-3]
      real*8, ALLOCATABLE,DIMENSION(:,:,:,:)  :: CCN_TOMAS
!@var TOMAS_QEXT/TOMAS_QSCA/TOMAS_QABS/TOMAS_GSCA : size-dependant radiative properties 
!@+      lookup tables based on Mie theory
      REAL*8, DIMENSION(124,101,91)      :: TOMAS_QEXT, TOMAS_QSCA,
     &     TOMAS_QABS,TOMAS_GSCA !,TOMAS_QBACK 

!@param TOMAS_DIAG_FC : Flag used in aerosol radiation calculation
!@+  2=external mixing (=icomp-2) radiation calls  |
!@+  1=internal mixing (but AECOB is externally mixed) (ANUM_01) radiation call
!@+  TOMAS_DIAG_FC=2 is only available now.

      INTEGER                            :: TOMAS_DIAG_FC = 2 

!@var number of supersaturations (0.1/0.2/0.3) 
      integer, parameter :: nsmax=3
      real*8, parameter :: Smax(nsmax)=(/0.1,0.2,0.3/)

  
      END MODULE TOMAS_AEROSOL


      SUBROUTINE TOMAS_DRV 

      USE DOMAIN_DECOMP_ATM, only : GRID, GET, write_parallel
     &     ,am_i_root
      USE TOMAS_AEROSOL 
      USE TRACER_COM

      USE TRDIAG_COM, only : taijs=>taijs_loc,taijls=>taijls_loc
     *     ,ijts_TOMAS,itcon_TOMAS
      USE FLUXES, only: tr3Dsource

      USE MODEL_COM, only : im,jm,lm     ! dimensions
     $                     ,t            ! potential temperature (C)
     $                     ,q            ! saturated pressure
     $                     ,dtsrc
      USE GEOM, only: axyp,imaxj,BYAXYP
      USE CONSTANT,   only:  lhe,mair,gasc   
      USE DYNAMICS,   only: pmid,pk,byam,gz, am   ! midpoint pressure in hPa (mb)
!                                           and pk is t mess up factor
!                                           BYAM  1/Air mass (m^2/kg)
      IMPLICIT NONE

C-----VARIABLE DECLARATIONS------------------------------------------

      INTEGER J_0, J_1, I_0, I_1

      integer i,j,l,n,jc,mt,k,np  !counters
      integer mpnum       !microphysical process id #
      real adt            !aerosol microphysics time step (seconds)
      real*8 qsat         !used in RH calculation
      integer tracnum
      integer flag
      real*8 frac
      
      real*8 Nkout(NBINS), Mkout(NBINS,icomp)
      real*8 Gcout(icomp-1)
!@var tot_aam : total aerosol ammonia per grid cell across all bins
      real*8 tot_aam  
!@var Gcavg : average h2so4 concentration during timstep
      real*8 Gcavg 
             
!@var H2SO4rate_o : H2SO4rate for the specific gridcell
!@var SOAmass : SOA mass in a gridcell
      real*8 H2SO4rate_o, SOAmass 
      integer num_iter
!@var fn : nucleation rate of clusters cm-3 s-1
      real*8 fn                 
!@var fn1: formation rate of particles to first size bin cm-3 s-1
      real*8 fn1               
      real*8 tot_n_1, tot_n_1a, tot_n_2, tot_n_i ! used for nitrogen mass checks
      real*8 tot_s_1, tot_s_1a,tot_s_1b, tot_s_2 ! used for sulfur mass checks
      real*8 Nknuc(ibins), Mknuc(ibins, icomp)
      real*8 Nkcond(ibins),Mkcond(ibins,icomp)
      real*8 INIT_Nk(ibins),INIT_Mk(ibins,icomp)
      real*8 INIT_H2SO4,INIT_NH3,INIT_NH4,INIT_SOA
      real*8 TSUM(2)

c$$$      real*4, dimension(GRID%J_STRT_HALO:GRID%J_STOP_HALO,lm) :: 
c$$$     &     nucrate,nucrate1

C-----CODE-----------------------------------------------------------    
C****
C**** Extract useful local domain parameters from "grid"
C****
      CALL GET(grid, J_STRT=J_0,       J_STOP=J_1)
      I_0 = grid%I_STRT
      I_1 = grid%I_STOP

!debug      nucrate(J_0:J_1,I_0:I_1)        = 0.d0  !DIAG: nucleation rate diagnotics (Jnuc)
!debug      nucrate1(J_0:J_1,I_0:I_1)       = 0.d0  !DIAG: new particle formation rate at lowest boundary in TOMAS


C     Loop over all grid cells
      DO L=1,LM                            
         DO J=J_0,J_1                          
            DO I=I_0,IMAXJ(J)


               temp = pk(l,i,j)*t(i,j,l) !should be in [K]
               rh = MIN(1.,q(i,j,l)/QSAT(temp,lhe,pmid(l,i,j))) ! rH [0-1]
               pres= pmid(l,i,j)*100. ! pmid in [hPa]
               boxmass=am(l,i,j)*axyp(i,j) !kg of air
               boxvol=boxmass/mair*1000.d0
     &              *gasc*temp/pres*1e6 !cm3
                          
Cjrp  Initialize all components condensible gas values to zero      
Cjrp  Gc(srtso4) will remain zero until within cond_nuc where the
Cjrp  pseudo steady state H2SO4 concentration will be put in this place.

               Gc(:)=0.d0
               
C     Swap T0M into Nk, Mk, Gc arrays

               do n=1,ibins
                  Nk(n)=TRM(i,j,l,IDTNUMD-1+n)
                  Mk(n,srtso4)=TRM(i,j,l,IDTSO4-1+n)
                  Mk(n,srtna) =TRM(i,j,l,IDTNA -1+n)
                  MK(n,srtecob)=TRM(i,j,l,IDTECOB -1+n)
                  MK(n,srtecil)=TRM(i,j,l,IDTECIL -1+n)
                  MK(n,srtocob)=TRM(i,j,l,IDTOCOB -1+n)
                  MK(n,srtocil)=TRM(i,j,l,IDTOCIL -1+n)      
                  Mk(n,srtdust)=TRM(i,j,l,IDTDUST -1+n)            
                  Mk(n,srth2o)=TRM(i,j,l,IDTH2O-1+n)
                  Mk(n,srtnh4)=0.
               enddo

               INIT_NK(:) = NK(:)
               INIT_MK(:,:)=MK(:,:)
               INIT_H2SO4 = H2SO4_chem(I,J,L)*dtsrc
               INIT_NH3=TRM(I,J,L,n_NH3)
               INIT_NH4=TRM(I,J,L,n_NH4)
               INIT_SOA=TRM(I,J,L,n_SOAgas)

! swap NH3 from giss to tomas               
               tot_n_i = TRM(i,j,l,n_NH3)*14.d0/17.d0 + 
     &              TRM(i,j,l,n_NH4)*14.d0/18.d0

               call NH3_GISStoTOMAS(TRM(i,j,l,n_NH3), 
     &              TRM(i,j,l,n_NH4),Gc,Mk)
              
                                ! nitrogen and sulfur mass checks
                                ! get the total mass of N
               tot_n_1 = Gc(srtnh4)*14.d0/17.d0
               do k=1,ibins
                  tot_n_1 = tot_n_1 + Mk(k,srtnh4)*14.d0/18.d0
               enddo
               
               H2SO4rate_o = H2SO4_chem(i,j,l) !kg of h2so4/sec  (from SO2+OH)
               SOAmass=TRM(i,j,l,n_SOAgas) !kg of SOA 

! get the total mass of S
               tot_s_1 = H2SO4rate_o*dtsrc*32.d0/98.d0
               do k=1,ibins
                  tot_s_1 = tot_s_1 + Mk(k,srtso4)*32.d0/96.d0
#ifdef TOMAS_HETCHEM
     &              +sum(H2SO4_hetchem(i,j,l,1:ibins))*dtsrc*32.d0/96.d0
#endif
                enddo
! timestep is shorten to be 10 mins 

               adt=dtsrc/3. 

               do mt=1,3        ! 3 * 10 min inside

                                !calculate SOA to condense
               SOArate = SOAmass*(1.d0-
     &              exp(-adt/(tau_soa*3600.*24.)))/adt

               SOAmass=SOAmass-SOArate*adt ! updated SOA mass within the time step loop

               if(H2SO4rate_o.le.0.) THEN
                  if(am_i_root())
     &                 print*,'problem in soa_amp',i,j,l,
     &                 h2so4rate_o,SOArate

                  soa_amp=0. !no SOA condensing! 
               else

               soa_amp = SOArate/H2SO4rate_o  ! TOMAS- floating invalid due to zero division??
               endif

C     ****************
C     Aerosol dynamics
C     ****************                                  


!MAR2013 - CALL MNFIX BEFORE COND_NUC??!
               call ezwatereqm(Mk)
               call storenm()
               call mnfix(Nk,Mk)
               mpnum=6 
               call aerodiag(mpnum,i,j,l)

               Gcavg = 0.0

               Gc(srtso4)=h2so4rate_o*adt ! this is for condensation diagnostics 

               call storenm()

               Gc(srtso4)=0.

C If any Nk are zero, then set them to a small value to avoid division by zero
               call cond_nuc(Nk,Mk,Gc,Nkout,Mkout,Gcout,fn,fn1,
     &             H2SO4rate_o,adt,num_iter,Nknuc,Mknuc,Nkcond,Mkcond,l)
                                !get nucleation diagnostic
 
               Mk(:,:)=Mknuc(:,:)
               Nk(:)=Nknuc(:)

               Gc(srtso4)=h2so4rate_o*adt !to make zero nucleation diag
            
               mpnum=3 
               call aerodiag(mpnum,i,j,l)

               Mk(:,:)=Mkcond(:,:)
               Nk(:)=Nkcond(:)

               
               Gc(srtnh4)=Gcout(srtnh4)
               Gc(srtso4)=Gcout(srtso4)

ccc               TRM(I,J,L,n_H2SO4)=Gc(srtso4)

               mpnum=1
               call aerodiag(mpnum,i,j,l)
               Mk(:,:)=Mkout(:,:)
               Nk(:)=Nkout(:)

!               nucrate(j,l)=nucrate(j,l)+fn
!               nucrate1(j,l)=nucrate1(j,l)+fn1
               
                                ! accumulate nucleation rate diagnostics
                                ! first sum for JL
!               TSUM(1)=TSUM(1)+fn*boxvol*adt ! number of particles generated per kg of air in timestep
!               TSUM(2)=TSUM(2)+fn1*boxvol*adt
                                ! IJ
c$$$               T3DC(I,J,L,1)=T3DC(I,J,L,1)+fn*boxvol*adt
c$$$               T3DC(I,J,L,2)=T3DC(I,J,L,2)+fn1*boxvol*adt      
               
                                ! nitrogen and sulfur mass checks
                                ! get the total mass of N
               tot_n_1a = Gc(srtnh4)*14.d0/17.d0
               do k=1,ibins
                  tot_n_1a = tot_n_1a + Mk(k,srtnh4)*14.d0/18.d0
               enddo
               
                                ! get the total mass of S
               tot_s_1b = 0.d0
               do k=1,ibins
                  tot_s_1b = tot_s_1b + Mk(k,srtso4)*32.d0/96.d0
               enddo               
        
               call storenm()
               call mnfix(Nk,Mk)
               mpnum=6 
               call aerodiag(mpnum,i,j,l)


C YHL - This should be called for next multicoag.  This does not need to be called in last timestep but do it anyway. 

                                !Coagulation
               call storenm()
               call multicoag(adt)

               mpnum=2
               call aerodiag(mpnum,i,j,l)

C     Do water eqm at appropriate times
               call eznh3eqm(Gc,Mk)
               call ezwatereqm(Mk)
                  
C     ***********************
C     End of aerosol dynamics
C     ***********************
              

             enddo              ! timestep


#ifdef TOMAS_HETCHEM
               call storenm()
               call mnfix(Nk,Mk)
               mpnum=6 
               call aerodiag(mpnum,i,j,l)

c             if(sum(h2so4_hetchem(i,j,l,1:ibins)).gt.0.and.i.eq.30)then
c               do k=1,ibins
c               print*,'debug bef hetchem',k,H2SO4_HETCHEM(I,J,L,k)
c     &               ,Nk(k) ,Mk(k,1),Mk(k,2)
c               enddo
c             endif

               call storenm()
                  call cond_hetchem(i,j,l) ! SO4 hetchem 


c             if(sum(h2so4_hetchem(i,j,l,1:ibins)).gt.0.and.i.eq.30)then
c              do k=1,ibins
c               print*,'debug aft hetchem',k,H2SO4_HETCHEM(I,J,L,k)
c     &                ,Nk(k),Mk(k,1),Mk(k,2)
c               enddo
c             endif

               mpnum=8 
               call aerodiag(mpnum,i,j,l)
#endif

! Do water eqm at appropriate times
               call ezwatereqm(Mk)

               call storenm()
               call mnfix(Nk,Mk)
               mpnum=6 
               call aerodiag(mpnum,i,j,l)
!in-cloud oxidation is not affected by time step.  
!! in-cloud oxidation! move from clouds2.f to here        

               call storenm()
                  call aqoxid(i,j,l,.TRUE.) ! Moist Convective clouds

               mpnum=4 
               call aerodiag(mpnum,i,j,l)

! Do water eqm at appropriate times
               call ezwatereqm(Mk)

               call storenm()
               call mnfix(Nk,Mk)
               mpnum=6 
               call aerodiag(mpnum,i,j,l)

               call storenm()
                  call aqoxid(i,j,l,.false.) ! Large scale clouds

               mpnum=5 
               call aerodiag(mpnum,i,j,l)

! Do water eqm at appropriate times
               call ezwatereqm(Mk)

               call storenm()
               call mnfix(Nk,Mk)
               mpnum=6 
               call aerodiag(mpnum,i,j,l)

                                ! do nitrogen and sulfur mass checks
                                ! get the total mass of N
               tot_n_2 = Gc(srtnh4)*14.d0/17.d0
               do k=1,ibins
                  tot_n_2 = tot_n_2 + Mk(k,srtnh4)*14.d0/18.d0
               enddo
               
                                ! get the total mass of S
               tot_s_2 = 0.0 !Gc(srtso4)
               do k=1,ibins
                  tot_s_2 = tot_s_2 + Mk(k,srtso4)*32.d0/96.d0
               enddo
               
               if(am_i_root())then

                  if (abs(tot_n_2-tot_n_1)/tot_n_1.gt.1.0D-4)then
                     print*,'Nitrogen not conserved in aerophys'
! 1                   print*,'i',i,'j',j,'l',l
!                     print*,'Init,Init1,Intm,Final',tot_n_i,tot_n_1,
!     *                    tot_n_1a,tot_n_2
                  endif
                  
                  if (abs(tot_s_2-tot_s_1)/tot_s_1.gt.1.0D-4)then !TOMAS - increase from 1.0D-4 
                     print*,'Sulfur not conserved in aerophys'
                     print*,'i',i,'j',j,'l',l,'Init,Final',tot_s_1,
     *                    tot_s_2
                  endif
               endif
              
C     Check for negative tracer problems
               flag=0
               do n=1,NBINS
                  if (Nk(n) .lt. 0.0) then
                     write(*,*) 'Nk < 0 for i,j,l,bin:',i,j,l,n
                     flag=1
                  endif
                  do jc=1,icomp
                     if (Mk(n,jc) .lt. 0.0) then
                    write(*,*) 'Mk < 0 for i,j,l,bin,comp:',i,j,l,n,jc
                        flag=1
                     endif
                  enddo
               enddo
               if (flag .eq. 1) 
     &              CALL STOP_MODEL('- tracer in TOMAS_DRV',255)   
                                ! regular bins
          call getCCN_kappa(i,j,l,Nk,Mk,Temp,boxmass,boxvol) 

!Save diagnostics! 
               do n=1,ibins       
!     Aerosol number             
                  tracnum=IDTNUMD-1+n 
                  tr3Dsource(i,j,l,nOther,tracnum)=
     &                 (NK(N)-INIT_NK(N))/dtsrc
                  
                  do np=1,ptype
                     if (ijts_TOMAS(np,tracnum).gt.0) 
     &                taijs(i,j,ijts_TOMAS(np,tracnum)) 
     &                    =taijs(i,j,ijts_TOMAS(np,tracnum))
     &                    +AEROD(i,j,l,tracnum,np) ! /adt
                     if (itcon_TOMAS(np,tracnum).gt.0) 
     &                    call inc_diagtcb(i,j,AEROD(i,j,l,tracnum,np),
     &                    itcon_TOMAS(np,tracnum),tracnum)
                  enddo

                  do jc=1,icomp-idiag
                     tracnum=IDTSO4-1+n+ibins*(jc-1)
                     tr3Dsource(i,j,l,nOther,tracnum)=
     &                    (MK(n,jc)-INIT_Mk(n,jc))/dtsrc

                  do np=1,ptype
                     if (ijts_TOMAS(np,tracnum).gt.0) 
     &                   taijs(i,j,ijts_TOMAS(np,tracnum)) 
     &                    =taijs(i,j,ijts_TOMAS(np,tracnum))
     &                    +AEROD(i,j,l,tracnum,np) ! /adt

                     if (itcon_TOMAS(np,tracnum).gt.0) 
     &                    call inc_diagtcb(i,j,AEROD(i,j,l,tracnum,np),
     &                    itcon_TOMAS(np,tracnum),tracnum)
                  enddo

                  enddo  
c$$$                  
                  tracnum=IDTH2O-1+n 
                  tr3Dsource(i,j,l,nOther,tracnum)=
     &                 (MK(N,SRTH2O)-INIT_MK(N,SRTH2O))/dtsrc
               enddo
               

!               tr3Dsource(i,j,l,nOther,n_H2SO4) =
!     *              (Gc(srtSO4)-INIT_H2SO4)/dtsrc
               TRM(I,J,L,n_H2SO4)=Gc(srtso4)              
                  do np=1,ptype
                     if (ijts_TOMAS(np,n_H2SO4).gt.0) 
     &                taijs(i,j,ijts_TOMAS(np,n_H2SO4)) 
     &                    =taijs(i,j,ijts_TOMAS(np,n_H2SO4))
     &                    +AEROD(i,j,l,n_H2SO4,np) ! /adt
                     if (itcon_TOMAS(np,n_H2SO4).gt.0) 
     &                    call inc_diagtcb(i,j,AEROD(i,j,l,n_H2SO4,np),
     &                    itcon_TOMAS(np,n_H2SO4),n_H2SO4)
                  enddo               


               tr3Dsource(i,j,l,nChemistry,n_NH3)=
     *              (Gc(srtNH4)-INIT_NH3)/dtsrc

                                ! aerosol ammonia
               tot_aam = 0.d0
               do n=1,NBINS
                  tot_aam = tot_aam + Mk(n,srtnh4)
               enddo
               
               tr3Dsource(i,j,l,nChemistry,n_NH4)=
     *              (tot_aam-INIT_NH4)/dtsrc

               tr3Dsource(i,j,l,nChemistry,n_SOAgas)=
     *              (SOAmass-INIT_SOA)/dtsrc !total SOA rate in 30 min

!               if(SOArate*dtsrc.GT.1)then
!              print*,'SOA_cond_AEROD',(SOAmass-INIT_SOA),
!     &                SUM(AEROD(i,j,l,IDTOCIL:IDTOCIL+NBINS-1,1)),
!     &                ijts_TOMAS(1,IDTOCIL), sum(taijs
!     &(i,j,ijts_TOMAS(1,IDTOCIL):ijts_TOMAS(1,IDTOCIL+NBINS-1)))
!               endif

               
               AEROD(i,j,l,:,:)=0.0 !for aeroupdate

C     End of loop over grid cells
            enddo               !I loop
c$$$  do K=1,2
c$$$  TAJLS(J,L,K+9) = TAJLS(J,L,K+9) + TSUM(K)
c$$$            enddo
         enddo                  !J loop
      enddo                     !L loop
      return
      END SUBROUTINE TOMAS_DRV                     !of main
      

      logical function is_nan(value)
      real*8 value
      if (abs(value).ge.0) then
         is_nan=.false.
      else
         is_nan=.true.
         endif
      return
      end function is_nan

      subroutine nanstop(value, line, var1, var2)
      real*8 value
      integer line, var1, var2

      if (abs(value).ge.0) then
      else
      write (*,*) 'line',line, var1, var2, value
      call flush(6)
      call stop_model('nan in nanstop',255)
      endif
      return
      end subroutine nanstop


!@sum  dep_getdp:  calculates the average diameter of aerosol
!@+     particles in a given GCM grid cell and size bin.
!@auth  Yunha Lee
!@ver   1.0

      subroutine dep_getdp(i,j,l,getdp,size_density)                                            
      USE TRACER_COM, only : nbins,IDTSO4,IDTNA,IDTECIL,
     &     IDTECOB,IDTOCIL,IDTOCOB,IDTDUST,IDTH2O,
     &     IDTNUMD,ntm,xk,trm
      USE CONSTANT,   only : pi,lhe,mair,gasc    
      USE MODEL_COM, only : t            ! potential temperature (C)
     $                     ,q            ! saturated pressure
      USE DYNAMICS,   only: pmid,pk    ! midpoint pressure in hPa (mb)
!                                           and pk is t mess up factor
      USE TOMAS_AEROSOL

      IMPLICIT NONE

C-----ARGUMENT DECLARATIONS------------------------------------------

      integer i,j,l  !coordinate of GCM grid cell
      integer n      !tracer index

C-----VARIABLE DECLARATIONS------------------------------------------

      integer k     !size bin index
      real density                 !density (kg/m3) of current size bin
      real mso4, mh2o, mno3, mnh4  !mass of each component (kg/grid box)
      real mecil,mecob,mocil,mocob
      real mdust,mtot,mnacl             
      real*8 mp          !particle mass (kg)
      real*8 mu          !air viscosity (kg/m s)
      real*8 qsat
      real aerodens
      real,intent(out),DIMENSION(nbins) :: getdp,size_density
C-----VARIABLE COMMENTS----------------------------------------------

C-----ADJUSTABLE PARAMETERS------------------------------------------
      real*8 Neps  !a small number of particles (#/box)
      parameter (Neps=1.d-20)

C-----CODE-----------------------------------------------------------

!     Compute particle diameter for each bins - YUNHA LEE 
      do k=1,nbins

C     Swap GCM variables into aerosol algorithm variables
        Nk(k)=trm(i,j,l,IDTNUMD-1+k)
        Mk(k,srtso4)=trm(i,j,l,IDTSO4-1+k)
        Mk(k,srtna )=trm(i,j,l,IDTNA -1+k)
        Mk(k,srtnh4)=0.1875*Mk(k,srtso4) ! artificial for now.. 0.0!t0m(i,j,l,IDTNH4-1+n)
        MK(k,srtecob)=trm(i,j,l,IDTECOB -1+k)
        MK(k,srtecil)=trm(i,j,l,IDTECIL -1+k)
        MK(k,srtocob)=trm(i,j,l,IDTOCOB -1+k)
        MK(k,srtocil)=trm(i,j,l,IDTOCIL -1+k) 
        MK(k,srtdust)=trm(i,j,l,IDTDUST -1+k) 
        Mk(k,srth2o)= trm(i,j,l,IDTH2O-1+k) !I don't think this is necessary!
      enddo

      temp = pk(l,i,j)*t(i,j,l) !should be in [K]
      rh = MIN(1.,q(i,j,l)/QSAT(temp,lhe,pmid(l,i,j))) ! rH [0-1]
            
      call mnfix(Nk,Mk)  
      call ezwatereqm(Mk)

      do k=1,nbins 
         if (Nk(k) .eq. 0.0) then
            if (Mk(k,srtso4) .gt. 1.) then
               print*, 'ERROR in getdp - # = but mass > 0',i,j
               print*, 'bin=',k
               print*, 'TRM(#)=',Nk(k)
               print*, 'TRM(SO4)=',Mk(k,srtso4)
               print*, 'TRM(NACL)=',Mk(k,srtna)
               print*, 'TRM(OCIL)=',Mk(k,srtocil)
               call stop_model('ERROR IN getdp',255)
            endif
         endif

         mso4=Mk(k,srtso4) 
         mnacl=Mk(k,srtna)
         mno3=0.e0
         if ((mso4+mno3) .lt. 1.e-8) mso4=1.e-8
         mnh4=0.1875*Mk(k,srtso4)  !assume ammonium bisulfate
         mecob=Mk(k,srtecob)
         mecil=Mk(k,srtecil)
         mocil=Mk(k,srtocil)
         mocob=Mk(k,srtocob)
         mdust=Mk(k,srtdust)          
         mh2o=Mk(k,srth2o)   

!in CLOUDS2.f - some tracers goes negative..
!So, set to zero to prevent a problem in density calculation
         if(mnacl.lt.0) mnacl=0.
         if(mecob.lt.0) mecob=0.
         if(mecil.lt.0) mecil=0.
         if(mocob.lt.0) mocob=0.
         if(mocil.lt.0) mocil=0.
         if(mdust.lt.0) mdust=0.
         if(mh2o.lt.0) mh2o=0.

         density=aerodens(mso4,mno3,mnh4 !mno3 taken off!
     *        ,mnacl,mecil,mecob,mocil,mocob,mdust,mh2o) !assume bisulfate   

         mtot= 1.1875*Mk(k,srtso4)+mnacl+mecil+mecob+
     *        mocil+mocob+mdust+mh2o
   
         size_density(k)=density      

         if (Nk(k) .gt. Neps.and.mtot.gt.0.) then
            mp=mtot/Nk(k)
         else
            mp=sqrt(xk(k+1)*xk(k))
            if(Nk(k) .gt. Neps) 
     &           print*,'Warning in getdp:#>Neps but mtot=0',
     &           k,mtot,Nk(k)
         endif
         
!     fix unrealistically large mp for low aerosol conc.
         if (mp .gt. 1.d3*xk(NBINS+1)) then
            
            if ((Nk(k) .lt. 1.d5) .and. !negligible amount of aerosol - fudge mp
     &           (Mk(k,srtso4) .lt. 3.)) then
               mp=sqrt(xk(k+1)*xk(k))
            else
               if (Nk(k) .gt. 1.d12) then
!MODELE-TOMAS: during CONDSE, TM(H2O) is so large that causes too big mp. 
!MODELE-TOMAS: So, if dry mass is less than the max size boundary, just take the max mp. 
                  if((mtot-Mk(k,srth2o)).lt.
     &                 1.d1*xk(nbins+1)*Nk(k))then
                     print*,'Fudge mp in getdp: large mp by AH2O'
                     mp=sqrt(xk(nbins+1)*xk(nbins))                     
                  else
                  print*,'ERROR in getdp: mp too large'
                  print*, 'bin=',k
                  print*, 'TM(#)=', Nk(k)
                  print*, 'TM(SO4)=', mso4, mh2o
                  print*, 'TM(NACL)=', mnacl, mdust
                  print*, 'TM(OC)=',mocob,mocil
                  print*, 'TM(EC)=',mecob,mecil
                  call stop_model('mp too large getdp',255)
                  endif
               endif
            endif
         endif
         getdp(k)=(6.d0*mp/(pi*size_density(k)))**(1.d0/3.d0)
      enddo

      RETURN
      END subroutine dep_getdp

!@sum  readfraction: read fraction lookup table for wet deposition
!@auth Peter Adams
!@ver  1.0
      subroutine readfraction(infile,fraction2)
      
      implicit none
      
#ifdef TOMAS_12_10NM
	character*17 infile
#endif
#ifdef TOMAS_12_3NM
	character*21 infile
#endif

      integer innum, ii, jj, kk
      
      real*8,intent(out),dimension(101,101,101):: fraction2
      parameter (innum=580)
 1    format(f6.5)
      open(unit=innum,file=infile,FORM='FORMATTED',STATUS='OLD')
      do ii=1,101
        do jj=1,101
          do kk=1,101
            read(innum,1) fraction2(kk,jj,ii)
            if (fraction2(kk,jj,ii).gt.1.) fraction2(kk,jj,ii)=0.
          enddo
        enddo
      enddo
!     print*,'fraction last',fraction(101,101,101)
      close(innum)
      return 
      end subroutine readfraction
      
!     @sum  readbinact: read binact lookup table for wet deposition
!@auth Peter Adams
!@ver  1.0
      subroutine readbinact(infile,binact)
      
      USE TOMAS_AEROSOL, ONLY : ibins
      
      implicit none
      
#ifdef TOMAS_12_10NM
	character*15 infile
#endif
#ifdef TOMAS_12_3NM
	character*19 infile
#endif
      integer innum, ii, jj, kk
      integer,intent(out),dimension(101,101,101):: binact
      parameter (innum=590)
 1    format(I2)
      open(unit=innum,file=infile,FORM='FORMATTED',STATUS='OLD')
      do ii=1,101
        do jj=1,101
          do kk=1,101
            read(innum,1) binact(kk,jj,ii)
            if (binact(kk,jj,ii).eq.0) binact(kk,jj,ii)=ibins+1
          enddo
        enddo
      enddo
      close(innum)
      return
      end subroutine readbinact
      
!     @sum  getfraction: compute composition fraction for wet deposition
!@auth Peter Adams/Yunha Lee 
!@ver  1.0
      subroutine getfraction(tr_conv,tm,fract)
      USE TOMAS_AEROSOL, ONLY : binact02,binact10,
     &     fraction02,fraction10 
      USE TRACER_COM, only : nbins,ntm,IDTECIL,
     &     IDTOCIL,IDTOCOB,IDTSO4,IDTNA,IDTDUST,
     &     IDTECOB
      
      implicit none
      
      real mecil, mocil, mocob, mso4, mnacl,mdust, mtot
      real xocil, xso4, xnacl
      integer iso4, inacl, iocil,k
      integer getbinact
      real*8,dimension(nbins) ::  fract
      REAL*8, DIMENSION(ntm) :: TM
      LOGICAL TR_CONV
      
      do k=1, nbins
        mecil=TM(IDTECIL-1+k)
	mocil=TM(IDTOCIL-1+k)
	mocob=TM(IDTOCOB-1+k)
	mso4=TM(IDTSO4-1+k)*1.2 !account for ammonium sulfate
	mnacl=TM(IDTNA-1+k)
	mdust=TM(IDTDUST-1+k)
        
CCC   (IMPORTANT) This is temporal treatment for negative tracers. 
CCC   In AR5_v2_branch,sometimes tracer mass become negative. 
CCC   Until solution is available, turn off the stop_model and set zero for fraction.
        if(mso4.lt.0) mso4=0.
        if(mnacl.lt.0) mnacl=0.
        if(mocil.lt.0) mocil=0.
        if(mocob.lt.0) mocob=0.
        if(mecil.lt.0) mecil=0.
        if(mdust.lt.0) mdust=0.
        
CCC   (IMPORTANT) END here 
        
	mtot=mecil+mocil+mocob+mso4+mnacl+mdust+1.e-20
	xocil=mocil/mtot
	xso4=mso4/mtot
	xnacl=mnacl/mtot
	iso4=min(101,int(xso4*100)+1)
	inacl=min(101,int(xnacl*100)+1)
	iocil=min(101,int(xocil*100)+1)

        if(xso4.lt.0.or.xnacl.lt.0.or. xocil.lt.0)then
           print*,'wrong getfraction'
           print*,'mass',mso4,mnacl,mecil,mocob
     &          ,mocil,mdust,tr_conv
             call stop_model('wrong getfraction',255)
        endif

        if (tr_conv)then        !convective clouds
          getbinact=binact10(iso4,inacl,iocil)
          
          if(binact10(iso4,inacl,iocil).lt.0)then
            print*,'wrong binact',binact10(iso4,inacl,iocil)
            print*,'iso4',iso4,inacl,iocil
            print*,'mass',mso4,mnacl,mocil,mtot
            call stop_model('wrong binact',255)
          endif
          
          if (getbinact.gt.k) then
            fract(k)=0.         !not activated
          else if (getbinact.eq.k) then
            
            fract(k)=fraction10(iso4,inacl,iocil) !partly activated
          else
            fract(k)=1.         !all sizebin activated
          endif
          if(getbinact.le.2)
     &         print*,'CONV CLD',getbinact,k,iso4,inacl,iocil
        else                    !large-scale
          getbinact=binact02(iso4,inacl,iocil)
          
          if(binact02(iso4,inacl,iocil).lt.0)then
            print*,'wrong binact',binact02(iso4,inacl,iocil)
            print*,'iso4',iso4,inacl,iocil
            print*,'mass',mso4,mnacl,mocil,mtot
          endif
          
          if (getbinact.gt.k) then
            fract(k)=0.         !not activated
          else if (getbinact.eq.k) then
            
            fract(k)=fraction02(iso4,inacl,iocil) !partly activated
          else
            fract(k)=1.         !all sizebin activated
            
            if(getbinact.le.4)
     &           print*,'STRAT CLD',getbinact,k,iso4,inacl,iocil
            
          endif
        endif
      enddo                     !k
      return
      end subroutine getfraction
      
!@sum  aqoxid: takes an amount of SO4 produced via in-cloud
!@+   oxidation and condenses it onto an existing aerosol size
!@+   distribution.  It assumes that only particles larger than the
!@+   critical activation diameter activate and that all of these have
!@+   grown to roughly the same size.  Therefore, the mass of SO4 
!@+   produced by oxidation is partitioned to the various size bins
!@+   according to the number of particles in that size bin.  
!@+   Values of tau are calculated for each size bin accordingly and
!@+   the cond subroutine is called to update Nk and Mk.
!@auth Peter Adams and Yunha Lee 
!@ver  1.0

C-----INPUTS------------------------------------------------------------

C-----OUTPUTS-----------------------------------------------------------

      SUBROUTINE aqoxid(i,j,l,tr_conv)

C-----INCLUDE FILES-----------------------------------------------------

      USE TOMAS_AEROSOL
      USE TRACER_COM, only : ntm, IDTECIL,
     &       IDTOCIL,IDTOCOB,IDTSO4,IDTNA,IDTDUST,
     &       IDTECOB,IDTH2O,xk,nbins

      IMPLICIT NONE
C-----ARGUMENT DECLARATIONS---------------------------------------------

      real*8 moxid !mass of new sulfate from in-cloud oxid.
      real*8, dimension(ibins) :: fraction    !fraction activated for every sizebin

c      real*8 tracc(NTM)

C-----VARIABLE DECLARATIONS---------------------------------------------

      real*8 Nact, Mact  !#/mass of activated particles
      real*8 mpo   !initial particle mass (kg)
      real*8 mpw   !initial particle wet mass (kg)
      real*8 aqtau(ibins)
      integer k,mpnum,n,tracnum,i,j,l
      real*8 Nko(ibins), Mko(ibins, icomp) !input to cond routine
      real*8 Nkf(ibins), Mkf(ibins, icomp) !output from cond routine
      real*8 tdt      !the value 2/3
      real*8,parameter :: eps=1.d-40
      integer jc
      real*8 frac      
      real*8 WR                ! wet ratio = total mass/ dry mass (win, 5/15/06)
      real*8 mox(ibins) !mass of new sulfate per particle in each bin
      real*8 tot_aam ! total aerosol ammonia
      real*8 TM(ntm) ! total aerosol ammonia

      LOGICAL TR_CONV

C-----CODE--------------------------------------------------------------

              
      if (tr_conv) then
         moxid=AQSO4oxid_mc(i,j,l) 
      else
         moxid=AQSO4oxid_ls(i,j,l) 
      endif

      if (moxid.eq.0.d0) return

      tdt=2.d0/3.d0

      TM(:)=0.0
!only mass needed for getfraction
      do n=1,IBINS
         do jc=1,icomp-idiag
            tracnum=IDTSO4-1+n+ibins*(jc-1)
            TM(tracnum)=Mk(n,jc)               
         enddo
         tracnum=IDTH2O-1+n
         TM(tracnum)=Mk(n,srth2o)
      enddo

      if (tr_conv) then
         CALL getfraction (.true.,TM,FRACTION) !1% supersaturation assumption
      else
         CALL getfraction(.false.,TM,FRACTION) !0.2% supersaturation assumption
      endif

      Nact=0.0
      Mact=0.0
      do k=1,ibins
         Nact=Nact+Nk(k)*fraction(k)
         do jc=1,icomp-idiag
            Mact=Mact+Mk(k,jc)*fraction(k)
         enddo
      enddo

      if ((Mact+moxid)/(Nact+eps) .gt. xk(ibins)) then !YHL- I change xk(ibins-1) to xk(ibins)
!            if (TAU .gt. 8350.) then
cdebug               write(*,*) 'ERROR in aqoxidcc: Ave size out of bounds'
c$$$               write(*,*) 'Nact: ',Nact
c$$$               write(*,*) 'moxid/Mact: ',moxid,Mact
c$$$               do k=1,ibins
c$$$                  write(*,*) 'k, N, MSO4, MH2O: ',k,Nk(k),
c$$$     &                  Mk(k,srtso4),Mk(k,srth2o)
c$$$               enddo
               goto 20
c$$$            else
c$$$               !don't worry about the first two weeks
c$$$               goto 20
c$$$            endif
       endif

C Calculate tau for each size bin
      moxid=moxid/(Nact+eps) !now kg H2SO4 per activated particle
      do k=1,ibins
         mox(k)=fraction(k)*moxid
         if (fraction(k) .eq. 0.) then
            !too small to activate - no sulfate for this particle
            aqtau(k)=0.0
         else
            !activated particle - calculate appropriate tau
            mpo=0.0
            mpw=0.0
            !WIN'S CODE MODIFICATION 6/19/06
            !THIS MUST CHANGED WITH THE NEW dmdt_int.f
            do jc=1,icomp-idiag
               mpo = mpo+Mk(k,jc)  !accumulate dry mass
            enddo
            do jc=1,icomp
               mpw = mpw+Mk(k,jc)  ! have wet mass include amso4
            enddo
            WR = mpw/mpo  !WR = wet ratio = total mass/dry mass
            if (Nk(k) .gt. 0.d0) then
               mpw=mpw/Nk(k)
               aqtau(k)=1.5d0*((mpw+mox(k)*WR)**tdt-mpw**tdt)  !added WR to moxid term (win, 5/15/06)
            else
               !nothing in this bin - set tau to zero
               aqtau(k)=0.0
               mox(k)=0.d0
            endif
         endif
      enddo

      Nko(:)=Nk(:)
      Mko(:,:)=Mk(:,:)

      call tmcond(aqtau,xk,Mko,Nko,Mkf,Nkf,srtso4,mox)

      Nk(:)=Nkf(:)
      Mk(:,:)=Mkf(:,:)   


 20      continue   !go here if process is skipped  

      RETURN
      END SUBROUTINE aqoxid

#ifdef TOMAS_HETCHEM

!@sum cond_hetchem: takes an amount of SO4 produced via HETCHEM
!@+   on dust particles and update number/size distribution. 
!@+   The amount of SO4 formed each size bin is determined from HETCHEM. 
!@+   
!@+   CURRENT (July 2014)- it works for H2SO4 condensation (srtso4 hardcoded)
!@+   To use Nitrate condensation, it must be modified!!! 
!@+ 
!@auth Yunha Lee 
!@ver  1.0

C-----INPUTS------------------------------------------------------------

C-----OUTPUTS-----------------------------------------------------------

      SUBROUTINE cond_hetchem(i,j,l)

C-----INCLUDE FILES-----------------------------------------------------

      USE TOMAS_AEROSOL
      USE TRACER_COM, only : ntm, IDTECIL,
     &       IDTOCIL,IDTOCOB,IDTSO4,IDTNA,IDTDUST,
     &       IDTECOB,IDTH2O,xk,nbins
      USE MODEL_COM, only : dtsrc

      IMPLICIT NONE

C-----VARIABLE DECLARATIONS---------------------------------------------

      real*8 mpo   !initial particle mass (kg)
      real*8 mpw   !initial particle wet mass (kg)
      real*8 hctau(ibins)
      integer k,mpnum,n,tracnum,i,j,l
      real*8 Nko(ibins), Mko(ibins, icomp) !input to cond routine
      real*8 Nkf(ibins), Mkf(ibins, icomp) !output from cond routine
      real*8 tdt      !the value 2/3
      real*8,parameter :: eps=1.d-40
      integer jc     
      real*8 WR                ! wet ratio = total mass/ dry mass (win, 5/15/06)
      real*8 mox(ibins) !mass of new sulfate per particle in each bin


      LOGICAL TR_CONV

C-----CODE--------------------------------------------------------------

      mox(:)=H2SO4_HETCHEM(i,j,l,:)*dtsrc 
C if hetchem doesn't occur, return. 
      if (sum(mox(1:ibins)).le.0.d0) return

      tdt=2.d0/3.d0

C Calculate tau for each size bin
      do k=1,ibins

c        mox(k)=mox(k)/Nk(k)

            mpo=0.0
            mpw=0.0
            !WIN'S CODE MODIFICATION 6/19/06
            !THIS MUST CHANGED WITH THE NEW dmdt_int.f
            do jc=1,icomp-idiag
               mpo = mpo+Mk(k,jc)  !accumulate dry mass
            enddo
            do jc=1,icomp
               mpw = mpw+Mk(k,jc)  ! have wet mass include amso4
            enddo
            WR = mpw/mpo  !WR = wet ratio = total mass/dry mass
            if (Nk(k) .gt. 0.d0) then
               mox(k)=mox(k)/Nk(k)
               mpw=mpw/Nk(k)
               hctau(k)=1.5d0*((mpw+mox(k)*WR)**tdt-mpw**tdt)  !added WR to moxid term (win, 5/15/06)
            else
               !nothing in this bin - set tau to zero
               hctau(k)=0.0
               mox(k)=0.d0
            endif
      enddo

      Nko(:)=Nk(:)
      Mko(:,:)=Mk(:,:)

      call tmcond(hctau,xk,Mko,Nko,Mkf,Nkf,srtso4,mox)

      Nk(:)=Nkf(:)
      Mk(:,:)=Mkf(:,:)   

c      if(sum(h2so4_hetchem(i,j,l,1:ibins)).gt.0.and.i.eq.30)
c     * print*,'COND_HET',sum(mox(1:ibins)),sum(Mko(1:ibins,1)),
c     *     sum(Mkf(1:ibins,1)),sum(Mko(1:ibins,2)),sum(Mkf(1:ibins,2))

      RETURN
      END SUBROUTINE cond_hetchem

#endif

!@sum  aerodens : this function calculates the aerosol density (kg/m3)  
!@auth Peter Adams
Ckpc  Jan.,2002 - extended to include carbonaceous aerosols
!@ver  1.0

      real FUNCTION aerodens(mso4,mno3,mnh4,mnacl,mecil,
     & mecob,mocil,mocob,mdust,mh2o)

      USE TRACER_COM, only : trpdens,IDTECOB,IDTOCOB,IDTDUST,N_AECOB
      IMPLICIT NONE
!@var mso4, mno3, mnh4, mh2o, mnacl - These are the masses of each aerosol
!@+   component.  Since the density is an intensive property,
!@+   these may be input in a variety of units (ug/m3, mass/cell, etc.).

      real mso4, mno3, mnh4, mnacl, mecil,mecob,mocil,mocob,mdust,mh2o
      real inodens, idensity !, dec,doc,ddust
!     external inodens
!!      parameter(dec=2200., doc=1400., ddust=2650.)

      idensity=inodens(mso4, mno3,mnh4, mnacl, mh2o)

      aerodens=(idensity*(mso4+mno3+mnh4+mnacl+mh2o) !mno3 taken out! 
!!     &  +dec*(mecil+mecob)+doc*(mocil+mocob)
!!     &     +ddust*mdust)
     &  +trpdens(idtecob)*(mecil+mecob)+trpdens(idtocob)*(mocil+mocob)
     &     +trpdens(idtdust)*mdust)
     &  /(mso4+mno3+mnh4+mnacl+mh2o+mecil+mecob+mocil+mdust+mocob)

      RETURN
      END FUNCTION aerodens


!@sum  inodens :this function calculates the density (kg/m3) of a sulfate-
!@+   nitrate-ammonium-nacl-water mixture that is assumed to be internally
!@+   mixed.  

!@auth Peter Adams, May 1999
!@+   November, 2001 - extended to include NaCl and bug fixed
!@+                    the bug was that species densities (dan, ds0,
!@+                    etc...) are supposed to be calculated based on
!@+                    *total* solute concentration, not each species
!@+                    contribution as it had been.
C-----Literature cited--------------------------------------------------
!@+   I. N. Tang and H. R. Munkelwitz, Water activities, densities, and
!@+     refractive indices of aqueous sulfates and sodium nitrate droplets
!@+     of atmospheric importance, JGR, 99, 18,801-18,808, 1994
!@+   Ignatius N. Tang, Chemical and size effects of hygroscopic aerosols
!@+     on light scattering coefficients, JGR, 101, 19,245-19,250, 1996
!@+   Ignatius N. Tang, Thermodynamic and optical properties of mixed-salt
!@+     aerosols of atmospheric importance, JGR, 102, 1883-1893, 1997

!@var   mso4, mno3, mnh4, mh2o, mnacl - These are the masses of each aerosol
!@+   component.  Since the density is an intensive property,
!@+   these may be input in a variety of units (ug/m3, mass/cell, etc.).

      real FUNCTION inodens(mso4, mno3, mnh4, mnacl, mh2o)

      IMPLICIT NONE
      real mso4, mno3, mnh4, mnacl, mh2o
      real so4temp, no3temp, nh4temp, nacltemp, h2otemp  !store initial values
      real mwso4, mwno3, mwnh4, mwnacl, mwh2o            !molecular weights
      real ntot, mtot, drytot                      !total number of moles, mass
      real nso4, nno3, nnh4, nnacl, nh2o       !moles of each species
      real xso4, xno3, xnh4, xnacl, xh2o       !mole fractions
      real rso4, rno3, rnh4, rnacl, rh2o       !partial molar refractions
      real ran, rs0, rs1, rs15, rs2       !same, but for solute species
      real asr                            !ammonium/sulfate molar ratio
      real nan, ns0, ns1, ns15, ns2, nss  !moles of dry solutes (nss = sea salt)
      real xan, xs0, xs1, xs15, xs2, xss  !mass % of dry solutes - Tang (1997) eq. 10
      real dan, ds0, ds1, ds15, ds2, dss  !binary solution densities - Tang (1997) eq. 10
      real mwan, mws0, mws1, mws15, mws2  !molecular weights
      real yan, ys0, ys1, ys15, ys2, yss  !mole fractions of dry solutes
      real yh2o
      real d                              !mixture density
      real xtot

C     In the lines above, "an" refers to ammonium nitrate, "s0" to 
C     sulfuric acid, "s1" to ammonium bisulfate, and "s2" to ammonium sulfate.
C     "nacl" or "ss" is sea salt.

      parameter(mwso4=96., mwno3=62., mwnh4=18., mwh2o=18., 
     &          mwnacl=58.45)
      parameter(mwan=mwnh4+mwno3, mws0=mwso4+2., mws1=mwso4+1.+mwnh4,
     &          mws2=2*mwnh4+mwso4)

C Save initial component masses to restore later 
!      mno3=0.d0
      so4temp=mso4
      no3temp=mno3
      nh4temp=mnh4
      h2otemp=mh2o
      nacltemp=mnacl
C Calculate mole fractions
      mtot = mso4+mno3+mnh4+mnacl+mh2o
      drytot = mso4+mno3+mnh4+mnacl
      if (drytot .lt. 1.e-15) then
      inodens=1000.
      return
      endif
      nso4 = mso4/mwso4
      nno3 = mno3/mwno3  !nno3 =zero 
      nnh4 = mnh4/mwnh4
      nnacl = mnacl/mwnacl
      nh2o = mh2o/mwh2o
      ntot = nso4+nno3+nnh4+nnacl+nh2o
      xso4 = nso4/ntot
      xno3 = nno3/ntot
      xnh4 = nnh4/ntot
      xnacl = nnacl/ntot
      xh2o = nh2o/ntot

C If there are more moles of nitrate than ammonium, treat unneutralized
C HNO3 as H2SO4
      if (nno3 .gt. nnh4) then  !will never occur as no3 is always zero
         !make the switch
         nso4=nso4+(nno3-nnh4)
         nno3=nnh4
         mso4=nso4*mwso4
         mno3=nno3*mwno3

         !recalculate quantities
         mtot = mso4+mno3+mnh4+mnacl+mh2o
         nso4 = mso4/mwso4
         nno3 = mno3/mwno3
         nnh4 = mnh4/mwnh4
         nnacl = mnacl/mwnacl
         nh2o = mh2o/mwh2o
         ntot = nso4+nno3+nnh4+nnacl+nh2o
         xso4 = nso4/ntot
         xno3 = nno3/ntot
         xnh4 = nnh4/ntot
         xnacl = nnacl/ntot
         xh2o = nh2o/ntot

      endif

C Calculate the mixture density
C Assume that nitrate exists as ammonium nitrate and that other ammonium
C contributes to neutralizing sulfate
      nan=nno3
      if (nnh4 .gt. nno3) then 
         !extra ammonium
         asr=(nnh4-nno3)/nso4 
      else
         !less ammonium than nitrate - all sulfate is sulfuric acid
         asr=0.0                !if nnh4=0, then asr=0 
      endif
      if (asr .ge. 2.) asr=2.0
      if (asr .ge. 1.) then
         !assume NH4HSO4 and (NH4)2(SO4) mixture
         !NH4HSO4
         ns1=nso4*(2.-asr)
         !(NH4)2SO4
         ns2=nso4*(asr-1.)
         ns0=0.0
      else
         !assume H2SO4 and NH4HSO4 mixture
         !NH4HSO4
         ns1=nso4*asr
         !H2SO4
         ns0=nso4*(1.-asr)
         ns2=0.0
      endif

      !Calculate weight percent of solutes
      xan=nan*mwan/mtot*100.
      xs0=ns0*mws0/mtot*100.
      xs1=ns1*mws1/mtot*100.
      xs2=ns2*mws2/mtot*100.
      xnacl=nnacl*mwnacl/mtot*100.
      xtot=xan+xs0+xs1+xs2+xnacl
!      call nanstop(xtot,153,0,0)
      !Calculate binary mixture densities (Tang, eqn 9)
      dan=0.9971 +4.05e-3*xtot +9.0e-6*xtot**2.
      ds0=0.9971 +7.367e-3*xtot -4.934e-5*xtot**2. +1.754e-6*xtot**3.
     &       -1.104e-8*xtot**4.
      ds1=0.9971 +5.87e-3*xtot -1.89e-6*xtot**2. +1.763e-7*xtot**3.
      ds2=0.9971 +5.92e-3*xtot -5.036e-6*xtot**2. +1.024e-8*xtot**3.
      dss=0.9971 +7.41e-3*xtot -3.741e-5*xtot**2. +2.252e-6*xtot**3.
     &       -2.06e-8*xtot**4.

      !Convert x's (weight percent of solutes) to fraction of dry solute (scale to 1)
      xtot=xan+xs0+xs1+xs2+xnacl
      xan=xan/xtot
      xs0=xs0/xtot
      xs1=xs1/xtot
      xs2=xs2/xtot
      xnacl=xnacl/xtot

      !Calculate mixture density
      d=1./(xan/dan+xs0/ds0+xs1/ds1+xs2/ds2+xnacl/dss)  !Tang, eq. 10
c      call nanstop(d,173,0,0)
      if (abs(d).ge.0) then
      else
         write(*,*) d,xtot,xan,xs0,xs1,xs2,xnacl,dan,ds0,ds1,ds2,dss
         write(*,*) 'woo',asr,mtot,mso4,mno3,mnh4,mnacl,mh2o
         write(*,*) 'woo',so4temp,no3temp,nh4temp,nacltemp,h2otemp
         call stop_model('nan in inodens',255)
      endif
      if ((d .gt. 2.) .or. (d .lt. 0.997)) then
         write(*,*) 'ERROR in inodens'
         write(*,*) mso4,mno3,mnh4,mnacl,mh2o
         call stop_model('ERROR in inodens',255)
      endif

C Restore masses passed
      mso4=so4temp
      mno3=no3temp
      mnh4=nh4temp
      mnacl=nacltemp
      mh2o=h2otemp
C Return the density
      inodens=1000.*d    !Convert g/cm3 to kg/m3

      RETURN
      END FUNCTION inodens


!@sum ezwatereqm : uses the current RH to calculate how much water is 
!@+   in equilibrium with the aerosol.  Aerosol water concentrations 
!@+   are assumed to be in equilibrium at all times and the array of 
!@+   concentrations is updated accordingly.

!@+   This version of the routine works for sulfate, sea salt, organic
!@+   particles.  They are assumed to be externally mixed and their
!@+   associated water is added up to get total aerosol water.
!@+   wr is the ratio of wet mass to dry mass of a particle.  Instead
!@+   of calling a thermodynamic equilibrium code, this routine uses a
!@+   simple curve fits to estimate wr based on the current humidity.
!@+   The curve fit is based on ISORROPIA results for ammonium bisulfate
!@+   at 273 K and sea salt at 273 K.
!@auth Peter Adams, March 2000

      SUBROUTINE ezwatereqm(Mke)

      USE TOMAS_AEROSOL
      IMPLICIT NONE

      real*8 Mke(ibins,icomp)
      integer k
      real*8 so4mass, naclmass, ocilmass
      real*8 wrso4, wrnacl, wrocil
      real*8 rhe
      real*8 waterso4, waternacl, waterocil

      rhe=100.d0*rh
      if (rhe .gt. 99.d0) rhe=99.d0
      if (rhe .lt. 1.d0) rhe=1.d0

      do k=1,ibins

         so4mass=Mke(k,srtso4)*1.1875  !converts kg so4 to kg nh4hso4
         naclmass=Mke(k,srtna)      !already as kg nacl - no conv necessary
         ocilmass=MKe(k,srtocil)    !already as kg ocil

         wrso4=waterso4(rhe)
         wrnacl=waternacl(rhe)
         wrocil=waterocil(rhe)

         Mke(k,srth2o)=so4mass*(wrso4-1.d0)+naclmass*(wrnacl-1.d0)
     &                 +ocilmass*(wrocil-1.d0)

      enddo

      RETURN
      END SUBROUTINE ezwatereqm


!@sum eznh3eqm :  puts ammonia to the particle phase until 
!@+   there is 2 moles of ammonium per mole of sulfate and the remainder
!@+   of ammonia is left in the gas phase.

!@auth Jeff Pierce, April 2007

      SUBROUTINE eznh3eqm(Gce,Mke)

      USE TOMAS_AEROSOL
      IMPLICIT NONE
      real*8 Gce(icomp)
      real*8 Mke(ibins,icomp)
      integer k
      real*8 tot_nh3  !total kmoles of ammonia
      real*8 tot_so4  !total kmoles of so4
      real*8 sfrac    !fraction of sulfate that is in that bin


      ! get the total number of kmol nh3
      tot_nh3 = Gce(srtnh4)/17.d0
      do k=1,ibins
         tot_nh3 = tot_nh3 + Mke(k,srtnh4)/18.d0
      enddo

      ! get the total number of kmol so4
      tot_so4 = 0.d0
      do k=1,ibins
         tot_so4 = tot_so4 + Mke(k,srtso4)/96.d0
      enddo

      ! see if there is free ammonia
      if (tot_nh3/2.d0.lt.tot_so4)then  ! no free ammonia
         Gce(srtnh4) = 0.d0 ! no gas phase ammonia
         do k=1,ibins
            sfrac = Mke(k,srtso4)/96.d0/tot_so4
            Mke(k,srtnh4) = sfrac*tot_nh3*18.d0 ! put the ammonia where the sulfate is
         enddo
      else ! free ammonia
         do k=1,ibins
            Mke(k,srtnh4) = Mke(k,srtso4)/96.d0*2.d0*18.d0 ! fill the particle phase
         enddo
         Gce(srtnh4) = (tot_nh3 - tot_so4*2.d0)*17.d0 ! put whats left over in the gas phase
      endif

      RETURN
      END 


!@sum waterso4 : uses the current RH to calculate how much water is 
!@+   in equilibrium with the sulfate.  Aerosol water concentrations 
!@+   are assumed to be in equilibrium at all times and the array of 
!@+   concentrations is updated accordingly.

!@+   waterso4 is the ratio of wet mass to dry mass of a particle.  Instead
!@+   of calling a thermodynamic equilibrium code, this routine uses a
!@+   simple curve fit to estimate wr based on the current humidity.
!@+   The curve fit is based on ISORROPIA results for ammonium bisulfate
!@+   at 273 K.

!@auth Peter Adams, March 2000

      real*8 FUNCTION waterso4(rhe)
      IMPLICIT NONE

      real*8 rhe   !relative humidity (0-100 scale)

!@var   waterso4 is the ratio of wet mass to dry mass of a particle.

C-----CODE--------------------------------------------------------------

      if (rhe .gt. 99.) rhe=99.
      if (rhe .lt. 1.) rhe=1.

         if (rhe .gt. 96.) then
            waterso4=
     &      0.7540688*rhe**3-218.5647*rhe**2+21118.19*rhe-6.801999e5
         else
         if (rhe .gt. 91.) then
            waterso4=8.517e-2*rhe**2 -15.388*rhe +698.25
         else
         if (rhe .gt. 81.) then
            waterso4=8.2696e-3*rhe**2 -1.3076*rhe +53.697
         else
         if (rhe .gt. 61.) then
            waterso4=9.3562e-4*rhe**2 -0.10427*rhe +4.3155
         else
         if (rhe .gt. 41.) then
            waterso4=1.9149e-4*rhe**2 -8.8619e-3*rhe +1.2535
         else
            waterso4=5.1337e-5*rhe**2 +2.6266e-3*rhe +1.0149
         endif
         endif
         endif
         endif
         endif

         !check for error
         if (waterso4 .gt. 30.) then
            write(*,*) 'ERROR in waterso4'
            write(*,*) rhe,waterso4
            call stop_model('ERROR in waterso4',255)
         endif

      RETURN
      END  FUNCTION waterso4




!@sum waternacl: This function uses the current RH to calculate how much water is 
!@+   in equilibrium with the seasalt.  Aerosol water concentrations 
!@+   are assumed to be in equilibrium at all times and the array of 
!@+   concentrations is updated accordingly.

!@+   waternacl is the ratio of wet mass to dry mass of a particle.  Instead
!@+   of calling a thermodynamic equilibrium code, this routine uses a
!@+   simple curve fit to estimate waternacl based on the current humidity.
!@+   The curve fit is based on ISORROPIA results for sodium sulfate
!@+   at 273 K.

!@auth Peter Adams, November 2001

      real*8 FUNCTION waternacl(rhe)
      IMPLICIT NONE

      real*8 rhe   !relative humidity (0-100 scale)
!@var   waternacl is the ratio of wet mass to dry mass of a particle.

      if (rhe .gt. 99.) rhe=99.
      if (rhe .lt. 1.) rhe=1.

         if (rhe .gt. 90.) then
            waternacl=5.1667642e-2*rhe**3-14.153121*rhe**2
     &               +1292.8377*rhe-3.9373536e4
         else
         if (rhe .gt. 80.) then
            waternacl=
     &      1.0629e-3*rhe**3-0.25281*rhe**2+20.171*rhe-5.3558e2
         else
         if (rhe .gt. 50.) then
            waternacl=
     &      4.2967e-5*rhe**3-7.3654e-3*rhe**2+.46312*rhe-7.5731
         else
         if (rhe .gt. 20.) then
            waternacl=
     &      2.9443e-5*rhe**3-2.4739e-3*rhe**2+7.3430e-2*rhe+1.3727
         else
            waternacl=1.17
         endif
         endif
         endif
         endif

         !check for error
         if (waternacl .gt. 45.) then
            write(*,*) 'ERROR in waternacl'
            write(*,*) rhe,waternacl
            call stop_model('ERROR in waternacl',255)
         endif

      RETURN
      END  FUNCTION waternacl


!@sum waterocil: This function uses the current RH to calculate how much water is 
!@+   in equilibrium with the hydrophillic OA.  Aerosol water concentrations 
!@+   are assumed to be in equilibrium at all times and the array of 
!@+   concentrations is updated accordingly.

!@+   waterocil is the ratio of wet mass to dry mass of a particle.  Instead
!@+   of calling a thermodynamic equilibrium code, this routine uses a
!@+   simple curve fit to estimate waterocil based on the current humidity.
!@+   The curve fit is based on observations of Dick et al. JGR D1 1471-1479

!@auth YUNHA LEE, AUG 2006

      real*8 FUNCTION waterocil(rhe)

      IMPLICIT NONE
      real*8 rhe   !relative humidity (0-100 scale)
      real*8 a, b, c, d, e, f, prefactor, activcoef
      parameter(a=1.0034, b=0.1614, c=1.1693,d=-3.1,
     & e=6.0)


      if (rhe .gt. 99.) rhe=99.
      if (rhe .lt. 1.) rhe=1.

      if (rhe .gt. 85.) then
         waterocil=d+e*(rhe/100) 
cyhl Growth factor above RH 85% is not available, so it assumes linear growth 
cyhl at above 85%.  
      else
         waterocil=a+b*(rhe/100)+c*(rhe/100)**2. 
cyhl This eq is based on the extrapolation curve obtained from  
cyhl Dick et al 2000 figure 5.(High organic,density=1400g/cm3)
      endif
      
         !check for error
      if (waterocil .gt. 10.) then
         write(*,*) 'ERROR in waterocil'
         write(*,*) rhe,waterocil
         call stop_model('ERROR in waterocil',255)
      endif

      RETURN
      END  FUNCTION waterocil

!@sum aeroupdate: update the aerosol water concentrations 
!@+      , fix size distributions, and check negative tracers
!@auth Peter Adams/YUNHA LEE

      SUBROUTINE aeroupdate

      USE DOMAIN_DECOMP_ATM, only : GRID, GET, write_parallel,
     &     am_i_root
      USE TOMAS_AEROSOL 
      USE GEOM, only: imaxj
      USE TRACER_COM, only : IDTSO4, IDTNA, IDTOCIL,IDTH2O,NBINS
     &     ,trm,IDTECOB,IDTECIL,IDTOCOB,IDTDUST,IDTNUMD,TRNAME
     *     ,ntm,ntm_TOMAS
      USE TRDIAG_COM, only : taijs=>taijs_loc !,taijls=>taijls_loc
      USE MODEL_COM, only : im,jm,lm     ! dimensions
     $                     ,q            ! saturated pressure
     $                     ,t
     $                     ,dtsrc
      USE CONSTANT,   only:  lhe
      USE DYNAMICS,   only: pmid,pk ! midpoint pressure in hPa (mb)

      IMPLICIT NONE

      INTEGER J_0, J_1, I_0, I_1

      integer i,j,l,k,n,jc,mpnum  !counters
      real*8 qsat         !used in RH calculation
      integer tracnum
      
      real*8 frac, Nkout(iBINS),Mkout(iBINS,icomp),Gcout(icomp-1)
      real*8 rhe
      real*8 waterso4, waternacl, waterocil

C****
C**** Extract useful local domain parameters from "grid"
C****
      CALL GET(grid, J_STRT=J_0,       J_STOP=J_1)
      I_0 = grid%I_STRT
      I_1 = grid%I_STOP


C     Loop over all grid cells
      DO L=1,LM                            
      DO J=J_0,J_1                          
      DO I=I_0,IMAXJ(J)

        temp = pk(l,i,j)*t(i,j,l) !should be in [K]
        rh = MIN(1.,q(i,j,l)/QSAT(temp,lhe,pmid(l,i,j))) ! rH [0-100%]
C     Swap GCM variables into aerosol algorithm variables
        do n=1,NBINS
          Nk(n)=trm(i,j,l,IDTNUMD-1+n)
          Mk(n,srtso4)=trm(i,j,l,IDTSO4-1+n)
          Mk(n,srtna )=trm(i,j,l,IDTNA -1+n)
          Mk(n,srtnh4)=0.1875*Mk(n,srtso4) ! artificial for now.. 0.0!t0m(i,j,l,IDTNH4-1+n)
          MK(n,srtecob)=trm(i,j,l,IDTECOB -1+n)
          MK(n,srtecil)=trm(i,j,l,IDTECIL -1+n)
          MK(n,srtocob)=trm(i,j,l,IDTOCOB -1+n)
          MK(n,srtocil)=trm(i,j,l,IDTOCIL -1+n) 
          MK(n,srtdust)=trm(i,j,l,IDTDUST -1+n) 
          Mk(n,srth2o)= trm(i,j,l,IDTH2O-1+n) !I don't think this is necessary!
        enddo

!     Do water eqm at appropriate times
        
        call storenm()
        call mnfix(Nk,Mk) 
        mpnum=7
        call aerodiag(mpnum,i,j,l)
        call ezwatereqm(Mk)
        
C     Swap Nk, Mk, and Gc arrays back to T0M
        do n=1,NBINS
          tracnum=IDTNUMD-1+n
          if (Nk(n) .ge. TRM(i,j,l,tracnum)) then
            TRM(i,j,l,tracnum)=Nk(n)
          else
            frac=Nk(n)/TRM(i,j,l,tracnum)
            call scalemom(i,j,l,tracnum,frac)
          endif
          do jc=1,icomp-idiag
            tracnum=IDTSO4-1+n+ibins*(jc-1)
            if (Mk(n,jc) .ge. TRM(i,j,l,tracnum)) then
              TRM(i,j,l,tracnum)=Mk(n,jc)
            else
              frac=Mk(n,jc)/TRM(i,j,l,tracnum)
              call scalemom(i,j,l,tracnum,frac)
            endif
          enddo
          tracnum=IDTH2O-1+n
          if (Mk(n,srth2o) .ge. TRM(i,j,l,tracnum)) then
            TRM(i,j,l,tracnum)=Mk(n,srth2o)
          else
            frac=Mk(n,srth2o)/TRM(i,j,l,tracnum)
            call scalemom(i,j,l,tracnum,frac)
          endif    
        enddo
        
C     Check for negative tracer problems
        do n=1,ntm_TOMAS
          if (TRM(i,j,l,IDTSO4+n-1) .lt. 0.0) then
            if (abs(TRM(i,j,l,IDTSO4+n-1)) .gt. 1.e-10) then
!serious problem - report error
              write(*,*) 'ERROR: Tracer ',trname(IDTSO4+n-1),
     &             trm(i,j,l,IDTSO4+n-1)
              write(*,*) ' < 0 in box ', i,j,l
              call stop_model('TRM<0 in aeroupdate',255)
c$$$  else
c$$$  !numerical problem - set to zero
c$$$  TRM(i,j,l,IDTSO4+n-1)=0.0!1.d-42 !5??
            endif
          endif
        enddo
        
      enddo
      enddo
      enddo
      
      RETURN
      END SUBROUTINE aeroupdate


!@sum scalemom: When a tracer concentration decreases, call this routine to
!@+   decrease T0M and all higher order moments in proportion.
!@auth Peter Adams, September 2000

      SUBROUTINE scalemom(i,j,l,tn,f)

      USE QUSDEF, only : nmom
      USE TRACER_COM, only : trm, trmom

      integer i,j,l       !grid box
      integer tn          !tracer id number
      real*8 f  !factor (0-1) by which to decrease all moments

      TRM(i,j,l,tn) = TRM(i,j,l,tn)*f

      do n=1,nmom
         trmom(n,i,j,l,tn)=f*trmom(n,i,j,l,tn)
      enddo
      
      RETURN
      END SUBROUTINE scalemom


!@sum storenm: Stores values of Nk and Mk into Nkd and Mkd for diagnostic
!@+   purposes.  Also do gas phase concentrations.
!@auth Peter Adams, June 2000

      SUBROUTINE storenm()

      USE TOMAS_AEROSOL
      IMPLICIT NONE
      integer j,k

      do j=1,icomp-1
         Gcd(j)=Gc(j)
      enddo
      do k=1,ibins
         Nkd(k)=Nk(k)
         do j=1,icomp
            Mkd(k,j)=Mk(k,j)
         enddo
      enddo

      RETURN
      END SUBROUTINE storenm


!@sum aerodiag : accumulates diagnostics on aerosol microphysical processes.
!@auth Peter Adams/Yunha Lee

      SUBROUTINE aerodiag(pt,i,j,l)

      USE TOMAS_AEROSOL
      USE TRACER_COM, only : IDTSO4,IDTNUMD,n_H2SO4
      implicit none 
      integer pt, i, j, l, jc, n
      integer k,kk,tracnum


      if(pt.eq.4.or.pt.eq.5.or.pt.eq.7 .or.pt.eq.8)then ! Aqoxid_mc,Aqoxid_ls,Aeroupdate,HETCHEM
C     Bulk species
        AEROD(i,j,l,n_H2SO4,pt)=
     &       (Gc(srtso4)-Gcd(srtso4))
        
        do n=1,ibins       
          
!     Aerosol number
          tracnum=IDTNUMD-1+n
          AEROD(i,j,l,tracnum,pt)= 
     &         (Nk(n)-Nkd(n))
!     Aerosol mass
          do jc=1,icomp-idiag
            tracnum=IDTSO4-1+n+ibins*(jc-1)
            AEROD(i,j,l,tracnum,pt)= 
     &           (Mk(n,jc)-Mkd(n,jc))
            
          enddo   
        enddo
        
      else
        
C Bulk species
        AEROD(i,j,l,n_H2SO4,pt)=AEROD(i,j,l,n_H2SO4,pt)
     &       +(Gc(srtso4)-Gcd(srtso4))
        
        do n=1,ibins  
!     Aerosol number
          tracnum=IDTNUMD-1+n
          
          AEROD(i,j,l,tracnum,pt)= 
     &         AEROD(i,j,l,tracnum,pt) 
     &         +(Nk(n)-Nkd(n))  
          
!     Aerosol mass
          do jc=1,icomp-idiag
            tracnum=IDTSO4-1+n+ibins*(jc-1)
            AEROD(i,j,l,tracnum,pt)= 
     &           AEROD(i,j,l,tracnum,pt) 
     &           + (Mk(n,jc)-Mkd(n,jc))
            
          enddo 
        enddo
      endif      
      RETURN
      END SUBROUTINE aerodiag


!@sum subgridcoag_drv: subgrid coagulation for 3D emissions  
!@+   No moments updated here because aerosol emission are positive! 
!@auth Jeff Pierce/Yunha Lee, July 2011 

      SUBROUTINE subgridcoag_drv(dtstep)


      USE DOMAIN_DECOMP_ATM, only : GRID, GET, write_parallel
     &     ,am_i_root
      USE TOMAS_AEROSOL 
      USE TRDIAG_COM, only : taijs=>taijs_loc,ijts_subcoag,itcon_subcoag
      USE FLUXES, only: tr3Dsource
      USE MODEL_COM, only : lm     ! dimensions
     $                     ,t            ! potential temperature (C)
     $                     ,q            ! saturated pressure
     $                     ,dtsrc
      USE GEOM, only: axyp,imaxj,BYAXYP,axyp
      USE CONSTANT,   only: pi,gasc,mair   
      USE DYNAMICS,   only: pmid,pk, am   ! midpoint pressure in hPa (mb)
!                                           and pk is t mess up factor
!                                           BYAM  1/Air mass (m^2/kg)
      USE TRACER_COM, only : nbins,xk,ntm,trm,trmom,ntsurfsrc,
     &     IDTSO4,IDTNA,IDTECOB,IDTECIL,IDTOCOB,
     &     IDTOCIL,IDTDUST,IDTNUMD,n_SO2,IDTH2O
      IMPLICIT NONE

      integer :: J_1, J_0, I_1, I_0
      INTEGER :: L,I,J

      REAL*8, INTENT(IN) :: dtstep
      INTEGER n,ns,c,k,tot_src,tracnum
      INTEGER tomas_ntsurf !same as ntsurfsrc
      real*8 ndistinit(nbins) !the number of particles being added to the gridbox before subgrid coag
      real*8,dimension(nbins) ::  ndist, ndist2, ndist0 !the number of particles in the box
      real*8,dimension(nbins,icomp) :: mdist,mdist2,mdist0 ! the mass of each component in the box. (kg)
      real tscale ! the scale time for mixing (s)
      real*8 ndistfinal(nbins),tot_ndistinit(nbins) !the number of particles being added to the gridbox after subgrid coag
      real*8 maddfinal(nbins) !the mass that should be added to each bin due to coagulation (kg)


      CALL GET(grid, J_STRT=J_0, J_STOP=J_1)
      I_0 = grid%I_STRT
      I_1 = grid%I_STOP

      DO L=1,LM; DO J=J_0,J_1; DO I=I_0,imaxj(j)
        
        tscale=5.*3600.
        temp = pk(l,i,j)*t(i,j,l) !should be in [K]
        pres= pmid(l,i,j)*100.  ! pmid in [hPa]
        boxvol=am(l,i,j)*axyp(i,j)/mair*1000.d0
     &       *gasc*temp/pres*1e6 !cm3
        
        do k=1,nbins
          ndist0(k)=TRM_EMIS(I,J,L,IDTNUMD+k-1)
          do c=1,icomp-idiag
            mdist0(k,c)=TRM_EMIS(I,J,L,IDTSO4+(c-1)*nbins+k-1)
          enddo
          mdist0(k,srtnh4)=0.0
          mdist0(k,srth2o)=TRM_EMIS(I,J,L,IDTH2O+k-1)
          ndistfinal(k)=0
          maddfinal(k)=0
        enddo
        
        ndist(:)=0.
        mdist(:,:)=0.
        ndistinit(:)=0.
        
        DO NS=1,3
          
          do k=1,nbins
            if(ns.lt.3) 
     &           ndistinit(k)=tr3Dsource(i,j,l,ns,IDTNUMD+K-1)*dtstep
            if(ns.eq.3) 
     &           ndistinit(k)=tr3Dsource(i,j,l,ns+1,IDTNUMD+K-1)*dtstep 
          enddo
          
          if(sum(ndistinit(1:nbins)).gt.0.)then
!     only when there is emission! 
            call subgridcoag(ndistinit,ndist0,mdist0,boxvol,
     &           tscale,ndistfinal,maddfinal) ! account for subgrid coagulation
            
            do k=1,nbins
              ndist(k)=ndist(k)+ndistfinal(k)
              
              IF(NS.EQ.1)THEN
!SO4
                mdist(k,srtso4)=
     &               ndistfinal(k)*(sqrt(xk(k)*xk(k+1)))+
     &               maddfinal(k)
                
              ELSEIF(NS.EQ.2)THEN
!EC
                mdist(k,srtecil)=
     &               ndistfinal(k)*0.2*(sqrt(xk(k)*xk(k+1)))+
     &               maddfinal(k)*0.2
                mdist(k,srtecob)=
     &               ndistfinal(k)*0.8*(sqrt(xk(k)*xk(k+1)))+
     &               maddfinal(k)*0.8   
                
              ELSEIF(NS.EQ.3)THEN
!OC
                mdist(k,srtocil)=
     &               ndistfinal(k)*0.5*(sqrt(xk(k)*xk(k+1)))+
     &               maddfinal(k)*0.5
                mdist(k,srtocob)=
     &               ndistfinal(k)*0.5*(sqrt(xk(k)*xk(k+1)))+
     &               maddfinal(k)*0.5   
              ENDIF
            enddo               !k
          endif
        enddo                   !ns
        
!     fix the inconsistancies in the distribution
        do k=1,nbins
          ndist2(k)=ndist(k)+ndist0(k)
          do c=1,icomp
            mdist2(k,c)=mdist(k,c)+mdist0(k,c)
          enddo
        enddo
        
        if(sum(ndist(1:nbins)).gt.0.) call mnfix(ndist2,mdist2)
              
        do k=1,nbins            
          tracnum=IDTNUMD-1+k  
          N_subgridcg(i,j,l,k,2)=(ndist2(k)- !this is emission after subgrid
     &         trm(i,j,l,tracnum)) 

          if(l.eq.1)then
            N_subgridcg(i,j,l,k,2)=N_subgridcg(i,j,l,k,2)+ 
     &           N_subgridcg(i,j,l,k,1) !from 2-d emission subgrid coagulation
          endif
          
          trm(i,j,l,tracnum)=ndist2(k)          
          taijs(i,j,ijts_subcoag(tracnum)) 
     &         =taijs(i,j,ijts_subcoag(tracnum))
     &         +N_subgridcg(i,j,l,k,2) ! /adt
          
          if (itcon_subcoag(tracnum).gt.0) 
     &         call inc_diagtcb(i,j,N_subgridcg(i,j,l,k,2) ,
     &         itcon_subcoag(tracnum),tracnum)
          
          do c=1,icomp-idiag            
            tracnum=IDTSO4-1+k+nbins*(c-1) 
            M_subgridcg(i,j,l,k,c,2)=mdist2(k,c)- !trm + emission after subgrid 
     &           trm(i,j,l,tracnum) !trm + emission before subgrid (which is computed in apply_tracer3d)

          if(l.eq.1)then
            M_subgridcg(i,j,l,k,c,2)=M_subgridcg(i,j,l,k,c,2)+
     &           M_subgridcg(i,j,l,k,c,1)
          endif

            trm(i,j,l,tracnum)=mdist2(k,c)
            taijs(i,j,ijts_subcoag(tracnum)) 
     &           =taijs(i,j,ijts_subcoag(tracnum))
     &           +M_subgridcg(i,j,l,k,c,2) ! /adt

            if (itcon_subcoag(tracnum).gt.0) 
     &           call inc_diagtcb(i,j,M_subgridcg(i,j,l,k,c,2),
     &           itcon_subcoag(tracnum),tracnum)
                       
          enddo
        enddo

        M_subgridcg(i,j,l,:,:,:)=0.0
        N_subgridcg(i,j,l,:,:)=0.0
        
       enddo; enddo; enddo
        
       return
       end subroutine subgridcoag_drv
    

!@sum subgridcoag_drv_2D: subgrid coagulation for 2D emissions  
!@+   No moments updated here because aerosol emission are positive! 
!@auth Jeff Pierce/Yunha Lee, July 2011 

      SUBROUTINE subgridcoag_drv_2D(dtstep)

      USE DOMAIN_DECOMP_ATM, only : GRID, GET, write_parallel
     &     ,am_i_root
      USE TOMAS_AEROSOL 
      USE MODEL_COM, only : lm     ! dimensions
     $                     ,t            ! potential temperature (C)
     $                     ,q            ! saturated pressure
     $                     ,dtsrc
      USE GEOM, only: axyp,imaxj,BYAXYP,axyp
      USE CONSTANT,   only: pi,gasc,mair   
      USE DYNAMICS,   only: pmid,pk, am   ! midpoint pressure in hPa (mb)
!                                           and pk is t mess up factor
!                                           BYAM  1/Air mass (m^2/kg)
      USE TRACER_COM, only : nbins,xk,ntm,trm,trmom,ntsurfsrc,
     &     IDTSO4,IDTNA,IDTECOB,IDTECIL,IDTOCOB,
     &     IDTOCIL,IDTDUST,IDTNUMD,n_SO2,IDTH2O
      USE FLUXES, only : trsource,trflux1, trsrfflx
      USE TRDIAG_COM, only : taijs=>taijs_loc
      IMPLICIT NONE

      integer :: J_1, J_0, I_1, I_0,L,I,J
      REAL*8, INTENT(IN) :: dtstep
      INTEGER n,ns,c,k,tot_src,tracnum
      INTEGER tomas_ntsurf !same as ntsurfsrc
      real*8 ndistinit(nbins) !the number of particles being added to the gridbox before subgrid coag
      real*8,dimension(nbins) ::  ndist,ndist0 !the number of particles in the box
      real*8,dimension(nbins,icomp) :: mdist,mdist0 ! the mass of each component in the box. (kg)
      real tscale ! the scale time for mixing (s)
      real*8 ndistfinal(nbins),tot_ndistinit(nbins) !the number of particles being added to the gridbox after subgrid coag
      real*8 maddfinal(nbins) !the mass that should be added to each bin due to coagulation (kg)


      CALL GET(grid, J_STRT=J_0, J_STOP=J_1)
      I_0 = grid%I_STRT
      I_1 = grid%I_STOP
      l=1 !2-D emission is at 1st layer

      DO J=J_0,J_1; DO I=I_0,imaxj(j)
          
!     subgrid timescale and met conditions
        tscale=5.*3600.
        temp = pk(l,i,j)*t(i,j,l) !should be in [K]
        pres= pmid(l,i,j)*100.  ! pmid in [hPa]
        boxvol=am(l,i,j)*axyp(i,j)/mair*1000.d0
     &       *gasc*temp/pres*1e6 !cm3
        
!     Amount of tracer before emission is applied.         
        do k=1,nbins
          ndist0(k)=TRM(I,J,L,IDTNUMD+k-1)
          do c=1,icomp-idiag
            mdist0(k,c)=TRM(I,J,L,IDTSO4+(c-1)*nbins+k-1)
          enddo
          mdist0(k,srtnh4)=0.0
          mdist0(k,srth2o)=TRM(I,J,L,IDTH2O+k-1)
          ndistfinal(k)=0
          maddfinal(k)=0
        enddo
              
!     Only initialize when 2-D emission starts! 
        tot_ndistinit(:)=0.
        ndist(:)=0.
        mdist(:,:)=0.
        ndistinit(:)=0.
                
        DO ns=1,ntsurfsrc(idtnumd)         
!     ns=1 for so4; ns=2 for ec; ns=3 for oc          
          do k=1,nbins            
            ndistinit(k)=trsource(i,j,NS,IDTNUMD+K-1)*dtstep           
            tot_ndistinit(k)=tot_ndistinit(k)+ndistinit(k) !sum of number emission for SO4, EC, and OC
          enddo
          
!     only when there is emission!   
          if(sum(ndistinit(1:nbins)).gt.0.)then        
            call subgridcoag(ndistinit,ndist0,mdist0,boxvol,
     &           tscale,ndistfinal,maddfinal) ! account for subgrid coagulation
            
            do k=1,nbins
              ndist(k)=ndistfinal(k)+ ndist(k) !sum of all emission type after subgrid coag              
              IF(NS.EQ.1)       !SO4
     &             mdist(k,srtso4)= 
     &             ndistfinal(k)*(sqrt(xk(k)*xk(k+1)))+
     &             maddfinal(k)
              
              IF(NS.EQ.2)       !EC
     &             mdist(k,srtecil)= 
     &             ndistfinal(k)*0.2*(sqrt(xk(k)*xk(k+1)))+
     &             maddfinal(k)*0.2                
              IF(NS.EQ.2)       !EC
     &             mdist(k,srtecob)= 
     &             ndistfinal(k)*0.8*(sqrt(xk(k)*xk(k+1)))+
     &             maddfinal(k)*0.8   
              
              IF(NS.EQ.3)       !OC
     &             mdist(k,srtocil)=
     &             ndistfinal(k)*0.5*(sqrt(xk(k)*xk(k+1)))+
     &             maddfinal(k)*0.5
              IF(NS.EQ.3)       !OC
     &             mdist(k,srtocob)=
     &             ndistfinal(k)*0.5*(sqrt(xk(k)*xk(k+1)))+
     &             maddfinal(k)*0.5 

            enddo
          endif !positive emission         
        enddo !ntsurfsrc    
        
        if(sum(ndist(1:nbins)).gt.0.) call mnfix(ndist,mdist)
        
!     DIAGNOSTICS!       
        do k=1,nbins            
          tracnum=IDTNUMD-1+k           
          N_subgridcg(i,j,l,k,1)=N_subgridcg(i,j,l,k,1)+ndist(k)- ! emission after subgrid
     &         tot_ndistinit(k) ! emission before subgrid
          
          trflux1(i,j,tracnum)=trflux1(i,j,tracnum)+
     &         (ndist(k)- tot_ndistinit(k))/dtstep !kg/sec
          
          do c=1,icomp-idiag            
            tracnum=IDTSO4-1+k+nbins*(c-1)
            
            if(c.eq.2.or.c.eq.7)then !no subgrid coagulation
              M_subgridcg(i,j,l,k,c,1) =0.
            else              
              if(c.eq.1) tomas_ntsurf=ntsurfsrc(n_SO2)
              if(c.eq.3.or.c.eq.4) tomas_ntsurf=ntsurfsrc(IDTECOB) !ecob
              if(c.eq.5.or.c.eq.6)  tomas_ntsurf=ntsurfsrc(IDTOCOB) !ecob              
              M_subgridcg(i,j,l,k,c,1)=M_subgridcg(i,j,l,k,c,1)
     &             + mdist(k,c)-
     &             (sum(trsource(i,j,1:tomas_ntsurf,tracnum))*dtstep) ! emission before subgrid              
            endif
            trflux1(i,j,tracnum)=trflux1(i,j,tracnum)+
     &           (mdist(k,c)-
     &             (sum(trsource(i,j,1:tomas_ntsurf,tracnum))*dtstep))
     &          /dtstep 
            
          enddo !c
        enddo ! k

      enddo; enddo ! i,j
      
      return
      end subroutine subgridcoag_drv_2D
     

!@sum subgridcoag: determine how much of each size of freshly emitted aerosol will 
!@+   be scavenged by coagulation prior to being completely mixed in the gridbox and will
!@+   give the new emissions size distribution along with where the mass of coagulated
!@+   particles should be added.

!@auth Jeff Pierce, Dec 2006

      SUBROUTINE subgridcoag(ndistinit,ndist,mdist,boxvolume,
     & tscale,ndistfinal,maddfinal)

      USE TRACER_COM, only : nbins,xk
      USE TOMAS_AEROSOL
      USE CONSTANT, ONLY : pi,gasc,mair
      IMPLICIT NONE

      INTEGER n,k,c,kk
      real*8, intent(in) :: ndistinit(nbins) !the number of particles being added to the gridbox before subgrid coag
      real*8, intent(in) :: ndist(nbins) !the number of particles in the box
      real*8, intent(in) :: mdist(nbins,icomp) ! the mass of each component in the box. (kg)
      real*8, intent(in) :: boxvolume  ! volume of box in cm3
      real, intent(in) :: tscale ! the scale time for mixing (s)
      real*8, intent(out) :: ndistfinal(nbins) !the number of particles being added to the gridbox after subgrid coag
      real*8, intent(out) :: maddfinal(nbins) !the mass that should be added to each bin due to coagulation (kg)
      real*8 mp ! mass of the particle (kg)
      real density                !density (kg/m3) of particles
      real*8 diameter(nbins) ! diamter of the particle (m)
      real*8 diaml(nbins) ! total diamter of particles larger (m/cm3)
      real*8 fracdiaml(nbins,nbins) ! fraction of coagulation that occurs with each bin larger
      real*8 kcoag(nbins) ! the coagulation rate for the particles in each bin (s^-1)
      real aerodens
      real mso4, mh2o, mno3, mnh4  !mass of each component (kg/grid box)
      real mecil,mecob,mocil,mocob
      real mdust,mnacl   
      real*8 v1,v2,v3  !for coag rate calculation
      parameter(v1=8.5708E-13,
     &          v2=-1.4174,
     &           v3=4.3047E-4)

C     get the wet diameter of particles in each size bin
      do k=1,nbins
         mp=0.1875*mdist(k,srtso4)
         do c=1,icomp
            mp = mp + mdist(k,c)
         enddo
         if (ndist(k).eq.0.)then
            mp=sqrt(xk(k)*xk(k+1))
         else
            mp = mp / ndist(k)
         endif
         if((mdist(k,srtso4)+mdist(k,srtna)+mdist(k,srtocil)+
     &        mdist(k,srtdust)).eq.0)then
            density=1400.
         else
         mso4=mdist(k,srtso4) 
         mnacl=mdist(k,srtna)
         mno3=0.e0
         if ((mso4+mno3) .lt. 1.e-8) mso4=1.e-8
         mnh4=0.1875*mso4  !assume ammonium bisulfate
         mecob=mdist(k,srtecob)
         mecil=mdist(k,srtecil)
         mocil=mdist(k,srtocil)
         mocob=mdist(k,srtocob)
         mdust=mdist(k,srtdust)          
         mh2o=mdist(k,srth2o)   

         density=aerodens(mso4,mno3,mnh4 !mno3 taken off!
     *        ,mnacl,mecil,mecob,mocil,mocob,mdust,mh2o) !assume bisulfate 
         endif
         diameter(k)=2.*(3./4./pi*mp/density)**(1./3.) ! m
      enddo

C     get the total diameter of particles larger than each size bin
      diaml(nbins)=0. !no diameter larger than largest bin
      do kk=1,nbins-1
         k=nbins-kk
         diaml(k) = diaml(k+1) + diameter(k+1)*ndist(k+1)/boxvolume ! m/cm3
      enddo
      
C     get the fraction of the diameter larger that comes from each bin larger
      do k=1,nbins
         do kk=1,nbins
            fracdiaml(k,kk)=0.
         enddo
      enddo
      do k=1,nbins-1
         do kk=k+1,nbins
            if (diaml(k).ne.0.0)then
               fracdiaml(k,kk)=diameter(kk)*ndist(kk)/boxvolume/diaml(k)
            else
               fracdiaml(k,kk)=0.0
            endif
         enddo
      enddo

C     determine the coagulation rate for each size bin
      do k=1,nbins
         if (diameter(k).gt.0.d0)then
            kcoag(k) = (v1*diameter(k)**(v2)+v3)*diaml(k)
         else
            kcoag(k) = 0.d0
         endif
      enddo

C     determine the number of new particles left after coagulation
      do k=1,nbins
         ndistfinal(k)=ndistinit(k)*exp(-kcoag(k)*tscale)
      enddo

C     determine the mass added to each bin coagulation
      do k=1,nbins
         maddfinal(k)=0.
      enddo
      do k=1,nbins-1
         do kk=k+1,nbins
            maddfinal(kk)=maddfinal(kk) + (ndistinit(k)-ndistfinal(k))*
     &           fracdiaml(k,kk)*sqrt(xk(k)*xk(k+1))
         enddo
      enddo

      return
      end SUBROUTINE subgridcoag

!@sum  alloc_tracer_TOMAS_com : alllocate arrays whose sizes need 
!@+    to be determined at run-time
!@auth Yunha Lee
      subroutine alloc_tracer_TOMAS_com(grid)

      use domain_decomp_atm, only : dist_grid, get
      use model_com, only     : lm
      use TOMAS_aerosol

      IMPLICIT NONE

      type (dist_grid), intent(in) :: grid
      integer :: ier, J_1H, J_0H, I_1H, I_0H
      logical :: init = .false.

      if(init)return
      init=.true.
    
      call get( grid , J_STRT_HALO=J_0H, J_STOP_HALO=J_1H )
      I_0H=GRID%I_STRT_HALO
      I_1H=GRID%I_STOP_HALO 

      allocate(  AQSO4oxid_mc(I_0H:I_1H,J_0H:J_1H,LM)   )
      allocate(  AQSO4oxid_ls(I_0H:I_1H,J_0H:J_1H,LM)   )
      allocate(  H2SO4_chem(I_0H:I_1H,J_0H:J_1H,LM)  )

#ifdef TOMAS_HETCHEM
      allocate(  H2SO4_hetchem(I_0H:I_1H,J_0H:J_1H,LM,NBINS)  )
#endif
      allocate( AEROD(I_0H:I_1H,J_0H:J_1H,LM,NTM,ptype) )
     
      allocate(  N_subgridcg(I_0H:I_1H,J_0H:J_1H,LM,IBINS,2) )
      allocate(  M_subgridcg(I_0H:I_1H,J_0H:J_1H,LM,IBINS,
     *     ICOMP-IDIAG,2))
     
      allocate(  TRM_EMIS(I_0H:I_1H,J_0H:J_1H,LM,NTM) )

      allocate(  CCN_TOMAS(I_0H:I_1H,J_0H:J_1H,LM,NSMAX) )

      return
      end subroutine alloc_tracer_TOMAS_com
