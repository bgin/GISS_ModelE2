#include "rundeck_opts.h"

      module AerParam_mod
!@sum This module reads, time-interpolates, and stores fields needed
!@+   by the radiation code in the prescribed-aerosol configuration
!@+   of modelE.  Subroutine updateAerosol2 provides the fields
!@+   used to calculate the direct radiative effects of aerosols.
!@+   Subroutine dCDNC_est provides a parameterized estimate of
!@+   changes of lower-tropospheric CDNC relative to 1850, as input
!@+   to prescriptions of aerosol indirect effects on cloud cover
!@+   and optical depth.
!@+   Dust aerosols are currently ingested via a separate module.
!@auth D. Koch, R. Ruedy
!@auth T. Clune reduced memory usage, created first ver. of readTable()
!@auth M. Kelley updated for grid flexibility, added comments, netcdf input

#ifdef NEW_IO_AERINP
      use timestream_mod, only : timestream
#endif
      implicit none
      save
      private

      public :: dCDNC_est
      public :: updateAerosol
      public :: updateAerosol2
      public :: DRYM2G, aermix
      public :: lma

      real*8, allocatable :: anssdd(:,:)
      real*8, allocatable :: mdpi(:,:,:)
      real*8, allocatable :: mdcur(:,:,:)
      real*8, allocatable :: md1850(:,:,:,:)

      character(len=3), dimension(6), parameter :: aernames=(/
!**** Sulfate
     &     'SUL',
!**** Sea Salt
     &     'SSA',
!**** Nitrate
     &     'NIT',
!**** Organic Carbon
     &     'OCA',
!**** Black Carbon from Fossil and bio fuel
     &     'BCA',
!**** Black Carbon from Biomass burning
     &     'BCB'
     &     /)

      real*8 :: DRYM2G(8) =
     &     (/4.667, 0.866, 4.448, 5.017, 9.000, 9.000, 1.000,1.000/)

C     Layer  1    2    3    4    5    6    7    8    9
      INTEGER :: La720=3 ! top low cloud level (aerosol-grid).
                         ! =3 for orig 9-level model
      REAL*8 , PARAMETER ::
     &     Za720=2635.                ! depth of low cloud region (m)
     &    ,byz_cm3 = 1.d-6 / Za720    ! 1d-6/depth in m (+conversion /m3 -> /cm3)
     &    ,byz_gcm3 = 1.d-3 * byz_cm3 ! g vs kg

      REAL*8, dimension(13) :: AERMIX=(/
C      Pre-Industrial+Natural 1850 Level  Industrial Process  BioMBurn
C      ---------------------------------  ------------------  --------
C       1    2    3    4    5    6    7    8    9   10   11   12   13
C      SNP  SBP  SSP  ANP  ONP  OBP  BBP  SUI  ANI  OCI  BCI  OCB  BCB
     + 1.0, 1.0, 1.0, 1.0, 2.5, 2.5, 1.9, 1.0, 1.0, 2.5, 1.9, 2.5, 1.9/)

      integer :: ima, jma, lma

#ifdef NEW_IO_AERINP
!@var A6streams interface for reading and time-interpolating AERO files
!@var BCdepstream interface for reading and time-interpolating BC_dep file
!@+   See usage notes in timestream_mod
      type(timestream), dimension(6) :: A6streams
      type(timestream), public :: BCdepstream
!@dbparam A6yr data year to use for each aerosol.  Rundecks should set
!@+       this array via XXX_yr where XXX is an element of aernames(:).
!@+       Elements of A6yr not specified in the rundeck will default to
!@+       the (possibly time-varying) value of jYearA passed to
!@+       updateaerosol2.
      integer, dimension(6) :: A6yr
#else
      REAL*4, allocatable :: A6YEAR2(:,:,:,:,:)
      REAL*4, dimension(:,:,:,:,:), allocatable ::
     &     SULDD,NITDD,OCADD,BCADD,BCBDD
      REAL*4, dimension(:,:,:,:), allocatable :: SSADD
      REAL*8, allocatable :: anfix(:,:,:)
#endif

!@var depoBC,depoBC_1990 prescribed black carbon deposition (curr,1990)
!@+   for parameterization of the BC effect on snow albedo
      REAL*8, ALLOCATABLE, DIMENSION(:,:), public :: depoBC,depoBC_1990

      contains

      subroutine dCDNC_EST(i,j,pland, dCDNC) !, table)
!@sum  finds change in cloud droplet number concentration since 1850
!@auth R. Ruedy
!@ver  1.0
      USE CONSTANT, only : pi
      implicit none
      integer, intent(in)  :: i,j ! grid indices
      real*8 , intent(in)  :: pland ! land fraction
      real*8 , intent(out) :: dCDNC ! CDNC(cur)-CDNC(1850)

      real*8, parameter, dimension(5) ::
C                TROPOSPHERIC AEROSOL PARAMETERS
C                  SO4     NO3    OCX    BCB   BCI
     &  f_act=(/ 1.0d0,  1.0d0, 0.8d0, 0.6d0, .8d0/), ! soluble fraction
     &  dens =(/1769d0, 1700d0,  1.d3,  1.d3, 1.d3/)  ! density

      real*8, parameter, dimension(2) ::
C                    Ocean         Land      ! r**3: r=.085,.052 microns
     &  radto3 =(/ 614.125d-24, 140.608d-24/),  ! used for SO4,NO3,OC,BC
     &  scl    =(/     162d0,       298d0/),  ! for Gultepe's formula
     &  offset =(/     273d0,       595d0/)   ! for Gultepe's formula

      integer it, n
      real*8  An,An0,cdnc(2),cdnc0(2),fbymass1

      do it=1,2  ! ocean, land
        An0 = anssdd(i,j)  !  aerosol number of sea salt and dust
        An  = An0          !  aerosol number of sea salt and dust
        do n=1,4
          fbymass1 =  F_act(n)*(.75d0/pi)/(dens(n)*radto3(it))
          An0 = An0 + mdpi (n,i,j)*fbymass1   ! +fact*tot_mass/part_mass
          An  = An  + mdcur(n,i,j)*fbymass1
        end do
        fbymass1 =  F_act(5)*(.75d0/pi)/(dens(5)*radto3(it))
        An  = An  + mdcur(5,i,j)*fbymass1

        if(An0.lt.1.) An0=1.
        if(An .lt.1.) An =1.
        cdnc0(it) = max( 20d0, scl(it)*log10(AN0)-offset(it))
        cdnc (it) = max( 20d0, scl(it)*log10(AN )-offset(it))
      end do

      dCDNC = (1-pland)*(cdnc(1)-cdnc0(1))+pland *(cdnc(2)-cdnc0(2))
      return
      end subroutine dCDNC_EST

      SUBROUTINE updateAerosol(JYEARA,JJDAYA, a6jday, plbaer)

cc    INCLUDE 'rad00def.radCOMMON.f'
C     ------------------------------------------------------------------
C     Reads: sep2003_XXX_Koch_kg_m2_72x46x9_1850-1990 aerosol kg/m2 data
C     for SUI,OCI,BCI, and PRE (PRE=SNP,SBP,SSP,ANP,ONP,OBP,ANI,OCB,BCB)
C
C     Makes: A6YEAR(72,46,9,0:12,6), A6JDAY(9,6,72,46) (dry aerosol Tau)
C     ------------------------------------------------------------------

      USE FILEMANAGER, only : openunit,closeunit
      implicit none

      INTEGER, intent(in) :: jyeara,jjdaya
      real*8, pointer :: a6jday(:,:,:,:)
      real*8, dimension(:), pointer ::  plbaer

      real*4, allocatable, dimension(:,:,:,:,:) :: A6YEAR,
     &     PREDD, SUIDD, OCIDD, BCIDD
      REAL*8  md1850(4,72,46,0:12),anfix(72,46,0:12)
      save A6YEAR,PREDD,SUIDD,OCIDD,BCIDD,md1850,anfix ! ,mddust

      CHARACTER*80 XTITLE
      CHARACTER*40 :: RDFILE(5) = (/                !  Input file names
     1            'sep2003_PRE_Koch_kg_m2_ChinSEA_Liao_1850'
     2           ,'sep2003_SUI_Koch_kg_m2_72x46x9_1875-1990'
     3           ,'sep2003_OCI_Koch_kg_m2_72x46x9_1875-1990'
     4           ,'sep2003_BCI_Koch_kg_m2_72x46x9_1875-1990'
     5           ,'low.dust.72x46.monthly.bin              '/)

      CHARACTER*40 :: RDFGEN(5) = (/                ! generic names
     * 'TAero_PRE','TAero_SUI','TAero_OCI','TAero_BCI','M_LowDust'/)

C                TROPOSPHERIC AEROSOL COMPOSITIONAL/TYPE PARAMETERS
C                   SO4    SEA    ANT    OCX    BCI    BCB   *BCB  *BCB
C     DATA REFDRY/0.200, 1.000, 0.300, 0.300, 0.100, 0.100, 0.200,0.050/
C
C     DATA REFWET/0.272, 1.808, 0.398, 0.318, 0.100, 0.100, 0.200,0.050/
C
C     DATA DRYM2G/4.667, 0.866, 4.448, 5.018, 9.000, 9.000, 5.521,8.169/
C
CKoch DATA DRYM2G/5.000, 2.866, 8.000, 8.000, 9.000, 9.000, 5.521,8.169/
C
C     DATA RHTMAG/1.788, 3.310, 1.756, 1.163, 1.000, 1.000, 1.000,1.000/
C
CRH70 DATA WETM2G/8.345, 2.866, 7.811, 5.836, 9.000, 9.000, 5.521,8.169/
C
C     DATA Q55DRY/2.191, 2.499, 3.069, 3.010, 1.560, 1.560, 1.914,0.708/
C
C     DATA DENAER/1.760, 2.165, 1.725, 1.500, 1.300, 1.300, 1.300,1.300/
C
C     ------------------------------------------------------------------
C          DRYM2G(I) = 0.75/DENAER(I)*Q55DRY(I)/REFDRY(I)
C          WETM2G(I) = DRYM2G(I)*RHTMAG(I)
C          RHTMAG(I) = Rel Humidity TAU Magnification factor  at RH=0.70
C          REFWET(I) = Rel Humidity REFDRY Magnification      at RH=0.70
C     ------------------------------------------------------------------

C     TROP AEROSOL 1850 BACKGROUND, INDUSTRIAL & BIO-BURNING PARAMETERS
C     DATA AERMIX/
C       Pre-Industrial+Natural 1850 Level  Industrial Process  BioMBurn
C       ---------------------------------  ------------------  --------
C        1    2    3    4    5    6    7    8    9   10   11   12   13
C       SNP  SBP  SSP  ANP  ONP  OBP  BBP  SUI  ANI  OCI  BCI  OCB  BCB
C    +  1.0, 1.0, 1.0, 1.0, 2.5, 2.5, 1.9, 1.0, 1.0, 2.5, 1.9, 2.5, 1.9/

C      A6YEAR          PRE                  SUI         OCI        BCI
C     ------------------------------------------------------------------
C     NAER=1=SO4 = SNP*1+SBP*2         +  SUI*I,J
C          2=SEA = SSP*3
C          3=ANT = ANP*4+ANI*0,8
C          4=OCX = ONP*5+OBP*6+OCB*0,9             +  OCI*I,J
C          5=BCI =                                               BCI*I,J
C          6=BCB = BBP*7,BCB*0,10
C     ------------------------------------------------------------------
C     Aerosol input data is from designated source files PRE SUI OCI BCI
C     Aerosol output is accumulated for 6-A6YEAR designated compositions
C           SNP*1 represents AERMIX(1)*PRE(I,J,L,M,1) = 1850 Natural SO4
C           SBP*2 represents AERMIX(2)*PRE(I,J,L,M,2) = 1850 BioBurn SO4
C         SUI*I,J represents AERMIX(8)*(SUI(I,J,L,M,I)intSUI(I,J,L,M,J))
C           SSP*1 represents AERMIX(3)*PRE(I,J,L,M,3) = 1850 SeaSalt
C         BCB*0,10represents AERMIX(11)*(0_interpol_PRE(I,J,L,M,10)) BCB
C         (which is interpolated linearly in time from 0 amount in 1850)
C
C         1850 Background   Sulfate  SO4 = 0.870 Natural + 0.130 BioBurn
C         1850 Background  Sea Salt  SEA =  all Natural-Mean SSP SeaSalt
C         1850 Background AmmNitrate ANT = ANP=(1.26/5.27)*1990(ANP+ANI)
C         1850 Background Org Carbon OCX = 0.162 Natural + 0.838 BioBurn
C         1850 Background Blk Carbon BCI = 0  (No Industrial BC in 1850)
C         1850 Background Blk Carbon BCB = BBP,all of 1850 BC is BioBurn
C     ------------------------------------------------------------------
      logical qexist
      INTEGER, save :: IFILE=11, IFIRST=1, JYRNOW=0

      INTEGER ia,idd,ndd,m,mi,mj,i,j,l,n,jyearx,iys,jys,iyc,jyc
      REAL*8 WTANI,WTOCB,WTBCB,wt75,swti,swtj,cwti,cwtj,xmi,wtmi,wtmj
      !REAL*8 , PARAMETER :: Za720=2635. ! depth of low cloud region (m)
      REAL*8 xsslt,byz ! ,xdust
      IF(IFIRST==1) THEN

        allocate(A6YEAR(72,46,9,0:12,6))
        allocate(PREDD(72,46,9,12,10))
        allocate(SUIDD(72,46,9,12,8))
        allocate(OCIDD(72,46,9,12,8))
        allocate(BCIDD(72,46,9,12,8))
        allocate(A6JDAY(9,6,72,46))

        allocate(anssdd(72,46))
        allocate(mdpi(4,72,46))
        allocate(mdcur(5,72,46))

C                                       READ Input PRE,SUI,OCI,BCI Files
C                                       --------------------------------
      inquire (file=RDFGEN(1),exist=qexist) ! decide whether specific or
      if(qexist) RDFILE=RDFGEN              !     generic names are used
      inquire (file=RDFILE(1),exist=qexist) !     stop if neither exist
      if(.not.qexist)
     &     call stop_model('updateAerosol: no TropAero files',255)

!**** Pre-industrial data
      call openunit (RDFILE(1),ifile,.true.,.true.)    ! unformatted,old
      DO 101 IDD=1,10
      DO 101 M=1,12
  101 READ (IFILE) XTITLE,PREDD(:,:,:,M,IDD)
      call closeunit (ifile)
!**** Industrial Sulfates
      call openunit (RDFILE(2),ifile,.true.,.true.)
      DO 102 IDD=1,8
      DO 102 M=1,12
  102 READ (IFILE) XTITLE,SUIDD(:,:,:,M,IDD)
      call closeunit (ifile)
!**** Industrial Organic Carbons
      call openunit (RDFILE(3),ifile,.true.,.true.)
      DO 103 IDD=1,8
      DO 103 M=1,12
  103 READ (IFILE) XTITLE,OCIDD(:,:,:,M,IDD)
      call closeunit (ifile)
!**** Industrial Black Carbons
      call openunit (RDFILE(4),ifile,.true.,.true.)
      DO 104 IDD=1,8
      DO 104 M=1,12
  104 READ (IFILE) XTITLE,BCIDD(:,:,:,M,IDD)
      call closeunit (ifile)

C**** Prepare for aerosol indirect effect parameterization:
C     - Collect the monthly aerosol number densities (an) for the time
C       independent aerosols (desert dust and sea salt)       an:  /cm^3
C     - Save the monthly 1850 mass densities (md) for the time dependent
C       aerosols (Sulfates,Nitrates,Organic & Black Carbons)  md: kg/cm3

!!!   call openunit (RDFILE(5),ifile,.true.,.true.) !neglect desert dust
!!!   xdust=.33/(2000.*4.1888*(.40d-6)**3)     ! f/[rho*4pi/3*r^3] (/kg)
      xsslt=aermix(3)/(2000.*4.1888*(.44d-6)**3) ! x/particle-mass (/kg)
      byz = 1d-6/za720 ! 1d-6/depth in m (+conversion /m3 -> /cm3)
      DO M=1,12
!!!     READ (IFILE) XTITLE,mddust
      DO J=1,46
      DO I=1,72
        anfix(i,j,m) = 0. !!! xdust*mddust(i,j) ! aerosol number (/cm^3)
     +               +    byz * SUM(PREDD(I,J,1:la720,M,3)) * Xsslt
C****   md1850(1:4,i,j,m)  !  mass density (kg/cm^3): SO4, NO3, OC, BCB
        md1850(1,i,j,m) = byz * SUM(AERMIX(1)*PREDD(I,J,1:La720,M,1) +
     +                              AERMIX(2)*PREDD(I,J,1:La720,M,2))
        md1850(2,i,j,m) = byz * SUM(AERMIX(4)*PREDD(I,J,1:La720,M,4))
        md1850(3,i,j,m) = byz * SUM(AERMIX(5)*PREDD(I,J,1:La720,M,5) +
     +                              AERMIX(6)*PREDD(I,J,1:La720,M,6))
        md1850(4,i,j,m) = byz * SUM(AERMIX(7)*PREDD(I,J,1:La720,M,7))
      end do
      end do
      end do
      anfix(:,:,0) = anfix(:,:,12) ; md1850(:,:,:,0) = md1850(:,:,:,12)
!!!   call closeunit (ifile)

      IFIRST=0
      ENDIF


C     To time input data READs, JYEARX is set ahead of JYEARA by 15 days
C     ------------------------------------------------------------------
      if(JYEARA<0) then
        JYEARX = -JYEARA
      else
        JYEARX=MIN(JYEARA+(JJDAYA+15)/366,2050)
      end if

      IF(JYEARX==JYRNOW) GO TO 500    ! Get A6JDAY from current A6YEAR

C     Begin current A6YEAR  with 1850 Background SO4,SEA,ANT,OCX,BCI,BCB
      DO 114 M=1,12
      A6YEAR(:,:,:,M,1) = AERMIX(1)*PREDD(:,:,:,M,1)*1000*DRYM2G(1)
     +                   +AERMIX(2)*PREDD(:,:,:,M,2)*1000*DRYM2G(1)
      A6YEAR(:,:,:,M,2) = AERMIX(3)*PREDD(:,:,:,M,3)*1000*DRYM2G(2)
      A6YEAR(:,:,:,M,3) = AERMIX(4)*PREDD(:,:,:,M,4)*1000*DRYM2G(3)
      A6YEAR(:,:,:,M,4) = AERMIX(5)*PREDD(:,:,:,M,5)*1000*DRYM2G(4)
     +                   +AERMIX(6)*PREDD(:,:,:,M,6)*1000*DRYM2G(4)
      A6YEAR(:,:,:,M,5) = 0
      A6YEAR(:,:,:,M,6) = AERMIX(7)*PREDD(:,:,:,M,7)*1000*DRYM2G(6)
  114 CONTINUE
!****                                   Define 1849 Background  Dec data
      DO N=1,6
        A6YEAR(:,:,:,0,N)=A6YEAR(:,:,:,12,N)
      END DO

      IF(JYEARX > 1850) THEN                           !   (JYEAR>1850)
        WTANI=GLOPOP(JYEARX)
        WTOCB=min( 1d0 , (JYEARX-1850)/140.D0 )
        WTBCB=min( 1d0 , (JYEARX-1850)/140.D0 )
        DO M=1,12            !  Add time dependent JYEAR ANI,OCB,BCB
          A6YEAR(:,:,:,M,3) = A6YEAR(:,:,:,M,3)+
     +        AERMIX( 9)*WTANI*PREDD(:,:,:,M, 8)*1000*DRYM2G(3)
          A6YEAR(:,:,:,M,4) = A6YEAR(:,:,:,M,4)+
     +        AERMIX(12)*WTOCB*PREDD(:,:,:,M, 9)*1000*DRYM2G(4)
          A6YEAR(:,:,:,M,6) = A6YEAR(:,:,:,M,6)+
     +        AERMIX(13)*WTBCB*PREDD(:,:,:,M,10)*1000*DRYM2G(6)
        END DO
        WTANI=GLOPOP(JYEARX-1)
        WTOCB=min( 139/140d0 , (JYEARX-1851)/140.D0 )
        WTBCB=min( 139/140d0 , (JYEARX-1851)/140.D0 )
        M=12        !  Add time dependent JYEAR-1 ANI,OCB,BCB Dec data
        A6YEAR(:,:,:,0,3) = A6YEAR(:,:,:,0,3)+
     +      AERMIX( 9)*WTANI*PREDD(:,:,:,M, 8)*1000*DRYM2G(3)
        A6YEAR(:,:,:,0,4) = A6YEAR(:,:,:,0,4)+
     +      AERMIX(12)*WTOCB*PREDD(:,:,:,M, 9)*1000*DRYM2G(4)
        A6YEAR(:,:,:,0,6) = A6YEAR(:,:,:,0,6)+
     +      AERMIX(13)*WTBCB*PREDD(:,:,:,M,10)*1000*DRYM2G(6)
      ENDIF

      IF(JYEARX > 1850.and.JYEARX < 1876) THEN   !   (1850<JYEAR<1876)
       WT75=(JYEARX-1850)/25.D0
       DO M=1,12          !    Add time dependent JYEAR SUI,OCI,BCI
       A6YEAR(:,:,:,M,1) = A6YEAR(:,:,:,M,1)+
     +                 WT75*SUIDD(:,:,:,M,1)*AERMIX( 8)*1000*DRYM2G(1)
       A6YEAR(:,:,:,M,4) = A6YEAR(:,:,:,M,4)+
     +                 WT75*OCIDD(:,:,:,M,1)*AERMIX(10)*1000*DRYM2G(4)
       A6YEAR(:,:,:,M,5) = A6YEAR(:,:,:,M,5)+
     +                 WT75*BCIDD(:,:,:,M,1)*AERMIX(11)*1000*DRYM2G(5)
       END DO

       WT75=(JYEARX-1851)/25.D0
       M=12          !  Add time dependent JYEAR-1 SUI,OCI,BCI Dec data
       A6YEAR(:,:,:,0,1) = A6YEAR(:,:,:,0,1)+
     +                 WT75*SUIDD(:,:,:,M,1)*AERMIX( 8)*1000*DRYM2G(1)
       A6YEAR(:,:,:,0,4) = A6YEAR(:,:,:,0,4)+
     +                 WT75*OCIDD(:,:,:,M,1)*AERMIX(10)*1000*DRYM2G(4)
       A6YEAR(:,:,:,0,5) = A6YEAR(:,:,:,0,5)+
     +                 WT75*BCIDD(:,:,:,M,1)*AERMIX(11)*1000*DRYM2G(5)
      ENDIF

      IF(JYEARX > 1875) THEN                         !     (JYEAR>1875)
      CALL STREND(JYEARX,IYS,JYS,SWTI,SWTJ)
      CALL CTREND(JYEARX,IYC,JYC,CWTI,CWTJ)
      DO 141 M=1,12            !    Add time dependent JYEAR SUI,OCI,BCI
      DO 141 L=1,9
      DO 141 J=1,46
      DO 141 I=1,72
      A6YEAR(I,J,L,M,1)=A6YEAR(I,J,L,M,1)+AERMIX( 8)*1000.D0*DRYM2G(1)*
     +                 (SWTI*SUIDD(I,J,L,M,IYS)+SWTJ*SUIDD(I,J,L,M,JYS))
      A6YEAR(I,J,L,M,4)=A6YEAR(I,J,L,M,4)+AERMIX(10)*1000.D0*DRYM2G(4)*
     +                 (CWTI*OCIDD(I,J,L,M,IYC)+CWTJ*OCIDD(I,J,L,M,JYC))
      A6YEAR(I,J,L,M,5)=A6YEAR(I,J,L,M,5)+AERMIX(11)*1000.D0*DRYM2G(5)*
     +                 (CWTI*BCIDD(I,J,L,M,IYC)+CWTJ*BCIDD(I,J,L,M,JYC))
      IF(A6YEAR(I,J,L,M,1) < 0.) A6YEAR(I,J,L,M,1)=0.
      IF(A6YEAR(I,J,L,M,4) < 0.) A6YEAR(I,J,L,M,4)=0.
      IF(A6YEAR(I,J,L,M,5) < 0.) A6YEAR(I,J,L,M,5)=0.
  141 CONTINUE

      CALL STREND(JYEARX-1,IYS,JYS,SWTI,SWTJ)
      CALL CTREND(JYEARX-1,IYC,JYC,CWTI,CWTJ)
      M=12            !  Add time dependent JYEAR-1 SUI,OCI,BCI Dec data
      DO 145 L=1,9
      DO 145 J=1,46
      DO 145 I=1,72
      A6YEAR(I,J,L,0,1)=A6YEAR(I,J,L,0,1)+AERMIX( 8)*1000.D0*DRYM2G(1)*
     +                 (SWTI*SUIDD(I,J,L,M,IYS)+SWTJ*SUIDD(I,J,L,M,JYS))
      A6YEAR(I,J,L,0,4)=A6YEAR(I,J,L,0,4)+AERMIX(10)*1000.D0*DRYM2G(4)*
     +                 (CWTI*OCIDD(I,J,L,M,IYC)+CWTJ*OCIDD(I,J,L,M,JYC))
      A6YEAR(I,J,L,0,5)=A6YEAR(I,J,L,0,5)+AERMIX(11)*1000.D0*DRYM2G(5)*
     +                 (CWTI*BCIDD(I,J,L,M,IYC)+CWTJ*BCIDD(I,J,L,M,JYC))
      IF(A6YEAR(I,J,L,0,1) < 0.) A6YEAR(I,J,L,0,1)=0.
      IF(A6YEAR(I,J,L,0,4) < 0.) A6YEAR(I,J,L,0,4)=0.
      IF(A6YEAR(I,J,L,0,5) < 0.) A6YEAR(I,J,L,0,5)=0.
  145 CONTINUE
      ENDIF
      JYRNOW=JYEARX
      if(jyeara<0) then  ! cyclic case
        DO N=1,6
          A6YEAR(:,:,:,0,N)=A6YEAR(:,:,:,12,N)
        END DO
      end if

C      A6JDAY is interpolated daily from A6YEAR seasonal data via JJDAYA
C      -----------------------------------------------------------------

  500 CONTINUE
      XMI=(JJDAYA+JJDAYA+31-(JJDAYA+15)/61+(JJDAYA+14)/61)/61.D0
      MI=XMI
      WTMJ=XMI-MI       !   Intra-year interpolation is linear in JJDAYA
      WTMI=1.D0-WTMJ
      IF(MI > 11) MI=0
      MJ=MI+1
      DO 510 J=1,46
      DO 510 I=1,72
      DO 510 N=1,6
      DO 510 L=1,9
      A6JDAY(L,N,I,J)=WTMI*A6YEAR(I,J,L,MI,N)+WTMJ*A6YEAR(I,J,L,MJ,N)
  510 CONTINUE

C**** Needed for aerosol indirect effect parameterization in GCM
      byz=1d-9/za720
      do j=1,46
      do i=1,72
C**** sea salt, desert dust
         anssdd(i,j) = WTMI*anfix(i,j,mi)+WTMJ*anfix(i,j,mj)
C**** SU4,NO3,OCX,BCB,BCI (reordered: no sea salt, no pre-ind BCI)
         mdpi(:,i,j) =
     &        WTMI*md1850(:,i,j,mi) + WTMJ*md1850(:,i,j,mj) !1:4
         mdcur(1,i,j) = SUM(A6JDAY(1:La720,1,I,J)) * byz/drym2g(1)
         mdcur(2,i,j) = SUM(A6JDAY(1:La720,3,I,J)) * byz/drym2g(3)
         mdcur(3,i,j) = SUM(A6JDAY(1:La720,4,I,J)) * byz/drym2g(4)
         mdcur(4,i,j) = SUM(A6JDAY(1:La720,6,I,J)) * byz/drym2g(6)
         mdcur(5,i,j) = SUM(A6JDAY(1:La720,5,I,J)) * byz/drym2g(5)
      end do
      end do

      RETURN        !  A6JDAY(9,6,72,46) is used in GETAER via ILON,JLAT
      END SUBROUTINE updateAerosol

#ifndef NEW_IO_AERINP
      subroutine updateAerosol2(jYearA, jjDaya, a6jday, plbaer)

C     ------------------------------------------------------------------
C     Reads: XXX_Koch2008_kg_m2_72x46x20_1880-2000 aerosol kg/m2 data
C     for SUL,NIT,OCA,BCA,BCB and
C            SSA_Koch2008_kg_m2_72x46x20
C
C    Makes: A6YEAR2(72,46,12,0:12,6),A6JDAY(20,6,72,46) (dry aerosol Tau)
C     ------------------------------------------------------------------

      USE FILEMANAGER, only : openunit,closeunit
      use domain_decomp_atm, only : grid
      implicit none

      INTEGER, intent(in) :: jyeara,jjdaya
      real*8, pointer :: a6jday(:,:,:,:)
      real*8, dimension(:), pointer ::  plbaer

      character*80 :: aertitle ! aerosol info
      integer, save :: ndeca, fdeca, ldeca ! dimensions needed for aerosols

      CHARACTER*40 :: RDFILE(7) = (/                !  Input file names
     1            'SUL_Koch2008_kg_m2_72x46x20_1890-2000h  '
     2           ,'SSA_Koch2008_kg_m2_72x46x20h            '
     3           ,'NIT_Bauer2008_kg_m2_72x46x20_1890-2000h '
     4           ,'OCA_Koch2008_kg_m2_72x46x20_1890-2000h  '
     5           ,'BCA_Koch2008_kg_m2_72x46x20_1890-2000h  '
     6           ,'BCB_Koch2008_kg_m2_72x46x20_1890-2000h  '
     7           ,'low.dust.72x46.monthly.bin              '/)

      CHARACTER*40 :: RDFGEN(7) = (/                ! generic names
     * 'TAero_SUL','TAero_SSA','TAero_NIT','TAero_OCA','TAero_BCA',
     * 'TAero_BCB','M_LowDust'/)

C                TROPOSPHERIC AEROSOL COMPOSITIONAL/TYPE PARAMETERS
C                   SO4    SEA    ANT    OCX    BCI    BCB   *BCB  *BCB
C     DATA REFDRY/0.200, 1.000, 0.300, 0.300, 0.100, 0.100, 0.200,0.050/
C
C     DATA REFWET/0.272, 1.808, 0.398, 0.318, 0.100, 0.100, 0.200,0.050/
C
C     DATA DRYM2G/4.667, 0.866, 4.448, 5.018, 9.000, 9.000, 5.521,8.169/
C
CKoch DATA DRYM2G/5.000, 2.866, 8.000, 8.000, 9.000, 9.000, 5.521,8.169/
C
C     DATA RHTMAG/1.788, 3.310, 1.756, 1.163, 1.000, 1.000, 1.000,1.000/
C
CRH70 DATA WETM2G/8.345, 2.866, 7.811, 5.836, 9.000, 9.000, 5.521,8.169/
C
C     DATA Q55DRY/2.191, 2.499, 3.069, 3.010, 1.560, 1.560, 1.914,0.708/
C
C     DATA DENAER/1.760, 2.165, 1.725, 1.500, 1.300, 1.300, 1.300,1.300/
C
C     ------------------------------------------------------------------
C          DRYM2G(I) = 0.75/DENAER(I)*Q55DRY(I)/REFDRY(I)
C          WETM2G(I) = DRYM2G(I)*RHTMAG(I)
C          RHTMAG(I) = Rel Humidity TAU Magnification factor  at RH=0.70
C          REFWET(I) = Rel Humidity REFDRY Magnification      at RH=0.70
C     ------------------------------------------------------------------
      integer :: ifile
      INTEGER, save :: JYRNOW=0

      INTEGER ia,idd,ndd,m,mi,mj,i,j,l,n,jyearx,iy,jy,iyc,jyc,iyp,jyp
      REAL*8 wti,wtj,cwti,cwtj,pwti,pwtj,xmi,wtmi,wtmj
      REAL*8 xsslt ! ,xdust
      real*8 :: dp,mindp

      logical, save :: init = .false.
      integer :: idxMonth0, idxMonth1
      integer :: idxDecade0, idxDecade1
      logical :: updateTable

      integer, save :: decadeNow = -9999
      integer, save :: yearNow = -9999
      integer, save :: monthNow = -1

      integer :: im, id, mm

      integer :: i_0,i_1,j_0,j_1

      i_0 = grid%i_strt
      i_1 = grid%i_stop
      j_0 = grid%j_strt
      j_1 = grid%j_stop
      
      if (.not. init) then
        init = .true.
        RDFILE=RDFGEN           !     generic names are used
c read table sizes then close
        call openunit (RDFILE(1),ifile,.true.,.true.) ! unformatted,old
        read(ifile) aertitle, ima, jma, lma, ndeca, fdeca, ldeca

!!    check whether the newer input files are being used
        if(aertitle(1:26)=='                          ') call
     *       stop_model('updateAerosol2: use newer input files',255)

        if (.not. associated(plbaer)) allocate( plbaer(lma+1) )

        allocate(
     &        suldd(i_0:i_1,j_0:j_1,lma,2,2)
     &       ,nitdd(i_0:i_1,j_0:j_1,lma,2,2)
     &       ,ocadd(i_0:i_1,j_0:j_1,lma,2,2)
     &       ,bcadd(i_0:i_1,j_0:j_1,lma,2,2)
     &       ,bcbdd(i_0:i_1,j_0:j_1,lma,2,2)
     &       ,ssadd(i_0:i_1,j_0:j_1,lma,2)
     &       )
        if (.not.associated(A6JDAY))
     *       allocate(A6JDAY(lma,6,i_0:i_1,j_0:j_1))
        allocate( A6YEAR2 (i_0:i_1,j_0:j_1,lma,2,6) )
        allocate( md1850 (4,i_0:i_1,j_0:j_1,0:12) )
        allocate( anfix (i_0:i_1,j_0:j_1,2) )

        allocate(anssdd(i_0:i_1,j_0:j_1))
        allocate(mdpi(4,i_0:i_1,j_0:j_1))
        allocate(mdcur(5,i_0:i_1,j_0:j_1))

        read(ifile) plbaer
        call closeUnit(ifile)

!**** For parameterized AIE, find level whose top is closest to 720 mb
        mindp = 1d30
        do La720=1,lma
          dp = abs(plbaer(La720)-720d0)
          if(dp > mindp) exit
          mindp = dp
        enddo
        La720 = La720 - 2

!**** Pre-industrial mass densities for parameterized AIE
        do m = 1, 12
          call readTable(RDFILE(1), SULDD(:,:,:,1,1), month=m,decade=1)
          md1850(1,:,:,m) = byz_cm3 * SUM(SULDD(:,:,1:La720,1,1), DIM=3)
          call readTable(RDFILE(3), NITDD(:,:,:,1,1), month=m,decade=1)
          md1850(2,:,:,m) = byz_cm3 * SUM(NITDD(:,:,1:La720,1,1), DIM=3)
          call readTable(RDFILE(4), OCADD(:,:,:,1,1), month=m,decade=1)
          md1850(3,:,:,m) = byz_cm3 * SUM(OCADD(:,:,1:La720,1,1), DIM=3)
          call readTable(RDFILE(5), BCBDD(:,:,:,1,1), month=m,decade=1)
          call readTable(RDFILE(6), BCADD(:,:,:,1,1), month=m,decade=1)
          md1850(4,:,:,m) = byz_cm3 * (
     &         SUM(BCBDD(:,:,1:La720,1,1), DIM=3) +
     &         SUM(BCADD(:,:,1:La720,1,1), DIM=3) )
        end do
        md1850(:,:,:,0) = md1850(:,:,:,12)

      end if

      call getRefMonth(jjdaya, idxMonth0)
      call getRefDecade(jjdaya, jyeara, idxMonth0, ndeca, fdeca, ldeca,
     &     idxDecade0)

      updateTable = (idxMonth0 /= monthNow .or.
     &     idxDecade0 /= decadeNow )

      if (updateTable) then
        monthNow = idxMonth0
        decadeNow = idxDecade0

        do im = 1, 2
          select case (im)
          case (1)
            mm = idxMonth0
            m = mm
            if (m == 0) m = 12
            call getRefDecade(jjdaya, jyeara, idxMonth0,
     &         ndeca, fdeca, ldeca,
     &         idxDecade0)
          case (2)
            mm = idxMonth0 + 1
            m = mm
            call getRefDecade(jjdaya, jyeara, idxMonth0+1,
     &         ndeca, fdeca, ldeca,
     &         idxDecade0)
          end select

!**** Sea salt
          call readTable(RDFILE(2), SSADD(:,:,:,im), month=m, decade=1)

          do id = 1, 2
            select case (id)
            case (1)
              idd = idxDecade0
            case (2)
              idd = min(idxDecade0+1, ndeca)
            end select

!**** Sulfate
            call readTable(RDFILE(1),SULDD(:,:,:,im,id),
     &           month=m, decade=idd)
!**** Nitrate
            call readTable(RDFILE(3),NITDD(:,:,:,im,id),
     &           month=m, decade=idd)
!**** Organic Carbon
            call readTable(RDFILE(4),OCADD(:,:,:,im,id),
     &           month=m, decade=idd)
!**** Black Carbon from Fossil and bio fuel
            call readTable(RDFILE(5),BCADD(:,:,:,im,id),
     &           month=m,decade=idd)
!**** Black Carbon from Biomass burning
            call readTable(RDFILE(6),BCBDD(:,:,:,im,id),
     &           month=m,decade=idd)
          end do

C**** Prepare for aerosol indirect effect parameterization:
C     - Collect the monthly aerosol number densities (an) for the time
C       independent aerosols (desert dust and sea salt)       an:  /cm^3
C     - Save the monthly 1850 mass densities (md) for the time dependent
C       aerosols (Sulfates,Nitrates,Organic & Black Carbons)  md: kg/cm3

!!!   xdust=.33/(2000.*4.1888*(.40d-6)**3)     ! f/[rho*4pi/3*r^3] (/kg)
        xsslt=1.d0/(2000.*4.1888*(.44d-6)**3) ! x/particle-mass (/kg)

        do J=J_0,J_1
          do I=I_0,I_1
c SUM to L=5 for low clouds only
c Using 1890 not 1850 values here
            anfix(i,j,im) = 0.   !!! xdust*mddust(i,j) ! aerosol number (/cm^3)
     +           +    byz_cm3 * SUM(SSADD(I,J,1:La720,im)) * Xsslt
          end do
        end do

      end do

      end if

C                                                   0            12
C     Collect 13 months of data for time period Dec/15/(yr-1)-Dec/15/yr:
C     To time input data READs, JYEARX is set ahead of JYEARA by 15 days
C     ------------------------------------------------------------------
      if(JYEARA<0) then
        JYEARX = -JYEARA
      else
!TODO - hardwired for 365 days?   Year 2050?
        !JYEARX=MIN(JYEARA+(JJDAYA+15)/366,2050)
        JYEARX=JYEARA+(JJDAYA+15)/366
      end if

      DO iM=1,2
      A6YEAR2(:,:,:,im,2) = SSADD(:,:,:,im)*1000*DRYM2G(2) !*AERMIX(3)
      END DO
!****

      IF(JYEARX < fdeca .or. JYEARX > ldeca) THEN   !   (use first/last decadal mean)
        DO iM=1,2               !    Set time dependent JYEAR SU,NIT,OC,BC
          A6YEAR2(:,:,:,im,1) = SULDD(:,:,:,im,1)*1000*DRYM2G(1) !*AERMIX( 8)
          A6YEAR2(:,:,:,im,3) = NITDD(:,:,:,im,1)*1000*DRYM2G(3) !*AERMIX( 4)
          A6YEAR2(:,:,:,im,4) = OCADD(:,:,:,im,1)*1000*DRYM2G(4) !*AERMIX(10)
          A6YEAR2(:,:,:,im,5) = BCADD(:,:,:,im,1)*1000*DRYM2G(5) !*AERMIX(11)
          A6YEAR2(:,:,:,im,6) = BCBDD(:,:,:,im,1)*1000*DRYM2G(6) !*AERMIX(13)
        END DO

      ELSE  !  IF(JYEARX.ge.fdeca.and.JYEARX.LE.ldeca) THEN
        iyc=INT((JYEARX-fdeca)/10.d0)+1 ! current year
        jyc=min(ndeca,iyc+1)    ! have only ndeca decades
        cwtj=(JYEARX-fdeca)/10.d0-INT((JYEARX-fdeca)/10.d0)
        cwti=1.d0-cwtj
        iyp=INT((JYEARX-(fdeca+1))/10.d0)+1 ! previous year
        jyp=iyp+1
        pwtj=(JYEARX-(fdeca+1))/10.d0-INT((JYEARX-(fdeca+1))/10.d0)
        if(JYEARX==fdeca) then
          iyp=1 ; jyp=1 ; pwtj=0
        end if
        pwti=1.d0-pwtj
        DO 141 im=1,2           !    Set time dependent JYEAR SU,NIT,OC,BC
          wti=cwti ; if (idxMonth0 == 0 .and. im == 1) wti=pwti
          wtj=cwtj ; if (idxMonth0 == 0 .and. im == 1) wtj=pwtj

          iy = 1
          jy = 2

          DO 141 L=1,lma        !   AERMIX scalings are removed
            DO 141 J=J_0,J_1
              DO 141 I=I_0,I_1
                A6YEAR2(I,J,L,im,1) = 1000.D0*DRYM2G(1)* ! AERMIX( 8)*
     +               (WTI*SULDD(I,J,L,im,IY)+WTJ*SULDD(I,J,L,im,JY))
                A6YEAR2(I,J,L,im,3) = 1000.D0*DRYM2G(3)* ! AERMIX( 4)*
     +               (WTI*NITDD(I,J,L,im,IY)+WTJ*NITDD(I,J,L,im,JY))
                A6YEAR2(I,J,L,im,4) = 1000.D0*DRYM2G(4)* ! AERMIX(10)*
     +               (WTI*OCADD(I,J,L,im,IY)+WTJ*OCADD(I,J,L,im,JY))
                A6YEAR2(I,J,L,im,5) = 1000.D0*DRYM2G(5)* ! AERMIX(11)*
     +               (WTI*BCADD(I,J,L,im,IY)+WTJ*BCADD(I,J,L,im,JY))
                A6YEAR2(I,J,L,im,6) = 1000.D0*DRYM2G(6)* ! AERMIX(13)*
     +               (WTI*BCBDD(I,J,L,im,IY)+WTJ*BCBDD(I,J,L,im,JY))
 141          CONTINUE
      ENDIF
      JYRNOW=JYEARX

C      A6JDAY is interpolated daily from A6YEAR2 seasonal data via JJDAYA
C      -----------------------------------------------------------------

  500 CONTINUE
      XMI=(JJDAYA+JJDAYA+31-(JJDAYA+15)/61+(JJDAYA+14)/61)/61.D0
      MI=XMI
      WTMJ=XMI-MI       !   Intra-year interpolation is linear in JJDAYA
      WTMI=1.D0-WTMJ
      IF(MI > 11) MI=0
      MJ=MI+1

      DO 510 J=J_0,J_1
      DO 510 I=I_0,I_1
      DO 510 N=1,6
      DO 510 L=1,lma
      A6JDAY(L,N,I,J)=WTMI*A6YEAR2(I,J,L,1,N)+WTMJ*A6YEAR2(I,J,L,2,N)
  510 CONTINUE

C**** Needed for aerosol indirect effect parameterization in GCM
      do j=J_0,J_1
      do i=I_0,I_1
C**** sea salt, desert dust
        anssdd(i,j) = WTMI*anfix(i,j,1)+WTMJ*anfix(i,j,2)
C**** SU4,NO3,OCX,BCB,BCI (reordered: no sea salt, no pre-ind BCI)
        mdpi(:,i,j) =
     &       WTMI*md1850(:,i,j,mi) + WTMJ*md1850(:,i,j,mj) !1:4
        mdcur(1,i,j) = SUM (A6JDAY(1:La720,1,I,J))*
     &       byz_gcm3/drym2g(1)
        mdcur(2,i,j) = SUM (A6JDAY(1:La720,3,I,J))*
     &       byz_gcm3/drym2g(3)
        mdcur(3,i,j) = SUM (A6JDAY(1:La720,4,I,J))*
     &       byz_gcm3/drym2g(4)
        mdcur(4,i,j) = SUM (A6JDAY(1:La720,6,I,J))*
     &       byz_gcm3/drym2g(6)
        mdcur(5,i,j) = SUM (A6JDAY(1:La720,5,I,J))*
     &       byz_gcm3/drym2g(5)
      end do
      end do

      return        !  A6JDAY(9,6,72,46) is used in GETAER via ILON,JLAT
      end subroutine updateAerosol2

      subroutine getRefMonth(day, idxMonth0)
      integer, intent(in) :: day
      integer, intent(out) :: idxMonth0

      idxMonth0 = (day + day + 31 - (day+15)/61 + (day+14)/61) / 61.d+0
      if (idxMonth0 > 11) idxMonth0 = 0

      end subroutine getRefMonth

      subroutine getRefDecade(day, year, month, ndeca, fdeca, ldeca,
     &     idxDecade0)
      integer, intent(in) :: day
      integer, intent(in) :: year
      integer, intent(in) :: month
      integer, intent(in) :: ndeca, fdeca, ldeca
      integer, intent(out) :: idxDecade0

      integer :: refYear
      integer :: iyc, jyc, iyp, jyp

      if (year < 0) then
        refYear = -year
      else
        !refYear = min(year + (day+15)/366,2050)
        refYear = year + (day+15)/366
      end if

      if (refYear < fdeca) then
        idxDecade0 = 1
      else if (refYear > ldeca) then
        idxDecade0 = ndeca
      else
        iyc = int((refYear - fdeca)/10.d0) + 1
        jyc = min(ndeca,idxDecade0 + 1)

        iyp = int((refYear - (fdeca+1))/10.d0)+1
        jyp = iyp + 1

        if (month == 0) then
          idxDecade0 = iyp
        else
          idxDecade0 = iyc
        end if

      end if

      end subroutine getRefDecade

      subroutine readTable(fileName, array, month, decade)
      use FileManager, only: openUnit, closeUnit
      use domain_decomp_atm, only : grid,dread_parallel
      character(len=*), intent(in) :: fileName
      real*4, intent(inout) :: array(grid%i_strt:,grid%j_strt:,:)
      integer, intent(in) :: month
      integer, intent(in) :: decade

      character*80 :: title
      integer :: ifile
      integer :: ima, jma, lma, ndeca, fdeca, ldeca
      integer :: m, nskip
      integer :: i,j,l
      integer :: i_0,i_1,j_0,j_1

      real*4, allocatable :: array_glob(:,:,:)
      real*8, allocatable :: r8array(:,:,:)
      real*8, dimension(:,:), allocatable :: r8glob,r8loc

      i_0 = grid%i_strt
      i_1 = grid%i_stop
      j_0 = grid%j_strt
      j_1 = grid%j_stop

      call openunit (fileName, ifile, .true., .true.)    ! unformatted,old
      read(ifile) title, ima, jma, lma, ndeca, fdeca, ldeca
      read(ifile) ! skip plbaer

      ! skip   earlier decades + earlier months in this decade
      nskip =    12*(decade-1) +           (month-1)
      do m = 1,nskip
        read(ifile)
      enddo

      ! read actual month
      if(grid%im_world == ima .and. grid%jm_world == jma) then
        allocate(r8array(grid%i_strt_halo:grid%i_stop_halo,
     &                   grid%j_strt_halo:grid%j_stop_halo,lma))
        call dread_parallel(grid,ifile,'AEROFILE',r8array)
        array = r8array(grid%i_strt:grid%i_stop,
     &                  grid%j_strt:grid%j_stop,:)
        deallocate(r8array)
      else
        allocate(array_glob(ima,jma,lma))
        allocate(r8glob(ima,jma))
        allocate(r8loc(grid%i_strt_halo:grid%i_stop_halo,
     &                 grid%j_strt_halo:grid%j_stop_halo))
        read (ifile) array_glob
        do l=1,lma
          r8glob = array_glob(:,:,l)
          call sample_latlon(ima,jma,r8glob,r8loc)
          array(i_0:i_1,j_0:j_1,l) = r8loc(i_0:i_1,j_0:j_1)
        enddo
        deallocate(array_glob,r8glob,r8loc)
      endif

      call closeunit (ifile)

      end subroutine readTable

#else /* NEW_IO_AERINP */

      subroutine updateAerosol2(jYearA, jjDaya, a6jday, plbaer)
!@sum updateAerosol2 reads aerosol file(s) and calculates A6JDAY(lma,6,:,:)
!@+   (dry aerosol Tau) for current day, year.  On startup, it allocates
!@+   a6jday and plbaer, and reads plbaer.  Note that jYearA may be
!@+   negative, which is the Model E method for indicating that the
!@+   data for abs(jYearA) is to be used for all years (this matters
!@+   when time-interpolating between December and January).
! This is the netcdf version.  For additional information, see the
! Fortran-IO version above.
! This version gives different results than the Fortran-IO version, since
! the latter calculates the inter-month time interpolation weight using
!   XMI=(JJDAYA+JJDAYA+31-(JJDAYA+15)/61+(JJDAYA+14)/61)/61.D0
! and read_stream calculates it based on the midpoint dates of the
! current and following months.
! As currently coded, the time-varying aerosols are all assumed to
! reside in a single input file TAERO.  Pass a per-aerosol path to
! init_stream to change this.  Multiple aerosol input files would
! permit multiple time axes/resolutions for example.  Note that
! the radiation code still assumes that all aerosols are on the same
! vertical grid.  A per-aerosol plbaer in calls to REPART in the
! radiation code would allow more flexibility.
      use domain_decomp_atm, only : grid
      use model_com, only : jdmidofm ! for md1850 month interp
      use timestream_mod, only : init_stream,read_stream
     &     ,reset_stream_properties,get_by_index
      use pario, only : par_open,par_close,read_dist_data,read_data
     &     ,get_dimlen,get_dimlens
      use param
      implicit none
c
!@var jyeara current year; negative if the same year is to be repeated
!@var jjdaya current day
      INTEGER, intent(in) :: jyeara,jjdaya
!@var a6jday optical depth for 6 aerosol types
      real*8, pointer :: a6jday(:,:,:,:)
!@var plbaer pressures of layer interfaces in aerosol in datafiles
      real*8, dimension(:), pointer ::  plbaer
c
      INTEGER m,mi,mj,i,j,l,n,A6yrx(6)
      REAL*8 wtmi,wtmj
      REAL*8 xsslt ! ,xdust
      real*8 :: dp,mindp
      logical, save :: init = .false.
      logical :: cyclic
      integer :: dlens(7)
      character(len=16) :: A6yrstr
      integer :: i_0,i_1,j_0,j_1

      integer :: fid
      real*8, allocatable :: aerarr(:,:,:),arr12(:,:,:,:)

      i_0 = grid%i_strt
      i_1 = grid%i_stop
      j_0 = grid%j_strt
      j_1 = grid%j_stop

      if (.not. init) then
        init = .true.

        allocate( md1850 (4,i_0:i_1,j_0:j_1,0:12) )
        allocate(anssdd(i_0:i_1,j_0:j_1))
        allocate(mdpi(4,i_0:i_1,j_0:j_1))
        allocate(mdcur(5,i_0:i_1,j_0:j_1))

        fid = par_open(grid,'TAero_SSA','read')

        !lma = get_dimlen(grid,fid,'lev')
        call get_dimlens(grid,fid,'plbaer',n,dlens)
        lma = dlens(1)-1

        if (.not.associated(A6JDAY))
     *       allocate(A6JDAY(lma,6,i_0:i_1,j_0:j_1))

        if (.not. associated(plbaer)) allocate( plbaer(lma+1) )
        call read_data(grid,fid,'plbaer',plbaer,bcast_all=.true.)

        call par_close(grid,fid)

!**** For parameterized AIE, find level whose top is closest to 720 mb
        mindp = 1d30
        do La720=1,lma
          dp = abs(plbaer(La720)-720d0)
          if(dp > mindp) exit
          mindp = dp
        enddo
        La720 = La720 - 2

        allocate(arr12(grid%i_strt_halo:grid%i_stop_halo,
     &                 grid%j_strt_halo:grid%j_stop_halo,lma,12))
        allocate(aerarr(i_0:i_1,j_0:j_1,12))

        do n=1,6
          if(n.eq.2) cycle ! skip sea salt
          ! check whether this aerosol has a fixed year
          A6yrstr = trim(aernames(n))//'_yr'
          if(is_set_param(A6yrstr)) then
            call get_param(A6yrstr,A6yr(n))
            A6yr(n) = abs(A6yr(n)) ! just in case
            cyclic = .true.
          else
            A6yr(n) = -9999
            cyclic = jyeara < 0
          endif
          ! Initialize the stream to the year 1850 to extract
          ! the monthly climatology of 1850 aerosols.
          ! The read_stream call below will jump to the current year.
          call init_stream(grid,A6streams(n),
     &         'TAero_'//trim(aernames(n)),trim(aernames(n)),
     &         0d0,1d30,'linm2m',1850,1,cyclic=.true.)
          do m=1,12
            call get_by_index(grid,A6streams(n),m,arr12(:,:,:,m))
          enddo
          aerarr = byz_cm3 *
     &         SUM(arr12(i_0:i_1,j_0:j_1,1:La720,:), DIM=3)
          select case (n)
          case (1)
            md1850(1,:,:,1:12) = aerarr
          case (3,4,5)
            md1850(n-1,:,:,1:12) = aerarr
          case (6)
            md1850(4,:,:,1:12) = md1850(4,:,:,1:12) + aerarr
          end select
          ! Need this call to allow year jumps for cyclic case
          call reset_stream_properties(grid,A6streams(n),cyclic=cyclic)
        enddo

        ! No interannual variation of sea salt
        n = 2
        A6yr(n) = 1850 ! nominal year is always 1850
        call init_stream(grid,A6streams(n),
     &       'TAero_SSA',trim(aernames(n)),
     &       0d0,1d30,'linm2m',A6yr(n),1,cyclic=.true.)
        
        md1850(:,:,:,0) = md1850(:,:,:,12)
        deallocate(arr12,aerarr)

      endif ! end init

C**** read and time-interpolate
      where(A6yr.gt.0)
        A6yrx = A6yr
      elsewhere
        A6yrx = abs(jyeara)
      end where
      allocate(aerarr(grid%i_strt_halo:grid%i_stop_halo,
     &                grid%j_strt_halo:grid%j_stop_halo,lma))
      do n=1,6
        call read_stream(grid,A6streams(n),A6yrx(n),jjdaya,aerarr)
        do j=j_0,j_1
        do i=i_0,i_1
        do l=1,lma
          a6jday(l,n,i,j)=aerarr(i,j,l)
        enddo
        enddo
        enddo
      enddo
      deallocate(aerarr)

! remove AERMIX scalings
      DO J=J_0,J_1
      DO I=I_0,I_1
      DO N=1,6
      DO L=1,lma
        A6JDAY(L,N,I,J)=(1000.D0*DRYM2G(N))*A6JDAY(L,N,I,J)
      ENDDO
      ENDDO
      ENDDO
      ENDDO

C**** Calculate terms for aerosol indirect effect parameterization

      do i=1,13
        if(jjdaya.le.jdmidofm(i)) then
          wtmi = real(jdmidofm(i)-jjdaya,kind=8)/
     &               (jdmidofm(i)-jdmidofm(i-1))
          wtmj = 1d0-wtmi
          mi = i-1
          if(mi > 11) mi=0
          mj = mi+1
          exit
        endif
      enddo


!!!   xdust=.33/(2000.*4.1888*(.40d-6)**3)     ! f/[rho*4pi/3*r^3] (/kg)
      xsslt=1.d0/(2000.*4.1888*(.44d-6)**3) ! x/particle-mass (/kg)

      do j=J_0,J_1
      do i=I_0,I_1

C**** sea salt contribution to anssdd
c SUM to L=5 for low clouds only
        anssdd(i,j) = SUM(A6JDAY(1:La720,2,I,J))/(1000.D0*DRYM2G(2))
     &       *byz_cm3 * Xsslt

C**** SU4,NO3,OCX,BCB,BCI (reordered: no sea salt, no pre-ind BCI)
        mdpi(:,i,j) =
     &       WTMI*md1850(:,i,j,mi) + WTMJ*md1850(:,i,j,mj)
        mdcur(1,i,j) = SUM (A6JDAY(1:La720,1,I,J))*
     &       byz_gcm3/drym2g(1)
        mdcur(2,i,j) = SUM (A6JDAY(1:La720,3,I,J))*
     &       byz_gcm3/drym2g(3)
        mdcur(3,i,j) = SUM (A6JDAY(1:La720,4,I,J))*
     &       byz_gcm3/drym2g(4)
        mdcur(4,i,j) = SUM (A6JDAY(1:La720,6,I,J))*
     &       byz_gcm3/drym2g(6)
        mdcur(5,i,j) = SUM (A6JDAY(1:La720,5,I,J))*
     &       byz_gcm3/drym2g(5)
      end do
      end do

      return
      end subroutine updateAerosol2

#endif /* NEW_IO_AERINP */

      REAL*8 FUNCTION GLOPOP(JYEAR)
      IMPLICIT none

C     ----------------------------------------------------------------
C     GLOPOP = normalized global population trend set to unity in 1990
C              based on UN statistics & population projections to 2050
C
C     GLOPOP = 0.000 for 1850 and earlier
C            = 1.000 for 1990
C            = 1.658 for 2050 and later
C     ----------------------------------------------------------------

      INTEGER, intent(in) :: jyear
      REAL*8 :: GPNORM = 5.27-1.26 ,DNGPOP(21), GPOP(21) = (/
C               1850                     1900                     1950
     A          1.26,1.33,1.41,1.49,1.57,1.65,1.75,1.86,2.07,2.30,2.52
C                                        2000                     2050
     B              ,3.02,3.70,4.44,5.27,6.06,6.79,7.50,8.11,8.58,8.91/)
      INTEGER i,iy
      REAL*8 xy,dy

      DO 110 I=1,21
      DNGPOP(I)=(GPOP(I)-GPOP(1))/GPNORM
  110 CONTINUE
      XY=(JYEAR-1840)/10.D0
      IY=XY
      DY=XY-IY
      if (IY <  1) then ; IY =  1 ; DY = 0 ; endif
      if (IY > 20) then ; IY = 20 ; DY = 1 ; endif
      GLOPOP=DNGPOP(IY)+DY*(DNGPOP(IY+1)-DNGPOP(IY))
      RETURN
      END FUNCTION GLOPOP
!!!   GLOPOP and STREND,CTREND(in RAD_UTILS.f) are only needed by updateAerosol

      end module AerParam_mod

#ifndef NEW_IO_AERINP
      subroutine updBCd (year)
!@sum  reads appropriate Black Carbon deposition data if necessary
!@auth R. Ruedy
!@ver  1.0
      USE FILEMANAGER
      use AerParam_mod, only: depoBC,depoBC_1990
      USE DOMAIN_DECOMP_ATM, only: AM_I_ROOT,GRID,GET
      USE RESOLUTION, ONLY : im,jm
      implicit none
      integer, intent(in)   :: year

      integer :: imr,jmr
      integer :: i_0h,i_1h,j_0h,j_1h
      real*8  wt

      integer :: iu,year1,year2
      integer, save :: year0=0,yearL=-2,year_old=-1

      character*80 title
      real*8, dimension(:,:), allocatable :: BCdep,BCdep1,BCdep2
      real*4, dimension(:,:), allocatable :: BCdep4
      integer :: nrecs,reclens(1)

      if(.not.allocated(depoBC)) then
        CALL GET(grid, I_STRT_HALO=I_0H, I_STOP_HALO=I_1H,
     &                 J_STRT_HALO=J_0H, J_STOP_HALO=J_1H)
        ALLOCATE(depoBC     (I_0H:I_1H, J_0H:J_1H),
     &           depoBC_1990(I_0H:I_1H, J_0H:J_1H))
      endif

C**** check whether update is needed
      if (year.eq.year_old) return
      if (year_old.eq.yearL.and.year.gt.yearL) return
      if (year_old.eq.year0.and.year.lt.year0) return

      ! try to guess the resolution of the input file
      ! from the length of the first record, assuming
      ! 80-byte title followed by im*jm 4-byte values
      nrecs = 1 ! only checking the first record
      call get_recordlengths('BC_dep',nrecs,reclens)
      select case (reclens(1))
      case ( 80 + 4*72*46 )
        imr = 72
        jmr = 46
      case ( 80 + 4*144*90 )
        imr = 144
        jmr = 90
      case default ! try model native grid resolution
        if(reclens(1) .eq. 80 + 4*im*jm) then
          imr = im
          jmr = jm
        else
          call stop_model('unrecognized resolution in updBCd',255)
        endif
      end select

      allocate(BCdep(imr,jmr),
     &     BCdep1(imr,jmr),BCdep2(imr,jmr),BCdep4(imr,jmr))

      call openunit('BC_dep',iu,.true.,.true.)

      if (year_old.lt.0) then
C****   read whole input file and find range: year0->yearL
   10   read(iu,end=20) title
        read(title,*) yearL
        go to 10
      end if

   20 rewind (iu)
      read(iu) title,BCdep4
      read(title,*) year0
      BCdep1=BCdep4 ; BCdep2=BCdep4 ; year2=year0 ; year1=year0
      if (year.le.year1)              year2=year+1

      do while (year2.lt.year .and. year2.ne.yearL)
         year1 = year2 ; BCdep1 = BCdep2
         read (iu) title,BCdep4
         read(title,*) year2
         BCdep2 = BCdep4
      end do

      if(year.le.year1) then
        wt = 0.
      else if (year.ge.yearL) then
        wt = 1.
      else
        wt = (year-year1)/(real(year2-year1,kind=8))
      end if

      if (AM_I_ROOT())  write(6,*)
     &     'Using BCdep data from year',year1+wt*(year2-year1)
      call closeunit(iu)

C**** Set the Black Carbon deposition array
      BCdep(:,:) = BCdep1(:,:) + wt*(BCdep2(:,:)-BCdep1(:,:))

      call sample_latlon(imr,jmr,BCdep,depoBC)

      year_old = year

      deallocate(BCdep,BCdep1,BCdep2,BCdep4)

      return
      end subroutine updBCd
#else

      subroutine updBCd(year)
!@sum updBCd reads timeseries file for black carbon deposition
!@+   and interpolates depoBC to requested year.  This is the
!@+   netcdf version based on the traditional-I/O version above.
      use domain_decomp_atm, only : grid,get
      use timestream_mod, only : init_stream,read_stream
      use AerParam_mod, only: BCdepstream,depoBC,depoBC_1990
      implicit none
      integer, intent(in) :: year
c
      logical, save :: init = .false.
      integer :: i_0h,i_1h,j_0h,j_1h
      integer :: day

      day = 1 ! to pass a required argument

      if (.not. init) then
        init = .true.

        call get(grid, i_strt_halo=i_0h, i_stop_halo=i_1h,
     &                 j_strt_halo=j_0h, j_stop_halo=j_1h)
        allocate(depoBC     (i_0h:i_1h, j_0h:j_1h),
     &           depoBC_1990(i_0h:i_1h, j_0h:j_1h))

        call init_stream(grid,BCdepstream,'BC_dep','BC_dep',
     &       0d0,1d30,'none',year,day)
      endif

      call read_stream(grid,BCdepstream,year,day,depoBC)

      end subroutine updBCd

#endif

      module DustParam_mod
      ! persistent data remains in RADPAR.  read_alloc_dust lives
      ! in a module to allow pointers to be passed to it
      ! and allocated within.
      contains
#ifndef NEW_IO_AERINP
      subroutine read_alloc_dust(nrfu,imd,jmd,lmd,nsized,nmond,
     &     plbdust,rodust,redust,
     &     tdust
     &     )
!@sum This routine reads dust aerosol fields needed by the radiation
!@+   code in the prescribed-dust configuration of modelE.
!@auth R. Miller
!@auth M. Kelley extracted from RCOMP1, updated for grid flexibility,
!@+              added comments
      use DOMAIN_DECOMP_ATM, only: GRID,AM_I_ROOT,DREAD_PARALLEL
      use filemanager, only : nameunit,get_recordlengths
      implicit none
      integer :: nrfu,imd,jmd,lmd,nsized,nmond
      real*8, dimension(:), pointer :: plbdust, redust, rodust
      real*4, dimension(:,:,:,:,:), pointer :: tdust

      REAL*4, dimension(:,:,:,:,:), allocatable :: TDUST_glob ! just for read-in
      real*8, dimension(:,:,:), allocatable :: r8array ! just for read-in
      real*8, dimension(:,:), allocatable :: r8glob,r8loc ! just for read-in
      integer :: k3,k4,k5,n,ios
      character*80 :: dtitle ! dust info                            !ron
      logical :: metadata_is_separated
      integer :: im_gcm,jm_gcm ! gcm domain size for checking input file dims
      integer :: i_0,i_1,j_0,j_1 ! gcm domain bounds
      integer :: nrecs

      im_gcm = grid%im_world
      jm_gcm = grid%jm_world
      i_0 = grid%i_strt
      i_1 = grid%i_stop
      j_0 = grid%j_strt
      j_1 = grid%j_stop

      ! Count whether the dust file has at least 3 records.  This is
      ! to check whether plbdust, redust, rodust were placed in the
      ! same fortran record as gridded data tdust.  Mixtures of
      ! metadata and fields in fortran records prevent the use of
      ! parallel read routines.
      nrecs = 0
      do while(nrecs.lt.3)
        read(nrfu,iostat=ios)
        if(ios.ne.0) exit
        nrecs = nrecs + 1
      enddo
      rewind nrfu
      metadata_is_separated = nrecs==3

      if(am_i_root())
     & write(*,*) 'READ Offline DUST Distribution:'                   !ron
      read(nrfu) dtitle, imd, jmd, lmd, nsized, nmond ! offline dims  !ron
      if(am_i_root()) then
      write(*,*) trim(dtitle)                                         !ron
      write(*,*) 'Offline Dust Dims: imd, jmd, lmd, nsized, nmond: ', !ron
     *     imd, jmd, lmd, nsized, nmond                               !ron
      endif
      if (imd.eq.0)
     *     call stop_model('Please update dust file RADN6', 255)

      allocate( plbdust(lmd+1) )                                      !ron
      allocate( redust(nsized), rodust(nsized) )                      !ron
      allocate( tdust(i_0:i_1,j_0:j_1,lmd,nsized,nmond) )

      if(metadata_is_separated .and.
     &     im_gcm == imd .and. jm_gcm == jmd) then
        ! newer-format file allowing use of parallel read routine
        read(nrfu) plbdust, redust, rodust
        allocate(r8array(grid%i_strt_halo:grid%i_stop_halo,
     &                   grid%j_strt_halo:grid%j_stop_halo,
     &                   lmd*nsized*nmond))
        call dread_parallel(grid,nrfu,'DUSTFILE',r8array)
        n = 0
        do k5=1,nmond; do k4=1,nsized; do k3=1,lmd
          n = n + 1
          tdust(i_0:i_1,j_0:j_1,k3,k4,k5) = r8array(i_0:i_1,j_0:j_1,n)
        enddo; enddo; enddo
        deallocate(r8array)
      else
        allocate( tdust_glob(imd,jmd,lmd,nsized,nmond) )
        allocate(r8glob(imd,jmd))
        allocate(r8loc(grid%i_strt_halo:grid%i_stop_halo,
     &                 grid%j_strt_halo:grid%j_stop_halo))
        if(metadata_is_separated) then
          read(nrfu) plbdust, redust, rodust
          read(nrfu) tdust_glob
        else
          read(nrfu) plbdust, redust, rodust, tdust_glob
        endif
        do k5=1,nmond; do k4=1,nsized; do k3=1,lmd
          r8glob = tdust_glob(:,:,k3,k4,k5)
          call sample_latlon(imd,jmd,r8glob,r8loc)
          tdust(i_0:i_1,j_0:j_1,k3,k4,k5) = r8loc(i_0:i_1,j_0:j_1)
        enddo; enddo; enddo
        deallocate(tdust_glob,r8glob,r8loc)
      endif
      return
      end subroutine read_alloc_dust
#else
      subroutine read_alloc_dust(nrfu,imd,jmd,lmd,nsized,nmond,
     &     plbdust,rodust,redust,
     &     tdust
     &     )
!@sum This routine reads dust aerosol fields needed by the radiation
!@+   code in the prescribed-dust configuration of modelE.  This
!@+   is the netcdf version based on the traditional-I/O version
!@+   above.
!@auth R. Miller
!@auth M. Kelley extracted from RCOMP1, updated for netcdf
      use domain_decomp_atm, only: grid
      use pario, only : par_open,par_close
     &     ,get_dimlens,read_data,read_dist_data
      implicit none
      integer :: nrfu,imd,jmd,lmd,nsized,nmond
      real*8, dimension(:), pointer :: plbdust, redust, rodust
      real*4, dimension(:,:,:,:,:), pointer :: tdust
c
      integer :: fid,ndims,dlens(7)
      real*8, dimension(:,:,:,:,:), allocatable :: r8arr ! just for read-in
      integer :: im_gcm,jm_gcm ! gcm domain size for checking input file dims
      integer :: i_0,i_1,j_0,j_1 ! gcm domain bounds
      integer :: i_0h,i_1h,j_0h,j_1h
c
      im_gcm = grid%im_world
      jm_gcm = grid%jm_world
      i_0 = grid%i_strt
      i_1 = grid%i_stop
      j_0 = grid%j_strt
      j_1 = grid%j_stop
      i_0h = grid%i_strt_halo
      i_1h = grid%i_stop_halo
      j_0h = grid%j_strt_halo
      j_1h = grid%j_stop_halo

      fid = par_open(grid,'DUSTaer','read')
      call get_dimlens(grid,fid,'DUST',ndims,dlens)
      if(ndims.ne.5) call stop_model(
     &     'incorrect ndims for DUSTaer variable DUST',255)
      imd = dlens(1)
      jmd = dlens(2)
      if(imd.ne.im_gcm .or. jmd.ne.jm_gcm) call stop_model(
     &       'incorrect im,jm in DUSTaer',255)
      lmd    = dlens(3)
      nsized = dlens(4)
      nmond  = dlens(5)

      allocate( plbdust(lmd+1) )
      allocate( redust(nsized), rodust(nsized) )
      allocate( tdust(i_0:i_1,j_0:j_1,lmd,nsized,nmond) )
      allocate( r8arr(i_0h:i_1h,j_0h:j_1h,lmd,nsized,nmond) )

      call read_data(grid,fid,'plbdust',plbdust,bcast_all=.true.)
      call read_data(grid,fid,'redust',redust,bcast_all=.true.)
      call read_data(grid,fid,'rodust',rodust,bcast_all=.true.)

      call read_dist_data(grid,fid,'DUST',r8arr)
      tdust = r8arr(i_0:i_1,j_0:j_1,:,:,:)
      deallocate(r8arr)

      call par_close(grid,fid)

      return
      end subroutine read_alloc_dust
#endif
      end module DustParam_mod
