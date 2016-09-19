#include "rundeck_opts.h"

      MODULE GM_COM
!@sum  GM_COM variables related to GM isopycnal and Redi fluxes
!@auth Gavin Schmidt/Dan Collins
!@ver  2009/02/13
      USE CONSTANT, only : radius,omega,grav
      USE OCEAN, only : im,jm,lmo,lmm,lmu,lmv,dts,cospo,sinpo,ze,dxypo
     *     ,mo,dypo,dyvo,dxpo,dzo
      USE KPP_COM, only : kpl
      USE OCEAN_DYN, only  : dh
#ifdef OCN_GISS_MESO
      USE ODIAG, only : oij=>oij_loc,ij_gmsc
#endif

!      USE DOMAIN_DECOMP_1D, ONLY : grid, GET, HALO_UPDATE, NORTH, SOUTH,
      USE DOMAIN_DECOMP_1D, ONLY : GET, HALO_UPDATE, NORTH, SOUTH,
     *                          PACK_DATA, AM_I_ROOT,
     *                          UNPACK_DATA, GLOBALSUM
#ifndef OCN_GISS_MESO
     *                         ,ESMF_BCAST
#endif
      USE OCEANR_DIM, only : grid=>ogrid
      IMPLICIT NONE

      SAVE

!@var AI0,AI1,AI2,AI3 Cmponents of GM mixing coeff = F(isopycnal slopes)
!@var SIX0-SIX3,SIY0-SIY3: Slopes calculated from 4 triads of density.
!     REAL*8, DIMENSION(IM,JM,LMO) ::
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) ::
     *     ASX0,ASX1,ASX2,ASX3,AIX0,AIX1,AIX2,AIX3,
     *     ASY0,ASY1,ASY2,ASY3,AIY0,AIY1,AIY2,AIY3,
     *     S2X0,S2X1,S2X2,S2X3,S2Y0,S2Y1,S2Y2,S2Y3

!     REAL*8, DIMENSION(IM,JM,LMO) ::
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) ::
     *     DXZ, BXZ, CXZ,CDXZ, EXZ, DYZ,BYZ,
     *     CYZ,CDYZ,EYZ,CEXZ,CEYZ

!     REAL*8, DIMENSION(IM,JM,LMO) ::
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) ::
     *     BXX, BYY, BZZ, AZX, BZX, CZX,AEZX,
     *     EZX,CEZX, AZY, BZY,  CZY,AEZY, EZY,CEZY

!     REAL*8, DIMENSION(IM,JM,LMO) ::
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) ::
     *     BYDZV,BYDH,DZV

!     REAL*8, DIMENSION(IM,JM,LMO) ::
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) ::
     *     RHOX,RHOY,RHOMZ,BYRHOZ

!@var AINV Calculated Isopycnal thickness diffusion (m^2/s)
!@var ARIV Calculated Redi diffusion (m^2/s)
!     REAL*8, DIMENSION(IM,JM) ::
#ifdef OCN_GISS_MESO
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) ::
     *     AINV,ARIV
#else
      REAL*8, ALLOCATABLE, DIMENSION(:,:) ::
     *     AINV,ARIV
#endif

#ifndef OCN_GISS_MESO
!@var ARAI Scaling for Redi diffusion terms to GM diffusion term (1)
!@var QCROSS true if cross terms should be calculated
      REAL*8, PARAMETER :: ARAI = 1d0
      LOGICAL, PARAMETER :: QCROSS = .NOT. (ARAI.eq.1d0) ! i.e..FALSE.
#else
!@var RGMI ratio of GM (thickness) diffusivity to isoneutral (Redi) diffusivity
      REAL*8 :: RGMI
!@var QCROSS true if cross terms should be calculated
      LOGICAL :: QCROSS
#endif
!@var AMU = Visbeck scheme scaling parameter (1)
      REAL*8, PARAMETER :: AMU = 0.13d0
!@var SLIM = Upper limit of isopycnal slopes (stability parameter)
      REAL*8, PARAMETER :: SLIM=2d-3, BYSLIM=1./SLIM

      REAL*8, PARAMETER :: eps=TINY(1.D0)

!     REAL*8, DIMENSION(JM) ::
      REAL*8, ALLOCATABLE, DIMENSION(:) ::
     *     BYDYP,BYDXP,BYDYV

!@var VBAR specific volume (ref to mid point pressure)
!@var dVBARdZ specific volume vertical difference (ref to lower point pressure)
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: VBAR, dVBARdZ,G3D,S3D,P3D
!@var RHOZ1K density gradient over top 1km
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: RHOZ1K

!@var LUP level corresponding to 1km depth
      INTEGER :: LUP

      contains

      SUBROUTINE ALLOC_GM_COM
        implicit none
c**** allocate arrays
        allocate( ASX0(IM,grid%j_strt_halo:grid%j_stop_halo,LMO) )
        allocate( ASX1(IM,grid%j_strt_halo:grid%j_stop_halo,LMO) )
        allocate( ASX2(IM,grid%j_strt_halo:grid%j_stop_halo,LMO) )
        allocate( ASX3(IM,grid%j_strt_halo:grid%j_stop_halo,LMO) )
        allocate( AIX0(IM,grid%j_strt_halo:grid%j_stop_halo,LMO) )
        allocate( AIX1(IM,grid%j_strt_halo:grid%j_stop_halo,LMO) )
        allocate( AIX2(IM,grid%j_strt_halo:grid%j_stop_halo,LMO) )
        allocate( AIX3(IM,grid%j_strt_halo:grid%j_stop_halo,LMO) )
        allocate( ASY0(IM,grid%j_strt_halo:grid%j_stop_halo,LMO) )
        allocate( ASY1(IM,grid%j_strt_halo:grid%j_stop_halo,LMO) )
        allocate( ASY2(IM,grid%j_strt_halo:grid%j_stop_halo,LMO) )
        allocate( ASY3(IM,grid%j_strt_halo:grid%j_stop_halo,LMO) )
        allocate( AIY0(IM,grid%j_strt_halo:grid%j_stop_halo,LMO) )
        allocate( AIY1(IM,grid%j_strt_halo:grid%j_stop_halo,LMO) )
        allocate( AIY2(IM,grid%j_strt_halo:grid%j_stop_halo,LMO) )
        allocate( AIY3(IM,grid%j_strt_halo:grid%j_stop_halo,LMO) )
        allocate( S2X0(IM,grid%j_strt_halo:grid%j_stop_halo,LMO) )
        allocate( S2X1(IM,grid%j_strt_halo:grid%j_stop_halo,LMO) )
        allocate( S2X2(IM,grid%j_strt_halo:grid%j_stop_halo,LMO) )
        allocate( S2X3(IM,grid%j_strt_halo:grid%j_stop_halo,LMO) )
        allocate( S2Y0(IM,grid%j_strt_halo:grid%j_stop_halo,LMO) )
        allocate( S2Y1(IM,grid%j_strt_halo:grid%j_stop_halo,LMO) )
        allocate( S2Y2(IM,grid%j_strt_halo:grid%j_stop_halo,LMO) )
        allocate( S2Y3(IM,grid%j_strt_halo:grid%j_stop_halo,LMO) )

        allocate( DXZ(IM,grid%j_strt_halo:grid%j_stop_halo,LMO) )
        allocate( BXZ(IM,grid%j_strt_halo:grid%j_stop_halo,LMO) )
        allocate( CXZ(IM,grid%j_strt_halo:grid%j_stop_halo,LMO) )
        allocate( CDXZ(IM,grid%j_strt_halo:grid%j_stop_halo,LMO) )
        allocate( EXZ(IM,grid%j_strt_halo:grid%j_stop_halo,LMO) )
        allocate( DYZ(IM,grid%j_strt_halo:grid%j_stop_halo,LMO) )
        allocate( BYZ(IM,grid%j_strt_halo:grid%j_stop_halo,LMO) )
        allocate( CYZ(IM,grid%j_strt_halo:grid%j_stop_halo,LMO) )
        allocate( CDYZ(IM,grid%j_strt_halo:grid%j_stop_halo,LMO) )
        allocate( EYZ(IM,grid%j_strt_halo:grid%j_stop_halo,LMO) )
        allocate( CEXZ(IM,grid%j_strt_halo:grid%j_stop_halo,LMO) )
        allocate( CEYZ(IM,grid%j_strt_halo:grid%j_stop_halo,LMO) )

        allocate( BXX (IM,grid%j_strt_halo:grid%j_stop_halo,LMO) )
        allocate( BYY (IM,grid%j_strt_halo:grid%j_stop_halo,LMO) )
        allocate( BZZ (IM,grid%j_strt_halo:grid%j_stop_halo,LMO) )
        allocate( AZX (IM,grid%j_strt_halo:grid%j_stop_halo,LMO) )
        allocate( BZX (IM,grid%j_strt_halo:grid%j_stop_halo,LMO) )
        allocate( CZX (IM,grid%j_strt_halo:grid%j_stop_halo,LMO) )
        allocate( AEZX(IM,grid%j_strt_halo:grid%j_stop_halo,LMO) )
        allocate( EZX (IM,grid%j_strt_halo:grid%j_stop_halo,LMO) )
        allocate( CEZX(IM,grid%j_strt_halo:grid%j_stop_halo,LMO) )
        allocate( AZY (IM,grid%j_strt_halo:grid%j_stop_halo,LMO) )
        allocate( BZY (IM,grid%j_strt_halo:grid%j_stop_halo,LMO) )
        allocate( CZY (IM,grid%j_strt_halo:grid%j_stop_halo,LMO) )
        allocate( AEZY(IM,grid%j_strt_halo:grid%j_stop_halo,LMO) )
        allocate( EZY (IM,grid%j_strt_halo:grid%j_stop_halo,LMO) )
        allocate( CEZY(IM,grid%j_strt_halo:grid%j_stop_halo,LMO) )

        allocate( BYDZV(IM,grid%j_strt_halo:grid%j_stop_halo,LMO) )
        allocate( BYDH (IM,grid%j_strt_halo:grid%j_stop_halo,LMO) )
        allocate( DZV  (IM,grid%j_strt_halo:grid%j_stop_halo,LMO) )

        allocate( VBAR   (IM,grid%j_strt_halo:grid%j_stop_halo,LMO) )
        allocate( DVBARDZ(IM,grid%j_strt_halo:grid%j_stop_halo,LMO) )
        allocate( G3D   (IM,grid%j_strt_halo:grid%j_stop_halo,LMO) )
        allocate( S3D   (IM,grid%j_strt_halo:grid%j_stop_halo,LMO) )
        allocate( P3D   (IM,grid%j_strt_halo:grid%j_stop_halo,LMO) )
        allocate( RHOZ1K (IM,grid%j_strt_halo:grid%j_stop_halo) )

        allocate( RHOX   (IM,grid%j_strt_halo:grid%j_stop_halo,LMO) )
        allocate( RHOY   (IM,grid%j_strt_halo:grid%j_stop_halo,LMO) )
        allocate( RHOMZ  (IM,grid%j_strt_halo:grid%j_stop_halo,LMO) )
        allocate( BYRHOZ (IM,grid%j_strt_halo:grid%j_stop_halo,LMO) )

#ifdef OCN_GISS_MESO
        allocate( AINV (IM,grid%j_strt_halo:grid%j_stop_halo,LMO) )
        allocate( ARIV (IM,grid%j_strt_halo:grid%j_stop_halo,LMO) )
#else
        allocate( AINV (IM,grid%j_strt_halo:grid%j_stop_halo) )
        allocate( ARIV (IM,grid%j_strt_halo:grid%j_stop_halo) )
#endif

        allocate( BYDYP (grid%j_strt_halo:grid%j_stop_halo) )
        allocate( BYDXP (grid%j_strt_halo:grid%j_stop_halo) )
        allocate( BYDYV (grid%j_strt_halo:grid%j_stop_halo) )

      END SUBROUTINE ALLOC_GM_COM

      END MODULE GM_COM

      SUBROUTINE GMKDIF
#ifdef OCN_GISS_MESO
     * (RGMI_in)
#endif
!@sum GMKDIF calculates density gradients and tracer operators
!@+   for Redi and GM isopycnal skew fluxes
!@auth Dan Collins/Gavin Schmidt
      USE GM_COM
#ifdef OCN_GISS_MESO
      USE Dictionary_mod
#endif
      IMPLICIT NONE
#ifdef OCN_GISS_MESO
      REAL*8, INTENT(IN) :: RGMI_in
#endif
      INTEGER I,J,L,IM1

      INTEGER :: J_0, J_1, J_0S, J_1S, J_0STG, J_1STG, J_0H, J_1H
#ifndef OCN_GISS_MESO
      LOGICAL :: HAVE_NORTH_POLE
#endif
      INTEGER, SAVE :: IFIRST = 1

      if ( IFIRST.ne.0 ) then
        IFIRST = 0
        call ALLOC_GM_COM
      endif

#ifdef OCN_GISS_MESO
      RGMI = RGMI_in
      QCROSS = .NOT. (RGMI.eq.1d0)
#endif
c**** Extract domain decomposition info
      CALL GET(grid, J_STRT = J_0, J_STOP = J_1,
     &               J_STRT_SKP  = J_0S,   J_STOP_SKP  = J_1S,
     &               J_STRT_STGR=J_0STG, J_STOP_STGR=J_1STG,
     &               J_STRT_HALO = J_0H,   J_STOP_HALO = J_1H
#ifdef OCN_GISS_MESO
     &               )
#else
     &               ,HAVE_NORTH_POLE = HAVE_NORTH_POLE)
#endif

C**** Initialize SLOPES common block of coefficients
!$OMP PARALLEL DO  PRIVATE(L)
      DO L=1,LMO
        ASX0(:,:,L)=0. ;ASX1(:,:,L)=0. ; ASX2(:,:,L)=0. ; ASX3(:,:,L)=0.
        AIX0(:,:,L)=0. ;AIX1(:,:,L)=0. ; AIX2(:,:,L)=0. ; AIX3(:,:,L)=0.
        ASY0(:,:,L)=0. ;ASY1(:,:,L)=0. ; ASY2(:,:,L)=0. ; ASY3(:,:,L)=0.
        AIY0(:,:,L)=0. ;AIY1(:,:,L)=0. ; AIY2(:,:,L)=0. ; AIY3(:,:,L)=0.
        S2X0(:,:,L)=0. ;S2X1(:,:,L)=0. ; S2X2(:,:,L)=0. ; S2X3(:,:,L)=0.
        S2Y0(:,:,L)=0. ;S2Y1(:,:,L)=0. ; S2Y2(:,:,L)=0. ; S2Y3(:,:,L)=0.
         BXX(:,:,L)=0. ; BYY(:,:,L)=0. ;  BZZ(:,:,L)=0.
      END DO
!$OMP END PARALLEL DO

C**** Calculate all diffusivities
      CALL ISOSLOPE4

      CALL HALO_UPDATE(grid,AIY0 (:,grid%j_strt_halo:grid%j_stop_halo,
     *           :) , FROM=SOUTH)
      CALL HALO_UPDATE(grid,AIY1 (:,grid%j_strt_halo:grid%j_stop_halo,
     *           :) , FROM=SOUTH)

C**** Surface layer is calculated directly with appropriate AI,S = 0
!$OMP PARALLEL DO  PRIVATE(L,J,IM1,I)
      DO L=1,LMO       !Calculating all layers!
      DO J=J_0STG,J_1STG
      IM1 = IM
      DO I=1,IM
      IF(LMM(I,J).lt.L) GO TO 110
C**** FXZPR  Tracer flux at center point in x-direction from tracer
C**** z-gradients by right-side slopes (F=A*S* -gradT)
C**** FXZPR(I,J,L) = AI0(I,J,L)*SIX0(I,J,L)*(TRM(I,J,LP1)/MO(I,J,LP1) -
C****                                        TRM(I,J,L)  /MO(I,J,L))
C****              + AI1(I,J,L)*SIX1(I,J,L)*(TRM(I,J,L)  /MO(I,J,L)   -
C****                                        TRM(I,J,LM1)/MO(I,J,LM1))
C**** FXZPL(I,J,L) = AI2(I,J,L)*SIX2(I,J,L)*
C****                              (TRM(I,J,LP1)/MO(I,J,LP1) -
C****                               TRM(I,J,L)  /MO(I,J,L))
C****                + AI3(I,J,L)*SIX3(I,J,L)*
C****                              (TRM(I,J,L)  /MO(I,J,L)   -
C****                               TRM(I,J,LM1)/MO(I,J,LM1))
C**** FXZPL X-direction tracer flux at center point from z-gradients
C**** by left side slopes.
C**** FXZ   X-direction tracer flux from tracer z-gradients, at center
C**** of right grid wall.
C**** FXZ(I,J,L) = (FXZPR(I,J,L) + FXZPL(IP1,J,L))*1/4 (Four slopes)
C**** TRMXZ+(I,J,L) = TRMXZ-(I,J,L) - (FXZ(I,J,L)-FXZ(IM1,J,L))*DTS/MO
C****
C**** Calculate coefficients of Tracer Operator
C**** ASX0,ASX1,ASX2,ASX3 (ASY0..) A*S for 4 triads in x (y) direction
C**** FXZ   coefficients
      IF (QCROSS) THEN
        DXZ(I,J,L)     = -ASX1(I,J,L) * BYDH(I,J,L)
        CDXZ(IM1,J,L)  = -ASX3(I,J,L) * BYDH(I,J,L)
        BXZ(I,J,L)     = (ASX1(I,J,L) - ASX0(I,J,L))*BYDH(I,J,L)
        CXZ(IM1,J,L)   = (ASX3(I,J,L) - ASX2(I,J,L))*BYDH(I,J,L)
        EXZ(I,J,L)     =  ASX0(I,J,L) * BYDH(I,J,L)
        CEXZ(IM1,J,L)  =  ASX2(I,J,L) * BYDH(I,J,L)
C**** Flux in Y-direction
C**** ASY1(I,JM-1,L), ASY0(I,JM-1,LM1) and ASY3(I,2,L), ASY2(I,2,LM1)??
        DYZ(I,J,L)     = -ASY1(I,J,L) * BYDH(I,J,L)
        BYZ(I,J,L)     = (ASY1(I,J,L) - ASY0(I,J,L))*BYDH(I,J,L)
        CDYZ(I,J-1,L)  = -ASY3(I,J,L) * BYDH(I,J,L)
        CYZ(I,J-1,L)   = (ASY3(I,J,L) - ASY2(I,J,L))*BYDH(I,J,L)
        CEYZ(I,J-1,L)  =  ASY2(I,J,L) * BYDH(I,J,L)
        EYZ(I,J,L)     =  ASY0(I,J,L) * BYDH(I,J,L)
      END IF
C**** Diagonal terms of Horizontal fluxes are decoupled
C**** from downgradients (-ve gradT)  !!!Sign!!!
      BXX(IM1,J,L) = (AIX2(I,J,L) + AIX0(IM1,J  ,L) +
     *                AIX3(I,J,L) + AIX1(IM1,J  ,L))
      BYY(I,J-1,L) = (AIY2(I,J,L) + AIY0(I  ,J-1,L) +
     *                AIY3(I,J,L) + AIY1(I  ,J-1,L))
#ifndef OCN_GISS_MESO
      IF(KPL(I,J).gt.L) GO TO 110
#else
C**** Z-direction fluxes
C**** Diagonal (Includes AI and DYV/DYP)
      IF (L.gt.KPL(I,J)) THEN
#endif
      IF (L.gt.1) BZZ(I,J,L-1) = (S2X1(I,J,L  ) + S2X3(I,J,L  ) +
     *                            S2X0(I,J,L-1) + S2X2(I,J,L-1) +
     *                            S2Y1(I,J,L  ) + S2Y3(I,J,L  ) +
     *                            S2Y0(I,J,L-1) + S2Y2(I,J,L-1))
C****
C**** Off-diagonal
C**** FZXbot(I,J,L) = ASX0(I,J,L)*(TR(I  ,J,L) - TR(IP1,J,L))
C****               + ASX2(I,J,L)*(TR(IM1,J,L) - TR(I  ,J,L))
C**** FZXtop(I,J,L) = ASX1(I,J,L)*(TR(I  ,J,L) - TR(IP1,J,L))
C****               + ASX3(I,J,L)*(TR(IM1,J,L) - TR(I  ,J,L))
C**** Coeficients for (I,J,L) are multiplied by T(IP1,J,L) and T(IM1)
      AZX(I,J,L)   =  ASX2(I,J,L)
      BZX(I,J,L)   =  ASX0(I,J,L) - ASX2(I,J,L)
      CZX(I,J,L)   = -ASX0(I,J,L)
      AEZX(I,J,L)  =  ASX3(I,J,L)
      EZX(I,J,L)   =  ASX1(I,J,L) - ASX3(I,J,L)
      CEZX(I,J,L)  = -ASX1(I,J,L)
C**** y-coefficients *DYV(J-1)/DYV(J-1) or *DYV(J)/DYV(J)
      AZY(I,J,L)   =  ASY2(I,J,L)
      BZY(I,J,L)   =  ASY0(I,J,L) - ASY2(I,J,L)
      CZY(I,J,L)   = -ASY0(I,J,L)
      AEZY(I,J,L)  =  ASY3(I,J,L)
      EZY(I,J,L)   =  ASY1(I,J,L) - ASY3(I,J,L)
      CEZY(I,J,L)  = -ASY1(I,J,L)
#ifdef OCN_GISS_MESO
      ELSE
      IF (L.gt.1) BZZ(I,J,L-1) = 0.
      AZX(I,J,L)   = 0.
      BZX(I,J,L)   = 0.
      CZX(I,J,L)   = 0.
      AEZX(I,J,L)  = 0.
      EZX(I,J,L)   = 0.
      CEZX(I,J,L)  = 0.
      AZY(I,J,L)   = 0.
      BZY(I,J,L)   = 0.
      CZY(I,J,L)   = 0.
      AEZY(I,J,L)  = 0.
      EZY(I,J,L)   = 0.
      CEZY(I,J,L)  = 0.
      END IF
#endif
  110 IM1 = I
      END DO
      END DO
      END DO
!$OMP END PARALLEL DO

!     calculate CDYZ,CYZ,CEYZ,BYY at J_1 from J_1 + 1
      IF (QCROSS) THEN
        CALL HALO_UPDATE(grid,ASY3 (:,grid%j_strt_halo:grid%j_stop_halo,
     *           :) , FROM=NORTH)
        CALL HALO_UPDATE(grid,BYDH (:,grid%j_strt_halo:grid%j_stop_halo,
     *           :) , FROM=NORTH)
        CALL HALO_UPDATE(grid,ASY2 (:,grid%j_strt_halo:grid%j_stop_halo,
     *           :) , FROM=NORTH)
      END IF
      CALL HALO_UPDATE(grid,AIY2 (:,grid%j_strt_halo:grid%j_stop_halo,
     *           :) , FROM=NORTH)
      CALL HALO_UPDATE(grid,AIY3 (:,grid%j_strt_halo:grid%j_stop_halo,
     *           :) , FROM=NORTH)
      J=J_1STG+1
      if(J.lt.JM) then
      DO L=1,LMO       !Calculating for all layers & all I's
        IM1 = IM
        DO I=1,IM
          IF(LMM(I,J).lt.L) GO TO 111
          IF (QCROSS) THEN
            CDYZ(I,J-1,L)  = -ASY3(I,J,L) * BYDH(I,J,L)
            CYZ(I,J-1,L)   = (ASY3(I,J,L) - ASY2(I,J,L))*BYDH(I,J,L)
            CEYZ(I,J-1,L)  =  ASY2(I,J,L) * BYDH(I,J,L)
          END IF
          BYY(I,J-1,L) = (AIY2(I,J,L) + AIY0(I  ,J-1,L) +
     *                    AIY3(I,J,L) + AIY1(I  ,J-1,L))
  111     IM1 = I
        END DO
      END DO
      endif

C**** End of Main Loop of GMKDIF
      RETURN
      END SUBROUTINE GMKDIF
C****
      SUBROUTINE GMFEXP (TRM, TXM, TYM, TZM, QLIMIT, GIJL)
!@sum  GMFEXP apply GM fluxes to tracer quantities
!@auth Gavin Schmidt/Dan Collins
!@ver  1.0
      USE GM_COM
      IMPLICIT NONE
!     REAL*8, DIMENSION(IM,JM,LMO), INTENT(INOUT) :: TRM,TXM,TYM,TZM
      REAL*8, DIMENSION(IM,grid%j_strt_halo:grid%j_stop_halo,LMO),
     *        INTENT(INOUT) :: TRM,TXM,TYM,TZM
      LOGICAL, INTENT(IN) :: QLIMIT
!@var GIJL Tracer flux
!     REAL*8, DIMENSION(IM,JM,LMO,3), INTENT(INOUT) :: GIJL
      REAL*8, DIMENSION(IM,grid%j_strt_halo:grid%j_stop_halo,LMO,3),
     *         INTENT(INOUT) :: GIJL
!     REAL*8, DIMENSION(IM,JM,LMO) :: TR
      REAL*8, DIMENSION(IM,grid%j_strt_halo:grid%j_stop_halo,LMO) :: TR
!     REAL*8, DIMENSION(IM,JM,LMO) :: FXX,FXZ,FYY,FYZ,FZZ,FZX,FZY
      REAL*8, DIMENSION(IM,grid%j_strt_halo:grid%j_stop_halo,LMO) ::
     *         FXX,FXZ,FYY,FYZ,FZZ,FZX,FZY
      REAL*8 MOFX, MOFY, MOFZ, RFXT, RFYT, RFZT, DT4, DT4DY, DT4DX,
     *     STRNP
      REAL*8, DIMENSION(IM,grid%j_strt_halo:grid%j_stop_halo,LMO) ::
     *         flux_x, flux_y, flux_z
      INTEGER I,J,L,IM1,IP1

      INTEGER :: J_0, J_1, J_0S, J_1S, J_0STG, J_1STG, J_0H, J_1H
#ifndef OCN_GISS_MESO
      LOGICAL :: HAVE_NORTH_POLE
#endif

c**** Extract domain decomposition info
      CALL GET(grid, J_STRT = J_0, J_STOP = J_1,
     &               J_STRT_SKP  = J_0S,   J_STOP_SKP  = J_1S,
     &               J_STRT_STGR=J_0STG, J_STOP_STGR=J_1STG,
     &               J_STRT_HALO = J_0H,   J_STOP_HALO = J_1H
#ifdef OCN_GISS_MESO
     &               )
#else
     &              ,HAVE_NORTH_POLE = HAVE_NORTH_POLE)
#endif

      DT4 = 0.25*DTS
C**** Explicit calculation of "fluxes" (R(T))
C**** Calculate tracer concentration
!$OMP PARALLEL DO  PRIVATE(L,J,I)
      DO L = 1,LMO
        DO J = J_0,J_1
          DO I = 1,IM
            IF(L.le.LMM(I,J)) THEN
              TR(I,J,L)=TRM(I,J,L)/(DXYPO(J)*MO(I,J,L))
            ELSE
              TR(I,J,L)=0.0
            END IF
            FXX(I,J,L) = 0. ; FXZ(I,J,L) = 0.
            FYY(I,J,L) = 0. ; FYZ(I,J,L) = 0.
            FZZ(I,J,L) = 0. ; FZX(I,J,L) = 0. ; FZY(I,J,L) = 0.
          END DO
        END DO
      END DO
!$OMP  END PARALLEL DO

      CALL HALO_UPDATE(grid,TR(:,grid%j_strt_halo:grid%j_stop_halo,:) ,
     *                 FROM=NORTH+SOUTH)

!$OMP PARALLEL DO  PRIVATE(L,J,DT4DY,DT4DX,IM1,I,IP1)
      DO L=1,LMO
      DO J=J_0S,J_1S
C**** 1/DYPO(DXP) from Y(X) gradient of T for FZY(FZX)
      DT4DY = 0.25*DTS*BYDYP(J)
      DT4DX = 0.25*DTS*BYDXP(J)
      IM1 = IM - 1
      I   = IM
      DO IP1=1,IM
      IF(LMM(I,J).le.0) GO TO 530
C**** Loop for Fluxes in X-direction
      IF(LMU(IM1,J).lt.L) GO TO 510
C**** Calculate fluxes
C**** Diagonal (Includes 1/DXP for gradTX, Not divFX)
      FXX(IM1,J,L) = (DT4DX * BXX(IM1,J,L)) * (TR(IM1,J,L) - TR(I,J,L))
C**** Off-diagonal
      IF (QCROSS) THEN
        IF(L.gt.1) FXZ(IM1,J,L) = FXZ(IM1,J,L) +
     *       DT4 * (DXZ(IM1,J,L) * TR(IM1,J,L-1) +
     *             CDXZ(IM1,J,L) * TR(I,J,L-1))
        FXZ(IM1,J,L) = FXZ(IM1,J,L) +
     *       DT4 * (BXZ(IM1,J,L) * TR(IM1,J,L) +
     *              CXZ(IM1,J,L)  * TR(I,J,L))
C**** Skip for L+1 greater than LMM(IM1,J)
        IF(LMM(IM1,J).gt.L) FXZ(IM1,J,L) =  FXZ(IM1,J,L) +
     *       DT4 * EXZ(IM1,J,L) * TR(IM1,J,L+1)
C**** Skip for L+1 greater than LMM(I,J)
        IF(LMM(I,J).gt.L) FXZ(IM1,J,L) =  FXZ(IM1,J,L) +
     *       DT4 * CEXZ(IM1,J,L) * TR(I,J,L+1)
        
#ifdef OCN_GISS_MESO
        FXZ(IM1,J,L) =  -FXZ(IM1,J,L)*(RGMI-1d0)
#else
        FXZ(IM1,J,L) =  FXZ(IM1,J,L)*(1d0-ARAI)
#endif
      END IF
  510 CONTINUE
C**** END of FX
C**** Loop for Fluxes in Y-direction
      IF(LMV(I,J).lt.L) GO TO 520
C**** Calculate fluxes in Y-direction
C**** Diagonal:  Must be divided by DYPO(J) for divF!
      FYY(I,J,L) = DT4*BYY(I,J,L)*(TR(I,J,L) - TR(I,J+1,L))*BYDYV(J)
C**** Off-diagonal:   Divide by DYPO(J) for divF!
C**** YZ terms are divided by appropriate DZV for Z-gradient of T
      IF (QCROSS) THEN
        IF(L.gt.1) FYZ(I,J,L) = FYZ(I,J,L) +
     *       DT4 * (DYZ(I,J,L) * TR(I,J  ,L-1) +
     *             CDYZ(I,J,L) * TR(I,J+1,L-1))
        FYZ(I,J,L) = FYZ(I,J,L) +
     *       DT4 * (BYZ(I,J,L) * TR(I,J  ,L) +
     *              CYZ(I,J,L) * TR(I,J+1,L))
C**** Skip for L+1 greater than LMM(I,J)
        IF(LMM(I,J).gt.L) FYZ(I,J,L) =  FYZ(I,J,L) +
     *       DT4 * EYZ(I,J,L) * TR(I,J,L+1)
C**** Skip for L+1 greater than LMM(I,J+1)
        IF(LMM(I,J+1).gt.L) FYZ(I,J,L) =  FYZ(I,J,L) +
     *       DT4 * CEYZ(I,J,L) * TR(I,J+1,L+1)
#ifdef OCN_GISS_MESO
        FYZ(I,J,L) =  -FYZ(I,J,L) *(RGMI-1d0)
#else
        FYZ(I,J,L) =  FYZ(I,J,L) *(1d0-ARAI)
#endif
      END IF
  520 CONTINUE
C**** END of FY
C**** Loop for Fluxes in Z-direction
#ifndef OCN_GISS_MESO
      IF(KPL(I,J).gt.L) GO TO 530
#endif
      IF(LMM(I,J).le.L) GO TO 530
C**** Calculate fluxes in Z-direction
C**** Diagonal      :  Need BYDH for divF!
#ifdef OCN_GISS_MESO
      IF(KPL(I,J).le. L) then
#endif
      FZZ(I,J,L) = DT4*BZZ(I,J,L)*(TR(I,J,L+1)-TR(I,J,L))*BYDZV(I,J,L)
#ifdef OCN_GISS_MESO
      endif
#endif
C**** Off-diagonal X:  May need to use IM2,IM1,I and F(IM1)
      FZX(I,J,L) = DT4DX * (BZX(I,J,L) * TR(I,J,L) +
     *     AZX(I,J,L) * TR(IM1,J,L) + CZX(I,J,L) * TR(IP1,J,L) +
     *     EZX(I,J,L+1) * TR(I,J,L+1)) !EZX=0 for LMU(IorIM1)=0
C**** Skip for L+1 greater than LMM(I,J)
      IF(LMM(IM1,J).gt.L) FZX(I,J,L) =  FZX(I,J,L) +
     *     DT4DX * AEZX(I,J,L+1) * TR(IM1,J,L+1)
C**** Skip for L+1 greater than LMM(I,J+1)
      IF(LMM(IP1,J).gt.L) FZX(I,J,L) =  FZX(I,J,L) +
     *     DT4DX * CEZX(I,J,L+1) * TR(IP1,J,L+1)
C**** Off-diagonal Y:  May need to use JM2,J-1,J and F(J-1)
      FZY(I,J,L) = DT4DY * (BZY(I,J,L) * TR(I,J,L) +
     *     AZY(I,J,L) * TR(I,J-1,L) + CZY(I,J,L) * TR(I,J+1,L) +
     *     EZY(I,J,L+1) * TR(I,J,L+1)) !EZY=0 for LMV(JorJ-1)=0
C**** Skip for L+1 greater than LMM(I,J)
      IF(LMM(I,J-1).gt.L) FZY(I,J,L) =  FZY(I,J,L) +
     *     DT4DY * AEZY(I,J,L+1) * TR(I,J-1,L+1)
C**** Skip for L+1 greater than LMM(I,J+1)
      IF(LMM(I,J+1).gt.L) FZY(I,J,L) =  FZY(I,J,L) +
     *     DT4DY * CEZY(I,J,L+1) * TR(I,J+1,L+1)
  530 IM1 = I
      I = IP1
      END DO
      END DO
      END DO
!$OMP  END PARALLEL DO

!     HALO data used in computeFluxes 
      CALL HALO_UPDATE(grid,MO(:,grid%j_strt_halo:grid%j_stop_halo
     *           ,:) , FROM=SOUTH+NORTH)
#ifndef OCN_GISS_MESO
      CALL HALO_UPDATE(grid,DXYPO(grid%j_strt_halo:grid%j_stop_halo)
     *               , FROM=SOUTH+NORTH)
#endif
      CALL HALO_UPDATE(grid,BYY(:,grid%j_strt_halo:grid%j_stop_halo
     *           ,:) , FROM=SOUTH)

      CALL HALO_UPDATE(grid,FYY(:,grid%j_strt_halo:grid%j_stop_halo
     *           ,:) , FROM=SOUTH)
      CALL HALO_UPDATE(grid,FYZ(:,grid%j_strt_halo:grid%j_stop_halo
     *           ,:) , FROM=SOUTH)

      call computeFluxes(DT4, flux_x, flux_y, flux_z, TXM, TYM,
     &                   TZM, FXX, FXZ, FYY, FYZ, FZZ, FZX, FZY)


      IF (QLIMIT) THEN
        !serial version to preserve original algorithm for TRM update 
        call wrapAdjustFluxes(TRM, flux_x, flux_y, flux_z, GIJL)
      ELSE
        CALL HALO_UPDATE(grid,TRM(:,grid%j_strt_halo:grid%j_stop_halo
     *           ,:) , FROM=SOUTH)
      CALL HALO_UPDATE(grid,flux_x(:,grid%j_strt_halo:grid%j_stop_halo
     *           ,:) , FROM=NORTH+SOUTH)
      CALL HALO_UPDATE(grid,flux_y(:,grid%j_strt_halo:grid%j_stop_halo
     *           ,:) , FROM=NORTH+SOUTH)
        call addFluxes (TRM, flux_x, flux_y, flux_z, GIJL)
      END IF

      RETURN
      END SUBROUTINE GMFEXP
C****
      SUBROUTINE computeFluxes(DT4, flux_x, flux_y, flux_z, TXM, TYM,
     &                         TZM, FXX, FXZ, FYY, FYZ, FZZ, FZX, FZY)
!@sum  computes GM fluxes for tracer quantities
!@ver  1.0
      USE GM_COM, ONLY: grid,GET,IM,JM,LMO,LMM,LMU,LMV,
     &                  DXYPO,MO, kpl,
#ifdef OCN_GISS_MESO
     &                  BXX, BYY, BZZ, BYDH, RGMI, BYDXP, BYDYP
#else
     &                  BXX, BYY, BZZ, BYDH, ARAI, BYDXP, BYDYP
#endif
      IMPLICIT NONE
      REAL*8, INTENT(IN) :: DT4
      REAL*8, DIMENSION(IM,grid%j_strt_halo:grid%j_stop_halo,LMO),
     &             INTENT(INOUT) :: flux_x, flux_y, flux_z
      REAL*8, DIMENSION(IM,grid%j_strt_halo:grid%j_stop_halo,LMO),
     &             INTENT(INOUT) :: TXM,TYM,TZM
      REAL*8, DIMENSION(IM,grid%j_strt_halo:grid%j_stop_halo,LMO),
     &             INTENT(IN) :: FXX,FXZ,FYY,FYZ,FZZ,FZX,FZY
      REAL*8 MOFX, MOFY, MOFZ, RFXT, RFYT, RFZT
      INTEGER I,J,L,IM1

      INTEGER :: J_0S, J_1S
      LOGICAL :: HAVE_NORTH_POLE
#ifdef OCN_GISS_MESO
      LOGICAL :: HAVE_SOUTH_POLE
#endif

c**** Extract domain decomposition info
      CALL GET(grid, J_STRT_SKP  = J_0S,   J_STOP_SKP  = J_1S
     &              ,HAVE_NORTH_POLE = HAVE_NORTH_POLE
#ifndef OCN_GISS_MESO
     &               )
#else
     &              ,HAVE_SOUTH_POLE = HAVE_SOUTH_POLE)
#endif

C**** Summation of explicit fluxes to TR, (R(T))
      flux_x = 0.0; flux_y = 0.0; flux_z = 0.0;

!$OMP PARALLEL DO  PRIVATE(L,J,IM1,I,MOFX,RFXT,MOFY,RFYT,MOFZ,RFZT)
      DO L=1,LMO
C**** Non-Polar boxes
      DO J=J_0S,J_1S
      IM1 = IM
      DO I=1,IM
      IF(LMM(I,J).le.0) GO TO 610
C**** Loop for Fluxes in X-direction
      IF(LMU(IM1,J).ge.L) THEN
        MOFX = (MO(IM1,J,L) +  MO(I  ,J,L)) * DXYPO(J) *BYDXP(J) *0.5
        RFXT =(FXX(IM1,J,L) + FXZ(IM1,J,L))*MOFX

        flux_x(I,J,L) =RFXT
      endif
C**** Gradient fluxes in X direction affected by diagonal terms
      IF (L.le.LMM(I,J)) TXM(I,J,L) = (TXM(I,J,L) -
     *     3.*(FXX(IM1,J,L)+FXX(I,J,L))*MO(I,J,L)*DXYPO(J)*BYDXP(J))/
     *     (1.+6.*DT4*(BXX(IM1,J,L)+BXX(I,J,L))*BYDXP(J)**2)

C**** Loop for Fluxes in Y-direction
      IF(LMV(I,J-1).ge.L) THEN
        MOFY =((MO(I,J-1,L)*BYDYP(J-1)*DXYPO(J-1)) +
     *         (MO(I,J  ,L)*BYDYP(J  )*DXYPO(J  )))*0.5
        RFYT =(FYY(I,J-1,L) + FYZ(I,J-1,L))*MOFY

        flux_y(I,J,L) =RFYT
      endif
C**** Gradient fluxes in Y direction affected by diagonal terms
      IF (L.le.LMM(I,J)) TYM(I,J,L) = (TYM(I,J,L) -
     *     3.*(FYY(I,J-1,L)+FYY(I,J,L))*MO(I,J,L)*DXYPO(J)*BYDYP(J))/
     *     (1.+6.*DT4*(BYY(I,J-1,L)+BYY(I,J,L))*BYDYP(J)**2)

C**** Summation of explicit fluxes to TR, (R(T)) --- Z-direction
C**** Loop for Fluxes in Z-direction
#ifndef OCN_GISS_MESO
      IF(KPL(I,J).gt.L) GO TO 610
#endif
      IF(LMM(I,J).lt.L) GO TO 610
      IF(LMM(I,J).gt.L) THEN
C**** Calculate new tracer/salinity/enthalpy
        MOFZ =((MO(I,J,L+1)*BYDH(I,J,L+1)) +
     *         (MO(I,J,L  )*BYDH(I,J,L  ))) * DXYPO(J) *0.5
#ifndef OCN_GISS_MESO
        RFZT =(FZZ(I,J,L) +(FZX(I,J,L)+FZY(I,J,L))*(1.d0+ARAI))*MOFZ
#else
        RFZT =(FZZ(I,J,L) +(FZX(I,J,L)+FZY(I,J,L))*(1.d0+RGMI))*MOFZ
#endif

        flux_z(I,J,L) =RFZT
      endif
C**** Gradient fluxes in Z direction affected by diagonal terms
c      IF (L.gt.1) THEN
c        TZM(I,J,L) = (TZM(I,J,L) - 3.*(FZZ(I,J,L)+FZZ(I,J,L-1))*
c     *       MO(I,J,L)*DXYPO(J)*BYDH(I,J,L))/(1.+6.*DT4*(BZZ(I,J,L)
c     *       +BZZ(I,J,L-1))*BYDH(I,J,L)**2)
c      ELSE
c        TZM(I,J,L) = (TZM(I,J,L) - 3.*FZZ(I,J,L)*MO(I,J,L)*DXYPO(J)
c     *       *BYDH(I,J,L))/(1.+12.*DT4*BZZ(I,J,L)*BYDH(I,J,L)**2)
c      END IF

C**** END of I and J loops
  610 IM1 = I
      END DO
      END DO

      IF(HAVE_NORTH_POLE) THEN

C****   North Polar box
C****   Fluxes in Y-direction
        DO I=1,IM
          IF(LMV(I,JM-1).ge.L) THEN
            MOFY =((MO(I,JM-1,L)*BYDYP(JM-1)*DXYPO(JM-1)) +
     *             (MO(1,JM  ,L)*BYDYP(JM  )*DXYPO(JM  )))*0.5
            RFYT =(FYY(I,JM-1,L) + FYZ(I,JM-1,L))*MOFY
            flux_y(I,JM,L) = RFYT
          END IF
        END DO

C****   Loop for Fluxes in Z-direction
#ifndef OCN_GISS_MESO
        IF(KPL(1,JM).gt.L) CYCLE  !GO TO 620
        IF(LMM(1,JM).lt.L) CYCLE  !GO TO 620
#else
        IF(LMM(1,JM).ge.L) THEN
#endif
        IF(LMM(1,JM).gt.L) THEN
C****     Calculate new tracer/salinity/enthalpy
          MOFZ =((MO(1,JM,L+1)*BYDH(1,JM,L+1)) +
     *           (MO(1,JM,L  )*BYDH(1,JM,L  ))) * DXYPO(JM) *0.5
#ifndef OCN_GISS_MESO
          RFZT =(FZZ(1,JM,L)+FZY(1,JM,L)*(1d0+ARAI))*MOFZ
#else
          RFZT =(FZZ(1,JM,L)+FZY(1,JM,L)*(1d0+RGMI))*MOFZ
#endif
          flux_z(1,JM,L) = RFZT
        END IF
C****   Gradient fluxes in Z direction affected by diagonal terms
c        IF (L.gt.1) THEN
c          TZM(1,JM,L) = (TZM(1,JM,L) - 3.*(FZZ(1,JM,L)+FZZ(1,JM,L-1))*
c     *      MO(1,JM,L)*DXYPO(JM)*BYDH(1,JM,L))/(1.+6.*DT4*(BZZ(1,JM,L)
c     *      +BZZ(1,JM,L-1))*BYDH(1,JM,L)**2)
c        ELSE
c          TZM(1,JM,L) = (TZM(1,JM,L) - 3.*FZZ(1,JM,L)*MO(1,JM,L)*
c     &      DXYPO(JM) *BYDH(1,JM,L))/(1.+12.*DT4*BZZ(1,JM,L)*
c     &      BYDH(1,JM,L)**2)
c        END IF

#ifdef OCN_GISS_MESO
        END IF
#endif
      END IF   ! HAVE_NORTH_POLE

#ifdef OCN_GISS_MESO
 615  IF(HAVE_SOUTH_POLE) THEN

C****   South Polar box
C****   Fluxes in Y-direction
        DO I=1,IM
          IF(LMV(I,2).ge.L) THEN
            MOFY =((MO(I,2,L)*BYDYP(2)*DXYPO(2)) +
     *             (MO(1,1,L)*BYDYP(1)*DXYPO(1)))*0.5
            RFYT =(FYY(I,2,L) + FYZ(I,2,L))*MOFY
            flux_y(I,1,L) = RFYT
          END IF
        END DO

C****   Loop for Fluxes in Z-direction
        IF(LMM(1,1).ge.L) THEN
        IF(LMM(1,1).gt.L) THEN
C****     Calculate new tracer/salinity/enthalpy
          MOFZ =((MO(1,1,L+1)*BYDH(1,1,L+1)) +
     *           (MO(1,1,L  )*BYDH(1,1,L  ))) * DXYPO(1) *0.5
#ifndef OCN_GISS_MESO
          RFZT =(FZZ(1,1,L)+FZY(1,1,L)*(1d0+ARAI))*MOFZ
#else
          RFZT =(FZZ(1,1,L)+FZY(1,1,L)*(1d0+RGMI))*MOFZ
#endif
          flux_z(1,1,L) = RFZT
        END IF
C****   Gradient fluxes in Z direction affected by diagonal terms
c        IF (L.gt.1) THEN
c          TZM(1,1,L) = (TZM(1,1,L) - 3.*(FZZ(1,1,L)+FZZ(1,1,L-1))*
c     *      MO(1,1,L)*DXYPO(1)*BYDH(1,1,L))/(1.+6.*DT4*(BZZ(1,1,L)
c     *      +BZZ(1,1,L-1))*BYDH(1,1,L)**2)
c        ELSE
c          TZM(1,1,L) = (TZM(1,1,L) - 3.*FZZ(1,1,L)*MO(1,1,L)*
c     &      DXYPO(1) *BYDH(1,1,L))/(1.+12.*DT4*BZZ(1,1,L)*
c     &      BYDH(1,1,L)**2)
c        END IF
        END IF
      END IF   ! HAVE_SOUTH_POLE
#endif

 620  END DO   ! L loop


      RETURN
      END SUBROUTINE computeFluxes

      SUBROUTINE wrapAdjustFluxes(TRM, flux_x, flux_y, flux_z, GIJL)
!@sum Applies limiter to GM flux divergence to prevent tracer quantities
!@+   from becoming negative.
!@ver  1.0
      USE GM_COM
      IMPLICIT NONE
      REAL*8, DIMENSION(IM,grid%j_strt_halo:grid%j_stop_halo,LMO),
     &        INTENT(INOUT) :: TRM, flux_x, flux_y, flux_z
      REAL*8, DIMENSION(IM,grid%j_strt_halo:grid%j_stop_halo,LMO,3),
     &        INTENT(INOUT) :: GIJL

      REAL*8, DIMENSION(IM,grid%j_strt_halo:grid%j_stop_halo,LMO)
     &        :: conv
      real*8, dimension(grid%j_strt_halo:grid%j_stop_halo) ::
     &     convadj_j,convpos_j
      real*8 :: sumadj,sumpos,posadj
      integer :: i,j,l

      INTEGER :: J_0S, J_1S, J_0, J_1
      LOGICAL :: HAVE_NORTH_POLE
#ifdef OCN_GISS_MESO
      LOGICAL :: HAVE_SOUTH_POLE
#endif

c**** Extract domain decomposition info
      CALL GET(grid,
     &     J_STRT=J_0, J_STOP=J_1,
     &     J_STRT_SKP=J_0S, J_STOP_SKP=J_1S,
     &     HAVE_NORTH_POLE = HAVE_NORTH_POLE
#ifndef OCN_GISS_MESO
     &     )
#else
     &    ,HAVE_SOUTH_POLE = HAVE_SOUTH_POLE)
#endif


      call halo_update(grid,flux_y,from=north)

c
c If any flux convergences are adjusted, these flux diagnostics
c will be wrong.  Since the adjustments are typically rare in
c time and space, the diagnostics error will be small.
c Shift x/y indices by 1 for consistency with ocean dynamics
c conventions.
c
      do l=1,lmo
        do j=j_0,j_1
          do i=1,im-1
            gijl(i,j,l,1) = gijl(i,j,l,1) + flux_x(i+1,j,l)
          enddo
          i=im
            gijl(i,j,l,1) = gijl(i,j,l,1) + flux_x(  1,j,l)
          do i=1,im
            gijl(i,j,l,3) = gijl(i,j,l,3) + flux_z(i,j,l)
          enddo
        enddo
        do j=j_0,min(j_1,jm-1)
          do i=1,im
            gijl(i,j,l,2) = gijl(i,j,l,2) + flux_y(i,j+1,l)
          enddo
        enddo
      enddo

c
c Compute flux convergence
c
      do l=lmo,1,-1
        do j=j_0s,j_1s
          do i=1,im-1
            conv(i,j,l) = flux_x(i,j,l)-flux_x(i+1,j,l)
     &                   +flux_y(i,j,l)-flux_y(i,j+1,l)
          enddo
          i=im
            conv(i,j,l) = flux_x(i,j,l)-flux_x(1,j,l)
     &                   +flux_y(i,j,l)-flux_y(i,j+1,l)
        enddo
        if(have_north_pole) conv(1,jm,l) = sum(flux_y(:,jm,l))/im
#ifdef OCN_GISS_MESO
        if(have_south_pole) conv(1,1,l) = sum(flux_y(:,1,l))/im
#endif
        if(l.lt.lmo) then
          do j=j_0s,j_1s
          do i=1,im
            conv(i,j,l  ) = conv(i,j,l  ) + flux_z(i,j,l)
            conv(i,j,l+1) = conv(i,j,l+1) - flux_z(i,j,l)
          enddo
          enddo
          if(have_north_pole) then
            conv(1,jm,l  ) = conv(1,jm,l  ) + flux_z(1,jm,l)
            conv(1,jm,l+1) = conv(1,jm,l+1) - flux_z(1,jm,l)
          endif
#ifdef OCN_GISS_MESO
          if(have_south_pole) then
            conv(1,1,l  ) = conv(1,1,l  ) + flux_z(1,1,l)
            conv(1,1,l+1) = conv(1,1,l+1) - flux_z(1,1,l)
          endif
#endif
        endif
      enddo
      if(have_north_pole) then
        do l=1,lmo
          trm(2:im,jm,l) = trm(1,jm,l)
          conv(2:im,jm,l) = conv(1,jm,l)
        enddo
      endif
#ifdef OCN_GISS_MESO
      if(have_south_pole) then
        do l=1,lmo
          trm(2:im,1,l) = trm(1,1,l)
          conv(2:im,1,l) = conv(1,1,l)
        enddo
      endif
#endif


c
c Reduce the magnitude of flux convergences where they cause
c TRM to become negative.  Keep a tally of the convergence
c adjustments in convadj_j(:)
c
      convadj_j(:) = 0d0
      convpos_j(:) = 0d0
      do l=1,lmo
      do j=j_0,j_1
      do i=1,im
        if(lmm(i,j).lt.l) cycle
        if(conv(i,j,l).gt.0.) then
          convpos_j(j) = convpos_j(j) +conv(i,j,l)
        elseif(conv(i,j,l).lt.-trm(i,j,l)) then
          convadj_j(j) = convadj_j(j) -trm(i,j,l)-conv(i,j,l)
          conv(i,j,l) = -trm(i,j,l)
        endif
      enddo
      enddo
      enddo

c
c To conserve the global sum of TRM, subtract the adjustments to
c negative convergences from the positive convergences.
c Use a multiplicative factor rather than straight subtraction.
c
      call globalsum(grid,convpos_j,sumpos,all=.true.)
      call globalsum(grid,convadj_j,sumadj,all=.true.)
      if(sumpos.gt.0.) then
        posadj = 1.-sumadj/sumpos
      else
        posadj = 1.
      endif

c      if(am_i_root() .and. posadj.ne.1.)
c     &     write(6,*) 'posadj ',posadj

c
c Apply the corrected flux convergence
c
      do l=1,lmo
      do j=j_0,j_1
      do i=1,im
        if(lmm(i,j).lt.l) cycle
        if(conv(i,j,l).gt.0.) conv(i,j,l) = conv(i,j,l)*posadj
        trm(i,j,l) = max(0d0,trm(i,j,l)+conv(i,j,l))
      enddo
      enddo
      enddo

      RETURN
      END SUBROUTINE wrapAdjustFluxes

cC****
c      SUBROUTINE limitedValue(RFXT,TRM1,TRM)
c      USE GM_COM, only : eps
c      IMPLICIT NONE
c
c      REAL*8, INTENT(INOUT) :: RFXT
c      REAL*8, INTENT( IN  ) :: TRM1, TRM
c
c      IF ( RFXT >  0.95d0*TRM1.and.ABS( RFXT ) >   eps ) THEN
c        RFXT = 0.95d0*TRM1
c      ELSEIF ( RFXT < -0.95d0*TRM.and.ABS( RFXT ) > eps) THEN
c        RFXT = -0.95d0*TRM
c      END IF
c
c      RETURN
c      END SUBROUTINE limitedValue
cC****
c      SUBROUTINE wrapAdjustFluxes(TRM, flux_x, flux_y, flux_z, GIJL)
c!@wrapper to set up local storage and call adjustAndAddFluxes
c      USE GM_COM, only: grid, IM, JM, LMO, KPL, KPL_GLOB,
c     &                  PACK_DATA, ESMF_BCAST, UNPACK_DATA 
c      IMPLICIT NONE
c      REAL*8, DIMENSION(IM,grid%j_strt_halo:grid%j_stop_halo,LMO),
c     &        INTENT(INOUT) :: TRM, flux_x, flux_y, flux_z
c      REAL*8, DIMENSION(IM,grid%j_strt_halo:grid%j_stop_halo,LMO,3),
c     &        INTENT(INOUT) :: GIJL
c!***  local storage for global variables
c      REAL*8, DIMENSION(IM,JM,LMO) ::
c     &            TRM_GLOB, flux_x_GLOB, flux_y_GLOB, flux_z_GLOB
c      REAL*8, DIMENSION(IM,JM,LMO,3) :: GIJL_GLOB
c
c      CALL PACK_DATA(grid, TRM , TRM_GLOB)
c      CALL PACK_DATA(grid, flux_x , flux_x_GLOB)
c      CALL PACK_DATA(grid, flux_y , flux_y_GLOB)
c      CALL PACK_DATA(grid, flux_z , flux_z_GLOB)
c      CALL PACK_DATA(grid, GIJL , GIJL_GLOB)
c      CALL PACK_DATA(grid, KPL , KPL_GLOB)
c
c      CALL ESMF_BCAST(grid, TRM_GLOB)
c      CALL ESMF_BCAST(grid, flux_x_GLOB)
c      CALL ESMF_BCAST(grid, flux_y_GLOB)
c      CALL ESMF_BCAST(grid, flux_z_GLOB)
c      CALL ESMF_BCAST(grid, GIJL_GLOB)
c      CALL ESMF_BCAST(grid, KPL_GLOB)
c
c      call adjustAndAddFluxes (TRM_GLOB, flux_x_GLOB, flux_y_GLOB,
c     &                         flux_z_GLOB, GIJL_GLOB)
c      call UNPACK_DATA (grid, TRM_GLOB, TRM)
c      call UNPACK_DATA (grid, GIJL_GLOB, GIJL)
c
c      RETURN
c      END SUBROUTINE wrapAdjustFluxes
cC****
c      SUBROUTINE adjustAndAddFluxes (TRM, flux_x, flux_y, flux_z, GIJL)
c!@sum applies limiting process to  GM fluxes for tracer quantities
c!@ver  1.0
c! this routine is serialized to preserve original algorithm for
c! TRM update
c      USE GM_COM, only: IM, JM, LMO, LMM, LMU, LMV, KPL_GLOB
c      IMPLICIT NONE
c      REAL*8, DIMENSION(IM,JM,LMO),
c     &        INTENT(INOUT) :: TRM, flux_x, flux_y, flux_z
c      REAL*8, DIMENSION(IM,JM,LMO,3),
c     &        INTENT(INOUT) :: GIJL
c
c      REAL*8  RFXT, RFYT, RFZT, STRNP
c      INTEGER I,J,L,IM1
c
c!updates TRM and GIJL with limit adjusted flux values
c
c      DO L=1,LMO
cC**** Non-Polar boxes
c      DO J=2,JM-1
c      IM1 = IM
c      DO I=1,IM
c      IF(LMM(I,J).le.0) GO TO 611
c
c!adjustments in x direction
c      IF(LMU(IM1,J).ge.L) THEN
cC**** If fluxes are more than 95% of tracer amount, limit fluxes
c        RFXT = flux_x(I,J,L)
c        call limitedValue(RFXT,TRM(IM1,J,L),TRM(I,J,L))
cC**** Add and Subtract horizontal X flux
c        TRM(I  ,J,L) = TRM(I  ,J,L) + RFXT
c        TRM(IM1,J,L) = TRM(IM1,J,L) - RFXT
cC**** Save Diagnostic, GIJL(1) = RFXT
c        GIJL(I,J,L,1) = GIJL(I,J,L,1) + RFXT
c        !update flux limited values
c        flux_x(I,J,L) =  RFXT
c      END IF
c
c!adjustments in y direction
cC**** Loop for Fluxes in Y-direction
c      IF(LMV(I,J-1).ge.L) THEN
cC**** If fluxes are more than 95% of tracer amount, limit fluxes
c        RFYT = flux_y(I,J,L)
c        call limitedValue(RFYT,TRM(I,J-1,L),TRM(I,J,L))
cC**** Add and Subtract horizontal Y fluxes
c        TRM(I,J  ,L) = TRM(I,J  ,L) + RFYT
c        TRM(I,J-1,L) = TRM(I,J-1,L) - RFYT
cC**** Save Diagnostic, GIJL(2) = RFYT
c        GIJL(I,J,L,2) = GIJL(I,J,L,2) + RFYT
c        !update flux limited values
c        flux_y(I, J, L) = RFYT
c      END IF
c
cC**** END of I and J loops
c  611 IM1 = I
c      END DO
c      END DO
c
cC**** North Polar box
c      STRNP=0.
cC**** Fluxes in Y-direction
c      DO I=1,IM
c      IF(LMV(I,JM-1).ge.L) THEN
c        RFYT = flux_y(I,JM,L)
cC**** If fluxes are more than 95% of tracer amount, limit fluxes
c        call limitedValue(RFYT,TRM(I,JM-1,L),TRM(1,JM,L))
cC**** Add and Subtract horizontal Y fluxes
c        STRNP= STRNP + RFYT
c        TRM(I,JM-1,L) = TRM(I,JM-1,L) - RFYT
c        flux_y(I,JM,L) = RFYT
c      END IF
c      END DO
cC**** adjust polar box
c      TRM(1,JM,L)=TRM(1,JM,L) + STRNP/IM
cC**** Save Diagnostic, GIJL(2) = RFYT
c      GIJL(1,JM,L,2) = GIJL(1,JM,L,2) + STRNP
c
c
c 621  END DO    !L loop
c
cC**** Summation of explicit fluxes to TR, (R(T)) --- Z-direction
cC**** Extracted from preceding L loop to allow parallelization of that
cC**** loop; this loop cannot be
c      DO L=1,LMO
cC**** Non-Polar boxes
c      DO J=2,JM-1
c      DO I=1,IM
c      IF(LMM(I,J).le.0) GO TO 710
cC**** Loop for Fluxes in Z-direction
c      IF(KPL_GLOB(I,J).gt.L) GO TO 710
c      IF(LMM(I,J).lt.L) GO TO 710
c      IF(LMM(I,J).gt.L) THEN
cC****   tracer/salinity/enthalpy
c        RFZT = flux_z(I,J,L)
cC**** If fluxes are more than 95% of tracer amount, limit fluxes
c        call limitedValue(RFZT,TRM(I,J,L+1),TRM(I,J,L))
cC**** Add and Subtract vertical flux. Note +ve upward flux
c        TRM(I,J,L  ) = TRM(I,J,L  ) + RFZT
c        TRM(I,J,L+1) = TRM(I,J,L+1) - RFZT
cC**** Save Diagnostic, GIJL(3) = RFZT
c        GIJL(I,J,L,3) = GIJL(I,J,L,3) + RFZT
c      END IF
c
cC**** END of I and J loops
c 710  CONTINUE
c      END DO
c      END DO
c
c
cC****   North Polar box
cC****   Loop for Fluxes in Z-direction
c        IF(KPL_GLOB(1,JM).gt.L) GO TO 720
c        IF(LMM(1,JM).lt.L) GO TO 720
c        IF(LMM(1,JM).gt.L) THEN
cC****     Calculate new tracer/salinity/enthalpy
c            RFZT = flux_z(1,JM,L)
cC****     If fluxes are more than 95% of tracer amount, limit fluxes
c          call limitedValue(RFZT,TRM(1,JM,L+1),TRM(1,JM,L))
cC****     Add and Subtract vertical flux. Note +ve upward flux
c          TRM(1,JM,L  ) = TRM(1,JM,L  ) + RFZT
c          TRM(1,JM,L+1) = TRM(1,JM,L+1) - RFZT
cC****     Save Diagnostic, GIJL(3) = RFZT
c          GIJL(1,JM,L,3) = GIJL(1,JM,L,3) + RFZT
c        END IF
c 720    CONTINUE
c
c      END DO
c
c      RETURN
c      END SUBROUTINE adjustAndAddFluxes
C****
      SUBROUTINE addFluxes (TRM, flux_x, flux_y, flux_z, GIJL)
!@sum applies GM fluxes to tracer quantities
!@ver  1.0
      USE GM_COM, only: grid, GET, IM, JM, LMO
#ifdef OCN_GISS_MESO
     &                       ,DXYPO,MO
#endif
     &                       ,LMM, LMU, LMV, KPL

      IMPLICIT NONE
      REAL*8, DIMENSION(IM,grid%j_strt_halo:grid%j_stop_halo,LMO),
     &         INTENT(INOUT) :: TRM
      REAL*8, DIMENSION(IM,grid%j_strt_halo:grid%j_stop_halo,LMO),
     &          INTENT(IN) :: flux_x, flux_y, flux_z
      REAL*8, DIMENSION(IM,grid%j_strt_halo:grid%j_stop_halo,LMO,3),
     &          INTENT(INOUT) :: GIJL
      REAL*8 STRNP, RFZT
#ifdef OCN_GISS_MESO
     &      ,STRSP
#endif
      INTEGER :: I, J, L, IM1

      INTEGER :: J_0, J_1, J_0S, J_1S, J_0STG, J_1STG, J_0H, J_1H
      LOGICAL :: HAVE_NORTH_POLE
#ifdef OCN_GISS_MESO
      LOGICAL :: HAVE_SOUTH_POLE
#endif

c**** Extract domain decomposition info
      CALL GET(grid, J_STRT = J_0, J_STOP = J_1,
     &               J_STRT_SKP  = J_0S,   J_STOP_SKP  = J_1S,
     &               J_STRT_STGR=J_0STG, J_STOP_STGR=J_1STG,
     &               J_STRT_HALO = J_0H,   J_STOP_HALO = J_1H,
     &               HAVE_NORTH_POLE = HAVE_NORTH_POLE
#ifndef OCN_GISS_MESO
     &               )
#else
     &              ,HAVE_SOUTH_POLE = HAVE_SOUTH_POLE)
#endif

!updates TRM and GIJL using flux values
      DO L=1,LMO
C**** Non-Polar boxes
      DO J=J_0S,J_1S
      IM1 = IM
      DO I=1,IM

      IF(LMM(I,J).le.0) GO TO 613

        IF(LMU(IM1,J).ge.L) THEN
          TRM(I  ,J,L) = TRM(I  ,J,L) + flux_x(I,J,L)
          TRM(IM1,J,L) = TRM(IM1,J,L) - flux_x(I,J,L)
#ifndef OCN_GISS_MESO
C**** Save Diagnostic, GIJL(1) = flux_x(I,J,L)
          GIJL(I,J,L,1) = GIJL(I,J,L,1) + flux_x(I,J,L)
#endif
        END IF

        IF(LMV(I,J-1).ge.L) THEN
          TRM(I,J  ,L) = TRM(I,J  ,L) + flux_y(I,J,L)
          TRM(I,J-1,L) = TRM(I,J-1,L) - flux_y(I,J,L)
#ifndef OCN_GISS_MESO
C**** Save Diagnostic, GIJL(2) = flux_y(I,J,L)
          GIJL(I,J,L,2) = GIJL(I,J,L,2) + flux_y(I,J,L)
#endif
        END IF

C**** END of I and J loops
  613 IM1 = I
      END DO
      END DO

      IF(HAVE_NORTH_POLE) THEN
C****   North Polar box
        STRNP=0.
C****   Fluxes in Y-direction
        DO I=1,IM
          IF(LMV(I,JM-1).ge.L) THEN
C****       Add and Subtract horizontal Y fluxes
            STRNP= STRNP + flux_y(I,JM,L)
            TRM(I,JM-1,L) = TRM(I,JM-1,L) - flux_y(I,JM,L)
          END IF
        END DO
C**** adjust polar box
        TRM(1,JM,L)=TRM(1,JM,L) + STRNP/IM
#ifndef OCN_GISS_MESO
C**** Save Diagnostic, GIJL(2) = STRNP
        GIJL(1,JM,L,2) = GIJL(1,JM,L,2) + STRNP
#endif
      END IF

#ifdef OCN_GISS_MESO
      IF(HAVE_SOUTH_POLE) THEN
C****   South Polar box
        STRSP=0.
C****   Fluxes in Y-direction
        DO I=1,IM
          IF(LMV(I,2).ge.L) THEN
C****       Add and Subtract horizontal Y fluxes
            STRSP= STRSP + flux_y(I,1,L)
            TRM(I,2,L) = TRM(I,2,L) - flux_y(I,1,L)
          END IF
        END DO
C**** adjust polar box
        TRM(1,1,L)=TRM(1,1,L) + STRSP/IM
      END IF
#endif

 622  END DO    !L loop


      J = J_1S + 1
      if (J.lt.JM) then

      DO L=1,LMO
C**** Non-Polar boxes
      IM1 = IM
      DO I=1,IM

      IF(LMM(I,J).le.0) GO TO 614

        IF(LMV(I,J-1).ge.L) THEN
          TRM(I,J-1,L) = TRM(I,J-1,L) - flux_y(I,J,L)
        END IF

C**** END of L and I loops
  614 IM1 = I
      END DO
      END DO

      endif

C**** Summation of explicit fluxes to TR, (R(T)) --- Z-direction
C**** Extracted from preceding L loop to allow parallelization of that
C**** loop; this loop can't be
      DO L=1,LMO
C**** Non-Polar boxes
      DO J=J_0S,J_1S
      DO I=1,IM
      IF(LMM(I,J).le.0) GO TO 710
C**** Loop for Fluxes in Z-direction
#ifndef OCN_GISS_MESO
      IF(KPL(I,J).gt.L) GO TO 710
#endif
      IF(LMM(I,J).lt.L) GO TO 710
      IF(LMM(I,J).gt.L) THEN
C****   tracer/salinity/enthalpy
        RFZT = flux_z(I,J,L)
C**** Add and Subtract vertical flux. Note +ve upward flux
        TRM(I,J,L  ) = TRM(I,J,L  ) + RFZT
        TRM(I,J,L+1) = TRM(I,J,L+1) - RFZT
#ifndef OCN_GISS_MESO
C**** Save Diagnostic, GIJL(3) = RFZT
        GIJL(I,J,L,3) = GIJL(I,J,L,3) + RFZT
#endif
      END IF

C**** END of I and J loops
 710  CONTINUE
      END DO
      END DO


      IF(HAVE_NORTH_POLE) THEN
C**** North Polar box
C**** Loop for Fluxes in Z-direction
#ifndef OCN_GISS_MESO
        IF(KPL(1,JM).gt.L) GO TO 720
#endif
        IF(LMM(1,JM).lt.L) GO TO 720
        IF(LMM(1,JM).gt.L) THEN
C**** Calculate new tracer/salinity/enthalpy
          RFZT = flux_z(1,JM,L)
C**** Add and Subtract vertical flux. Note +ve upward flux
          TRM(1,JM,L  ) = TRM(1,JM,L  ) + RFZT
          TRM(1,JM,L+1) = TRM(1,JM,L+1) - RFZT
#ifndef OCN_GISS_MESO
C**** Save Diagnostic, GIJL(3) = RFZT
          GIJL(1,JM,L,3) = GIJL(1,JM,L,3) + RFZT
#endif
        END IF
 720    CONTINUE
      END IF     ! HAVE_NORTH_POLE
#ifndef OCN_GISS_MESO
      END DO
#endif


#ifdef OCN_GISS_MESO
      IF(HAVE_SOUTH_POLE) THEN
C**** South Polar box
C**** Loop for Fluxes in Z-direction
        IF(LMM(1,1).lt.L) GO TO 730
        IF(LMM(1,1).gt.L) THEN
C**** Calculate new tracer/salinity/enthalpy
          RFZT = flux_z(1,1,L)
C**** Add and Subtract vertical flux. Note +ve upward flux
          TRM(1,1,L  ) = TRM(1,1,L  ) + RFZT
          TRM(1,1,L+1) = TRM(1,1,L+1) - RFZT
        END IF
 730    CONTINUE
      END IF     ! HAVE_SOUTH_POLE
      END DO

c
c Save flux diagnostics.
c Shift x/y indices by 1 for consistency with ocean dynamics
c conventions.
c
      do l=1,lmo
        do j=j_0,j_1
          do i=1,im-1
            gijl(i,j,l,1) = gijl(i,j,l,1) + flux_x(i+1,j,l)
          enddo
          i=im
            gijl(i,j,l,1) = gijl(i,j,l,1) + flux_x(  1,j,l)
          do i=1,im
            gijl(i,j,l,3) = gijl(i,j,l,3) + flux_z(i,j,l)
          enddo
        enddo
        do j=j_0,min(j_1,jm-1)
          do i=1,im
            gijl(i,j,l,2) = gijl(i,j,l,2) + flux_y(i,j+1,l)
          enddo
        enddo
      enddo
#endif

      RETURN
      END SUBROUTINE addFluxes
C****
      SUBROUTINE ISOSLOPE4
!@sum  ISOSLOPE4 calculates the isopycnal slopes from density triads
!@auth Gavin Schmidt/Dan Collins
!@ver  2009/02/13
      USE GM_COM
      IMPLICIT NONE
      INTEGER I,J,L,IM1
      REAL*8 :: AIX0ST,AIX2ST,AIY0ST,AIY2ST,SIX0,SIX2,SIY0,SIY2,
     *          AIX1ST,AIX3ST,AIY1ST,AIY3ST,SIX1,SIX3,SIY1,SIY3
      Real*8 :: byAIDT, DSX0sq,DSX1sq,DSX2sq,DSX3sq,
     *                  DSY0sq,DSY1sq,DSY2sq,DSY3sq
      INTEGER :: J_0, J_1, J_0S, J_1S, J_0STG, J_1STG, J_0H, J_1H
#ifndef OCN_GISS_MESO
      LOGICAL :: HAVE_NORTH_POLE
#endif
      LOGICAL :: dothis

c**** Extract domain decomposition info
      CALL GET(grid, J_STRT = J_0, J_STOP = J_1,
     &               J_STRT_SKP  = J_0S,   J_STOP_SKP  = J_1S,
     &               J_STRT_STGR=J_0STG, J_STOP_STGR=J_1STG,
     &               J_STRT_HALO = J_0H,   J_STOP_HALO = J_1H
#ifdef OCN_GISS_MESO
     &               )
#else
     &              ,HAVE_NORTH_POLE = HAVE_NORTH_POLE)
#endif

C**** Calculate horizontal and vertical density gradients.
      CALL DENSGRAD
C****
C**** Three grid boxes (triads) are used for each slope. The GM skew
C**** diffusion coefficient is calculated for each triad as well.
C**** The diffusion coefficient is taken from Visbeck et al (1997)
C****
C**** Main Loop over I,J and L
!$OMP PARALLEL DO  PRIVATE(L,J,IM1,I, byAIDT,
!$omp&  AIX0ST,AIX1ST,AIX2ST,AIX3ST, AIY0ST,AIY1ST,AIY2ST,AIY3ST,
!$OMP&  SIX0  ,SIX1  ,SIX2  ,SIX3  , SIY0  ,SIY1  ,SIY2  ,SIY3  ,
!$OMP&  DSX0sq,DSX1sq,DSX2sq,DSX3sq, DSY0sq,DSY1sq,DSY2sq,DSY3sq)
      DO L=1,LMO
      DO J=J_0STG,J_1STG
      IM1=IM
      DO I=1,IM
      IF(LMM(I,J).lt.L) GO TO 800
#ifndef OCN_GISS_MESO
      AIX0(I,J,L) = AINV(I,J)
      AIX2(I,J,L) = AINV(I,J)
      AIY0(I,J,L) = AINV(I,J)
      AIY2(I,J,L) = AINV(I,J)
      AIX1(I,J,L) = AINV(I,J)
      AIX3(I,J,L) = AINV(I,J)
      AIY1(I,J,L) = AINV(I,J)
      AIY3(I,J,L) = AINV(I,J)
#endif
C**** SIX0, SIY0, SIX2, SIY2: four slopes that use RHOMZ(L)
      IF(L.EQ.LMM(I,J) .or. RHOMZ(I,J,L).eq.0.) THEN
        AIX0ST = 0.
        AIX2ST = 0.
        AIY0ST = 0.
        AIY2ST = 0.
        SIX0 = 0.
        SIX2 = 0.
        SIY0 = 0.
        SIY2 = 0.
#ifdef OCN_GISS_MESO
        AIX0(I,J,L) = 0.
        AIX2(I,J,L) = 0.
        AIY0(I,J,L) = 0.
        AIY2(I,J,L) = 0.
#endif
      ELSE
#ifdef OCN_GISS_MESO
        AIX0ST = ARIV(I,J,L)
        AIX2ST = ARIV(I,J,L)
        AIY0ST = ARIV(I,J,L)
        AIY2ST = ARIV(I,J,L)
#else
        AIX0ST = AINV(I,J)
        AIX2ST = AINV(I,J)
        AIY0ST = AINV(I,J)
        AIY2ST = AINV(I,J)
#endif
        SIX0 = RHOX(I  ,J,L) * BYRHOZ(I,J,L)
        SIX2 = RHOX(IM1,J,L) * BYRHOZ(I,J,L)
        SIY2 = RHOY(I,J-1,L) * BYRHOZ(I,J,L)
        SIY0 = RHOY(I,J  ,L) * BYRHOZ(I,J,L)
#ifdef OCN_GISS_MESO
        IF (ARIV(I,J,L).gt.0.) THEN ! limit slopes <ML
          byAIDT = 1 / (4*DTS*(AINV(I,J,L)+ARIV(I,J,L)))
#else
        IF (L.ge.KPL(I,J).and.AINV(I,J).gt.0.) THEN ! limit slopes <ML
          byAIDT = 1 / (4*DTS*(AINV(I,J)+ARIV(I,J)))
#endif
          DSX0sq = DZV(I,J,L)**2 * byAIDT
          DSX2sq = DZV(I,J,L)**2 * byAIDT
          DSY0sq = DZV(I,J,L)**2 * byAIDT
          DSY2sq = DZV(I,J,L)**2 * byAIDT
          If (SIX0**2 > DSX0sq)  AIX0ST = AIX0ST * DSX0sq / SIX0**2
          If (SIX2**2 > DSX2sq)  AIX2ST = AIX2ST * DSX2sq / SIX2**2
          If (SIY0**2 > DSY0sq)  AIY0ST = AIY0ST * DSY0sq / SIY0**2
          If (SIY2**2 > DSY2sq)  AIY2ST = AIY2ST * DSY2sq / SIY2**2
        END IF
#ifdef OCN_GISS_MESO
C**** AI are always * layer thickness for vertical gradient in FXX, FYY
        AIX0(I,J,L) = AIX0ST * DZV(I,J,L) * BYDH(I,J,L)
        AIX2(I,J,L) = AIX2ST * DZV(I,J,L) * BYDH(I,J,L)
        AIY0(I,J,L) = AIY0ST * DZV(I,J,L) * BYDH(I,J,L)
        AIY2(I,J,L) = AIY2ST * DZV(I,J,L) * BYDH(I,J,L)
#endif
      ENDIF
C**** SIX1, SIY1, SIX3, SIY3: four slopes that use RHOMZ(L-1)
      IF(L.EQ.1.) THEN
        AIX1ST = 0. ; AIX3ST = 0. ; AIY1ST = 0. ;  AIY3ST = 0.
        SIX1 = 0.   ; SIX3 = 0.   ; SIY1 = 0.   ;  SIY3 = 0.
#ifdef OCN_GISS_MESO
        AIX1(I,J,L) = 0. ; AIX3(I,J,L) = 0.
        AIY1(I,J,L) = 0. ; AIY3(I,J,L) = 0.
#endif
      ELSEIF (RHOMZ(I,J,L-1).eq.0.) THEN
        AIX1ST = 0. ; AIX3ST = 0. ; AIY1ST = 0. ;  AIY3ST = 0.
        SIX1 = 0.   ; SIX3 = 0.   ; SIY1 = 0.   ;  SIY3 = 0.
#ifdef OCN_GISS_MESO
        AIX1(I,J,L) = 0. ; AIX3(I,J,L) = 0.
        AIY1(I,J,L) = 0. ; AIY3(I,J,L) = 0.
#endif
      ELSE
#ifdef OCN_GISS_MESO
        AIX1ST = ARIV(I,J,L)
        AIX3ST = ARIV(I,J,L)
        AIY1ST = ARIV(I,J,L)
        AIY3ST = ARIV(I,J,L)
#else
        AIX1ST = AINV(I,J)
        AIX3ST = AINV(I,J)
        AIY1ST = AINV(I,J)
        AIY3ST = AINV(I,J)
#endif
        SIX1 = RHOX(I  ,J,L) * BYRHOZ(I,J,L-1)
        SIX3 = RHOX(IM1,J,L) * BYRHOZ(I,J,L-1)
        SIY1 = RHOY(I,J  ,L) * BYRHOZ(I,J,L-1)
        SIY3 = RHOY(I,J-1,L) * BYRHOZ(I,J,L-1)
#ifdef OCN_GISS_MESO
        IF (ARIV(I,J,L).gt.0.) THEN ! limit slopes <ML
          byAIDT = 1 / (4*DTS*(AINV(I,J,L)+ARIV(I,J,L)))
#else
        IF (L.ge.KPL(I,J)+1.and.AINV(I,J).gt.0.) THEN ! limit slopes <ML
          byAIDT = 1 / (4*DTS*(AINV(I,J)+ARIV(I,J)))
#endif
          DSX1sq = DZV(I,J,L-1)**2 * byAIDT
          DSX3sq = DZV(I,J,L-1)**2 * byAIDT
          DSY1sq = DZV(I,J,L-1)**2 * byAIDT
          DSY3sq = DZV(I,J,L-1)**2 * byAIDT
          If (SIX1**2 > DSX1sq)  AIX1ST = AIX1ST * DSX1sq / SIX1**2
          If (SIX3**2 > DSX3sq)  AIX3ST = AIX3ST * DSX3sq / SIX3**2
          If (SIY1**2 > DSY1sq)  AIY1ST = AIY1ST * DSY1sq / SIY1**2
          If (SIY3**2 > DSY3sq)  AIY3ST = AIY3ST * DSY3sq / SIY3**2
        END IF
#ifdef OCN_GISS_MESO
C**** AI are always * layer thickness for vertical gradient in FXX, FYY
        AIX1(I,J,L) = AIX1ST * DZV(I,J,L-1) * BYDH(I,J,L)
        AIX3(I,J,L) = AIX3ST * DZV(I,J,L-1) * BYDH(I,J,L)
        AIY1(I,J,L) = AIY1ST * DZV(I,J,L-1) * BYDH(I,J,L)
        AIY3(I,J,L) = AIY3ST * DZV(I,J,L-1) * BYDH(I,J,L)
#endif
      ENDIF
C**** AIX0...AIX3, AIY0...AIY3
      ASX0(I,J,L) = AIX0ST * SIX0
      ASX1(I,J,L) = AIX1ST * SIX1
      ASX2(I,J,L) = AIX2ST * SIX2
      ASX3(I,J,L) = AIX3ST * SIX3
      ASY0(I,J,L) = AIY0ST * SIY0
      ASY1(I,J,L) = AIY1ST * SIY1
      ASY2(I,J,L) = AIY2ST * SIY2
      ASY3(I,J,L) = AIY3ST * SIY3
C**** S2X0...S2X3, S2Y0...S2Y3
      S2X0(I,J,L) = AIX0ST * SIX0 * SIX0
      S2X1(I,J,L) = AIX1ST * SIX1 * SIX1
      S2X2(I,J,L) = AIX2ST * SIX2 * SIX2
      S2X3(I,J,L) = AIX3ST * SIX3 * SIX3
      S2Y0(I,J,L) = AIY0ST * SIY0 * SIY0 * BYDYP(J) * DYVO(J)
      S2Y1(I,J,L) = AIY1ST * SIY1 * SIY1 * BYDYP(J) * DYVO(J)
      S2Y2(I,J,L) = AIY2ST * SIY2 * SIY2 * BYDYP(J) * DYVO(J-1)
      S2Y3(I,J,L) = AIY3ST * SIY3 * SIY3 * BYDYP(J) * DYVO(J-1)
#ifdef OCN_GISS_MESO
      AIX0(I,J,L) = ARIV(I,J,L)
      AIX2(I,J,L) = ARIV(I,J,L)
      AIY0(I,J,L) = ARIV(I,J,L)
      AIY2(I,J,L) = ARIV(I,J,L)
      AIX1(I,J,L) = ARIV(I,J,L)
      AIX3(I,J,L) = ARIV(I,J,L)
      AIY1(I,J,L) = ARIV(I,J,L)
      AIY3(I,J,L) = ARIV(I,J,L)
#endif
  800 IM1 = I
      END DO
      END DO
      END DO
!$OMP END PARALLEL DO
      RETURN
      END SUBROUTINE ISOSLOPE4
C****
      SUBROUTINE DENSGRAD
!@sum  DENSGRAD calculates all horizontal and vertical density gradients
!@auth Gavin Schmidt/Dan Collins
!@ver  1.0
      USE GM_COM
      USE FILEMANAGER
#ifdef OCN_GISS_MESO
      use odiag, only : oijl=>oijl_loc,ijl_mfvb
      use ocean, only : dzo,dxvo
      use ocean, only : nbyzm,i1yzm,i2yzm,lmm
      use ocean, only : nbyzu,i1yzu,i2yzu,lmu
      use ocean, only : nbyzv,i1yzv,i2yzv
      use constant, only : rhows
#else
      USE ODIAG, only : oij=>oij_loc,ij_gmsc
#endif

      IMPLICIT NONE

      !REAL*8, DIMENSION(IM,JM,LMO) :: RHO
      REAL*8, DIMENSION(IM,grid%j_strt_halo:grid%j_stop_halo,LMO) :: RHO
#ifdef OCN_GISS_MESO
     *  ,KBYRHOZ,kappam3d
      REAL*8, DIMENSION(IM,grid%j_strt_halo:grid%j_stop_halo,0:LMO) ::
     &     PSIY
#endif

      REAL*8  BYRHO,DZVLM1,CORI,BETA,ARHO,ARHOX,ARHOY,ARHOZ,AN,RD
     *     ,BYTEADY,DH0,DZSUMX,DZSUMY,R1,R2,P12
      REAL*8, SAVE :: HUP
      INTEGER I,J,L,IM1,LAV,iu_ODIFF
      INTEGER, SAVE :: IFIRST = 1
      CHARACTER TITLE*80
      Real*8,External   :: VOLGSP
      REAL*8, DIMENSION(IM,JM) ::  AINV_glob

      INTEGER :: J_0, J_1, J_0S, J_1S, J_0STG, J_1STG, J_0H, J_1H
      INTEGER :: J_1HR
      LOGICAL :: HAVE_NORTH_POLE
      logical :: dothis
#ifdef OCN_GISS_MESO
      LOGICAL :: HAVE_SOUTH_POLE
      real*8 :: k1d(lmo)
      real*8 :: wtup,wtdn,rhomid,rhoy_,sly,krat,kmax
      integer :: n
#endif

#ifdef SIMPLE_MESODIFF
!@var snavg vertical average of the product of isoneutral slope times
!@+         N (sqrt of B-V freq).  snavg = 1/t_eady
      real*8 :: snavg,snsum,delz,ksurf,arhoh
      real*8, parameter ::
!@param maxslope maximum slope used in t_eady calculation
     &         maxslope=4d-4
!@param mindepth (m) minimum averaging depth in t_eady calculation.  If
!@+              the local depth is less than mindepth, the calculation
!@+              is performed as if the depth were mindepth and isoneutral
!@+              slopes were zero between the local depth and mindepth.
     &        ,mindepth=400d0
!@param kfac the value of AMU used for SIMPLE_MESODIFF
     &        ,kfac=12d-3
!@param lscale (m) fixed length scale for calculating mesoscale diffusivity
     &        ,lscale=2.5d5
!@param mink,maxk (m2/s) lower, upper bounds for mesoscale diffusivity
     &        ,mink=100d0,maxk=6000d0
!@param minsinlat minimum for 1/sin(lat) scaling of mesoscale diffusivity
     &        ,minsinlat=.1d0 ! roughly sin(6 degrees)
#endif

c**** Extract domain decomposition info
      CALL GET(grid, J_STRT = J_0, J_STOP = J_1,
     &               J_STRT_SKP  = J_0S,   J_STOP_SKP  = J_1S,
     &               J_STRT_STGR=J_0STG, J_STOP_STGR=J_1STG,
     &               J_STRT_HALO = J_0H,   J_STOP_HALO = J_1H,
     &               HAVE_NORTH_POLE = HAVE_NORTH_POLE
#ifndef OCN_GISS_MESO
     &               )
#else
     &              ,HAVE_SOUTH_POLE = HAVE_SOUTH_POLE)
#endif

C**** set up geometry needed
      IF (IFIRST.eq.1) THEN
        !DO J=1,JM
        DO J=J_0,J_1
          BYDYP(J)=1d0/DYPO(J)
          BYDXP(J)=1d0/DXPO(J)
        END DO
        !DO J=1,JM-1
        DO J=J_0,J_1S
          BYDYV(J)=1d0/DYVO(J)
        END DO
C**** Calculate level at 1km depth
        LUP=0
   10   LUP=LUP + 1
        IF (ZE(LUP+1).lt.1d3) GOTO 10
        HUP = ZE(LUP)
#ifdef OCN_GISS_MESO
        IFIRST = 0.
#endif

        CALL HALO_UPDATE(grid,BYDYP (grid%j_strt_halo:grid%j_stop_halo),
     *                 FROM=SOUTH+NORTH)
        CALL HALO_UPDATE(grid,BYDYV (grid%j_strt_halo:grid%j_stop_halo),
     *                 FROM=SOUTH+NORTH)
        CALL HALO_UPDATE(grid,BYDXP (grid%j_strt_halo:grid%j_stop_halo),
     *                 FROM=SOUTH+NORTH)

      END IF
C****
      call vbar_gm0

      !initialize
      RHO = -0.0; BYDH = -0.0;
      DZV = -0.0; BYDZV = -0.0;

!$OMP PARALLEL DO  PRIVATE(L,J,I,BYRHO,DH0,DZVLM1)
      DO L=1,LMO
      !DO J=1,JM
      DO J=J_0,J_1
      DO I=1,IM

#ifndef OCN_GISS_MESO
C**** Bottom layer LMM(I,J)
        RHOMZ(I,J,L) = 0d0
C**** Initialize RHOX and RHOY as 0.
        RHOX(I,J,L) = 0d0
        RHOY(I,J,L) = 0d0
#endif
C**** Skip non-ocean grid points
        IF (L.le.LMM(I,J)) THEN
C**** RHO(I,J,L)  Density=1/specific volume
          BYRHO = VBAR(I,J,L)
            DH0 =   DH(I,J,L)
#ifdef OCN_GISS_MESO
          IF (J.eq.1)  THEN         ! South pole
            BYRHO = VBAR(1,1,L)
              DH0 =   DH(1,1,L)
          END IF
#endif
          IF (J.eq.JM) THEN         ! North pole
            BYRHO = VBAR(1,JM,L)
              DH0 =   DH(1,JM,L)
          END IF
          RHO(I,J,L)  = 1d0/BYRHO
          BYDH(I,J,L) = 1d0/DH0
          IF(L.gt.1) THEN
            DZVLM1     = 0.5* (DH(I,J,L) + DH(I,J,L-1))
            DZV(I,J,L-1) = DZVLM1
            BYDZV(I,J,L-1) = 1d0/DZVLM1
          END IF
        END IF
 911  END DO
      END DO
      END DO
!$OMP  END PARALLEL DO

      CALL HALO_UPDATE(grid,RHO (:,grid%j_strt_halo:grid%j_stop_halo,:),
     *                 FROM=SOUTH+NORTH)

      CALL HALO_UPDATE(grid,G3D)
      CALL HALO_UPDATE(grid,S3D)
      CALL HALO_UPDATE(grid,P3D)

C**** Calculate density gradients

      !initialize
      RHOMZ = -0.0; BYRHOZ = -0.0;
      RHOY = -0.0; RHOX = -0.0;

      J_1HR = min(J_1STG+1, JM)
!$OMP PARALLEL DO  PRIVATE(L,J,IM1,I)
      DO L=1,LMO
       !DO J=2,JM
        DO J=J_0STG,J_1HR
          IM1 = IM
          DO I=1,IM
C**** Skip non-ocean grid points
            IF(LMM(I,J).lt.L) GO TO 931
C**** minus vertical gradient
            IF(L.gt.1) THEN
              RHOMZ(I,J,L-1)=MAX(0d0,
     *           -dVBARdZ(I,J,L)*BYDZV(I,J,L-1)/VBAR(I,J,L-1)**2)
              IF(RHOMZ(I,J,L-1).ne.0.)
     *             BYRHOZ(I,J,L-1)=1./RHOMZ(I,J,L-1)
            END IF
C**** Calculate horizontal gradients
            IF(LMV(I,J-1).ge.L) THEN
              p12 = .5d0*(p3d(i,j-1,l)+p3d(i,j,l))
              r1 = 1d0/volgsp(g3d(i,j-1,l),s3d(i,j-1,l),p12)
              r2 = 1d0/volgsp(g3d(i,j  ,l),s3d(i,j  ,l),p12)
              RHOY(I,J-1,L) =
C     *           (RHO(I,J,L) - RHO(I,J-1,L))*BYDYV(J-1)
     *           (r2 - r1)*BYDYV(J-1)
            ENDIF
            IF(LMU(IM1,J).ge.L) THEN
              p12 = .5d0*(p3d(im1,j,l)+p3d(i,j,l))
              r1 = 1d0/volgsp(g3d(im1,j,l),s3d(im1,j,l),p12)
              r2 = 1d0/volgsp(g3d(i  ,j,l),s3d(i  ,j,l),p12)
              RHOX(IM1,J,L) =
C     *           (RHO(I,J,L) - RHO(IM1,J,L))*BYDXP(J)
     *           (r2 - r1)*BYDXP(J)
            ENDIF
  931       IM1 = I
          END DO
        END DO
      END DO
#ifdef OCN_GISS_MESO 
      CALL OCN_mesosc(kappam3d)
#endif
!$OMP  END PARALLEL DO
      AINV = 0.
      ARIV = 0.
#ifdef SIMPLE_MESODIFF
! A prescription of mesoscale K a la "Visbeck" that is being used for
! a GISS-MIT MIP (excepting the latitude dependence).
! It is like the default method in that
! K = factor*L*L*[upper-ocean average of 1/t_eady=N*isopycnal_slope]
! but with the following differences:
! 1. The length scale L is a constant rather than the deformation
!    radius R_d=NH/f (or beta-plane R_d=sqrt(NH/beta) at low latitudes)
! 2. Off the equator, the explicit latitude dependence is 1/sin(lat)
!    rather than the 1/sin**2(lat) arising from use of L=R_d.  Near the
!    equator, the 1/sin(lat) dependence is capped, while the use of
!    beta-plane R_d in the default method makes K independent of N.
! 3. As N goes to zero, K goes to zero via imposition of a maximum
!    isopycnal slope rather than via R_d.  As N goes to infinity
!    but 1/t_eady goes to zero, K goes to zero, while the N**2 in
!    R_d in the default method permits large K (off the equator)
! 4. N*isopycnal_slope is computed locally and averaged over the upper
!    ocean rather than computed using averages of density gradients.
!    Where the ocean is shallower than an imposed limit, the
!    upper-ocean average of 1/t_eady is computed as per the
!    comments above for the mindepth parameter.
! 5. Minimum and maximum K values are imposed.
! 
      do j=j_0s,j_1s
      do i=1,im
        if(lmm(i,j).eq.0) cycle
        im1 = i-1
        if(i.eq.1) im1=im
        arho = 0.
        lav = min(lup,lmm(i,j))
        snsum = 0.
        do l=1,lav
          arho  = arho  + rho(i,j,l)
          arhoz = .5d0*(rhomz(i,j,max(1,l-1))
     &                + rhomz(i,j,min(l,lmm(i,j)-1)) )
          if(arhoz.le.0.) cycle
          arhox = .5d0*(rhox(im1,j  ,l)+rhox(i,j,l))
          arhoy = .5d0*(rhoy(i  ,j-1,l)+rhoy(i,j,l))
          arhoh = sqrt(arhox*arhox+arhoy*arhoy)
          snsum = snsum + dzo(l)*sqrt(arhoz)*min(maxslope, arhoh/arhoz)
        enddo
        arho  = arho / real(lav,kind=8)
        delz = max(mindepth, .5*(ze(lav)+ze(lav-1)-ze(1)))
        snavg = snsum * sqrt(grav/arho) / delz
        ksurf = (kfac * (lscale**2)) *snavg
        ksurf = ksurf/max(abs(sinpo(j)),minsinlat)
        ksurf = max(mink, min(maxk, ksurf))
#ifdef OCN_GISS_MESO
        ainv(i,j,l) = ksurf
#else
        ainv(i,j) = ksurf
#endif
      enddo
      enddo
#else
C**** Calculate VMHS diffusion = amu* min(NH/f,equ.rad)^2 /Teady
!$OMP PARALLEL DO  PRIVATE(J,CORI,BETA,IM1,I,ARHO,ARHOX,ARHOY,
!$OMP&  DZSUMX,DZSUMY,LAV,L,ARHOZ,AN,RD,BYTEADY)
      !DO J=2,JM-1
      DO J=J_0S,J_1S
        CORI = ABS(2d0*OMEGA*SINPO(J))
        BETA = ABS(2d0*OMEGA*COSPO(J)/RADIUS)
        IM1=IM
        DO I=1,IM
          IF (LMM(I,J).gt.0) THEN
C**** Calculate average density + gradients over [1,LUP]
            ARHO  = 0.
            ARHOX = 0.
            ARHOY = 0.
            DZSUMX = 0.
            DZSUMY = 0.
            LAV = MIN(LUP,LMM(I,J))
            DO L=1,LAV
              ARHO  = ARHO  + RHO(I,J,L)
              IF(LMU(IM1,J).ge.L) THEN
                ARHOX = ARHOX + RHOX(IM1,J,L)*DZO(L)
                DZSUMX = DZSUMX + DZO(L)
              END IF
              IF(LMU(I  ,J).ge.L) THEN
                ARHOX = ARHOX + RHOX(I  ,J,L)*DZO(L)
                DZSUMX = DZSUMX + DZO(L)
              END IF
              IF(LMV(I,J-1).ge.L) THEN
                ARHOY = ARHOY + RHOY(I,J-1,L)*DZO(L)
                DZSUMY = DZSUMY + DZO(L)
              END IF
              IF(LMV(I,J  ).ge.L) THEN
                ARHOY = ARHOY + RHOY(I,J  ,L)*DZO(L)
                DZSUMY = DZSUMY + DZO(L)
              END IF
            END DO
            ARHO  = ARHO / REAL(LAV,KIND=8)
            IF (DZSUMX.gt.0.) ARHOX = ARHOX / DZSUMX
            IF (DZSUMY.gt.0.) ARHOY = ARHOY / DZSUMY
            IF (LAV.gt.1) THEN
#ifdef OCN_GISS_MESO
              ARHOZ=2*RHOZ1K(I,J)/(ZE(LAV)+ZE(LAV-1)-ZE(1))
#else
              ARHOZ = 2.*(RHO(I,J,LAV)-RHO(I,J,1))/
     *             (ZE(LAV)+ZE(LAV-1)-ZE(1))
#endif
            ELSE
              ARHOZ = 0.
            END IF
C**** avoid occasional inversions. IF ARHOZ<=0 then GM is pure vertical
C**** so keep at zero, and let KPP do the work.
            IF (ARHOZ.gt.0) THEN
#ifdef USE_1D_DIFFUSIVITY
              ARIV(I,J,:) = k1d(:)
#else
#ifdef OCN_GISS_MESO
              ARIV(I,J,:) = kappam3d(I,J,:)
#else
              AN = SQRT(GRAV * ARHOZ / ARHO)
#ifdef CONSTANT_MESO_LENSCALE
              RD = meso_lenscale_const
#else
              RD = AN * HUP / CORI
              IF (RD.gt.ABS(J-.5*(JM+1))*DYPO(J)) RD=SQRT(AN*HUP/BETA)
#endif
              BYTEADY = GRAV * SQRT(ARHOX*ARHOX + ARHOY*ARHOY) / (AN
     *             *ARHO)
              AINV(I,J) = AMU * RD**2 * BYTEADY ! was = AIN
#endif
#endif
            END IF
#ifdef OCN_GISS_MESO
            AINV(I,J,:) = RGMI * ARIV(I,J,:)
#endif
          END IF
          IM1=I
#ifdef OCN_GISS_MESO
C**** Set diagnostics
          OIJ(I,J,IJ_GMSC) = OIJ(I,J,IJ_GMSC) + AINV(I,J,1) ! GM-scaling
#endif
        END DO
      END DO
!$OMP  END PARALLEL DO
#endif

C**** North pole
      if ( HAVE_NORTH_POLE ) then
        IF (LMM(1,JM).gt.0) THEN
C**** Calculate average density + gradients over [1,LUP]
          ARHO  = 0. ; ARHOY = 0. ;  DZSUMY = 0.
          LAV = MIN(LUP,LMM(1,JM))
          DO L=1,LAV
            ARHO  = ARHO  + RHO(1,JM,L)
            DO I=1,IM
              IF(LMV(I,JM-1).ge.L) THEN
! take abs to get a non-directional scale
                ARHOY = ARHOY + ABS(RHOY(I,JM-1,L))*DZO(L)
                DZSUMY = DZSUMY + DZO(L)
              END IF
            END DO
          END DO
          ARHO  = ARHO / REAL(LAV,KIND=8)
          IF (DZSUMY.gt.0.) ARHOY = ARHOY / DZSUMY
          IF (LAV.gt.1) THEN
            ARHOZ=2*RHOZ1K(1,JM)/(ZE(LAV)+ZE(LAV-1)-ZE(1))
          ELSE
            ARHOZ = 0.
          END IF
C**** avoid occasional inversions. IF ARHOZ<=0 then GM is pure vertical
C**** so keep at zero, and let KPP do the work.
          IF (ARHOZ.gt.0) THEN
#ifdef USE_1D_DIFFUSIVITY
            ARIV(1,JM,:) = k1d(:)
#else
#ifdef OCN_GISS_MESO
            ARIV(1,JM,:) = kappam3d(1,JM,:)
#else
            AN = SQRT(GRAV * ARHOZ / ARHO)
            CORI = ABS(2d0*OMEGA*SINPO(JM))
#ifdef CONSTANT_MESO_LENSCALE
            RD = meso_lenscale_const
#else
            RD = AN * HUP / CORI
#endif
            BYTEADY = GRAV * ARHOY / (AN*ARHO)
            AINV(1,JM) = AMU * RD**2 * BYTEADY ! was = AIN
#endif
#endif
          END IF
#ifdef OCN_GISS_MESO
          AINV(1,JM,:) = RGMI * ARIV(1,JM,:)
#endif
        END IF
#ifdef OCN_GISS_MESO
        DO L=1,LMO
          AINV(2:IM,JM,L)=AINV(1,JM,L)
          ARIV(2:IM,JM,L)=ARIV(1,JM,L)
        ENDDO
C**** Set diagnostics
        OIJ(1,JM,IJ_GMSC) = OIJ(1,JM,IJ_GMSC) + AINV(1,JM,1) ! GM-scaling
#else
        AINV(2:IM,JM)=AINV(1,JM)
#endif
      endif

#ifndef OCN_GISS_MESO
      do j=j_0,j_1
      do i=1,im
        if(lmm(i,j).eq.0) cycle
        ! set isoneutral diffusivity from GM diffusivity
        ariv(i,j) = arai*ainv(i,j)
        ! set diagnostics
        oij(i,j,ij_gmsc) = oij(i,j,ij_gmsc) + ainv(i,j)
      enddo
      enddo
#endif

#ifdef OCN_GISS_MESO
C**** South pole
      if ( HAVE_SOUTH_POLE ) then
        IF (LMM(1,1).gt.0) THEN
C**** Calculate average density + gradients over [1,LUP]
          ARHO  = 0. ; ARHOY = 0. ;  DZSUMY = 0.
          LAV = MIN(LUP,LMM(1,1))
          DO L=1,LAV
            ARHO  = ARHO  + RHO(1,1,L)
            DO I=1,IM
              IF(LMV(I,2).ge.L) THEN
! take abs to get a non-directional scale
                ARHOY = ARHOY + ABS(RHOY(I,2,L))*DZO(L)
                DZSUMY = DZSUMY + DZO(L)
              END IF
            END DO
          END DO
          ARHO  = ARHO / REAL(LAV,KIND=8)
          IF (DZSUMY.gt.0.) ARHOY = ARHOY / DZSUMY
          IF (LAV.gt.1) THEN
            ARHOZ=2*RHOZ1K(1,1)/(ZE(LAV)+ZE(LAV-1)-ZE(1))
          ELSE
            ARHOZ = 0.
          END IF
C**** avoid occasional inversions. IF ARHOZ<=0 then GM is pure vertical
C**** so keep at zero, and let KPP do the work.
          IF (ARHOZ.gt.0) THEN
#ifdef USE_1D_DIFFUSIVITY
            ARIV(1,1,:) = k1d(:)
#else
#ifdef OCN_GISS_MESO
            ARIV(1,1,:) = kappam3d(1,1,:)
#else
            AN = SQRT(GRAV * ARHOZ / ARHO)
            CORI = ABS(2d0*OMEGA*SINPO(JM))
#ifdef CONSTANT_MESO_LENSCALE
            RD = meso_lenscale_const
#else
            RD = AN * HUP / CORI
#endif
            BYTEADY = GRAV * ARHOY / (AN*ARHO)
            ARIV(1,1,:) = AMU * RD**2 * BYTEADY ! was = AIN
#endif
#endif
          END IF
          !ARIV(1,1,:) = ARAI * AINV(1,1,:) ! was = ARI
          AINV(1,1,:) = RGMI * ARIV(1,1,:)
        END IF
c       AINV(2:IM,1)=AINV(1,1)
c       ARIV(2:IM,1)=ARIV(1,1)
        DO L=1,LMO
          AINV(2:IM,1,L)=AINV(1,1,L)
          ARIV(2:IM,1,L)=ARIV(1,1,L)
        ENDDO
C**** Set diagnostics
        OIJ(1,1,IJ_GMSC) = OIJ(1,1,IJ_GMSC) + AINV(1,1,1) ! GM-scaling
      endif
C****

c
c Calculate bolus velocity diagnostics and store them in units of
c accumulated mass flux (kg) for consistency with the diagnostics
c for resolved velocities. For convenience, a constant reference
c density is used to convert the eddy-induced streamfunction to
c mass-flux units.  Local densities can be used if the O(1%) discrepancy
c proves to be a confounding factor in analyses of this diagnostic.
c
c Since the bolus velocity does not appear explicitly in the skew-flux
c form of GM, the attempt has not been made here to exactly reproduce
c the effective interpolation stencil for the diffusivities and density
c gradients used in ISOSLOPE4.  However, a future commit will attempt
c to insert code into ISOSLOPE4 to store the diffusivities in effect
c after modification by CFL criteria and ML exclusion (rather than
c mimicking that logic here).
c 
c For now, only the north-south component is being calculated/stored.
c
      kbyrhoz = 0.
      do l=1,lmo-1
        wtup = dzo(l+1)/(dzo(l)+dzo(l+1))
        wtdn = 1.-wtup
        do j=j_0,j_1
        do n=1,nbyzm(j,l+1)
        do i=i1yzm(n,j,l+1),i2yzm(n,j,l+1)
          if(l.gt.kpl(i,j)) then ! GM excludes ML currently
            ! interpolate K/(drho/dz) to layer edges
            kbyrhoz(i,j,l) = (wtup*ainv(i,j,l)+wtdn*ainv(i,j,l+1))*
     &           byrhoz(i,j,l)
          endif
        enddo
        enddo
        enddo
      enddo
      call halo_update(grid,kbyrhoz,from=north)
      call halo_update(grid,byrhoz,from=north) ! needed?
      psiy = 0.
      rhomid = rhows ! for now
      do l=1,lmo-1
        kmax = ( .25*(dzo(l)+dzo(l+1))**2 )/dts
        wtup = dzo(l+1)/(dzo(l)+dzo(l+1))
        wtdn = 1.-wtup
        do j=j_0s,j_1s
        do n=1,nbyzv(j,l+1)
        do i=i1yzv(n,j,l+1),i2yzv(n,j,l+1)
          ! verticall interpolate rhoy to layer edges
          rhoy_ = (wtup*rhoy(i,j,l)+wtdn*rhoy(i,j,l+1))
          ! calculate y-slope
          sly = rhoy_*(byrhoz(i,j,l)+byrhoz(i,j+1,l))*.5
          ! streamfunc = K
          psiy(i,j,l) = rhoy_ *(kbyrhoz(i,j,l)+kbyrhoz(i,j+1,l))*.5
          ! mimic the "CFL tapering" in ISOSLOPE4: locally reduce
          ! K so that the vertical component of Redi diffusion does
          ! not violate CFL conditions for numerical stability
          krat = psiy(i,j,l)*sly/kmax
          if(krat .gt. 1d0) psiy(i,j,l) = psiy(i,j,l)/krat
          psiy(i,j,l) = psiy(i,j,l)*rhomid
        enddo
        enddo
        enddo
      enddo
      do l=1,lmo
        do j=j_0s,j_1s
        do n=1,nbyzv(j,l)
        do i=i1yzv(n,j,l),i2yzv(n,j,l)
          oijl(i,j,l,ijl_mfvb) = oijl(i,j,l,ijl_mfvb) +
     &         (psiy(i,j,l)-psiy(i,j,l-1))*dxvo(j)*dts
        enddo
        enddo
        enddo
      enddo

#endif

      RETURN
C****
      END SUBROUTINE DENSGRAD

      Subroutine VBAR_GM0
!@sum VBAR_GM0 calculates specific volume and vertical gradients for GM
!@auth Gary Russell/Gavin Schmidt
      Use CONSTANT, only: grav
      Use GM_COM, only: vbar,dvbardz,rhoz1k,lup,g3d,s3d,p3d
      Use OCEAN, Only: IM,JM,LMO, LMOM=>LMM,DXYPO, MO, ZE,
     *                 G0M,GZM=>GZMO, S0M,SZM=>SZMO, OPRESS, FOCEAN
      Use DOMAIN_DECOMP_1D, Only: HALO_UPDATE, NORTH
      USE OCEANR_DIM, only : grid=>ogrid

      Implicit None
      Real*8,Parameter :: z12eH=.28867513d0  !  z12eH = 1/SQRT(12)
      Integer*4 I,J,L,IMAX, J1,JN,JNH, LAV,LAVM
      Real*8, dimension(lmo) :: gup,gdn,sup,sdn,pm,gmd,smd
      Real*8 :: vup,vdn,vupu,vdnu,pe,bym
      Logical*4 QSP,QNP
      Real*8,External   :: VOLGSP,temgs

C****
C**** Extract domain decomposition band parameters
C****                          Band1  Band2  BandM
      J1  = GRID%J_STRT     !    1      5     JM-3   Band minimum
      JN  = GRID%J_STOP     !    4      8     JM     Band maximum
      JNH = Min(JN+1,JM)    !    5      9     JM     Halo maximum
      QSP = J1==1           !    T      F      F
      QNP = JN==JM          !    F      F      T
C****
      Call HALO_UPDATE (GRID,OPRESS,FROM=NORTH)
      Call HALO_UPDATE (GRID,  MO, FROM=NORTH)
      Call HALO_UPDATE (GRID, G0M, FROM=NORTH)
      Call HALO_UPDATE (GRID, GZM, FROM=NORTH)
      Call HALO_UPDATE (GRID, S0M, FROM=NORTH)
      Call HALO_UPDATE (GRID, SZM, FROM=NORTH)

      Do J=J1,JNH
        IMAX=IM  ;  If(J==1.or.J==JM) IMAX=1
        Do I=1,IMAX
          If (FOCEAN(I,J) == 0)  CYCLE

          Do L=1,LMOM(I,J)
            BYM=1d0/(MO(I,J,L)*DXYPO(J))
            GMD(L)= G0M(I,J,L)*BYM
            GUP(L)=(G0M(I,J,L)-2*z12eH*GZM(I,J,L))*BYM
            GDN(L)=(G0M(I,J,L)+2*z12eH*GZM(I,J,L))*BYM
            SMD(L)= S0M(I,J,L)*BYM
            SUP(L)=(S0M(I,J,L)-2*z12eH*SZM(I,J,L))*BYM
            SDN(L)=(S0M(I,J,L)+2*z12eH*SZM(I,J,L))*BYM
          EndDo
          
C**** Calculate pressure by integrating from the top down
          PE = OPRESS(I,J)
          Do L=1,LMOM(I,J)
            PM(L) = PE + MO(I,J,L)*GRAV*.5
            PE    = PE + MO(I,J,L)*GRAV
            G3D(I,J,L) = GMD(L)
            S3D(I,J,L) = SMD(L)
            P3D(I,J,L) = PM(L)
          EndDo
C**** Calculate potential specific volume (ref to mid-point pr)
          Do L=LMOM(I,J),1,-1
            VUP = VOLGSP (GUP(L),SUP(L),PM(L))
            VDN = VOLGSP (GDN(L),SDN(L),PM(L))
            VBAR(I,J,L) = (VUP + VDN)*.5
C**** Vertical gradient calculated using lower box mid-point pr
            IF (L.gt.1) then 
              VUPU = VOLGSP (GUP(L-1),SUP(L-1),PM(L))
              VDNU = VOLGSP (GDN(L-1),SDN(L-1),PM(L))
              dVBARdZ(I,J,L) = .5* (VUP + VDN - VUPU - VDNU)
            end if
          EndDo
C**** Vertical potential gradient in top 1km
          LAV = MIN(LUP,LMOM(I,J))
          LAVM = MAX(LAV/2,1)   ! mid depth
          VUP = VOLGSP (GMD(1),SMD(1),PM(LAVM))
          VDN = VOLGSP (GMD(LAV),SMD(LAV),PM(LAVM))
          RHOZ1K(I,J) = (VUP - VDN)/VBAR(I,J,LAVM)**2
        EndDo
      EndDo
C**** Copy VBAR to all longitudes at poles
      If (QNP) Then
        Do L=1,LMOM(1,JM)
          VBAR(2:IM,JM,L) = VBAR(1,JM,L)
          G3D(2:IM,JM,L) = G3D(1,JM,L)
          S3D(2:IM,JM,L) = S3D(1,JM,L)
          P3D(2:IM,JM,L) = P3D(1,JM,L)
          dVBARdZ(2:IM,JM,L) = dVBARdZ(1,JM,L)
        EndDo
      EndIf
      If (QSP) Then
        Do L=1,LMOM(1,1)
          VBAR(2:IM,1,L) = VBAR(1,1,L)
          G3D(2:IM,1,L) = G3D(1,1,L)
          S3D(2:IM,1,L) = S3D(1,1,L)
          P3D(2:IM,1,L) = P3D(1,1,L)
          dVBARdZ(2:IM,1,L) = dVBARdZ(1,1,L)
        EndDo
      EndIf

      Return
      End Subroutine VBAR_GM0
