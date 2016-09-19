#ifdef USE_ESMF
      module domain_decomp
      implicit none

#include "mpi_defs.h"
#include "mpif.h"

      integer npes, mype
      integer, allocatable :: jbounds(:)
      integer j0,j1

      contains

      subroutine app_init(jm,j0_out,j1_out)
      integer, intent(in) :: jm
      integer, intent(out) :: j0_out,j1_out
      !----
      integer ierr !, j
cddd      real*8 dj
cddd      integer :: sumj(46) = (/
cddd     & 0,    0,    0,    0,    0,    0,
cddd     & 0,    0,    0,    1,    1,    2,
cddd     & 3,    5,    8,   15,   16,   18,
cddd     &17,   15,   15,   18,   17,   16,
cddd     &17,   19,   17,   21,   24,   29,
cddd     &32,   29,   32,   35,   39,   41,
cddd     &44,   38,   49,   54,   38,   15,
cddd     & 7,    2,    0,    0 /)


      call MPI_INIT(ierr)
      call MPI_COMM_RANK(MPI_COMM_WORLD, mype, ierr)
      call MPI_COMM_SIZE(MPI_COMM_WORLD, npes, ierr)
      !print *,"status = ", ierr

      !print *,"npes, mype= ", npes, mype

      allocate( jbounds(0:npes+1) )

cddd      dj = dble(jm-12)/dble(npes)
      jbounds(0) = 0
cddd      do j=1,npes-1
cddd        jbounds(j) = dj*j + 8
cddd      enddo
      jbounds(npes) = jm

cddd      dj = sum( sumj ) / dble(npes)
cddd      do j=1,jm
cddd        jbounds( int(sum( sumj(1:j) )/dj) + 1 ) = j
cddd      enddo
cddd      jbounds(npes) = jm

      if ( npes == 4 ) then
        !jbounds(1) = jbounds(1) - 3 ! works for npes=4
        jbounds(1) = 23
        jbounds(2) = 33
        jbounds(3) = 38
        jbounds(4) = 46
      endif

cddd      if ( npes == 8 ) then
cddd        jbounds(1) = 18
cddd        jbounds(2) = 22
cddd        jbounds(3) = 28
cddd        jbounds(4) = 32
cddd        jbounds(5) = 35
cddd        jbounds(6) = 37
cddd        jbounds(7) = 39
cddd        jbounds(8) = 46
cddd      endif

      if ( npes == 8 ) then
        jbounds(1) = 18
        jbounds(2) = 23
        jbounds(3) = 29
        jbounds(4) = 33
        jbounds(5) = 35
        jbounds(6) = 37
        jbounds(7) = 39
        jbounds(8) = 46
      endif

      if ( npes == 16 ) then
        jbounds(1) = 16
        jbounds(2) = 18
        jbounds(3) = 21
        jbounds(4) = 23
        jbounds(5) = 26
        jbounds(6) = 29
        jbounds(7) = 31
        jbounds(8) = 33
        jbounds(9) = 34
        jbounds(10) = 35
        jbounds(11) = 36
        jbounds(12) = 37
        jbounds(13) = 38
        jbounds(14) = 39
        jbounds(15) = 40
        jbounds(16) = 46
      endif

! hack to deal with jm != 46
      if ( jm==90 ) then
        jbounds(:) = jbounds(:)*2
        jbounds(npes) = jm
      endif

      j0 = jbounds(mype) + 1
      j1 = jbounds(mype+1)
      j0_out = j0
      j1_out = j1

      end subroutine app_init


      subroutine array_gather( a )
      real*8 a(:,:)
      integer :: status
      integer :: sendcount, recvcounts(npes), displs(npes)
      integer idim, jdim, n
      real*8, allocatable :: b(:,:)

      idim = size(a,1)
      jdim = size(a,2)

      allocate( b(idim,jdim) )

      do n=1,npes
        recvcounts(n) = idim*(jbounds(n) - jbounds(n-1))
        displs(n) = idim*jbounds(n-1)
!        print *,"mype",mype,"n",n,"recvcounts(n)",recvcounts(n)
!     &       ,"displs(n)",displs(n),"jbounds(n)",jbounds(n)
      enddo
      sendcount = recvcounts(mype+1)

!      print *,"recvcounts = ", recvcounts
!      print *,"displs = ", displs
!      print *,"sendcount = ", sendcount

      call MPI_Barrier(MPI_COMM_WORLD, status)
      
      call MPI_GatherV(a(:,j0:j1), sendcount, MPI_DOUBLE_PRECISION,
     &     b, recvcounts, displs,
     &     MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, status)

!      print *,"status = ", status

      a = b

      deallocate( b )

      end subroutine array_gather

      subroutine array_bcast_r4( a )
      real*4 a(:,:)
      integer :: status

      call MPI_Barrier(MPI_COMM_WORLD, status)

      call MPI_Bcast ( a, size(a), MPI_REAL, 0, 
     &     MPI_COMM_WORLD, status )

      end subroutine array_bcast_r4

      end module domain_decomp
#else
      module domain_decomp
      implicit none

      integer npes, mype
      integer, allocatable :: jbounds(:)
      !integer j0,j1

      contains

      subroutine app_init(jm,j0_out,j1_out)
      integer, intent(in) :: jm
      integer, intent(out) :: j0_out,j1_out
      !----

      npes = 1
      mype = 0
      allocate( jbounds(0:npes+1) )
      jbounds(0) = 0
      jbounds(1) = jm
      jbounds(2) = jm ! not needed, just in case...

      j0_out = 1
      j1_out = jm
      end subroutine app_init

      subroutine array_gather( a )
      real*8 a(:,:)
      ! do nothing
      end subroutine array_gather

      subroutine array_bcast_r4( a )
      real*4 a(:,:)
      ! do nothing
      end subroutine array_bcast_r4

      end module domain_decomp
#endif

