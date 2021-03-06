module testGeometry_mod
!@sum testGeometry_mod contains tests for the various convenience
!@+   parameters and arrays that are associated with the lat-lon grid.
!@+   Essentially these are test routines for GEOM_B.f
!@auth Tom Clune <thomas.l.clune@nasa.gov>
  use pFUnit
  use Geometry_mod
  implicit none
  private

  public :: testGetIndicesFromLatLon
  public :: testGetIndices4x5
  public :: testGetIndices8x10

contains

  subroutine testGetIndicesFromLatLon()
!@sum Is the correct latitude index determined from latitude (in degrees)
    
    integer, parameter :: IM = 3, JM = 4, KM = 2 ! coarse grid
    type (Geometry_type) :: geometry

    type Case_type
      real*8  :: degreesLongitude
      real*8  :: degreesLatitude
      integer :: expectedIndices(2)
    end type Case_type

    integer, parameter :: numCases = 20
    type (Case_type) :: cases(numCases) = [ &
         & Case_type(0., -91., [2,1]), &
         & Case_type(0., -90., [2,1]), &
         & Case_type(0., -46., [2,1]), &
         & Case_type(0., -44., [2,2]), &
         & Case_type(0., -01., [2,2]), &
         & Case_type(0., +01., [2,3]), &
         & Case_type(0., +44., [2,3]), &
         & Case_type(0., +46., [2,4]), &
         & Case_type(0., +90., [2,4]), &
         & Case_type(0., +91., [2,4]), &

         & Case_type(-181., 20., [1,3]), &
         & Case_type(-180., 20., [1,3]), &
         & Case_type(-179., 20., [1,3]), &
         & Case_type( -61., 20., [1,3]), &
         & Case_type( -59., 20., [2,3]), &
         & Case_type( +59., 20., [2,3]), &
         & Case_type( +61., 20., [3,3]), &
         & Case_type(+179., 20., [3,3]), &
         & Case_type(+180., 20., [3,3]), &
         & Case_type(+181., 20., [3,3])  &
         & ]

    integer :: i

    geometry = newGeometry(IM, JM, KM)

    do i = 1, numCases
      call checkCase(geometry, cases(i))
    end do
    
  contains

    subroutine checkCase(geometry, case)
      type (Geometry_type), intent(in) :: geometry
      type (Case_type), intent(in) :: case

      call assertEqual(case%expectedIndices, &
           & getIndicesFromLatLon(geometry, case%degreesLongitude, case%degreesLatitude))
    end subroutine checkCase
    
  end subroutine testGetIndicesFromLatLon

  subroutine testGetIndices4x5()
!@sum Ensure that box at poles on modelE 4x5 grid is half size of others.
    integer, parameter :: IM = 72
    integer, parameter :: JM = 46
    integer, parameter :: KM = 2
    type (Geometry_type) :: geometry

    real*8, parameter :: lonSpacing = 5 ! degrees
    real*8, parameter :: latSpacing = 4 ! degrees
    real*8, parameter :: southPoleBoxCenter = -89. ! degrees

    real*8 :: degreesLat
    real*8 :: degreesLon
    integer :: i, j

    geometry = newGeometry4x5(KM)

    do j = 1, JM
      degreesLat = southPoleBoxCenter + (j - 1)*latSpacing
      do i = 1, IM
        degreesLon = -180. + (i - 0.5)*lonSpacing
        call assertEqual([i,j], getIndicesFromLatLon(geometry, degreesLon, degreesLat))
      end do
    end do

  end subroutine testGetIndices4x5

  subroutine testGetIndices8x10()
!@sum Ensure that box at poles on modelE 8x10 grid is quarter size of others.
    integer, parameter :: IM = 36
    integer, parameter :: JM = 24
    integer, parameter :: KM = 2
    type (Geometry_type) :: geometry

    real*8, parameter :: lonSpacing = 10 ! degrees
    real*8, parameter :: latSpacing = 8 ! degrees
    real*8, parameter :: southPoleBoxCenter = -89. ! degrees

    real*8 :: degreesLat
    real*8 :: degreesLon
    integer :: i, j

    geometry = newGeometry8x10(KM)

    do j = 1, JM
      degreesLat = southPoleBoxCenter + (j - 1)*latSpacing
      do i = 1, IM
        degreesLon = -180. + (i - 0.5)*lonSpacing
        call assertEqual([i,j], getIndicesFromLatLon(geometry, degreesLon, degreesLat))
      end do
    end do

  end subroutine testGetIndices8x10
  
end module testGeometry_mod
