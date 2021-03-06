module test_ReportColumn_mod
   use pFUnit, only: assertEqual
   use ReportColumn_mod
   use Timer_mod
   use TimeFormatUtilities_mod
   implicit none
   private

   public :: test_getFieldWidthNameColumn
   public :: test_getHeaderNameColumnA
   public :: test_getHeaderNameColumnB
   public :: test_getFieldNameColumn

   public :: test_getFieldWidthInclusiveTime
   public :: test_getHeaderInclusiveTime
   public :: test_getFieldInclusiveTime

   public :: test_getFieldWidthExclusiveTime
   public :: test_getHeaderExclusiveTime
   public :: test_getFieldExclusiveTime

   public :: test_getHeaderTripCounts
   public :: test_getFieldTripCounts

   public :: test_getHeaderMaximumTime
   public :: test_getFieldMaximumTime

   public :: test_getHeaderMinimumTime
   public :: test_getFieldMinimumTime

   public :: test_getHeaderAverageTime
   public :: test_getFieldAverageTime

   public :: test_setScale
   public :: test_setPrecision
   public :: test_getUnits

contains

   subroutine test_getFieldWidthNameColumn()
      type (ReportColumn_type) :: column
      integer :: fieldWidth
      fieldWidth = 5
      column = newColumn(NAME_COLUMN, fieldWidth)
      call assertEqual(fieldWidth, getFieldWidth(column))
   end subroutine test_getFieldWidthNameColumn

   subroutine test_getHeaderNameColumnA()
      type (ReportColumn_type) :: column
      integer :: fieldWidth
      fieldWidth = 6
      column = newColumn(NAME_COLUMN, fieldWidth)
      call assertEqual(fieldWidth, len(getHeader(column)))
      call assertEqual(' Name ', getHeader(column))
   end subroutine test_getHeaderNameColumnA

   subroutine test_getHeaderNameColumnB()
      type (ReportColumn_type) :: column
      integer :: fieldWidth
      fieldWidth = 3
      column = newColumn(NAME_COLUMN, fieldWidth)
      call assertEqual(fieldWidth, len(getHeader(column)))
      call assertEqual('Nam', getHeader(column))
   end subroutine test_getHeaderNameColumnB

   subroutine test_getFieldNameColumn()
      type (ReportColumn_type) :: column
      type (Timer_type) :: timer
      integer :: fieldWidth
      fieldWidth = 6
      column = newColumn(NAME_COLUMN, fieldWidth)
      call assertEqual(fieldWidth, len(getField(column, name = 'abcdefghij', timer=timer)))
      call assertEqual('abcdef', getField(column, name = 'abcdefghij', timer=timer))
   end subroutine test_getFieldNameColumn

   subroutine test_getFieldWidthInclusiveTime()
      type (ReportColumn_type) :: column
      integer :: fieldWidth
      fieldWidth = 5
      column = newColumn(INCLUSIVE_TIME_COLUMN, fieldWidth)
      call assertEqual(fieldWidth, getFieldWidth(column))
   end subroutine test_getFieldWidthInclusiveTime

   subroutine test_getHeaderInclusiveTime()
      type (ReportColumn_type) :: column
      integer :: fieldWidth
      fieldWidth = 6
      column = newColumn(INCLUSIVE_TIME_COLUMN, fieldWidth)
      call assertEqual(fieldWidth, len(getHeader(column)))
      call assertEqual('Inclus', getHeader(column))
   end subroutine test_getHeaderInclusiveTime

   subroutine test_getFieldInclusiveTime()
      type (ReportColumn_type) :: column
      type (Timer_type) :: timer
      integer :: fieldWidth
      real(kind=r64) :: time

      fieldWidth = 6
      column = newColumn(INCLUSIVE_TIME_COLUMN, fieldWidth)
      call assertEqual(fieldWidth, len(getField(column, name = 'abcdefghij', timer=timer)))

      time = 12.34
      call start(timer, 0._r64)
      call stop(timer, time)

      call assertEqual(formatSeconds(time,3,2), getField(column, name = 'abcdefghij', timer=timer))
   end subroutine test_getFieldInclusiveTime

   subroutine test_getFieldWidthExclusiveTime()
      type (ReportColumn_type) :: column
      integer :: fieldWidth
      fieldWidth = 5
      column = newColumn(EXCLUSIVE_TIME_COLUMN, fieldWidth)
      call assertEqual(fieldWidth, getFieldWidth(column))
   end subroutine test_getFieldWidthExclusiveTime

   subroutine test_getHeaderExclusiveTime()
      type (ReportColumn_type) :: column
      integer :: fieldWidth
      fieldWidth = 6
      column = newColumn(EXCLUSIVE_TIME_COLUMN, fieldWidth)
      call assertEqual(fieldWidth, len(getHeader(column)))
      call assertEqual('Exclus', getHeader(column))
   end subroutine test_getHeaderExclusiveTime

   subroutine test_getFieldExclusiveTime()
      type (ReportColumn_type) :: column
      type (Timer_type) :: timer
      type (Timer_type) :: otherTimer
      integer :: fieldWidth
      real(kind=r64) :: time

      fieldWidth = 6
      column = newColumn(EXCLUSIVE_TIME_COLUMN, fieldWidth)
      call assertEqual(fieldWidth, len(getField(column, name = 'abcdefghij', timer=timer)))

      time = 12.34
      call start(timer, 0._r64)
      call start(otherTimer, 0._r64)
      call stop(otherTimer, 1._r64)
      call stop(timer, time)

      call assertEqual(formatSeconds(time-1,3,2), getField(column, name = 'abcdefghij', timer=timer))
   end subroutine test_getFieldExclusiveTime

   subroutine test_getHeaderTripCounts()
      type (ReportColumn_type) :: column
      integer :: fieldWidth
      fieldWidth = 12
      column = newColumn(TRIP_COUNTS_COLUMN, fieldWidth)
      call assertEqual(fieldWidth, len(getHeader(column)))
      call assertEqual('   Trips ', getHeader(column))
   end subroutine test_getHeaderTripCounts

   subroutine test_getFieldTripCounts()
      type (ReportColumn_type) :: column
      type (Timer_type) :: timer
      integer :: fieldWidth
      real(kind=r64) :: time

      fieldWidth = 6
      column = newColumn(TRIP_COUNTS_COLUMN, fieldWidth)
      call assertEqual(fieldWidth, len(getField(column, name = 'abcdefghij', timer=timer)))

      call start(timer, 0._r64)
      call stop(timer, time)

      call assertEqual('     1', getField(column, name = 'abcdefghij', timer=timer))
   end subroutine test_getFieldTripCounts

   subroutine test_getHeaderMaximumTime()
      type (ReportColumn_type) :: column
      integer :: fieldWidth
      fieldWidth = 12
      column = newColumn(MAXIMUM_TIME_COLUMN, fieldWidth)
      call assertEqual(fieldWidth, len(getHeader(column)))
      call assertEqual('  Maximum   ', getHeader(column))
   end subroutine test_getHeaderMaximumTime

   subroutine test_getFieldMaximumTime()
      type (ReportColumn_type) :: column
      type (Timer_type) :: timer
      integer :: fieldWidth
      real(kind=r64) :: time

      fieldWidth = 6
      column = newColumn(MAXIMUM_TIME_COLUMN, fieldWidth)
      call assertEqual(fieldWidth, len(getField(column, name = 'abcdefghij', timer=timer)))

      time = 12.34
      call start(timer, 0._r64)
      call stop(timer, time)
      call start(timer, time)
      call stop(timer, 3*time)

      call assertEqual(formatSeconds(2*time,3,2), getField(column, name = 'abcdefghij', timer=timer))
   end subroutine test_getFieldMaximumTime

   subroutine test_getHeaderMinimumTime()
      type (ReportColumn_type) :: column
      integer :: fieldWidth
      fieldWidth = 12
      column = newColumn(MINIMUM_TIME_COLUMN, fieldWidth)
      call assertEqual(fieldWidth, len(getHeader(column)))
      call assertEqual('  Minimum   ', getHeader(column))
   end subroutine test_getHeaderMinimumTime

   subroutine test_getFieldMinimumTime()
      type (ReportColumn_type) :: column
      type (Timer_type) :: timer
      integer :: fieldWidth
      real(kind=r64) :: time

      fieldWidth = 6
      column = newColumn(MINIMUM_TIME_COLUMN, fieldWidth)
      call assertEqual(fieldWidth, len(getField(column, name = 'abcdefghij', timer=timer)))

      time = 12.34
      call start(timer, 0._r64)
      call stop(timer, time)
      call start(timer, time)
      call stop(timer, 3*time)

      call assertEqual(formatSeconds(1*time,3,2), getField(column, name = 'abcdefghij', timer=timer))
   end subroutine test_getFieldMinimumTime

   subroutine test_getHeaderAverageTime()
      type (ReportColumn_type) :: column
      integer :: fieldWidth
      fieldWidth = 12
      column = newColumn(AVERAGE_TIME_COLUMN, fieldWidth)
      call assertEqual(fieldWidth, len(getHeader(column)))
      call assertEqual('  Average   ', getHeader(column))
   end subroutine test_getHeaderAverageTime

   subroutine test_getFieldAverageTime()
      type (ReportColumn_type) :: column
      type (Timer_type) :: timer
      integer :: fieldWidth
      real(kind=r64) :: time

      fieldWidth = 6
      column = newColumn(AVERAGE_TIME_COLUMN, fieldWidth)
      call assertEqual(fieldWidth, len(getField(column, name = 'abcdefghij', timer=timer)))

      time = 12.34
      call start(timer, 0._r64)
      call stop(timer, time)
      call start(timer, time)
      call stop(timer, 3*time)

      call assertEqual(formatSeconds(1.5*time,3,2), getField(column, name = 'abcdefghij', timer=timer))
   end subroutine test_getFieldAverageTime

   subroutine test_setScale()
      type (ReportColumn_type) :: column
      type (Timer_type) :: timer
      type (Timer_type) :: otherTimer
      real(kind=r64) :: time, scale

      time = 12.34
      call start(timer, 0._r64)
      call stop(timer, time)

      column = newColumn(INCLUSIVE_TIME_COLUMN, 6)
      scale = 0.1_r64
      call setScale(column, scale, 'units')

      call assertEqual('units', getUnits(column))
      call assertEqual(formatSeconds(time*scale,3,2), getField(column, name = 'abcdefghij', timer=timer))

   end subroutine test_setScale

   subroutine test_setPrecision()
      type (ReportColumn_type) :: column
      type (Timer_type) :: timer
      type (Timer_type) :: otherTimer
      real(kind=r64) :: time, scale

      time = 12.34
      call start(timer, 0._r64)
      call stop(timer, time)

      column = newColumn(INCLUSIVE_TIME_COLUMN, 8)
      call setPrecision(column, 4)
      call assertEqual(formatSeconds(time,3,4), getField(column, name = 'abcdefghij', timer=timer))
   end subroutine test_setPrecision

   subroutine test_getUnits()
      type (ReportColumn_type) :: column
      real(kind=r64) :: scale

      call assertEqual('', getUnits(newColumn(NAME_COLUMN, 1)))
      call assertEqual('', getUnits(newColumn(TRIP_COUNTS_COLUMN, 1)))
      call assertEqual('seconds ', getUnits(newColumn(INCLUSIVE_TIME_COLUMN, 8)))
      call assertEqual('seconds ', getUnits(newColumn(EXCLUSIVE_TIME_COLUMN, 8)))
      call assertEqual('seconds ', getUnits(newColumn(MAXIMUM_TIME_COLUMN, 8)))
      call assertEqual('seconds ', getUnits(newColumn(MINIMUM_TIME_COLUMN, 8)))
      call assertEqual('seconds ', getUnits(newColumn(AVERAGE_TIME_COLUMN, 8)))

      scale = 1.
      column = newColumn(INCLUSIVE_TIME_COLUMN, 4)
      call setScale(column, scale, 'units')
      call assertEqual('unit', getUnits(column))

   end subroutine test_getUnits

end module test_ReportColumn_mod
