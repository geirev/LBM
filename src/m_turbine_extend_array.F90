!==============================================================
!  m_turbine_extend.F90
!  Small utility to grow allocatable arrays of point_t
!==============================================================
module m_turbine_extend_array
   use mod_turbines, only : point_t
   implicit none
contains

!--------------------------------------------------------------
!  subroutine turbine_extend_array
!
!  PURPOSE:
!    Grow an allocatable array of point_t to size n.
!
!  NOTE:
!    Copies existing data and move_allocs to original.
!--------------------------------------------------------------
subroutine turbine_extend_array(arr, n)
   type(point_t), allocatable, intent(inout) :: arr(:)
   integer,           intent(in)            :: n

   type(point_t), allocatable :: tmp(:)
   integer :: old

   if (.not. allocated(arr)) then
      allocate(arr(n))
   else
      allocate(tmp(n))
      old = size(arr)
      tmp(1:old) = arr(1:old)
      call move_alloc(tmp, arr)
   end if
end subroutine turbine_extend_array

end module m_turbine_extend_array
