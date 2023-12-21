module m_wtime
   real start,finish,cpu0,cpu1
   real :: cputime(10)=0.0
   real :: waltime(10)=0.0
   character(len=9) :: cpuname(1:8) = ['drift    ',&
                                       'defbnd   ',&
                                       'rho      ',&
                                       'velocity ',&
                                       'feq      ',&
                                       'collision',&
                                       'bndapply ',&
                                       'printing ']
contains
function wtime ( )

!*****************************************************************************80
!
!! WTIME returns a reading of the wall clock time.
!
!  Discussion:
!
!    To get the elapsed wall clock time, call WTIME before and after a given
!    operation, and subtract the first reading from the second.
!
!    This function is meant to suggest the similar routines:
!
!      "omp_get_wtime ( )" in OpenMP,
!      "MPI_Wtime ( )" in MPI,
!      and "tic" and "toc" in MATLAB.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 April 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = rk ) WTIME, the wall clock reading, in seconds.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer clock_max
  integer clock_rate
  integer clock_reading
  real ( kind = rk ) wtime

  call system_clock ( clock_reading, clock_rate, clock_max )

  wtime = real(clock_reading, kind = rk) / real(clock_rate, kind = rk)

end function
end module
