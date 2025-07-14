module m_wtime
   real start,finish,cpu0,cpu1
   integer, parameter :: nrtimes=13
   real :: cputime(1:nrtimes)=0.0
   real :: waltime(1:nrtimes)=0.0
   character(len=9) :: cpuname(1:nrtimes) = ['rho      ',&
                                       'velocity ',&
                                       'printing ',&
                                       'fequil   ',&
                                       'turbine  ',&
                                       'collision',&
                                       'applyturb',&
                                       'boundary ',&
                                       'solids   ',&
                                       'drift    ',&
                                       'fequil3  ',&
                                       'fregular ',&
                                       'vreman   ']
contains

subroutine cpustart()
   cpu0=wtime()
   call cpu_time(start)
end

subroutine cpufinish(icpu)
   integer, intent(in) :: icpu
   cpu1=wtime()
   call cpu_time(finish)
   cputime(icpu)=cputime(icpu)+finish-start
   waltime(icpu)=waltime(icpu)+cpu1-cpu0
end

subroutine cpuprint()
   implicit none
   integer l
   print '(tr22,3a)','cputime    ','walltime    ','speedup     '
   do l=1,nrtimes
      print '(tr10,a9,3f12.4)',cpuname(l),cputime(l),waltime(l),cputime(l)/(waltime(l)+tiny(cpu1))
   enddo
   print '(tr10,a9,3f12.4)','summary  ',sum(cputime(1:nrtimes)),sum(waltime(1:nrtimes)),&
                                        sum(cputime(1:nrtimes))/sum(waltime(1:nrtimes))
end subroutine

function wtime ( )

!*****************************************************************************
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
