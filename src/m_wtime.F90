module m_wtime
!@cuf use cudafor
   real(kind(1.0D+00)) :: xt0,xt1
   integer, parameter  :: nrtimes=16
   real(kind(1.0D+00)) :: walltime(1:nrtimes)=0.0
   integer istat
contains

subroutine cpustart()
!@cuf istat = cudaDeviceSynchronize()
   xt0 = wtime()
end

subroutine cpufinish(icpu)
!@cuf istat = cudaDeviceSynchronize()
   xt1=wtime();  walltime(icpu)=walltime(icpu)+xt1-xt0
end

subroutine cpuprint()
   implicit none
   integer l
   print '(a,1x,2f13.7)','initialization     time =',walltime(1)
   print '(a,1x,2f13.7)','turbine forcing    time =',walltime(2)
   print '(a,1x,2f13.7)','turbulence forcing time =',walltime(3)
   print '(a,1x,2f13.7)','equil              time =',walltime(4)
   print '(a,1x,2f13.7)','regularization     time =',walltime(5)
   print '(a,1x,2f13.7)','vreman             time =',walltime(6)
   print '(a,1x,2f13.7)','collisions         time =',walltime(7)
   print '(a,1x,2f13.7)','applyturbines      time =',walltime(8)
   print '(a,1x,2f13.7)','applyturbulence    time =',walltime(9)
   print '(a,1x,2f13.7)','solids             time =',walltime(10)
   print '(a,1x,2f13.7)','boundarycond       time =',walltime(11)
   print '(a,1x,2f13.7)','drift              time =',walltime(12)
   print '(a,1x,2f13.7)','macrovars          time =',walltime(13)
   print '(a,1x,2f13.7)','diag               time =',walltime(14)
   print '(a,1x,2f13.7)','averaging and turb time =',walltime(15)
   print '(a,1x,2f13.7)','final stuff        time =',walltime(16)
   print '(a,1x,2f13.7)','Total wall time    time =',sum(walltime(1:nrtimes))
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
