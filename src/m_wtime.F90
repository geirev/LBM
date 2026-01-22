module m_wtime
#ifdef _CUDA
   use cudafor
#endif
   implicit none
   integer, parameter  :: nrtimes = 30
   real(kind(1.0D+00)) :: walltime(1:nrtimes) = 0.0D+00
   real(kind(1.0D+00)) :: tstart = 0.0D+00
   integer :: istat
   integer :: clock_max
   integer :: clock_rate
   integer :: clock_reading

contains

subroutine cpustart()
   use m_readinfile, only : ltiming
   implicit none
   if (.not.ltiming) return
   tstart = wtime()
end subroutine

subroutine cpufinish(icpu)
   use m_readinfile, only : ltiming
   implicit none
   integer, intent(in) :: icpu
   real(kind(1.0D+00)) :: tend
   real(kind(1.0D+00)) :: elapsed
   if (.not.ltiming) return

   if (icpu < 1 .or. icpu > nrtimes) then
      print *, "Error: Invalid timer index icpu = ", icpu
      stop
   endif
   tend=wtime()
   if (tend >= tstart) then
      elapsed = real(tend - tstart) / real(clock_rate)
   else
      elapsed = real((clock_max - tstart + tend + 1)) / real(clock_rate)
   end if
   walltime(icpu) = walltime(icpu) + elapsed
!   if (icpu == 4) print '(a,i3,3f13.5)','icpu=',icpu,tstart/real(clock_rate),tend/real(clock_rate),walltime(icpu)
end subroutine

subroutine cpuprint()
   use m_readinfile, only : ltiming
   implicit none
   if (.not.ltiming) return
   print '(a,1x,f13.5)','initialization     time =',walltime(1)
   print '(a,1x,f13.5)','turbine forcing    time =',walltime(2)
   print '(a,1x,f13.5)','turbulence forcing time =',walltime(3)
   print '(a,1x,f13.5)','equil              time =',walltime(4)
   print '(a,1x,f13.5)','regularization     time =',walltime(5)
   print '(a,1x,f13.5)','vreman             time =',walltime(6)
   print '(a,1x,f13.5)','forcings_apply     time =',walltime(7)
   print '(a,1x,f13.5)','buoyancy forcing   time =',walltime(8)
   print '(a,1x,f13.5)','applyturbulence    time =',walltime(9)
   print '(a,1x,f13.5)','solids             time =',walltime(10)
   print '(a,1x,f13.5)','boundarycond       time =',walltime(11)
   print '(a,1x,f13.5)','drift              time =',walltime(12)
   print '(a,1x,f13.5)','macrovars          time =',walltime(13)
   print '(a,1x,f13.5)','diag               time =',walltime(14)
   print '(a,1x,f13.5)','averaging and turb time =',walltime(15)
   print '(a,1x,f13.5)','final stuff        time =',walltime(16)
   print '(a,1x,f13.5)','compute_f(neq)     time =',walltime(17)
   print '(a,1x,f13.5)','post coll          time =',walltime(18)
   print '(a,1x,f13.5)','21                 time =',walltime(21)
   print '(a,1x,f13.5)','22                 time =',walltime(22)
   print '(a,1x,f13.5)','23                 time =',walltime(23)
   print '(a,1x,f13.5)','turbine_deposit    time =',walltime(24)
   print '(a,1x,f13.5)','MPI halo communica time =',walltime(25)
   print '(a,1x,f13.5)','buoyancy force     time =',walltime(28)
   print '(a,1x,f13.5)','tracer             time =',walltime(30)
   print '(a,1x,f13.5)','Total wall time    time =',sum(walltime(1:nrtimes))
   print '(a,1x,f13.5)','Total TS routines  time =',sum(walltime(1:nrtimes))-walltime(1)-walltime(14)-walltime(16)
end subroutine

function wtime()
   implicit none
   integer, parameter :: rk = kind(1.0D+00)
   real(kind=rk) :: wtime
#ifdef _CUDA
   integer loc_istat
   loc_istat = cudaDeviceSynchronize()
#endif
   call system_clock(clock_reading, clock_rate, clock_max)
   wtime = real(clock_reading, kind = rk)! / real(clock_rate, kind = rk)
end function

end module

