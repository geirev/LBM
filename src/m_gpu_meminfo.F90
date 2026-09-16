module m_gpu_meminfo
#ifdef _CUDA
  use cudafor
#endif
  implicit none
contains
#ifdef _CUDA
   subroutine gpu_meminfo(msg)
   use mod_dimensions
   implicit none
   character(len=*), intent(in) :: msg
   integer(8) :: free_mem, total_mem
   real(8) :: free_gb, total_gb, used_gb
   integer istat

   istat=cudaMemGetInfo(free_mem, total_mem)

   free_gb = real(free_mem) / (1024.0d0**3)
   total_gb = real(total_mem) / (1024.0d0**3)
   used_gb = total_gb - free_gb

   print *, "-------------------------------"
   print '(a,a)',          "GPU Memory Info - ", trim(msg)

   print '(a,3i6,f10.2,a)', &
     "Global grid : ", nx,nyg,nz, &
     real(nx,8)*real(nyg,8)*real(nz,8)/1.0d6, " Million"

   print '(a,3i6,f10.2,a)', &
     "Local tile  : ", nx,ny,nz, &
     real(nx,8)*real(ny,8)*real(nz,8)/1.0d6, " Million"

   print '(a,f8.2)',       "Total GPU Memory (GB): ", total_gb
   print '(a,f8.2)',       "Used  GPU Memory (GB): ", used_gb
   print '(a,f8.2)',       "Free  GPU Memory (GB): ", free_gb
   print '(a,3i5,f8.2,a)', "Number of cells      : ", nx,nyg,nz,real(nx*nyg*nz)/1000000.0,' Million'
   print *, "-------------------------------"
end subroutine
#endif
end module
