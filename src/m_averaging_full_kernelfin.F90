module m_averaging_full_kernelfin
contains
#ifdef _CUDA
   attributes(global)&
#endif
   subroutine averaging_full_kernelfin(uave, vave, wave, uave2, vave2, wave2, Ti, uini, iave)
#ifdef _CUDA
   use cudafor
#endif
   use mod_dimensions, only : nx,ny,nz
   implicit none
   integer, value      :: iave
   real,    value      :: uini

   real, intent(inout) ::  uave(nx,ny,nz)
   real, intent(inout) ::  vave(nx,ny,nz)
   real, intent(inout) ::  wave(nx,ny,nz)

   real, intent(inout) :: uave2(nx,ny,nz)
   real, intent(inout) :: vave2(nx,ny,nz)
   real, intent(inout) :: wave2(nx,ny,nz)

   real, intent(out)   ::    Ti(nx,ny,nz)

   integer :: i, j, k

#ifdef _CUDA
   i = threadIdx%x + (blockIdx%x - 1) * blockDim%x
   j = threadIdx%y + (blockIdx%y - 1) * blockDim%y
   k = threadIdx%z + (blockIdx%z - 1) * blockDim%z
   if (k > nz) return
   if (j > ny) return
   if (i > nx) return
#else
!$OMP PARALLEL DO PRIVATE(i,j,k) SHARED(uave,vave,wave,uave2,vave2,wave2,Ti,iave,uini)
   do k=1,nz
   do j=1,ny
   do i=1,nx
#endif
      uave(i,j,k)=uave(i,j,k)/real(iave)
      vave(i,j,k)=vave(i,j,k)/real(iave)
      wave(i,j,k)=wave(i,j,k)/real(iave)

      uave2(i,j,k)=uave2(i,j,k)/real(iave)
      vave2(i,j,k)=vave2(i,j,k)/real(iave)
      wave2(i,j,k)=wave2(i,j,k)/real(iave)

      Ti(i,j,k)=uave2(i,j,k)-uave(i,j,k)**2 + vave2(i,j,k)-vave(i,j,k)**2 + wave2(i,j,k)-wave(i,j,k)**2
      Ti(i,j,k)=sqrt(Ti(i,j,k)/3.0)  !/uini

      !uave(i,j,k)=uave(i,j,k)/uini
      !vave(i,j,k)=vave(i,j,k)/uini
      !wave(i,j,k)=wave(i,j,k)/uini
#ifndef _CUDA
    enddo
    enddo
    enddo
!$OMP END PARALLEL DO
#endif

end subroutine
end module
