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

   real, intent(inout) ::  uave(0:nx+1,0:ny+1,0:nz+1)
   real, intent(inout) ::  vave(0:nx+1,0:ny+1,0:nz+1)
   real, intent(inout) ::  wave(0:nx+1,0:ny+1,0:nz+1)

   real, intent(inout) :: uave2(0:nx+1,0:ny+1,0:nz+1)
   real, intent(inout) :: vave2(0:nx+1,0:ny+1,0:nz+1)
   real, intent(inout) :: wave2(0:nx+1,0:ny+1,0:nz+1)

   real, intent(inout) ::    Ti(0:nx+1,0:ny+1,0:nz+1)
   real Ti2

   integer :: i, j, k

   if (iave == 0) then
      print *,'iave=0'
      return
   endif

#ifdef _CUDA
   i = threadIdx%x + (blockIdx%x - 1) * blockDim%x
   j = threadIdx%y + (blockIdx%y - 1) * blockDim%y
   k = threadIdx%z + (blockIdx%z - 1) * blockDim%z
   if (k > nz) return
   if (j > ny) return
   if (i > nx) return
#else
!$OMP PARALLEL DO PRIVATE(i,j,k,Ti2) SHARED(uave,vave,wave,uave2,vave2,wave2,Ti,iave,uini)
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

!      Ti(i,j,k)=uave2(i,j,k)-uave(i,j,k)**2 + vave2(i,j,k)-vave(i,j,k)**2 + wave2(i,j,k)-wave(i,j,k)**2
!      Ti(i,j,k)=sqrt(Ti(i,j,k)/3.0)  !/uini

      Ti2 = ( &
          max(uave2(i,j,k) - uave(i,j,k)**2, 0.0) + &
          max(vave2(i,j,k) - vave(i,j,k)**2, 0.0) + &
          max(wave2(i,j,k) - wave(i,j,k)**2, 0.0) ) / 3.0

         Ti(i,j,k) = sqrt(max(Ti2, 0.0))

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
