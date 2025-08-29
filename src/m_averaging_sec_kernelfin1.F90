module m_averaging_sec_kernelfin1
contains
#ifdef _CUDA
   attributes(global)&
#endif
   subroutine averaging_sec_kernelfin1(jdim, kdim, nrsec, iave,  uave, vave, wave, uave2, vave2, wave2, Ti, uini )
#ifdef _CUDA
   use cudafor
#endif
   implicit none
   integer, value      :: jdim,kdim,nrsec,iave
   real,    value      :: uini

   real, intent(inout) ::  uave(nrsec,jdim,kdim)
   real, intent(inout) ::  vave(nrsec,jdim,kdim)
   real, intent(inout) ::  wave(nrsec,jdim,kdim)

   real, intent(inout) :: uave2(nrsec,jdim,kdim)
   real, intent(inout) :: vave2(nrsec,jdim,kdim)
   real, intent(inout) :: wave2(nrsec,jdim,kdim)

   real, intent(out)   ::    Ti(nrsec,jdim,kdim)

   integer :: isec, j, k

#ifdef _CUDA
   isec = threadIdx%x + (blockIdx%x - 1) * blockDim%x
   j = threadIdx%y + (blockIdx%y - 1) * blockDim%y
   k = threadIdx%z + (blockIdx%z - 1) * blockDim%z
   if (k > kdim) return
   if (j > jdim) return
   if (isec > nrsec) return
#else
!$OMP PARALLEL DO PRIVATE(isec,j,k) SHARED(jdim,kdim,nrsec,uave,vave,wave,uave2,vave2,wave2,Ti,iave,uini)
   do k=1,kdim
   do j=1,jdim
   do isec=1,nrsec
#endif
      uave(isec,j,k)=uave(isec,j,k)/real(iave)
      vave(isec,j,k)=vave(isec,j,k)/real(iave)
      wave(isec,j,k)=wave(isec,j,k)/real(iave)

      uave2(isec,j,k)=uave2(isec,j,k)/real(iave)
      vave2(isec,j,k)=vave2(isec,j,k)/real(iave)
      wave2(isec,j,k)=wave2(isec,j,k)/real(iave)

      Ti(isec,j,k)=uave2(isec,j,k)-uave(isec,j,k)**2 + vave2(isec,j,k)-vave(isec,j,k)**2 + wave2(isec,j,k)-wave(isec,j,k)**2
      Ti(isec,j,k)=sqrt(Ti(isec,j,k)/3.0)/uini

      uave(isec,j,k)=uave(isec,j,k)/uini
      vave(isec,j,k)=vave(isec,j,k)/uini
      wave(isec,j,k)=wave(isec,j,k)/uini
#ifndef _CUDA
    enddo
    enddo
    enddo
!$OMP END PARALLEL DO
#endif

end subroutine
end module
