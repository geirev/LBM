module m_averaging_sec_kernel1
contains
#ifdef _CUDA
   attributes(global)&
#endif
   subroutine averaging_sec_kernel1(jdim, kdim, nrsec, ja, ka, iseci, u, v, w, uave, vave, wave, uave2, vave2, wave2 )
                            
#ifdef _CUDA
   use cudafor
#endif
   use mod_dimensions, only : nx,ny,nz
   implicit none
   integer, value      :: jdim,kdim,ja, ka, nrsec
   real, intent(in)    :: u(nx,ny,nz)
   real, intent(in)    :: v(nx,ny,nz)
   real, intent(in)    :: w(nx,ny,nz)


   real, intent(inout) :: uave(nrsec,jdim,kdim)
   real, intent(inout) :: vave(nrsec,jdim,kdim)
   real, intent(inout) :: wave(nrsec,jdim,kdim)

   real, intent(inout) :: uave2(nrsec,jdim,kdim)
   real, intent(inout) :: vave2(nrsec,jdim,kdim)
   real, intent(inout) :: wave2(nrsec,jdim,kdim)

   integer, intent(in) :: iseci(nrsec)

   integer :: isec, j, k, ii

#ifdef _CUDA
   isec = threadIdx%x + (blockIdx%x - 1) * blockDim%x
   j = threadIdx%y + (blockIdx%y - 1) * blockDim%y
   k = threadIdx%z + (blockIdx%z - 1) * blockDim%z
   if (k > kdim) return
   if (j > jdim) return
   if (isec > nrsec) return
#else
!$OMP PARALLEL DO PRIVATE(isec,j,k,ii) SHARED(jdim,kdim,ja, ka, nrsec,uave,vave,wave,uave2,vave2,wave2,u,v,w,iseci)
   do k=1,kdim
   do j=1,jdim
   do isec=1,nrsec
#endif
      ii=iseci(isec)

      uave(isec,j,k)=uave(isec,j,k)+u(ii,ja-1+j,ka-1+k)
      vave(isec,j,k)=vave(isec,j,k)+v(ii,ja-1+j,ka-1+k)
      wave(isec,j,k)=wave(isec,j,k)+w(ii,ja-1+j,ka-1+k)

      uave2(isec,j,k)=uave2(isec,j,k)+u(ii,ja-1+j,ka-1+k)**2
      vave2(isec,j,k)=vave2(isec,j,k)+v(ii,ja-1+j,ka-1+k)**2
      wave2(isec,j,k)=wave2(isec,j,k)+w(ii,ja-1+j,ka-1+k)**2
#ifndef _CUDA
    enddo
    enddo
    enddo
!$OMP END PARALLEL DO
#endif

end subroutine
end module
