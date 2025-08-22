module m_averaging_kernelfin2
contains
#ifdef _CUDA
   attributes(global)&
#endif
   subroutine averaging_kernelfin2(nx, ny, iave, uxave, vxave, wxave, uxave2, vxave2, wxave2, Tix, uini )
#ifdef _CUDA
   use cudafor
#endif
   implicit none
   integer, value   :: nx,ny,iave
   real,    value   :: uini

   real, intent(inout) :: uxave(nx,ny)
   real, intent(inout) :: vxave(nx,ny)
   real, intent(inout) :: wxave(nx,ny)

   real, intent(inout) :: uxave2(nx,ny)
   real, intent(inout) :: vxave2(nx,ny)
   real, intent(inout) :: wxave2(nx,ny)

   real, intent(inout) :: Tix(nx,ny)

   integer ::  i,j
#ifdef _CUDA
   integer k
#endif

#ifdef _CUDA
   i = threadIdx%x + (blockIdx%x - 1) * blockDim%x
   j = threadIdx%y + (blockIdx%y - 1) * blockDim%y
   k=1
   if ((j > ny).or.(i > nx)) return
#else
!$OMP PARALLEL DO PRIVATE(i,j) SHARED(nx, ny, uxave, vxave,wxave, uxave2, vxave2, wxave2, Tix, iave, uini )
   do j=1,ny
   do i=1,nx
#endif
      uxave(i,j)=uxave(i,j)/real(iave)
      vxave(i,j)=vxave(i,j)/real(iave)
      wxave(i,j)=wxave(i,j)/real(iave)

      uxave2(i,j)=uxave2(i,j)/real(iave)
      vxave2(i,j)=vxave2(i,j)/real(iave)
      wxave2(i,j)=wxave2(i,j)/real(iave)


      Tix(i,j)=uxave2(i,j)-uxave(i,j)**2 + vxave2(i,j)-vxave(i,j)**2 + wxave2(i,j)-wxave(i,j)**2
      Tix(i,j)=sqrt(Tix(i,j)/3.0)/uini

#ifndef _CUDA
    enddo
    enddo
!$OMP END PARALLEL DO
#endif

end subroutine
end module
