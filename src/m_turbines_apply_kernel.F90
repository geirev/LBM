module m_turbines_apply_kernel
contains
#ifdef _CUDA
   attributes(global) &
#endif
subroutine turbines_apply_kernel(f,df,tau,ip,n,iforce,nturbines)
   use mod_dimensions,  only : nx,ny,nz
   use mod_D3Q27setup,  only : nl
   use m_turbines_init, only : ieps
   implicit none
   integer, value      :: nturbines
   real, intent(inout) :: f(nl,0:nx+1,0:ny+1,0:nz+1)        ! distribution
   real, intent(in)    :: df(nl,-ieps:ieps,ny,nz,nturbines) ! forcing distributions
   real, intent(in)    :: tau(0:nx+1,0:ny+1,0:nz+1)         ! Tau
   integer, value      :: ip
   integer, value      :: n
   integer, value      :: iforce
   integer i,j,k,l,gx

   gx=1

#ifdef _CUDA
   gx = threadIdx%x + (blockIdx%x - 1) * blockDim%x
   i  = (ip-ieps) + gx -1
   j = threadIdx%y + (blockIdx%y - 1) * blockDim%y
   k = threadIdx%z + (blockIdx%z - 1) * blockDim%z
   if (i > nx .or. j > ny .or. k > nz) return
   if (i < ip-ieps  .or. i > ip+ieps) return
#else
!$OMP PARALLEL DO collapse(3) DEFAULT(NONE)&
!$OMP PRIVATE(i,j,k,l,n)&
!$OMP SHARED(f,df,ip)
      do k=1,nz
      do j=1,ny
      do i=ip-ieps,ip+ieps
#endif
!         if (iforce == 8) fac=(1.0-0.5/tau(i,j,k))
         do l=1,nl
            !f(l,i,j,k) = f(l,i,j,k) + fac*df(l,i-ip,j,k,n)
            f(l,i,j,k) = f(l,i,j,k) + df(l,i-ip,j,k,n)
         enddo
#ifndef _CUDA
      enddo
      enddo
      enddo
#endif
end subroutine
end module
