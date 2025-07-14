module m_drift
contains
subroutine drift(f,feq)
! f enter routine in feq following collisions
! f is returned in f
   use mod_dimensions
   use mod_D3Q27setup
   use m_wtime
   implicit none
   real, intent(out) :: f(nl,0:nx+1,0:ny+1,0:nz+1)
   real, intent(in)  :: feq(nl,0:nx+1,0:ny+1,0:nz+1)
#ifdef _CUDA
   attributes(device) :: f
   attributes(device) :: feq
#endif
   integer i,j,k
   integer, parameter :: icpu=12
   call cpustart()

#ifdef _CUDA
!$cuf kernel do(3) <<<*,*>>>
#endif
   do k=1,nz
   do j=1,ny
#ifndef _CUDA
!$OMP PARALLEL DO COLLAPSE(1) PRIVATE(i) SHARED(f, feq,j,k)
#endif
   do i=1,nx
!      do l=1,nl
!         f(l,i,j,k) = feq(l,i-cxs(l),j-cys(l),k-czs(l))
!      enddo
         f( 1,i,j,k) = feq( 1,i  ,j  ,k  )
         f( 2,i,j,k) = feq( 2,i-1,j  ,k  )
         f( 3,i,j,k) = feq( 3,i+1,j  ,k  )
         f( 4,i,j,k) = feq( 4,i  ,j-1,k  )
         f( 5,i,j,k) = feq( 5,i  ,j+1,k  )
         f( 6,i,j,k) = feq( 6,i  ,j  ,k+1)
         f( 7,i,j,k) = feq( 7,i  ,j  ,k-1)
         f( 8,i,j,k) = feq( 8,i-1,j-1,k  )
         f( 9,i,j,k) = feq( 9,i+1,j+1,k  )
         f(10,i,j,k) = feq(10,i-1,j+1,k  )
         f(11,i,j,k) = feq(11,i+1,j-1,k  )
         f(12,i,j,k) = feq(12,i+1,j  ,k+1)
         f(13,i,j,k) = feq(13,i-1,j  ,k-1)
         f(14,i,j,k) = feq(14,i  ,j-1,k-1)
         f(15,i,j,k) = feq(15,i  ,j+1,k+1)
         f(16,i,j,k) = feq(16,i+1,j  ,k-1)
         f(17,i,j,k) = feq(17,i-1,j  ,k+1)
         f(18,i,j,k) = feq(18,i  ,j+1,k-1)
         f(19,i,j,k) = feq(19,i  ,j-1,k+1)
         f(20,i,j,k) = feq(20,i+1,j-1,k-1)
         f(21,i,j,k) = feq(21,i-1,j+1,k+1)
         f(22,i,j,k) = feq(22,i+1,j+1,k+1)
         f(23,i,j,k) = feq(23,i-1,j-1,k-1)
         f(24,i,j,k) = feq(24,i-1,j-1,k+1)
         f(25,i,j,k) = feq(25,i+1,j+1,k-1)
         f(26,i,j,k) = feq(26,i+1,j-1,k+1)
         f(27,i,j,k) = feq(27,i-1,j+1,k-1)
   enddo
#ifndef _CUDA
!$OMP END PARALLEL DO
#endif
   enddo
   enddo

    call cpufinish(icpu)

end subroutine
end module
