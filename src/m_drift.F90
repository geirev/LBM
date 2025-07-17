module m_drift
contains
subroutine drift(f,feq,it)
! f enter routine in feq following collisions
! f is returned in f
   use mod_dimensions
   use mod_D3Q27setup
   use m_wtime
   use m_drift_kernel
   implicit none
   real, intent(out) :: f(nl,0:nx+1,0:ny+1,0:nz+1)
   real, intent(in)  :: feq(nl,0:nx+1,0:ny+1,0:nz+1)
   integer it
#ifdef _CUDA
   attributes(device) :: f
   attributes(device) :: feq
#endif
   integer i,j,k,l
   integer, parameter :: icpu=12
   integer :: tx, ty, tz, bx, by, bz
   call cpustart()
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!@cuf istat = cudaDeviceSynchronize()
   t0 = wallclock()


#ifdef _CUDA
   tx=8; bx=(nx+2+tx-1)/tx
   ty=8; by=(ny+2+ty-1)/ty
   tz=8; bz=(nz+2+tz-1)/tz
#endif
   call drift_kernel&
#ifdef _CUDA
        &<<<dim3(bx,by,bz), dim3(tx,ty,tz)>>>&
#endif
        &(f,feq, nx+2, ny+2, nz+2, nl, cxs, cys, czs)





!!  #ifdef _CUDA
!!  !$cuf kernel do(3) <<<*,*>>>
!!  #else
!!  !$OMP PARALLEL DO PRIVATE(i,j,k,l) SHARED(f, feq, cxs,cys,czs)
!!  #endif
!!     do k=1,nz
!!     do j=1,ny
!!     do i=1,nx
!!        do l=1,nl
!!           f(l,i,j,k) = feq(l,i-cxs(l),j-cys(l),k-czs(l))
!!        enddo
!!      enddo
!!      enddo
!!      enddo
#ifndef _CUDA
!$OMP END PARALLEL DO
#endif


!@cuf istat = cudaDeviceSynchronize()
   t1 = wallclock(); walltimelocal(21)=walltimelocal(21)+t1-t0







   call cpufinish(icpu)
   if (it==999) then
      print *
      do j=21,21
         print '(a,i3,g13.5)','drift       :',j,walltimelocal(j)
      enddo
   endif

end subroutine
end module
