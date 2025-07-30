module m_hermite_polynomials
contains
subroutine hermite_polynomials()
   use mod_dimensions
   use mod_D3Q27setup
#ifdef _CUDA
   use cudafor
#endif
   implicit none
   integer :: l, p, q, r,ierr
   real :: c_h(3,nl)
   real :: delta_h(3,3)

   real, allocatable :: c(:,:)
   real, allocatable :: delta(:,:)
#ifdef _CUDA
   attributes(device) :: c
   attributes(device) :: delta
#endif


   allocate(cxs(nl))
   allocate(cys(nl))
   allocate(czs(nl))
   allocate(bounce(nl))
   allocate(weights(1:nl))

!!  #ifdef _CUDA
!!      ierr = cudaMemcpy(cxs, cxs_h, nl*4, cudaMemcpyHostToDevice)
!!      ierr = cudaMemcpy(cys, cys_h, nl*4, cudaMemcpyHostToDevice)
!!      ierr = cudaMemcpy(czs, czs_h, nl*4, cudaMemcpyHostToDevice)
!!      ierr = cudaMemcpy(bounce, bounce_h, nl*4, cudaMemcpyHostToDevice)
!!  !   ierr = cudaMemcpy(weights, weights_h, nl*sz, cudaMemcpyHostToDevice)
!!      weights=weights_h
!!  #else
    cxs=cxs_h
    cys=cys_h
    czs=czs_h
    bounce=bounce_h
    weights=weights_h
!! #endif


! Build lattice velocity matrix
   allocate(c(3,nl))
   c_h(1,:) = real(cxs_h(:))
   c_h(2,:) = real(cys_h(:))
   c_h(3,:) = real(czs_h(:))
#ifdef _CUDA
   ! ierr = cudaMemcpy(c, c_h, 3*nl*sz, cudaMemcpyHostToDevice)
#endif
   c=c_h

! Kronecker delta
   allocate(delta(3,3))
   delta_h = 0.0
   delta_h(1,1) = 1.0
   delta_h(2,2) = 1.0
   delta_h(3,3) = 1.0
#ifdef _CUDA
    !ierr = cudaMemcpy(delta, delta_h, 3*3*sz, cudaMemcpyHostToDevice); print *,'cpdelta:ierr',ierr
#endif
   delta=delta_h

   allocate(H2(3,3,1:nl))
   allocate(H3(3,3,3,1:nl))
! Compute second-order Hermite polynomials
#ifdef _CUDA
!$cuf kernel do(1) <<<*,*>>>
#endif
   do l = 1, nl
      do p = 1, 3
         do q = 1, 3
            H2(p,q,l) = c(p,l)*c(q,l) - cs2 * delta(p,q)
         enddo
      enddo
   enddo

   ! Compute third-order Hermite polynomials
#ifdef _CUDA
!$cuf kernel do(1) <<<*,*>>>
#endif
   do l = 1, nl
      do p = 1, 3
         do q = 1, 3
            do r = 1, 3
               H3(p,q,r,l) = c(p,l)*c(q,l)*c(r,l) - cs2 * ( &
                               c(p,l)*delta(q,r) + &
                               c(q,l)*delta(p,r) + &
                               c(r,l)*delta(p,q))
            enddo
         enddo
      enddo
   enddo

end subroutine
end module

!   open(10,file='H2n.dat')
!       do l=1,nl
!          write(10,'(i3,a,9f10.5)')l,':',H2(:,:,l)
!       enddo
!   close(10)
!
!
!   open(10,file='H3n.dat')
!      do l=1,nl
!         write(10,'(i3,a,27g12.5)')l,':',H3(:,:,:,l)
!      enddo
!   close(10)

