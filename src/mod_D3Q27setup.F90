module mod_D3Q27setup
   use mod_dimensions
#ifdef _CUDA
   use cudafor
#endif
   implicit none

   real, parameter :: cs2=1/3.0
   real, parameter :: cs4=1/9.0
   real, parameter :: cs6=1/27.0


#ifdef D3Q19
! D3Q19 weights: [rest, 6 faces, 12 edges]
   integer, parameter :: cxs_h(1:nl)     = [0, 1,-1, 0, 0, 0, 0, 1,-1, 1,-1,-1, 1, 0, 0,-1, 1, 0, 0]
   integer, parameter :: cys_h(1:nl)     = [0, 0, 0, 1,-1, 0, 0, 1,-1,-1, 1, 0, 0, 1,-1, 0, 0,-1, 1]
   integer, parameter :: czs_h(1:nl)     = [0, 0, 0, 0, 0,-1, 1, 0, 0, 0, 0,-1, 1, 1,-1, 1,-1, 1,-1]
   integer, parameter :: bounce_h(1:nl)  = [1, 3, 2, 5, 4, 7, 6, 9, 8,11,10,13,12,15,14,17,16,19,18]
   real,    parameter :: weights_h(1:nl) = [1.0/3.0, &
                                            1.0/18.0, 1.0/18.0, 1.0/18.0, 1.0/18.0, 1.0/18.0, 1.0/18.0, &
                                            1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0, &
                                            1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0 ]

#else
! Define velocity vectors
!                                           1  2  3  4  5  6  7  8  9  0  1  2  3  4  5  6  7  8  9  0  1  2  3  4  5  6  7
   integer, parameter :: cxs_h(1:nl)     = [0, 1,-1, 0, 0, 0, 0, 1,-1, 1,-1,-1, 1, 0, 0,-1, 1, 0, 0,-1, 1,-1, 1, 1,-1,-1, 1]
   integer, parameter :: cys_h(1:nl)     = [0, 0, 0, 1,-1, 0, 0, 1,-1,-1, 1, 0, 0, 1,-1, 0, 0,-1, 1, 1,-1,-1, 1, 1,-1, 1,-1]
   integer, parameter :: czs_h(1:nl)     = [0, 0, 0, 0, 0,-1, 1, 0, 0, 0, 0,-1, 1, 1,-1, 1,-1, 1,-1, 1,-1,-1, 1,-1, 1,-1, 1]
   integer, parameter :: bounce_h(1:nl)  = [1, 3, 2, 5, 4, 7, 6, 9, 8,11,10,13,12,15,14,17,16,19,18,21,20,23,22,25,24,27,26]
   real,    parameter :: weights_h(1:nl) = [8.0/27.0, &
                                            2.0/27.0, 2.0/27.0, 2.0/27.0, 2.0/27.0, 2.0/27.0, 2.0/27.0, &
                                            1.0/54.0, 1.0/54.0, 1.0/54.0, 1.0/54.0, 1.0/54.0, 1.0/54.0, &
                                            1.0/54.0, 1.0/54.0, 1.0/54.0, 1.0/54.0, 1.0/54.0, 1.0/54.0, &
                                            1.0/216.0, 1.0/216.0, 1.0/216.0, 1.0/216.0, 1.0/216.0, 1.0/216.0, 1.0/216.0, 1.0/216.0]
#endif


   integer, allocatable :: cxs(:)
   integer, allocatable :: cys(:)
   integer, allocatable :: czs(:)
   integer, allocatable :: bounce(:)
   real,  dimension(:)      , allocatable :: weights
   real,  dimension(:,:,:)  , allocatable :: H2
   real,  dimension(:,:,:,:), allocatable :: H3
#ifdef _CUDA
   attributes(device) :: cxs
   attributes(device) :: cys
   attributes(device) :: czs
   attributes(device) :: bounce
   attributes(device) :: weights
   attributes(device) :: H2
   attributes(device) :: H3
#endif

contains

subroutine hermite_polynomials()
!   use mod_dimensions
!   use mod_D3Q27setup
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



end module
