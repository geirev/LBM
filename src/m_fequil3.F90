module m_fequil3
contains

subroutine fequil3(feq, rho, u, v, w, A0_2, A0_3, vel)
   use mod_dimensions
   use mod_D3Q27setup
   use m_readinfile
   use m_wtime

   implicit none
   real, intent(in)      :: rho(nx,ny,nz)
   real, intent(in)      :: u(nx,ny,nz)
   real, intent(in)      :: v(nx,ny,nz)
   real, intent(in)      :: w(nx,ny,nz)
   real, intent(out)     :: feq(nl,0:nx+1,0:ny+1,0:nz+1)
#ifdef _CUDA
   attributes(device) :: rho
   attributes(device) :: u
   attributes(device) :: v
   attributes(device) :: w
   attributes(device) :: feq
#endif

   real, intent(out)   :: A0_2(3,3,nx,ny,nz)
   real, intent(out)   :: A0_3(3,3,3,nx,ny,nz)
   real, intent(out)   :: vel(3,nx,ny,nz)
#ifdef _CUDA
   attributes(device) :: A0_2
   attributes(device) :: A0_3
   attributes(device) :: vel
#endif

   integer :: i, j, k, l, p, q, r, ia

   real, parameter :: inv1cs2 = 1.0/(cs2)
   real, parameter :: inv2cs4 = 1.0/(2.0*cs4)
   real, parameter :: inv2cs6 = 1.0/(2.0*cs6)
   real, parameter :: inv6cs6 = 1.0/(6.0*cs6)
   integer, parameter :: icpu=4
   call cpustart()

!  Compute A0_2, A0_3, and vel
! Loop over grid
#ifdef _CUDA
!$cuf kernel do(3) <<<*,*>>>
#else
!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(i,j,k,l,p,q,r)  SHARED(rho,u,v,w,vel,A0_2,A0_3)
#endif
   do k=1,nz
      do j=1,ny
         do i=1,nx

            vel(1,i,j,k)=u(i,j,k)
            vel(2,i,j,k)=v(i,j,k)
            vel(3,i,j,k)=w(i,j,k)

! A0_2 and A0_3 from \citet{fen21a} (following Eq. 32)
            do q=1,3
            do p=1,3
               A0_2(p,q,i,j,k)=rho(i,j,k)*vel(p,i,j,k)*vel(q,i,j,k)
            enddo
            enddo

            do r=1,3
            do q=1,3
            do p=1,3
               A0_3(p,q,r,i,j,k)=rho(i,j,k)*vel(p,i,j,k)*vel(q,i,j,k)*vel(r,i,j,k)
            enddo
            enddo
            enddo
         enddo
      enddo
   enddo
#ifndef _CUDA
!$OMP END PARALLEL DO
#endif

!  Compute 1st order eqilibrium distribution
#ifdef _CUDA
!$cuf kernel do(3) <<<*,*>>>
#else
!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(i,j,k) SHARED(feq, rho, vel)
#endif
   do k=1,nz
      do j=1,ny
         do i=1,nx
! Equilibrium distribution \citet{fen21a} Eq. (32) or jac18a eq (27)
            feq( 1,i,j,k) = rho(i,j,k) * (cs2                                                ) / cs2
            feq( 2,i,j,k) = rho(i,j,k) * (cs2  + vel(1,i,j,k)                                ) / cs2
            feq( 3,i,j,k) = rho(i,j,k) * (cs2  - vel(1,i,j,k)                                ) / cs2
            feq( 4,i,j,k) = rho(i,j,k) * (cs2                  + vel(2,i,j,k)                ) / cs2
            feq( 5,i,j,k) = rho(i,j,k) * (cs2                  - vel(2,i,j,k)                ) / cs2
            feq( 6,i,j,k) = rho(i,j,k) * (cs2                                  - vel(3,i,j,k)) / cs2
            feq( 7,i,j,k) = rho(i,j,k) * (cs2                                  + vel(3,i,j,k)) / cs2
            feq( 8,i,j,k) = rho(i,j,k) * (cs2  + vel(1,i,j,k)  + vel(2,i,j,k)                ) / cs2
            feq( 9,i,j,k) = rho(i,j,k) * (cs2  - vel(1,i,j,k)  - vel(2,i,j,k)                ) / cs2
            feq(10,i,j,k) = rho(i,j,k) * (cs2  + vel(1,i,j,k)  - vel(2,i,j,k)                ) / cs2
            feq(11,i,j,k) = rho(i,j,k) * (cs2  - vel(1,i,j,k)  + vel(2,i,j,k)                ) / cs2
            feq(12,i,j,k) = rho(i,j,k) * (cs2  - vel(1,i,j,k)                  - vel(3,i,j,k)) / cs2
            feq(13,i,j,k) = rho(i,j,k) * (cs2  + vel(1,i,j,k)                  + vel(3,i,j,k)) / cs2
            feq(14,i,j,k) = rho(i,j,k) * (cs2                  + vel(2,i,j,k)  + vel(3,i,j,k)) / cs2
            feq(15,i,j,k) = rho(i,j,k) * (cs2                  - vel(2,i,j,k)  - vel(3,i,j,k)) / cs2
            feq(16,i,j,k) = rho(i,j,k) * (cs2  - vel(1,i,j,k)                  + vel(3,i,j,k)) / cs2
            feq(17,i,j,k) = rho(i,j,k) * (cs2  + vel(1,i,j,k)                  - vel(3,i,j,k)) / cs2
            feq(18,i,j,k) = rho(i,j,k) * (cs2                  - vel(2,i,j,k)  + vel(3,i,j,k)) / cs2
            feq(19,i,j,k) = rho(i,j,k) * (cs2                  + vel(2,i,j,k)  - vel(3,i,j,k)) / cs2
            feq(20,i,j,k) = rho(i,j,k) * (cs2  - vel(1,i,j,k)  + vel(2,i,j,k)  + vel(3,i,j,k)) / cs2
            feq(21,i,j,k) = rho(i,j,k) * (cs2  + vel(1,i,j,k)  - vel(2,i,j,k)  - vel(3,i,j,k)) / cs2
            feq(22,i,j,k) = rho(i,j,k) * (cs2  - vel(1,i,j,k)  - vel(2,i,j,k)  - vel(3,i,j,k)) / cs2
            feq(23,i,j,k) = rho(i,j,k) * (cs2  + vel(1,i,j,k)  + vel(2,i,j,k)  + vel(3,i,j,k)) / cs2
            feq(24,i,j,k) = rho(i,j,k) * (cs2  + vel(1,i,j,k)  + vel(2,i,j,k)  - vel(3,i,j,k)) / cs2
            feq(25,i,j,k) = rho(i,j,k) * (cs2  - vel(1,i,j,k)  - vel(2,i,j,k)  + vel(3,i,j,k)) / cs2
            feq(26,i,j,k) = rho(i,j,k) * (cs2  - vel(1,i,j,k)  + vel(2,i,j,k)  - vel(3,i,j,k)) / cs2
            feq(27,i,j,k) = rho(i,j,k) * (cs2  + vel(1,i,j,k)  - vel(2,i,j,k)  + vel(3,i,j,k)) / cs2
         enddo
      enddo
   enddo
#ifndef _CUDA
!$OMP END PARALLEL DO
#endif

!  Compute 2st order eqilibrium distribution
#ifdef _CUDA
!$cuf kernel do(2) <<<*,*>>>
#else
!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(i, j, k, l, p, q) SHARED(feq, H2, A0_2)
#endif
   do k=1,nz
      do j=1,ny
         do i=1,nx
            do l=1,nl
               do p=1,3
               do q=1,3
                  feq(l,i,j,k)=feq(l,i,j,k) + H2(p,q,l)*A0_2(p,q,i,j,k)*inv2cs4
               enddo
               enddo
            enddo
         enddo
      enddo
   enddo
#ifndef _CUDA
!$OMP END PARALLEL DO
#endif

!  the above identically recovers the BGK equilibrium, now we add third order contributions
   if (ibgk == 3) then
#ifdef _CUDA
!$cuf kernel do(2) <<<*,*>>>
#else
!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(i,j,k,l,p,q,r)  SHARED(feq,H3,A0_3)
#endif
      do k=1,nz
         do j=1,ny
            do i=1,nx
               do l=1,nl
                  do p=1,3
                  do q=1,3
                  do r=1,3
                     feq(l,i,j,k)=feq(l,i,j,k) + H3(p,q,r,l)*A0_3(p,q,r,i,j,k)*inv6cs6
                  enddo
                  enddo
                  enddo
               enddo
            enddo
         enddo
      enddo
#ifndef _CUDA
!$OMP END PARALLEL DO
#endif
   endif

!  Scale with the weights
#ifdef _CUDA
!$cuf kernel do(2) <<<*,*>>>
#else
!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(i,j,k,l) SHARED(feq, weights)
#endif
      do k=1,nz
         do j=1,ny
            do i=1,nx
               do l=1,nl
                  feq(l,i,j,k)= weights(l)*feq(l,i,j,k)
               enddo
            enddo
         enddo
      enddo
#ifndef _CUDA
!$OMP END PARALLEL DO
#endif

   call cpufinish(icpu)


end subroutine

end module
