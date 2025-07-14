module m_regularization
contains

subroutine regularization(f, feq, u, v, w, A1_2, A1_3, vel)
   use mod_dimensions
   use mod_D3Q27setup
   use m_ablim
   use m_readinfile
!   use m_regnoneq
!   use m_reggetA1
!   use m_reg2ord
!   use m_reg3ord
!   use m_fequilscale
   use m_wtime
   use cublas
   implicit none
   real, intent(in)       :: u(nx,ny,nz)
   real, intent(in)       :: v(nx,ny,nz)
   real, intent(in)       :: w(nx,ny,nz)
   real, intent(in)       :: feq(nl,0:nx+1,0:ny+1,0:nz+1)
   real, intent(inout)    :: f(nl,0:nx+1,0:ny+1,0:nz+1)
#ifdef _CUDA
   attributes(device) :: u
   attributes(device) :: v
   attributes(device) :: w
   attributes(device) :: f
   attributes(device) :: feq
#endif

   real, intent(out)   :: A1_2(3,3,nx,ny,nz)
   real, intent(out)   :: A1_3(3,3,3,nx,ny,nz)
   real, intent(out)   :: vel(1:3,nx,ny,nz)
#ifdef _CUDA
   attributes(device) :: A1_2
   attributes(device) :: A1_3
   attributes(device) :: vel
#endif

   integer :: i, j, k, l, p, q, r

   real, parameter :: inv2cs4 = 1.0/(2.0*cs4)
   real, parameter :: inv2cs6 = 1.0/(2.0*cs6)
   real, parameter :: inv6cs6 = 1.0/(6.0*cs6)
   integer, parameter :: icpu=5
   call cpustart()


! Computing non-equilibrium distribution defined in \citet{fen21a} between Eqs (32) and (33)
! scaling by the weights
#ifdef _CUDA
!$cuf kernel do(2) <<<*,*>>>
#else
!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(i,j,k,l) SHARED(feq, f)
#endif
      do k=1,nz
         do j=1,ny
            do i=1,nx
               do l=1,nl
                  f(l,i,j,k)=f(l,i,j,k)-feq(l,i,j,k)
               enddo
            enddo
         enddo
      enddo
#ifndef _CUDA
!$OMP END PARALLEL DO
#endif



!   if (ihrr == 1) then
! projecting non-equilibrium distribution on the Hermitian polynomials for regularization
#ifdef _CUDA
!$cuf kernel do(2) <<<*,*>>>
#else
!$OMP PARALLEL DO collapse(1) DEFAULT(NONE) PRIVATE(i, j, k, l, p, q, r ) &
!$OMP&                           SHARED(f, u, v, w,  H2, H3, vel, A1_2, A1_3)
#endif
      do k=1,nz
      do j=1,ny
      do i=1,nx
         vel(1,i,j,k)=u(i,j,k)
         vel(2,i,j,k)=v(i,j,k)
         vel(3,i,j,k)=w(i,j,k)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Eq (11) from  Jacob 2018 is identical to the 33a from Feng (2021)
!! Used for regularization and turbulence calculation
!              call dgemv('n', 9,27,1.0,H2, 9,f(1,i,j,k),1,0.0,A1_2,1)

         do q=1,3
         do p=1,3
            A1_2(p,q,i,j,k)=0.0
         enddo
         enddo
         do l=1,nl
            do q=1,3
            do p=1,3
               A1_2(p,q,i,j,k) = A1_2(p,q,i,j,k) + H2(p,q,l)*f(l,i,j,k)
            enddo
            enddo
         enddo

! A1_3(HRR) from \citet{fen21a}, as defined after Eq. (34). Proof by Malaspinas 2015, Appendix B
         do r=1,3
         do q=1,3
         do p=1,3
            A1_3(p,q,r,i,j,k)=vel(p,i,j,k)*A1_2(q,r,i,j,k) +vel(q,i,j,k)*A1_2(r,p,i,j,k) +  vel(r,i,j,k)*A1_2(p,q,i,j,k)
         enddo
         enddo
         enddo

      enddo
      enddo
      enddo
#ifndef _CUDA
!$OMP END PARALLEL DO
#endif
! Rfneq from \citet{fen21a}, as defined in Eq. (34)
!              call dgemv('t', 9,27,inv2cs4,H2, 9,A1_2,1,0.0,f(1,i,j,k),1)
!              call dgemv('t',27,27,inv6cs6,H3,27,A1_3,1,1.0,f(1,i,j,k),1)
!              f(:,i,j,k)=weights(:)*f(:,i,j,k)


#ifdef _CUDA
!$cuf kernel do(2) <<<*,*>>>
#else
!$OMP PARALLEL DO collapse(3) DEFAULT(NONE) PRIVATE(i, j, k, l, p, q, r ) SHARED(f, H2, A1_2)
#endif
      do k=1,nz
      do j=1,ny
      do i=1,nx
         do l=1,nl
            f(l,i,j,k)=0.0
            do q=1,3
            do p=1,3
               f(l,i,j,k)=f(l,i,j,k) + H2(p,q,l)*A1_2(p,q,i,j,k)*inv2cs4
            enddo
            enddo
         enddo
      enddo
      enddo
      enddo
#ifndef _CUDA
!$OMP END PARALLEL DO
#endif

#ifdef _CUDA
!$cuf kernel do(2) <<<*,*>>>
#else
!$OMP PARALLEL DO collapse(1) DEFAULT(NONE) PRIVATE(i, j, k, l, p, q, r ) SHARED(f, H3, A1_3)
#endif
      do k=1,nz
      do j=1,ny
      do i=1,nx
         do l=1,nl
            do r=1,3
            do q=1,3
            do p=1,3
               f(l,i,j,k)=f(l,i,j,k) + H3(p,q,r,l)*A1_3(p,q,r,i,j,k)*inv6cs6
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

! scaling by the weights
#ifdef _CUDA
!$cuf kernel do(2) <<<*,*>>>
#else
!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(i,j,k,l) SHARED(f, weights)
#endif
      do k=1,nz
         do j=1,ny
            do i=1,nx
               do l=1,nl
                  f(l,i,j,k)= weights(l)*f(l,i,j,k)
               enddo
            enddo
         enddo
      enddo
#ifndef _CUDA
!$OMP END PARALLEL DO
#endif
!   endif

   call cpufinish(icpu)

end subroutine

end module
