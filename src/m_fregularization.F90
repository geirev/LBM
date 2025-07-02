module m_fregularization
contains

subroutine fregularization(f, feq, u, v, w)
   use mod_dimensions
   use mod_D3Q27setup
   use m_ablim
   use m_readinfile
   use m_wtime
   implicit none
   real, intent(in)      :: u(nx,ny,nz)
   real, intent(in)      :: v(nx,ny,nz)
   real, intent(in)      :: w(nx,ny,nz)
   real, intent(in)      :: feq(nl,0:nx+1,0:ny+1,0:nz+1)

   real, intent(inout)   :: f(nl,0:nx+1,0:ny+1,0:nz+1)

   logical, save         :: lfirst=.true.
   real, save            :: H2(3,3,nl)      ! Second order Hermite polynomial
   real, save            :: H3(3,3,3,nl)    ! Third order Hermite polynomial
   real, save            :: c(3,nl)         ! Array storage of cxs, cys, and czs

   real                  :: A1_2(3,3)
   real                  :: A1_3(3,3,3)

   real                  :: delta(1:3, 1:3) = reshape([1.0, 0.0, 0.0, &
                                                       0.0, 1.0, 0.0, &
                                                       0.0, 0.0, 1.0], [3, 3])

   real                  :: vel(1:3)

   integer :: i, j, k, l, p, q, r

   real, parameter :: inv2cs4 = 1.0/(2.0*cs4)
   real, parameter :: inv2cs6 = 1.0/(2.0*cs6)
   real, parameter :: inv6cs6 = 1.0/(6.0*cs6)
   integer, parameter :: icpu=12
   call cpustart()

   if (lfirst) then
      c(1,:)=real(cxs(:))
      c(2,:)=real(cys(:))
      c(3,:)=real(czs(:))

! Hermitian polynomials of 2nd order H2
      do l=1,nl
         do q=1,3
         do p=1,3
            H2(p,q,l)=c(p,l)*c(q,l) - cs2*delta(p,q)
         enddo
         enddo
      enddo

!      open(10,file='H2.dat')
!         do l=1,nl
!            write(10,'(i3,a,9f10.5)')l,':',H2(:,:,l)
!         enddo
!      close(10)

! Hermitian polynomials of 3nd order H3 with
! [c_\alpha \delta]_{ijk} = c_{\alpha,i} \delta_{jk} + c_{\alpha,j} \delta_{ik} +c_{\alpha,k} \delta_{ij}
      do l=1,nl
         do r=1,3
         do q=1,3
         do p=1,3
            H3(p,q,r,l)=c(p,l)*c(q,l)*c(r,l) - cs2*(c(p,l)*delta(q,r) + c(q,l)*delta(p,r) +  c(r,l)*delta(p,q))
         enddo
         enddo
         enddo
      enddo

!      open(10,file='H3.dat')
!         do l=1,nl
!            write(10,'(i3,a,27g12.5)')l,':',H3(:,:,:,l)
!         enddo
!      close(10)

      lfirst=.false.
   endif

! Computing non-equilibrium distribution defined in \citet{fen21a} between Eqs (32) and (33)
!$OMP PARALLEL DO collapse(2) DEFAULT(NONE) PRIVATE(i, j, k) &
!$OMP&                           SHARED(feq, f)
   do k=1,nz
      do j=1,ny
         do i=1,nx
            f(:,i,j,k)=f(:,i,j,k)-feq(:,i,j,k)
         enddo
      enddo
   enddo
!$OMP END PARALLEL DO

   if (ihrr == 1) then
! projecting non-equilibrium distribution on the Hermitian polynomials for regularization
      do k=1,nz
         do j=1,ny
!$OMP PARALLEL DO collapse(1) DEFAULT(NONE) PRIVATE(i, l, p, q, r, vel, A1_2, A1_3 ) &
!$OMP&                           SHARED( j, k, f, u, v, w,  weights, H2, c , H3)
            do i=1,nx
                  vel(1)=u(i,j,k)
                  vel(2)=v(i,j,k)
                  vel(3)=w(i,j,k)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Eq (11) from  Jacob 2018 is identical to the 33a from Feng (2021)
! Used for regularization and turbulence calculation
                  call dgemv('n', 9,27,1.0,H2, 9,f(1,i,j,k),1,0.0,A1_2,1)
!                  A1_2=0.0
!                  do l=1,nl
!                     do q=1,3
!                     do p=1,3
!                        A1_2(p,q) = A1_2(p,q) + H2(p,q,l)*f(l,i,j,k)
!!                        A1_2(p,q) = A1_2(p,q) + c(p,l)*c(q,l)*f(l,i,j,k)
!                     enddo
!                     enddo
!                  enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! A1_3(HRR) from \citet{fen21a}, as defined after Eq. (34). Proof by Malaspinas 2015, Appendix B
                    do r=1,3
                    do q=1,3
                    do p=1,3
                       A1_3(p,q,r)=vel(p)*A1_2(q,r) +vel(q)*A1_2(r,p) +  vel(r)*A1_2(p,q)
                    enddo
                    enddo
                    enddo

!! Unstable (as explained by Malaspinas 2015 ???)
!                    call dgemv('n',27,27,1.0,H3,27,f(1,i,j,k),1,0.0,A1_3,1)
!                    write(*,'(a,3i4,27g13.5)')'A1_2:',i,j,k,A1_3(:,:,:)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Rfneq from \citet{fen21a}, as defined in Eq. (34)
                  call dgemv('t', 9,27,inv2cs4,H2, 9,A1_2,1,0.0,f(1,i,j,k),1)
                  call dgemv('t',27,27,inv6cs6,H3,27,A1_3,1,1.0,f(1,i,j,k),1)
                  f(:,i,j,k)=weights(:)*f(:,i,j,k)

!                  do l=1,nl
!                     f(l,i,j,k)=0.0
!
!                       do q=1,3
!                       do p=1,3
!                          f(l,i,j,k)=f(l,i,j,k) + H2(p,q,l)*A1_2(p,q)*inv2cs4
!                       enddo
!                       enddo
!
!                       do r=1,3
!                       do q=1,3
!                       do p=1,3
!                          f(l,i,j,k)=f(l,i,j,k) + H3(p,q,r,l)*A1_3(p,q,r)*inv6cs6
!                       enddo
!                       enddo
!                       enddo
!
!                      f(l,i,j,k)=weights(l)*f(l,i,j,k)
!                  enddo

            enddo
!$OMP END PARALLEL DO
         enddo
      enddo
  endif

   call cpufinish(icpu)


end subroutine

end module
