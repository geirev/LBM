module m_fregularization
contains

subroutine fregularization(f, feq, rho, u, v, w)
   use mod_dimensions
   use mod_D3Q27setup
   use m_ablim
   use m_readinfile
   use m_wtime
   implicit none
   real, intent(in)      :: rho(nx,ny,nz)
   real, intent(in)      :: u(nx,ny,nz)
   real, intent(in)      :: v(nx,ny,nz)
   real, intent(in)      :: w(nx,ny,nz)
   real, intent(in)      :: feq(nl,0:nx+1,0:ny+1,0:nz+1)
   real, intent(inout)   :: f(nl,0:nx+1,0:ny+1,0:nz+1)

   logical, save         :: lfirst=.true.

   real, save            :: H2(3,3,27)      ! Second order Hermite polynomial
   real, save            :: H3(3,3,3,27)    ! Third order Hermite polynomial
   real, save            :: c(3,nl)         ! Array storage of cxs, cys, and czs

   real                  :: A0_2(3,3)
   real                  :: A1_2(3,3)
   real                  :: A1_3(3,3,3)

   real                  :: lfneq(nl)       ! Local non-equilibrium distribution

   real                  :: delta(1:3, 1:3) = reshape([1.0, 0.0, 0.0, &
                                                       0.0, 1.0, 0.0, &
                                                       0.0, 0.0, 1.0], [3, 3])

   real                  :: vel(1:3),dens

   integer :: i, j, k, l, p, q, r

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

      lfirst=.false.
   endif

! lfneq is defined in \citet{fen21a} between Eqs (32) and (33)
      do k=1,nz
         do j=1,ny
            do i=1,nx
               do l=1,nl
                  f(l,i,j,k)=f(l,i,j,k)-feq(l,i,j,k)
               enddo
            enddo
         enddo
      enddo

   if (ihrr == 1) then
! projecting non-equilibrium distribution on the Hermitian polynomials for regularization
!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(i, j, k, l, p, q, r, vel, dens, lfneq, A0_2, A1_2, A1_3 ) &
!$OMP&                           SHARED(feq, f, rho, u, v, w,  weights, H2, H3)
      do k=1,nz
         do j=1,ny
            do i=1,nx

               vel(1)=u(i,j,k)
               vel(2)=v(i,j,k)
               vel(3)=w(i,j,k)
               dens=rho(i,j,k)
               lfneq(:)=f(:,i,j,k)

! A0_2 from \citet{fen21a} (following Eq. 32)
               do q=1,3
               do p=1,3
                  A0_2(p,q)=dens*vel(p)*vel(q)
               enddo
               enddo

! Eq (11) from  Jacob 2018 is identical to the 33a from Feng (2021)
! Used for regularization and turbulence calculation
               A1_2=0.0
               do l=1,nl
                  do q=1,3
                  do p=1,3
                     A1_2(p,q) = A1_2(p,q) + H2(p,q,l)*lfneq(l)
                  enddo
                  enddo
               enddo

! A1_3(HRR) from \citet{fen21a}, as defined after Eq. (34)
                do r=1,3
                do q=1,3
                do p=1,3
                   A1_3(p,q,r)=vel(p)*A1_2(q,r) +vel(q)*A1_2(r,p) +  vel(r)*A1_2(p,q)
                enddo
                enddo
                enddo

! Rfneq from \citet{fen21a}, as defined in Eq. (34)
               do l=1,nl
                  lfneq(l)=0.0

                  do p=1,3
                  do q=1,3
                     lfneq(l)=lfneq(l) + H2(p,q,l)*A1_2(p,q)/(2.0*cs4)
                  enddo
                  enddo

                  lfneq(l)=lfneq(l)   &
                      + ( H3(1,1,2,l) + H3(2,3,3,l) ) * ( A1_3(1,1,2) + A1_3(2,3,3) )/(2.0*cs6) &
                      + ( H3(1,3,3,l) + H3(1,2,2,l) ) * ( A1_3(1,3,3) + A1_3(1,2,2) )/(2.0*cs6) &
                      + ( H3(2,2,3,l) + H3(1,1,3,l) ) * ( A1_3(2,2,3) + A1_3(1,1,3) )/(2.0*cs6) &
                      + ( H3(1,1,2,l) - H3(2,3,3,l) ) * ( A1_3(1,1,2) - A1_3(2,3,3) )/(6.0*cs6) &
                      + ( H3(1,3,3,l) - H3(1,2,2,l) ) * ( A1_3(1,3,3) - A1_3(1,2,2) )/(6.0*cs6) &
                      + ( H3(2,2,3,l) - H3(1,1,3,l) ) * ( A1_3(2,2,3) - A1_3(1,1,3) )/(6.0*cs6)

                  lfneq(l)=weights(l)*lfneq(l)
               enddo

               f(:,i,j,k) = lfneq(:)

            enddo
         enddo
      enddo
!$OMP END PARALLEL DO
   endif

   call cpufinish(icpu)


end subroutine

end module
