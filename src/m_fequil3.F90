module m_fequil3
contains

subroutine fequil3(feq, rho, u, v, w)
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
   real, intent(out)     :: feq(0:nx+1,0:ny+1,0:nz+1,nl)

   logical, save         :: lfirst=.true.
   real, save            :: H2(3,3,27)      ! Second order Hermite polynomial
   real, save            :: H3(3,3,3,27)    ! Third order Hermite polynomial
   real, save            :: c(3,nl)         ! Array storage of cxs, cys, and czs

   real                  :: lfeq(nl)

   real                  :: A0_2(3,3)
   real                  :: A0_3(3,3,3)
   real                  :: delta(1:3, 1:3) = reshape([1.0, 0.0, 0.0, &
                                                       0.0, 1.0, 0.0, &
                                                       0.0, 0.0, 1.0], [3, 3])

   real                  :: vel(1:3)

   integer :: i, j, k, l, p, q, r

   integer, parameter :: icpu=11
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


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Loop over grid
!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(i, j, k, l, p, q, r, vel,A0_2, A0_3, lfeq    ) &
!$OMP&                          SHARED(feq, rho, u, v, w, weights, c, H2, H3, ibgk)
   do k=1,nz
      do j=1,ny
         do i=1,nx

            vel(1)=u(i,j,k)
            vel(2)=v(i,j,k)
            vel(3)=w(i,j,k)

! A0_2 and A0_3 from \citet{fen21a} (following Eq. 32)
            do q=1,3
            do p=1,3
               A0_2(p,q)=rho(i,j,k)*vel(p)*vel(q)
            enddo
            enddo

            do r=1,3
            do q=1,3
            do p=1,3
               A0_3(p,q,r)=rho(i,j,k)*vel(p)*vel(q)*vel(r)
            enddo
            enddo
            enddo


! Equilibrium distribution \citet{fen21a} Eq. (32) or jac18a eq (27)
            do l=1,nl
               lfeq(l)=rho(i,j,k)

               lfeq(l)=lfeq(l) + rho(i,j,k)*( c(1,l)*vel(1) + c(2,l)*vel(2) + c(3,l)*vel(3) )/cs2

               do p=1,3
               do q=1,3
                  lfeq(l)=lfeq(l) + H2(p,q,l)*A0_2(p,q)/(2.0*cs4)
               enddo
               enddo

! the above identically recovers the BGK equilibrium, now we add third order contributions
               if (ibgk == 3) then
                  lfeq(l)=lfeq(l)   &
                      + ( H3(1,1,2,l) + H3(2,3,3,l) ) * ( A0_3(1,1,2) + A0_3(2,3,3) )/(2.0*cs6) &
                      + ( H3(1,3,3,l) + H3(1,2,2,l) ) * ( A0_3(1,3,3) + A0_3(1,2,2) )/(2.0*cs6) &
                      + ( H3(2,2,3,l) + H3(1,1,3,l) ) * ( A0_3(2,2,3) + A0_3(1,1,3) )/(2.0*cs6) &
                      + ( H3(1,1,2,l) - H3(2,3,3,l) ) * ( A0_3(1,1,2) - A0_3(2,3,3) )/(6.0*cs6) &
                      + ( H3(1,3,3,l) - H3(1,2,2,l) ) * ( A0_3(1,3,3) - A0_3(1,2,2) )/(6.0*cs6) &
                      + ( H3(2,2,3,l) - H3(1,1,3,l) ) * ( A0_3(2,2,3) - A0_3(1,1,3) )/(6.0*cs6)
               endif

               lfeq(l)= weights(l)*lfeq(l)

            enddo
            feq(i,j,k,:)=lfeq(:)
         enddo
      enddo
   enddo
!$OMP END PARALLEL DO

   call cpufinish(icpu)

end subroutine
end module
