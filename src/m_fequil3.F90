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
   real, intent(out)     :: feq(nl,0:nx+1,0:ny+1,0:nz+1)


   logical, save         :: lfirst=.true.
   real, save            :: H2(nl,3,3)      ! Second order Hermite polynomial
   real, save            :: H3(nl,3,3,3)    ! Third order Hermite polynomial

   real, save            :: c(3,nl)         ! Array storage of cxs, cys, and czs

   real                  :: A0_2(3,3)

   real                  :: A0_3(3,3,3)

   real                  :: delta(1:3, 1:3) = reshape([1.0, 0.0, 0.0, &
                                                       0.0, 1.0, 0.0, &
                                                       0.0, 0.0, 1.0], [3, 3])

   real                  :: vel(1:3)

   integer :: i, j, k, l, p, q, r, ia

   real, parameter :: inv1cs2 = 1.0/(cs2)
   real, parameter :: inv2cs4 = 1.0/(2.0*cs4)
   real, parameter :: inv2cs6 = 1.0/(2.0*cs6)
   real, parameter :: inv6cs6 = 1.0/(6.0*cs6)
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
            H2(l,p,q)=c(p,l)*c(q,l) - cs2*delta(p,q)
         enddo
         enddo
      enddo

! Hermitian polynomials of 3nd order H3 with
! [c_\alpha \delta]_{ijk} = c_{\alpha,i} \delta_{jk} + c_{\alpha,j} \delta_{ik} +c_{\alpha,k} \delta_{ij}
      do l=1,nl
         do r=1,3
         do q=1,3
         do p=1,3
            H3(l,p,q,r)=c(p,l)*c(q,l)*c(r,l) - cs2*(c(p,l)*delta(q,r) + c(q,l)*delta(p,r) +  c(r,l)*delta(p,q))
         enddo
         enddo
         enddo
      enddo

      lfirst=.false.
   endif


! Loop over grid
!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(i, ia, j, k, l, p, q, r, vel, A0_2, A0_3 )      &
!$OMP&                          SHARED(feq, rho, u, v, w, weights, c, H2, ibgk, H3)
   do k=1,nz
      do j=1,ny
         do i=1,nx
            ia=max(i,1)

            vel(1)=u(ia,j,k)
            vel(2)=v(ia,j,k)
            vel(3)=w(ia,j,k)

! A0_2 and A0_3 from \citet{fen21a} (following Eq. 32)
            do q=1,3
            do p=1,3
               A0_2(p,q)=rho(ia,j,k)*vel(p)*vel(q)
            enddo
            enddo

            do r=1,3
            do q=1,3
            do p=1,3
               A0_3(p,q,r)=rho(ia,j,k)*vel(p)*vel(q)*vel(r)
            enddo
            enddo
            enddo


! Equilibrium distribution \citet{fen21a} Eq. (32) or jac18a eq (27)
            feq( 1,i,j,k) = rho(ia,j,k) * (cs2                              ) / cs2
            feq( 2,i,j,k) = rho(ia,j,k) * (cs2  + vel(1)                    ) / cs2
            feq( 3,i,j,k) = rho(ia,j,k) * (cs2  - vel(1)                    ) / cs2
            feq( 4,i,j,k) = rho(ia,j,k) * (cs2            + vel(2)          ) / cs2
            feq( 5,i,j,k) = rho(ia,j,k) * (cs2            - vel(2)          ) / cs2
            feq( 6,i,j,k) = rho(ia,j,k) * (cs2                      - vel(3)) / cs2
            feq( 7,i,j,k) = rho(ia,j,k) * (cs2                      + vel(3)) / cs2
            feq( 8,i,j,k) = rho(ia,j,k) * (cs2  + vel(1)  + vel(2)          ) / cs2
            feq( 9,i,j,k) = rho(ia,j,k) * (cs2  - vel(1)  - vel(2)          ) / cs2
            feq(10,i,j,k) = rho(ia,j,k) * (cs2  + vel(1)  - vel(2)          ) / cs2
            feq(11,i,j,k) = rho(ia,j,k) * (cs2  - vel(1)  + vel(2)          ) / cs2
            feq(12,i,j,k) = rho(ia,j,k) * (cs2  - vel(1)            - vel(3)) / cs2
            feq(13,i,j,k) = rho(ia,j,k) * (cs2  + vel(1)            + vel(3)) / cs2
            feq(14,i,j,k) = rho(ia,j,k) * (cs2            + vel(2)  + vel(3)) / cs2
            feq(15,i,j,k) = rho(ia,j,k) * (cs2            - vel(2)  - vel(3)) / cs2
            feq(16,i,j,k) = rho(ia,j,k) * (cs2  - vel(1)            + vel(3)) / cs2
            feq(17,i,j,k) = rho(ia,j,k) * (cs2  + vel(1)            - vel(3)) / cs2
            feq(18,i,j,k) = rho(ia,j,k) * (cs2            - vel(2)  + vel(3)) / cs2
            feq(19,i,j,k) = rho(ia,j,k) * (cs2            + vel(2)  - vel(3)) / cs2
            feq(20,i,j,k) = rho(ia,j,k) * (cs2  - vel(1)  + vel(2)  + vel(3)) / cs2
            feq(21,i,j,k) = rho(ia,j,k) * (cs2  + vel(1)  - vel(2)  - vel(3)) / cs2
            feq(22,i,j,k) = rho(ia,j,k) * (cs2  - vel(1)  - vel(2)  - vel(3)) / cs2
            feq(23,i,j,k) = rho(ia,j,k) * (cs2  + vel(1)  + vel(2)  + vel(3)) / cs2
            feq(24,i,j,k) = rho(ia,j,k) * (cs2  + vel(1)  + vel(2)  - vel(3)) / cs2
            feq(25,i,j,k) = rho(ia,j,k) * (cs2  - vel(1)  - vel(2)  + vel(3)) / cs2
            feq(26,i,j,k) = rho(ia,j,k) * (cs2  - vel(1)  + vel(2)  - vel(3)) / cs2
            feq(27,i,j,k) = rho(ia,j,k) * (cs2  + vel(1)  - vel(2)  + vel(3)) / cs2


!            do p=1,3
!            do q=1,3
!               do l=1,nl
!                  feq(l,i,j,k)=feq(l,i,j,k) + H2(l,p,q)*A0_2(p,q)/(2.0*cs4)
!               enddo
!            enddo
!            enddo
            call dgemv('n', 27,9,inv2cs4,H2, 27,A0_2,1,1.0,feq(1,i,j,k),1)

! the above identically recovers the BGK equilibrium, now we add third order contributions
            if (ibgk == 3) then
!               do l=1,nl
!                     do p=1,3
!                     do q=1,3
!                     do r=1,3
!                        feq(l,i,j,k)=feq(l,i,j,k) + H3(l,p,q,r)*A0_3(p,q,r)*inv6cs6
!                     enddo
!                     enddo
!                     enddo
!               enddo
               call dgemv('n',27,27,inv6cs6,H3,27,A0_3,1,1.0,feq(1,i,j,k),1)
            endif

! scaling by the weights
            do l=1,nl
               feq(l,i,j,k)= weights(l)*feq(l,i,j,k)
            enddo

         enddo
      enddo
   enddo
!$OMP END PARALLEL DO

   call cpufinish(icpu)


end subroutine

end module
