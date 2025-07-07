module m_hermite_polynomials
contains
subroutine hermite_polynomials()
   use mod_dimensions
   use mod_D3Q27setup
   implicit none
   integer :: l, p, q, r
   real :: c(3,nl), delta(3,3)

   allocate(weights(1:nl))
   allocate(H2(3,3,1:nl))
   allocate(H3(3,3,3,1:nl))

! Build lattice velocity matrix
   c(1,:) = real(cxs(:))
   c(2,:) = real(cys(:))
   c(3,:) = real(czs(:))

! Kronecker delta
   delta = 0.0
   delta(1,1) = 1.0
   delta(2,2) = 1.0
   delta(3,3) = 1.0
! Set weights
   weights = [8.0/27.0, &
              2.0/27.0, 2.0/27.0, 2.0/27.0, 2.0/27.0, 2.0/27.0, 2.0/27.0, &
              1.0/54.0, 1.0/54.0, 1.0/54.0, 1.0/54.0, 1.0/54.0, 1.0/54.0, &
              1.0/54.0, 1.0/54.0, 1.0/54.0, 1.0/54.0, 1.0/54.0, 1.0/54.0, &
              1.0/216.0, 1.0/216.0, 1.0/216.0, 1.0/216.0, 1.0/216.0, 1.0/216.0, 1.0/216.0, 1.0/216.0]

! Compute second-order Hermite polynomials
   do l = 1, nl
      do p = 1, 3
         do q = 1, 3
            H2(p,q,l) = c(p,l)*c(q,l) - cs2 * delta(p,q)
         enddo
      enddo
   enddo

   ! Compute third-order Hermite polynomials
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


end subroutine
end module
