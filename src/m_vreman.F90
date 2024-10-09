module m_vreman
!  Vreman (2004) subgridscale turbulence model
contains

subroutine vreman(f,tau)
   use mod_dimensions
   use mod_D3Q27setup
   use m_readinfile
   use m_wtime
   implicit none
   real, intent(in)      :: f(0:nx+1,0:ny+1,0:nz+1,nl) ! Nonequilibrium f as input
   real, intent(out)     :: tau(nx,ny,nz)              ! Tau including subgrid scale mixing

   logical, save         :: lfirst=.true.

   real, save            :: H2(3,3,27)      ! Second order Hermite polynomial
   real, save            :: c(3,nl)         ! Array storage of cxs, cys, and czs
   real, save            :: dx              ! length scale lattice to physical
   real, save            :: const           ! c in Vreman 2004 Eq (5)

   real                  :: lfneq(nl)       ! Local equilibrium distribution

   real                  :: delta(1:3, 1:3) = reshape([1.0, 0.0, 0.0, &
                                                       0.0, 1.0, 0.0, &
                                                       0.0, 0.0, 1.0], [3, 3])


   integer :: i, j, k, l, m, p, q

   real eddyvisc  ! nu in Vreman 2004 Eq (5)
   real Bbeta     ! B_beta in Vreman 2004 Eq (5)
   real alpha(3,3)
   real beta(3,3)
   real alphamag

   integer, parameter :: icpu=13
   call cpustart()

   if (ivreman /= 1) then
      tau = 3.0*kinevisc + 0.5
      return
   endif

   if (lfirst) then
      const=2.5*smagorinsky**2
      dx=p2l%length

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

      lfirst=.false.
   endif


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Loop over grid
!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(i, j, k, l, p, q, lfneq, eddyvisc, Bbeta, alpha, beta, alphamag) &
!$OMP&                           SHARED(f, tau, weights, c, H2, tauin, kinevisc, dx, const, ivreman)
   do k=1,nz
      do j=1,ny
         do i=1,nx

            lfneq(:)=f(i,j,k,:)

! Eq (11) from Jacob 2018 is identical to the 33a from Feng (2021)
            alpha=0.0
            do l=1,nl
               do q=1,3
               do p=1,3
                  alpha(p,q) = alpha(p,q) + H2(p,q,l)*lfneq(l)
               enddo
               enddo
            enddo

! alphamag
            alphamag=0.00001
            do q=1,3
            do p=1,3
               alphamag=alphamag+alpha(p,q)*alpha(p,q)
            enddo
            enddo

! beta = del^2 * alpha' * alpha
            beta=0.00001
            do q=1,3
            do p=1,3
               do m=1,3
                  beta(p,q)=beta(p,q)+alpha(m,p)*alpha(m,q)
               enddo
            enddo
            enddo
            beta=dx**2*beta

! Vreman 2004 Eq (8)
            Bbeta=beta(1,1)*beta(2,2) - beta(1,2)**2  &
                 +beta(1,1)*beta(3,3) - beta(1,3)**2  &
                 +beta(2,2)*beta(3,3) - beta(2,3)**2

! Vreman 2004 Eq (5)
            eddyvisc=const*sqrt(Bbeta/alphamag)

            tau(i,j,k) = 3.0*(kinevisc + eddyvisc) + 0.5

         enddo
      enddo
   enddo
!$OMP END PARALLEL DO

   call cpufinish(icpu)

end subroutine
end module

