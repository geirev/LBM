module m_HRRequil
! NOTE:  f^coll = f - (1/tau) * (f - f^eq)
!               = f^eq + (1-1/tau) * f^neq       # f^neq= f-f^eq
!               ~ f^eq + (1-1/tau) * R(f^neq)
contains

subroutine HRRequil(feq, f, rho, u, v, w, tau)
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
   real, intent(inout)   :: f(0:nx+1,0:ny+1,0:nz+1,nl)
   real, intent(out)     :: tau(nx,ny,nz)
   !                        ihrr      ! (1) full third order HRR scheme returing feq and Rfneq
                                      ! (2) second order standard BGK returing feq and fneq
                                      ! (3) third order BGK that returns fneq instead of Rfneq,
                                      ! (0) Initialization by feq

   real, parameter       :: cs2=1/3.0
   real, parameter       :: cs4=1/9.0
   real, parameter       :: cs6=1/27.0

   logical, save         :: lfirst=.true.

   real, save            :: H2(3,3,27)      ! Second order Hermite polynomial
   real, save            :: H3(3,3,3,27)    ! Third order Hermite polynomial
   real, save            :: c(3,nl)         ! Array storage of cxs, cys, and czs

   real, save            :: uscale          ! velocity scale lattice to physical
   real, save            :: dx              ! length scale lattice to physical
   real, save            :: dt              ! time scale lattice to physical
   real, save            :: const           ! c in Vreman 2004 Eq (5)

   real                  :: A0_2(3,3)
   real                  :: A0_3(3,3,3)
   real                  :: A1_2(3,3)
   real                  :: A1_2FD(3,3)
   real                  :: A1_2HRR(3,3)
   real                  :: A1_3HRR(3,3,3)

   real                  :: Rfneq(nl)       ! Local axroximate nonequilibrium distribution
   real                  :: lfneq(nl)       ! Local equilibrium distribution
   real                  :: lfeq(nl)        ! Local equilibrium distribution
   real                  :: lf(nl)          ! Local predicted distribution


   real                  :: delta(1:3, 1:3) = reshape([1, 0, 0, &
                                                       0, 1, 0, &
                                                       0, 0, 1], [3, 3])

   real dudx, dudy, dudz, dvdx, dvdy, dvdz, dwdx, dwdy, dwdz
   real dxfac, dyfac, dzfac
   real                  :: vel(1:3),dens

   integer :: i, j, k, l, m, p, q, r, ia, ib, ja, jb, ka , kb


   real, parameter :: sigma=1.00
   real, parameter :: smagorinsky=0.10      !0.15      !0.18   Smagorinsky 0.065 from abk18a

   real eddyvisc  ! nu in Vreman 2004 Eq (5)
   real Bbeta     ! B_beta in Vreman 2004 Eq (5)
   real alpha(3,3)
   real beta(3,3)
   real alphamag






   integer, parameter :: icpu=4
   call cpustart()


   if (lfirst) then
      const=2.5*smagorinsky**2
      uscale=p2l%vel
      dx=p2l%length
      dt=p2l%time


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
!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(i, ia, ib, j, ja, jb, k, ka, kb, l, p, q, r,            &
!$OMP&                                  dxfac, dyfac, dzfac, vel, dens, lf, lfeq, lfneq, Rfneq, &
!$OMP&                                  A0_2, A0_3, A1_2, A1_2FD, A1_2HRR, A1_3HRR,             &
!$OMP&                                  dudx, dudy, dudz, dvdx, dvdy, dvdz, dwdx, dwdy, dwdz,   &
!$OMP&                                  eddyvisc, Bbeta, alpha, beta, alphamag) &
!$OMP&                           SHARED(feq, f, rho, u, v, w, tau, weights, c, H2, H3, tauin,   &
!$OMP&                                  kinevisc, uscale, dx, dt, const, ihrr)
   do k=1,nz
      call ablim(k,nz,dzfac,ka,kb)
      do j=1,ny
         call ablim(j,ny,dyfac,ja,jb)
         do i=1,nx
            call ablim(i,nx,dxfac,ia,ib)

            vel(1)=u(i,j,k)
            vel(2)=v(i,j,k)
            vel(3)=w(i,j,k)
            dens=rho(i,j,k)
            lf(:)=f(i,j,k,:)

! A0_2 and A0_3 from \citet{fen21a} (following Eq. 32)
            do q=1,3
            do p=1,3
               A0_2(p,q)=dens*vel(p)*vel(q)
            enddo
            enddo

            do r=1,3
            do q=1,3
            do p=1,3
               A0_3(p,q,r)=dens*vel(p)*vel(q)*vel(r)
            enddo
            enddo
            enddo

! Equilibrium distribution \citet{fen21a} Eq. (32) or jac18a eq (27)
            do l=1,nl
               lfeq(l)=dens

               lfeq(l)=lfeq(l) + dens*( c(1,l)*vel(1) + c(2,l)*vel(2) + c(3,l)*vel(3) )/cs2

               do p=1,3
               do q=1,3
                  lfeq(l)=lfeq(l) + H2(p,q,l)*A0_2(p,q)/(2.0*cs4)
               enddo
               enddo
               ! the above identically recovers the BGK equilibrium, below we add third order contributions
               ! Add third order correction
               if (ihrr /= 2) then
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

! lfneq is defined in \citet{fen21a} between Eqs (32) and (33)
            lfneq(:)=lf(:)-lfeq(:)

            if (ihrr == 1) then
! A1_2 from \citet{fen21a} Eq. (33a)
!              A1_2=0.0
!              do l=1,nl
!                 do q=1,3
!                 do p=1,3
!                    A1_2(p,q) = A1_2(p,q) + c(p,l)*c(q,l)*lfneq(l)
!                 enddo
!                 enddo
!              enddo
! Eq (11) from  Jacob 2018 is idetical to the 33a from Feng (2021)
               A1_2=0.0
               do l=1,nl
                  do q=1,3
                  do p=1,3
                     A1_2(p,q) = A1_2(p,q) + H2(p,q,l)*lfneq(l)
                  enddo
                  enddo
               enddo

               if (sigma /= 1.0) then
! A1_2FD from \citet{fen21a} Eq. (33b)
                  dudx=(u(ib,j,k)-u(ia,j,k))/(dxfac*dx)
                  dudy=(u(i,jb,k)-u(i,ja,k))/(dyfac*dx)
                  dudz=(u(i,j,kb)-u(i,j,ka))/(dzfac*dx)

                  dvdx=(v(ib,j,k)-v(ia,j,k))/(dxfac*dx)
                  dvdy=(v(i,jb,k)-v(i,ja,k))/(dyfac*dx)
                  dvdz=(v(i,j,kb)-v(i,j,ka))/(dzfac*dx)

                  dwdx=(w(ib,j,k)-w(ia,j,k))/(dxfac*dx)
                  dwdy=(w(i,jb,k)-w(i,ja,k))/(dyfac*dx)
                  dwdz=(w(i,j,kb)-w(i,j,ka))/(dzfac*dx)

                  A1_2FD(1,1) = dudx + dudx - (2.0/3.0)* (dudx + dvdy + dwdz)
                  A1_2FD(2,1) = dudy + dvdx
                  A1_2FD(3,1) = dudz + dwdx
                  A1_2FD(1,2) = A1_2FD(2,1)
                  A1_2FD(2,2) = dvdy + dvdy - (2.0/3.0)* (dudx + dvdy + dwdz)
                  A1_2FD(3,2) = dvdz + dwdy
                  A1_2FD(1,3) = A1_2FD(3,1)
                  A1_2FD(2,3) = A1_2FD(3,2)
                  A1_2FD(3,3) = dwdz + dwdz - (2.0/3.0)* (dudx + dvdy + dwdz)

                  A1_2FD = -2.0 * uscale * tauin * dens * cs2 * A1_2FD

! A1_2HRR from \citet{fen21a}, as defined after Eq. (34)
                  A1_2HRR = sigma*A1_2 + (1.0-sigma)*A1_2FD
               else
                  A1_2HRR=A1_2
               endif

! A1_3HRR from \citet{fen21a}, as defined after Eq. (34)
               do r=1,3
               do q=1,3
               do p=1,3
                  A1_3HRR(p,q,r)=vel(p)*A1_2HRR(q,r) +vel(q)*A1_2HRR(r,p) +  vel(r)*A1_2HRR(p,q)
               enddo
               enddo
               enddo

! Rfneq from \citet{fen21a}, as defined in Eq. (34)
               do l=1,nl
                  Rfneq(l)=0.0

                  do p=1,3
                  do q=1,3
                     Rfneq(l)=Rfneq(l) + H2(p,q,l)*A1_2HRR(p,q)/(2.0*cs4)
                  enddo
                  enddo

                  Rfneq(l)=Rfneq(l)   &
                      + ( H3(1,1,2,l) + H3(2,3,3,l) ) * ( A1_3HRR(1,1,2) + A1_3HRR(2,3,3) )/(2.0*cs6) &
                      + ( H3(1,3,3,l) + H3(1,2,2,l) ) * ( A1_3HRR(1,3,3) + A1_3HRR(1,2,2) )/(2.0*cs6) &
                      + ( H3(2,2,3,l) + H3(1,1,3,l) ) * ( A1_3HRR(2,2,3) + A1_3HRR(1,1,3) )/(2.0*cs6) &
                      + ( H3(1,1,2,l) - H3(2,3,3,l) ) * ( A1_3HRR(1,1,2) - A1_3HRR(2,3,3) )/(6.0*cs6) &
                      + ( H3(1,3,3,l) - H3(1,2,2,l) ) * ( A1_3HRR(1,3,3) - A1_3HRR(1,2,2) )/(6.0*cs6) &
                      + ( H3(2,2,3,l) - H3(1,1,3,l) ) * ( A1_3HRR(2,2,3) - A1_3HRR(1,1,3) )/(6.0*cs6)

                  Rfneq(l)=weights(l)*Rfneq(l)
               enddo
            else ! (ihrr /= 1)
! Third order BGK without expansion for fneq
               RFneq(:)=lfneq(:)
            endif

!            if ((i==nx/2).and.(j==ny/2).and.(k==nz/2)) then
!               print '(a,i2)','A1_2 ihrr=',ihrr
!               print '(3e13.5)',A1_2
!               print '(a)','A1_2FD'
!               print '(3e13.5)',A1_2FD
!               print '(a)','A1_2 - A1_2FD'
!               print '(3e13.5)',A1_2-A1_2FD
!               print *
!               print '(a)','lfeq:'
!               print '(10e13.5)',lfeq(:)
!               print '(a)','lfneq:'
!               print '(10e13.5)',lfneq(:)
!               print '(a)','Rfneq'
!               print '(10e13.5)',Rfneq(:)
!               print '(a)','Rfneq-lfneq'
!               print '(10e13.5)',Rfneq(:)-lfneq(:)
!               print *
!            endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  Vreman (2004) subgridscale turbulence model
            eddyvisc=0.0
            if (ihrr == 1) then
! S_ij
               alpha=A1_2

! alphamag
               alphamag=0.00001
               do q=1,3
               do p=1,3
                  alphamag=alphamag+alpha(p,q)*alpha(p,q)
               enddo
               enddo
               !print *,'HRR alphamag:',alphamag

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
               !print *,'HRR Bbeta:',Bbeta,Bbeta/alphamag

! Vreman 2004 Eq (5)
               eddyvisc=const*sqrt(Bbeta/alphamag)
               !print *,'HRR eddyvisc:',eddyvisc

            endif

            tau(i,j,k) = 3.0*(kinevisc + eddyvisc) + 0.5
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

            f(i,j,k,:) = Rfneq(:)
            feq(i,j,k,:)= lfeq(:)
           !if ((i==nx/2).and.(j==ny/2).and.(k==nz/2)) then
           !   print '(a,3e13.5)','tau=',tau(i,j,k),kinevisc,eddyvisc
           !endif

         enddo
      enddo
   enddo
!$OMP END PARALLEL DO

   call cpufinish(icpu)


end subroutine

end module
