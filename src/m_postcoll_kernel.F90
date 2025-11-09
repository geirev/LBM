! Equilibrium distribution \citet{fen21a} Eq. (32) or jac18a eq (27)
module m_postcoll_kernel
contains
#ifdef _CUDA
   attributes(global)&
#endif
   subroutine postcoll_kernel(f, tau, rho, u, v, w, inv1cs2, inv2cs4, inv6cs6, eps, kinevisc, const, ibgk, ihrr, ivreman)
#ifdef _CUDA
   use cudafor
#endif
   use mod_dimensions, only : nx,ny,nz
   use mod_D3Q27setup, only : nl,cxs,cys,czs,cs2, weights, H2, H3
   implicit none
   real, intent(inout)  :: f(nl,0:nx+1,0:ny+1,0:nz+1)
   real, intent(inout)  :: tau(0:nx+1,0:ny+1,0:nz+1)
   real, intent(in)      :: rho(0:nx+1,0:ny+1,0:nz+1)
   real, intent(in)      ::   u(0:nx+1,0:ny+1,0:nz+1)
   real, intent(in)      ::   v(0:nx+1,0:ny+1,0:nz+1)
   real, intent(in)      ::   w(0:nx+1,0:ny+1,0:nz+1)
   real, value          :: inv1cs2
   real, value          :: inv2cs4
   real, value          :: inv6cs6
   integer, value       :: ibgk
   integer, value       :: ihrr
   integer, value       :: ivreman
   real, value, intent(in) :: eps
   real, value, intent(in) :: kinevisc
   real, value, intent(in) :: const


! fequil and regularization variables
   real :: A2(3,3)
   real :: A3(3,3,3)
   real :: vel(3)
   real :: ratio
   real :: vratio
   real :: cu
   real :: tmpeq(nl)
   real :: tmpneq(nl)
   real :: dens

! vreman variables
   real   :: alpha(3,3)
   real   :: alphamag
   real   :: beta(3,3)
   real   :: Bbeta
   real   :: eddyvisc
   real   :: tautmp

   real   :: fac

   integer :: i, j, k, l, m, p, q, r
#ifdef _CUDA
   i = threadIdx%x + (blockIdx%x - 1) * blockDim%x
   j = threadIdx%y + (blockIdx%y - 1) * blockDim%y
   k = threadIdx%z + (blockIdx%z - 1) * blockDim%z
   if (i > nx .or. j > ny .or. k > nz) return
   ratio = inv6cs6 / inv2cs4
#else
   ratio = inv6cs6 / inv2cs4
!$OMP PARALLEL DO collapse(3) DEFAULT(NONE) &
!$OMP  PRIVATE(i,j,k,l,p,q,r,vel,dens,cu,tmpeq,tmpneq,A2,A3,vratio,alpha,alphamag,beta,Bbeta,eddyvisc,tautmp,fac)&
!$OMP  SHARED(f,rho,u,v,w,H2,H3,cxs,cys,czs,weights,inv1cs2,inv2cs4,inv6cs6,ratio,ibgk,ihrr,ivreman,eps,kinevisc,const,tau)
   do k=1,nz
   do j=1,ny
   do i=1,nx
#endif

! Copy u,v,w to vel(1:3)
      vel(1)=u(i,j,k)
      vel(2)=v(i,j,k)
      vel(3)=w(i,j,k)
      dens=rho(i,j,k)
! A2
      do q=1,3
      do p=1,3
         A2(p,q)=dens*vel(p)*vel(q)*inv2cs4
      enddo
      enddo

     if (ibgk == 3) then
! compute A3 using A2 to save arithmetic
         do r=1,3
         vratio = vel(r) * ratio        ! combine vel(r) with ratio to 1 value
         do q=1,3
         do p=1,3
            A3(p,q,r) = A2(p,q) * vratio
         enddo
         enddo
         enddo
      endif


! 1. order
      do l=1,nl
         cu = real(cxs(l))*vel(1) + real(cys(l))*vel(2) + real(czs(l))*vel(3)
         tmpeq(l) = dens * (1.0 + cu*inv1cs2)

! 2nd order equilibrium distribution \citet{fen21a} Eq. (32) or jac18a eq (27)
         do q=1,3
         do p=1,3
            tmpeq(l)=tmpeq(l) + H2(p,q,l)*A2(p,q)
         enddo
         enddo

         if (ibgk == 3 .and. l > 1) then
            do r=1,3
            do q=1,3
            do p=1,3
               tmpeq(l)=tmpeq(l) + H3(p,q,r,l)*A3(p,q,r)
            enddo
            enddo
            enddo
         endif
         tmpeq(l)=weights(l)*tmpeq(l)
      enddo


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Regulariation starts here

      if (ihrr==1) then
! computing A2
         A2(:,:)=0.0
         do l=1,nl
            tmpneq(l)=f(l,i,j,k)-tmpeq(l)
            do q=1,3
            do p=1,3
               A2(p,q) = A2(p,q) + h2(p,q,l)*tmpneq(l)
            enddo
            enddo
         enddo

! computing A3 using Maspalinas (2015, appedix B)
         do r=1,3
         do q=1,3
         do p=1,3
            A3(p,q,r)=vel(p)*A2(q,r) + vel(q)*A2(r,p) +  vel(r)*A2(p,q)
         enddo
         enddo
         enddo

! scale A2 and A3
         do q=1,3
         do p=1,3
            A2(p,q) = A2(p,q)*inv2cs4
         enddo
         enddo

         do r=1,3
         do q=1,3
         do p=1,3
            A3(p,q,r) = A3(p,q,r)*inv6cs6
         enddo
         enddo
         enddo

! computing fneq
         do l=1,nl
            tmpneq(l) = 0.0
            do q=1,3
            do p=1,3
               tmpneq(l) = tmpneq(l) + h2(p,q,l)*A2(p,q)
            end do
            end do

            do r=1,3
            do q=1,3
            do p=1,3
               tmpneq(l) = tmpneq(l) + h3(p,q,r,l)*A3(p,q,r)
            end do
            end do
            end do

            tmpneq(l) = weights(l) * tmpneq(l)
            !f(l,i,j,k) = tmpneq(l)
         end do
      else
         do l=1,nl
            tmpneq(l) =  f(l,i,j,k) - tmpeq(l)
            !f(l,i,j,k) = tmpneq(l)
            !f(l,i,j,k) = f(l,i,j,k) -  tmpeq(l)
         enddo
      endif

! vreman stuff
! Eq (11) from Jacob 2018 is identical to the 33a from Feng (2021)
      if (ivreman == 1) then
         do q=1,3
         do p=1,3
            alpha(p,q) =  0.0
         enddo
         enddo

         do l=1,nl
            do q=1,3
            do p=1,3
               alpha(p,q) = alpha(p,q) + H2(p,q,l)*tmpneq(l)
            enddo
            enddo
         enddo

         alphamag=0.0
         do q=1,3
         do p=1,3
            alphamag=alphamag+alpha(p,q)*alpha(p,q)
         enddo
         enddo
         alphamag = max(alphamag, tiny(alphamag))

   !! beta = del^2 * alpha' * alpha
         do q=1,3
         do p=1,3
            beta(p,q)=eps
            do m=1,3
               beta(p,q)=beta(p,q)+alpha(m,p)*alpha(m,q)
            enddo
         enddo
         enddo


         Bbeta=beta(1,1)*beta(2,2) - beta(1,2)**2  &
              +beta(1,1)*beta(3,3) - beta(1,3)**2  &
              +beta(2,2)*beta(3,3) - beta(2,3)**2


         tautmp=Bbeta/alphamag
         if (tautmp > eps) then
            eddyvisc=const*sqrt(tautmp)
         else
            eddyvisc=0.0
         endif
         tau(i,j,k) = 3.0*(kinevisc + eddyvisc) + 0.5
      endif


      fac = 1.0 - 1.0/tau(i,j,k)
      do l=1,nl
         f(l,i,j,k) = tmpeq(l) + fac * tmpneq(l)
      enddo

#ifndef _CUDA
    enddo
    enddo
    enddo
!$OMP END PARALLEL DO
#endif

end subroutine
end module
