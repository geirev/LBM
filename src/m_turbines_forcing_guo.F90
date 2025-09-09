module m_turbines_forcing_guo
contains

subroutine turbines_forcing_guo(df,du,dv,dw,rho,u,v,w,ip,iradius,nturbines,n)
! (8) Guo 2002
!     function [U,S]=SchemeVII(A,dt,tau,f,Rho,U)
!        U= U +  dt*A./2
!        S= Rho.*W’.*(Cdot(A)- sum(A.*(U+du),1))./ Cs2 *(1-1/(2*tau))
!        S=S+Rho.*W’.*(Cdot(A) .* Cdot(U+du) )./ Cs2^2*(1-1/(2*tau))
!     end
   use mod_dimensions
   use m_turbines_init, only : ieps
   use m_wtime
   use mod_D3Q27setup
#ifdef _CUDA
   use cudafor
#endif
   implicit none

! Computing the S_i term returned in df
   integer, intent(in) :: nturbines
   real, intent(inout) :: df(nl,-ieps:ieps,ny,nz,nturbines)
   real, intent(in)    :: du(-ieps:ieps,ny,nz)
   real, intent(in)    :: dv(-ieps:ieps,ny,nz)
   real, intent(in)    :: dw(-ieps:ieps,ny,nz)

   real, intent(in)    :: rho(nx,ny,nz)
   real, intent(inout) :: u(nx,ny,nz)
   real, intent(inout) :: v(nx,ny,nz)
   real, intent(inout) :: w(nx,ny,nz)


   integer, intent(in) :: ip, iradius, n
   integer :: l
!   real cx,cy,cz,cdotuu,cminx,cminy,cminz
   real cdota,udota,cdotu

#ifdef _CUDA
   attributes(device)  :: df
   attributes(device)  :: du
   attributes(device)  :: dv
   attributes(device)  :: dw
   attributes(device)  :: rho
   attributes(device)  :: u
   attributes(device)  :: v
   attributes(device)  :: w
#endif

   integer i,j,k



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! updating the forced equilibrium velocities ueq=u + F/(2rho)

#ifdef _CUDA
!$cuf kernel do(2) <<<*,*>>>
#else
!!  !$OMP PARALLEL DO PRIVATE(i,j,k,cdotuu,cminx,cminy,cminz,cx,cy,cz) &
!$OMP PARALLEL DO PRIVATE(i,j,k,l,cdota,cdotu,udota) &
!$OMP             SHARED(n, df, rho, weights, u, v, w, ip, du, dv, dw ,cxs, cys, czs)
#endif
   do k=1,nz
   do j=1,ny
   do i=-ieps,ieps
! updating the forced equilibrium velocities ueq=u + F/(2rho)
      u(ip+i,j,k)=u(ip+i,j,k)+0.5*du(i,j,k)
      v(ip+i,j,k)=v(ip+i,j,k)+0.5*dv(i,j,k)
      w(ip+i,j,k)=w(ip+i,j,k)+0.5*dw(i,j,k)

! Computing the S_i term returned in df
      do l=1,nl
         cdota =real(cxs(l))*du(i,j,k) + real(cys(l))*dv(i,j,k) + real(czs(l))*dw(i,j,k)
         cdotu =real(cxs(l))*u(ip+i,j,k) + real(cys(l))*v(ip+i,j,k) + real(czs(l))*w(ip+i,j,k)
         udota =u(ip+i,j,k)*du(i,j,k) + v(ip+i,j,k)*dv(i,j,k) + w(ip+i,j,k)*dw(i,j,k)
         df(l,i,j,k,n)= rho(ip+i,j,k) * weights(l) * ( (cdota - udota)/cs2 + (cdota*cdotu )/cs4 )
      enddo

!      do l=1,nl
!         cdotuu=(real(cxs(l))*u(ip+i,j,k) + real(cys(l))*v(ip+i,j,k) + real(czs(l))*w(ip+i,j,k))/cs4
!         cminx=(real(cxs(l))-u(ip+i,j,k))/cs2
!         cminy=(real(cys(l))-v(ip+i,j,k))/cs2
!         cminz=(real(czs(l))-w(ip+i,j,k))/cs2
!
!         cx=(cminx + cdotuu*real(cxs(l)))*du(i,j,k)
!         cy=(cminy + cdotuu*real(cys(l)))*dv(i,j,k)
!         cz=(cminz + cdotuu*real(czs(l)))*dw(i,j,k)
!
!         df(l,i,j,k,n)= rho(ip+i,j,k) * weights(l) * (cx + cy + cz)
!      enddo
   enddo
   enddo
   enddo
#ifndef _CUDA
!$OMP END PARALLEL DO
#endif

end subroutine
end module
