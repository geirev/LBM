module m_fequilscalar
contains
#ifdef _CUDA
!attributes(device) &
#endif
function fequilscalar(rho,u,v,w,weights,cxs,cys,czs,H2,H3) result(feq)
!function fequilscalar(rho, u, v, w) result(feq)
   use mod_dimensions
   implicit none
   real,    intent(in) :: rho
   real,    intent(in) :: u
   real,    intent(in) :: v
   real,    intent(in) :: w

   real, intent(in) :: cxs(nl)
   real, intent(in) :: cys(nl)
   real, intent(in) :: czs(nl)
   real,    intent(in) :: weights(nl)
   real,    intent(in) :: H2(3,3,nl)
   real,    intent(in) :: H3(3,3,3,nl)
#ifdef _CUDA
   attributes(device) :: rho
   attributes(device) :: u
   attributes(device) :: v
   attributes(device) :: w
   attributes(device) :: cxs
   attributes(device) :: cys
   attributes(device) :: czs
   attributes(device) :: weights
   attributes(device) :: H2
   attributes(device) :: H3
#endif

   real :: feq(nl)

   real :: A0_2(3,3)
   real :: A0_3(3,3,3)
   real :: vel(3)

   real dens

   integer l, p, q, r

   real, parameter :: cs2=1/3.0
   real, parameter :: cs4=1/9.0
   real, parameter :: cs6=1/27.0
   real, parameter :: inv6cs6 = 1.0/(6.0*cs6)


   vel(1)=u
   vel(2)=v
   vel(3)=w
   dens=rho

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

#ifndef _CUDA
! Equilibrium distribution \citet{fen21a} Eq. (32) or jac18a eq (27)
   do l=1,nl
      feq(l)=dens
      !feq(l)=feq(l) + dens*( cc(1,l)*vel(1) + cc(2,l)*vel(2) + cc(3,l)*vel(3) )/cs2
      feq(l)=feq(l) + dens*( cxs(l)*vel(1) +cys(l)*vel(2) + czs(l)*vel(3) )/cs2

      do p=1,3
      do q=1,3
         feq(l)=feq(l) + H2(p,q,l)*A0_2(p,q)/(2.0*cs4)
      enddo
      enddo

      ! the above identically recovers the BGK equilibrium, below we add third order contributions
      do p=1,3
      do q=1,3
      do r=1,3
         feq(l)=feq(l) + H3(p,q,r,l)*A0_3(p,q,r)*inv6cs6
      enddo
      enddo
      enddo

      feq(l)= weights(l)*feq(l)

   enddo
#endif

end function
end module
