module m_fequilscalar
contains
#ifdef _CUDA
attributes(device) &
#endif
function fequilscalar(rho, u, v, w) result(feq)
   use mod_dimensions
   use mod_D3Q27setup
   use m_readinfile, only : ibgk
   implicit none
   real,    intent(in) :: rho
   real,    intent(in) :: u
   real,    intent(in) :: v
   real,    intent(in) :: w
   real :: feq(nl)


   real  :: cc(3,nl)         ! Array storage of cxs, cys, and czs
   real  :: A0_2(3,3)
   real  :: A0_3(3,3,3)
   real  :: vel(1:3),dens
   real tmp

   integer l, p, q, r

   real, parameter :: inv6cs6 = 1.0/(6.0*cs6)

   cc(1,1:nl) = [0., 1.,-1., 0., 0., 0., 0., 1.,-1., 1.,-1.,-1., 1., 0., 0.,-1., 1., 0., 0.,-1., 1.,-1., 1., 1.,-1.,-1., 1.]
   cc(2,1:nl) = [0., 0., 0., 1.,-1., 0., 0., 1.,-1.,-1., 1., 0., 0., 1.,-1., 0., 0.,-1., 1., 1.,-1.,-1., 1., 1.,-1., 1.,-1.]
   cc(3,1:nl) = [0., 0., 0., 0., 0.,-1., 1., 0., 0., 0., 0.,-1., 1., 1.,-1., 1.,-1., 1.,-1., 1.,-1.,-1., 1.,-1., 1.,-1., 1.]


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


! Equilibrium distribution \citet{fen21a} Eq. (32) or jac18a eq (27)
   do l=1,nl
      feq(l)=dens
      feq(l)=feq(l) + dens*( cc(1,l)*vel(1) + cc(2,l)*vel(2) + cc(3,l)*vel(3) )/cs2

      do p=1,3
      do q=1,3
         feq(l)=feq(l) + H2(p,q,l)*A0_2(p,q)/(2.0*cs4)
      enddo
      enddo
      ! the above identically recovers the BGK equilibrium, below we add third order contributions
!      if (ibgk == 3) then
         do p=1,3
         do q=1,3
         do r=1,3
            feq(l)=feq(l) + H3(l,p,q,r)*A0_3(p,q,r)*inv6cs6
         enddo
         enddo
         enddo
!      endif

      feq(l)= weights(l)*feq(l)

   enddo

end function
end module
