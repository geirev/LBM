module m_feqscalar
contains
subroutine feqscalar(feq,rho,u,v,w)
   use mod_dimensions
   use mod_D3Q27setup
   implicit none
   real,    intent(in) :: rho
   real,    intent(in) :: u
   real,    intent(in) :: v
   real,    intent(in) :: w
   real,    intent(inout):: feq(nl)
   integer l

   do l=1,nl
      feq(l) = rho*weights(l) *(1.0                                                  &
                     +3.0*(real(cxs(l))*u + real(cys(l))*v + real(czs(l))*w)         &
                     +9.0*(real(cxs(l))*u + real(cys(l))*v + real(czs(l))*w)**2/2.0  &
                     -3.0*(u**2 + v**2 + w**2)/2.0)
   enddo

end subroutine
end module
