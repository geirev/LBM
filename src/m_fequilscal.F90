module m_fequilscal
contains
#ifdef _CUDA
attributes(global) &
#endif
subroutine fequilscal(feq,rho,u,v,w,weights,cxs,cys,czs,H2,H3)
   use mod_dimensions
   implicit none
   real, intent(out)   :: feq(nl,ny,nz)
   real,    intent(in) :: rho(ny,nz)
   real,    intent(in) :: u(ny,nz)
   real,    intent(in) :: v(ny,nz)
   real,    intent(in) :: w(ny,nz)

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


   real :: A0_2(3,3)
   real :: A0_3(3,3,3)
   real :: vel(3)

   real dens
   real tmp

   integer j, k, l, p, q, r

   real, parameter :: cs2=1/3.0
   real, parameter :: cs4=1/9.0
   real, parameter :: cs6=1/27.0
   real, parameter :: inv6cs6 = 1.0/(6.0*cs6)

#ifdef _CUDA
   j= blockDim%x*(blockIdx%x-1)+threadIdx%x
   k= blockDim%y*(blockIdx%y-1)+threadIdx%y
#endif
   if ((j>ny).or.(k>nz)) return

   vel(1)=u(j,k)
   vel(2)=v(j,k)
   vel(3)=w(j,k)
   dens=rho(j,k)

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
      feq(l,j,k)=dens
      !feq(l)=feq(l) + dens*( cc(1,l)*vel(1) + cc(2,l)*vel(2) + cc(3,l)*vel(3) )/cs2
      feq(l,j,k)=feq(l,j,k) + dens*( (cxs(l))*vel(1) +(cys(l))*vel(2) + (czs(l))*vel(3) )/cs2

      do p=1,3
      do q=1,3
         feq(l,j,k)=feq(l,j,k) + H2(p,q,l)*A0_2(p,q)/(2.0*cs4)
      enddo
      enddo

      ! the above identically recovers the BGK equilibrium, below we add third order contributions
      do p=1,3
      do q=1,3
      do r=1,3
         feq(l,j,k)=feq(l,j,k) + H3(l,p,q,r)*A0_3(p,q,r)*inv6cs6
      enddo
      enddo
      enddo

      feq(l,j,k)= weights(l)*feq(l,j,k)

   enddo

end subroutine
end module
