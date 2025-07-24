module m_fequilscal
contains
#ifdef _CUDA
attributes(global) &
#endif
subroutine fequilscal(feq,rho,vel,weights,cxs,cys,czs,H2,H3,ii)
   use mod_dimensions
#ifdef _CUDA
   use cudafor
#endif
   implicit none
   integer, value      :: ii
   real, intent(out)   :: feq(nl,ii,ny,nz)
   real,    intent(in) :: rho(ii,ny,nz)
   real,    intent(in) :: vel(3,ii,ny,nz)

   real, intent(in) :: cxs(nl)
   real, intent(in) :: cys(nl)
   real, intent(in) :: czs(nl)
   real,    intent(in) :: weights(nl)
   real,    intent(in) :: H2(3,3,nl)
   real,    intent(in) :: H3(3,3,3,nl)
#ifdef _CUDA
   attributes(device) :: rho
   attributes(device) :: vel
   attributes(device) :: cxs
   attributes(device) :: cys
   attributes(device) :: czs
   attributes(device) :: weights
   attributes(device) :: H2
   attributes(device) :: H3
#endif


   real :: A0_2(3,3)
   real :: A0_3(3,3,3)

   integer i, j, k, l, p, q, r

   real, parameter :: cs2=1/3.0
   real, parameter :: cs4=1/9.0
   real, parameter :: cs6=1/27.0
   real, parameter :: inv6cs6 = 1.0/(6.0*cs6)

#ifdef _CUDA
   i = threadIdx%x + (blockIdx%x - 1) * blockDim%x
   j = threadIdx%y + (blockIdx%y - 1) * blockDim%y
   k = threadIdx%z + (blockIdx%z - 1) * blockDim%z

   if (i > ii .or. j > ny .or. k > nz) return
   !print *,i,j,k
#endif

! A0_2 and A0_3 from \citet{fen21a} (following Eq. 32)
   do q=1,3
   do p=1,3
      A0_2(p,q)=rho(i,j,k)*vel(p,i,j,k)*vel(q,i,j,k)
   enddo
   enddo

   do r=1,3
   do q=1,3
   do p=1,3
      A0_3(p,q,r)=rho(i,j,k)*vel(p,i,j,k)*vel(q,i,j,k)*vel(r,i,j,k)
   enddo
   enddo
   enddo


! Equilibrium distribution \citet{fen21a} Eq. (32) or jac18a eq (27)
   do l=1,nl
      feq(l,i,j,k)=rho(i,j,k)
      feq(l,i,j,k)=feq(l,i,j,k) + rho(i,j,k)*((cxs(l))*vel(1,i,j,k)+(cys(l))*vel(2,i,j,k)+(czs(l))*vel(3,i,j,k))/cs2

      do p=1,3
      do q=1,3
         feq(l,i,j,k)=feq(l,i,j,k) + H2(p,q,l)*A0_2(p,q)/(2.0*cs4)
      enddo
      enddo

      ! the above identically recovers the BGK equilibrium, below we add third order contributions
      do p=1,3
      do q=1,3
      do r=1,3
         feq(l,i,j,k)=feq(l,i,j,k) + H3(l,p,q,r)*A0_3(p,q,r)*inv6cs6
      enddo
      enddo
      enddo

      feq(l,i,j,k)= weights(l)*feq(l,i,j,k)

   enddo

end subroutine
end module
