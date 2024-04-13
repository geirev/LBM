module m_stress
! check Eq 3.6
contains
   subroutine stress(f,feq,cxs,cys,czs,sigma,tau,dx)
   use mod_dimensions
   implicit none
   real,    intent(in) :: f(0:nx+1,0:ny+1,0:nz+1,nl)
   real,    intent(in) :: feq(0:nx+1,0:ny+1,0:nz+1,nl)
   real,    intent(out):: sigma(3,3,nx,ny,nz)
   integer, intent(in) :: cxs(nl)
   integer, intent(in) :: cys(nl)
   integer, intent(in) :: czs(nl)
   real,    intent(in) :: tau
   real,    intent(in) :: dx
   integer i,j,k,l
   real scaling,sigmaave
   real, parameter :: smagorinsky=0.08

   scaling=-(1.0 - 1.0/(2.0*tau))
   print '(a,2f12.5)','scaling=',scaling,tau

   sigma=0.0
   do l=1,nl
   do k=1,nz
   do j=1,ny
   do i=1,nx
      sigma(1,1,i,j,k)=sigma(1,1,i,j,k)+cxs(l)*cxs(l)*(f(i,j,k,l)-feq(i,j,k,l))
      sigma(2,1,i,j,k)=sigma(2,1,i,j,k)+cys(l)*cxs(l)*(f(i,j,k,l)-feq(i,j,k,l))
      sigma(3,1,i,j,k)=sigma(3,1,i,j,k)+czs(l)*cxs(l)*(f(i,j,k,l)-feq(i,j,k,l))
      sigma(2,2,i,j,k)=sigma(2,2,i,j,k)+cys(l)*cys(l)*(f(i,j,k,l)-feq(i,j,k,l))
      sigma(3,2,i,j,k)=sigma(3,2,i,j,k)+czs(l)*cys(l)*(f(i,j,k,l)-feq(i,j,k,l))
      sigma(3,3,i,j,k)=sigma(3,3,i,j,k)+czs(l)*cys(l)*(f(i,j,k,l)-feq(i,j,k,l))
   enddo
   enddo
   enddo
   enddo
   sigma(1,2,:,:,:)=sigma(2,1,:,:,:)
   sigma(1,3,:,:,:)=sigma(3,1,:,:,:)
   sigma(2,3,:,:,:)=sigma(3,2,:,:,:)

   sigma=sigma*scaling

   print *,'Stress tensor:'
   print '(3f13.6)',sigma(:,:,19,50,25)
   print '(3f13.6)',sigma(:,:,nx/2,ny/2,nz/2)
   sigmaave=0.0
   do j=1,3
   do i=1,3
      sigmaave=sigmaave+sigma(i,j,19,50,25)*sigma(i,j,19,50,25)
   enddo
   enddo
   sigmaave=sqrt(2.0*sigmaave)
   print '(2(a,f12.5))','sigave=',sigmaave,'  eddyvisc=',(smagorinsky*dx)**2*sigmaave

   end subroutine
end module
