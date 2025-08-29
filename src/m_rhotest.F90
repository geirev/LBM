module m_rhotest
contains
subroutine rhotest(f,rho,string)
   use mod_dimensions
   use mod_D3Q27setup, only : nl
   implicit none
   real,   intent(in)  :: f(nl,0:nx+1,0:ny+1,0:nz+1)
   real,   intent(out) :: rho(nx,ny,nz)
   character(len=*),   intent(in)   :: string
   real, allocatable   :: rho_h(:,:,:)
   integer i,j,k,l
   real x
#ifdef _CUDA
   attributes(device) :: f
   attributes(device) :: rho
#endif
   print '(a)',trim(string)

#ifdef _CUDA
!$cuf kernel do(3) <<<*,*>>>
#else
!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(i,j,k,l) SHARED(rho, f)
#endif
   do k=1,nz
   do j=1,ny
   do i=1,nx
      rho(i,j,k)=f(1,i,j,k)
      do l = 2, nl
         rho(i,j,k)=rho(i,j,k)+f(l,i,j,k)
      enddo
   enddo
   enddo
   enddo
#ifndef _CUDA
!$OMP END PARALLEL DO
#endif

!      if (runtest) write(98,'(a,100g15.7)')'f:',feq(1:10,ip,jp,kp)

   allocate(rho_h(nx,ny,nz))
   rho_h=rho
   do k=1,nz
   do j=1,ny
   do i=1,nx
      if (.not. ( (rho_h(i,j,k)==rho_h(i,j,k)) .and. (abs(rho_h(i,j,k)) < huge(x)) )) then
         print *, "Bad value at", i,j,k,rho_h(i,j,k)
         stop
      endif
      if (rho_h(i,j,k) < 0.0) then
         print '(2a,3I4,g13.5)',trim(string),', rho at:',i,j,k,rho_h(i,j,k)
         stop
      endif
   enddo
   enddo
   enddo
   deallocate(rho_h)
end subroutine
end module
