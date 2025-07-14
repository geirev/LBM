module m_inipert
contains
subroutine inipert(rho,u,v,w,uvel)
   use mod_dimensions
   use m_readinfile, only : rho0,linipert,udir
   implicit none
   real, intent(inout)  :: rho(nx,ny,nz)
   real, intent(inout)  :: u(nx,ny,nz)
   real, intent(inout)  :: v(nx,ny,nz)
   real, intent(inout)  :: w(nx,ny,nz)
   real, intent(in)     :: uvel(nz)
#ifdef _CUDA
   attributes(device) :: rho
   attributes(device) :: u
   attributes(device) :: v
   attributes(device) :: w
   attributes(device) :: uvel
#endif

   real :: stddev=0.000010  ! value running stable in uniform flow for 2000 timesteps
   real, parameter :: pi=3.1415927410125732
   integer i,j,k
   real, allocatable :: rho_h(:,:,:)

   allocate(rho_h(nx, ny, nz))
   call random_number(rho_h)
   rho = rho_h  ! Copy to managed or device array

#ifdef _CUDA
!$cuf kernel do(3) <<<*,*>>>
#else
!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(i, j, k) SHARED(rho, u, v, w, rho0, uvel, stddev, udir)
#endif
   do k=1,nz
      do j=1,ny
         do i=1,nx
            u(i,j,k)=uvel(k)*cos(udir*pi/180.0)
            v(i,j,k)=uvel(k)*sin(udir*pi/180.0)
            w(i,j,k)=0.0
            rho(i,j,k)=rho0 + stddev*rho(i,j,k)
         enddo
      enddo
   enddo


end subroutine
end module
