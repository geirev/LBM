module m_boundary_mass_monitor_kernel
contains

#ifdef _CUDA
attributes(global) &
#endif
subroutine boundary_mass_monitor_kernel(rho,u,v,diag, &
                                        j0_is_phys,jN_is_phys)

!-----------------------------------------------------------------------
! Monitor total mass and mass fluxes through the four i-j boundaries.
!
! diag(1) = total mass in the local computational domain
!
! diag(2) = outward mass flux through i = 0
! diag(3) = outward mass flux through i = nx+1
!
! diag(4) = outward mass flux through j = 0
! diag(5) = outward mass flux through j = ny+1
!
! Sign convention:
!
!       positive = mass leaving the computational domain
!
! Therefore:
!
!       i=0      outward normal = (-1,0,0)
!       i=nx+1   outward normal = (+1,0,0)
!
!       j=0      outward normal = (0,-1,0)
!       j=ny+1   outward normal = (0,+1,0)
!
! giving
!
!       Qi0 = -rho*u
!       QiN = +rho*u
!
!       Qj0 = -rho*v
!       QjN = +rho*v
!
! The fluxes are evaluated at the first interior lattice planes:
!
!       i = 1, nx
!       j = 1, ny
!
! rather than in the ghost cells.
!
! In MPI decomposition along j:
!
!   - every rank contributes to i-face fluxes;
!   - only the physical south/north ranks contribute to j-face fluxes.
!
! After this kernel, perform MPI_SUM on diag(:) if running with MPI.
!
!-----------------------------------------------------------------------

#ifdef _CUDA
   use cudafor
#endif

   use mod_dimensions, only : nx,ny,nz

   implicit none

   real, intent(in) :: rho(0:nx+1,0:ny+1,0:nz+1)
   real, intent(in) :: u  (0:nx+1,0:ny+1,0:nz+1)
   real, intent(in) :: v  (0:nx+1,0:ny+1,0:nz+1)

   real, intent(inout) :: diag(5)

   logical, value :: j0_is_phys
   logical, value :: jN_is_phys

   integer :: i,j,k

   real :: rr
   real :: flux
#ifdef _CUDA
   real :: old
#endif


!=======================================================================
! CUDA
!=======================================================================

#ifdef _CUDA

   i = threadIdx%x + (blockIdx%x-1)*blockDim%x
   j = threadIdx%y + (blockIdx%y-1)*blockDim%y
   k = threadIdx%z + (blockIdx%z-1)*blockDim%z

   if (i < 1 .or. i > nx) return
   if (j < 1 .or. j > ny) return
   if (k < 1 .or. k > nz) return


!-----------------------------------------------------------------------
! Total mass.
!-----------------------------------------------------------------------

   rr = rho(i,j,k)

   old = atomicadd(diag(1),rr)


!-----------------------------------------------------------------------
! i = 0 face.
!
! Each j,k point must contribute only once, so use i=1 threads.
!-----------------------------------------------------------------------

   if (i == 1) then

      flux = -rho(1,j,k)*u(1,j,k)

      old = atomicadd(diag(2),flux)

   endif


!-----------------------------------------------------------------------
! i = nx+1 face.
!-----------------------------------------------------------------------

   if (i == nx) then

      flux = rho(nx,j,k)*u(nx,j,k)

      old = atomicadd(diag(3),flux)

   endif


!-----------------------------------------------------------------------
! j = 0 physical face.
!
! Each i,k point contributes once using j=1 threads.
!-----------------------------------------------------------------------

   if (j == 1 .and. j0_is_phys) then

      flux = -rho(i,1,k)*v(i,1,k)

      old = atomicadd(diag(4),flux)

   endif


!-----------------------------------------------------------------------
! j = ny+1 physical face.
!-----------------------------------------------------------------------

   if (j == ny .and. jN_is_phys) then

      flux = rho(i,ny,k)*v(i,ny,k)

      old = atomicadd(diag(5),flux)

   endif


!=======================================================================
! CPU/OpenMP
!=======================================================================

#else

!$OMP PARALLEL DO COLLAPSE(3) DEFAULT(NONE)                         &
!$OMP& PRIVATE(i,j,k,rr,flux)                                       &
!$OMP& SHARED(rho,u,v,diag,j0_is_phys,jN_is_phys,nx,ny,nz)

   do k = 1,nz
   do j = 1,ny
   do i = 1,nx

      rr = rho(i,j,k)

!$OMP ATOMIC UPDATE
      diag(1) = diag(1) + rr


      if (i == 1) then

         flux = -rho(1,j,k)*u(1,j,k)

!$OMP ATOMIC UPDATE
         diag(2) = diag(2) + flux

      endif


      if (i == nx) then

         flux = rho(nx,j,k)*u(nx,j,k)

!$OMP ATOMIC UPDATE
         diag(3) = diag(3) + flux

      endif


      if (j == 1 .and. j0_is_phys) then

         flux = -rho(i,1,k)*v(i,1,k)

!$OMP ATOMIC UPDATE
         diag(4) = diag(4) + flux

      endif


      if (j == ny .and. jN_is_phys) then

         flux = rho(i,ny,k)*v(i,ny,k)

!$OMP ATOMIC UPDATE
         diag(5) = diag(5) + flux

      endif

   enddo
   enddo
   enddo

!$OMP END PARALLEL DO

#endif


end subroutine boundary_mass_monitor_kernel

end module m_boundary_mass_monitor_kernel
