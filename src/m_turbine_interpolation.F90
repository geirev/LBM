!==============================================================
!  m_turbine_interpolation.F90
!  Velocity & density interpolation (CPU/GPU unified)
!==============================================================
module m_turbine_interpolation
   use mod_dimensions, only : nx, ny, nz
   use mod_turbines,   only : pi
   implicit none
contains

!--------------------------------------------------------------
!  function turbine_trilinear
!
!  PURPOSE:
!    Trilinear interpolation of a scalar given its 8 corner
!    values and local fractional offsets fx, fy, fz.
!
!  NOTE:
!    Marked host/device so it can be called from both CPU and
!    CUDA kernels.
!--------------------------------------------------------------
#ifdef _CUDA
attributes(host,device) &
#endif
real function turbine_trilinear(c000,c100,c010,c110,c001,c101,c011,c111,fx,fy,fz)
   implicit none
   real, intent(in) :: c000,c100,c010,c110,c001,c101,c011,c111
   real, intent(in) :: fx,fy,fz

   turbine_trilinear = c000*(1-fx)*(1-fy)*(1-fz) + c100*fx*(1-fy)*(1-fz) + &
                       c010*(1-fx)*fy*(1-fz)     + c110*fx*fy*(1-fz)     + &
                       c001*(1-fx)*(1-fy)*fz     + c101*fx*(1-fy)*fz     + &
                       c011*(1-fx)*fy*fz         + c111*fx*fy*fz
end function turbine_trilinear


!--------------------------------------------------------------
!  subroutine turbine_interpolate_velocity
!
!  PURPOSE:
!    Trilinear interpolation of u, v, w, rho at a global
!    actuator point (xg, yg, zg). Works on both CPU and GPU.
!
!  ARGUMENTS:
!    u,v,w,rho : staggered fields on tile [0:nx+1,0:ny+1,0:nz+1]
!    xg,yg,zg  : global coordinates (in grid units)
!    j_start   : global j-index of local tile j=1
!    ux,uy,uz  : interpolated velocity components
!    dens      : interpolated density
!
!  NOTE:
!    - On CPU, j_start is from m_mpi_decomp_init or 1 if no MPI.
!    - On GPU, caller must pass the correct j_start.
!--------------------------------------------------------------
#ifdef _CUDA
attributes(host,device) &
#endif
subroutine turbine_interpolate_velocity(u,v,w,rho, xg,yg,zg, j_start, ux,uy,uz,dens)
   use mod_dimensions, only : nx, ny, nz
   implicit none

   real, intent(in) :: u  (0:nx+1,0:ny+1,0:nz+1)
   real, intent(in) :: v  (0:nx+1,0:ny+1,0:nz+1)
   real, intent(in) :: w  (0:nx+1,0:ny+1,0:nz+1)
   real, intent(in) :: rho(0:nx+1,0:ny+1,0:nz+1)

   real, intent(in)  :: xg, yg, zg
   integer, value    :: j_start
   real, intent(out) :: ux, uy, uz, dens

   integer :: ig, jg, kg
   integer :: i, j, k
   real    :: fx, fy, fz
   real    :: c000,c100,c010,c110,c001,c101,c011,c111

   ! Global integer indices
   ig = floor(xg)
   jg = floor(yg)
   kg = floor(zg)

   ! Fractional offsets
   fx = xg - real(ig)
   fy = yg - real(jg)
   fz = zg - real(kg)

   ! Global -> local tile index in j
   j = jg - j_start + 1
   i = ig
   k = kg

   ! Clamp into interior region
   i = max(1, min(nx-1, i))
   j = max(1, min(ny-1, j))
   k = max(1, min(nz-1, k))

   ! u
   c000=u(i,j,k);   c100=u(i+1,j,k);   c010=u(i,j+1,k);   c110=u(i+1,j+1,k)
   c001=u(i,j,k+1); c101=u(i+1,j,k+1); c011=u(i,j+1,k+1); c111=u(i+1,j+1,k+1)
   ux = turbine_trilinear(c000,c100,c010,c110,c001,c101,c011,c111,fx,fy,fz)

   ! v
   c000=v(i,j,k);   c100=v(i+1,j,k);   c010=v(i,j+1,k);   c110=v(i+1,j+1,k)
   c001=v(i,j,k+1); c101=v(i+1,j,k+1); c011=v(i,j+1,k+1); c111=v(i+1,j+1,k+1)
   uy = turbine_trilinear(c000,c100,c010,c110,c001,c101,c011,c111,fx,fy,fz)

   ! w
   c000=w(i,j,k);   c100=w(i+1,j,k);   c010=w(i,j+1,k);   c110=w(i+1,j+1,k)
   c001=w(i,j,k+1); c101=w(i+1,j,k+1); c011=w(i,j+1,k+1); c111=w(i+1,j+1,k+1)
   uz = turbine_trilinear(c000,c100,c010,c110,c001,c101,c011,c111,fx,fy,fz)

   ! rho
   c000=rho(i,j,k);   c100=rho(i+1,j,k);   c010=rho(i,j+1,k);   c110=rho(i+1,j+1,k)
   c001=rho(i,j,k+1); c101=rho(i+1,j,k+1); c011=rho(i,j+1,k+1); c111=rho(i+1,j+1,k+1)
   dens = turbine_trilinear(c000,c100,c010,c110,c001,c101,c011,c111,fx,fy,fz)
end subroutine turbine_interpolate_velocity

end module m_turbine_interpolation
