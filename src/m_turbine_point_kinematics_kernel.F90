module m_turbine_point_kinematics_kernel
contains
#ifdef _CUDA
attributes(global) &
#endif
subroutine turbine_point_kinematics_kernel(points,np,rho,u,v,w,j_start,j_end,kin)
   use mod_dimensions, only : nx,ny,nz
   use mod_turbines, only : point_t
   use m_turbine_interpolate_velocity, only : turbine_interpolate_velocity
   use m_turbine_rotor_basis, only : turbine_rotor_basis
   implicit none

   type(point_t) :: points(:)
   integer, value :: np,j_start,j_end
   real :: rho(0:nx+1,0:ny+1,0:nz+1)
   real :: u  (0:nx+1,0:ny+1,0:nz+1)
   real :: v  (0:nx+1,0:ny+1,0:nz+1)
   real :: w  (0:nx+1,0:ny+1,0:nz+1)
   real :: kin(3,np)

   integer :: p,ig,jg,kg
   real :: ux,uy,uz,dens
   real :: e_axis(3),e1(3),e2(3),e_tan(3)
   real :: ct,st

#ifdef _CUDA
   p=(blockIdx%x-1)*blockDim%x+threadIdx%x
   if (p<1 .or. p>np) return
#else
   do p=1,np
#endif
      kin(:,p)=0.0

      ig=floor(points(p)%xg)
      jg=floor(points(p)%yg)
      kg=floor(points(p)%zg)

      if (ig<1 .or. ig>nx .or. kg<1 .or. kg>nz .or. &
          jg<j_start .or. jg>j_end) then
#ifdef _CUDA
         return
#else
         cycle
#endif
      endif

      call turbine_interpolate_velocity(u,v,w,rho, &
           points(p)%xg,points(p)%yg,points(p)%zg,j_start, &
           ux,uy,uz,dens)

      call turbine_rotor_basis(points(p)%yaw,points(p)%tilt,e_axis,e1,e2)
      ct=cos(points(p)%theta)
      st=sin(points(p)%theta)
      e_tan=-st*e1+ct*e2

      kin(1,p)=ux*e_axis(1)+uy*e_axis(2)+uz*e_axis(3)
      kin(2,p)=points(p)%omegand*points(p)%relm - &
               (ux*e_tan(1)+uy*e_tan(2)+uz*e_tan(3))
      kin(3,p)=dens
#ifndef _CUDA
   enddo
#endif
end subroutine turbine_point_kinematics_kernel
end module m_turbine_point_kinematics_kernel
