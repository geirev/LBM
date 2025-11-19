module m_turbines_bounding_box
contains
subroutine turbines_bounding_box(points, np, krad)
   use mod_dimensions, only : nx,ny,nz
   use mod_turbines
   implicit none
   type(point_t), intent(in) :: points(:)
   integer, intent(in)       :: np, krad

   integer :: p, i0,j0,k0

   t_imin = nx+1;  t_imax = -1
   t_jmin = ny+1;  t_jmax = -1
   t_kmin = nz+1;  t_kmax = -1

   do p=1,np
      i0 = floor(points(p)%xg)
      j0 = floor(points(p)%yg)
      k0 = floor(points(p)%zg)

      t_imin = min(t_imin, i0-krad)
      t_imax = max(t_imax, i0+krad)

      t_jmin = min(t_jmin, j0-krad)
      t_jmax = max(t_jmax, j0+krad)

      t_kmin = min(t_kmin, k0-krad)
      t_kmax = max(t_kmax, k0+krad)
   end do

   ! clip to interior
   t_imin = max(1, t_imin)
   t_imax = min(nx, t_imax)

   t_jmin = 1  !max(1, t_jmin)
   t_jmax = ny !min(nyg, t_jmax)

   t_kmin = max(1, t_kmin)
   t_kmax = min(nz, t_kmax)
   !print *,'bounding box' ,t_imin,t_imax,t_jmin,t_jmax,t_kmin,t_kmax
end subroutine
end module
