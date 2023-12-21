module m_density
! Computes density as a sum over particles with different velocities
contains
function density(f) result(dens)
   use mod_dimensions
   real,    intent(in) :: f(nx,ny,nl)
   real dens(nx,ny)
   dens(:,:)=sum(f(:,:,:),dim=3)
end function
end module
