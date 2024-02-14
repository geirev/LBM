module m_velocity
! Computes density as a sum over particles with different velocities
contains
function velocity(f,rho,cs,blanking) result(vel)
   use mod_dimensions
   real,    intent(in) :: f(0:nx+1,ny,nz,nl)
   real,    intent(in) :: rho(nx,ny,nz)
   integer, intent(in) :: cs(nl)
   logical, intent(in) :: blanking(nx,ny,nz)
   real vel(nx,ny,nz)
   integer i,j,k,l

   vel=0.0
   do l=1,nl
      if (cs(l) /=0) then
         do k=1,nz
!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(i,j) SHARED(k, l, vel, f, cs)
         do j=1,ny
         do i=1,nx
            vel(i,j,k) =  vel(i,j,k) + f(i,j,k,l)*real(cs(l))
         enddo
         enddo
!$OMP END PARALLEL DO
         enddo
      endif
   enddo

   do k=1,nz
!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(i,j) SHARED(k, vel, cs, rho, blanking)
      do j=1,ny
      do i=1,nx
         if (blanking(i,j,k)) then
            vel(i,j,k)=0.0
         else
            vel(i,j,k) =  vel(i,j,k)/rho(i,j,k)
         endif
      enddo
      enddo
!$OMP END PARALLEL DO
   enddo

end function
end module
