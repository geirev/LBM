module m_heatflux
contains
subroutine heatflux(tempout,tempin)
! f enter routine in feq following collisions
! f is returned in f
   use mod_dimensions, only : nx, ny, nz
#ifdef _CUDA
   use cudafor
   use m_readinfile, only : ntx,nty,ntz
#endif
   use m_readinfile, only : iablvisc,istable,p2l
   use m_theta_fluxbc_kernel
   implicit none
   real, intent(inout) :: tempin (0:nx+1,0:ny+1,0:nz+1)
   real, intent(inout) :: tempout(0:nx+1,0:ny+1,0:nz+1)

#ifdef _CUDA
   attributes(device) :: tempin
   attributes(device) :: tempout
   integer :: tx, ty, tz, bx, by, bz
#endif
   integer i,j,k
   integer :: nsel
   real :: q
   real :: kappa
   real :: heating        ! Heating parameter for lower boundary in unstable case


   select case (istable)
   case (-1)
      q= 0.1       ! K m/s   0.05-0.1 for moderate, 0.2 -0.3 for strong
      !consider making the wall heating depend on local κ/tt to get a truly constant physical flux.
      kappa= 5.0  ! m^2/s (eddy diffusivity)
      q=q/p2l%vel
      kappa=kappa/(p2l%vel*p2l%length)
      heating=q/kappa
   case default
      return
   end select


! Apply vertical heat-flux BC on θ via the ghost cell k=0
   if (iablvisc == 2 .and. istable==-1) then
#ifdef _CUDA
      tx=ntx; bx=(nx+tx-1)/tx
      ty=nty; by=(ny+ty-1)/ty
      tz=1; bz=1
#endif
      call theta_fluxbc_kernel&
#ifdef _CUDA
        &<<<dim3(bx,by,bz), dim3(tx,ty,tz)>>>&
#endif
        &(tempin, heating)

#ifdef _CUDA
      tx=ntx; bx=(nx+tx-1)/tx
      ty=nty; by=(ny+ty-1)/ty
      tz=1; bz=1
#endif
      call theta_fluxbc_kernel&
#ifdef _CUDA
        &<<<dim3(bx,by,bz), dim3(tx,ty,tz)>>>&
#endif
        &(tempout, heating)
   endif

end subroutine
end module

