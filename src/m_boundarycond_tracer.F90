module m_boundarycond_tracer
contains
subroutine boundarycond_tracer(tracer)
   use mod_dimensions
   use m_readinfile,   only : ibnd,jbnd,kbnd,rho0,udir,iablvisc,istable
#ifdef _CUDA
   use m_readinfile, only : ntx,nty,ntz
#endif
   use m_wtime

   use m_boundary_periodic_i_kernel
   use m_boundary_periodic_j_kernel
   use m_boundary_periodic_k_kernel

   implicit none
   real, intent(inout):: tracer(ntracer,0:nx+1,0:ny+1,0:nz+1)

   real, parameter   :: pi=acos(-1.0)
#ifdef _CUDA
   attributes(device) :: tracer
#endif
   integer, parameter :: icpu=11
#ifdef _CUDA
   integer :: tx, ty, tz, bx, by, bz
#endif

   if (ntracer == 0) return
   call cpustart()

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Periodic boundary conditions in i-direction
   if (ibnd==0) then
#ifdef _CUDA
         tx=1;  bx=1
         ty=32; by=(ny+2+ty-1)/ty
         tz=8;  bz=(nz+2+tz-1)/tz
#endif
         call boundary_i_periodic_kernel&
#ifdef _CUDA
         &<<<dim3(bx,by,bz), dim3(tx,ty,tz)>>>&
#endif
         &(tracer,ntracer)
   endif

! Periodic boundary conditions in j-direction.
   if (jbnd==0) then
#ifdef _CUDA
      tx=256; bx=(nx+2+tx-1)/tx
      ty=1;   by=1
      tz=1;   bz=(nz+2+tz-1)/tz
#endif
      call boundary_j_periodic_kernel&
#ifdef _CUDA
      &<<<dim3(bx,by,bz), dim3(tx,ty,tz)>>>&
#endif
      &(tracer,ntracer)
   endif

! Periodic boundary conditions in k-direction.
   if (kbnd==0) then
#ifdef _CUDA
      tx=256; bx=(nx+2+tx-1)/tx
      ty=1;   by=1
      tz=1;   bz=(nz+2+tz-1)/tz
#endif
      call boundary_k_periodic_kernel&
#ifdef _CUDA
      &<<<dim3(bx,by,bz), dim3(tx,ty,tz)>>>&
#endif
      &(tracer,ntracer)
   endif


   call cpufinish(icpu)

end subroutine
end module
