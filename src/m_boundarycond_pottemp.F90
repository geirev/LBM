module m_boundarycond_pottemp
contains
subroutine boundarycond_pottemp(pottemp)
   use mod_dimensions
   use m_readinfile,   only : ibnd,jbnd,kbnd,iablvisc,istable
#ifdef _CUDA
   use m_readinfile, only : ntx,nty,ntz
#endif
   use m_wtime
   use m_boundary_i_periodic_kernel
   use m_boundary_j_periodic_kernel
!   use m_boundary_k_periodic_kernel

   implicit none
   real, intent(inout):: pottemp(0:nx+1,0:ny+1,0:nz+1)

   real, parameter   :: pi=3.1415927410125732
#ifdef _CUDA
   attributes(device) :: pottemp
#endif
   integer, parameter :: icpu=11
#ifdef _CUDA
   integer :: tx, ty, tz, bx, by, bz
#endif

   if (iablvisc /= 2) return

   call cpustart()
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Update sides for inflow conditions for periodic boundary conditions in i-direction
   if (ibnd==0) then
#ifdef _CUDA
         tx=1  ; bx=1
         ty=32 ; by=(ny+2+ty-1)/ty
         tz=8  ; bz=(nz+2+tz-1)/tz
#endif
         call boundary_i_periodic_kernel&
#ifdef _CUDA
         &<<<dim3(bx,by,bz), dim3(tx,ty,tz)>>>&
#endif
         &(pottemp,1)
   endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
         &(pottemp,1)
   endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! i-inflow boundary condition
   if (ibnd==1) then
! Periodic bnd cond for pottemp in unstable case
      if (istable == -1) then
#ifdef _CUDA
         tx=1  ; bx=1
         ty=32 ; by=(ny+2+ty-1)/ty
         tz=8  ; bz=(nz+2+tz-1)/tz
#endif
         call boundary_i_periodic_kernel&
#ifdef _CUDA
         &<<<dim3(bx,by,bz), dim3(tx,ty,tz)>>>&
#endif
         &(pottemp,1)
      endif
   endif

! j-inflow boundary condition
   if (jbnd==1) then
! Periodic bnd cond for pottemp in unstable case
      if (istable == -1) then
#ifdef _CUDA
         tx=1  ; bx=1
         ty=32 ; by=(ny+2+ty-1)/ty
         tz=8  ; bz=(nz+2+tz-1)/tz
#endif
         call boundary_j_periodic_kernel&
#ifdef _CUDA
         &<<<dim3(bx,by,bz), dim3(tx,ty,tz)>>>&
#endif
         &(pottemp,1)
      endif
   endif

! Update edges for inflow conditions for periodic boundary conditions in i-direction
   if (ibnd==0) then
#ifdef _CUDA
         tx=1  ; bx=1
         ty=32 ; by=(ny+2+ty-1)/ty
         tz=8  ; bz=(nz+2+tz-1)/tz
#endif
         call boundary_i_periodic_kernel&
#ifdef _CUDA
         &<<<dim3(bx,by,bz), dim3(tx,ty,tz)>>>&
#endif
         &(pottemp,1)
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
         &(pottemp,1)
   endif

   call cpufinish(icpu)

end subroutine
end module
