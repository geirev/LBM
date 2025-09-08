module m_averaging_full
   real,    dimension(:,:,:), allocatable, private, save :: uave, vave, wave
   real,    dimension(:,:,:), allocatable, private, save :: uave2, vave2, wave2
   real,    dimension(:,:,:), allocatable, private, save :: Ti


#ifdef _CUDA
   attributes(device) uave, vave, wave
   attributes(device) uave2, vave2, wave2
   attributes(device) Ti
#endif

contains
subroutine averaging_full(u,v,w,rho,lblanking,lfinal)
   use mod_dimensions
#ifdef _CUDA
   use m_readinfile, only : ntx,nty,ntz
#endif
   use m_readinfile, only : uini
   use m_averaging_full_kernel
   use m_averaging_full_kernelfin
   use m_diag
#ifdef _CUDA
   use cudafor
#endif
   implicit none
   real, intent(in)    :: u(nx,ny,nz)        ! x component of fluid velocity
   real, intent(in)    :: v(nx,ny,nz)        ! y component of fluid velocity
   real, intent(in)    :: w(nx,ny,nz)        ! z component of fluid velocity
   real, intent(in)    :: rho(nx,ny,nz)        ! z component of fluid velocity
   logical, intent(in) :: lblanking(0:nx+1,0:ny+1,0:nz+1)  ! z component of fluid velocity
#ifdef _CUDA
   attributes(device) :: u
   attributes(device) :: v
   attributes(device) :: w
   attributes(device) :: lblanking,rho
#endif
   logical, intent(in) :: lfinal

   integer, save :: ifirst=1
   integer, save :: iave=0


#ifdef _CUDA
   integer :: tx, ty, tz, bx, by, bz
#endif

   if (ifirst == 1) then
      print '(a)','Starting averaging.'
      allocate(  uave(nx,ny,nz) )
      allocate(  vave(nx,ny,nz) )
      allocate(  wave(nx,ny,nz) )
      allocate( uave2(nx,ny,nz) )
      allocate( vave2(nx,ny,nz) )
      allocate( wave2(nx,ny,nz) )
      allocate(    Ti(nx,ny,nz) )

      uave  =0.0; vave  =0.0; wave  =0.0; uave2 =0.0; vave2 =0.0; wave2 =0.0; Ti=0.0
      ifirst=0
   endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   iave=iave+1
#ifdef _CUDA
   tx=ntx; bx=(nx+tx-1)/tx
   ty=nty; by=(ny+ty-1)/ty
   tz=ntz; bz=(nz+tz-1)/tz
#endif
   call averaging_full_kernel&
#ifdef _CUDA
          &<<<dim3(bx,by,bz), dim3(tx,ty,tz)>>>&
#endif
          &(nx, ny, nz, u, v, w, uave, vave, wave, uave2, vave2, wave2 )


   if (lfinal) then
#ifdef _CUDA
      tx=ntx; bx=(nx+tx-1)/tx
      ty=nty; by=(ny+ty-1)/ty
      tz=ntz; bz=(nz+tz-1)/tz
#endif
      call averaging_full_kernelfin&
#ifdef _CUDA
          &<<<dim3(bx,by,bz), dim3(tx,ty,tz)>>>&
#endif
          &(nx, ny, nz, uave, vave, wave, uave2, vave2, wave2, Ti, uini, iave)

      call diag(2,0,rho,uave,vave,wave,lblanking,Ti)

      deallocate( uave , vave , wave , uave2 , vave2 , wave2 , Ti)
      print '(a)','Done with averaging.'
   endif

end subroutine

end module
