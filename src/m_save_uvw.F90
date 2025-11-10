module m_save_uvw
contains
subroutine save_uvw(it,u,v,w)
   use mod_dimensions, only : nx,ny,nz
   use m_readinfile,   only : lmeasurements
   implicit none
   integer, intent(in)  :: it
   real,    intent(in)  :: u(0:nx+1,0:ny+1,0:nz+1)
   real,    intent(in)  :: v(0:nx+1,0:ny+1,0:nz+1)
   real,    intent(in)  :: w(0:nx+1,0:ny+1,0:nz+1)
#ifdef _CUDA
   attributes(device) :: u
   attributes(device) :: v
   attributes(device) :: w
#endif
   real(kind=4) :: u_h(nx,ny,nz)
   real(kind=4) :: v_h(nx,ny,nz)
   real(kind=4) :: w_h(nx,ny,nz)
   character(len=6) cit
   integer iunit
   if (.not. lmeasurements) return

   u_h(1:nx,1:ny,1:nz)=u(1:nx,1:ny,1:nz)
   v_h(1:nx,1:ny,1:nz)=v(1:nx,1:ny,1:nz)
   w_h(1:nx,1:ny,1:nz)=w(1:nx,1:ny,1:nz)

   write(cit,'(i6.6)')it-1
   open(newunit=iunit,file='uvw'//cit//'.uf',form="unformatted", status='unknown')
         write(iunit)u_h,v_h,w_h
   close(iunit)

end subroutine
end module

