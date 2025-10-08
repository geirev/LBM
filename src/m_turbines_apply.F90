module m_turbines_apply
contains
subroutine turbines_apply(f,df,tau)
   use mod_dimensions,  only : nx,ny,nz
   use mod_D3Q27setup,  only : nl
   use m_turbines_init, only : ieps
   use m_readinfile   , only : ipos,nturbines,iforce
#ifdef _CUDA
   use m_readinfile,    only : ntx,nty,ntz
#endif
   use m_turbines_apply_kernel
   use m_wtime
   implicit none
   real, intent(inout) :: f(nl,0:nx+1,0:ny+1,0:nz+1)        ! distribution
   real, intent(in)    :: df(nl,-ieps:ieps,ny,nz,nturbines) ! forcing distributions
   real, intent(in)    :: tau(0:nx+1,0:ny+1,0:nz+1)         ! Tau
#ifdef _CUDA
   attributes(device) :: f
   attributes(device) :: df
   attributes(device) :: tau
#endif
   integer n,ip

#ifdef _CUDA
   integer :: tx, ty, tz, bx, by, bz
#endif
   integer, parameter :: icpu=8
   call cpustart()
#ifdef _CUDA
   tx=16; bx=((2*ieps+1)+tx-1)/tx
   ty=8;  by=(ny+ty-1)/ty
   tz=1;  bz=(nz+tz-1)/tz
#endif
   do n=1,nturbines
      ip=ipos(n)
      call turbines_apply_kernel&
#ifdef _CUDA
          &<<<dim3(bx,by,bz), dim3(tx,ty,tz)>>>&
#endif
          &(f,df,tau,ip,n,iforce,nturbines)
   enddo
   call cpufinish(icpu)
end subroutine
end module
