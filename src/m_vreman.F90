module m_vreman
!  Vreman (2004) subgridscale turbulence model
contains

subroutine vreman(f, tau, eddyvisc ,Bbeta ,alphamag ,alpha ,beta, it, nt1)
   use mod_dimensions
   use mod_D3Q27setup
   use m_readinfile, only : ivreman,kinevisc,p2l,smagorinsky,tauin
   use m_wtime
   use m_vreman_alpha_kernel
   use m_vreman_alphamag_kernel
   use m_vreman_beta_kernel
   use m_vreman_Bbeta_kernel
   use m_vreman_tau_kernel
   implicit none
   real, intent(in)      :: f(nl,0:nx+1,0:ny+1,0:nz+1) ! Nonequilibrium f as input
   real, intent(out)     :: tau(nx,ny,nz)              ! Tau including subgrid scale mixing
   integer, intent(in)   :: it
   integer, intent(in)   :: nt1
#ifdef _CUDA
   attributes(device) :: f
   attributes(device) :: tau
#endif

   real          :: dx              ! length scale lattice to physical
   real          :: const           ! c in Vreman 2004 Eq (5)
   real :: tmp
   real tau_h,eddyvisc_h,Bbeta_h,alphamag_h


   integer :: i, j, k, l, m, p, q

   real, intent(out) :: eddyvisc(nx,ny,nz)  ! nu in Vreman 2004 Eq (5)
   real, intent(out) :: Bbeta(nx,ny,nz)     ! B_beta in Vreman 2004 Eq (5)
   real, intent(out) :: alphamag(nx,ny,nz)
   real, intent(out) :: alpha(3,3,nx,ny,nz)
   real, intent(out) :: beta(3,3,nx,ny,nz)
#ifdef _CUDA
   attributes(device) :: alpha
   attributes(device) :: beta
   attributes(device) :: eddyvisc
   attributes(device) :: Bbeta
   attributes(device) :: alphamag
#endif

   integer, parameter :: icpu=6
   integer :: tx, ty, tz, bx, by, bz
   call cpustart()

   if (ivreman /= 1) then
#ifdef _CUDA
!$cuf kernel do(3) <<<*,*>>>
#else
!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(i, j, k) SHARED(tau, kinevisc)
#endif
      do k=1,nz
         do j=1,ny
            do i=1,nx
               tau(i,j,k) = 3.0*kinevisc + 0.5
            enddo
         enddo
      enddo
#ifndef _CUDA
!$OMP END PARALLEL DO
#endif
      return
   endif

   const=2.5*smagorinsky**2
   dx=p2l%length
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Compute alpha
!@cuf istat = cudaDeviceSynchronize()
      t0 = wtime()
#ifdef _CUDA
      tx=ntx; bx=(nx+2+tx-1)/tx
      ty=nty; by=(ny+2+ty-1)/ty
      tz=ntz; bz=(nz+2+tz-1)/tz
#endif
      call vreman_alpha_kernel&
#ifdef _CUDA
        &<<<dim3(bx,by,bz), dim3(tx,ty,tz)>>>&
#endif
        &(f, H2, alpha, nx+2, ny+2, nz+2, nl)
!@cuf istat = cudaDeviceSynchronize()
   t1 = wtime(); walltimelocal(41)=walltimelocal(41)+t1-t0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Compute alphamag
!@cuf istat = cudaDeviceSynchronize()
      t0 = wtime()
#ifdef _CUDA
      tx=ntx; bx=(nx+tx-1)/tx
      ty=nty; by=(ny+ty-1)/ty
      tz=ntz; bz=(nz+tz-1)/tz
#endif
      call vreman_alphamag_kernel&
#ifdef _CUDA
        &<<<dim3(bx,by,bz), dim3(tx,ty,tz)>>>&
#endif
        &(alphamag, alpha, nx, ny, nz, nl)
!@cuf istat = cudaDeviceSynchronize()
   t1 = wtime(); walltimelocal(42)=walltimelocal(42)+t1-t0



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Compute beta
!@cuf istat = cudaDeviceSynchronize()
      t0 = wtime()
#ifdef _CUDA
      tx=ntx; bx=(nx+tx-1)/tx
      ty=nty; by=(ny+ty-1)/ty
      tz=ntz; bz=(nz+tz-1)/tz
#endif
      call vreman_beta_kernel&
#ifdef _CUDA
        &<<<dim3(bx,by,bz), dim3(tx,ty,tz)>>>&
#endif
        &(beta, alpha, nx, ny, nz, nl)
!@cuf istat = cudaDeviceSynchronize()
   t1 = wtime(); walltimelocal(43)=walltimelocal(43)+t1-t0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Compute Bbeta
! Vreman 2004 Eq (8)

!@cuf istat = cudaDeviceSynchronize()
      t0 = wtime()
#ifdef _CUDA
      tx=ntx; bx=(nx+tx-1)/tx
      ty=nty; by=(ny+ty-1)/ty
      tz=ntz; bz=(nz+tz-1)/tz
#endif
      call vreman_Bbeta_kernel&
#ifdef _CUDA
        &<<<dim3(bx,by,bz), dim3(tx,ty,tz)>>>&
#endif
        &(Bbeta, beta, nx, ny, nz, nl)
!@cuf istat = cudaDeviceSynchronize()
   t1 = wtime(); walltimelocal(44)=walltimelocal(44)+t1-t0


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Vreman 2004 Eq (5)
!@cuf istat = cudaDeviceSynchronize()
      t0 = wtime()
#ifdef _CUDA
      tx=ntx; bx=(nx+tx-1)/tx
      ty=nty; by=(ny+ty-1)/ty
      tz=ntz; bz=(nz+tz-1)/tz
#endif
      call vreman_tau_kernel&
#ifdef _CUDA
        &<<<dim3(bx,by,bz), dim3(tx,ty,tz)>>>&
#endif
        &(tau, eddyvisc, Bbeta, alphamag, kinevisc, const, nx, ny, nz)
!@cuf istat = cudaDeviceSynchronize()
   t1 = wtime(); walltimelocal(45)=walltimelocal(45)+t1-t0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   call cpufinish(icpu)
   if (it==nt1) then
      do j=41,45
         print '(a24,i3,g13.5)','vreman:',j,walltimelocal(j)
      enddo
      print '(a24,g13.5)',      'vreman:',sum(walltimelocal(41:45))
   endif

end subroutine
end module

