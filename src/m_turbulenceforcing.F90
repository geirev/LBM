module m_turbulenceforcing
   integer, parameter     :: iturb_pos=10
   integer, parameter     :: iturb_radius=2

 ! Persistent device arrays
   real, allocatable :: utmp(:,:), vtmp(:,:), wtmp(:,:), rtmp(:,:)
   real, allocatable :: dfeq1(:,:,:), dfeq2(:,:,:)
   real, allocatable :: cx(:), cy(:), cz(:)
   real, allocatable :: turb_df(:,:,:)

#ifdef _CUDA
   attributes(device) :: utmp, vtmp, wtmp, rtmp
   attributes(device) :: dfeq1, dfeq2
   attributes(device) :: cx, cy, cz
   attributes(device) :: turb_df
#endif

contains

subroutine init_turbulenceforcing
   use mod_dimensions
   use mod_D3Q27setup
   implicit none
   integer l

   if (.not. allocated(turb_df)) allocate(turb_df(1:nl,1:ny,1:nz))
   if (.not. allocated(utmp))    allocate(utmp(ny,nz))
   if (.not. allocated(vtmp))    allocate(vtmp(ny,nz))
   if (.not. allocated(wtmp))    allocate(wtmp(ny,nz))
   if (.not. allocated(rtmp))    allocate(rtmp(ny,nz))
   if (.not. allocated(dfeq1))   allocate(dfeq1(nl,ny,nz))
   if (.not. allocated(dfeq2))   allocate(dfeq2(nl,ny,nz))
   if (.not. allocated(cx))      allocate(cx(nl))
   if (.not. allocated(cy))      allocate(cy(nl))
   if (.not. allocated(cz))      allocate(cz(nl))

   do l=1,nl
     cx(l)=real(cxs_h(l))
     cy(l)=real(cys_h(l))
     cz(l)=real(czs_h(l))
   enddo
end subroutine init_turbulenceforcing


subroutine turbulenceforcing(turb_df,rho,u,v,w,uu,vv,ww,it)
   use mod_dimensions
   use m_fequilscal
   use m_fequilscalar
#ifdef _CUDA
   use cudafor
#endif
   use mod_D3Q27setup
   use m_wtime

   real, intent(out)      :: turb_df(nl,ny,nz)                 ! forcing distributions
   real, intent(inout)    :: rho(nx,ny,nz)                     ! density
   real, intent(inout)    :: u(nx,ny,nz)                       ! velocity
   real, intent(inout)    :: v(nx,ny,nz)                       ! velocity
   real, intent(inout)    :: w(nx,ny,nz)                       ! velocity
   real, intent(in)       :: uu(ny,nz,0:nrturb)
   real, intent(in)       :: vv(ny,nz,0:nrturb)
   real, intent(in)       :: ww(ny,nz,0:nrturb)
   integer, intent(in)    :: it
#ifdef _CUDA
   attributes(device) :: turb_df
   attributes(device) :: rho
   attributes(device) :: u
   attributes(device) :: v
   attributes(device) :: w
   attributes(device) :: uu
   attributes(device) :: vv
   attributes(device) :: ww
#endif

   integer lit,i,j,k,ip,l
   integer tblock
#ifdef _CUDA
   type(dim3) :: grid
#endif

   integer, parameter :: icpu=3
   call cpustart()

!@cuf istat = cudaDeviceSynchronize()
   t0 = wallclock()

   ip=iturb_pos
   lit=mod(it,nrturb)
   if (lit==0) lit=nrturb

#ifdef _CUDA
   tBlock = 256
   grid%x=ny
   grid%y=nz
   grid%z=1
#endif
!@cuf istat = cudaDeviceSynchronize()
   t1 = wallclock(); walltimelocal(50)=walltimelocal(50)+t1-t0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Computing the S_i term returned in df
!@cuf istat = cudaDeviceSynchronize()
      t0 = wallclock()
#ifdef _CUDA
!$cuf kernel do(2) <<<*,*>>>
#else
!$OMP PARALLEL DO PRIVATE(j,k) SHARED(rho, u, v, w, ip, rtmp, utmp, vtmp, wtmp)
#endif
   do k=1,nz
   do j=1,ny
      rtmp(j,k)=rho(ip,j,k)
      utmp(j,k)=u(ip,j,k)
      vtmp(j,k)=v(ip,j,k)
      wtmp(j,k)=w(ip,j,k)
   enddo
   enddo
#ifndef _CUDA
!$OMP END PARALLEL DO
#endif
!@cuf istat = cudaDeviceSynchronize()
   t1 = wallclock(); walltimelocal(51)=walltimelocal(51)+t1-t0


!@cuf istat = cudaDeviceSynchronize()
   t0 = wallclock()
#ifdef _CUDA
   call fequilscal<<<grid,tblock>>>(dfeq1,rtmp, utmp, vtmp, wtmp, weights, cx, cy, cz, H2, H3)
#else
!$OMP PARALLEL DO PRIVATE(j,k) SHARED(dfeq1, rtmp, utmp, vtmp, wtmp, weights, cx, cy, cz, H2, H3 )
   do k=1,nz
   do j=1,ny
      dfeq1(:,j,k)=fequilscalar(rtmp(j,k), utmp(j,k), vtmp(j,k), wtmp(j,k), weights, cx, cy, cz, H2, H3)
   enddo
   enddo
!$OMP END PARALLEL DO
#endif
!@cuf istat = cudaDeviceSynchronize()
    t1 = wallclock(); walltimelocal(52)=walltimelocal(52)+t1-t0


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!@cuf istat = cudaDeviceSynchronize()
   t0 = wallclock()
#ifdef _CUDA
!$cuf kernel do(2) <<<*,*>>>
#else
!$OMP PARALLEL DO PRIVATE(j,k) SHARED(rho, u, v, w, ip, rtmp, utmp, vtmp, wtmp, uu, vv, ww, lit)
#endif
   do k=1,nz
   do j=1,ny
      rtmp(j,k)=rho(ip,j,k)
      utmp(j,k)=utmp(j,k)+0.001*uu(j,k,lit)
      vtmp(j,k)=vtmp(j,k)+0.001*vv(j,k,lit)
      wtmp(j,k)=wtmp(j,k)+0.001*ww(j,k,lit)
   enddo
   enddo
#ifndef _CUDA
!$OMP END PARALLEL DO
#endif
!@cuf istat = cudaDeviceSynchronize()
   t1 = wallclock(); walltimelocal(53)=walltimelocal(53)+t1-t0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!@cuf istat = cudaDeviceSynchronize()
   t0 = wallclock()
#ifdef _CUDA
   call fequilscal<<<grid,tblock>>>(dfeq2,rtmp, utmp, vtmp, wtmp, weights, cx, cy, cz, H2, H3)
#else
!$OMP PARALLEL DO PRIVATE(j,k) SHARED(dfeq2, rtmp, utmp, vtmp, wtmp, weights, cx, cy, cz, H2, H3 )
   do k=1,nz
   do j=1,ny
      dfeq2(:,j,k)=fequilscalar(rtmp(j,k), utmp(j,k), vtmp(j,k), wtmp(j,k), weights, cx, cy, cz, H2, H3)
   enddo
   enddo
!$OMP END PARALLEL DO
#endif
!@cuf istat = cudaDeviceSynchronize()
   t1 = wallclock(); walltimelocal(54)=walltimelocal(54)+t1-t0


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!@cuf istat = cudaDeviceSynchronize()
   t0 = wallclock()
#ifdef _CUDA
!$cuf kernel do(2) <<<*,*>>>
#else
!$OMP PARALLEL DO PRIVATE(j,k,l) SHARED(turb_df, dfeq1,dfeq2)
#endif
   do k=1,nz
   do j=1,ny
      do l=1,nl
         turb_df(l,j,k)=dfeq2(l,j,k)-dfeq1(l,j,k)
      enddo
   enddo
   enddo
#ifndef _CUDA
!$OMP END PARALLEL DO
#endif
!@cuf istat = cudaDeviceSynchronize()
   t1 = wallclock(); walltimelocal(55)=walltimelocal(55)+t1-t0

   call cpufinish(icpu)
   if (it==999) then
      do j=50,55
         print '(a24,i3,g13.5)','turbulenceforcing:',j,walltimelocal(j)
      enddo
      print '(a24,g13.5)',      'turbulenceforcing:',sum(walltimelocal(50:55))
   endif

end subroutine
end module
