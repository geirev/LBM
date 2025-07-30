module m_turbulenceforcing
   integer, parameter, public :: iturb_pos=10
   integer, parameter, public :: iturb_radius=2
   real, allocatable,  public :: turbulence_df(:,:,:)

! Persistent local device arrays
   real, allocatable, private :: vel(:,:,:,:), rtmp(:,:,:)
   real, allocatable, private :: dfeq1(:,:,:,:), dfeq2(:,:,:,:)
   real, allocatable, private :: cx(:), cy(:), cz(:)

#ifdef _CUDA
   attributes(device) :: turbulence_df
   attributes(device) :: vel, rtmp
   attributes(device) :: dfeq1, dfeq2
   attributes(device) :: cx, cy, cz
#endif

contains

subroutine init_turbulenceforcing
   use mod_dimensions
   use mod_D3Q27setup
   implicit none
   integer l

   if (.not. allocated(turbulence_df)) allocate(turbulence_df(1:nl,1:ny,1:nz))
   if (.not. allocated(vel))           allocate(vel(3,1,ny,nz))
   if (.not. allocated(rtmp))          allocate(rtmp(1,ny,nz))
   if (.not. allocated(dfeq1))         allocate(dfeq1(nl,1,ny,nz))
   if (.not. allocated(dfeq2))         allocate(dfeq2(nl,1,ny,nz))
   if (.not. allocated(cx))            allocate(cx(nl))
   if (.not. allocated(cy))            allocate(cy(nl))
   if (.not. allocated(cz))            allocate(cz(nl))

   do l=1,nl
     cx(l)=real(cxs_h(l))
     cy(l)=real(cys_h(l))
     cz(l)=real(czs_h(l))
   enddo
end subroutine init_turbulenceforcing


subroutine turbulenceforcing(rho,u,v,w,uu,vv,ww,ampl,it,nt1)
   use mod_dimensions
   use m_fequilscal
   use m_fequilscalar
#ifdef _CUDA
   use cudafor
#endif
   use mod_D3Q27setup
   use m_wtime

   real, intent(inout)    :: rho(nx,ny,nz)                     ! density
   real, intent(inout)    :: u(nx,ny,nz)                       ! velocity
   real, intent(inout)    :: v(nx,ny,nz)                       ! velocity
   real, intent(inout)    :: w(nx,ny,nz)                       ! velocity
   real, intent(in)       :: uu(ny,nz,0:nrturb)
   real, intent(in)       :: vv(ny,nz,0:nrturb)
   real, intent(in)       :: ww(ny,nz,0:nrturb)
   integer, intent(in)    :: it
   integer, intent(in)    :: nt1
   real, intent(in)       :: ampl
#ifdef _CUDA
   attributes(device) :: rho
   attributes(device) :: u
   attributes(device) :: v
   attributes(device) :: w
   attributes(device) :: uu
   attributes(device) :: vv
   attributes(device) :: ww
#endif

   integer lit,i,j,k,ip,l,i1,i2
#ifdef _CUDA
   type(dim3) :: grid,tblock
#endif

   integer, parameter :: icpu=3
   real :: dff(1:nl)
   call cpustart()

!@cuf istat = cudaDeviceSynchronize()
   t0 = wallclock()

   ip=iturb_pos
   lit=mod(it,nrturb)
   if (lit==0) lit=nrturb

#ifdef _CUDA
   ii = 1  ! i is fixed at 1

   tBlock%x = 1         ! only 1 thread in x (i-direction)
   tBlock%y = 8         ! 8 threads in y-direction
   tBlock%z = 8         ! 8 threads in z-direction

   grid%x = 1           ! only one block in x (i-direction)
   grid%y = (ny + tBlock%y - 1) / tBlock%y
   grid%z = (nz + tBlock%z - 1) / tBlock%z
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
!$OMP PARALLEL DO PRIVATE(j,k) SHARED(rho, u, v, w, ip, rtmp, vel)
#endif
   do k=1,nz
   do j=1,ny
      rtmp(1,j,k)=rho(ip,j,k)
      vel(1,1,j,k)=u(ip,j,k)
      vel(2,1,j,k)=v(ip,j,k)
      vel(3,1,j,k)=w(ip,j,k)
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
   call fequilscal<<<grid,tBlock>>>(dfeq1, rtmp, vel, weights, cx, cy, cz, H2, H3, ii)
#else
!$OMP PARALLEL DO PRIVATE(j,k) SHARED(dfeq1, rtmp, vel, weights, cx, cy, cz, H2, H3 )
   do k=1,nz
   do j=1,ny
      dfeq1(:,1,j,k)=fequilscalar(rtmp(1,j,k),vel(1,1,j,k), vel(2,1,j,k), vel(3,1,j,k), weights, cx, cy, cz, H2, H3)
   enddo
   enddo
!$OMP END PARALLEL DO
#endif
!   dff(1:10)=dfeq1(1:10,1,48,48)
!   print '(a,10g13.5)','dfeq1          :',dff(1:10)
!@cuf istat = cudaDeviceSynchronize()
    t1 = wallclock(); walltimelocal(52)=walltimelocal(52)+t1-t0


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!@cuf istat = cudaDeviceSynchronize()
   t0 = wallclock()
#ifdef _CUDA
!$cuf kernel do(2) <<<*,*>>>
#else
!$OMP PARALLEL DO PRIVATE(j,k) SHARED(rho, u, v, w, ip, rtmp, vel,  uu, vv, ww, lit)
#endif
   do k=1,nz
   do j=1,ny
      rtmp(1,j,k)=rho(ip,j,k)
      vel(1,1,j,k)=vel(1,1,j,k)+ampl*uu(j,k,lit)
      vel(2,1,j,k)=vel(2,1,j,k)+ampl*vv(j,k,lit)
      vel(3,1,j,k)=vel(3,1,j,k)+ampl*ww(j,k,lit)
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
   call fequilscal<<<grid,tBlock>>>(dfeq2, rtmp, vel, weights, cx, cy, cz, H2, H3, ii)
#else
!$OMP PARALLEL DO PRIVATE(j,k) SHARED(dfeq2, rtmp, vel, weights, cx, cy, cz, H2, H3 )
   do k=1,nz
   do j=1,ny
      dfeq2(:,1,j,k)=fequilscalar(rtmp(1,j,k),vel(1,1,j,k), vel(2,1,j,k), vel(3,1,j,k),&
                     weights, cx, cy, cz, H2, H3)
   enddo
   enddo
!$OMP END PARALLEL DO
#endif
!   dff(1:10)=dfeq2(1:10,1,48,48)
!   print '(a,10g13.5)','dfeq2          :',dff(1:10)
!@cuf istat = cudaDeviceSynchronize()
   t1 = wallclock(); walltimelocal(54)=walltimelocal(54)+t1-t0


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!@cuf istat = cudaDeviceSynchronize()
   t0 = wallclock()
#ifdef _CUDA
!$cuf kernel do(2) <<<*,*>>>
#else
!$OMP PARALLEL DO PRIVATE(j,k,l) SHARED(turbulence_df, dfeq1,dfeq2)
#endif
   do k=1,nz
   do j=1,ny
      do l=1,nl
         turbulence_df(l,j,k)=dfeq2(l,1,j,k)-dfeq1(l,1,j,k)
      enddo
   enddo
   enddo
#ifndef _CUDA
!$OMP END PARALLEL DO
#endif
!@cuf istat = cudaDeviceSynchronize()
   t1 = wallclock(); walltimelocal(55)=walltimelocal(55)+t1-t0

   call cpufinish(icpu)
   if (it==nt1) then
      do j=50,55
         print '(a24,i3,g13.5)','turbulenceforcing:',j,walltimelocal(j)
      enddo
      print '(a24,g13.5)',      'turbulenceforcing:',sum(walltimelocal(50:55))
   endif

end subroutine
end module
