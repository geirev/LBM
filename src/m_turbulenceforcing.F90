module m_turbulenceforcing
   integer, parameter     :: iturb_pos=10
   integer, parameter     :: iturb_radius=2
contains
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

   real, allocatable :: cx(:)
   real, allocatable :: cy(:)
   real, allocatable :: cz(:)
   real, allocatable :: utmp(:,:)
   real, allocatable :: vtmp(:,:)
   real, allocatable :: wtmp(:,:)
   real, allocatable :: rtmp(:,:)
   real, allocatable :: dfeq1(:,:,:)
   real, allocatable :: dfeq2(:,:,:)
#ifdef _CUDA
   attributes(device) :: cx
   attributes(device) :: cy
   attributes(device) :: cz
   attributes(device) :: utmp
   attributes(device) :: vtmp
   attributes(device) :: wtmp
   attributes(device) :: rtmp
   attributes(device) :: dfeq1
   attributes(device) :: dfeq2
#endif

   real ra,ua,va,wa
#ifdef _CUDA
   attributes(device) :: ua
   attributes(device) :: va
   attributes(device) :: wa
   attributes(device) :: ra
#endif
   integer lit,i,j,k,ip,l
   integer tblock
#ifdef _CUDA
   type(dim3) :: grid
#endif

   integer, parameter :: icpu=3
   call cpustart()

!   t0 = wallclock()
   allocate(utmp(ny,nz))
   allocate(vtmp(ny,nz))
   allocate(wtmp(ny,nz))
   allocate(rtmp(ny,nz))
   allocate(dfeq1(nl,ny,nz))
   allocate(dfeq2(nl,ny,nz))
   allocate(cx(nl))
   allocate(cy(nl))
   allocate(cz(nl))
!   t1 = wallclock(); walltime(1)=walltime(1)+t1-t0

   do l=1,nl
     cx(l)=real(cxs_h(l))
     cy(l)=real(cys_h(l))
     cz(l)=real(czs_h(l))
   enddo

   ip=iturb_pos
   lit=mod(it,nrturb)
   if (lit==0) lit=nrturb

#ifdef _CUDA
   tBlock = 256
   grid%x=ny
   grid%y=nz
   grid%z=1
#endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Computing the S_i term returned in df
!   t0 = wallclock()
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
!   t1 = wallclock(); walltime(2)=walltime(2)+t1-t0


!   t0 = wallclock()
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
!   t1 = wallclock(); walltime(3)=walltime(3)+t1-t0


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   t0 = wallclock()
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
!   t1 = wallclock(); walltime(4)=walltime(4)+t1-t0

!   t0 = wallclock()
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
!   t1 = wallclock(); walltime(5)=walltime(5)+t1-t0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   t0 = wallclock()
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
!   t1 = wallclock(); walltime(6)=walltime(6)+t1-t0

!   t0 = wallclock()
   deallocate(utmp)
   deallocate(vtmp)
   deallocate(wtmp)
   deallocate(rtmp)
   deallocate(dfeq1)
   deallocate(dfeq2)
   deallocate(cx)
   deallocate(cy)
   deallocate(cz)
!   t1 = wallclock(); walltime(7)=walltime(7)+t1-t0
!   if (it==999) then
!      do j=1,7
!         print '(a,i3,g13.5)','    turbtime:',j,walltime(j)
!      enddo
!      print '(a,g13.5)',   '    turbtotal  :',sum(walltime(1:7))
!   endif
   call cpufinish(icpu)

end subroutine
end module
