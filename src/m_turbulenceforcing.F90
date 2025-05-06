module m_turbulenceforcing
   integer, parameter     :: iturb_pos=10
   integer, parameter     :: iturb_radius=2
contains
subroutine turbulenceforcing(turb_df,rho,u,v,w,uu,vv,ww,it)
   use mod_dimensions
   use m_fequilscalar
   use m_wtime
   use m_readinfile, only : lturb

   real, intent(out)      :: turb_df(nl,-ieps:ieps,ny,nz)      ! forcing distributions
   real, intent(inout)    :: rho(nx,ny,nz)                     ! density
   real, intent(inout)    :: u(nx,ny,nz)                       ! velocity
   real, intent(inout)    :: v(nx,ny,nz)                       ! velocity
   real, intent(inout)    :: w(nx,ny,nz)                       ! velocity
   real, intent(in)       :: uu(ny,nz,0:nrturb)
   real, intent(in)       :: vv(ny,nz,0:nrturb)
   real, intent(in)       :: ww(ny,nz,0:nrturb)
   integer, intent(in)    :: it


   real dfeq(nl)
   real :: utmp(ny,nz)
   real :: vtmp(ny,nz)
   real :: wtmp(ny,nz)

   integer lit,i,j,k,ip

   if (.not.lturb) then
      turb_df=0.0
      return
   endif

   ip=iturb_pos

   lit=mod(it,nrturb)
   if (lit==0) lit=nrturb
   do k=1,nz
   do j=1,ny
      utmp(j,k)=0.00001*uu(j,k,lit)
      vtmp(j,k)=0.00001*vv(j,k,lit)
      wtmp(j,k)=0.00001*ww(j,k,lit)
   enddo
   enddo

! Computing the S_i term returned in df
!$OMP PARALLEL DO PRIVATE(i,j,k,dfeq) SHARED(turb_df, rho, u, v, w, ip )
   do k=1,nz
   do j=1,ny
      do i=-iturb_radius,iturb_radius
         dfeq(:)         =fequilscalar(rho(ip+i,j,k), u(ip+i,j,k),           v(ip+i,j,k),           w(ip+i,j,k))
         turb_df(:,i,j,k)=fequilscalar(rho(ip+i,j,k), u(ip+i,j,k)+utmp(j,k), v(ip+i,j,k)+vtmp(j,k), w(ip+i,j,k)+wtmp(j,k))
         turb_df(:,i,j,k)=turb_df(:,i,j,k)-dfeq(:)
      enddo
   enddo
   enddo
!$OMP END PARALLEL DO

end subroutine
end module
