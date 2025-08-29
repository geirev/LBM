module m_inflow_turbulence_compute
contains
subroutine inflow_turbulence_compute(uu,vv,ww,rr,lfirst,nrturb)
   use mod_dimensions
   use m_pseudo2D
   implicit none
   integer, intent(in)  :: nrturb
   real, intent(inout)  :: uu(ny,nz,0:nrturb)
   real, intent(inout)  :: vv(ny,nz,0:nrturb)
   real, intent(inout)  :: ww(ny,nz,0:nrturb)
   real, intent(inout)  :: rr(ny,nz,0:nrturb)
#ifdef _CUDA
   attributes(device) :: uu
   attributes(device) :: vv
   attributes(device) :: ww
   attributes(device) :: rr
#endif
   real, allocatable, dimension(:,:,:) :: uu_h
   real, allocatable, dimension(:,:,:) :: vv_h
   real, allocatable, dimension(:,:,:) :: ww_h
   real, allocatable, dimension(:,:,:) :: rr_h
   logical, intent(in)  :: lfirst

   real cor1,cor2,dx,dy,dir
   integer(kind=4) n1,n2
   integer i,j,k
   logical :: verbose=.false.
   real aveu,avev,avew,aver
   real varu,varv,varw,varr

   real :: timecor=0.98

   allocate( uu_h(ny,nz,0:nrturb) )
   allocate( vv_h(ny,nz,0:nrturb) )
   allocate( ww_h(ny,nz,0:nrturb) )
   allocate( rr_h(ny,nz,0:nrturb) )

! Simulating a time series of inflow boundary perturbations for u
   print '(a)','inflow_turbulence_compute: Simulating inflow turbulence forcing'
   if (lfirst) then
      uu_h(:,:,0)=0.0
      vv_h(:,:,0)=0.0
      ww_h(:,:,0)=0.0
      rr_h(:,:,0)=0.0
   else
      uu_h(:,:,0)=uu(:,:,nrturb)
      vv_h(:,:,0)=vv(:,:,nrturb)
      ww_h(:,:,0)=ww(:,:,nrturb)
      rr_h(:,:,0)=rr(:,:,nrturb)
   endif

   cor1=10.0/sqrt(3.0)
   cor2=10.0/sqrt(3.0)
   dir=0.0
   dx=1.0
   dy=1.0
   n1=ny
   n2=nz
   call pseudo2d(uu_h(:,:,0:nrturb),ny,nz,nrturb+1,cor1,cor2,dx,dy,n1,n2,dir,verbose)
   call pseudo2d(vv_h(:,:,0:nrturb),ny,nz,nrturb+1,cor1,cor2,dx,dy,n1,n2,dir,verbose)
   call pseudo2d(ww_h(:,:,0:nrturb),ny,nz,nrturb+1,cor1,cor2,dx,dy,n1,n2,dir,verbose)
   call pseudo2d(rr_h(:,:,0:nrturb),ny,nz,nrturb+1,cor1,cor2,dx,dy,n1,n2,dir,verbose)

! Imposing time correlations
   do i=1,nrturb
      uu_h(:,:,i)=timecor*uu_h(:,:,i-1)+sqrt(1.0-timecor**2)*uu_h(:,:,i)
      vv_h(:,:,i)=timecor*vv_h(:,:,i-1)+sqrt(1.0-timecor**2)*vv_h(:,:,i)
      ww_h(:,:,i)=timecor*ww_h(:,:,i-1)+sqrt(1.0-timecor**2)*ww_h(:,:,i)
      rr_h(:,:,i)=timecor*rr_h(:,:,i-1)+sqrt(1.0-timecor**2)*rr_h(:,:,i)
   enddo

! Ensure zero mean and variance equal to one
!   print *,'removing nonzero average and scaling variance'
   do i=1,nrturb
      aveu=0.0
      avev=0.0
      avew=0.0
      aver=0.0
      do k=1,nz
      do j=1,ny
         aveu = aveu + uu_h(j,k,i)
         avev = avev + vv_h(j,k,i)
         avew = avew + ww_h(j,k,i)
         aver = aver + rr_h(j,k,i)
      enddo
      enddo
      aveu=aveu/real(ny*nz)
      avev=avev/real(ny*nz)
      avew=avew/real(ny*nz)
      aver=aver/real(ny*nz)
      uu_h(:,:,i)=uu_h(:,:,i)-aveu
      vv_h(:,:,i)=vv_h(:,:,i)-avev
      ww_h(:,:,i)=ww_h(:,:,i)-avew
      rr_h(:,:,i)=rr_h(:,:,i)-aver

      varu=0.0
      varv=0.0
      varw=0.0
      varr=0.0
      do k=1,nz
      do j=1,ny
         varu = varu + uu_h(j,k,i)**2
         varv = varv + vv_h(j,k,i)**2
         varw = varw + ww_h(j,k,i)**2
         varr = varr + rr_h(j,k,i)**2
      enddo
      enddo
      varu=sqrt(varu/real(ny*nz-1))
      varv=sqrt(varv/real(ny*nz-1))
      varw=sqrt(varw/real(ny*nz-1))
      varr=sqrt(varr/real(ny*nz-1))
      uu_h(:,:,i)=(1.0/varu)*uu_h(:,:,i)
      vv_h(:,:,i)=(1.0/varv)*vv_h(:,:,i)
      ww_h(:,:,i)=(1.0/varw)*ww_h(:,:,i)
      rr_h(:,:,i)=(1.0/varr)*rr_h(:,:,i)

   enddo

!   print *,' copy to device '

! Copy to device
   uu=uu_h
   vv=vv_h
   ww=ww_h
   rr=rr_h

!   print *,' copied to device '

   if (allocated(uu_h)) deallocate(uu_h)
   if (allocated(vv_h)) deallocate(vv_h)
   if (allocated(ww_h)) deallocate(ww_h)
   if (allocated(rr_h)) deallocate(rr_h)

!   print *,' deallocated '

end subroutine
end module
