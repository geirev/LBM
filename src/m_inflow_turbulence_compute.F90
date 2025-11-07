module m_inflow_turbulence_compute
contains
subroutine inflow_turbulence_compute(uu,vv,ww,rr,ny_here,nz,nrturb,lfirst)
   use m_pseudo2D
   use m_tecfld
   implicit none
   integer, intent(in) :: ny_here,nz,nrturb
   logical, intent(in) :: lfirst
   real, intent(inout) :: uu(ny_here,nz,0:nrturb)
   real, intent(inout) :: vv(ny_here,nz,0:nrturb)
   real, intent(inout) :: ww(ny_here,nz,0:nrturb)
   real, intent(inout) :: rr(ny_here,nz,0:nrturb)
   real :: cor1, cor2
   real :: dx=1.0, dy=1.0, dir=0.0, timecor=0.98
   integer :: i,j,k
   real :: aveu,avev,avew,aver,varu,varv,varw,varr
   integer n1,n2

   cor1=10.0/sqrt(3.0)!*real(ny_here/nz)
   cor2=10.0/sqrt(3.0)


   print *,'compute_turbulence_field: generating pseudo-2D inflow forcing'

   if (lfirst) then
      uu(:,:,0)=0.0; vv(:,:,0)=0.0; ww(:,:,0)=0.0; rr(:,:,0)=0.0
   endif
   n1=ny_here
   n2=nz
   call pseudo2d(uu(:,:,0:nrturb),ny_here,nz,nrturb+1,cor1,cor2,dx,dy,n1,n2,dir,.false.)
   call pseudo2d(vv(:,:,0:nrturb),ny_here,nz,nrturb+1,cor1,cor2,dx,dy,n1,n2,dir,.false.)
   call pseudo2d(ww(:,:,0:nrturb),ny_here,nz,nrturb+1,cor1,cor2,dx,dy,n1,n2,dir,.false.)
   call pseudo2d(rr(:,:,0:nrturb),ny_here,nz,nrturb+1,cor1,cor2,dx,dy,n1,n2,dir,.false.)

   do i=1,nrturb
      uu(:,:,i)=timecor*uu(:,:,i-1)+sqrt(1.0-timecor**2)*uu(:,:,i)
      vv(:,:,i)=timecor*vv(:,:,i-1)+sqrt(1.0-timecor**2)*vv(:,:,i)
      ww(:,:,i)=timecor*ww(:,:,i-1)+sqrt(1.0-timecor**2)*ww(:,:,i)
      rr(:,:,i)=timecor*rr(:,:,i-1)+sqrt(1.0-timecor**2)*rr(:,:,i)

      aveu=sum(uu(:,:,i))/real(ny_here*nz)
      avev=sum(vv(:,:,i))/real(ny_here*nz)
      avew=sum(ww(:,:,i))/real(ny_here*nz)
      aver=sum(rr(:,:,i))/real(ny_here*nz)

      uu(:,:,i)=uu(:,:,i)-aveu
      vv(:,:,i)=vv(:,:,i)-avev
      ww(:,:,i)=ww(:,:,i)-avew
      rr(:,:,i)=rr(:,:,i)-aver

      varu=sqrt(sum(uu(:,:,i)**2)/real(ny_here*nz-1))
      varv=sqrt(sum(vv(:,:,i)**2)/real(ny_here*nz-1))
      varw=sqrt(sum(ww(:,:,i)**2)/real(ny_here*nz-1))
      varr=sqrt(sum(rr(:,:,i)**2)/real(ny_here*nz-1))

      uu(:,:,i)=uu(:,:,i)/varu
      vv(:,:,i)=vv(:,:,i)/varv
      ww(:,:,i)=ww(:,:,i)/varw
      rr(:,:,i)=rr(:,:,i)/varr
   enddo
!   call tecfld('turb',ny_here,nz,nrturb,uu)
end subroutine
end module

