module m_inflow_turbulence_compute
contains
subroutine inflow_turbulence_compute(uu,vv,ww,rr,ny_here,nz,nrturb,lfirst)
   use m_pseudo2D
   use m_tecfld
   use m_readinfile, only : p2l,timecor, turb_length
   implicit none
   integer, intent(in) :: ny_here,nz,nrturb
   logical, intent(in) :: lfirst
   real, intent(inout) :: uu(ny_here,nz,0:nrturb)
   real, intent(inout) :: vv(ny_here,nz,0:nrturb)
   real, intent(inout) :: ww(ny_here,nz,0:nrturb)
   real, intent(inout) :: rr(ny_here,nz,0:nrturb)
   real :: cor1, cor2
   real :: dx, dy, dir=0.0
   integer :: i,j,k
   real :: aveu,avev,avew,aver,varu,varv,varw,varr
   integer n1,n2
   integer n0,nn


   dx=p2l%length
   dy=p2l%length

   cor1=turb_length/sqrt(3.0)
   cor2=turb_length/sqrt(3.0)

   if (cor1 < 3.0*dx .or. cor1 > 100.0*dx) then
      print *,'WARNING: turbulence correlation length poorly scaled'
      print *,'  cor1 = ',cor1
      print *,'  dx   = ',dx
      print *,'  cor1/dx = ',cor1/dx
   endif

   print *,'compute_turbulence_field: generating pseudo-2D inflow forcing'


   if (lfirst) then
      n0=0
      nn=1
   else
      n0=1
      nn=0
      uu(:,:,0)=uu(:,:,nrturb)
      vv(:,:,0)=vv(:,:,nrturb)
      ww(:,:,0)=ww(:,:,nrturb)
      rr(:,:,0)=rr(:,:,nrturb)
   endif

   n1=ny_here
   n2=nz
   call pseudo2d(uu(:,:,n0:nrturb),ny_here,nz,nrturb+nn,cor1,cor2,dx,dy,n1,n2,dir,.false.)

   n1=ny_here; n2=nz
   call pseudo2d(vv(:,:,n0:nrturb),ny_here,nz,nrturb+nn,cor1,cor2,dx,dy,n1,n2,dir,.false.)

   n1=ny_here; n2=nz
   call pseudo2d(ww(:,:,n0:nrturb),ny_here,nz,nrturb+nn,cor1,cor2,dx,dy,n1,n2,dir,.false.)

   n1=ny_here; n2=nz
   call pseudo2d(rr(:,:,n0:nrturb),ny_here,nz,nrturb+nn,cor1,cor2,dx,dy,n1,n2,dir,.false.)

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
!   call tecfld('tec_turb',ny_here,nz,nrturb,uu)
end subroutine
end module

