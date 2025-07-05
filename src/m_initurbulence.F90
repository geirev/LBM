module m_initurbulence
contains
subroutine initurbulence(uu,vv,ww,rr,lfirst)
   use mod_dimensions
   use m_set_random_seed2
   use m_pseudo2D
   implicit none
   real, intent(inout)  :: uu(ny,nz,0:nrturb)
   real, intent(inout)  :: vv(ny,nz,0:nrturb)
   real, intent(inout)  :: ww(ny,nz,0:nrturb)
   real, intent(inout)  :: rr(ny,nz,0:nrturb)
   logical, intent(in)  :: lfirst

   real cor1,cor2,dx,dy,dir
   integer(kind=4) n1,n2
   integer i,j,k
   integer, save :: nyy,nzz,nt
   logical :: verbose=.false.
   real aveu,avev,avew,aver
   real varu,varv,varw,varr

   real :: timecor=0.98

! New seed
   call system('rm seed.dat')
   call set_random_seed2

! Simulating a time series of inflow boundary perturbations for u
   print '(a)','initurbulence: Simulating inflow turbulence forcing'
   if (lfirst) then
      print '(a,l1)','initurbulence: lfirst=',lfirst
      uu(:,:,0)=0.0
      vv(:,:,0)=0.0
      ww(:,:,0)=0.0
      rr(:,:,0)=0.0
   else
      print '(a,l1)','initurbulence: lfirst=',lfirst
      uu(:,:,0)=uu(:,:,nrturb)
      vv(:,:,0)=vv(:,:,nrturb)
      ww(:,:,0)=ww(:,:,nrturb)
      rr(:,:,0)=rr(:,:,nrturb)
   endif

   cor1=10.0/sqrt(3.0)
   cor2=10.0/sqrt(3.0)
   dir=0.0
   dx=1.0
   dy=1.0
   n1=int(ny,4)
   n2=int(nz,4)
   nyy=int(ny,4)
   nzz=int(nz,4)
   nt=int(nrturb)
   call pseudo2d(uu(:,:,1:nt),nyy,nzz,nt,cor1,cor2,dx,dy,n1,n2,dir,verbose)
   call pseudo2d(vv(:,:,1:nt),nyy,nzz,nt,cor1,cor2,dx,dy,n1,n2,dir,verbose)
   call pseudo2d(ww(:,:,1:nt),nyy,nzz,nt,cor1,cor2,dx,dy,n1,n2,dir,verbose)
   call pseudo2d(rr(:,:,1:nt),nyy,nzz,nt,cor1,cor2,dx,dy,n1,n2,dir,verbose)

! Imposing time correlations
   do i=1,nrturb
      uu(:,:,i)=timecor*uu(:,:,i-1)+sqrt(1.0-timecor**2)*uu(:,:,i)
      vv(:,:,i)=timecor*vv(:,:,i-1)+sqrt(1.0-timecor**2)*vv(:,:,i)
      ww(:,:,i)=timecor*ww(:,:,i-1)+sqrt(1.0-timecor**2)*ww(:,:,i)
      rr(:,:,i)=timecor*rr(:,:,i-1)+sqrt(1.0-timecor**2)*rr(:,:,i)
   enddo

! Ensure zero mean and variance equal to one
   print *,'removing nonzero average and scaling variance'
   do i=1,nrturb
      aveu=0.0
      avev=0.0
      avew=0.0
      aver=0.0
      do k=1,nz
      do j=1,ny
         aveu = aveu + uu(j,k,i)
         avev = avev + vv(j,k,i)
         avew = avew + ww(j,k,i)
         aver = aver + rr(j,k,i)
      enddo
      enddo
      aveu=aveu/real(ny*nz)
      avev=avev/real(ny*nz)
      avew=avew/real(ny*nz)
      aver=aver/real(ny*nz)
      uu(:,:,i)=uu(:,:,i)-aveu
      vv(:,:,i)=vv(:,:,i)-avev
      ww(:,:,i)=ww(:,:,i)-avew
      rr(:,:,i)=rr(:,:,i)-aver

      varu=0.0
      varv=0.0
      varw=0.0
      varr=0.0
      do k=1,nz
      do j=1,ny
         varu = varu + uu(j,k,i)**2
         varv = varv + vv(j,k,i)**2
         varw = varw + ww(j,k,i)**2
         varr = varr + rr(j,k,i)**2
      enddo
      enddo
      varu=sqrt(varu/real(ny*nz-1))
      varv=sqrt(varv/real(ny*nz-1))
      varw=sqrt(varw/real(ny*nz-1))
      varr=sqrt(varr/real(ny*nz-1))
      uu(:,:,i)=(1.0/varu)*uu(:,:,i)
      vv(:,:,i)=(1.0/varv)*vv(:,:,i)
      ww(:,:,i)=(1.0/varw)*ww(:,:,i)
      rr(:,:,i)=(1.0/varr)*rr(:,:,i)

   enddo

   print *,'checking average and variance'
   do i=1,nrturb
      aveu=0.0
      avev=0.0
      avew=0.0
      aver=0.0
      do k=1,nz
      do j=1,ny
         aveu = aveu + uu(j,k,i)
         avev = avev + vv(j,k,i)
         avew = avew + ww(j,k,i)
         aver = aver + rr(j,k,i)
      enddo
      enddo
      aveu=aveu/real(ny*nz)
      avev=avev/real(ny*nz)
      avew=avew/real(ny*nz)
      aver=aver/real(ny*nz)

      varu=0.0
      varv=0.0
      varw=0.0
      varr=0.0
      do k=1,nz
      do j=1,ny
         varu = varu + uu(j,k,i)**2
         varv = varv + vv(j,k,i)**2
         varw = varw + ww(j,k,i)**2
         varr = varr + rr(j,k,i)**2
      enddo
      enddo
      varu=varu/real(ny*nz-1)
      varv=varv/real(ny*nz-1)
      varw=varw/real(ny*nz-1)
      varr=varr/real(ny*nz-1)
  !    print '(a,i5,8g13.5)','check turbulence stat:',i,aveu,avev,avew,aver,varu,varv,varw,varr

   enddo

end subroutine
end module
