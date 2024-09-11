module m_averaging
contains
subroutine averaging(u,v,w,lfinal,iradius)
   use mod_dimensions
   use m_readinfile, only : ipos,jpos,kpos,uini,p2l
   use mod_nrel5mw, only : rotorradius,hubradius
   use m_tecfld
   real, intent(in)    :: u(nx,ny,nz)        ! x component of fluid velocity
   real, intent(in)    :: v(nx,ny,nz)        ! y component of fluid velocity
   real, intent(in)    :: w(nx,ny,nz)        ! z component of fluid velocity
   logical, intent(in) :: lfinal
   integer, intent(in) :: iradius

   integer, save :: ifirst=1
   integer, save :: iave=0
   real, save :: D=0.0

   real, save, allocatable :: uave(:,:,:), vave(:,:,:), wave(:,:,:)
   real, save, allocatable :: uave2(:,:,:), vave2(:,:,:), wave2(:,:,:)
   real, save, allocatable :: Ti(:,:,:)

   real, save, allocatable :: uxave(:), vxave(:), wxave(:)
   real, save, allocatable :: uxave2(:), vxave2(:), wxave2(:)
   real, save, allocatable :: Tix(:)

   integer, save, allocatable :: iseci(:)

   integer, save :: nrsec
   integer j,k,i,isec
   integer, save :: ja,jb,ka,kb
   integer, save :: jaa,jbb,kaa,kbb
   character(len=2) csec
   real x

   if (ifirst == 1) then
      D=2.0*(rotorradius+hubradius)
      print *,'rotor diameter',D*p2l%length
      nrsec=int((nx-ipos(1))/(2*iradius))+1
      nrsec=min(nrsec,nx-1)
      print *,'Number of diagnostic cross sections=',nrsec
      ja=max(1,jpos(1)-2*iradius)
      jb=min(ny,jpos(1)+2*iradius)
      ka=max(1,kpos(1)-2*iradius)
      kb=min(nz,kpos(1)+2*iradius)
      jaa=max(1,jpos(1)-10)
      jbb=min(ny,jpos(1)+10)
      kaa=max(1,kpos(1)-10)
      kbb=min(nz,kpos(1)+10)
      allocate( uave(0:nrsec-1,ja:jb,ka:kb) )
      allocate( vave(0:nrsec-1,ja:jb,ka:kb) )
      allocate( wave(0:nrsec-1,ja:jb,ka:kb) )
      allocate( uave2(0:nrsec-1,ja:jb,ka:kb) )
      allocate( vave2(0:nrsec-1,ja:jb,ka:kb) )
      allocate( wave2(0:nrsec-1,ja:jb,ka:kb) )
      allocate( Ti(0:nrsec-1,ja:jb,ka:kb) )

      allocate( uxave(1:nx),  vxave(1:nx),  wxave(1:nx) )
      allocate( uxave2(1:nx), vxave2(1:nx), wxave2(1:nx) )
      allocate( Tix(1:nx) )

      allocate( iseci(0:nrsec-1) )
      uave=0.0 ;  vave=0.0;   wave=0.0
      uave2=0.0;  vave2=0.0;  wave2=0.0
      uxave=0.0;  vxave=0.0;  wxave=0.0
      uxave2=0.0; vxave2=0.0; wxave2=0.0

      ifirst=0
      print *,'Diagnostic section i indices'
      isec=0
      do i=ipos(1),nx,2*iradius
         iseci(isec)=i
         print *,'i=',isec,iseci(isec)
         isec=isec+1
      enddo
   endif


   iave=iave+1
! section averages
   do isec=0,nrsec-1
      uave(isec,:,:)=uave(isec,:,:)+u(iseci(isec),:,:)
      vave(isec,:,:)=vave(isec,:,:)+v(iseci(isec),:,:)
      wave(isec,:,:)=wave(isec,:,:)+w(iseci(isec),:,:)

      uave2(isec,:,:)=uave2(isec,:,:)+u(iseci(isec),:,:)**2
      vave2(isec,:,:)=vave2(isec,:,:)+v(iseci(isec),:,:)**2
      wave2(isec,:,:)=wave2(isec,:,:)+w(iseci(isec),:,:)**2
   enddo

! Alongflow averages at center +- 10 grid points (r/D=0.625)
   uxave(:)=uxave(:) + u(:,jaa,kpos(1)) + u(:,jbb,kpos(1)) + u(:,jpos(1),kaa) + u(:,jpos(1),kbb)
   vxave(:)=vxave(:) + v(:,jaa,kpos(1)) + v(:,jbb,kpos(1)) + v(:,jpos(1),kaa) + v(:,jpos(1),kbb)
   wxave(:)=wxave(:) + w(:,jaa,kpos(1)) + w(:,jbb,kpos(1)) + w(:,jpos(1),kaa) + w(:,jpos(1),kbb)

   uxave2(:)=uxave2(:) + u(:,jaa,kpos(1))**2 + u(:,jbb,kpos(1))**2 + u(:,jpos(1),kaa)**2 + u(:,jpos(1),kbb)**2
   vxave2(:)=vxave2(:) + v(:,jaa,kpos(1))**2 + v(:,jbb,kpos(1))**2 + v(:,jpos(1),kaa)**2 + v(:,jpos(1),kbb)**2
   wxave2(:)=wxave2(:) + w(:,jaa,kpos(1))**2 + w(:,jbb,kpos(1))**2 + w(:,jpos(1),kaa)**2 + w(:,jpos(1),kbb)**2

   if (lfinal) then
      uave=uave/real(iave)
      vave=vave/real(iave)
      wave=wave/real(iave)

      uxave=uxave/real(4*iave)
      vxave=vxave/real(4*iave)
      wxave=wxave/real(4*iave)

      uave2=uave2/real(iave)
      vave2=vave2/real(iave)
      wave2=wave2/real(iave)

      uxave2=uxave2/real(4*iave)
      vxave2=vxave2/real(4*iave)
      wxave2=wxave2/real(4*iave)

      Ti=uave2-uave**2 + vave2-vave**2 + wave2-wave**2
      Ti=sqrt(Ti/3.0)/uini

      Tix=uxave2-uxave**2 + vxave2-vxave**2 + wxave2-wxave**2
      Tix=sqrt(Tix/3.0)/uini

      uave=uave/uini
      vave=vave/uini
      wave=wave/uini

      open(10,file='aveJ.dat')
         call writetecheader(nrsec,ja,jb)
         do j=ja,jb
            x=real(j-jpos(1))/D
            write(10,'(i3,100g13.5)')j,x,uave(:,j,kpos(1)),wave(:,j,kpos(1)),Ti(:,j,kpos(1))
         enddo
      close(10)

      open(10,file='aveK.dat')
         call writetecheader(nrsec,ka,kb)
         do k=ka,kb
            x=real(k-kpos(1))/D
            write(10,'(i3,100g13.5)')k,x,uave(:,jpos(1),k),vave(:,jpos(1),k),Ti(:,jpos(1),k)
         enddo
      close(10)

      open(10,file='aveI.dat')
         write(10,'(a)')'VARIABLES = "i" "x" "uave" "vave" "wave" "Ti"'
         write(10,'(a,i4,a)')'ZONE  T="Averages in i-direction" F=Point, I=  ',nx,', J=   1, K=1'
         do i=1,nx
            x=real(i-ipos(1))/D
            write(10,'(i3,5g13.5)')i,x,uxave(i),vxave(i),wxave(i),Tix(i)
         enddo
      close(10)

      do i=0,nrsec-1
         write(csec,'(i2.2)')i
         call tecfld('averages2dim'//csec//'D.dat',ja,jb,ka,kb,uave(i,:,:),vave(i,:,:),wave(i,:,:),Ti(i,:,:))
      enddo

      deallocate( uave , vave , wave , uave2 , vave2 , wave2 , Ti)
      deallocate( uxave , vxave , wxave , uxave2 , vxave2 , wxave2 , Tix)
      deallocate( iseci )
   endif

end subroutine

subroutine writetecheader(nrsec,ia,ib)
   integer, intent(in) :: nrsec,ia,ib
   character(len=2) csec
   integer i
   write(10,'(a)',advance='no')'VARIABLES = "i" "x" '
   do i=0,nrsec-1
      write(csec,'(i2.2)')i
      write(10,'(3a)',advance='no')' "u D',csec,'"'
   enddo
   do i=0,nrsec-1
      write(csec,'(i2.2)')i
      write(10,'(3a)',advance='no')' "v D',csec,'"'
   enddo
   do i=0,nrsec-1
      write(csec,'(i2.2)')i
      write(10,'(3a)',advance='no')' "Ti D',csec,'"'
   enddo
   write(10,*)
   write(10,'(a,i4,a)')'ZONE  T="Averages in j-direction" F=Point, I=  ',ib-ia+1,', J=   1, K=1'
end subroutine
end module
