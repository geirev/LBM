module m_averaging
contains
subroutine averaging(u,v,w,lfinal,iradius)
   use mod_dimensions
   use m_readinfile, only : ipos,jpos,kpos,uini,p2l
   use mod_nrel5mw, only : rotorradius,hubradius
   use m_tecfldsec
   real, intent(in)    :: u(nx,ny,nz)        ! x component of fluid velocity
   real, intent(in)    :: v(nx,ny,nz)        ! y component of fluid velocity
   real, intent(in)    :: w(nx,ny,nz)        ! z component of fluid velocity
#ifdef _CUDA
   attributes(managed) :: u
   attributes(managed) :: v
   attributes(managed) :: w
#endif
   logical, intent(in) :: lfinal
   integer, intent(in) :: iradius

   integer, save :: ifirst=1
   integer, save :: iave=0
   real,    save :: D=0.0

   real,    save, allocatable :: uave(:,:,:), vave(:,:,:), wave(:,:,:)
   real,    save, allocatable :: uave2(:,:,:), vave2(:,:,:), wave2(:,:,:)
   real,    save, allocatable :: Ti(:,:,:)

   real,    save, allocatable :: uxave(:), vxave(:), wxave(:)
   real,    save, allocatable :: uxave2(:), vxave2(:), wxave2(:)
   real,    save, allocatable :: Tix(:)

   integer, save, allocatable :: iseci(:)
   integer, save              :: nrsec
   integer, save              :: ja,jb,ka,kb
   integer, save              :: jaa,jbb,kaa,kbb

   character(len=2) csec
   real x
   integer j,k,i,isec

   if (ifirst == 1) then
      D=2.0*(rotorradius+hubradius)
      print *,'rotor diameter',D*p2l%length

      nrsec=int(real(nx-ipos(1))/real(2*iradius))+1
      if (ipos(1)+(nrsec-1)*2*iradius > nx-1) nrsec=nrsec-1
      print *,'Number of diagnostic cross sections=',nrsec
      if (nrsec < 1) stop 'averaging: nrsec less than one'

      print *,'Diagnostic section i indices'
      allocate( iseci(0:nrsec-1) )
      do isec=0,nrsec-1
         iseci(isec)=ipos(1)+isec*2*iradius
         print *,'i=',isec,iseci(isec)
      enddo

      ja=max(1  ,jpos(1)-2*iradius)
      jb=min(ny ,jpos(1)+2*iradius)
      ka=max(1  ,kpos(1)-2*iradius)
      kb=min(nz ,kpos(1)+2*iradius)
      print '(4(a,i3.3))','averaging ja:jb ka:kb ',ja,':',jb,' x ',ka,':',kb

      jaa=max(1 ,jpos(1)-10)
      jbb=min(ny,jpos(1)+10)
      kaa=max(1 ,kpos(1)-10)
      kbb=min(nz,kpos(1)+10)
      print '(4(a,i3.3))','averaging jaa=',jaa,', jbb=',jbb,', kaa=',kaa,', kbb=',kbb

      allocate(  uave(0:nrsec-1,ja:jb,ka:kb) )
      allocate(  vave(0:nrsec-1,ja:jb,ka:kb) )
      allocate(  wave(0:nrsec-1,ja:jb,ka:kb) )
      allocate( uave2(0:nrsec-1,ja:jb,ka:kb) )
      allocate( vave2(0:nrsec-1,ja:jb,ka:kb) )
      allocate( wave2(0:nrsec-1,ja:jb,ka:kb) )
      allocate(    Ti(0:nrsec-1,ja:jb,ka:kb) )

      allocate(  uxave(1:nx),  vxave(1:nx),  wxave(1:nx) )
      allocate( uxave2(1:nx), vxave2(1:nx), wxave2(1:nx) )
      allocate(    Tix(1:nx) )

      uave  =0.0; vave  =0.0; wave  =0.0
      uave2 =0.0; vave2 =0.0; wave2 =0.0
      uxave =0.0; vxave =0.0; wxave =0.0
      uxave2=0.0; vxave2=0.0; wxave2=0.0

      ifirst=0
   endif

   iave=iave+1
! section averages
   do k=ka,kb
   do j=ja,jb
   do isec=0,nrsec-1
      uave(isec,j,k)=uave(isec,j,k)+u(iseci(isec),j,k)
      vave(isec,j,k)=vave(isec,j,k)+v(iseci(isec),j,k)
      wave(isec,j,k)=wave(isec,j,k)+w(iseci(isec),j,k)

      uave2(isec,j,k)=uave2(isec,j,k)+u(iseci(isec),j,k)**2
      vave2(isec,j,k)=vave2(isec,j,k)+v(iseci(isec),j,k)**2
      wave2(isec,j,k)=wave2(isec,j,k)+w(iseci(isec),j,k)**2
   enddo
   enddo
   enddo

! Alongflow averages at center +- 10 grid points (r/D=0.625)
   do i=1,nx
      uxave(i)=uxave(i) + u(i,jaa,kpos(1)) + u(i,jbb,kpos(1)) + u(i,jpos(1),kaa) + u(i,jpos(1),kbb)
      vxave(i)=vxave(i) + v(i,jaa,kpos(1)) + v(i,jbb,kpos(1)) + v(i,jpos(1),kaa) + v(i,jpos(1),kbb)
      wxave(i)=wxave(i) + w(i,jaa,kpos(1)) + w(i,jbb,kpos(1)) + w(i,jpos(1),kaa) + w(i,jpos(1),kbb)

      uxave2(i)=uxave2(i) + u(i,jaa,kpos(1))**2 + u(i,jbb,kpos(1))**2 + u(i,jpos(1),kaa)**2 + u(i,jpos(1),kbb)**2
      vxave2(i)=vxave2(i) + v(i,jaa,kpos(1))**2 + v(i,jbb,kpos(1))**2 + v(i,jpos(1),kaa)**2 + v(i,jpos(1),kbb)**2
      wxave2(i)=wxave2(i) + w(i,jaa,kpos(1))**2 + w(i,jbb,kpos(1))**2 + w(i,jpos(1),kaa)**2 + w(i,jpos(1),kbb)**2
   enddo

   if (lfinal) then
      uave=uave/real(iave)
      vave=vave/real(iave)
      wave=wave/real(iave)

      uave2=uave2/real(iave)
      vave2=vave2/real(iave)
      wave2=wave2/real(iave)

      uxave=uxave/real(4*iave)
      vxave=vxave/real(4*iave)
      wxave=wxave/real(4*iave)

      uxave2=uxave2/real(4*iave)
      vxave2=vxave2/real(4*iave)
      wxave2=wxave2/real(4*iave)

      do k=ka,kb
      do j=ja,jb
      do isec=0,nrsec-1
         Ti(isec,j,k)=uave2(isec,j,k)-uave(isec,j,k)**2 + vave2(isec,j,k)-vave(isec,j,k)**2 + wave2(isec,j,k)-wave(isec,j,k)**2
         Ti(isec,j,k)=sqrt(Ti(isec,j,k)/3.0)/uini
      enddo
      enddo
      enddo

      do i=1,nx
         Tix(i)=uxave2(i)-uxave(i)**2 + vxave2(i)-vxave(i)**2 + wxave2(i)-wxave(i)**2
         Tix(i)=sqrt(Tix(i)/3.0)/uini
      enddo

      uave=uave/uini
      vave=vave/uini
      wave=wave/uini

      open(10,file='aveJ.dat')
         write(10,'(a)',advance='no')'VARIABLES = "j" "x" "u" "w" "Ti"'
         do isec=0,nrsec-1
            write(csec,'(i2.2)')isec
            write(10,'(3a,i4,a)')'ZONE  T="D',csec,'" F=Point, I=  ',jb-ja+1,', J=   1, K=1'
            do j=ja,jb
               x=real(j-jpos(1))/D
               write(10,'(i3,100g13.5)')j,x,uave(isec,j,kpos(1)),wave(isec,j,kpos(1)),Ti(isec,j,kpos(1))
            enddo
         enddo
      close(10)

      open(10,file='aveK.dat')
         write(10,'(a)',advance='no')'VARIABLES = "k" "x" "u" "v" "Ti"'
         do isec=0,nrsec-1
            write(csec,'(i2.2)')isec
            write(10,'(3a,i4,a)')'ZONE  T="D',csec,'" F=Point, I=  ',kb-ka+1,', J=   1, K=1'
            do k=ka,kb
               x=real(k-kpos(1))/D
               write(10,'(i3,100g13.5)')k,x,uave(isec,jpos(1),k),vave(isec,jpos(1),k),Ti(isec,jpos(1),k)
            enddo
         enddo
      close(10)

      open(10,file='aveI.dat')
         write(10,'(a)')'VARIABLES = "i" "x" "u" "v" "w" "Ti"'
         write(10,'(a,i4,a)')'ZONE  T="Averages in i-direction" F=Point, I=  ',nx,', J=   1, K=1'
         do i=1,nx
            x=real(i-ipos(1))/D
            write(10,'(i3,5g13.5)')i,x,uxave(i),vxave(i),wxave(i),Tix(i)
         enddo
      close(10)

      call tecfldsec('averages2dim.dat',nrsec,ja,jb,ka,kb,uave,vave,wave,Ti)

      deallocate( uave , vave , wave , uave2 , vave2 , wave2 , Ti)
      deallocate( uxave , vxave , wxave , uxave2 , vxave2 , wxave2 , Tix)
      deallocate( iseci )
      print '(a)','Done with averaging.'
   endif

end subroutine

!subroutine writetecheader(nrsec,ia,ib)
!   integer, intent(in) :: nrsec,ia,ib
!   character(len=2) csec
!   integer i
!   print *,'techeader'
!   write(10,'(a)',advance='no')'VARIABLES = "i" "x" '
!   do i=0,nrsec-1
!      write(csec,'(i2.2)')i
!      write(10,'(3a)',advance='no')' "u D',csec,'"'
!   enddo
!   do i=0,nrsec-1
!      write(csec,'(i2.2)')i
!      write(10,'(3a)',advance='no')' "v D',csec,'"'
!   enddo
!   do i=0,nrsec-1
!      write(csec,'(i2.2)')i
!      write(10,'(3a)',advance='no')' "Ti D',csec,'"'
!   enddo
!   write(10,*)
!   write(10,'(a,i4,a)')'ZONE  T="Averages in j-direction" F=Point, I=  ',ib-ia+1,', J=   1, K=1'
!end subroutine
end module
