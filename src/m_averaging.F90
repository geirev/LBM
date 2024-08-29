module m_averaging
contains
subroutine averaging(u,v,w,lfinal,iradius)
   use mod_dimensions
   use m_readinfile, only : ipos,jpos,kpos,uini,turbrad,p2l
   use mod_nrl5mw, only : rotorradius,hubradius
   use m_tecfld
   real, intent(in)    :: u(nx,ny,nz)        ! x component of fluid velocity
   real, intent(in)    :: v(nx,ny,nz)        ! y component of fluid velocity
   real, intent(in)    :: w(nx,ny,nz)        ! z component of fluid velocity
   logical, intent(in) :: lfinal
   integer, intent(in) :: iradius

   integer, save :: ifirst=1
   integer, save :: iave=0
   real, save :: D=0.0

   real, save, allocatable :: uave(:,:,:)
   real, save, allocatable :: vave(:,:,:)
   real, save, allocatable :: wave(:,:,:)
   integer, save, allocatable :: iseci(:)

   integer, save :: nrsec
   integer j,k,i,isec
   integer, save :: ja,jb,ka,kb
   character(len=2) csec
   real x

   if (ifirst == 1) then
      D=2.0*(rotorradius+hubradius)
      print *,'rotor diameter',D*p2l%length
      nrsec=int((nx-ipos(1))/(2*iradius))+1
      print *,'Number of diagnostic cross sections=',nrsec
      ja=max(1,jpos(1)-2*iradius)
      jb=min(ny,jpos(1)+2*iradius)
      ka=max(1,kpos(1)-2*iradius)
      kb=min(nz,kpos(1)+2*iradius)
      allocate( uave(0:nrsec-1,ja:jb,ka:kb) )
      allocate( vave(0:nrsec-1,ja:jb,ka:kb) )
      allocate( wave(0:nrsec-1,ja:jb,ka:kb) )
      allocate( iseci(0:nrsec-1) )
      uave=0.0
      vave=0.0
      wave=0.0
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
   do isec=0,nrsec-1
      uave(isec,ja:jb,ka:kb)=uave(isec,ja:jb,ka:kb)+u(iseci(isec),ja:jb,ka:kb)
      vave(isec,ja:jb,ka:kb)=vave(isec,ja:jb,ka:kb)+v(iseci(isec),ja:jb,ka:kb)
      wave(isec,ja:jb,ka:kb)=wave(isec,ja:jb,ka:kb)+w(iseci(isec),ja:jb,ka:kb)
   enddo

   if (lfinal) then
      uave=uave/real(iave)
      vave=vave/real(iave)
      wave=wave/real(iave)

      uave=uave/uini
      vave=vave/uini
      wave=wave/uini

      open(10,file='aveJ.dat')
         do j=ja,jb
            x=real(j-jpos(1))/D
            write(10,'(i3,60f10.3)')j,x,uave(:,j,kpos(1)),wave(:,j,kpos(1))
         enddo
      close(10)

      open(10,file='aveK.dat')
         do k=ka,kb
            x=real(k-kpos(1))/D
            write(10,'(i3,60f10.3)')k,x,uave(:,jpos(1),k),vave(:,jpos(1),k)
         enddo
      close(10)

      do i=0,nrsec-1
         write(csec,'(i2.2)')i
         call tecfld('averages2dim'//csec//'D.dat',ja,jb,ka,kb,uave(i,:,:),vave(i,:,:),wave(i,:,:))
      enddo

   endif


end subroutine
end module
