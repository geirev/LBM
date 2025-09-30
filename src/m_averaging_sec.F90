module m_averaging_sec
   real,    dimension(:,:,:), allocatable, private, save :: uave, vave, wave
   real,    dimension(:,:,:), allocatable, private, save :: uave2, vave2, wave2
   real,    dimension(:,:,:), allocatable, private, save :: Ti

   real,    dimension(:,:)  , allocatable, private, save :: uxave, vxave, wxave
   real,    dimension(:,:)  , allocatable, private, save :: uxave2, vxave2, wxave2
   real,    dimension(:,:)  , allocatable, private, save :: Tix

   integer, allocatable, private, save                   :: iseci(:)

#ifdef _CUDA
   attributes(device) uave, vave, wave
   attributes(device) uave2, vave2, wave2
   attributes(device) Ti
   attributes(device) uxave, vxave, wxave
   attributes(device) uxave2, vxave2, wxave2
   attributes(device) Tix
   attributes(device) iseci
#endif

contains
subroutine averaging_sec(u,v,w,lfinal,iradius)
   use mod_dimensions
   use m_readinfile, only : ipos,jpos,kpos,uini,p2l
#ifdef _CUDA
   use m_readinfile, only : ntx,nty,ntz
#endif
   use mod_nrel5mw,  only : rotorradius,hubradius
   use m_tecfldsec
   use m_averaging_sec_kernel1
   use m_averaging_sec_kernel2
   use m_averaging_sec_kernelfin1
   use m_averaging_sec_kernelfin2
#ifdef _CUDA
   use cudafor
#endif
   implicit none
   real, intent(in)    :: u(nx,ny,nz)        ! x component of fluid velocity
   real, intent(in)    :: v(nx,ny,nz)        ! y component of fluid velocity
   real, intent(in)    :: w(nx,ny,nz)        ! z component of fluid velocity
#ifdef _CUDA
   attributes(device) :: u
   attributes(device) :: v
   attributes(device) :: w
#endif
   logical, intent(in) :: lfinal
   integer, intent(in) :: iradius

   integer, save :: ifirst=1
   integer, save :: iave=0
   real,    save :: D=0.0


   real, allocatable, dimension(:,:)   :: uxave_h, vxave_h, wxave_h
   real, allocatable, dimension(:,:)   :: Tix_h
   real, allocatable, dimension(:,:)   :: ave_h

   integer, save              :: nrsec
   integer, save              :: ja,jb,ka,kb,jdim,kdim

   character(len=2) csec
   real x
   integer j,k,i,isec,itmp
   real utmp,vtmp,wtmp,Ttmp
#ifdef _CUDA
   integer :: tx, ty, tz, bx, by, bz
#endif

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
         itmp=iseci(isec)
         print *,'i=',isec,itmp
      enddo

      ja=max(1  ,jpos(1)-2*iradius)
      jb=min(ny ,jpos(1)+2*iradius)
      ka=max(1  ,kpos(1)-2*iradius)
      kb=min(nz ,kpos(1)+2*iradius)
      ja=1
      jb=ny
      ka=1
      kb=nz
      print '(4(a,i3.3))','averaging ja:jb ka:kb ',ja,':',jb,' x ',ka,':',kb


! cross wake sections
      allocate(  uave(0:nrsec-1,ja:jb,ka:kb) )
      allocate(  vave(0:nrsec-1,ja:jb,ka:kb) )
      allocate(  wave(0:nrsec-1,ja:jb,ka:kb) )
      allocate( uave2(0:nrsec-1,ja:jb,ka:kb) )
      allocate( vave2(0:nrsec-1,ja:jb,ka:kb) )
      allocate( wave2(0:nrsec-1,ja:jb,ka:kb) )
      allocate(    Ti(0:nrsec-1,ja:jb,ka:kb) )

! along wake section
      allocate(  uxave(1:nx,1:ny) )
      allocate(  vxave(1:nx,1:ny) )
      allocate(  wxave(1:nx,1:ny) )
      allocate( uxave2(1:nx,1:ny) )
      allocate( vxave2(1:nx,1:ny) )
      allocate( wxave2(1:nx,1:ny) )
      allocate(    Tix(1:nx,1:ny) )

      uave  =0.0; vave  =0.0; wave  =0.0; uave2 =0.0; vave2 =0.0; wave2 =0.0; Ti=0.0
      uxave =0.0; vxave =0.0; wxave =0.0; uxave2=0.0; vxave2=0.0; wxave2=0.0; Tix=0.0
      ifirst=0
   endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   iave=iave+1
   jdim=jb-ja+1
   kdim=kb-ka+1
   !print '(a,8i5)','computing j-k sections averages',ja, jb, ka, kb, ny, jb-ja+1, nz, kb-ka+1,iave
#ifdef _CUDA
   tx=ntx; bx=(nrsec+tx-1)/tx
   ty=nty; by=(jdim +ty-1)/ty
   tz=ntz; bz=(kdim +tz-1)/tz
#endif
   call averaging_sec_kernel1&
#ifdef _CUDA
          &<<<dim3(bx,by,bz), dim3(tx,ty,tz)>>>&
#endif
          &(jdim, kdim, nrsec, ja, ka, iseci, u, v, w, uave, vave, wave, uave2, vave2, wave2 )

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   print *,'computing i-j  section averages'
#ifdef _CUDA
   tx=ntx; bx=(nx+tx-1)/tx
   ty=nty; by=(ny+ty-1)/ty
   tz=1; bz=1
#endif
   call averaging_sec_kernel2&
#ifdef _CUDA
          &<<<dim3(bx,by,bz), dim3(tx,ty,tz)>>>&
#endif
          &(u, v, w, uxave, vxave, wxave, uxave2, vxave2, wxave2, kpos(1))
!   print *,'done computing along-flow section averages'



   if (lfinal) then
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#ifdef _CUDA
   tx=ntx; bx=(nrsec  +tx-1)/tx
   ty=nty; by=(jdim+ty-1)/ty
   tz=ntz; bz=(kdim+tz-1)/tz
#endif
   call averaging_sec_kernelfin1&
#ifdef _CUDA
          &<<<dim3(bx,by,bz), dim3(tx,ty,tz)>>>&
#endif
          &(jdim, kdim, nrsec, iave,  uave, vave, wave, uave2, vave2, wave2, Ti, uini )


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#ifdef _CUDA
   tx=ntx; bx=(nx+tx-1)/tx
   ty=nty; by=(ny+ty-1)/ty
   tz=1; bz=1
#endif
   call averaging_sec_kernelfin2&
#ifdef _CUDA
          &<<<dim3(bx,by,bz), dim3(tx,ty,tz)>>>&
#endif
          &(iave, uxave, vxave, wxave, uxave2, vxave2, wxave2, Tix, uini )

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      print *,'kpos=',kpos(1),jpos(1)
      open(10,file='aveJ.dat')
         write(10,'(a)')'VARIABLES = "j" "x" "u" "v" "Ti"'
         do isec=0,nrsec-1
            write(csec,'(i2.2)')isec
            write(10,'(3a,i4,a)')'ZONE  T="D',csec,'" F=Point, I=  ',jdim,', J=   1, K=1'
            do j=ja,jb
               x=real(j-jpos(1))/D
               utmp=uave(isec,j,kpos(1))
               wtmp=wave(isec,j,kpos(1))
               Ttmp=  Ti(isec,j,kpos(1))
               write(10,'(i3,100g13.5)')j,x,utmp,wtmp,Ttmp
            enddo
         enddo
      close(10)

      open(10,file='aveK.dat')
         write(10,'(a)')'VARIABLES = "j" "x" "u" "v" "Ti"'
         do isec=0,nrsec-1
            write(csec,'(i2.2)')isec
            write(10,'(3a,i4,a)')'ZONE  T="D',csec,'" F=Point, I=  ',kdim,', J=   1, K=1'
            do k=ka,kb
               x=real(k-kpos(1))/D
               utmp=uave(isec,jpos(1),k)
               vtmp=vave(isec,jpos(1),k)
               Ttmp=  Ti(isec,jpos(1),k)
               write(10,'(i3,100g13.5)')k,x,utmp,vtmp,Ttmp
            enddo
         enddo
      close(10)

      allocate(ave_h(ja:jb,ka:kb))
      open(10,file='aveJK.dat',status='unknown')
         write(10,*)'TITLE = "aveJK"'
         write(10,*)'VARIABLES = "j" "k" "u" "v" "w" "Ti"'
         do isec=0,nrsec-1
            write(csec,'(i2.2)')isec
            if (isec==0) then
               write(10,'(3a,i5,a,i5,a)')' ZONE T="D',csec,'" F=BLOCK, I=',jdim,', J=',kdim,', K=1'
               write(10,'(30I5)')((j,j=ja,jb),k=ka,kb)
               write(10,'(30I5)')((k,j=ja,jb),k=ka,kb)
            else
               write(10,'(3a,i5,a,i5,a)')' ZONE T="D',csec,'" F=BLOCK, I=',jdim,', J=',kdim,', K=1 VARSHARELIST=([1,2]=1)'
            endif
            ave_h(:,:)=uave(isec,:,:); write(10,900)((ave_h(j,k),j=ja,jb),k=ka,kb)
            ave_h(:,:)=vave(isec,:,:); write(10,900)((ave_h(j,k),j=ja,jb),k=ka,kb)
            ave_h(:,:)=wave(isec,:,:); write(10,900)((ave_h(j,k),j=ja,jb),k=ka,kb)
            ave_h(:,:)=Ti(isec,:,:); write(10,900)((ave_h(j,k),j=ja,jb),k=ka,kb)
         enddo
      close(10)
      deallocate(ave_h)

      allocate(uxave_h(nx,ny))
      allocate(vxave_h(nx,ny))
      allocate(wxave_h(nx,ny))
      allocate(Tix_h(nx,ny))

      uxave_h=uxave
      vxave_h=vxave
      wxave_h=wxave
      Tix_h=Tix
      open(10,file='aveIJ.dat')
         write(10,*)'TITLE = "aveIJ.dat"'
         write(10,'(a)')'VARIABLES = "i" "j" "x" "y" "u" "v" "w" "Ti"'
         write(10,'(a,2(i4,a))')'ZONE  T="Averages in i-direction" F=BLOCK, I=  ',nx,', J=  ',ny,', K=1'
         write(10,'(30I5)')((i,i=1,nx),j=1,ny)
         write(10,'(30I5)')((j,i=1,nx),j=1,ny)
         write(10,900)((real(i-ipos(1))/D,i=1,nx),j=1,ny)
         write(10,900)((real(j-jpos(1))/D,i=1,nx),j=1,ny)
         write(10,900)((uxave_h(i,j),i=1,nx),j=1,ny)
         write(10,900)((vxave_h(i,j),i=1,nx),j=1,ny)
         write(10,900)((wxave_h(i,j),i=1,nx),j=1,ny)
         write(10,900)((Tix_h(i,j),i=1,nx),j=1,ny)
      close(10)
 900 format(10(1x,e12.5))


      deallocate( uave , vave , wave , uave2 , vave2 , wave2 , Ti)
      deallocate( uxave , vxave , wxave , uxave2 , vxave2 , wxave2 , Tix)
      deallocate( iseci )
      deallocate(uxave_h, vxave_h, wxave_h)
      deallocate(Tix_h)

      print '(a)','Done with averaging_sec.'
   endif

end subroutine

end module
