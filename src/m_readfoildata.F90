module m_readfoildata

contains
subroutine readfoildata(cl,cd)
   use mod_nrl5mw
   implicit none

   type foildata
      real deg
      real cl
      real cd
      real cm
   end type

   real, intent(inout) :: cl(nrchords)
   real, intent(inout) :: cd(nrchords)

   integer, parameter :: nrc=8
   integer nrang(nrc),i,j,k
   real :: pitch=0.0

   real slope,angle

   type(foildata), allocatable :: foil(:,:)

   character(len=14) :: filenames(1:nrc)= [ 'Cylinder1.dat ', &
                                            'Cylinder2.dat ', &
                                            'DU40_A17.dat  ', &
                                            'DU35_A17.dat  ', &
                                            'DU30_A17.dat  ', &
                                            'DU25_A17.dat  ', &
                                            'DU21_A17.dat  ', &
                                            'NACA64_A17.dat' ]

   allocate(foil(150,8))
   foil(:,:)%deg=0.0
   foil(:,:)%cl=0.0
   foil(:,:)%cd=0.0
   foil(:,:)%cm=0.0

   do k=1,nrc
!      print *,'filename: ',filenames(k)
      open(10,file='./Airfoils/'//trim(filenames(k)))
      do i=1,50
         read(10,*)
      enddo
      read(10,*)nrang(k)
      read(10,*)
      read(10,*)
      do i=1,nrang(k)
         read(10,*)foil(i,k)
      enddo
      close(10)
   enddo

!   do k=1,nrc
!      print *,'filename: ',filenames(k),nrang(k)
!      do i=1,nrang(k)
!         write(*,'(i4,f12.3,3f8.4)')i,foil(i,k)
!      enddo
!      print *
!   enddo

   do j=1,nrchords
      angle=twist(j) + pitch
      k=nfoil(j)
      do i=1,nrang(k)
         if ( foil(i,k)%deg < angle .and. angle <= foil(i+1,k)%deg ) then
            slope=(foil(i+1,k)%cl-foil(i,k)%cl)/(foil(i+1,k)%deg-foil(i,k)%deg)
            cl(j)=foil(i,k)%cl + slope*(angle - foil(i,k)%deg)

            slope=(foil(i+1,k)%cd-foil(i,k)%cd)/(foil(i+1,k)%deg-foil(i,k)%deg)
            cd(j)=foil(i,k)%cd + slope*(angle - foil(i,k)%deg)

!            print '(2i3,3(a,3f10.4))',j,i , 'angle: ', foil(i,k)%deg , angle, foil(i+1,k)%deg, &
!                          ' -- ',  foil(i,k)%cl, cl(j), foil(i+1,k)%cl,&
!                          ' -- ',  foil(i,k)%cd, cd(j), foil(i+1,k)%cd
            exit
         endif
      enddo

   enddo

   deallocate(foil)

end subroutine
end module
