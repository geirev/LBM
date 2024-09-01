module m_nrelreadfoil
contains
subroutine nrelreadfoil()
   use mod_nrel5mw
   implicit none


   integer i,j,k


   character(len=14) :: filenames(1:nrc)= [ 'Cylinder1.dat ', &
                                            'Cylinder2.dat ', &
                                            'DU40_A17.dat  ', &
                                            'DU35_A17.dat  ', &
                                            'DU30_A17.dat  ', &
                                            'DU25_A17.dat  ', &
                                            'DU21_A17.dat  ', &
                                            'NACA64_A17.dat' ]

   foil(:,:)%deg=0.0
   foil(:,:)%cl=0.0
   foil(:,:)%cd=0.0
   foil(:,:)%cm=0.0

   do k=1,nrc
      print *,'filename: ',filenames(k)
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

   do k=1,nrc
      print *,'filename: ',filenames(k),nrang(k)
      do i=1,nrang(k)
         write(*,'(i4,f12.3,3f8.4)')i,foil(i,k)
      enddo
      print *
   enddo
end subroutine
end module
