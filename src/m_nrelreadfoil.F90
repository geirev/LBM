module m_nrelreadfoil
contains
subroutine nrelreadfoil()
   use mod_nrel5mw
   implicit none


   integer i,k


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

   open(11,file='tables.dat')
      write(11,'(a)')'VARIABLES = "i" "angle" "CL" "CD" "CM"'
      do k=1,nrc
         write(11,'(3a,i4,a)')'ZONE  T="',trim(filenames(k)),'" F=Point, I=  ',nrang(k),' J=   1, K=1'
         do i=1,nrang(k)
            write(11,'(i4,50f10.4)')i,foil(i,k)
         enddo
         write(11,*)
      enddo
   close(11)
end subroutine
end module
