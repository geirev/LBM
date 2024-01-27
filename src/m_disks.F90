module m_disks
contains
subroutine disks(blanking)
   use mod_dimensions
   logical, intent(inout) :: blanking(nx,ny,nz)
   integer i,j

   open(10,file='disk.dat')
   do i=100,400,50
   do j=20,25
      blanking(i,j,:)=.true.
      write(10,'(2i3)')i,j
   enddo
   do j=50,55
      blanking(i,j,:)=.true.
      write(10,'(2i3)')i,j
   enddo
   do j=75,80
      blanking(i,j,:)=.true.
      write(10,'(2i3)')i,j
   enddo
   enddo

   close(10)

   open(10,file='blanking.dat')
      do j=1,ny
      do i=1,nx
         if (blanking(i,j,1)) write(10,'(2i4,a)')i,j,' 1.0'
      enddo
      enddo
   close(10)
end subroutine
end module

