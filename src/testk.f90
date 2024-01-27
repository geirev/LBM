program testk
   integer, parameter :: nx = 10
   integer, parameter :: ny = 10
   integer, parameter :: nz = 10

   integer i,j,k,ia,ib,ja,jb,ka,kb

   do i=1,nx
      ia=mod(nx+i-1-1,nx)+1
      ib=mod(nx+i-1+1,nx)+1
      print *,'I',ia,i,ib
   enddo

end program
