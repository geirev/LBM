program testk
   implicit none
   integer, parameter :: nx = 10
   integer, parameter :: cxs(1:3) = [0, 1,-1]

   integer i,j,k,ia,ib,ja,jb,ka,kb,i1,i2,shift
   real a(nx),b(nx)
   call random_number(a)

   do j=1,3
   do i=1,nx
      i1=mod(nx+i-1-cxs(j),nx)+1
      b(i)=a(i1)
      print *,i,i1,i2,cxs(j)
   enddo
   enddo

end program
