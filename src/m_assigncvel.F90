module m_assigncvel
! Not used but cool ordering of the lattice to simplify boundary conditions in the future
   integer :: cu(-1:1,-1:1,-1:1)
   integer :: cv(-1:1,-1:1,-1:1)
   integer :: cw(-1:1,-1:1,-1:1)
   real    :: cweights(-1:1,-1:1,-1:1)
contains
subroutine assigncvel()
   integer i,j,k,itot

   do k=-1,1
   do j=-1,1
   do i=-1,1
      cu(i,j,k)=i
      cv(i,j,k)=j
      cw(i,j,k)=k

      itot=abs(i)+abs(j)+abs(k)
      if (itot==0) then
         cweights(i,j,k)=8.0/27.0
      elseif (itot==1) then
         cweights(i,j,k)=2.0/27.0
      elseif (itot==2) then
         cweights(i,j,k)=1.0/54.0
      elseif (itot==3) then
         cweights(i,j,k)=1.0/216.0
      else
         stop 'error in itot ifstatement'
      endif

      print '(3i3,f10.3)',i,j,k,cweights(i,j,k)
   enddo
   enddo
   enddo
   stop


end subroutine
end module
