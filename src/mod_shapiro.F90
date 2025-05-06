module mod_shapiro
contains

subroutine shfact(n,sh)
   implicit none
   integer n,i,j
   real sh(0:n),f2n,fj,f2nj,ff
!-------------------------------------------------------------------
!      This routine calculates the factors in the general formula
!      for the Shapiro filter, given as
!
!                                      (-1)^{n+j-1}(2n)!
!        y(i)=x(i) + \Sum_{j=0}^{2n} { -------------------- }
!                                       2^{2n} j!  (2n-j)!
!
!      which is calculated from equation (39) in
!      (The multi-dimensional Crowley advection scheme,
!      Piotr K. Smolarkiewicz, Monthley weater review,
!      volume 110, pages 1968-1983, december 1982)
!      Note the error in the paper.
!
!       n      is the order of the filter. n=2^j  j=0,1,2, ...
!       sh     is the vector of dimension sh(0:2n) which
!              returns the factors in the filter.
!
!      Normally this routine is called first to compute the factors
!      and followed by "shfilt" for the filtering.
!
!                  Made by Geir Evensen juni 1991,
!                          geir@nrsc.no
!
!-------------------------------------------------------------------
   if (.not.(n==1 .or. n==2 .or. n==4 .or. n== 8 )) then
      write(*,*)'Error in "shfact"'
      write(*,'(a,i3)')'Invalid order of filter (1, 2, 4, 8): n=',n
      stop
   endif

   ff=2.0**(2*n)


   f2n=1
   do j=1,2*n
      f2n=f2n*real(j)
   enddo

   fj=1
   sh(0)=((-1.0)**(n-1))/ff
   do j=1,n

!        calculates (2n-j)! for each j
         f2nj=1
         do i=1,2*n-j
            f2nj=f2nj*real(i)
         enddo

!        calculates j!
         fj=fj*real(j)
!        calculates the factors
         sh(j)=(-1.0)**(n+j-1) *real( f2n/(ff*fj*f2nj) )
   enddo
   return
end subroutine shfact


subroutine shfilt(n,sh,ndim,x,incx,y,incy,shdim)
   implicit logical (a-z)
   integer n,i,j,incx,incy,ndim,ix,iy,shdim,m
   real sh(0:shdim),x(ndim*incx),y(ndim*incy),sum
!-------------------------------------------------------------------
!      This routine calculates the filtered vector y(i) from
!      the vector x(i), using the shapirofilter of order n,
!      and the method described in "shfact"
!
!      This routine must be called after the routine "shfact", which
!      computes the factors used in the filter.
!
!       n      The order of the filter.
!       sh     The vector of dimension sh(0:2n) which
!              returns the factors in the filter.
!       ndim   The number of elements in the vectors x and y.
!       x      The vector containing the elements to be filtered.
!       incx   The increment between elements which should be filtered.
!       y      The vector containing the filtered result.
!       incy   The increment between filtered elements in y.
!
!
!                  Made by Geir Evensen juni 1991,
!                          geir@nrsc.no
!
!-------------------------------------------------------------------
   if (n.le.0) return

   y(1)=x(1)
   y(ndim)=x(ndim)

   if ((incx.eq.1).and.(incy.eq.1)) then
!        code for both increments equal to 1
      do i=2,n
         sum=0.0
         do j=0,n-1
            m=1-(i-n+j)
            if (m.ge.0) then
               sum=sum+sh(j)*(2.0*x(1)-x(1+m) + x(i+n-j))
            else
               sum=sum+sh(j)*(x(i+n-j)+x(i-n+j))
            endif
         enddo
         y(i)=(1.0+sh(n))*x(i)+sum
      enddo

      do i=n+1,ndim-n
         sum=0.0
         do j=0,n-1
            sum=sum+sh(j)*(x(i+n-j)+x(i-n+j))
         enddo
         y(i)=(1.0+sh(n))*x(i)+sum
      enddo

      do i=ndim-n+1,ndim-1
         sum=0.0
         do j=0,n-1
            m=(i+n-j)-ndim
            if (m.ge.0) then
               sum=sum+sh(j)*( 2.0*x(ndim)-x(ndim-m) + x(i-n+j))
            else
               sum=sum+sh(j)*(x(i+n-j)+x(i-n+j))
            endif
         enddo
         y(i)=(1.0+sh(n))*x(i)+sum
      enddo


   else
!     code for unequal increments or equal increments not equal to 1
      do i=n+1,ndim-n
         ix=(i-1)*incx + 1
         iy=(i-1)*incy + 1

         sum=0.0
         do j=0,n-1
            sum=sum+sh(j)*(x(ix+(n-j)*incx)+x(ix-(n-j)*incx))
         enddo

         y(iy)=(1.0+sh(n))*x(ix)+sum
      enddo
   endif
end subroutine shfilt

subroutine shfilt2D(ish,sh,shdim,field,nx,ny)
   implicit none
   integer, intent(in) :: nx,ny
   integer, intent(in) :: ish
   integer, intent(in) :: shdim
   real,    intent(in) :: sh(shdim)
   real,    intent(inout) :: field(nx,ny)
   real, allocatable :: x(:),y(:)
   integer i,j
   integer :: incx=1
   integer :: incy=1
   allocate(x(nx),y(nx))
   do j=1,ny
      x(:)=field(:,j)
      call shfilt(ish,sh,nx,x,incx,y,incy,shdim)
      field(:,j)=y(:)
   enddo
   deallocate(x,y)
   allocate(x(ny),y(ny))
   do i=1,nx
      x(:)=field(i,:)
      call shfilt(ish,sh,nx,x,incx,y,incy,shdim)
      field(i,:)=y(:)
   enddo
   deallocate(x,y)
end subroutine

end module
