program test_cshift
    integer, dimension(3,3) :: a,b
    integer :: cs=1
    a = reshape( (/ 1, 2, 3, 4, 5, 6, 7, 8, 9 /), (/ 3, 3 /))
    print '(3i3)', a(1,:)
    print '(3i3)', a(2,:)
    print '(3i3)', a(3,:) 
!    b = cshift(a, -1, 1)
    b = cshift(a, cs, 2)
    print *
    print '(3i3)', b(1,:)
    print '(3i3)', b(2,:)
    print '(3i3)', b(3,:)

    do j=1,3
       j1=mod(3+j-1-cs,3)+1
       j2=mod(3+j-1+cs,3)+1
       print *,'aaa',j1,j,j2
    do i=1,3
       b(i,j)=a(i,j1)
    enddo
    enddo

    print *
    print '(3i3)', b(1,:)
    print '(3i3)', b(2,:)
    print '(3i3)', b(3,:)

end program test_cshift
