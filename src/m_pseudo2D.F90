module m_pseudo2D
contains
subroutine pseudo2D(Amat,nx,ny,lde,rx,ry,dx,dy,n1,n2,theta,verbose,lmean)
! This routine calculates the pseudo random filds using
! the procedure outlined in Evensen (1994,2009).
! But with analytic solution for r1 and r2 (Evensen, 2020).

! n1 should visely be selected as n1=nx*rx/dx unless periodic fields
! For periodic fields n1=nx. If rx=n*dx and n > 0.1*nx there is a risk for negative
! long range correlations.  Likewise for n2.

! Setting lmean=.false. generates random fields with mean approx equal zero for each realization
! by setting the first fourier amplitude equal to zero.

   use mod_fftw3
   implicit none
   logical, intent(in) :: verbose
   integer, intent(in) :: nx,ny           ! horizontal dimensions
   integer, intent(in) :: lde             ! number of fields stored at the time
   real, intent(out)   :: Amat(nx,ny,lde) ! generated random fields
   real, intent(in)    :: rx,ry           ! Horizontal decorrelation lengths (rx in dir theta)
   real, intent(in)    :: dx,dy           ! grid spacing
   real, intent(in)    :: theta           ! rotation angle in deg (theta=0 is east, rotation anticlocwise)
   integer, intent(inout) :: n1,n2           ! horizontal dimensions in fft grid (even numbers)
   logical, intent(in), optional :: lmean ! random mean

   real r1,r2,c

   integer(kind=8) plan

   integer l,p,j,m,i
   real kappa2,lambda2,kappa,lambda
   real pi2,deltak,summ
   real a11tmp,a22tmp,a11,a22,a12,torad
   real, parameter :: r8 = 8.0

   real, allocatable    :: fampl(:,:,:)
   real, allocatable    :: phi(:,:)
   real, allocatable    :: y(:,:)   ! Physical field
   complex, allocatable :: x(:,:)   ! Fourier amplitudes

   real, parameter :: pi=3.141592653589
   logical lm
   real e

   if (lde < 1)    stop 'pseudo2D: error lde < 1'
   if (rx <= 0.0)  stop 'pseudo2D: error, rx <= 0.0'
   if (ry <= 0.0)  stop 'pseudo2D: error, ry <= 0.0'
   if (n1 < nx)    stop 'pseudo2D: n1 < nx'
   if (n2 < ny)    stop 'pseudo2D: n2 < ny'

   if (mod(n1,2) /= 0) n1=2*int(real(n1+1)/2.0) ! make n1 even
   if (mod(n2,2) /= 0) n2=2*int(real(n2+1)/2.0) ! make n2 even

   if (present(lmean)) then
      lm=lmean    ! optional override or confirm default
   else
      lm=.true.   ! default is true
   endif

   allocate(fampl(2,0:n1/2,-n2/2:n2/2))
   allocate(phi(0:n1/2,-n2/2:n2/2))
   allocate(y(0:n1+1,0:n2-1))
   allocate(x(0:n1/2,0:n2-1))

   pi2=2.0*pi
   deltak=pi2**2/(real(n1*n2)*dx*dy)
   kappa=pi2/(real(n1)*dx)
   kappa2=kappa**2
   lambda=pi2/(real(n2)*dy)
   lambda2=lambda**2

   if (allocated(y)) deallocate(y)
   allocate(y(0:n1-1,0:n2-1))
#ifdef DOUBLE_PRECISION
   call dfftw_plan_dft_c2r_2d(plan,n1,n2,x,y,FFTW_ESTIMATE)
#else
   call sfftw_plan_dft_c2r_2d(plan,n1,n2,x,y,FFTW_ESTIMATE)
#endif


! computing the coefficients r1, r2, and c

   r1 = sqrt(real(r8, kind=kind(rx))) / rx
   r2 = sqrt(real(r8, kind=kind(ry))) / rx
!   r1=sqrt(8.0_c_float)/rx
!   r2=sqrt(8.0_c_float)/ry

!   if (verbose) print '(a,2f12.5,2i6,4f10.2)','pseudo2D: Call newton with ',r1,r2,n1,n2,dx,dy,rx,ry
!   call newton2D(r1,r2,n1,n2,dx,dy,rx,ry,cnv,verbose)
!   if (.not.cnv) then
!      stop 'newton did not converge'
!   endif

   summ=0.0
   do p=-n2/2+1,n2/2
   do l=-n1/2+1,n1/2
      summ=summ+exp(-2.0*(kappa2*real(l*l)/r1**2 + lambda2*real(p*p)/r2**2))
   enddo
   enddo

   if (.not.lm) summ=summ-1.0   ! Subtract one if mean is not sampled.

   c=sqrt(1.0/(deltak*summ))

   if (verbose) then
      print *,'pseudo2D: r1=  ',r1
      print *,'pseudo2D: r2=  ',r2
      print *,'pseudo2D:  c=  ',c
   end if

! Rotation to angle theta
   a11tmp=1.0/r1**2
   a22tmp=1.0/r2**2
   torad=-pi/180.0
   a11=a11tmp*cos(theta*torad)**2 + a22tmp*sin(theta*torad)**2
   a22=a11tmp*sin(theta*torad)**2 + a22tmp*cos(theta*torad)**2
   a12=(a22tmp-a11tmp)*cos(theta*torad)*sin(theta*torad)

   do j=1,lde
      ! Calculating the random wave phases
      call random_number(phi)
      phi=pi2*phi



      ! Calculating the wave amplitues
      do p=-n2/2,n2/2
      do l=0,n1/2
!         e=exp(-(kappa2*real(l*l)/r1**2+lambda2*real(p*p)/r2**2))
         e=exp(-( a11*kappa2*real(l*l) + 2.0*a12*kappa*lambda*real(l*p) + a22*lambda2*real(p*p) ))
         fampl(1,l,p)=e*cos(phi(l,p))*sqrt(deltak)*c
         fampl(2,l,p)=e*sin(phi(l,p))*sqrt(deltak)*c
      enddo
      enddo
      if (.not.lm) fampl(1,0,0)=0.0
      fampl(2,0,0)=0.0

      do p=0,n2/2-1
         x(:,p)=cmplx(fampl(1,:,p),fampl(2,:,p))
      enddo

      do p=n2/2,n2-1
         x(:,p)=cmplx(fampl(1,:,-n2+p),fampl(2,:,-n2+p))
      enddo

#ifdef DOUBLE_PRECISION
      call dfftw_execute_dft_c2r(plan,x,y)
#else
      call sfftw_execute_dft_c2r(plan,x,y)
#endif
      do m=1,ny
      do i=1,nx
         Amat(i,m,j)=y(i-1,m-1)
      enddo
      enddo

   enddo

#ifdef DOUBLE_PRECISION
   call dfftw_destroy_plan(plan)
#else
   call sfftw_destroy_plan(plan)
#endif

   deallocate(fampl, phi, y, x)

end subroutine
end module
