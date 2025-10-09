module m_turbines_compute_force_kernel
contains
#ifdef _CUDA
  attributes(global) &
#endif
subroutine turbines_compute_force_kernel(force, forceN, forceT, thetain, jpos, kpos, iradius, relm)
   use mod_dimensions,  only : ny,nz
   use m_turbines_init, only : ieps
   use mod_nrel5mw,     only : nrchords
    implicit none
    ! device arrays / arguments
    real, intent(out)   :: force(0:ieps,ny,nz,3)
    real, intent(in)    :: forceN(nrchords,3)
    real, intent(in)    :: forceT(nrchords,3)
    real, intent(in)    :: relm(nrchords)
    real, intent(in), value     :: thetain
    integer, intent(in), value  :: jpos, kpos, iradius
#ifdef _CUDA
    attributes(device) :: force
    attributes(device) :: forceN
    attributes(device) :: forceT
    attributes(device) :: relm
#endif

    ! local and indexing
    integer :: i, j, k
    integer :: ja,jb,ka,kb
    integer :: iblade, ichord
    real :: y0, z0
    real :: yp, zp, costheta, sintheta
    real :: gauss
    real, parameter :: pi = 3.1415927410125732
    real, parameter :: pi2 = 2.0*pi
    real, parameter :: rad120 = pi2*120.0/360.0
    real, parameter :: eps = 2.5
    real :: f1, f2, f3
    real :: theta_local

    ja=max(1,jpos-iradius-2)
    jb=min(ny,jpos+iradius+2)
    ka=max(1,kpos-iradius-2)
    kb=min(nz,kpos+iradius+2)

    ! compute global thread indices (CUDA Fortran are 1-based)
#ifdef _CUDA
    i = threadIdx%x + (blockIdx%x - 1) * blockDim%x  -1  ! i=0...0eps
    j = threadIdx%y + (blockIdx%y - 1) * blockDim%y      ! 1..ny
    k = threadIdx%z + (blockIdx%z - 1) * blockDim%z      ! 1..nz

    ! simple bounds checks
    if (i < 0 .or. i > ieps) return
    if (j < ja .or. j > jb) return
    if (k < ka .or. k > kb) return
#else
!$OMP PARALLEL DO PRIVATE(i,j,k,f1,f2,f3,theta_local,costheta,sintheta,ichord,iblade,yp,zp,gauss)&
!$OMP  &SHARED(force, jpos, kpos, y0,z0, relm)
      do k=1,nz
      do j=1,ny
      do i=0,ieps
#endif

         ! set center coords and rotation
         y0 = real(jpos)
         z0 = real(kpos)

         ! initialize local accumulator for three components
         f1 = 0.0
         f2 = 0.0
         f3 = 0.0

         ! loop over three blades; each thread computes contributions of all chords for its (i,j,k)
         ! We follow the same rotation increments as CPU code: theta, theta+120deg, theta+240deg
         theta_local = thetain
         do iblade = 1, 3
            costheta = cos(theta_local)
            sintheta = sin(theta_local)

            ! loop chords
            do ichord = 1, nrchords
               yp = y0 + relm(ichord) * costheta
               zp = z0 + relm(ichord) * sintheta

               if (ieps == 0) then
                  gauss = (1.0/(pi*eps**2)) * exp( - ((real(j)-yp)**2 + (real(k)-zp)**2) / eps**2 )
               else
                  gauss = (1.0/(sqrt(pi**3)*eps**3)) * exp( - ((real(j)-yp)**2 + (real(k)-zp)**2 + real(i)**2) / eps**2 )
               endif

               f1 = f1 + forceN(ichord,iblade) * gauss
               f2 = f2 - forceT(ichord,iblade) * sintheta * gauss
               f3 = f3 + forceT(ichord,iblade) * costheta * gauss
            end do

            theta_local = theta_local + rad120
         end do

         force(i,j,k,1) =  f1
         force(i,j,k,2) =  f2
         force(i,j,k,3) =  f3
#ifndef _CUDA
      enddo
      enddo
      enddo
#endif


end subroutine
end module 
