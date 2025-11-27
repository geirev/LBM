module m_mechanical_ablvisc
   use mod_dimensions, only : nz
   real, save :: ablvisc(nz)
contains
subroutine mechanical_ablvisc(ir)
   use mod_dimensions, only : nz
   use m_readinfile,   only : uini,p2l,iablvisc,ablheight
   implicit none

   integer, intent(in) :: ir
   real, parameter   :: kappa=0.4       ! von Karman constant
   real, parameter   :: z=100.0         ! typical height
   real, parameter   :: z0=0.1          ! roughness z0

   real height,uref
   integer k


   if (iablvisc /= 1) then
      ablvisc(:)=0.0
      return
   endif

   uref=kappa*uini*p2l%vel/log(z/z0)

   do k=1,nz
      height=real(k-1)*p2l%length
      if (height < ablheight) then
         ablvisc(k)=kappa*uref*height*(1.0-height/ablheight)**2/p2l%visc
      else
         ablvisc(k)=0.0
      endif
   enddo

   if (ir==1) then
      open(10, file='ablvisc.dat')
         do k=1,nz
            height=real(k-1)*p2l%length
            write(10,'(i4,3g13.5)')k,height,ablvisc(k),ablvisc(k)*p2l%visc
         enddo
      close(10)
   endif
end subroutine
end module 
