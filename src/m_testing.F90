module m_testing
contains
subroutine testing(it,f,feq)
   use mod_dimensions, only : nx,ny,nz
   use mod_D3Q27setup, only : nl
   use m_readinfile,   only : ltesting
   implicit none
   integer, intent(in)  :: it
   integer, parameter   :: ntot=nl*(nx+2)*(ny+2)*(nz+2)
   real,    intent(in)  :: f(nl,0:nx+1,0:ny+1,0:nz+1)
   real,    intent(out) :: feq(nl,0:nx+1,0:ny+1,0:nz+1)
#ifdef _CUDA
   attributes(device) :: f
   attributes(device) :: feq
#endif
   real(kind=4) :: f_h(nl,0:nx+1,0:ny+1,0:nz+1)
   character(len=6) cit
   logical ex
   integer iunit,i,j,k,l
   real fsum,fmax,eps,diff
   if (.not. ltesting) return

   eps = sqrt(tiny(1.0))

   write(cit,'(i6.6)')it-1
   inquire(file='testing'//cit//'.uf',exist=ex)

   if (ex) then
      open(newunit=iunit,file='testing'//cit//'.uf',form="unformatted", status='old')
         read(iunit)f_h
      close(iunit)
      feq=f_h
#ifdef _CUDA
!$cuf kernel do(1) <<<*,*>>>
#endif
     fsum=0.0
     fmax=0.0
     do k=1,nz
     do j=1,ny
     do i=1,nx
     do l=1,nl
        diff=abs(f(l,i,j,k)-feq(l,i,j,k))
!        if (diff > 1.0E-4) then
!           print *,l,i,j,k,diff
!        endif
        fsum=fsum+diff
        fmax=max(diff,fmax)
     enddo
     enddo
     enddo
     enddo
     fsum=fsum/real(ntot)
     print *,'Total misfit: ',fsum,fmax
   else
      open(newunit=iunit,file='testing'//cit//'.uf',form="unformatted", status='replace')
         f_h=f
         write(iunit)f_h
      close(iunit)
   endif

end subroutine
end module

