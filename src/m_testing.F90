module m_testing
contains
subroutine testing(it,f,feq)
   use mod_dimensions, only : nx,ny,nz
   use mod_D3Q27setup, only : nl
   use m_readinfile,   only : ltesting
#ifdef MPI
   use m_mpi_decomp_init, only : mpi_rank
#endif
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
   logical ex
   integer iunit,i,j,k,l
   real fsum,fmax,eps,diff
   integer ir
   character(len=4) ctile
   character(len=6) cit
   character(len=3) ext
   character(len=5) suffix
   character(len=10) prefix
   character(len=10)  directory
   character(len=100) fname

   if (.not. ltesting) return
! File names
#ifdef MPI
   ir=mpi_rank
#else
   ir=0
#endif
   write(ctile,'(i4.4)') ir
   suffix = '_' // trim(ctile)
   ext='.uf'
   write(cit,'(i6.6)')it-1

   eps = sqrt(tiny(1.0))

   directory='testing/'
   call system('mkdir -p '//trim(directory))
   prefix='testing'

   fname = trim(directory) // trim(prefix) // trim(suffix) // '_' // trim(cit) // trim(ext)
   inquire(file=trim(fname),exist=ex)
   if (ex) then
      open(newunit=iunit,file=trim(fname),form="unformatted", status='old')
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
!        if (diff > 1.0E-5) then
!           print *,l,i,j,k,diff
!        endif
        fsum=fsum+diff
        fmax=max(diff,fmax)
     enddo
     enddo
     enddo
     enddo
     fsum=fsum/real(ntot)
     print '(a,i4,a,g12.5,a,g12.5)','Total misfit ',ir,': Mean abs error=',fsum,' max error=',fmax
   else
      open(newunit=iunit,file=trim(fname),form="unformatted", status='replace')
         print '(a,a)','Writing new:',trim(fname)
         f_h=f
         write(iunit)f_h
      close(iunit)
   endif

end subroutine
end module

