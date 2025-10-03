module m_testing
contains
subroutine testing(it,f,feq)
   use mod_dimensions, only : nx,ny,nz
   use mod_D3Q27setup, only : nl
   use m_readinfile,   only : ltesting
   implicit none
   integer, intent(in)  :: it
   integer, parameter   :: ntot=nl*(nx+2)*(ny+2)*(nz+2)
   real,    intent(in)  :: f(ntot)
   real,    intent(out) :: feq(ntot)
#ifdef _CUDA
   attributes(device) :: f
   attributes(device) :: feq
#endif
   real(kind=4) :: f_h(ntot)
   character(len=6) cit
   logical ex
   integer iunit,i,j,k,l,idx
   real fsum,fmax,eps,diff
  integer :: strideI, strideJ, strideK
  integer :: remk, remj
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
     do idx=1,ntot
        diff=abs(f(idx)-feq(idx))
        if (diff > 1.0E-3) then
           strideI = nl
           strideJ = nl*(nx+2)
           strideK = nl*(nx+2)*(ny+2)

           k = (idx-1) / strideK + 1
           remk = mod(idx-1, strideK)
           j = remk / strideJ + 1
           remj = mod(remk, strideJ)
           I = remj / strideI + 1
           l = mod(remj, strideI) + 1

           print *,idx,i,j,k,l,diff
        endif
        fsum=fsum+diff
        fmax=max(diff,fmax)
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

