module m_turbine_deposit_gpu_kernel
contains
#ifdef _CUDA
attributes(global) &
#endif
subroutine turbine_deposit_gpu_kernel( F_turb, xg, yg, zg, Fvec, np, nx, ny, nz, nyg, j_start, j_end, krad, sigma2)
  implicit none
#ifdef _CUDA
  integer, parameter :: tpb = 256
#endif
  integer, value :: np, nx, ny, nz, nyg, j_start, j_end, krad
  real,    value :: sigma2

#ifdef _CUDA
  real, device :: F_turb(3,0:nx+1,0:ny+1,0:nz+1)
  real, device :: xg(np), yg(np), zg(np)
  real, device :: Fvec(3,np)
#else
  real :: F_turb(3,0:nx+1,0:ny+1,0:nz+1)
  real :: xg(np), yg(np), zg(np)
  real :: Fvec(3,np)
#endif

  integer :: p, tid, idx, stride
  integer :: i0, j0, k0
  integer :: ii1, ii2, jj1g, jj2g, kk1, kk2
  integer :: ni, nj, nk, nTot
  integer :: ii, jj, kk, dj, di, dk
  integer :: jj_local
  real    :: xgp, ygp, zgp
  real    :: fp1, fp2, fp3
  real    :: dx, dy, dz, r2, w, arg
  real    :: sum_local, sumW
  real    :: old
  real, shared :: shsum(tpb)

  p   = blockIdx%x
  tid = threadIdx%x
  if (p < 1 .or. p > np) return


  fp1 = Fvec(1,p); fp2 = Fvec(2,p); fp3 = Fvec(3,p)
  if (fp1 == 0.0 .and. fp2 == 0.0 .and. fp3 == 0.0) return

  xgp = xg(p); ygp = yg(p); zgp = zg(p)

  i0 = int(floor(xgp))
  j0 = int(floor(ygp))
  k0 = int(floor(zgp))

  ii1  = max(1,   i0 - krad);  ii2  = min(nx,  i0 + krad)
  jj1g = max(1,   j0 - krad);  jj2g = min(nyg, j0 + krad)
  kk1  = max(1,   k0 - krad);  kk2  = min(nz,  k0 + krad)

  ni = ii2 - ii1 + 1
  nj = jj2g - jj1g + 1
  nk = kk2 - kk1 + 1
  nTot = ni*nj*nk
  if (nTot <= 0) return

  ! Pass 1: sumW
  sum_local = 0.0
  stride = blockDim%x
  idx = tid - 1
  do while (idx < nTot)
     di = mod(idx, ni)
     dj = mod(idx/ni, nj)
     dk = idx/(ni*nj)
     ii = ii1 + di
     jj = jj1g + dj
     kk = kk1 + dk

     dx = real(ii) - xgp
     dy = real(jj) - ygp
     dz = real(kk) - zgp
     r2 = dx*dx + dy*dy + dz*dz

     sum_local = sum_local + exp( -r2 / (2.0*sigma2) )
     idx = idx + stride
  end do

  shsum(tid) = sum_local
  call syncthreads()

  idx = stride/2
  do while (idx > 0)
     if (tid <= idx) shsum(tid) = shsum(tid) + shsum(tid+idx)
     call syncthreads()
     idx = idx/2
  end do

  sumW = shsum(1)
  if (sumW < 1.0e-6) return

  ! Pass 2: deposit
  idx = tid - 1
  do while (idx < nTot)
     di = mod(idx, ni)
     dj = mod(idx/ni, nj)
     dk = idx/(ni*nj)
     ii = ii1 + di
     jj = jj1g + dj
     kk = kk1 + dk

     if (jj >= j_start .and. jj <= j_end) then
        jj_local = jj - j_start + 1

        if (jj_local >= 1 .and. jj_local <= ny) then
           dx = real(ii) - xgp
           dy = real(jj) - ygp
           dz = real(kk) - zgp
           r2 = dx*dx + dy*dy + dz*dz

           arg = r2 / (2.0*sigma2)
           if (arg > 80.0) then
              w = 0.0
           else
              w = exp(-arg) / sumW
           end if

           old = atomicAdd(F_turb(1,ii,jj_local,kk), w * fp1)
           old = atomicAdd(F_turb(2,ii,jj_local,kk), w * fp2)
           old = atomicAdd(F_turb(3,ii,jj_local,kk), w * fp3)
        end if
     end if

     idx = idx + stride
  end do

end subroutine
end module

