!==============================================================
!  m_turbine_deposit.F90
!  Deposit smoothed actuator forces on the lattice
!==============================================================
module m_turbine_deposit
   use mod_turbines
   use mod_dimensions, only : nx, ny, nz, nyg
#ifdef MPI
   use m_mpi_decomp_init, only : j_start, j_end
#endif
   implicit none
contains

!--------------------------------------------------------------
!  subroutine turbine_deposit
!
!  PURPOSE:
!    Deposit global point forces Fvec_global onto the
!    tile-local forcing field F_turb using a Gaussian stencil.
!
!  CALL:
!    call turbine_deposit(F_turb, points_global, Fvec_global, np)
!--------------------------------------------------------------
subroutine turbine_deposit(F_turb, points_global, Fvec_global, np)
   real,          intent(inout) :: F_turb(3,0:nx+1,0:ny+1,0:nz+1)
   integer,       intent(in)    :: np
   type(point_t), intent(in)    :: points_global(np)
   real,          intent(in)    :: Fvec_global(3, np)

   integer :: p, i0, j0, k0
   integer :: ii, jj, kk
   integer :: ii1, ii2, jj1g, jj2g, kk1, kk2
   integer :: jj_local
   real    :: dx, dy, dz, weight, sumW
   real    :: sigma, epsilon
   real    :: xg, yg, zg
   real    :: Fp(3)
   integer :: krad

   epsilon = 2.0
   sigma   = epsilon / sqrt(2.0)
   krad    = ceiling(3.0 * sigma)

   do p = 1, np
      xg = points_global(p)%xg
      yg = points_global(p)%yg
      zg = points_global(p)%zg

      Fp(:) = Fvec_global(:,p)
      if (Fp(1) == 0.0 .and. Fp(2) == 0.0 .and. Fp(3) == 0.0) cycle

      i0 = floor(xg)
      j0 = floor(yg)
      k0 = floor(zg)

      ii1  = max(1,   i0 - krad)
      ii2  = min(nx,  i0 + krad)
      jj1g = max(1,   j0 - krad)
      jj2g = min(nyg, j0 + krad)
      kk1  = max(1,   k0 - krad)
      kk2  = min(nz,  k0 + krad)

      sumW = 0.0
      do kk = kk1, kk2
         dz = real(kk) - zg
         do jj = jj1g, jj2g
            dy = real(jj) - yg
            do ii = ii1, ii2
               dx = real(ii) - xg
               sumW = sumW + exp(-(dx*dx + dy*dy + dz*dz)/(2.0*sigma*sigma))
            end do
         end do
      end do

      if (sumW < 1.0e-6) cycle

      do kk = kk1, kk2
         dz = real(kk) - zg
         do jj = jj1g, jj2g
#ifdef MPI
            if (jj < j_start .or. jj > j_end) cycle
            jj_local = jj - j_start + 1
#else
            jj_local = jj
#endif
            if (jj_local < 1 .or. jj_local > ny) cycle

            dy = real(jj) - yg

            do ii = ii1, ii2
               dx = real(ii) - xg

               weight = exp(-(dx*dx + dy*dy + dz*dz)/(2.0*sigma*sigma)) / sumW

               F_turb(1,ii,jj_local,kk) = F_turb(1,ii,jj_local,kk) + weight*Fp(1)
               F_turb(2,ii,jj_local,kk) = F_turb(2,ii,jj_local,kk) + weight*Fp(2)
               F_turb(3,ii,jj_local,kk) = F_turb(3,ii,jj_local,kk) + weight*Fp(3)
            end do
         end do
      end do
   end do
end subroutine turbine_deposit

end module m_turbine_deposit
