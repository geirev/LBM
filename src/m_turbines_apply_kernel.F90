module m_turbines_apply_kernel
contains

#ifdef _CUDA
   attributes(global) &
#endif
   subroutine turbines_apply_kernel(f, rho, u, v, w, F_turb, inv1cs2, inv2cs4, ratio, ibgk,&
                                    t_imin,t_imax, t_jmin,t_jmax, t_kmin,t_kmax)
#ifdef _CUDA
   use cudafor
#endif
   use mod_dimensions,  only : nx, ny, nz
   use mod_D3Q27setup,  only : nl, cxs, cys, czs, weights, H2, H3
   implicit none

   !-----------------------------------------------------------------
   ! Arguments
   !-----------------------------------------------------------------
   real, intent(inout) :: f(nl,0:nx+1,0:ny+1,0:nz+1)
   real, intent(in)    :: rho(0:nx+1,0:ny+1,0:nz+1)
   real, intent(in)    ::   u(0:nx+1,0:ny+1,0:nz+1)
   real, intent(in)    ::   v(0:nx+1,0:ny+1,0:nz+1)
   real, intent(in)    ::   w(0:nx+1,0:ny+1,0:nz+1)
   real, intent(in)    :: F_turb(3,0:nx+1,0:ny+1,0:nz+1)

   real, value         :: inv1cs2, inv2cs4, ratio
   integer, value      :: ibgk
   integer, value      :: t_imin,t_imax, t_jmin,t_jmax, t_kmin,t_kmax

   !-----------------------------------------------------------------
   ! Locals
   !-----------------------------------------------------------------
   integer :: i, j, k, l
   integer :: p, q, r

   real :: rr
   real :: Fx, Fy, Fz
   real :: vel0(3), vel1(3)
   real :: A0_2(3,3), A1_2(3,3)
   real :: A0_3(3,3,3), A1_3(3,3,3)
   real :: cu, tmp1, tmp2, tmp3
   real :: vratio

   real :: feq0, feq1
   real :: du, dv, dw

   !-----------------------------------------------------------------
   ! Map threads to (i,j,k)
   !-----------------------------------------------------------------
#ifdef _CUDA
!   i = t_imin + (blockIdx%x - 1) * blockDim%x + threadIdx%x - 1
!   j = t_jmin + (blockIdx%y - 1) * blockDim%y + threadIdx%y - 1
!   k = t_kmin + (blockIdx%z - 1) * blockDim%z + threadIdx%z - 1
   i = (blockIdx%x - 1) * blockDim%x + threadIdx%x
   j = (blockIdx%y - 1) * blockDim%y + threadIdx%y
   k = (blockIdx%z - 1) * blockDim%z + threadIdx%z

   ! Interior points only (assumes ghost cells at 0 and nx+1 etc.)
!   if (i < t_imin .or. i > t_imax) return
!   if (j < t_jmin .or. j > t_jmax) return
!   if (k < t_kmin .or. k > t_kmax) return
   if (i < 1 .or. i > nx) return
   if (j < 1 .or. j > ny) return
   if (k < 1 .or. k > nz) return
#else
!   do k=t_kmin,t_kmax
!   do j=t_jmin,t_jmax
!   do i=t_imin,t_imax
   do k=1,nz
   do j=1,ny
   do i=1,nx
#endif



   !-----------------------------------------------------------------
   ! Read local force and skip if zero
   !-----------------------------------------------------------------
   Fx = F_turb(1,i,j,k)
   Fy = F_turb(2,i,j,k)
   Fz = F_turb(3,i,j,k)

   if (Fx == 0.0 .and. Fy == 0.0 .and. Fz == 0.0) return

   !-----------------------------------------------------------------
   ! Original density and velocity
   !-----------------------------------------------------------------
   rr      = rho(i,j,k)
   if (rr <= 0.0) return   ! safety

   vel0(1) = u(i,j,k)
   vel0(2) = v(i,j,k)
   vel0(3) = w(i,j,k)

   !-----------------------------------------------------------------
   ! 1) Equilibrium with original velocity: feq0
   !    Build A0_2 and (optionally) A0_3
   !-----------------------------------------------------------------
   do q = 1,3
      do p = 1,3
         A0_2(p,q) = rr * vel0(p) * vel0(q) * inv2cs4
      enddo
   enddo

   if (ibgk == 3) then
      do r = 1,3
         vratio = vel0(r) * ratio
         do q = 1,3
            do p = 1,3
               A0_3(p,q,r) = A0_2(p,q) * vratio
            enddo
         enddo
      enddo
   endif

   !-----------------------------------------------------------------
   ! 2) Compute forced velocity vel1 = vel0 + du, dv, dw
   !    Here: du = -F / rho (add dt or coefficients if needed)
   !-----------------------------------------------------------------
   du = -Fx / rr
   dv = -Fy / rr
   dw = -Fz / rr

   vel1(1) = vel0(1) + du
   vel1(2) = vel0(2) + dv
   vel1(3) = vel0(3) + dw

   !-----------------------------------------------------------------
   ! 3) Equilibrium with forced velocity: feq1
   !    Build A1_2 and (optionally) A1_3
   !-----------------------------------------------------------------
   do q = 1,3
      do p = 1,3
         A1_2(p,q) = rr * vel1(p) * vel1(q) * inv2cs4
      enddo
   enddo

   if (ibgk == 3) then
      do r = 1,3
         vratio = vel1(r) * ratio
         do q = 1,3
            do p = 1,3
               A1_3(p,q,r) = A1_2(p,q) * vratio
            enddo
         enddo
      enddo
   endif

   !-----------------------------------------------------------------
   ! 4) For each lattice direction: feq0, feq1 and apply Δf = feq1 - feq0
   !-----------------------------------------------------------------
   do l = 1, nl

      ! ----- feq0 -----
      cu   = real(cxs(l))*vel0(1) + real(cys(l))*vel0(2) + real(czs(l))*vel0(3)
      tmp1 = rr * (1.0 + cu*inv1cs2)

      tmp2 = 0.0
      do q = 1,3
         do p = 1,3
            tmp2 = tmp2 + H2(p,q,l) * A0_2(p,q)
         enddo
      enddo

      tmp3 = 0.0
      if (ibgk == 3) then
         do r = 1,3
            do q = 1,3
               do p = 1,3
                  tmp3 = tmp3 + H3(p,q,r,l) * A0_3(p,q,r)
               enddo
            enddo
         enddo
      endif

      feq0 = weights(l) * (tmp1 + tmp2 + tmp3)

      ! ----- feq1 -----
      cu   = real(cxs(l))*vel1(1) + real(cys(l))*vel1(2) + real(czs(l))*vel1(3)
      tmp1 = rr * (1.0 + cu*inv1cs2)

      tmp2 = 0.0
      do q = 1,3
         do p = 1,3
            tmp2 = tmp2 + H2(p,q,l) * A1_2(p,q)
         enddo
      enddo

      tmp3 = 0.0
      if (ibgk == 3) then
         do r = 1,3
            do q = 1,3
               do p = 1,3
                  tmp3 = tmp3 + H3(p,q,r,l) * A1_3(p,q,r)
               enddo
            enddo
         enddo
      endif

      feq1 = weights(l) * (tmp1 + tmp2 + tmp3)

      ! Apply turbine forcing directly to f
      f(l,i,j,k) = f(l,i,j,k) + (feq1 - feq0)

   enddo
#ifndef _CUDA
   enddo
   enddo
   enddo
#endif

   end subroutine

end module

