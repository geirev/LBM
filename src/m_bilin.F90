module m_bilin
contains
subroutine bilin(xb, yb, u, v, w, r, ic, jc, nx, ny, u_interp, v_interp, w_interp, r_interp)
    implicit none
    integer, intent(in) :: ic, jc, nx, ny
    real, intent(in) :: xb, yb
    real, intent(in) :: u(nx, ny), v(nx, ny), w(nx, ny), r(nx, ny)
    real, intent(out) :: u_interp, v_interp, w_interp, r_interp

    real :: t, u1, u2, v1, v2, w1, w2, r1 ,r2
    integer :: i1, i2, j1, j2

    ! Ensure indices are within bounds
    i1 = max(1, min(ic, nx-1))
    i2 = i1 + 1
    j1 = max(1, min(jc, ny-1))
    j2 = j1 + 1

    ! Interpolation factors
    t = (xb - real(i1)) / (real(i2) - real(i1))

    ! Bilinear interpolation for u
    u1 = u(i1, j1) + t * (u(i2, j1) - u(i1, j1))
    u2 = u(i1, j2) + t * (u(i2, j2) - u(i1, j2))
    u_interp = u1 + (yb - real(j1)) / (real(j2) - real(j1)) * (u2 - u1)

    ! Bilinear interpolation for v
    v1 = v(i1, j1) + t * (v(i2, j1) - v(i1, j1))
    v2 = v(i1, j2) + t * (v(i2, j2) - v(i1, j2))
    v_interp = v1 + (yb - real(j1)) / (real(j2) - real(j1)) * (v2 - v1)

    ! Bilinear interpolation for w
    w1 = w(i1, j1) + t * (w(i2, j1) - w(i1, j1))
    w2 = w(i1, j2) + t * (w(i2, j2) - w(i1, j2))
    w_interp = w1 + (yb - real(j1)) / (real(j2) - real(j1)) * (w2 - w1)

    ! Bilinear interpolation for r
    r1 = r(i1, j1) + t * (r(i2, j1) - r(i1, j1))
    r2 = r(i1, j2) + t * (r(i2, j2) - r(i1, j2))
    r_interp = r1 + (yb - real(j1)) / (real(j2) - real(j1)) * (r2 - r1)

end subroutine
end module
