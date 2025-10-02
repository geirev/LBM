#ifdef NETCDF
module m_netcdfout
  use netcdf
  use mod_dimensions
  implicit none

contains

  subroutine netcdfout(filename, it, variables_string, num_of_variables, lblanking, rho, u, v, w, Ti)
    implicit none
    character(len=*), intent(in) :: filename         ! Output filename
    integer,          intent(in) :: it               ! Timestep index
    character(len=*), intent(in) :: variables_string ! Variables to print separated by ,
    integer,          intent(in) :: num_of_variables ! Number of variables to print

    real,     intent(in)             :: rho(nx,ny,nz)       ! fluid density
    real,     intent(in)             :: u(nx,ny,nz)         ! x component of fluid velocity
    real,     intent(in)             :: v(nx,ny,nz)         ! y component of fluid velocity
    real,     intent(in)             :: w(nx,ny,nz)         ! z component of fluid velocity
    logical,  intent(in)             :: lblanking(0:nx+1,0:ny+1,0:nz+1) ! blanking
    real,     intent(in), optional   :: Ti(nx,ny,nz)        ! Turbulent kinetic energy

    ! NetCDF variables
    integer :: ncid, varid_rho, varid_u, varid_v, varid_w, varid_Ti, varid_blank
    integer :: dimids(3)
    integer :: ierr
    real :: blanking(nx,ny,nz)
    integer :: i, j, k

    ! Create blanking array (convert logical to real)
    blanking = 0.0
    do k=1,nz
      do j=1,ny
        do i=1,nx
          if (lblanking(i,j,k)) blanking(i,j,k)=1.0
        end do
      end do
    end do

    ! Create NetCDF file
    ierr = nf90_create(filename, nf90_clobber, ncid)
    if (ierr /= nf90_noerr) stop "Error creating NetCDF file"

    ! Define dimensions
    ierr = nf90_def_dim(ncid, "x", nx, dimids(1))
    ierr = nf90_def_dim(ncid, "y", ny, dimids(2))
    ierr = nf90_def_dim(ncid, "z", nz, dimids(3))

    ! Define variables
    ierr = nf90_def_var(ncid, "rho", nf90_real, dimids, varid_rho)
    ierr = nf90_def_var(ncid, "u", nf90_real, dimids, varid_u)
    ierr = nf90_def_var(ncid, "v", nf90_real, dimids, varid_v)
    ierr = nf90_def_var(ncid, "w", nf90_real, dimids, varid_w)
    ierr = nf90_def_var(ncid, "blanking", nf90_real, dimids, varid_blank)
    if (present(Ti)) then
      ierr = nf90_def_var(ncid, "Ti", nf90_real, dimids, varid_Ti)
    end if

    ! End define mode
    ierr = nf90_enddef(ncid)

    ! Write data
    ierr = nf90_put_var(ncid, varid_rho, rho)
    ierr = nf90_put_var(ncid, varid_u, u)
    ierr = nf90_put_var(ncid, varid_v, v)
    ierr = nf90_put_var(ncid, varid_w, w)
    ierr = nf90_put_var(ncid, varid_blank, blanking)
    if (present(Ti)) then
      ierr = nf90_put_var(ncid, varid_Ti, Ti)
    end if

    ! Close file
    ierr = nf90_close(ncid)

    print *, "netcdfout: Finished writing NetCDF file ", trim(filename)

  end subroutine netcdfout

end module m_netcdfout
#endif
