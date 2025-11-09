module m_diag
contains
subroutine diag(filetype,it,rho,u,v,w,lblanking,Ti)
   use mod_dimensions
   use m_readinfile, only : iout, iprt1, iprt2, dprt, nt1, lnodump
   use m_tecout
#ifdef NETCDF
   use m_netcdfout
#endif
#ifdef MPI
   use m_mpi_decomp_init, only : mpi_rank
#endif
   use m_wtime
   implicit none
   integer, intent(in)   :: filetype, it
   real,    intent(in)   :: rho(0:nx+1,0:ny+1,0:nz+1)
   real,    intent(in)   ::   u(0:nx+1,0:ny+1,0:nz+1)
   real,    intent(in)   ::   v(0:nx+1,0:ny+1,0:nz+1)
   real,    intent(in)   ::   w(0:nx+1,0:ny+1,0:nz+1)
   logical, intent(in)   :: lblanking(0:nx+1,0:ny+1,0:nz+1)
   real,    optional, intent(in) :: Ti(0:nx+1,0:ny+1,0:nz+1)
#ifdef _CUDA
   attributes(device) :: rho, u, v, w, lblanking, Ti
#endif

   ! Host scratch
   real,    allocatable :: u_h(:,:,:), v_h(:,:,:), w_h(:,:,:), rho_h(:,:,:), Ti_h(:,:,:)
   logical, allocatable :: lblanking_h(:,:,:)

   ! Output bookkeeping
   integer            :: num_of_vars
   character(len=200) :: variables
   character(len=10)  :: cit
   character(len=4)   :: ctile
   character(len=5)   :: suffix
   character(len=3)   :: prefix
   character(len=10)  :: directory
   character(len=4)   :: ext
   character(len=200) :: fname
   integer, parameter :: icpu = 14

   if (lnodump) return
   if (.not. (mod(it,iout)==0 .or. it==nt1 .or. ((it<=iprt1 .or. it>=iprt2) .and. mod(it,dprt)==0))) return

   call cpustart()

   ! -------------------------
   ! Decide "cit", variables, num_of_vars
   ! -------------------------
   select case (filetype)
   case (0,3)   ! Full field with indices
      if (present(Ti)) then
         cit         = 'F_AVERAGE'
         variables   = 'i,j,k,blanking,rho,u,v,w,Ti'
         num_of_vars = 9
      else
         write(cit,'(a1,i6.6)') 'F', it
         variables   = 'i,j,k,blanking,rho,u,v,w'
         num_of_vars = 8
      end if

   case (1)     ! Grid only
      cit         = 'GRID'
      variables   = 'i,j,k,blanking'
      num_of_vars = 4

   case (2)     ! Solution (no indices)
      if (present(Ti)) then
         cit         = '_AVERAGE'
         variables   = 'rho,u,v,w,Ti'
         num_of_vars = 5
      else
         write(cit,'(i6.6)') it
         variables   = 'rho,u,v,w'
         num_of_vars = 4
      end if

   case default
      print *, 'diag filetype:', filetype
      stop 'invalid filetype in diag'
   end select

   ! -------------------------
   ! Build filename pieces
   ! -------------------------
#ifdef MPI
   write(ctile,'(i4.4)') mpi_rank
   suffix = '_' // trim(ctile)
#else
   suffix = ''
#endif

#ifdef NETCDF
   prefix = 'out'
   ext    = '.nc'
#else
   prefix = 'tec'
   ext    = '.plt'
#endif
   directory='output/'

   call system('mkdir -p '//trim(directory))
   fname = trim(directory) // trim(prefix) // trim(suffix) // '_' // trim(cit) // trim(ext)

   ! -------------------------
   ! Host copies (needed by writers)
   ! -------------------------
   allocate(u_h(0:nx+1,0:ny+1,0:nz+1), v_h(0:nx+1,0:ny+1,0:nz+1), w_h(0:nx+1,0:ny+1,0:nz+1), rho_h(0:nx+1,0:ny+1,0:nz+1))
   allocate(lblanking_h(0:nx+1,0:ny+1,0:nz+1))

   u_h         = u
   v_h         = v
   w_h         = w
   rho_h       = rho
   lblanking_h = lblanking

   if (present(Ti)) then
      allocate(Ti_h(0:nx+1,0:ny+1,0:nz+1))
      Ti_h = Ti
   end if

   if (minval(rho_h) < 0.0) then
      print *,'iter=',it,'  minmaxrho=',minval(rho_h),' -- ',maxval(rho_h)
      print *,'iter=',it,'  minmaxloc=',minloc(rho_h),' -- ',maxloc(rho_h)
      stop 'Unstable simulation'
   end if

   ! -------------------------
   ! Single, non-duplicated writer dispatch
   ! -------------------------
   select case (filetype)

   case (0,1,2)
      if (present(Ti)) then
         call tecout(filetype, trim(fname), it, trim(variables), num_of_vars, &
                     lblanking_h, rho_h, u_h, v_h, w_h, Ti_h)
      else
         call tecout(filetype, trim(fname), it, trim(variables), num_of_vars, &
                     lblanking_h, rho_h, u_h, v_h, w_h)
      end if

   case (3)
#ifdef NETCDF
      if (present(Ti)) then
         call netcdfout(trim(fname), it, trim(variables), num_of_vars, &
                        lblanking_h, rho_h, u_h, v_h, w_h, Ti_h)
      else
         call netcdfout(trim(fname), it, trim(variables), num_of_vars, &
                        lblanking_h, rho_h, u_h, v_h, w_h)
      end if
#else
      print *, 'filetype=3 requires compiling with NETCDF'
      stop
#endif

   end select

   ! -------------------------
   ! Cleanup
   ! -------------------------
   deallocate(u_h, v_h, w_h, rho_h)
   if (allocated(lblanking_h)) deallocate(lblanking_h)
   if (allocated(Ti_h))        deallocate(Ti_h)

   call cpufinish(icpu)
end subroutine
end module

