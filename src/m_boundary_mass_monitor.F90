module m_boundary_mass_monitor
contains

subroutine boundary_mass_monitor(rho,u,v,it)

!-----------------------------------------------------------------------
! Driver for the mass / density / boundary-flux diagnostics.
!
! Call, for example, every 100 time steps:
!
!    if (mod(it,100) == 0) then
!       call boundary_mass_monitor( &
!            rho,u,v,it,j0_is_phys,jN_is_phys)
!    endif
!
!
! Quantities printed:
!
!   <rho>    global domain-mean density
!
!   Mass     total mass in the computational domain
!
!   Qi0      outward mass flux through i=0
!   QiN      outward mass flux through i=nx+1
!   Qj0      outward mass flux through j=0
!   QjN      outward mass flux through j=ny+1
!
!   Qnet     total outward mass flux
!
!            Qnet = Qi0 + QiN + Qj0 + QjN
!
!   dM/dt    measured change of total mass per lattice time step
!            since the previous call
!
!   residual mass-balance residual
!
!            residual = dM/dt + Qnet
!
! For exact conservation:
!
!            dM/dt = -Qnet
!
! and therefore residual should be approximately zero.
!
!
! Sign convention:
!
!   positive boundary flux = mass leaving domain
!   negative boundary flux = mass entering domain
!
!
! MPI:
!
! Local diagnostics from all MPI ranks are summed with MPI_Allreduce.
! Only rank 0 writes the diagnostics.
!
!-----------------------------------------------------------------------

   use mod_dimensions, only : nx,ny,nz

   use m_boundary_mass_monitor_kernel, &
       only : boundary_mass_monitor_kernel

#ifdef MPI
   use mpi
   use m_mpi_decomp_init, only : north, south
#endif

#ifdef _CUDA
   use cudafor
#endif

   implicit none


!-----------------------------------------------------------------------
! Arguments
!-----------------------------------------------------------------------

   real, intent(in) :: rho(0:nx+1,0:ny+1,0:nz+1)
   real, intent(in) :: u  (0:nx+1,0:ny+1,0:nz+1)
   real, intent(in) :: v  (0:nx+1,0:ny+1,0:nz+1)

   integer, intent(in) :: it

   logical :: j0_is_phys
   logical :: jN_is_phys


#ifdef _CUDA
   attributes(device) :: rho,u,v
#endif


!-----------------------------------------------------------------------
! Diagnostics
!
! local/global:
!
!   1 = mass
!   2 = Qi0
!   3 = QiN
!   4 = Qj0
!   5 = QjN
!-----------------------------------------------------------------------

   real :: diag_local(5)
   real :: diag_global(5)

#ifdef _CUDA
   real, device :: diag_device(5)
#endif


!-----------------------------------------------------------------------
! Derived diagnostics
!-----------------------------------------------------------------------

   real :: mass
   real :: rho_mean

   real :: Qi0
   real :: QiN
   real :: Qj0
   real :: QjN

   real :: Qnet

   real :: dMdt
   real :: residual

   real :: ncells_local
   real :: ncells_global


!-----------------------------------------------------------------------
! Previous diagnostic state.
!
! SAVE means these values survive between calls.
!-----------------------------------------------------------------------

   real, save :: mass_previous = 0.0

   integer, save :: it_previous = 0

   logical, save :: first_call = .true.


!-----------------------------------------------------------------------
! MPI
!-----------------------------------------------------------------------

#ifdef MPI
   integer :: ierr
   integer :: myrank
#endif


!-----------------------------------------------------------------------
! CUDA launch configuration
!-----------------------------------------------------------------------

#ifdef _CUDA

   integer :: tx,ty,tz
   integer :: bx,by,bz

#endif


!=======================================================================
! Number of cells represented by this MPI rank.
!
! This is summed over MPI rather than assuming anything about the
! global j dimension.
!=======================================================================

   ncells_local = real(nx)*real(ny)*real(nz)

#ifdef MPI
      j0_is_phys = (south == MPI_PROC_NULL)
      jN_is_phys = (north == MPI_PROC_NULL)
#else
      j0_is_phys = .true.
      jN_is_phys = .true.
#endif

#ifdef MPI

   call MPI_Allreduce( &
        ncells_local,ncells_global,1,MPI_REAL,MPI_SUM, &
        MPI_COMM_WORLD,ierr)

   call MPI_Comm_rank(MPI_COMM_WORLD,myrank,ierr)

#else

   ncells_global = ncells_local

#endif


!=======================================================================
! Compute local diagnostics
!=======================================================================

#ifdef _CUDA

!-----------------------------------------------------------------------
! Reset device diagnostic array.
!-----------------------------------------------------------------------

   diag_device = 0.0


!-----------------------------------------------------------------------
! Kernel launch configuration.
!-----------------------------------------------------------------------

   tx = 8
   ty = 8
   tz = 4

   bx = (nx + tx - 1)/tx
   by = (ny + ty - 1)/ty
   bz = (nz + tz - 1)/tz


!-----------------------------------------------------------------------
! Launch diagnostic kernel.
!-----------------------------------------------------------------------

   call boundary_mass_monitor_kernel &
        <<<dim3(bx,by,bz),dim3(tx,ty,tz)>>> &
        (rho,u,v,diag_device,j0_is_phys,jN_is_phys)


!-----------------------------------------------------------------------
! Copy result to host.
!-----------------------------------------------------------------------

   diag_local = diag_device


#else

!-----------------------------------------------------------------------
! CPU version.
!-----------------------------------------------------------------------

   diag_local = 0.0

   call boundary_mass_monitor_kernel( &
        rho,u,v,diag_local,j0_is_phys,jN_is_phys)

#endif


!=======================================================================
! Sum diagnostics over MPI ranks
!=======================================================================

#ifdef MPI

   call MPI_Allreduce( &
        diag_local,diag_global,5,MPI_REAL,MPI_SUM, &
        MPI_COMM_WORLD,ierr)

#else

   diag_global = diag_local

#endif


!=======================================================================
! Unpack diagnostics
!=======================================================================

   mass = diag_global(1)

   Qi0 = diag_global(2)
   QiN = diag_global(3)

   Qj0 = diag_global(4)
   QjN = diag_global(5)


!=======================================================================
! Mean density
!=======================================================================

   rho_mean = mass/ncells_global


!=======================================================================
! Net outward flux
!=======================================================================

   Qnet = Qi0 + QiN + Qj0 + QjN


!=======================================================================
! Actual mass change since previous diagnostic.
!
! LBM timestep is assumed to be one lattice time unit.
!=======================================================================

   if (first_call) then

      dMdt     = 0.0
      residual = 0.0

      first_call = .false.

   else

      if (it > it_previous) then

         dMdt = &
              (mass-mass_previous) / &
              real(it-it_previous)

      else

         dMdt = 0.0

      endif


!-----------------------------------------------------------------------
! Conservation requires
!
!       dM/dt + Qnet = 0
!-----------------------------------------------------------------------

      residual = dMdt + Qnet

   endif


!=======================================================================
! Write diagnostics
!=======================================================================

#ifdef MPI

   if (myrank == 0) then

#endif

      write(*,'(a)') &
         '----------------------------------------------------------------------'

      write(*,'(a,i12)') &
         'Boundary mass diagnostics, step = ',it

      write(*,'(a,es18.9)') &
         '  mean rho          = ',rho_mean

      write(*,'(a,es18.9)') &
         '  total mass        = ',mass

      write(*,'(a,es18.9)') &
         '  Qi0  outward      = ',Qi0

      write(*,'(a,es18.9)') &
         '  QiN  outward      = ',QiN

      write(*,'(a,es18.9)') &
         '  Qj0  outward      = ',Qj0

      write(*,'(a,es18.9)') &
         '  QjN  outward      = ',QjN

      write(*,'(a,es18.9)') &
         '  Qnet outward      = ',Qnet

      write(*,'(a,es18.9)') &
         '  dM/dt measured    = ',dMdt

      write(*,'(a,es18.9)') &
         '  dM/dt + Qnet      = ',residual

#ifdef MPI

   endif

#endif


!=======================================================================
! Save state for next call
!=======================================================================

   mass_previous = mass
   it_previous   = it


end subroutine boundary_mass_monitor

end module m_boundary_mass_monitor
