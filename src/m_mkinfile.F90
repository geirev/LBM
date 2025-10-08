module m_mkinfile
contains
subroutine mkinfile()

open(10,file='infile.in',status='new')
   write(10,'(a)')'# Development variables'
   write(10,'(a)')' F                ! ltiming       : CPU timing'
   write(10,'(a)')' F                ! ltesting      : testing final solution while developing code.'
   write(10,'(a)')' F                ! lnodump       : No saving of diagnostics or restarts while optimizing GPU code'
   write(10,'(a)')' 256 1 1          ! ntx, nty, ntz : number of threads per block in x, y, and z direction'
   write(10,'(a)')'# Experiment configuration'
   write(10,'(a)')' windfarm         ! experiment    : Cylinder, airfoil etc'
   write(10,'(a)')' 3                ! ibgk          : BGK order, ibgk=2 standard second order, ibgk=3 gives third order BGK'
   write(10,'(a)')' 1                ! ihrr          : Collision operator ihrr=1 includes regularization by 3rd ord Hermite'
   write(10,'(a)')' 1 0.15           ! ivreman smagor: Vreman subgridscale mixing (1) with Smagorinsky constant (0.15)'
   write(10,'(a)')'# Time stepping and outputs'
   write(10,'(a)')' 00000            ! nt0           : First timestep'
   write(10,'(a)')' 200              ! nt1           : Final timestep'
   write(10,'(a)')' 5000             ! iout          : Number of steps between diag outputs'
   write(10,'(a)')' 10000            ! irestarts     : Number of steps between restart outputs'
   write(10,'(a)')' 00 60000 1       ! iprt1 iprt2 x : Output every x timestep for it <= iprt1 and it >= iprt2'
   write(10,'(a)')' 2                ! tecout        : full tecplot files (0), only solution (2), netcdf (3)'
   write(10,'(a)')'# Boundary conditions'
   write(10,'(a)')' 1                ! ibnd          : 0-periodic, 1 in/out flow,'
   write(10,'(a)')' 0                ! jbnd          : 0-periodic, 11,12,21,22 no-slip bb(1), free-slip bb(2) for j=1 and j=ny'
   write(10,'(a)')' 0                ! kbnd          : 0-periodic, 11,12,21,22 no-slip bb(1), free-slip bb(2) for k=1 and k=nz'
   write(10,'(a)')'# Inflow variables'
   write(10,'(a)')' 8.0 0.0          ! uini, udir    : Inflow wind velocity [m/s], direction in degrees (-45:45)'
   write(10,'(a)')' F 0.00005  100   ! lturb amp nrtu: Add turbulence forcing on inflow, amplitude, number of prestored time ste'
   write(10,'(a)')'# Physical variables'
   write(10,'(a)')' 0.0000178        ! visckin       : Dimensional kinematic viscosity'
   write(10,'(a)')' 1.225            ! C_rho - Density of air at surface 15C and  101.325 kPa  [kg/m^3] Eq. (7.12)'
   write(10,'(a)')' 4.0              ! C_l   - Length of a lattice cell in meters [m]'
   write(10,'(a)')' 75.0             ! C_u   - Wind velocity conversion [m/s]   -> C_t=C_l/C_u'
   write(10,'(a)')'# Averaging variables'
   write(10,'(a)')' F F              ! lave lavesec  : Switch on/off full averaging and turbine section averaging'
   write(10,'(a)')' 20000            ! avestart      : Iteration to start computing section averages'
   write(10,'(a)')' 40000            ! avestop       : Iteration to save section averages'
   write(10,'(a)')'# Turbine-definitions'
   write(10,'(a)')' 1                ! nturbines     : Number of turbines'
   write(10,'(a)')' 0.0              ! pitchangle    : Imposed pitchangle (0 until u=11.4, see table 7.1 in NREL doc).'
   write(10,'(a)')' 8.95             ! turbrpm       : Turbine RPM for actuator line model (max 12.1 9.22 8.95 12.06'
   write(10,'(a)')' 0.00             ! tipspeed ratio: Tipspeed ratio (7.55) (if given will override the given turbine RPM'
   write(10,'(a)')' 0                ! itiploss      : Tiploss(0-none, 1-Prandl, 2-Shen)'
   write(10,'(a)')' 96               ! ipos          : i-location turbine one'
   write(10,'(a)')' 61               ! jpos          : j-location turbine one'
   write(10,'(a)')' 61               ! kpos          : k-location turbine one'
   write(10,'(a)')' 64               ! ipos          : i-location turbine two'
   write(10,'(a)')' 75               ! jpos          : j-location turbine two'
   write(10,'(a)')' 48               ! kpos          : k-location turbine two'
   write(10,'(a)')
   write(10,'(a)')
   write(10,'(a)')' kinematic viscosity of air is 1.78E-5'
   write(10,'(a)')' Reynolds number becomes Re= u D /nu = 10^1 * 10^2 / 10^(-5)  = 10^7 - 10^8'
close(10)

end subroutine
end module
