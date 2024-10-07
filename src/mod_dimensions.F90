module mod_dimensions
   integer, parameter :: nx = 100          ! 928 resolution x-dir (east)
   integer, parameter :: ny = 96           ! resolution y-dir (north)
   integer, parameter :: nz = 96           ! resolution z-dir (up)
   integer, parameter :: nl = 27           ! number of componets of D3Q27
   integer, parameter :: ieps = 5          ! number of gridcells for smoothing actuatorline forcing in i-dir
   integer, parameter :: nrturb = 10       ! number of precomputed batches if inflow turbulence for u,v,w, and rho
end module
