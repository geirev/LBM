module m_printdefines
contains
subroutine printdefines()

#ifdef _CUDA
   print*, "Running GPU version!"
#elif OPEN_MP
   print*, "Running OPEN_MP version!"
#else
   print*, "Running single CPU version!"
#endif

#ifdef D3Q19
   print*, "D3Q19 lattice"
#else
   print*, "D3Q27 lattice"
#endif

#ifdef DOUBLE_PRECISION
   print*, "Double precision code (-r8)"
#else
   print*, "Single precision code"
#endif

#ifdef NETCDF
   print*, "Netcdf libraries included"
#else
   print*, "No Netcdf"
#endif

end subroutine
end module
