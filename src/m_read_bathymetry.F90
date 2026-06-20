module m_read_bathymetry
contains
subroutine read_bathymetry(blanking,experiment)
   use mod_dimensions, only : nx, nyg, nz
   implicit none
   logical, intent(inout) :: blanking(0:nx+1,0:nyg+1,0:nz+1)
   character(len=*) experiment
   integer iunit,i,j,k

   open(newunit=iunit,file='bathymetry_'//trim(experiment)//'.uf',form='unformatted')
      read(iunit)i,j,k
      if (i == nx .and. j==nyg .and. k== nz) then
         read(iunit)blanking
      else
         print *,'Bathymetry dimension issue: ',i,nx,j,nyg,k,nz
      endif
   close(iunit)
end subroutine
end module
