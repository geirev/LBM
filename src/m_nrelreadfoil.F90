module m_nrelreadfoil
contains

subroutine nrelreadfoil()
   use mod_nrel5mw
#ifdef _CUDA
   use cudafor
#endif
   implicit none

   integer :: i,k
   integer :: istat

   character(len=14) :: filenames(1:nrc)= [ 'Cylinder1.dat ', &
                                            'Cylinder2.dat ', &
                                            'DU40_A17.dat  ', &
                                            'DU35_A17.dat  ', &
                                            'DU30_A17.dat  ', &
                                            'DU25_A17.dat  ', &
                                            'DU21_A17.dat  ', &
                                            'NACA64_A17.dat' ]

   ! Reset host foil table
   foil_h(:,:)%deg = 0.0
   foil_h(:,:)%cl  = 0.0
   foil_h(:,:)%cd  = 0.0
   foil_h(:,:)%cm  = 0.0
   nrang_h(:)      = 0

   ! ------------------------------
   ! READ FOIL TABLES TO HOST
   ! ------------------------------
   do k=1,nrc
      open(10,file='./Airfoils/'//trim(filenames(k)))

      do i=1,50
         read(10,*)
      enddo

      read(10,*) nrang_h(k)
      read(10,*)
      read(10,*)

      do i=1,nrang_h(k)
         read(10,*) foil_h(i,k)
      enddo

      close(10)
   enddo

   ! ------------------------------
   ! COPY HOST -> DEVICE TABLES
   ! ------------------------------
   nrang = nrang_h
   foil  = foil_h
#ifdef _CUDA
   istat = cudaDeviceSynchronize()
#endif

end subroutine nrelreadfoil
end module m_nrelreadfoil
