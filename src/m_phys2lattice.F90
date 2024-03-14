module m_phys2lattice

   type physconv
      real rho
      real length
      real time
      real vel
      real dynvisc
      real kinvisc
   end type
   real reynoldsnr
   real machnr

   type(physconv) p2l


contains

   subroutine phys2lattice()
      use m_readinfile

      real nu,newtau

      p2l%rho=1.225          ! Air density at 15C and  101.325 kPa  (kg/m^3)
      p2l%length=10.0        ! Length 10 m  (Rotor diameter is 200 m)
      p2l%vel=15.0           ! wind velocity (15 m/s)
      p2l%dynvisc=17.89E-06  ! Dynamic viscosity of air at 15C 17.89x10^-6  Ns/m^2   kg m/s^2  s/m^2 = kg /(m s)
      p2l%kinvisc=14.61E-06  ! Dynamic viscosity of air at 15C 14.61x10^-6  m^2/s


      p2l%time=p2l%length/p2l%vel

      reynoldsnr=(20.0*p2l%length)*p2l%vel/p2l%kinvisc

      print '(a,f12.4,a)','p2l_rho     =',p2l%rho      ,' [kg/m^2]'
      print '(a,f12.4,a)','p2l_length  =',p2l%length   ,' [m]'
      print '(a,f12.4,a)','p2l_vel     =',p2l%vel      ,' [m/s]'
      print '(a,f12.4,a)','p2l_time    =',p2l%time     ,' [s]'
      print '(a,f12.4,a)','p2l_dynvis  =',p2l%dynvisc  ,' [kg/(ms)]'
      print '(a,f12.4,a)','p2l_kinvis  =',p2l%kinvisc  ,' [m^2/s]'
      print '(a,i12  ,a)','Reynolds nr = ',nint(reynoldsnr)  ,' [ ]'

      nu=(1.0/3.0)*(tau - 0.5) * (p2l%length**2/p2l%time)
      print '(a,f12.4,a)','nu          =',nu           ,' [m^2/s]'

      newtau=0.5+ p2l%kinvisc * (p2l%time/p2l%length**2) / 3.0
      print '(a,f12.4,a)','newtau      =',newtau       ,' [ ]'

! We have 
!     nu* = cs*^2 (tau* -0.5)
!     nu* = (p2l%time/p2l%length**2) * nu

   end subroutine

end module

