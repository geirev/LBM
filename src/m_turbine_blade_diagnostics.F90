module m_turbine_blade_diagnostics
   implicit none
contains
!=======================================================================
! Write spanwise ALM blade diagnostics in Tecplot ASCII format.
!
! All velocities, radius and chord are SI (m/s, m).
! fn and ft are blade-element forces in LBM force units, as supplied
! by the caller; no physical-force conversion is performed here.
! Gamma = 0.5 * |Urel| * chord * CL [m^2/s], using corrected velocity.
! du = corrected velocity - uncorrected velocity [m/s].
!
! Call once per turbine/blade after final forces have been evaluated.
! Each call appends one ordered Tecplot zone to a blade-specific file.
! Set first=.true. for the first output to truncate an old file.
! All MPI ranks must NOT write the same files: call only on rank 0.
!=======================================================================
subroutine turbine_blade_diagnostics(iturb,iblade,itimestep,time,n, &
                                     radius,chord,twist,pitch, &
                                     ua,ut,uac,utc,cl,fn,ft,first)
   implicit none
   integer, intent(in) :: iturb,iblade,itimestep,n
   real, intent(in) :: time,pitch
   real, intent(in) :: radius(n),chord(n),twist(n)
   real, intent(in) :: ua(n),ut(n),uac(n),utc(n),cl(n),fn(n),ft(n)
   logical, intent(in) :: first
   integer :: iu,ic
   real :: rmax,phi0,phi,alpha0,alpha,urel,gamma
   character(len=128) :: filename

   if (n<1) return
   write(filename,'("output/blade_t",i3.3,"_b",i2.2,".dat")') iturb,iblade
   if (first) then
      open(newunit=iu,file=trim(filename),status='replace',action='write')
      write(iu,'(a)') 'TITLE="ALM blade diagnostics"'
      write(iu,'(a)') 'VARIABLES="r_m","r_over_R","alpha_deg",'// &
           '"alpha0_deg","CL","Gamma_m2s","FN_LBM","FT_LBM",'// &
           '"du_axial_ms","du_tangential_ms","ua_ms","ut_ms",'// &
           '"uac_ms","utc_ms","phi_deg","phi0_deg","Urel_ms"'
   else
      open(newunit=iu,file=trim(filename),status='unknown', &
           position='append',action='write')
   endif
   rmax=maxval(radius)
   write(iu,'(a,i0,a,i0,a,es16.8,a,i0,a)') &
        'ZONE T="step=',itimestep,' time=',nint(time), &
        '", SOLUTIONTIME=',time,', I=',n,', F=POINT'
   do ic=1,n
      phi0=atan2(ua(ic),ut(ic))
      phi=atan2(uac(ic),utc(ic))
      alpha0=phi0*180.0/acos(-1.0)-twist(ic)-pitch
      alpha=phi*180.0/acos(-1.0)-twist(ic)-pitch
      urel=sqrt(uac(ic)**2+utc(ic)**2)
      gamma=0.5*urel*chord(ic)*cl(ic)
      write(iu,'(17(1x,es18.9))') radius(ic),radius(ic)/max(rmax,tiny(rmax)), &
           alpha,alpha0,cl(ic),gamma,fn(ic),ft(ic), &
           uac(ic)-ua(ic),utc(ic)-ut(ic),ua(ic),ut(ic), &
           uac(ic),utc(ic),phi*180.0/acos(-1.0), &
           phi0*180.0/acos(-1.0),urel
   enddo
   close(iu)
end subroutine turbine_blade_diagnostics
end module m_turbine_blade_diagnostics
