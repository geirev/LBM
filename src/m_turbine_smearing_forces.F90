module m_turbine_smearing_forces
   use m_turbine_smearing_correction, only : turbine_smear_init, turbine_smear_correct_blade
   implicit none
   logical, save :: smear_initialized=.false.

   ! IMPORTANT: keep identical to turbine_deposit/turbine_deposit_gpu.
   real, parameter :: epsilon_lattice=2.0

contains
subroutine turbine_smearing_forces(turbines,points,kin,Fvec,np,dx,uconv,dt,itimestep, &
                                    use_smearing,use_prandtl)
   use mod_turbines, only : turbine_t,point_t,pi
   use mod_turbine_def, only : nrchords,relm,chord,nblades,radiusnd,nrang
   use m_readinfile, only : actuator_model,nazim
   use m_liftdrag, only : liftdrag_host
   use m_turbine_compute_bladeforce, only : turbine_compute_bladeforce
   use m_turbine_blade_diagnostics, only : turbine_blade_diagnostics
   use m_turbine_rotor_basis, only : turbine_rotor_basis
#ifdef MPI
   use mpi, only : MPI_Comm_rank, MPI_COMM_WORLD
#endif
   implicit none

   type(turbine_t), intent(in) :: turbines(:)
   type(point_t), intent(in) :: points(:)
   integer, intent(in) :: np,itimestep
   real, intent(in) :: kin(3,np),dx,uconv,dt
   logical, intent(in) :: use_smearing,use_prandtl
   real, intent(out) :: Fvec(3,np)

   integer :: it,ib,ic,p0,p,nbl,niter
   real :: eps_si,resid,phi,cl,cd
   real :: ua(nrchords),ut(nrchords),uac(nrchords),utc(nrchords)
   real :: cl_diag(nrchords)
   real :: fn_diag(nrchords),ft_diag(nrchords)
   real :: e_axis(3),e1(3),e2(3),e_tan(3)
   integer :: myrank
#ifdef MPI
   integer :: ierr
#endif
   logical, save :: first_diagnostic=.true.
   real :: r_si(nrchords),c_si(nrchords)
   real :: twist_b(nrchords)
   integer :: foil_b(nrchords)

   myrank=0
#ifdef MPI
   call MPI_Comm_rank(MPI_COMM_WORLD,myrank,ierr)
#endif

   Fvec=0.0
   eps_si=epsilon_lattice*dx
   r_si=relm*dx
   c_si=chord*dx

   if (use_smearing .and. .not.smear_initialized) then
      call turbine_smear_init(size(turbines),nblades,nrchords,r_si,eps_si)
      smear_initialized=.true.
   endif

   ! points are ordered turbine -> blade/azimuth -> chord.
   do it=1,size(turbines)
      if (actuator_model==0) then
         nbl=nblades
      else
         nbl=nazim
      endif

      do ib=1,nbl
         p0=0
         do p=1,np
            if (points(p)%iturb==it .and. points(p)%iblade==ib .and. &
                points(p)%ichord==1) then
               p0=p
               exit
            endif
         enddo
         if (p0==0) cycle

         do ic=1,nrchords
            p=p0+ic-1
            ua(ic)=kin(1,p)*uconv
            ut(ic)=kin(2,p)*uconv
            twist_b(ic)=points(p)%twist
            foil_b(ic)=points(p)%foil
         enddo

         if (actuator_model==0 .and. use_smearing) then
            call turbine_smear_correct_blade(it,ib,dt,eps_si,r_si,c_si, &
                 twist_b,turbines(it)%pitchangle,foil_b,ua,ut,uac,utc, &
                 niter,resid)
         else
            uac=ua
            utc=ut
         endif

         do ic=1,nrchords
            p=p0+ic-1
            ! Convert corrected SI velocities back to LBM units before
            ! using the existing force routine. Density remains LBM density.
            phi=atan2(uac(ic),utc(ic))
            if (foil_b(ic) < 1) then
               write(*,*) 'ERROR foil index: it,ib,ic,p,foil = ', &
                          it,ib,ic,p,foil_b(ic)
               stop
            endif
            call liftdrag_host(cl,cd,phi*180.0/pi-points(p)%twist-points(p)%pitch, &
                          points(p)%foil)
            ! Use the common blade-force routine for all combinations of
            ! smearing correction and Prandtl tip loss.
            call turbine_compute_bladeforce(Fvec(:,p),points(p), &
                 uac(ic)/uconv,utc(ic)/uconv,kin(3,p),cl,cd, &
                 radiusnd,use_prandtl)

            Fvec(:,p)=Fvec(:,p)*points(p)%force_scale
            ! Use exactly the rotor basis and tangential direction from
            ! turbine_compute_bladeforce. Project the FINAL force vector,
            ! including Prandtl and force_scale, without recomputing loads.
            call turbine_rotor_basis(points(p)%yaw,points(p)%tilt, &
                 e_axis,e1,e2)
            e_tan=-sin(points(p)%theta)*e1 + cos(points(p)%theta)*e2
            cl_diag(ic)=cl
            fn_diag(ic)=dot_product(Fvec(:,p),e_axis)
            ft_diag(ic)=dot_product(Fvec(:,p),e_tan)
         enddo
         if (actuator_model==0 .and. myrank==0 .and. mod(itimestep,100)==0) then
            call turbine_blade_diagnostics( &
                 it,ib,itimestep,real(itimestep)*dt,nrchords, &
                 r_si,c_si,twist_b,turbines(it)%pitchangle, &
                 ua,ut,uac,utc,cl_diag,fn_diag,ft_diag, &
                 first_diagnostic)
         endif
         if (use_smearing .and. itimestep==1 .and. it==1 .and. ib==1) then
            write(*,'(a,f8.3,a,i4,a,es12.4)') &
                 'NW smearing correction: epsilon=',eps_si,' m, iterations=', &
                 niter,', residual=',resid
         endif
      enddo
   enddo
   if (myrank==0 .and. mod(itimestep,100)==0) first_diagnostic=.false.
end subroutine turbine_smearing_forces
end module m_turbine_smearing_forces
