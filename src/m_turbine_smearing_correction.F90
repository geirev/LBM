module m_turbine_smearing_correction
   use m_liftdrag, only : liftdrag_host
!-----------------------------------------------------------------------
! Near-wake / Gaussian-smearing correction for the LBM actuator-line model.
!
! Adapted from the Meyer Forsting / Pirrung near-wake correction.
!
! IMPORTANT:
!   - Called once per complete blade, not once per actuator point.
!   - All arguments are dimensional SI quantities:
!       r,chord,eps [m], uaxis,utan,omega*r [m/s], dt [s].
!   - utan is the relative tangential velocity used by the LBM ALM:
!       utan = omega*r - u_cfd . e_tan
!   - The routine returns corrected uaxis and utan.  The existing LBM
!     liftdrag() and turbine_compute_bladeforce() can then be used unchanged.
!   - This first LBM port assumes straight, unswept, unconed blades.
!   - Do NOT combine this initially with the Prandtl tip-loss multiplier;
!     validate the near-wake correction separately.
!-----------------------------------------------------------------------
   implicit none
   private

   public :: turbine_smear_init
   public :: turbine_smear_correct_blade
   public :: turbine_smear_reset

   real, parameter :: pi = acos(-1.0)

   logical, save :: initialized = .false.
   integer, save :: nt_save=0, nb_save=0, ns_save=0, nv_save=0

   ! History, indexed by turbine, blade, section/vortex.
   real, allocatable, save :: vn_hist(:,:,:), vt_hist(:,:,:), w_hist(:,:,:)
   real, allocatable, save :: gamma_hist(:,:,:)
   real, allocatable, save :: dyn_hist(:,:,:,:)
   real, allocatable, save :: x_hist(:,:,:,:), y_hist(:,:,:,:)

   ! Geometry-dependent near-wake coefficients.
   real, allocatable, save :: nw_phi(:,:,:,:), nw_phi_s(:,:,:,:)
   real, allocatable, save :: nw_dx(:,:,:,:), nw_dy(:,:,:,:)
   real, allocatable, save :: nw_a(:,:,:,:,:)
   real, allocatable, save :: smear0(:,:,:,:)

   real, save :: nw_P(4,5), nw_N(4,5)

contains

subroutine turbine_smear_init(nturb, nblades, nsections, r, eps)
!-----------------------------------------------------------------------
! Initialize history and geometry-dependent quantities.
!
! r(nsections) : actuator/control-point radii [m], increasing outward.
! eps          : Gaussian smearing length [m].  Use the SAME physical
!                epsilon as in turbine_deposit.
!
! Vortex trailing points are placed at section boundaries.  For a blade
! represented by section centres r(:), the interior boundaries are the
! midpoints.  End boundaries are extrapolated by half a section spacing.
!-----------------------------------------------------------------------
   integer, intent(in) :: nturb, nblades, nsections
   real, intent(in)    :: r(nsections), eps

   integer :: it, ib, k, l
   real, allocatable :: rv(:)
   real :: h, hr

   if (nsections < 2) stop 'turbine_smear_init: nsections < 2'
   if (eps <= 0.0)    stop 'turbine_smear_init: eps <= 0'
   if (any(r(2:nsections) <= r(1:nsections-1))) &
      stop 'turbine_smear_init: radii must increase'

   call turbine_smear_reset()

   nt_save=nturb; nb_save=nblades; ns_save=nsections; nv_save=nsections+1

   allocate(vn_hist(nturb,nblades,nsections)); vn_hist=0.0
   allocate(vt_hist(nturb,nblades,nsections)); vt_hist=0.0
   allocate(w_hist(nturb,nblades,nsections)); w_hist=0.0
   allocate(gamma_hist(nturb,nblades,nsections)); gamma_hist=0.0
   allocate(dyn_hist(nturb,nblades,nsections,3)); dyn_hist=0.0
   allocate(x_hist(nturb,nblades,nsections+1,nsections)); x_hist=0.0
   allocate(y_hist(nturb,nblades,nsections+1,nsections)); y_hist=0.0

   allocate(nw_phi(nturb,nblades,nsections+1,nsections)); nw_phi=0.0
   allocate(nw_phi_s(nturb,nblades,nsections+1,nsections)); nw_phi_s=0.0
   allocate(nw_dx(nturb,nblades,nsections+1,nsections)); nw_dx=0.0
   allocate(nw_dy(nturb,nblades,nsections+1,nsections)); nw_dy=0.0
   allocate(nw_a(nturb,nblades,nsections+1,nsections,4)); nw_a=0.0
   allocate(smear0(nturb,nblades,nsections+1,nsections)); smear0=0.0

   call set_nw_constants()
   allocate(rv(nsections+1))
   call make_vortex_radii(r,rv)

   do it=1,nturb
      do ib=1,nblades
         do k=1,nsections+1
            do l=1,nsections
               h  = rv(k)-r(l)
               hr = h/rv(k)
               call hr_components(h,hr,nw_phi(it,ib,k,l), &
                    nw_phi_s(it,ib,k,l),nw_dx(it,ib,k,l), &
                    nw_dy(it,ib,k,l),nw_a(it,ib,k,l,:))
               smear0(it,ib,k,l)=exp(-(h*h)/(eps*eps))
            enddo
         enddo
      enddo
   enddo

   deallocate(rv)
   initialized=.true.
end subroutine turbine_smear_init


subroutine turbine_smear_correct_blade(iturb, iblade, dt, eps, r, chord, &
                                       twist, pitch, foil, uaxis_cfd, &
                                       utan_cfd, uaxis_corr, utan_corr, &
                                       niter, resid)
!-----------------------------------------------------------------------
! Correct the local AL velocities for the induction lost through Gaussian
! force smearing.  One call handles ONE COMPLETE BLADE.
!
! Input velocities are dimensional [m/s].  utan_cfd is already the
! relative tangential velocity omega*r-u.e_tan used by the current LBM ALM.
!
! twist,pitch [deg]; foil(:) are the existing airfoil table identifiers.
!
! The present straight-blade port maps the near-wake correction directly
! onto rotor axial/tangential components.  This is the appropriate first
! implementation for the current rigid SWT/NREL blade representation.
!-----------------------------------------------------------------------
   integer, intent(in) :: iturb, iblade
   real, intent(in)    :: dt, eps
   real, intent(in)    :: r(:), chord(:), twist(:), pitch
   integer, intent(in) :: foil(:)
   real, intent(in)    :: uaxis_cfd(:), utan_cfd(:)
   real, intent(out)   :: uaxis_corr(:), utan_corr(:)
   integer, intent(out), optional :: niter
   real, intent(out), optional    :: resid

   integer :: ns,nv,l,k,iter
   real :: resmax, relax, relax_safety
   real :: h,hr,phi_cfd,phi_star,dbeta
   real :: xnew,ynew,xexp,yexp
   real :: cl,cd,alpha,vrel
   real, allocatable :: rv(:), phi_v(:), gamma(:), gamma_steady(:)
   real, allocatable :: dgamma(:), dyn(:,:), wold(:), wnew(:)
   real, allocatable :: vn(:),vt(:),vrelv(:), smear(:,:)
   real, allocatable :: xnew_all(:,:),ynew_all(:,:)


   ns=size(r); nv=ns+1
   if (.not.initialized) stop 'turbine_smear_correct_blade: not initialized'
   if (ns /= ns_save) stop 'turbine_smear_correct_blade: wrong nsections'
   if (iturb<1 .or. iturb>nt_save) stop 'turbine_smear_correct_blade: bad iturb'
   if (iblade<1 .or. iblade>nb_save) stop 'turbine_smear_correct_blade: bad iblade'

   allocate(rv(nv),phi_v(nv),gamma(ns),gamma_steady(ns),dgamma(nv))
   allocate(dyn(ns,3),wold(ns),wnew(ns),vn(ns),vt(ns),vrelv(nv))
   allocate(smear(nv,ns),xnew_all(nv,ns),ynew_all(nv,ns))

   call make_vortex_radii(r,rv)

   ! Wake/helix angle from the uncorrected CFD relative velocity.
   do l=1,ns
      phi_cfd=atan2(uaxis_cfd(l),max(abs(utan_cfd(l)),1.0e-8))
      if (utan_cfd(l)<0.0) phi_cfd=atan2(uaxis_cfd(l),utan_cfd(l))
      gamma(l)=phi_cfd
   enddo
   call interp_sections_to_vortices(r,gamma,rv,phi_v)

   ! In the straight-blade implementation delta=0, hence the smearing
   ! factor reduces to exp[-(rv-r)^2/eps^2].  Recompute here so eps may
   ! be changed between runs without rebuilding the module.
   do k=1,nv
      do l=1,ns
         h=rv(k)-r(l)
         smear(k,l)=exp(-(h*h)/(eps*eps))
      enddo
   enddo

   vn=vn_hist(iturb,iblade,:)
   vt=vt_hist(iturb,iblade,:)
   wold=w_hist(iturb,iblade,:)
   dyn=dyn_hist(iturb,iblade,:,:)
   relax_safety=0.2
   relax=0.2
   resmax=huge(1.0)

   do iter=1,500
      ! Corrected local velocity -> AoA -> circulation.
      do l=1,ns
         uaxis_corr(l)=uaxis_cfd(l)+vn(l)
         utan_corr(l) =utan_cfd(l) -vt(l)
         vrel=sqrt(uaxis_corr(l)**2+utan_corr(l)**2)
         alpha=atan2(uaxis_corr(l),utan_corr(l))*180.0/pi-twist(l)-pitch
         call liftdrag_host(cl,cd,alpha,foil(l))
         gamma_steady(l)=0.5*vrel*chord(l)*cl
      enddo

      call dynamic_circulation(dt,chord, &
           sqrt(uaxis_corr*uaxis_corr+utan_corr*utan_corr), &
           gamma_steady,gamma_hist(iturb,iblade,:), &
           dyn_hist(iturb,iblade,:,:),gamma,dyn)

      dgamma(2:ns)=gamma(2:ns)-gamma(1:ns-1)
      dgamma(1)=gamma(1)
      dgamma(nv)=-gamma(ns)

      call interp_sections_to_vortices(r, &
           sqrt(uaxis_corr*uaxis_corr+utan_corr*utan_corr),rv,vrelv)

      wnew=0.0; xnew_all=0.0; ynew_all=0.0
      do k=1,nv
         dbeta=max(vrelv(k)*dt/rv(k),0.0)
         do l=1,ns
            h=rv(k)-r(l); hr=h/rv(k)
            call phi_components(hr,nw_a(iturb,iblade,k,l,:),phi_v(k), &
                 nw_phi(iturb,iblade,k,l),nw_phi_s(iturb,iblade,k,l),phi_star)
            phi_star=max(phi_star,1.0e-8)
            xexp=exp(-dbeta/phi_star)
            yexp=exp(-4.0*dbeta/phi_star)
            xnew=dgamma(k)*nw_dx(iturb,iblade,k,l)*phi_star*(1.0-xexp) + &
                 x_hist(iturb,iblade,k,l)*xexp
            ynew=dgamma(k)*nw_dy(iturb,iblade,k,l)*phi_star*(1.0-yexp) + &
                 y_hist(iturb,iblade,k,l)*yexp
            xnew_all(k,l)=xnew; ynew_all(k,l)=ynew
            wnew(l)=wnew(l)+smear(k,l)*(xnew+ynew)
         enddo
      enddo

      ! Conservative fixed relaxation for the first LBM integration.
      ! The original source computes an adaptive relaxation factor; once
      ! the port is validated that routine can be added unchanged.
      wnew=relax*wold+(1.0-relax)*wnew
      resmax=maxval(abs(wnew-wold))
      wold=wnew

      ! Resolve correction along the local wake helix into axial and
      ! tangential rotor components.
      do l=1,ns
         phi_cfd=atan2(uaxis_cfd(l)+vn(l),utan_cfd(l)-vt(l))
         vn(l)=wold(l)*cos(phi_cfd)
         vt(l)=wold(l)*sin(phi_cfd)
      enddo

      if (resmax < 1.0e-6) exit
      if (iter==100) relax=min(0.5,relax+relax_safety)
   enddo

   uaxis_corr=uaxis_cfd+vn
   utan_corr =utan_cfd -vt

   ! Save converged history for next LBM time step.
   vn_hist(iturb,iblade,:)=vn
   vt_hist(iturb,iblade,:)=vt
   w_hist(iturb,iblade,:)=wold
   gamma_hist(iturb,iblade,:)=gamma_steady
   dyn_hist(iturb,iblade,:,:)=dyn
   x_hist(iturb,iblade,:,:)=xnew_all
   y_hist(iturb,iblade,:,:)=ynew_all

   if (present(niter)) niter=iter
   if (present(resid)) resid=resmax

   deallocate(rv,phi_v,gamma,gamma_steady,dgamma,dyn,wold,wnew,vn,vt, &
              vrelv,smear,xnew_all,ynew_all)
end subroutine turbine_smear_correct_blade


subroutine dynamic_circulation(dt_global,chord,vrel,gsteady,ghist,dhist, &
                               gamma,dyn)
   real, intent(in)  :: dt_global,chord(:),vrel(:),gsteady(:),ghist(:)
   real, intent(in)  :: dhist(:,:)
   real, intent(out) :: gamma(:),dyn(:,:)
   real :: tau(size(chord)), q(size(chord))

   tau=max(chord/(2.0*max(vrel,1.0e-8)),1.0e-8)
   q=dt_global/tau
   dyn(:,1)=dhist(:,1)*exp(-0.3064*q) + &
            0.27735*(gsteady+ghist)*(1.0-exp(-0.3064*q))
   dyn(:,2)=dhist(:,2)*exp(-0.0439*q) + &
            0.0914*(gsteady+ghist)*(1.0-exp(-0.0439*q))
   dyn(:,3)=dhist(:,3)*exp(-3.227*q) + &
            0.13125*(gsteady+ghist)*(1.0-exp(-3.227*q))
   gamma=dyn(:,1)+dyn(:,2)+dyn(:,3)
end subroutine dynamic_circulation


subroutine make_vortex_radii(r,rv)
   real,intent(in)  :: r(:)
   real,intent(out) :: rv(size(r)+1)
   integer :: n
   n=size(r)
   rv(2:n)=0.5*(r(1:n-1)+r(2:n))
   rv(1)=max(1.0e-6,r(1)-0.5*(r(2)-r(1)))
   rv(n+1)=r(n)+0.5*(r(n)-r(n-1))
end subroutine make_vortex_radii


subroutine interp_sections_to_vortices(r,y,rv,yv)
   real,intent(in)  :: r(:),y(:),rv(:)
   real,intent(out) :: yv(:)
   integer :: n,k,l
   real :: a
   n=size(r)
   yv(1)=y(1); yv(n+1)=y(n)
   do k=2,n
      l=k-1
      a=(rv(k)-r(l))/(r(l+1)-r(l))
      yv(k)=y(l)+a*(y(l+1)-y(l))
   enddo
end subroutine interp_sections_to_vortices


subroutine hr_components(h,hr,phi,phi_s,dx,dy,a)
   real,intent(in) :: h,hr
   real,intent(out):: phi,phi_s,dx,dy,a(4)
   real :: dw,root_corr
   integer :: k

   if (abs(h)<1.0e-10 .or. abs(hr)<1.0e-10) then
      phi=1.0e-8; phi_s=1.0e-8; dx=0.0; dy=0.0; a=0.0
      return
   endif

   if (hr>0.0) then
      phi=pi/4.0*abs((1.0+0.5*hr)*log(max(1.0-hr,1.0e-8)))
   else
      phi=log(max(1.0-hr,1.0e-8))/(1.5+log(max(1.0-0.5*hr,1.0e-8)))
   endif
   phi=max(abs(phi),1.0e-8)

   dw=1.0/(4.0*pi*abs(h)*hr)
   dx=1.359*dw
   dy=-0.359/4.*dw
   root_corr = (dx*(1.0-exp(-pi/(2.0*phi))) + &
                dy*(1.0-exp(-2.0*pi/phi))) / (dx+dy)
   phi = max(abs(root_corr*phi),1.0e-8)

   phi_s=0.788*abs(hr)
   if (hr<0.0) then
      do k=1,4
         a(k)=nw_N(k,1)+nw_N(k,2)*exp(nw_N(k,3)*hr)+ &
              nw_N(k,4)*exp(nw_N(k,5)*hr)-nw_N(k,4)-nw_N(k,2)
      enddo
   else
      do k=1,4
         a(k)=nw_P(k,1)+nw_P(k,2)*hr+nw_P(k,3)*hr**2+nw_P(k,4)*hr**3
      enddo
   endif
end subroutine hr_components


subroutine phi_components(hr,a,phi_flow,phi,phi_s,phi_star)
   real,intent(in) :: hr,a(4),phi_flow,phi,phi_s
   real,intent(out):: phi_star
   real :: kk
   if (hr<0.0) then
      kk=a(1)+a(2)*exp(a(3)*(pi/2.0-phi_flow))+ &
         a(4)*exp(-8.0*(pi/2.0-phi_flow))-a(2)-a(4)
   else
      kk=a(1)+a(2)*phi_flow+a(3)*phi_flow**2+a(4)*phi_flow**3
   endif
   kk=min(1.0,max(0.0,kk))
   phi_star=kk*phi_s+(1.0-kk)*phi
end subroutine phi_components


subroutine set_nw_constants()
   nw_P(1,:)=(/-1.64636754, 8.14821475,-12.17848617, 5.02652655, 21.77131119/)
   nw_P(2,:)=(/-0.49900667, 6.08465269,-15.17120054,14.82541353, -2.42318681/)
   nw_P(3,:)=(/ 3.90835606,-18.76623305,39.12432824,-29.48700861,-60.29473093/)
   nw_P(4,:)=(/-1.60622844, 7.42952712,-15.85947771,11.68701640,-195.06086584/)
   nw_N(1,:)=(/ 1.01933206,-0.13566815,0.39552471,0.08018230,44.83475035/)
   nw_N(2,:)=(/12.98745129,49.99999999,0.00235345,11.31160750,3935.3432335/)
   nw_N(3,:)=(/-0.69015914,101.238776,-0.00154468,3.99520376,0.39454150/)
   nw_N(4,:)=(/-0.26925341,49.99999999,-0.00247976,0.40364252,1.16610463/)
end subroutine set_nw_constants


subroutine turbine_smear_reset()
   if (allocated(vn_hist)) deallocate(vn_hist,vt_hist,w_hist,gamma_hist,dyn_hist, &
        x_hist,y_hist,nw_phi,nw_phi_s,nw_dx,nw_dy,nw_a,smear0)
   initialized=.false.
   nt_save=0; nb_save=0; ns_save=0; nv_save=0
end subroutine turbine_smear_reset

end module m_turbine_smearing_correction
