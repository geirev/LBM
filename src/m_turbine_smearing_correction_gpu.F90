! CUDA Fortran batch near-wake solver.  One block per blade, one thread per section.
! This module intentionally retains the original CPU solver's state equations.
module m_turbine_smearing_correction_gpu
#ifdef _CUDA
 use cudafor
 use m_liftdrag, only: liftdrag
#endif
 implicit none
 private
 public :: smear_gpu_batch, smear_gpu_reset
 integer, parameter :: maxsec=32
 real, parameter :: pi=acos(-1.0)
#ifdef _CUDA
 real, allocatable, device, save :: vn_d(:,:),vt_d(:,:),w_d(:,:),gh_d(:,:),dyn_d(:,:,:)
 real, allocatable, device, save :: xh_d(:,:,:),yh_d(:,:,:)
 real, allocatable, device, save :: ph_d(:,:),ps_d(:,:),ddx_d(:,:),ddy_d(:,:),a_d(:,:,:),sm_d(:,:)
 real, allocatable, device, save :: r_d(:),rv_d(:),ch_d(:)
 real, allocatable, device, save :: ua_d(:,:),ut_d(:,:),tw_d(:,:),pitch_d(:),uac_d(:,:),utc_d(:,:)
 integer, allocatable, device, save :: foil_d(:,:),iters_d(:)
 real, allocatable, device, save :: residual_d(:)
#endif
 integer, save :: ns_save=0,nb_save=0
 real, save :: eps_save=-1., dx_save=-1.
 real, parameter :: npcoef(4,5)=reshape([ &
 -1.64636754,-0.49900667,3.90835606,-1.60622844, &
 8.14821475,6.08465269,-18.76623305,7.42952712, &
 -12.17848617,-15.17120054,39.12432824,-15.85947771, &
 5.02652655,14.82541353,-29.48700861,11.68701640, &
 21.77131119,-2.42318681,-60.29473093,-195.06086584],[4,5])
 real, parameter :: nncoef(4,5)=reshape([ &
 1.01933206,12.98745129,-0.69015914,-0.26925341, &
 -0.13566815,49.99999999,101.238776,49.99999999, &
 0.39552471,0.00235345,-0.00154468,-0.00247976, &
 0.08018230,11.31160750,3.99520376,0.40364252, &
 44.83475035,3935.3432335,0.39454150,1.16610463],[4,5])
contains
subroutine smear_gpu_reset()
#ifdef _CUDA
 if (allocated(vn_d)) then
   deallocate(vn_d,vt_d,w_d,gh_d,dyn_d,xh_d,yh_d,ph_d,ps_d,ddx_d,ddy_d,a_d,sm_d, &
              r_d,rv_d,ch_d,ua_d,ut_d,tw_d,pitch_d,uac_d,utc_d,foil_d,iters_d,residual_d)
 endif
#endif
 ns_save=0;nb_save=0;eps_save=-1.;dx_save=-1.
end subroutine

subroutine smear_gpu_batch(nsec,nblade,dt,eps,r,chord,ua,ut,twist,foil,pitch,uac,utc, &
                           iters,residual)
 integer,intent(in)::nsec,nblade,foil(nsec,nblade)
 real,intent(in)::dt,eps,r(nsec),chord(nsec),ua(nsec,nblade),ut(nsec,nblade), &
                  twist(nsec,nblade),pitch(nblade)
 real,intent(out)::uac(nsec,nblade),utc(nsec,nblade)
 integer,intent(out)::iters(nblade)
 real,intent(out)::residual(nblade)
#ifdef _CUDA
 integer::nv,k,l,j
 real::rv(nsec+1),h,hr,dw,rootcorr,phi,phis,ddx,ddy,aa(4),s
 real,allocatable:: ph(:,:),ps(:,:),dxh(:,:),dyh(:,:),ah(:,:,:),sm(:,:)
 type(dim3)::grid,block
 if(nsec<2.or.nsec>maxsec) error stop 'smear_gpu_batch: unsupported nsec (2..32)'
 if(nblade<1) return
 if(eps<=0.) error stop 'smear_gpu_batch: eps<=0'
 if(any(r(2:nsec)<=r(1:nsec-1))) error stop 'smear_gpu_batch: radii not increasing'
 nv=nsec+1
 if(ns_save/=nsec.or.nb_save/=nblade.or.eps_save/=eps.or. &
    (allocated(vn_d).eqv..false.)) then
   call smear_gpu_reset()
   allocate(vn_d(nsec,nblade),vt_d(nsec,nblade),w_d(nsec,nblade),gh_d(nsec,nblade))
   allocate(dyn_d(nsec,3,nblade),xh_d(nv,nsec,nblade),yh_d(nv,nsec,nblade))
   allocate(ph_d(nv,nsec),ps_d(nv,nsec),ddx_d(nv,nsec),ddy_d(nv,nsec),a_d(nv,nsec,4),sm_d(nv,nsec))
   allocate(r_d(nsec),rv_d(nv),ch_d(nsec))
   allocate(ua_d(nsec,nblade),ut_d(nsec,nblade),tw_d(nsec,nblade),pitch_d(nblade))
   allocate(uac_d(nsec,nblade),utc_d(nsec,nblade),foil_d(nsec,nblade),iters_d(nblade),residual_d(nblade))
   vn_d=0.;vt_d=0.;w_d=0.;gh_d=0.;dyn_d=0.;xh_d=0.;yh_d=0.
   rv(2:nsec)=0.5*(r(1:nsec-1)+r(2:nsec))
   rv(1)=max(1.e-6,r(1)-0.5*(r(2)-r(1)))
   rv(nv)=r(nsec)+0.5*(r(nsec)-r(nsec-1))
   allocate(ph(nv,nsec),ps(nv,nsec),dxh(nv,nsec),dyh(nv,nsec),ah(nv,nsec,4),sm(nv,nsec))
   do l=1,nsec
     do k=1,nv
       h=rv(k)-r(l);hr=h/rv(k)
       call geometry(h,hr,phi,phis,ddx,ddy,aa)
       ph(k,l)=phi;ps(k,l)=phis;dxh(k,l)=ddx;dyh(k,l)=ddy;ah(k,l,:)=aa
       sm(k,l)=exp(-(h*h)/(eps*eps))
     enddo
   enddo
   ph_d=ph;ps_d=ps;ddx_d=dxh;ddy_d=dyh;a_d=ah;sm_d=sm
   r_d=r;rv_d=rv;ch_d=chord
   ns_save=nsec;nb_save=nblade;eps_save=eps
 endif
 ua_d=ua;ut_d=ut;tw_d=twist;foil_d=foil;pitch_d=pitch
 block=dim3(32,1,1);grid=dim3(nblade,1,1)
 call smear_kernel<<<grid,block>>>(nsec,nblade,dt,r_d,rv_d,ch_d,ua_d,ut_d,tw_d,foil_d,pitch_d, &
         ph_d,ps_d,ddx_d,ddy_d,a_d,sm_d,vn_d,vt_d,w_d,gh_d,dyn_d,xh_d,yh_d, &
         uac_d,utc_d,iters_d,residual_d)
 uac=uac_d;utc=utc_d;iters=iters_d;residual=residual_d
#else
 error stop 'smear_gpu_batch requires compilation with -D_CUDA'
#endif
end subroutine

subroutine geometry(h,hr,phi,phis,dx,dy,a)
 real,intent(in)::h,hr
 real,intent(out)::phi,phis,dx,dy,a(4)
 real::dw,rootcorr
 integer::j
 if(abs(h)<1.e-10.or.abs(hr)<1.e-10) then
   phi=1.e-8;phis=1.e-8;dx=0.;dy=0.;a=0.;return
 endif
 if(hr>0.) then
   phi=pi/4.*abs((1.+0.5*hr)*log(max(1.-hr,1.e-8)))
 else
   phi=log(max(1.-hr,1.e-8))/(1.5+log(max(1.-0.5*hr,1.e-8)))
 endif
 phi=max(abs(phi),1.e-8)
 dw=1./(4.*pi*abs(h)*hr)
 dx=1.359*dw;dy=0.359/4.*dw
 rootcorr=(dx*(1.-exp(-pi/(2.*phi)))-dy*(1.-exp(-2.*pi/phi)))/(dx-dy)
 dx=rootcorr*dx;dy=rootcorr*dy;phi=max(abs(rootcorr*phi),1.e-8)
 phis=0.788*abs(hr)
 if(hr<0.) then
   do j=1,4
     a(j)=nncoef(j,1)+nncoef(j,2)*exp(nncoef(j,3)*hr)+ &
          nncoef(j,4)*exp(nncoef(j,5)*hr)-nncoef(j,4)-nncoef(j,2)
   enddo
 else
   do j=1,4
     a(j)=npcoef(j,1)+npcoef(j,2)*hr+npcoef(j,3)*hr**2+npcoef(j,4)*hr**3
   enddo
 endif
end subroutine

#ifdef _CUDA
attributes(global) subroutine smear_kernel(ns,nb,dt,r,rv,ch,ua,ut,tw,fo,pitch, &
            ph,ps,ddx,ddy,a,sm,vnh,vth,wh,gh,dh,xh,yh,uac,utc,niter,resid)
 use m_liftdrag, only: liftdrag
 implicit none
 integer,value::ns,nb
 real,value::dt
 real,device::r(ns),rv(ns+1),ch(ns),ua(ns,nb),ut(ns,nb),tw(ns,nb),pitch(nb)
 integer,device::fo(ns,nb)
 real,device::ph(ns+1,ns),ps(ns+1,ns),ddx(ns+1,ns),ddy(ns+1,ns),a(ns+1,ns,4),sm(ns+1,ns)
 real,device::vnh(ns,nb),vth(ns,nb),wh(ns,nb),gh(ns,nb),dh(ns,3,nb)
 real,device::xh(ns+1,ns,nb),yh(ns+1,ns,nb),uac(ns,nb),utc(ns,nb),resid(nb)
 integer,device::niter(nb)
 real,shared::vn(maxsec),vt(maxsec),wold(maxsec),wnew(maxsec)
 real,shared::uax(maxsec),utan(maxsec),vr(maxsec),gsteady(maxsec),gcur(maxsec)
 real,shared::gdiff(maxsec+1),phiv(maxsec+1),vrelv(maxsec+1)
 real,shared::dyn(maxsec,3),err(maxsec),maxerr,relax
 integer::l,b,k,j,iter,km
 real::alpha,cl,cd,vrel,q,e1,e2,e3,dbeta,hr,kk,phistar,xe,ye,xnew,ynew,sumw,angle,oldx,oldy
 l=threadIdx%x;b=blockIdx%x
 if(b>nb)return
 if(l<=ns) then
   vn(l)=vnh(l,b);vt(l)=vth(l,b);wold(l)=wh(l,b)
 endif
 if(l==1)relax=0.2
 call syncthreads()
 do k=l,ns+1,blockDim%x
   if(k==1) then
     angle=atan2(ua(1,b),max(abs(ut(1,b)),1.e-8))
     if(ut(1,b)<0.)angle=atan2(ua(1,b),ut(1,b))
   else if(k==ns+1)then
     angle=atan2(ua(ns,b),max(abs(ut(ns,b)),1.e-8))
     if(ut(ns,b)<0.)angle=atan2(ua(ns,b),ut(ns,b))
   else
     angle=atan2(ua(k-1,b),max(abs(ut(k-1,b)),1.e-8))
     if(ut(k-1,b)<0.)angle=atan2(ua(k-1,b),ut(k-1,b))
     alpha=atan2(ua(k,b),max(abs(ut(k,b)),1.e-8))
     if(ut(k,b)<0.)alpha=atan2(ua(k,b),ut(k,b))
     angle=angle+(rv(k)-r(k-1))/(r(k)-r(k-1))*(alpha-angle)
   endif
   phiv(k)=angle
 enddo
 call syncthreads()
 do iter=1,500
   if(l<=ns)then
     uax(l)=ua(l,b)+vn(l);utan(l)=ut(l,b)-vt(l)
     vr(l)=sqrt(uax(l)*uax(l)+utan(l)*utan(l))
     alpha=atan2(uax(l),utan(l))*180./pi-tw(l,b)-pitch(b)
     call liftdrag(cl,cd,alpha,fo(l,b))
     gsteady(l)=0.5*vr(l)*ch(l)*cl
     q=dt/max(ch(l)/(2.*max(vr(l),1.e-8)),1.e-8)
     e1=exp(-0.3064*q);e2=exp(-0.0439*q);e3=exp(-3.227*q)
     dyn(l,1)=dh(l,1,b)*e1+0.27735*(gsteady(l)+gh(l,b))*(1.-e1)
     dyn(l,2)=dh(l,2,b)*e2+0.0914*(gsteady(l)+gh(l,b))*(1.-e2)
     dyn(l,3)=dh(l,3,b)*e3+0.13125*(gsteady(l)+gh(l,b))*(1.-e3)
     gcur(l)=dyn(l,1)+dyn(l,2)+dyn(l,3)
   endif
   call syncthreads()
   do k=l,ns+1,blockDim%x
     if(k==1)then
       gdiff(k)=gcur(1);vrelv(k)=vr(1)
     else if(k==ns+1)then
       gdiff(k)=-gcur(ns);vrelv(k)=vr(ns)
     else
       gdiff(k)=gcur(k)-gcur(k-1)
       vrelv(k)=vr(k-1)+(rv(k)-r(k-1))/(r(k)-r(k-1))*(vr(k)-vr(k-1))
     endif
   enddo
   call syncthreads()
   if(l<=ns)then
     sumw=0.
     do k=1,ns+1
       dbeta=max(vrelv(k)*dt/rv(k),0.)
       hr=(rv(k)-r(l))/rv(k)
       if(hr<0.)then
         kk=a(k,l,1)+a(k,l,2)*exp(a(k,l,3)*(pi/2.-phiv(k)))+ &
            a(k,l,4)*exp(-8.*(pi/2.-phiv(k)))-a(k,l,2)-a(k,l,4)
       else
         kk=a(k,l,1)+a(k,l,2)*phiv(k)+a(k,l,3)*phiv(k)**2+a(k,l,4)*phiv(k)**3
       endif
       kk=min(1.,max(0.,kk))
       phistar=max(kk*ps(k,l)+(1.-kk)*ph(k,l),1.e-8)
       xe=exp(-dbeta/phistar);ye=exp(-4.*dbeta/phistar)
       oldx=xh(k,l,b);oldy=yh(k,l,b)
       xnew=gdiff(k)*ddx(k,l)*phistar*(1.-xe)+oldx*xe
       ynew=gdiff(k)*ddy(k,l)*phistar*(1.-ye)+oldy*ye
       ! Do not update xh/yh until convergence (CPU uses previous-step history in each iteration).
       ! Recompute converged x/y once below.
       sumw=sumw+sm(k,l)*(xnew+ynew)
     enddo
     wnew(l)=relax*wold(l)+(1.-relax)*sumw
     err(l)=abs(wnew(l)-wold(l))
   endif
   call syncthreads()
   if(l==1)then
     maxerr=0.
     do j=1,ns
       maxerr=max(maxerr,err(j))
     enddo
   endif
   call syncthreads()
   if(l<=ns)then
     wold(l)=wnew(l)
     angle=atan2(ua(l,b)+vn(l),ut(l,b)-vt(l))
     vn(l)=wold(l)*cos(angle)
     vt(l)=wold(l)*sin(angle)
   endif
   call syncthreads()
   if(maxerr<1.e-6)exit
   if(iter==100.and.l==1)relax=0.4
   call syncthreads()
 enddo
 if(l<=ns)then
   uac(l,b)=ua(l,b)+vn(l);utc(l,b)=ut(l,b)-vt(l)
   vnh(l,b)=vn(l);vth(l,b)=vt(l);wh(l,b)=wold(l)
   gh(l,b)=gsteady(l)
   do j=1,3
     dh(l,j,b)=dyn(l,j)
   enddo
   do k=1,ns+1
     dbeta=max(vrelv(k)*dt/rv(k),0.)
     hr=(rv(k)-r(l))/rv(k)
     if(hr<0.)then
       kk=a(k,l,1)+a(k,l,2)*exp(a(k,l,3)*(pi/2.-phiv(k)))+ &
          a(k,l,4)*exp(-8.*(pi/2.-phiv(k)))-a(k,l,2)-a(k,l,4)
     else
       kk=a(k,l,1)+a(k,l,2)*phiv(k)+a(k,l,3)*phiv(k)**2+a(k,l,4)*phiv(k)**3
     endif
     kk=min(1.,max(0.,kk))
     phistar=max(kk*ps(k,l)+(1.-kk)*ph(k,l),1.e-8)
     xe=exp(-dbeta/phistar);ye=exp(-4.*dbeta/phistar)
     xh(k,l,b)=gdiff(k)*ddx(k,l)*phistar*(1.-xe)+xh(k,l,b)*xe
     yh(k,l,b)=gdiff(k)*ddy(k,l)*phistar*(1.-ye)+yh(k,l,b)*ye
   enddo
 endif
 if(l==1)then
   niter(b)=iter;resid(b)=maxerr
 endif
end subroutine
#endif
end module
