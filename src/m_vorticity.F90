module m_vorticity
contains
subroutine vorticity(u,v,w,vortx,vorty,vortz,vort,blanking)
   use mod_dimensions
   real,    intent(in)  :: u(nx,ny,nz)
   real,    intent(in)  :: v(nx,ny,nz)
   real,    intent(in)  :: w(nx,ny,nz)
   real,    intent(out) :: vortx(nx,ny,nz)
   real,    intent(out) :: vorty(nx,ny,nz)
   real,    intent(out) :: vortz(nx,ny,nz)
   real,    intent(out) :: vort(nx,ny,nz)
   logical, intent(in)  :: blanking(nx,ny,nz)
   integer i,j,k,ia,ib,ja,jb,ka,kb

!   vortx = (cshift(w,1,2) - cshift(w,-1,2)) - (cshift(v,1,3) - cshift(v,-1,3))
!$OMP PARALLEL DO PRIVATE(i,j,k,ia,ib,ja,jb,ka,kb) SHARED(u,v,w)
      do k=1,nz
         ka=mod(nz+k-1-1,nz)+1
         kb=mod(nz+k-1+1,nz)+1
         do j=1,ny
            ja=mod(ny+j-1-1,ny)+1
            jb=mod(ny+j-1+1,ny)+1
            vortx(:,j,k) = (w(:,jb,k) - w(:,ja,k))   -   (v(:,j,kb) - v(:,j,ka))
         enddo
       enddo
!$OMP END PARALLEL DO

!   vorty = (cshift(u,1,3) - cshift(u,-1,3)) - (cshift(w,1,1) - cshift(w,-1,1))
!$OMP PARALLEL DO PRIVATE(i,j,k,ia,ib,ja,jb,ka,kb) SHARED(u,v,w)
       do k=1,nz
          ka=mod(nz+k-1-1,nz)+1
          kb=mod(nz+k-1+1,nz)+1
          do j=1,ny
             do i=1,nx
                ia=mod(nx+i-1-1,nx)+1
                ib=mod(nx+i-1+1,nx)+1
                vorty(i,j,k) = (u(i,j,kb) - u(i,j,ka))   -   (w(ib,j,k) - w(ia,j,k))
             enddo
          enddo
       enddo
!$OMP END PARALLEL DO

!   vortz = (cshift(v,1,1) - cshift(v,-1,1)) - (cshift(u,1,2) - cshift(u,-1,2))
!$OMP PARALLEL DO PRIVATE(i,j,k,ia,ib,ja,jb,ka,kb) SHARED(u,v,w)
       do k=1,nz
          do j=1,ny
             ja=mod(ny+j-1-1,ny)+1
             jb=mod(ny+j-1+1,ny)+1
             do i=1,nx
                ia=mod(nx+i-1-1,nx)+1
                ib=mod(nx+i-1+1,nx)+1
                vortz(i,j,k) = (v(ib,j,k) - v(ia,j,k))   -   (u(i,jb,k) - u(i,ja,k))
             enddo
          enddo
       enddo
!$OMP END PARALLEL DO

!$OMP PARALLEL DO PRIVATE(i,j,k) SHARED(vortx,vorty,vortz)
      do k=1,nz
         do j=1,ny
            do i=1,nx
               if (blanking(i,j,k)) then
                  vortx(i,j,k)=0.0
                  vorty(i,j,k)=0.0
                  vortz(i,j,k)=0.0
               endif
            enddo
         enddo
       enddo
!$OMP END PARALLEL DO

!$OMP PARALLEL DO PRIVATE(i,j,k) SHARED(vort,vortx,vorty,vortz)
      do k=1,nz
         do j=1,ny
            do i=1,nx
               vort(i,j,k) =sqrt(vortx(i,j,k)**2 + vorty(i,j,k)**2 + vortz(i,j,k)**2)
            enddo
         enddo
       enddo
!$OMP END PARALLEL DO

end subroutine
end module
