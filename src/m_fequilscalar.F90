module m_fequilscalar
contains
function fequilscalar(rho, u, v, w) result(feq)
   use mod_dimensions
   use mod_D3Q27setup
   use m_readinfile, only : ihrr
   implicit none
   real,    intent(in) :: rho
   real,    intent(in) :: u
   real,    intent(in) :: v
   real,    intent(in) :: w
   real feq(nl)


   real, save            :: H2(3,3,27)      ! Second order Hermite polynomial
   real, save            :: H3(3,3,3,27)    ! Third order Hermite polynomial
   real, save            :: c(3,nl)         ! Array storage of cxs, cys, and czs
   real                  :: A0_2(3,3)
   real                  :: A0_3(3,3,3)
   real                  :: delta(1:3, 1:3) = reshape([1, 0, 0, &
                                                       0, 1, 0, &
                                                       0, 0, 1], [3, 3])
   real                  :: vel(1:3),dens

   logical, save         :: lfirst=.true.
   integer l, p, q, r



   if (lfirst) then
      c(1,:)=real(cxs(:))
      c(2,:)=real(cys(:))
      c(3,:)=real(czs(:))

! Hermitian polynomials of 2nd order H2
      do l=1,nl
         do q=1,3
         do p=1,3
            H2(p,q,l)=c(p,l)*c(q,l) - cs2*delta(p,q)
         enddo
         enddo
      enddo

! Hermitian polynomials of 3nd order H3 with
! [c_\alpha \delta]_{ijk} = c_{\alpha,i} \delta_{jk} + c_{\alpha,j} \delta_{ik} +c_{\alpha,k} \delta_{ij}
      do l=1,nl
         do r=1,3
         do q=1,3
         do p=1,3
            H3(p,q,r,l)=c(p,l)*c(q,l)*c(r,l) - cs2*(c(p,l)*delta(q,r) + c(q,l)*delta(p,r) +  c(r,l)*delta(p,q))
         enddo
         enddo
         enddo
      enddo

      lfirst=.false.
   endif

   vel(1)=u
   vel(2)=v
   vel(3)=w
   dens=rho

! A0_2 and A0_3 from \citet{fen21a} (following Eq. 32)
   do q=1,3
   do p=1,3
      A0_2(p,q)=dens*vel(p)*vel(q)
   enddo
   enddo

   do r=1,3
   do q=1,3
   do p=1,3
      A0_3(p,q,r)=dens*vel(p)*vel(q)*vel(r)
   enddo
   enddo
   enddo

! Equilibrium distribution \citet{fen21a} Eq. (32) or jac18a eq (27)
   do l=1,nl
      feq(l)=dens

      feq(l)=feq(l) + dens*( c(1,l)*vel(1) + c(2,l)*vel(2) + c(3,l)*vel(3) )/cs2

      do p=1,3
      do q=1,3
         feq(l)=feq(l) + H2(p,q,l)*A0_2(p,q)/(2.0*cs4)
      enddo
      enddo
      ! the above identically recovers the BGK equilibrium, below we add third order contributions
      ! Add third order correction
      if (ihrr /= 2) then
         feq(l)=feq(l)   &
             + ( H3(1,1,2,l) + H3(2,3,3,l) ) * ( A0_3(1,1,2) + A0_3(2,3,3) )/(2.0*cs6) &
             + ( H3(1,3,3,l) + H3(1,2,2,l) ) * ( A0_3(1,3,3) + A0_3(1,2,2) )/(2.0*cs6) &
             + ( H3(2,2,3,l) + H3(1,1,3,l) ) * ( A0_3(2,2,3) + A0_3(1,1,3) )/(2.0*cs6) &
             + ( H3(1,1,2,l) - H3(2,3,3,l) ) * ( A0_3(1,1,2) - A0_3(2,3,3) )/(6.0*cs6) &
             + ( H3(1,3,3,l) - H3(1,2,2,l) ) * ( A0_3(1,3,3) - A0_3(1,2,2) )/(6.0*cs6) &
             + ( H3(2,2,3,l) - H3(1,1,3,l) ) * ( A0_3(2,2,3) - A0_3(1,1,3) )/(6.0*cs6)
      endif

      feq(l)= weights(l)*feq(l)

   enddo

end function
end module
