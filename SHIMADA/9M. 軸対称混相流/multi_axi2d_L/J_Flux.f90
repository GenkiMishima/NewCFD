		subroutine J_Flux
!
		include 'common.h'
		dimension dw_b(8), dw_f(8), w_c(8), w_R(8,nj), w_L(8,nj)
		dimension f_g(4,nj), f_p(4,nj), tmatrix(4,4)
! for van Leer
		dimension fl(4), fr(4)
! for Roe
!		dimension delta_q(4), xinv(4,4), alpha(4), aramda(4), dissip(4)
!		dimension xm(4,4), fl(4), fr(4)
!
		include 'function.h'
!
		do i=1, ni-1
			do j=1, nj-1
				w_c(1) = fnr(i,j)
				w_c(2) = fnu(i,j)
				w_c(3) = fnv(i,j)
				w_c(4) = fnp(i,j)
				w_c(5) = fnr_p(i,j)
				w_c(6) = fnu_p(i,j)
				w_c(7) = fnv_p(i,j)
				w_c(8) = fnt_p(i,j)
!
				dw_b(1) = w_c(1) - fnr(i,j-1)
				dw_b(2) = w_c(2) - fnu(i,j-1)
				dw_b(3) = w_c(3) - fnv(i,j-1)
				dw_b(4) = w_c(4) - fnp(i,j-1)
				dw_b(5) = w_c(5) - fnr_p(i,j-1)
				dw_b(6) = w_c(6) - fnu_p(i,j-1)
				dw_b(7) = w_c(7) - fnv_p(i,j-1)
				dw_b(8) = w_c(8) - fnt_p(i,j-1)
!
				dw_f(1) = fnr(i,j+1) - w_c(1)
				dw_f(2) = fnu(i,j+1) - w_c(2)
				dw_f(3) = fnv(i,j+1) - w_c(3)
				dw_f(4) = fnp(i,j+1) - w_c(4)
				dw_f(5) = fnr_p(i,j+1) - w_c(5)
				dw_f(6) = fnu_p(i,j+1) - w_c(6)
				dw_f(7) = fnv_p(i,j+1) - w_c(7)
				dw_f(8) = fnt_p(i,j+1) - w_c(8)
!
				do k=1, 8
					del_w_r=0.5*((1.+kappa)*dw_b(k)+(1.-kappa)*dw_f(k))
					del_w_l=0.5*((1.-kappa)*dw_b(k)+(1.+kappa)*dw_f(k))

					del_w_lim=0.
					if(dw_b(k).ne.0.) then
						theta=dw_f(k)/dw_b(k)
						if(theta.gt.0.) then
							dl=del_w_l*min(4.*theta/(1.-kappa+(1.+kappa)*theta),1.)
							dr=del_w_r*min(4./(1.+kappa+(1.-kappa)*theta),1.)
							del_w_lim=0.5*(sign(1.,dl)+sign(1.,dr))*min(abs(dl),abs(dr))
						end if
					end if

					w_R(k,j  )=w_c(k)-0.5*phi*del_w_lim
					w_L(k,j+1)=w_c(k)+0.5*phi*del_w_lim
				end do
			end do
!
! boundary values
!
			w_R(1,nj)=fnr(i,nj)
			w_R(2,nj)=fnu(i,nj)
			w_R(3,nj)=fnv(i,nj)
			w_R(4,nj)=fnp(i,nj)
			w_R(5,nj)=fnr_p(i,nj)
			w_R(6,nj)=fnu_p(i,nj)
			w_R(7,nj)=fnv_p(i,nj)
			w_R(8,nj)=fnt_p(i,nj)
!
			w_L(1,1)=fnr(i,0)
			w_L(2,1)=fnu(i,0)
			w_L(3,1)=fnv(i,0)
			w_L(4,1)=fnp(i,0)
			w_L(5,1)=fnr_p(i,0)
			w_L(6,1)=fnu_p(i,0)
			w_L(7,1)=fnv_p(i,0)
			w_L(8,1)=fnt_p(i,0)
!
! Numerical Flux Evaluation
!
			do j=1, nj
				tmatrix=0.
				tmatrix(1,1)=1.
				tmatrix(2,2)= norm_j(1,i,j)
				tmatrix(2,3)=-norm_j(2,i,j)
				tmatrix(3,2)= norm_j(2,i,j)
				tmatrix(3,3)= norm_j(1,i,j)
				tmatrix(4,4)=1.
! gas-phase
				rL =w_L(1,j)
				uxL=w_L(2,j)
				uyL=w_L(3,j)
				pL =w_L(4,j)

				rR =w_R(1,j)
				uxR=w_R(2,j)
				uyR=w_R(3,j)
				pR =w_R(4,j)
!
				uL= norm_j(1,i,j)*uxL+norm_j(2,i,j)*uyL
				vL=-norm_j(2,i,j)*uxL+norm_j(1,i,j)*uyL
				uR= norm_j(1,i,j)*uxR+norm_j(2,i,j)*uyR
				vR=-norm_j(2,i,j)*uxR+norm_j(1,i,j)*uyR

!				include 'roe2d.h'
				include 'vanLeer2d.h'
!
				f_g(1:4,j)=fL(1:4)+fR(1:4)
!
! particle-phase
				r_p_L =w_L(5,j)
				ux_p_L=w_L(6,j)
				uy_p_L=w_L(7,j)
				t_p_L =w_L(8,j)

				r_p_R =w_R(5,j)
				ux_p_R=w_R(6,j)
				uy_p_R=w_R(7,j)
				t_p_R =w_R(8,j)
!
				un_p_L= norm_j(1,i,j)*ux_p_L+norm_j(2,i,j)*uy_p_L
				un_p_R= norm_j(1,i,j)*ux_p_R+norm_j(2,i,j)*uy_p_R
				u_p_L=0.5*(un_p_L+abs(un_p_L))
				u_p_R=0.5*(un_p_R-abs(un_p_R))

				v_p_L=-norm_j(2,i,j)*ux_p_L+norm_j(1,i,j)*uy_p_L
				v_p_R=-norm_j(2,i,j)*ux_p_R+norm_j(1,i,j)*uy_p_R

				fL(1) = r_p_L * u_p_L
				fL(2) = r_p_L * u_p_L * u_p_L
				fL(3) = r_p_L * u_p_L * v_p_L
				fL(4) = r_p_L * u_p_L * (cc*t_p_L+0.5*(u_p_L**2+v_p_L**2))

				fR(1) = r_p_R * u_p_R
				fR(2) = r_p_R * u_p_R * u_p_R
				fR(3) = r_p_R * u_p_R * v_p_R
				fR(4) = r_p_R * u_p_R * (cc*t_p_R+0.5*(u_p_R**2+v_p_R**2))
!
				f_p(1:4,j)=fL(1:4)+fR(1:4)
!
				f_g(1:4,j)=ds_j(i,j)*rave_j(i,j)*matmul(tmatrix(1:4,1:4),f_g(1:4,j))
				f_p(1:4,j)=ds_j(i,j)*rave_j(i,j)*matmul(tmatrix(1:4,1:4),f_p(1:4,j))
	
			end do
!			
			do j=1, nj-1
				rhs_g(1:4,i,j)=rhs_g(1:4,i,j)-(f_g(1:4,j+1)-f_g(1:4,j))*iflag(i,j)
				rhs_p(1:4,i,j)=rhs_p(1:4,i,j)-(f_p(1:4,j+1)-f_p(1:4,j))*iflag(i,j)
			end do
		end do
!
		return
		end