		subroutine I_Flux
!
		include 'common.h'
		dimension dw_b(8), dw_f(8), w_c(8), w_R(8,ni), w_L(8,ni)
		dimension f_g(4,ni), f_p(4,ni), tmatrix(4,4)
! for van Leer
		dimension fl(4), fr(4)
! for Roe
!		dimension delta_q(4), xinv(4,4), alpha(4), aramda(4), dissip(4)
!		dimension xm(4,4), fl(4), fr(4)
!
		include 'function.h'
!
		do j=1, nj-1
			do i=1, ni-1
				w_c(1) = fnr(i,j)
				w_c(2) = fnu(i,j)
				w_c(3) = fnv(i,j)
				w_c(4) = fnp(i,j)
				w_c(5) = fnr_p(i,j)
				w_c(6) = fnu_p(i,j)
				w_c(7) = fnv_p(i,j)
				w_c(8) = fnt_p(i,j)
!
				dw_b(1) = w_c(1) - fnr(i-1,j)
				dw_b(2) = w_c(2) - fnu(i-1,j)
				dw_b(3) = w_c(3) - fnv(i-1,j)
				dw_b(4) = w_c(4) - fnp(i-1,j)
				dw_b(5) = w_c(5) - fnr_p(i-1,j)
				dw_b(6) = w_c(6) - fnu_p(i-1,j)
				dw_b(7) = w_c(7) - fnv_p(i-1,j)
				dw_b(8) = w_c(8) - fnt_p(i-1,j)
!
				dw_f(1) = fnr(i+1,j) - w_c(1)
				dw_f(2) = fnu(i+1,j) - w_c(2)
				dw_f(3) = fnv(i+1,j) - w_c(3)
				dw_f(4) = fnp(i+1,j) - w_c(4)
				dw_f(5) = fnr_p(i+1,j) - w_c(5)
				dw_f(6) = fnu_p(i+1,j) - w_c(6)
				dw_f(7) = fnv_p(i+1,j) - w_c(7)
				dw_f(8) = fnt_p(i+1,j) - w_c(8)
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

					w_R(k,i  )=w_c(k)-0.5*phi*del_w_lim
					w_L(k,i+1)=w_c(k)+0.5*phi*del_w_lim
				end do
			end do
!
! boundary values
!
			w_R(1,ni)=fnr(ni,j)
			w_R(2,ni)=fnu(ni,j)
			w_R(3,ni)=fnv(ni,j)
			w_R(4,ni)=fnp(ni,j)
			w_R(5,ni)=fnr_p(ni,j)
			w_R(6,ni)=fnu_p(ni,j)
			w_R(7,ni)=fnv_p(ni,j)
			w_R(8,ni)=fnt_p(ni,j)
!
			w_L(1,1)=fnr(0,j)
			w_L(2,1)=fnu(0,j)
			w_L(3,1)=fnv(0,j)
			w_L(4,1)=fnp(0,j)
			w_L(5,1)=fnr_p(0,j)
			w_L(6,1)=fnu_p(0,j)
			w_L(7,1)=fnv_p(0,j)
			w_L(8,1)=fnt_p(0,j)
!
! Numerical Flux Evaluation
!
			do i=1, ni
				tmatrix=0.
				tmatrix(1,1)=1.
				tmatrix(2,2)= norm_i(1,i,j)
				tmatrix(2,3)=-norm_i(2,i,j)
				tmatrix(3,2)= norm_i(2,i,j)
				tmatrix(3,3)= norm_i(1,i,j)
				tmatrix(4,4)=1.
! gas-phase
				rL =w_L(1,i)
				uxL=w_L(2,i)
				uyL=w_L(3,i)
				pL =w_L(4,i)

				rR =w_R(1,i)
				uxR=w_R(2,i)
				uyR=w_R(3,i)
				pR =w_R(4,i)
!
				uL= norm_i(1,i,j)*uxL+norm_i(2,i,j)*uyL
				vL=-norm_i(2,i,j)*uxL+norm_i(1,i,j)*uyL
				uR= norm_i(1,i,j)*uxR+norm_i(2,i,j)*uyR
				vR=-norm_i(2,i,j)*uxR+norm_i(1,i,j)*uyR

!				include 'roe2d.h'
				include 'vanLeer2d.h'
!
				f_g(1:4,i)=fL(1:4)+fR(1:4)
!
! particle-phase
				r_p_L =w_L(5,i)
				ux_p_L=w_L(6,i)
				uy_p_L=w_L(7,i)
				t_p_L =w_L(8,i)

				r_p_R =w_R(5,i)
				ux_p_R=w_R(6,i)
				uy_p_R=w_R(7,i)
				t_p_R =w_R(8,i)
!
				un_p_L= norm_i(1,i,j)*ux_p_L+norm_i(2,i,j)*uy_p_L
				un_p_R= norm_i(1,i,j)*ux_p_R+norm_i(2,i,j)*uy_p_R
				u_p_L=0.5*(un_p_L+abs(un_p_L))
				u_p_R=0.5*(un_p_R-abs(un_p_R))

				v_p_L=-norm_i(2,i,j)*ux_p_L+norm_i(1,i,j)*uy_p_L
				v_p_R=-norm_i(2,i,j)*ux_p_R+norm_i(1,i,j)*uy_p_R

				fL(1) = r_p_L * u_p_L
				fL(2) = r_p_L * u_p_L * u_p_L
				fL(3) = r_p_L * u_p_L * v_p_L
				fL(4) = r_p_L * u_p_L * (cc*t_p_L+0.5*(u_p_L**2+v_p_L**2))

				fR(1) = r_p_R * u_p_R
				fR(2) = r_p_R * u_p_R * u_p_R
				fR(3) = r_p_R * u_p_R * v_p_R
				fR(4) = r_p_R * u_p_R * (cc*t_p_R+0.5*(u_p_R**2+v_p_R**2))
!
				f_p(1:4,i)=fL(1:4)+fR(1:4)
!
				f_g(1:4,i)=ds_i(i,j)*rave_i(i,j)*matmul(tmatrix(1:4,1:4),f_g(1:4,i))
				f_p(1:4,i)=ds_i(i,j)*rave_i(i,j)*matmul(tmatrix(1:4,1:4),f_p(1:4,i))
	
			end do
!			
			do i=1, ni-1
				rhs_g(1:4,i,j)=-(f_g(1:4,i+1)-f_g(1:4,i))*iflag(i,j)
				rhs_p(1:4,i,j)=-(f_p(1:4,i+1)-f_p(1:4,i))*iflag(i,j)
			end do
		end do
!
		return
		end