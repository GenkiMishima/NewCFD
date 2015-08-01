		subroutine Point_Implicit
!
		include 'common.h'
		dimension tdq(8)
		include 'function.h'
!
		do j=1, nj-1
		do i=1, ni-1
			r_i = 0.5*(rave_i(i,j)+rave_i(i+1,j))
			r_j = 0.5*(rave_j(i,j)+rave_j(i,j+1))
			enx_i=0.5*(norm_i(1,i,j)+norm_i(1,i+1,j))
			enr_i=0.5*(norm_i(2,i,j)+norm_i(2,i+1,j))
			ddd=sqrt(enx_i**2+enr_i**2)
			enx_i=enx_i/ddd
			enr_i=enr_i/ddd
			enx_j=0.5*(norm_j(1,i,j)+norm_j(1,i,j+1))
			enr_j=0.5*(norm_j(2,i,j)+norm_j(2,i,j+1))
			ddd=sqrt(enx_j**2+enr_j**2)
			enx_j=enx_j/ddd
			enr_j=enr_j/ddd
			dsi=0.5*(ds_i(i,j)+ds_i(i+1,j))
			dsj=0.5*(ds_j(i,j)+ds_j(i,j+1))

			dratio=fnr_p(i,j) / fnr(i,j)
			temp=fnt(i,j)
			sos=fna(i,j)
			u=fnu(i,j)
			v=fnv(i,j)
			up=fnu_p(i,j)
			vp=fnv_p(i,j)
            un_i = enx_i * u + enr_i * v
			un_j = enx_j * u + enr_j * v
			unp_i= enx_i * up + enr_i * vp
			unp_j= enx_j * up + enr_j * vp

			umax_i_p=abs(unp_i)
			umax_j_p=abs(unp_j)
			umax_i_g=abs(un_i)+sos
			umax_j_g=abs(un_j)+sos

			un_p_i_1 = abs( norm_i(1,  i,  j)* up + norm_i(2,  i,  j) * vp ) 
			un_p_i_2 = abs( norm_i(1,i+1,  j)* up + norm_i(2,i+1,  j) * vp )
			un_p_j_1 = abs( norm_j(1,  i,  j)* up + norm_j(2,  i,  j) * vp )
			un_p_j_2 = abs( norm_j(1,  i,j+1)* up + norm_j(2,  i,j+1) * vp )
!
			dtmove = vol(i,j) / (ds_i(i,j)*rave_i(i,j)*un_p_i_1 + ds_i(i+1,j)*rave_i(i+1,j)*un_p_i_2 &
			    		        +ds_j(i,j)*rave_j(i,j)*un_p_j_1 + ds_j(i,j+1)*rave_j(i,j+1)*un_p_j_2)

			if(dtmove.gt.max(tau_v(i,j),tau_t(i,j))) then
				nu_a_g(i,j)=r_i*max(umax_i_p,umax_i_g)
				nu_b_g(i,j)=r_j*max(umax_j_p,umax_j_g)
				nu_a_p(i,j)=nu_a_g(i,j)
				nu_b_p(i,j)=nu_b_g(i,j)
			else
				nu_a_g(i,j)=r_i*umax_i_g
				nu_b_g(i,j)=r_j*umax_j_g
				nu_a_p(i,j)=r_i*umax_i_p
				nu_b_p(i,j)=r_j*umax_j_p
			end if
!
			alpha_g(i,j)=1.+dt(i,j)/vol(i,j)*(dsi*nu_a_g(i,j)+dsj*nu_b_g(i,j))
			alpha_p(i,j)=1.+dt(i,j)/vol(i,j)*(dsi*nu_a_p(i,j)+dsj*nu_b_p(i,j))
			rhs_g(1:4,i,j)=rhs_g(1:4,i,j)/alpha_g(i,j)
			rhs_p(1:4,i,j)=rhs_p(1:4,i,j)/alpha_p(i,j)
!
			rad=vol(i,j)/area(i,j)
			rhs_g(3,i,j)=1./(rad/dt(i,j)+gm*v)*(0.5*gm*(u**2+v**2)*rhs_g(1,i,j)+gm*u*rhs_g(2,i,j) &
						+rad/dt(i,j)*rhs_g(3,i,j)+gm*rhs_g(4,i,j))
!
			tau_vbydt=tau_v(i,j)/dt(i,j)
			fact=1./(1.+dratio+tau_vbydt)
!
			do k=1, 4
				tdq(k  )=rhs_g(k,i,j)
				tdq(k+4)=rhs_p(k,i,j)
			end do

			rhs_g(1,i,j)=tdq(1)
			rhs_g(2,i,j)=fact*(dratio*u*tdq(1)+(1.+tau_vbydt)*tdq(2) &
			            -u*tdq(5)+tdq(6))
			rhs_g(3,i,j)=fact*(dratio*v*tdq(1)+(1.+tau_vbydt)*tdq(3) &
			            -v*tdq(5)+tdq(7))
			rhs_p(1,i,j)=tdq(5)
			rhs_p(2,i,j)=fact*(-dratio*u*tdq(1)+dratio*tdq(2)+u*tdq(5)+(dratio+tau_vbydt)*tdq(6))
			rhs_p(3,i,j)=fact*(-dratio*v*tdq(1)+dratio*tdq(3)+v*tdq(5)+(dratio+tau_vbydt)*tdq(7))

			dtbytau_t=dt(i,j)/tau_t(i,j)
			dtbytau_v=dt(i,j)/tau_v(i,j)
			a4=dratio*( (cv*temp-0.5*(u**2+v**2))*beta*dtbytau_t + (u*up+v*vp)*dtbytau_v)*rhs_g(1,i,j) &
			  +dratio*(u*beta*dtbytau_t-up*dtbytau_v)*rhs_g(2,i,j) &
			  +dratio*(v*beta*dtbytau_t-vp*dtbytau_v)*rhs_g(3,i,j) &
			  -((beta*cv*temp-0.5*(up**2+vp**2))*dtbytau_t+(up**2+vp**2)*dtbytau_v)*rhs_p(1,i,j) &
			  -(up*dtbytau_t+(u-2.*up)*dtbytau_v)*rhs_p(2,i,j) &
			  -(vp*dtbytau_t+(v-2.*vp)*dtbytau_v)*rhs_p(3,i,j)
			a8=-a4
			a4=tdq(4) + a4
			a8=tdq(8) + a8

			fact=1./(1.+beta*dratio+tau_t(i,j)/dt(i,j))
			rhs_g(4,i,j)=fact*((1.+tau_t(i,j)/dt(i,j))*a4+a8)
			rhs_p(4,i,j)=fact*(beta*dratio*a4+(beta*dratio+tau_t(i,j)/dt(i,j))*a8)
		end do
		end do
		
		return
		end