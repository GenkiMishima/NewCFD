		subroutine I_boundary
!
		include 'common.h'
		include 'function.h'

! i = 0   :  subsonic inlet
		i=0
		un = un_in
		ut = 0.

		do j=1, nj-1
			u=norm_i(1,1,j)*un-norm_i(2,1,j)*ut
			v=norm_i(2,1,j)*un+norm_i(1,1,j)*ut
			q_g(1,i,j) = rho_chamber
			q_g(2,i,j) = rho_chamber * u
			q_g(3,i,j) = rho_chamber * v
			q_g(4,i,j) = pres_chamber /(gam_eff-1.) + 0.5*rho_chamber*(u**2+v**2)
			q_p(1,i,j) = rho_p_chamber
			q_p(2,i,j) = rho_p_chamber * u
			q_p(3,i,j) = rho_p_chamber * v
			q_p(4,i,j) = rho_p_chamber * (cc*t_p_chamber+0.5*(u**2+v**2))
		end do
!
! i=ni : out flow
		i=ni
		do j=1, nj-1
			ucheck=norm_i(1,i,j)*fnu(i-1,j)+norm_i(2,i,j)*fnv(i-1,j)
			sos =fna(i-1,j)
			amach=ucheck/sos
			if(amach.le.0.) then
				q_g(1,i,j)=rho_amb
				q_g(2,i,j)=0.
				q_g(3,i,j)=0.
				q_g(4,i,j)=pres_amb/gm
			else if(amach.ge.1.) then
				q_g(1:4,i,j)=q_g(1:4,i-1,j)
			else
				rho = fnr(i-1,j)
				xvel= fnu(i-1,j)
				yvel= fnv(i-1,j)
				pres= pres_amb
				q_g(1,i,j)=rho
				q_g(2,i,j)=rho*xvel
				q_g(3,i,j)=rho*yvel
				q_g(4,i,j)=pres/gm+0.5*rho*(xvel**2+yvel**2)
			end if
!
			u_p_check=norm_i(1,i,j)*fnu_p(i-1,j)+norm_i(2,i,j)*fnv_p(i-1,j)
			if(u_p_check.le.0.) then
				q_p(1,i,j)=rho_p_amb
				q_p(2,i,j)=0.
				q_p(3,i,j)=0.
				q_p(4,i,j)=rho_p_amb*cc*t_p_amb
			else
				q_p(1:4,i,j)=q_p(1:4,i-1,j)
			end if
		end do		
!
		return
		end
