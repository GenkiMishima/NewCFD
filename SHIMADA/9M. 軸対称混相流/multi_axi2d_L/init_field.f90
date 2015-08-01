		subroutine init_field
!
		include 'common.h'
		dimension ar(ni)
!
		i1d=1
!
		if(itype.eq.0) then
			do j=0, nj
			do i=0, ni
				q_g(1,i,j)=rho_amb
				q_g(2,i,j)=rho_amb * u_amb
				q_g(3,i,j)=rho_amb * v_amb
				q_g(4,i,j)=pres_amb/gm + 0.5*rho_amb*(u_amb**2+v_amb**2)
				q_p(1,i,j)=rho_p_amb
				q_p(2,i,j)=rho_p_amb * u_p_amb
				q_p(3,i,j)=rho_p_amb * u_p_amb
				q_p(4,i,j)=rho_p_amb * (cc * t_p_amb + 0.5*(u_p_amb**2+v_p_amb**2))
			end do
			end do
!
			if(i1d.eq.1) then

			gm_eff=gam_eff-1.

			do i=1, ni
				ar(i)=0.
				do j=1, nj-1	
					ar(i)=ar(i)+ds_i(i,j)*rave_i(i,j)
				end do
			end do
!
			at=ar(1)
			do i=1, ni
				if(at.gt.ar(i)) then
					it=i
					at=ar(i)
				end if
			end do
!
			do i=1, ni
				if(i.lt.it) then
					am=0.001
				else if(i.eq.it) then
					am=1.0
					goto 999
				else
					am=1.1
				end if
!
				do k=1, 500
					ff=((2.+am**2*gm_eff)/(gam_eff+1.))**((gam_eff+1.)/gm_eff)/am**2-(ar(i)/at)**2
					df=4.**(gam_eff/gm_eff)*(am**2-1.)*((2.+gm_eff*am**2)/(2.+2.*gam_eff))**(2./gm_eff) &
						/ (am**3*(gam_eff+1.))
					am=am-ff/df
				end do
999             continue

!
				pres=pres_chamber/(1.+.5*gm_eff*am**2)**(gam_eff/gm_eff)
				temp=t_chamber/(1.+.5*gm_eff*am**2)
				sos=sqrt(gam_eff*gasc_eff*temp)
				velo=sos*am
				rho=pres/(gasc_eff*temp)

				rho_p=rho * aloadratio
				temp_p = temp

				do j=1, nj-1
					vx=norm_i(1,i,j)*velo
					vy=norm_i(2,i,j)*velo
					vx_p=vx
					vy_p=vy
!
					q_g(1,i,j)=rho
					q_g(2,i,j)=rho*vx
					q_g(3,i,j)=rho*vy
					q_g(4,i,j)=pres/gm_eff+0.5*rho*(vx**2+vy**2)
					q_p(1,i,j)=rho_p
					q_p(2,i,j)=rho_p * vx_p
					q_p(3,i,j)=rho_p * vy_p
					q_p(4,i,j)=rho_p * (cc*temp_p + 0.5*(vx_p**2+vy_p**2))
				end do
			end do
			end if
!
			i=0
			q_g(1:4,i,1:nj-1)=q_g(1:4,i+1,1:nj-1)
			q_p(1:4,i,1:nj-1)=q_p(1:4,i+1,1:nj-1)

		else
			read(50) time, ic
			read(50) q_g
			read(50) q_p
		end if
!
		return
		end