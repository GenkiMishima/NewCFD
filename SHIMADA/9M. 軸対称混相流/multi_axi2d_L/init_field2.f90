		subroutine init_field
!
		include 'common.h'
		dimension ar(ni)
!
		if(itype.eq.0) then
			do j=0, nj
			do i=0, ni
				q(1,i,j)=rho_amb
				q(2,i,j)=rho_amb * u_amb
				q(3,i,j)=rho_amb * v_amb
				q(4,i,j)=pres_amb/gm + 0.5*rho_amb*(u_amb**2+v_amb**2)
			end do
			end do
!
			i=0
			q(1:4,i,1:nj-1)=q(1:4,i+1,1:nj-1)

		else
			read(50) time, ic
			read(50) q
		end if
!
		return
		end