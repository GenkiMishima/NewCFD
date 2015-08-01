		subroutine time_step
!
		include 'common.h'
		include 'function.h'
!
		do j=1, nj-1
			do i=1, ni-1
				if(iflag(i,j).eq.1) then
					sos=fna(i,j)
					u = fnu(i,j)
					v = fnv(i,j)
					u_p = fnu_p(i,j)
					v_p = fnv_p(i,j)
!
					un_i_1 = abs( norm_i(1,  i,  j)* u + norm_i(2,  i,  j) * v ) + sos 
					un_i_2 = abs( norm_i(1,i+1,  j)* u + norm_i(2,i+1,  j) * v ) + sos
					un_j_1 = abs( norm_j(1,  i,  j)* u + norm_j(2,  i,  j) * v ) + sos
					un_j_2 = abs( norm_j(1,  i,j+1)* u + norm_j(2,  i,j+1) * v ) + sos

					un_p_i_1 = abs( norm_i(1,  i,  j)* u_p + norm_i(2,  i,  j) * v_p ) 
					un_p_i_2 = abs( norm_i(1,i+1,  j)* u_p + norm_i(2,i+1,  j) * v_p )
					un_p_j_1 = abs( norm_j(1,  i,  j)* u_p + norm_j(2,  i,  j) * v_p )
					un_p_j_2 = abs( norm_j(1,  i,j+1)* u_p + norm_j(2,  i,j+1) * v_p )

!					unpi1=0.5*(un_p_i_1+abs(un_p_i_1))
!					unpi2=0.5*(un_p_i_2+abs(un_p_i_2))
!					unpj1=0.5*(un_p_j_1+abs(un_p_j_1))
!					unpj2=0.5*(un_p_j_2+abs(un_p_j_2))
!
					dt_gas = vol(i,j) / (ds_i(i,j)*rave_i(i,j)*un_i_1 + ds_i(i+1,j)*rave_i(i+1,j)*un_i_2 &
						    		    +ds_j(i,j)*rave_j(i,j)*un_j_1 + ds_j(i,j+1)*rave_j(i,j+1)*un_j_2)
					dt_p   = vol(i,j) / (ds_i(i,j)*rave_i(i,j)*un_p_i_1 + ds_i(i+1,j)*rave_i(i+1,j)*un_p_i_2 &
						    		    +ds_j(i,j)*rave_j(i,j)*un_p_j_1 + ds_j(i,j+1)*rave_j(i,j+1)*un_p_j_2)
!                    dt_p= vol(i,j) / &
!					    (ds_i(i,j)*rave_i(i,j)*unpi1 + ds_i(i+1,j)*rave_i(i+1,j)*unpi2 &
!						+ds_j(i,j)*rave_j(i,j)*unpj1 + ds_j(i,j+1)*rave_j(i,j+1)*unpj2) 
					dt(i,j) = cfl * min(dt_gas, dt_p)
				end if
			end do
		end do
!
		dtmin=1.e30
		if(ilocal.eq.0) then ! global time step is employed
			do j=1, nj-1
				do i=1, ni-1
					if(iflag(i,j).eq.1) then
						dtmin=min(dtmin,dt(i,j))
					end if
				end do
			end do
			do j=1, nj-1
				do i=1, ni-1
					dt(i,j)=dtmin
				end do
			end do
		end if
!
		return
		end