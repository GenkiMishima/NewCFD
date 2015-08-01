 		subroutine update
!
		include 'common.h'
		real rtime, ta(2)
		include 'function.h'
!

		 fact=1.
		 do j=1, nj-1
         do i=1, ni-1
			up=fnu_p(i,j)
			vp=fnv_p(i,j)
			tp=fnt_p(i,j)
			rp=fnr_p(i,j)
			fact_e=abs(rp*cc*tp) / ( &
			       abs( (0.5*(up**2+vp**2)-cc*tp)*rhs_p(1,i,j) ) &
			     + abs( up*rhs_p(2,i,j) ) &
				 + abs( vp*rhs_p(3,i,j) ) &
				 + abs( rhs_p(4,i,j) )  )
			fact_r=abs(rp)/ (abs(rhs_p(1,i,j))+1.e-30)
			if(tp.eq.0.) write(6,*) i,j,q_p(1:4,i,j)
			fact=min(fact,fact_e,fact_r)
		 end do
		 end do

		residual=0.
		ic=ic+1
		do j=1, nj-1
			do i=1, ni-1
				if(iflag(i,j).eq.1) then
					eng_old=q_g(4,i,j)
					ep_old=q_p(4,i,j)
					q_g(1:4,i,j)=q_g(1:4,i,j)+rhs_g(1:4,i,j) * fact
					q_p(1:4,i,j)=q_p(1:4,i,j)+rhs_p(1:4,i,j) * fact
!					residual=residual+(q_g(4,i,j)-eng_old)**2
					residual=residual+(q_g(4,i,j)-eng_old)**2+(q_p(4,i,j)-ep_old)**2
				end if
			end do
		end do
		residual = sqrt(residual)
!
		do j=1, nj-1
		do i=1, ni-1
			if(fnp(i,j).le.0.) write(6,*) i,j,q_p(1:4,i,j)
		end do
		end do
		return
		end