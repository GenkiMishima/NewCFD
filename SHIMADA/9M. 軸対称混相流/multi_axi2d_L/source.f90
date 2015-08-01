		subroutine source
		include 'common.h'
		dimension source_g(4), source_p(4)
		include 'function.h'
	    func_g(a$)=(1.+a$*(12.278+0.548*a$))/(1.+11.278*a$)
!
		do j=1, nj-1
		do i=1, ni-1
			rad=vol(i,j) / area(i,j)
			source_g=0. ; source_p=0.
			source_g(3)= fnp(i,j) / rad
!
			rho = fnr(i,j)
		    u =fnu(i,j)
			v =fnv(i,j)

			rhop = fnr_p(i,j)
			up=fnu_p(i,j)
			vp=fnv_p(i,j)

		    u_rel=u-up
			v_rel=v-vp
			rvelo=sqrt(u_rel**2+v_rel**2)
			reynolds=rho*rvelo*dp/amyu
			
			sos=sqrt(gam*fnp(i,j)/rho)
			amach=rvelo/sos

			if(reynolds.gt.0.) then
               if(reynolds.lt.1.e3) then
			      CD0=24./reynolds*(1.+reynolds**(2./3.)/6.)
			   else
			      CD0=0.4392
			   end if
			   cfactor=amach/(reynolds*Prandtl)
			   aNusselt0=2.+0.459*reynolds**0.55*Prandtl**0.33
			else
			   cfactor=alambda/(rho*dp*sos*cp)
			   aNusselt0=2.
			end if

			aNusselt=aNusselt0/(1.+3.42*cfactor*aNusselt0)
			tau_t(i,j)=cc*sigmap*dp**2/(6.*aNusselt*alambda)

            temp_g=fnt(i,j)
			temp_p=fnt_p(i,j)
			stemp=temp_g-temp_p
			source_g(4)=-rhop*cc*stemp/tau_t(i,j)

			if(rvelo.ne.0.) then
  			   tratio=temp_p / temp_g
			   func_h=5.6/(1.+amach)+1.7*sqrt(tratio)
               CD=2.+(CD0-2.)*Exp(-3.07*Sqrt(gam)*amach/reynolds*func_g(reynolds)) &
			     +func_h/(sqrt(gam)*amach)*exp(-0.5*reynolds/amach)
			   tau_v(i,j)=4./3.*sigmap*dp**2/(reynolds*amyu*CD)

			   source_g(2)=           -rhop*u_rel/tau_v(i,j)
			   source_g(3)=source_g(3)-rhop*v_rel/tau_v(i,j)
			   source_g(4)=source_g(4)-rhop*(u_rel*up+v_rel*vp)/tau_v(i,j)
			   source_p(2)=            rhop*u_rel/tau_v(i,j)
			   source_p(3)=            rhop*v_rel/tau_v(i,j)
			else
			   tau_v(i,j)=4./3.*sigmap*dp**2/(24.*amyu)
			end if

		    source_p(4)=-source_g(4)

!
! Right-Hand Side 
!
			rhs_g(1:4,i,j)=rhs_g(1:4,i,j)+Vol(i,j)*source_g(1:4)
			rhs_p(1:4,i,j)=rhs_p(1:4,i,j)+Vol(i,j)*source_p(1:4)

			rhs_g(1:4,i,j)=rhs_g(1:4,i,j) * dt(i,j)/vol(i,j)
			rhs_p(1:4,i,j)=rhs_p(1:4,i,j) * dt(i,j)/vol(i,j)

		end do
		end do
!
		return
		end
