!****************************************************************************
!
!  PROGRAM:    multi_q1d_LU-SGS.f90
!
! Quasi 1-D Euler Solver, multi-phase flow
! Roe's Scheme : Flux Difference Splliting Scheme for Euler Equations
!
!  LU-SGS time integration
!  higher-order space approximation
!  Limiting at most
!  Entropy Fix is done
!
!     phi=1.  :  Higher-Oder Space discritization
!     phi=0.  :  1st-Oder Space discritization
!	  kappa=-1    : full upwind 2nd order
!	  kappa=0.    : Fromm's 2nd order
!	  kappa=1./3. : 3rd order 
!	  kappa=1.	  : central 2nd order
!
! Copyright 2001, T.Shimada 01/3/2
!
!****************************************************************************
      program multi_q1d_lusgs

      implicit real*8 (a-h,o-z)
      real*8 kappa, nu_max_g, nu_max_p
      parameter(n=201)
!
      dimension x(n),q_g(3,-1:n+1), q_p(3,-1:n+1) ,f_g(3,n), f_p(3,n) 
	  dimension w_L(6,0:n+1),w_R(6,0:n+1)
	  dimension A(n), Vol(n)
	  dimension source_g(3), source_p(3)
      dimension dw_b(6), dw_g(6), w_c(6)
      dimension rhs_g(3),rhs_p(3), dq_g$(3,n), dq_p$(3,n), abar(n)
	  dimension dfdq_g(3,3,n), dfdq_p(3,3,n)
	  dimension nu_max_g(n), nu_max_p(n), factor_g(n), factor_p(n)
	  dimension ds(3), a_g$(3,3), b_g$(3,3), c_g$(3,3)
	  dimension a_p$(3,3), b_p$(3,3), c_p$(3,3)
!
      namelist /input/ CFL, maxiter, ipiter, phi, kappa, delta, imp
      namelist /gas/ rL,uL,pL,tL,rR,uR,pR,tR, amol, gam
	  namelist /particle/ rpL, upL, tpL, rpR, upR, tpR, dp, aloadratio
	  namelist /space/ xmin, xmax
	  namelist /boundary/ ibtype_in, ibtype_out
! 
! Definition of statement functions
!
      fnp(i)=gm*(q_g(3,i)-0.5*q_g(2,i)**2/q_g(1,i))
      fnu(i)=q_g(2,i)/q_g(1,i)
      fna(i)=sqrt(gam*fnp(i)/q_g(1,i))
	  fnt(i)=(q_g(3,i)/q_g(1,i)-0.5*fnu(i)**2)/cv
	  fnarea(xx)=1.+2.2*(xx-1.5)**2
	  fnia(xx)=xx+2.2/3.*(xx-1.5)**3
	  fnup(i)=q_p(2,i)/q_p(1,i)
	  fntp(i)=(q_p(3,i)/q_p(1,i)-0.5*fnup(i)**2)/cc
	  func_g(a$)=(1.+a$*(12.278+0.548*a$))/(1.+11.278*a$)
!
! Set Default values
!
      gam=1.211
	  amol=20.33					! g / mole
      amyu=7.e-5					! Pa-s
!
!     particle phase
!
      sigmap=3204.					! kg/m^3
	  dp = 1.e-6					! m
	  cc = 1380.					! J/kg-K
	  aloadratio=0.4

! types of inflow/outflow boundaries
!     ibtype_in =11  !  subsonic,  1st-order
      ibtype_in =12  !  subsonic,  2nd-order
!	  ibtype_in =20  !  supersonic
!
!     ibtype_out=11  !  subsonic,  1st-order
	  ibtype_out=12  !  subsonic,  2nd-order
!	  ibtype_out=21  !  supersonic, 1st-order
!	  ibtype_out=22  !  supersonic, 2nd-order
!

! read input data
!
      open(99,file='.\input_Mq1d.txt',status='unknown')
      read(99,input)
      read(99,gas)
      read(99,particle)
	  read(99,space)
	  read(99,boundary)
      close(99)


	  gm=gam-1.
	  gasc=8314.511212/amol			! J/kg-K
	  cp=gam/gm*gasc
	  cv=1./gm*gasc
	  Prandtl=4.*gam/(9.*gam-5.)
	  alambda=amyu*cp/Prandtl
      beta=cc/cv

	  rL=pL/(gasc*tL)
	  rR=pR/(gasc*tR)
	  rpL=aloadratio*rL
	  rpR=aloadratio*rR

	  p0=pL
	  r0=rL
	  t0=tL
	  rp0=rpL
	  tp0=tpL
	  a0=sqrt(gam*gasc*tL)
	  ae0=a0*sqrt((gam+aloadratio*beta)/gam/(1.+aloadratio*beta)/(1.+aloadratio))

      write(6,input)
      write(6,gas)
	  write(6,particle)
	  write(6,space)
	  write(6,boundary)
      write(60,input)
      write(60,gas)
	  write(60,particle)
	  write(60,space)
	  write(60,boundary)
!
! Cell locations
!
      dx=(xmax-xmin)/float(n-1)
      do i=1, n
         x(i)=xmin+dx*float(i-1)
		 A(i)=fnarea(x(i))
      end do
!
! Set Initial State
!
      do i=1, n-1
         dr$=(rR-rL)/float(n-1)
		 du$=(uR-uL)/float(n-1)
		 dp$=(pR-pL)/float(n-1)
         drp$=(rpR-rpL)/float(n-1)
		 dup$=(upR-upL)/float(n-1)
		 dtp$=(tpR-tpL)/float(n-1)
!
		 r$=rL+dr$*float(i-1)
		 p$=pL+dp$*float(i-1)
		 u$=uL+du$*float(i-1)
		 rp$=rpL+drp$*float(i-1)
		 tp$=tpL+dtp$*float(i-1)
		 up$=upL+dup$*float(i-1)

		 q_g(1,i)=r$
		 q_g(2,i)=r$*u$
		 q_g(3,i)=p$/gm+0.5*r$*u$**2
		 q_p(1,i)=rp$
		 q_p(2,i)=rp$*up$
		 q_p(3,i)=rp$*(cc*tp$+0.5*up$**2)

		 x1=x(i)
		 x2=x(i+1)
		 vol(i)=fnia(x2)-fnia(x1)
      end do

	  time=0.
      ic=0
      ip=10
	  ipout=ipiter
      write(ip,'(e16.8)') (0.5*(x(i)+x(i+1)),i=1,n-1)
	  ip=ip+1
	  write(ip,'(6e16.8)') (q_g(1,i)/r0,fnu(i)/ae0,fnt(i)/t0,q_p(1,i)/rp0,fnup(i)/ae0,fntp(i)/tp0,i=1,n-1)

	  open(51,file='residual.dat',status='unknown')
!
! Main Loop
!
      do while (ic.le.maxiter)

!Boudary condition
        if(ibtype_in.eq.11) then
			r_in=rL
			u_in=fnu(1)
			p_in=pL
		else if(ibtype_in.eq.12) then
			r_in=rL
			u_in=2.*fnu(1)-fnu(2)
			p_in=pL
		else
		   stop 'inflow boundary type does not match'
		end if
		rp_in=r_in*aloadratio
		up_in=u_in
		tp_in=p_in/(gasc*r_in)
!
		if(ibtype_out.eq.11) then
			r_out=q_g(1,n-1)
			u_out=fnu(n-1)
			p_out=pR
		else if(ibtype_out.eq.12) then
			r_out=2.*q_g(1,n-1)-q_g(1,n-2)
			u_out=2.*fnu(n-1)-fnu(n-2)
			p_out=pR
		else if(ibtype_out.eq.21) then
			r_out=q_g(1,n-1)
			u_out=fnu(n-1)
			p_out=fnp(n-1)
		else if(ibtype_out.eq.22) then
			r_out=2.*q_g(1,n-1)-q_g(1,n-2)
			u_out=2.*fnu(n-1)-fnu(n-2)
			p_out=2.*fnp(n-1)-fnp(n-2)
		else
			stop 'outflow boundary type does not match'
		end if

		rp_out=q_p(1,n-1)
		up_out=fnup(n-1)
		tp_out=fntp(n-1)

		q_g(1,0)=r_in
		q_g(2,0)=r_in*u_in
		q_g(3,0)=p_in/gm+0.5*r_in*u_in**2
		q_p(1,0)=rp_in
		q_p(2,0)=rp_in*up_in
		q_p(3,0)=rp_in*(cc*tp_in+0.5*up_in**2)

		q_g(1,n)=r_out
		q_g(2,n)=r_out*u_out
		q_g(3,n)=p_out/gm+0.5*r_out*u_out**2
		q_p(1,n)=rp_out
		q_p(2,n)=rp_out*up_out
		q_p(3,n)=rp_out*(cc*tp_out+0.5*up_out**2)
!
! Time step evaluation
         dtmin=1.e30
         do i=1,n-1
		    dtf=(x(i+1)-x(i))/max(abs(fnu(i))+fna(i),abs(fnup(i)))
            dtmin=min(dtmin,dtf)
         end do
         dt=cfl*dtmin
         time=time+dt
         ic=ic+1
!
! evaluate left and right values at cell interfaces
!
         do i=1, n-1
            dw_b(1) = q_g(1,i  ) - q_g(1,i-1)
            dw_b(2) = fnu(i  ) - fnu(i-1)
            dw_b(3) = fnp(i  ) - fnp(i-1)
			dw_b(4) = q_p(1,i  ) - q_p(1,i-1)
			dw_b(5) = fnup(i ) - fnup(i-1)
			dw_b(6) = fntp(i ) - fntp(i-1)

            dw_g(1) = q_g(1,i+1) - q_g(1,i  )
            dw_g(2) = fnu(i+1) - fnu(i  )
            dw_g(3) = fnp(i+1) - fnp(i  )
			dw_g(4) = q_p(1,i+1) - q_p(1,i  )
			dw_g(5) = fnup(i+1)- fnup(i )
			dw_g(6) = fntp(i+1)- fntp(i )

            w_c(1) = q_g(1,i)
            w_c(2) = fnu(i)
            w_c(3) = fnp(i)
			w_c(4) = q_p(1,i)
			w_c(5) = fnup(i)
			w_c(6) = fntp(i)

            do j=1, 6
               del_w_r=0.5*((1.+kappa)*dw_b(j)+(1.-kappa)*dw_g(j))
               del_w_l=0.5*((1.-kappa)*dw_b(j)+(1.+kappa)*dw_g(j))

               del_w_lim=0.
               if(dw_b(j).ne.0.) then
                  theta=dw_g(j)/dw_b(j)
                  if(theta.gt.0.) then
                     dl=del_w_l*min(4.*theta/(1.-kappa+(1.+kappa)*theta),1.)
                     dr=del_w_r*min(4./(1.+kappa+(1.-kappa)*theta),1.)
                     del_w_lim=0.5*(sign(1.,dl)+sign(1.,dr))*min(abs(dl),abs(dr))
                  end if
               end if

               w_R(j,i  )=w_c(j)-0.5*phi*del_w_lim
               w_L(j,i+1)=w_c(j)+0.5*phi*del_w_lim

            end do
         end do
!
		w_R(1,n)=q_g(1,n)
		w_R(2,n)=fnu(n)
		w_R(3,n)=fnp(n)
		w_R(4,n)=q_p(1,n)
		w_R(5,n)=fnup(n)
		w_R(6,n)=fntp(n)

		w_L(1,1)=q_g(1,0)
		w_L(2,1)=fnu(0)
		w_L(3,1)=fnp(0)
		w_L(4,1)=q_p(1,0)
		w_L(5,1)=fnup(0)
		w_L(6,1)=fntp(0)

!
! Numerical Flux Evaluation  : Roe's Scheme
!
         do i=1, n
            HL=gam/gm*w_L(3,i)/w_L(1,i)+0.5*w_L(2,i)**2
            HR=gam/gm*w_R(3,i)/w_R(1,i)+0.5*w_R(2,i)**2

            rrl=sqrt(w_L(1,i))
            rrr=sqrt(w_R(1,i))
            fctL=rrl/(rrl+rrr)
            fctR=rrr/(rrl+rrr)

            uH=fctL*w_L(2,i) + fctR*w_R(2,i)
            HH=fctL*HL       + fctR*HR

            aH=sqrt(gm*(HH-0.5*uH**2))

            dq1=w_R(1,i)-w_L(1,i)
            dq2=w_R(1,i)*w_R(2,i)-w_L(1,i)*w_L(2,i)
            dq3=w_R(3,i)/gm+0.5*w_R(1,i)*w_R(2,i)**2 &
			   -w_L(3,i)/gm-0.5*w_L(1,i)*w_L(2,i)**2

            alpha1=(1.-0.5*gm*(uH/aH)**2)*dq1 &
                  +gm/aH**2*uH           *dq2 &
                  -gm/aH**2              *dq3
            alpha2=0.25*uH*(2.*aH+gm*uH) *dq1 &
                  +0.5*(-aH-gm*uH)       *dq2 &
                  +0.5*gm                *dq3
            alpha3=0.25*uH*(-2.*aH+gm*uH)*dq1 &
                  +0.5*(aH-gm*uH)        *dq2 &
                  +0.5*gm                *dq3

            aramda1=abs(    uH)
            aramda2=abs(-aH+uH)
            aramda3=abs( aH+uH)

! Entropy Fix Considered
			If(aramda2.lt.delta*aH) aramda2=(aramda2**2+delta**2*aH**2)/(2.*delta*aH)
			If(aramda3.lt.delta*aH) aramda3=(aramda3**2+delta**2*aH**2)/(2.*delta*aH)
!
            d1=alpha1                                *aramda1 &
			  +alpha2 /aH**2                         *aramda2 &
			  +alpha3 /aH**2                         *aramda3
            d2=alpha1 *uH                            *aramda1 &
			  +alpha2 *(-aH+uH)/aH**2                *aramda2 &
			  +alpha3 *( aH+uH)/aH**2                *aramda3
            d3=alpha1 *0.5*uH**2                     *aramda1 &
			  +alpha2 *(-uH/aH+0.5*(uH/aH)**2+1./gm) *aramda2 &
              +alpha3 *( uH/aH+0.5*(uH/aH)**2+1./gm) *aramda3

            f1L=w_L(1,i)*w_L(2,i)
            f2L=w_L(1,i)*w_L(2,i)**2+w_L(3,i)
            f3L=w_L(1,i)*w_L(2,i)*(w_L(3,i)*gam/gm/w_L(1,i)+0.5*w_L(2,i)**2)

            f1R=w_R(1,i)*w_R(2,i)
            f2R=w_R(1,i)*w_R(2,i)**2+w_R(3,i)
            f3R=w_R(1,i)*w_R(2,i)*(w_R(3,i)*gam/gm/w_R(1,i)+0.5*w_R(2,i)**2)

            f_g(1,i)=0.5*(f1L+f1R-d1) * A(i)
            f_g(2,i)=0.5*(f2L+f2R-d2) * A(i)
            f_g(3,i)=0.5*(f3L+f3R-d3) * A(i)
!
! particle phase numerical flux
!
			upL$=0.5*(w_L(5,i)+abs(w_L(5,i)))
			upR$=0.5*(w_R(5,i)-abs(w_R(5,i)))
			rpL$=w_L(4,i)
			rpR$=w_R(4,i)
			epL$=cc*w_L(6,i)+0.5*w_L(5,i)**2
			epR$=cc*w_R(6,i)+0.5*w_R(5,i)**2

            f_p(1,i)=(rpL$*upL$ + rpR$*upR$) * A(i)
            f_p(2,i)=(rpL$*upL$**2 + rpR$*upR$**2) * A(i)
            f_p(3,i)=(rpL$*upL$*epL$ + rpR$*upR$*epR$) * A(i)
         end do
!
! Source term, 
!
         do i=1, n-1
			source_g(1:3)=0.
			source_p(1:3)=0.
			source_g(2)=dt*fnp(i)*(A(i+1)-A(i))/Vol(i)
!
		    u =fnu(i)
			up=fnup(i)
		    u_rel=u-up
			reynolds=q_g(1,i)*abs(u_rel)*dp/amyu
			
			sos=sqrt(gam*fnp(i)/q_g(1,i))
			amach=abs(u_rel)/sos

			if(reynolds.gt.0.) then
               if(reynolds.lt.1.e3) then
			      CD0=24./reynolds*(1.+reynolds**(2./3.)/6.)
			   else
			      CD0=0.4392
			   end if
			   cfactor=amach/(reynolds*Prandtl)
			   aNusselt0=2.+0.459*reynolds**0.55*Prandtl**0.33
			else
			   cfactor=alambda/(q_g(1,i)*dp*sos*cp)
			   aNusselt0=2.
			end if

			aNusselt=aNusselt0/(1.+3.42*cfactor*aNusselt0)
			tau_t=cc*sigmap*dp**2/(6.*aNusselt*alambda)
			dtbytau_t=dt/tau_t


            temp_g=fnt(i)
			temp_p=fntp(i)
			stemp=temp_g-temp_p
			source_g(3)=-q_p(1,i)*cc*stemp*dtbytau_t

			tau_v_eq=4./3.*sigmap*dp**2/(24.*amyu)
			if(abs(u_rel).gt.0.) then
  			   tratio=temp_p / temp_g
			   func_h=5.6/(1.+amach)+1.7*sqrt(tratio)
               CD=2.+(CD0-2.)*Exp(-3.07*Sqrt(gam)*amach/reynolds*func_g(reynolds)) &
			     +func_h/(sqrt(gam)*amach)*exp(-0.5*reynolds/amach)
			   tau_v=4./3.*sigmap*dp**2/(reynolds*amyu*CD)
			   dtbytau_v=dt/tau_v

			   source_g(2)=source_g(2)-q_p(1,i)*u_rel*dtbytau_v
			   source_g(3)=source_g(3)-q_p(1,i)*u_rel*up*dtbytau_v
			   source_p(2)=q_p(1,i)*u_rel*dtbytau_v
			else
			   tau_v=4./3.*sigmap*dp**2/(24.*amyu)
			   dtbytau_v=dt/tau_v
			end if

		    source_p(3)=-source_g(3)

!
! Right-Hand Side 
!
			rhs_g(1:3)=-f_g(1:3,i+1)+f_g(1:3,i)
			rhs_p(1:3)=-f_p(1:3,i+1)+f_p(1:3,i)
			rhs_g(1:3)=rhs_g(1:3)+Vol(i)*source_g(1:3)/dt
			rhs_p(1:3)=rhs_p(1:3)+Vol(i)*source_p(1:3)/dt

		
! 
! LU-SGS Scheme
!
! the 1st step :  point-implicit scheme
!
			dratio=q_p(1,i)/q_g(1,i)
			umax_p=abs(fnup(i))
			umax_g=abs(fnu(i))+fna(i)
			distance=umax_p*max(tau_v,tau_t)

			if(distance.lt.dx) then
				nu_max_g(i)=max(umax_p,umax_g)
				nu_max_p(i)=nu_max_g(i)
			else
				nu_max_g(i)=umax_g
				nu_max_p(i)=umax_p
			end if

			abar(i)=0.5*(A(i)+A(i+1))
			if(imp.ne.0) then
				factor_g(i)=Vol(i)/dt+abar(i)*nu_max_g(i)
				factor_p(i)=Vol(i)/dt+abar(i)*nu_max_p(i)
			else
				factor_g(i)=Vol(i)/dt
				factor_p(i)=Vol(i)/dt
			end if

			rhs_g(1:3)=rhs_g(1:3)/factor_g(i)
			rhs_p(1:3)=rhs_p(1:3)/factor_p(i)
!
			fact=(A(i+1)-A(i))/Vol(i)*gm
			u=fnu(i)
			ds(1)= fact*0.5*u*u
			ds(2)=-fact*u
			ds(3)= fact
!
			a_g$=0.
			a_p$=0.
			a_g$(1,1) = 1.
			a_g$(2,1) = dt*ds(1)/(1.-dt*ds(2))
			a_g$(2,2) = 1./(1.-dt*ds(2))
			a_g$(2,3) = dt*ds(3)/(1.-dt*ds(2))
			a_g$(3,3) = 1.
			a_p$(1,1) = 1.
			a_p$(2,2) = 1.
			a_p$(3,3) = 1.
!
			dq_g$(1:3,i)=MATMUL(a_g$,rhs_g(1:3))
			rhs_g(1:3)=dq_g$(1:3,i)
			dq_p$(1:3,i)=rhs_p(1:3)

			temp_g=fnt(i)
			up=fnup(i)

			fact=1./(dtbytau_v+1.+dtbytau_v*dratio)
			dq_g$(2,i)=fact * (rhs_g(1) * dtbytau_v * u * dratio &
			                  +rhs_g(2) * (dtbytau_v + 1.) &
							  -rhs_p(1) * dtbytau_v * u  &
							  +rhs_p(2) * dtbytau_v  )
							  
			dq_p$(2,i)=fact *(-rhs_g(1) * dtbytau_v * u * dratio &
			                  +rhs_g(2) * dtbytau_v * dratio &
							  +rhs_p(1) * dtbytau_v * u  &
							  +rhs_p(2) * (dtbytau_v * dratio + 1.) )

			dsdq31=-((cv*temp_g - u**2/2.)*dratio*beta*dtbytau_t + u*up*dratio*dtbytau_v)
			dsdq32=-(u*beta*dratio*dtbytau_t- up*dratio*dtbytau_v)
			dsdq33=1.+beta*dratio*dtbytau_t
			dsdq34=-((+up**2/2. - cv*temp_g*beta)*dtbytau_t - up**2*dtbytau_v)
			dsdq35=-(-up*dtbytau_t + (-u + 2.*up)*dtbytau_v)
			dsdq36=-dtbytau_t

			dsdq61= (cv*temp_g - u**2/2.)*dratio*beta*dtbytau_t + u*up*dratio*dtbytau_v
			dsdq62= u*beta*dratio*dtbytau_t - up*dratio*dtbytau_v
			dsdq63= -beta*dratio*dtbytau_t
			dsdq64= (+up**2/2. - cv*temp_g*beta)*dtbytau_t - up**2*dtbytau_v
			dsdq65= -up*dtbytau_t + (-u + 2.*up)*dtbytau_v
			dsdq66=1.+dtbytau_t

			rhs_g(3)=rhs_g(3)-(dsdq31*dq_g$(1,i)+dsdq32*dq_g$(2,i)+dsdq34*dq_p$(1,i)+dsdq35*dq_p$(2,i))
			rhs_p(3)=rhs_p(3)-(dsdq61*dq_g$(1,i)+dsdq62*dq_g$(2,i)+dsdq64*dq_p$(1,i)+dsdq65*dq_p$(2,i))

			dq_g$(3,i)= (dsdq66*rhs_g(3) - dsdq36*rhs_p(3))/ (-dsdq36*dsdq63 + dsdq33*dsdq66)
			dq_p$(3,i)= (dsdq63*rhs_g(3) - dsdq33*rhs_p(3))/ ( dsdq36*dsdq63 - dsdq33*dsdq66)
		end do

		if(imp.ne.0) then
!
! the 2nd step : forward substitution
!
		do i=1, n-1
			u=fnu(i)
			H=0.5*u**2+fna(i)**2/gm
			dfdq_g(1:3,1:3,i)=0.

			dfdq_g(1,1,i)=0.
			dfdq_g(1,2,i)=1.
			dfdq_g(1,3,i)=0.
			dfdq_g(2,1,i)=0.5*(gam-3.)*u*u
			dfdq_g(2,2,i)=(3.-gam)*u
			dfdq_g(2,3,i)=gm
			dfdq_g(3,1,i)=0.5*u*(-2.*H+gm*u*u)
			dfdq_g(3,2,i)=H-gm*u*u
			dfdq_g(3,3,i)=gam*u

			dfdq_p(1:3,1:3,i)=0.
			up=fnup(i)
			ep=q_p(3,i)/q_p(1,i)
			dfdq_p(1,1,i)=0.
			dfdq_p(1,2,i)=1.
			dfdq_p(1,3,i)=0.
			dfdq_p(2,1,i)=-up**2
			dfdq_p(2,2,i)=2.*up
			dfdq_p(2,3,i)=0.
			dfdq_p(3,1,i)=-ep*up
			dfdq_p(3,2,i)=ep
			dfdq_p(3,3,i)=up
		end do
!
		do i=2, n-1
			b_g$ = 0.5*dfdq_g(1:3,1:3,i-1)
			b_g$(1,1)=b_g$(1,1)+0.5*nu_max_g(i-1)
			b_g$(2,2)=b_g$(2,2)+0.5*nu_max_g(i-1)
			b_g$(3,3)=b_g$(3,3)+0.5*nu_max_g(i-1)
			b_p$ = 0.5*dfdq_p(1:3,1:3,i-1)
			b_p$(1,1)=b_p$(1,1)+0.5*nu_max_p(i-1)
			b_p$(2,2)=b_p$(2,2)+0.5*nu_max_p(i-1)
			b_p$(3,3)=b_p$(3,3)+0.5*nu_max_p(i-1)
!
			fact_g=abar(i-1)/factor_g(i)
			fact_p=abar(i-1)/factor_p(i)
			dq_g$(1:3,i)=dq_g$(1:3,i)+fact_g*MATMUL(b_g$,dq_g$(1:3,i-1))
			dq_p$(1:3,i)=dq_p$(1:3,i)+fact_p*MATMUL(b_p$,dq_p$(1:3,i-1))
		end do
!
! the 3rd step : backword substitution
!
		do i=n-2, 1, -1
			c_g$ = 0.5*dfdq_g(1:3,1:3,i+1)
			c_g$(1,1)=c_g$(1,1)-0.5*nu_max_g(i+1)
			c_g$(2,2)=c_g$(2,2)-0.5*nu_max_g(i+1)
			c_g$(3,3)=c_g$(3,3)-0.5*nu_max_g(i+1)

			c_p$ = 0.5*dfdq_p(1:3,1:3,i+1)
			c_p$(1,1)=c_p$(1,1)-0.5*nu_max_p(i+1)
			c_p$(2,2)=c_p$(2,2)-0.5*nu_max_p(i+1)
			c_p$(3,3)=c_p$(3,3)-0.5*nu_max_p(i+1)
!
			fact_g=abar(i+1)/factor_g(i)
			fact_p=abar(i+1)/factor_p(i)
			dq_g$(1:3,i)=dq_g$(1:3,i)-fact_g*MATMUL(c_g$,dq_g$(1:3,i+1))
			dq_p$(1:3,i)=dq_p$(1:3,i)-fact_p*MATMUL(c_p$,dq_p$(1:3,i+1))
		end do

		end if
			
!
! Update field
!
		 residual=0.

		 fact=1.
         do i=1, n-1
			fact_e=abs(q_p(1,i)*cc*fntp(i)) / (abs((0.5*fnup(i)**2-cc*fntp(i))*dq_p$(1,i)) &
				+abs(fnup(i)*dq_p$(2,i))+abs(dq_p$(3,i)))
			fact_r=abs(q_p(1,i))/ (abs(dq_p$(1,i))+1.e-30)
			fact=min(fact,fact_e,fact_r)
		 end do

		imin=n
		tpmin=fntp(0)
         do i=1, n-1
			eng_old=q_g(3,i)
			q_g(1:3,i)=q_g(1:3,i)+dq_g$(1:3,i) * fact
			q_p(1:3,i)=q_p(1:3,i)+dq_p$(1:3,i) * fact

			if(tpmin.gt.fntp(i)) then
				imin=i
				tpmin=fntp(i)
			end if
			residual=residual+(q_g(3,i)-eng_old)**2
         end do
		 write(6,*) 'imin=',imin, '  tpmin=',tpmin
		 residual=sqrt(residual)
!
! output
!
         write(51,'(i6,e16.8)') ic,log10(residual+1.e-20)

         if(ic.ge.ipout) then
            ipout=ipout+ipiter
            ip=ip+1
            write(6,*) 'Time=',time,'    residual=',log10(residual+1.e-20)
			
    		write(ip,'(6e16.8)') (q_g(1,i)/r0,fnu(i)/ae0,fnt(i)/t0,q_p(1,i)/rp0,fnup(i)/ae0,fntp(i)/tp0,i=1,n-1)
         end if
      end do

      stop
	end