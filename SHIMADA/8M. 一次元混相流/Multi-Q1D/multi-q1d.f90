!****************************************************************************
!
!  PROGRAM: multi_q1d.f90
!
! Quasi 1-D Multi-phase flow Solver
! Roe's Scheme for gas phase
! FVS Scheme for particle phase
!
!  first-order explicit time integration
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
! T.Shimada 01/1/25
!
!****************************************************************************
      program multi_q1d

      implicit real*8 (a-h,o-z)
      real*8 kappa
      parameter(n=201)
!
      dimension x(n),q(6,0:n),f(6,n),w_L(6,n),w_R(6,n), A(n), Vol(n), &
				source(6), dq(6,n), rhs(6)
      dimension dw_b(6), dw_f(6), w_c(6)
	  dimension dfdq(6,6)
      namelist /input/ CFL, maxiter, ipiter, phi, kappa, delta, imp
      namelist /gas/ rL,uL,pL,tL,rR,uR,pR,tR, amol, gam
	  namelist /particle/ rpL,upL,tpL,rpR,upR,tpR, dp, aloadratio
	  namelist /space/ xmin, xmax
	  namelist /boundary/ ibtype_in, ibtype_out

! 
! Definition of statement functions
!
      fnp(i)=gm*(q(3,i)-0.5*q(2,i)**2/q(1,i))
      fnu(i)=q(2,i)/q(1,i)
      fna(i)=sqrt(gam*fnp(i)/q(1,i))
	  fnt(i)=(q(3,i)/q(1,i)-0.5*fnu(i)**2)/cv

	  fnarea(xx)=1.+2.2*(xx-1.5)**2
	  fnia(xx)=xx+2.2/3.*(xx-1.5)**3

	  fnup(i)=q(5,i)/q(4,i)
	  fntp(i)=(q(6,i)/q(4,i)-0.5*fnup(i)**2)/cc
	  func_g(a$)=(1.+a$*(12.278+0.548*a$))/(1.+11.278*a$)

! Set Default values
!!
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

		 q(1,i)=r$
		 q(2,i)=r$*u$
		 q(3,i)=p$/gm+0.5*r$*u$**2
		 q(4,i)=rp$
		 q(5,i)=rp$*up$
		 q(6,i)=rp$*(cc*tp$+0.5*up$**2)

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
	  write(ip,'(6e16.8)') (q(1,i)/r0,fnu(i)/ae0,fnt(i)/t0,q(4,i)/rp0,fnup(i)/ae0,fntp(i)/tp0,i=1,n-1)

	  open(51,file='residual.dat',status='unknown')
!
! Main Loop
!
      do while (ic.le.maxiter)

!Boudary condition
        if(ibtype_in.eq.11) then
		   u_in=fnu(1)
		   r_in=rL
		   p_in=pL
		else if(ibtype_in.eq.12) then
		   u_in=2.0*fnu(1)-fnu(2)
		   r_in=rL
		   p_in=pL
		else if(ibtype_in.eq.20) then
		   u_in=uL
		   r_in=rL
		   p_in=pL
		else
		   stop 'inflow boundary type does not match'
		end if
!
		rp_in=r_in*aloadratio
		up_in=u_in
		tp_in=p_in/(gasc*r_in)
!
		if(ibtype_out.eq.11) then
		   u_out=fnu(n-1)
		   r_out=q(1,n-1)
		   p_out=pR
		else if(ibtype_out.eq.12) then
		   u_out=2.0*fnu(n-1)-fnu(n-2)
		   r_out=2.0*q(1,n-1)-q(1,n-2)
		   p_out=pR
		else if(ibtype_out.eq.21) then
		   u_out=fnu(n-1)
		   r_out=q(1,n-1)
		   p_out=fnp(n-1)
		else if(ibtype_out.eq.22) then
		   u_out=2.0*fnu(n-1)-fnu(n-2)
		   r_out=2.0*q(1,n-1)-q(1,n-2)
		   p_out=2.0*fnp(n-1)-fnp(n-2)
		else
			stop 'outflow boundary type does not match'
		end if

		rp_out=q(4,n-1)
		up_out=fnup(n-1)
		tp_out=fntp(n-1)
!
!		rp_out=1.5*q(4,n-1) -0.5*q(4,n-2)
!		up_out=1.5*fnup(n-1)-0.5*fnup(n-2)
!		tp_out=1.5*fntp(n-1)-0.5*fntp(n-2)
!
		q(1,0)=r_in
		q(2,0)=r_in*u_in
		q(3,0)=p_in/gm+0.5*r_in*u_in**2
		q(4,0)=rp_in
		q(5,0)=rp_in*up_in
		q(6,0)=rp_in*(cc*tp_in+0.5*up_in**2)
!
		q(1,n)=r_out
		q(2,n)=r_out*u_out
		q(3,n)=p_out/gm+0.5*r_out*u_out**2
		q(4,n)=rp_out
		q(5,n)=rp_out*up_out
		q(6,n)=rp_out*(cc*tp_out+0.5*up_out**2)
!
! Time step evaluation
         dtmin=1.e30
         do i=1,n-1
		    dti=(x(i+1)-x(i))/max(abs(fnu(i))+fna(i),abs(fnup(i)))
            dtmin=min(dtmin,dti)
         end do
         dt=cfl*dtmin

         time=time+dt
         ic=ic+1
!
! evaluate left and right values at cell interfaces
!
         do i=1, n-1
            dw_b(1) = q(1,i  ) - q(1,i-1)
            dw_b(2) = fnu(i  ) - fnu(i-1)
            dw_b(3) = fnp(i  ) - fnp(i-1)
			dw_b(4) = q(4,i  ) - q(4,i-1)
			dw_b(5) = fnup(i ) - fnup(i-1)
			dw_b(6) = fntp(i ) - fntp(i-1)

            dw_f(1) = q(1,i+1) - q(1,i  )
            dw_f(2) = fnu(i+1) - fnu(i  )
            dw_f(3) = fnp(i+1) - fnp(i  )
			dw_f(4) = q(4,i+1) - q(4,i  )
			dw_f(5) = fnup(i+1)- fnup(i )
			dw_f(6) = fntp(i+1)- fntp(i )

            w_c(1) = q(1,i)
            w_c(2) = fnu(i)
            w_c(3) = fnp(i)
			w_c(4) = q(4,i)
			w_c(5) = fnup(i)
			w_c(6) = fntp(i)

            do j=1, 6
               del_w_r=0.5*((1.+kappa)*dw_b(j)+(1.-kappa)*dw_f(j))
               del_w_l=0.5*((1.-kappa)*dw_b(j)+(1.+kappa)*dw_f(j))

               del_w_lim=0.
               if(dw_b(j).ne.0.) then
                  theta=dw_f(j)/dw_b(j)
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
! boundary conditions
!
         w_R(1,n)=q(1,n)
         w_R(2,n)=fnu(n)
         w_R(3,n)=fnp(n)
		 w_R(4,n)=q(4,n)
		 w_R(5,n)=fnup(n)
		 w_R(6,n)=fntp(n)

         w_L(1,1)=q(1,0)
         w_L(2,1)=fnu(0)
         w_L(3,1)=fnp(0)
		 w_L(4,1)=q(4,0)
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

            f(1,i)=0.5*(f1L+f1R-d1) * A(i)
            f(2,i)=0.5*(f2L+f2R-d2) * A(i)
            f(3,i)=0.5*(f3L+f3R-d3) * A(i)
!!
! particle phase numerical flux
!
			upL$=0.5*(w_L(5,i)+abs(w_L(5,i)))
			upR$=0.5*(w_R(5,i)-abs(w_R(5,i)))
			rpL$=w_L(4,i)
			rpR$=w_R(4,i)
			epL$=cc*w_L(6,i)+0.5*w_L(5,i)**2
			epR$=cc*w_R(6,i)+0.5*w_R(5,i)**2

            f(4,i)=(rpL$*upL$ + rpR$*upR$) * A(i)
            f(5,i)=(rpL$*upL$**2 + rpR$*upR$**2) * A(i)
            f(6,i)=(rpL$*upL$*epL$ + rpR$*upR$*epR$) * A(i)
         end do
!
! Source term, point implicit time integration
!
		 residual=0.

         do i=1, n-1
		    source(1:6)=0.
			source(2)=dt*fnp(i)*(A(i+1)-A(i))/Vol(i)
!
		    u =fnu(i)
			up=fnup(i)
		    u_rel=u-up
			reynolds=q(1,i)*abs(u_rel)*dp/amyu
			
			sos=sqrt(gam*fnp(i)/q(1,i))
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
			   cfactor=alambda/(q(1,i)*dp*sos*cp)
			   aNusselt0=2.
			end if

			aNusselt=aNusselt0/(1.+3.42*cfactor*aNusselt0)
			tau_t=cc*sigmap*dp**2/(6.*aNusselt*alambda)
			dtbytau_t=dt/tau_t

            temp_g=fnt(i)
			temp_p=fntp(i)
			stemp=temp_g-temp_p
			source(3)=-q(4,i)*cc*stemp*dtbytau_t

			if(u_rel.ne.0.) then
   			   tratio=temp_p / temp_g
			   func_h=5.6/(1.+amach)+1.7*sqrt(tratio)
               CD=2.+(CD0-2.)*Exp(-3.07*Sqrt(gam)*amach/reynolds*func_g(reynolds)) &
			     +func_h/(sqrt(gam)*amach)*exp(-0.5*reynolds/amach)
			   tau_v=4./3.*sigmap*dp**2/(reynolds*amyu*CD)
			   dtbytau_v=dt/tau_v

			   source(2)=source(2)-q(4,i)*u_rel*dtbytau_v
			   source(3)=source(3)-q(4,i)*u_rel*up*dtbytau_v
			   source(5)=q(4,i)*u_rel*dtbytau_v
			else
			   tau_v=4./3.*sigmap*dp**2/(24.*amyu)
			   dtbytau_v=dt/tau_v
			end if

 		    source(6)=-source(3)
!
! point implicit
!
			rhs(1:6)=-dt/vol(i)*(f(1:6,i+1)-f(1:6,i))+source(1:6)
			dq(1:6,i)=rhs(1:6)
			fact=(A(i+1)-A(i))*gm/vol(i)

			dq(2,i)=(dt*fact*u**2*rhs(1)+2.*rhs(2)+2.*dt*fact*rhs(3))/(1.+dt*fact*u)/2.
			rhs(2)=dq(2,i)

			dratio=q(4,i)/q(1,i)
			dq(2,i)=(rhs(2)*(dtbytau_v+1.)+dtbytau_v*(rhs(5)-rhs(4)*u+rhs(1)*u*dratio)) / &
			        (dtbytau_v+1.+dtbytau_v*dratio)
			dq(5,i)=(rhs(4)*u*dtbytau_v+rhs(5)+(rhs(2)+rhs(5)-rhs(1)*u)*dtbytau_v*dratio) / &
			        (dtbytau_v+1.+dtbytau_v*dratio)
			dsdq31=-(((cv*temp_g - u**2/2.)*dratio*beta)*dtbytau_t + (u*up*dratio)*dtbytau_v)
			dsdq32=-((u*beta*dratio)*dtbytau_t - (up*dratio)*dtbytau_v)
			dsdq33=1.-(-(beta*dratio)*dtbytau_t)
			dsdq34=-(-((-up**2/2. + cv*temp_g*beta)*dtbytau_t) - up**2*dtbytau_v)
			dsdq35=-(-(up*dtbytau_t) + (-u + 2*up)*dtbytau_v)
			dsdq36=-dtbytau_t

			dsdq61= (((cv*temp_g - u**2/2.)*dratio*beta)*dtbytau_t + (u*up*dratio)*dtbytau_v)
			dsdq62= ((u*beta*dratio)*dtbytau_t - (up*dratio)*dtbytau_v)
			dsdq63= (-(beta*dratio)*dtbytau_t)
			dsdq64= (-((-up**2/2. + cv*temp_g*beta)*dtbytau_t) - up**2*dtbytau_v)
			dsdq65= (-(up*dtbytau_t) + (-u + 2*up)*dtbytau_v)
			dsdq66=1.+dtbytau_t

			rhs(3)=rhs(3)-(dsdq31*dq(1,i)+dsdq32*dq(2,i)+dsdq34*dq(4,i)+dsdq35*dq(5,i))
			rhs(6)=rhs(6)-(dsdq61*dq(1,i)+dsdq62*dq(2,i)+dsdq64*dq(4,i)+dsdq65*dq(5,i))

			dq(3,i)= (dsdq66*rhs(3) - dsdq36*rhs(6))/ (-dsdq36*dsdq63 + dsdq33*dsdq66)
			dq(6,i)= (dsdq63*rhs(3) - dsdq33*rhs(6))/ ( dsdq36*dsdq63 - dsdq33*dsdq66)
		end do
!
! update field
!
		 fact=1.
         do i=1, n-1
			fact_e=abs(q(4,i)*cc*fntp(i)) / (abs((0.5*fnup(i)**2-cc*fntp(i))*dq(4,i)) &
				+abs(fnup(i)*dq(5,i))+abs(dq(6,i)))
			fact_r=abs(q(4,i))/ (abs(dq(4,i))+1.e-30)
			fact=min(fact,fact_e,fact_r)
		 end do

		imin=n
		tpmin=fntp(0)
         do i=1, n-1
			eng_old=q(3,i)
			q(1:6,i)=q(1:6,i)+dq(1:6,i) * fact
			if(tpmin.gt.fntp(i)) then
				imin=i
				tpmin=fntp(i)
			end if
			residual=residual+(q(3,i)-eng_old)**2
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
			write(ip,'(6e16.8)') (q(1,i)/r0,fnu(i)/ae0,fnt(i)/t0,q(4,i)/rp0,fnup(i)/ae0,fntp(i)/tp0,i=1,n-1)
         end if
      end do

      stop
	end program multi_q1d