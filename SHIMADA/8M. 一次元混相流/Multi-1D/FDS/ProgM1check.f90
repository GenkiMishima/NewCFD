!****************************************************************************
!
!  PROGRAM: prog-M1
!
! 1-D Multi-Phase Flow Solver
! Scheme : Flux Difference Splliting Scheme
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
! T.Shimada 03/07/01
!
!****************************************************************************
      program progM1

      implicit real*8 (a-h,o-z)
      real*8 kappa
!      implicit real*4 (a-h,o-z)
!      real*4 kappa
      parameter(n=501)
!
      dimension x(n),q(6,0:n),f(6,n),w_L(6,n),w_R(6,n)
	  dimension source(6,n), rhs(6), dq(6)
      dimension dw_b(6), dw_f(6), w_c(6)
      namelist /input/ CFL, tend, toutput, phi, kappa, delta, delta_p, imp
      namelist /flow/ rL,uL,pL,rR,uR,pR
	  namelist /particle/ rpL,upL,tpL,rpR,upR,tpR
	  namelist /space/ xmin, xmax, xdia
! 
! Definition of statement functions
!
      fnp(i)=gm*(q(3,i)-0.5*q(2,i)**2/q(1,i))
      fnu(i)=q(2,i)/q(1,i)
      fna(i)=sqrt(gam*fnp(i)/q(1,i))
	  fnt(i)=fnp(i)/(q(1,i)*gasc)
!	  fnM(i)=fnu(i)/fna(i)
	  fnup(i)=q(5,i)/q(4,i)
	  fntp(i)=(q(6,i)/q(4,i)-0.5*fnup(i)**2)/cc
	  func_g(a$)=(1.+a$*(12.278+0.548*a$))/(1.+11.278*a$)

!
! Set constants
!
!     gas phase
!
      gam=1.211
      gm=gam-1.
	  amol=20.33					! g / mole
	  gasc=8314.511212/amol			! J/kg-K
	  cp=gam/gm*gasc
	  cv=1./gm*gasc
	  Prandtl=4.*gam/(9.*gam-5.)
      amyu=7.e-5					! Pa-s
	  alambda=amyu*cp/Prandtl
!
!     particle phase
!
      sigmap=3204.					! kg/m^3
	  dp = 1.e-6					! m
	  cc = 1380.					! J/kg-K

      akappa=cc/cv

! Set Default values
!
      tend=0.5
      toutput=0.1
      cfl=0.1
      kappa=1./3.
      phi=1. 

      rL=5.
      uL=0.
      pL=5.

      rR=1.
      uR=0.
      pR=1.

	  rpL=rL*0.1
	  upL=uL
	  tpL=pL/(gasc*rL)

	  rpR=rR*0.1
	  upR=uR
	  tpR=pR/(gasc*rR)

	  xmin=-1.
	  xmax=1.
	  xdia=0.
!
! read input data
!
      open(99,file='input_M1.txt',status='unknown')
      read(99,input)
      read(99,flow)
	  read(99,particle)
	  read(99,space)
      close(99)

      tpL=pL/(gasc*rL)
	  tpR=pR/(gasc*rR)

      write(6,input)
      write(6,flow)
	  write(6,particle)
	  write(6,space)
      write(60,input)
      write(60,flow)
	  write(60,particle)
	  write(60,space)


!
! Cell locations
!
      dx=(xmax-xmin)/float(n-1)
      do i=1, n
         x(i)=xmin+dx*float(i-1)
      end do
!
! Set Initial State
!
      do i=1, n-1
         xh=0.5*(x(i)+x(i+1))
         if(xh.le.xdia) then
            q(1,i)=rL
            q(2,i)=rL*uL
            q(3,i)=pL/gm+0.5*rL*uL**2
			q(4,i)=rpL
			q(5,i)=rpL*upL
			q(6,i)=rpL*(cc*tpL+0.5*upL**2)
         else
            q(1,i)=rR
            q(2,i)=rR*uR
            q(3,i)=pR/gm+0.5*rR*uR**2
			q(4,i)=rpR
			q(5,i)=rpR*upR
			q(6,i)=rpR*(cc*tpR+0.5*upR**2)
         end if
      end do

      q(1:6,0)=q(1:6,1)
      q(1:6,n)=q(1:6,n-1)
        
      time=0.
      tout=0.
      ic=0
      ip=10
      write(ip,'(e16.8)') (0.5*(x(i)+x(i+1)),i=1,n-1)
	  ip=ip+1
	  write(ip,'(6e16.8)') (q(1,i)/rL,fnu(i)/sqrt(gam*pL/rL),fnp(i)/pL, &
	      q(4,i)/rpL,fnup(i)/sqrt(gam*pL/rL),fntp(i)/tpL,i=1,n-1)	  
	  tout=toutput
!
! Main Loop
!
      do while(time.lt.tend)
         umax=0.
         do i=1,n-1
            umax=max(umax,abs(fnu(i))+fna(i),abs(fnup(i)))
         end do
         dt=cfl*dx/umax

         time=time+dt
         ic=ic+1

!		adjust the output time
		 if(time.ge.tout) then
			dt=dt-(time-tout)
			time=tout
		 end if
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
         w_R(1,n)=rR
         w_R(2,n)=uR
         w_R(3,n)=pR
		 w_R(4,n)=rpR
		 w_R(5,n)=upR
		 w_R(6,n)=tpR

         w_L(1,1)=rL
         w_L(2,1)=uL
         w_L(3,1)=pL
         w_L(4,1)=rpL
         w_L(5,1)=upL
         w_L(6,1)=tpL
!
! Numerical Flux Evaluation  : FDS Scheme
!
         do i=1, n
! gas phase numerical flux
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

            f(1,i)=0.5*(f1L+f1R-d1)
            f(2,i)=0.5*(f2L+f2R-d2)
            f(3,i)=0.5*(f3L+f3R-d3)
!
! particle phase numerical flux
!
			upL$=max(0.,w_L(5,i))
			upR$=min(0.,w_R(5,i))
			rpL$=w_L(4,i)
			rpR$=w_R(4,i)
			epL$=cc*w_L(6,i)+0.5*w_L(5,i)**2
			epR$=cc*w_R(6,i)+0.5*w_R(5,i)**2

!			if(upL.gt.0.and.upR.gt.0.and.upL.lt.upR) then
!				rpH=rpL
!				upH=upL
!			else if(upL.lt.0.and.upR.gt.0.) then
!				rpH=0.
!			else if(upL.lt.0.and.upR.lt.0.and.upL.lt.upR) then
!				rpH=rpR
!				upH=upR
!			else if(upL.gt.0..and.upR.gt.0..and.upL.gt.upR) then
!				rpH=rpL
!				upH=upL
!			else if(upL.gt.0..and.upR.lt.0.) then
!				rpH=rpL+rpR
!				upH=(rpL*upL+rpR*upR)/(rpL+rpR)
!			else
!				rpH=rpR
!				upH=upR
!			end if

			if(upL$.ne.0..or.upR$.ne.0.) then
				rpH=(rpL$*upL$+rpR$*upR$)**2/(rpL$*upL$**2+rpR$*upR$**2)
				upH=(rpL$*upL$**2+rpR$*upR$**2)/(rpL$*upL$+rpR$*upR$)
				epH=(epL$*rpL$*upL$+epR$*rpR$*upR$)/(rpL$*upL$+rpR$*upR$)
			else
				rpH=0.5*(rpL$+rpR$)
				upH=0.
				epH=0.5*(epL$+epR$)
			end if

            f(4,i)=rpH*upH
            f(5,i)=rpH*upH**2
            f(6,i)=rpH*upH*epH

         end do
!
! Source Term
!
		 do i=1, n-1
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

			source(1:6,i)=0.

			source(3,i)=-q(4,i)*cc*(fnt(i)-fntp(i))/tau_t

            temp_g=fnt(i)
			temp_p=fntp(i)
			if(u_rel.ne.0.) then
   			   tratio=temp_p / temp_g
!			   write(6,*) tratio
			   func_h=5.6/(1.+amach)+1.7*sqrt(tratio)

               CD=2.+(CD0-2.)*Exp(-3.07*Sqrt(gam)*amach/reynolds*func_g(reynolds)) &
			     +func_h/(sqrt(gam)*amach)*exp(-0.5*reynolds/amach)
			   tau_v=4./3.*sigmap*dp**2/(reynolds*amyu*CD)

			   source(2,i)=-q(4,i)*u_rel/tau_v
			   source(3,i)=source(3,i)-q(4,i)*u_rel*up/tau_v
			   source(5,i)=-source(2,i)
			else
			   tau_v=4./3.*sigmap*dp**2/(24.*amyu)
			end if

 		    source(6,i)=-source(3,i)

			akai=q(4,i)/q(1,i)
			rhs(1:6)=-dt/dx*(f(1:6,i+1)-f(1:6,i))+dt*source(1:6,i)
			dq(1:6)=rhs(1:6)

			if(imp.eq.1) then
			dq(2)=(rhs(2)*(dt+tau_v)+dt*(rhs(5)-rhs(4)*u+rhs(1)*u*akai)) / &
			      (dt+tau_v+dt*akai)
			dq(5)=(rhs(4)*u*dt+rhs(5)*tau_v+(rhs(2)+rhs(5)-rhs(1)*u)*dt*akai) / &
			      (dt+tau_v+dt*akai)
			dsdq31=-dt*(((cv*temp_g - u**2/2.)*akai*akappa)/tau_t + (u*up*akai)/tau_v)
			dsdq32=-dt*((u*akappa*akai)/tau_t - (up*akai)/tau_v)
			dsdq33=1.-dt*(-(akappa*akai)/tau_t)
			dsdq34=-dt*(-((-up**2/2. + cv*temp_g*akappa)/tau_t) - up**2/tau_v)
			dsdq35=-dt*(-(up/tau_t) + (-u + 2*up)/tau_v)
			dsdq36=-dt/tau_t

			dsdq61= dt*(((cv*temp_g - u**2/2.)*akai*akappa)/tau_t + (u*up*akai)/tau_v)
			dsdq62= dt*((u*akappa*akai)/tau_t - (up*akai)/tau_v)
			dsdq63= dt*(-(akappa*akai)/tau_t)
			dsdq64= dt*(-((-up**2/2. + cv*temp_g*akappa)/tau_t) - up**2/tau_v)
			dsdq65= dt*(-(up/tau_t) + (-u + 2*up)/tau_v)
			dsdq66=1.+dt/tau_t

			rhs(3)=rhs(3)-(dsdq31*dq(1)+dsdq32*dq(2)+dsdq34*dq(4)+dsdq35*dq(5))
			rhs(6)=rhs(6)-(dsdq61*dq(1)+dsdq62*dq(2)+dsdq64*dq(4)+dsdq65*dq(5))

			dq(3)= (dsdq66*rhs(3) - dsdq36*rhs(6))/ (-dsdq36*dsdq63 + dsdq33*dsdq66)
			dq(6)= (dsdq63*rhs(3) - dsdq33*rhs(6))/ ( dsdq36*dsdq63 - dsdq33*dsdq66)
		    end if
			q(1:6,i)=q(1:6,i)+sign(1.,dq(1:6))*max(1.e-80,abs(dq(1:6)))
!			q(1:6,i)=q(1:6,i)+dq(1:6)
			
		end do

		q(1:6,0)=q(1:6,1)
		q(1:6,n)=q(1:6,n-1)
!!
! output
!
         if(time.ge.tout) then
            tout=tout+toutput
            ip=ip+1
            write(6,*) 'Time=',time
	        write(ip,'(6e16.8)') (q(1,i)/rL,fnu(i)/sqrt(gam*pL/rL),fnp(i)/pL, &
			    q(4,i)/rpL,fnup(i)/sqrt(gam*pL/rL),fntp(i)/tpL,i=1,n-1)
         end if
      end do

      stop
	end program ProgM1
