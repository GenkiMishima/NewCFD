!****************************************************************************
!
!  PROGRAM: prog7
!
! 1-D Euler Solver
! Roe's Scheme : Flux Difference Splliting Scheme for Euler Equations
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
! T.Shimada 00/12/11
!
!****************************************************************************
      program prog7

      implicit real*8 (a-h,o-z)
      real*8 kappa
      parameter(n=501)
!
      dimension x(n),q(3,0:n),f(3,n),w_L(3,n),w_R(3,n)
      dimension dw_b(3), dw_f(3), w_c(3)
      namelist /input/ CFL, tend, toutput, phi, kappa, delta
      namelist /flow/ rL,uL,pL,rR,uR,pR
	  namelist /space/ xmin, xmax, xdia
! 
! Definition of statement functions
!
      fnp(i)=gm*(q(3,i)-0.5*q(2,i)**2/q(1,i))
      fnu(i)=q(2,i)/q(1,i)
      fna(i)=sqrt(gam*fnp(i)/q(1,i))
	  fnM(i)=fnu(i)/fna(i)
!
! Set constants
!
      gam=1.4
      gm=gam-1.

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

	  xmin=-1.
	  xmax=1.
	  xdia=0.
!
! read input data
!
      open(99,file='input_prog7.txt',status='unknown')
      read(99,input)
      read(99,flow)
	  read(99,space)
      close(99)

      write(6,input)
      write(6,flow)
	  write(6,space)
      write(60,input)
      write(60,flow)
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
            q(3,i)=pL/(gam-1.)+0.5*rL*uL**2
         else
            q(1,i)=rR
            q(2,i)=rR*uR
            q(3,i)=pR/(gam-1.)+0.5*rR*uR**2
         end if
      end do

      q(1,0)=q(1,1)
      q(2,0)=q(2,1)
      q(3,0)=q(3,1)
      q(1,n)=q(1,n-1)
      q(2,n)=q(2,n-1)
      q(3,n)=q(3,n-1)
        
      time=0.
      tout=0.
      ic=0
      ip=10
      write(ip,'(e16.8)') (0.5*(x(i)+x(i+1)),i=1,n-1)
	  ip=ip+1
	  write(ip,'(4e16.8)') (q(1,i),fnu(i),fnp(i),fnM(i),i=1,n-1)
	  tout=toutput
!
! Main Loop
!
      do while(time.lt.tend)
         umax=0.
         do i=1,n-1
            umax=max(umax,abs(fnu(i))+fna(i))
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
            dw_f(1) = q(1,i+1) - q(1,i  )
            dw_f(2) = fnu(i+1) - fnu(i  )
            dw_f(3) = fnp(i+1) - fnp(i  )
            w_c(1) = q(1,i)
            w_c(2) = fnu(i)
            w_c(3) = fnp(i)

            do j=1, 3
               del_w_r=0.5*((1.+kappa)*dw_b(j)+(1.-kappa)*dw_f(j))
               del_w_l=0.5*((1.-kappa)*dw_b(j)+(1.+kappa)*dw_f(j))

               del_w_lim=0.
               if(dw_b(j).ne.0.) then
                  theta=dw_f(j)/dw_b(j)
                  if(theta.gt.0.) then
                     dl=del_w_l*min(4.*theta/(1.-kappa+(1.+kappa)*theta),1.)
                     dr=del_w_r*min(4./(1.+kappa+(1.-kappa)*theta),1.)
                     del_w_lim=0.5*(sign(1.,dl)+sign(1.,dr))*min(abs(dl),abs(dr))
					 del_w_lim=0.5*(dl+dr)
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

         w_L(1,1)=rL
         w_L(2,1)=uL
         w_L(3,1)=pL
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

            f(1,i)=0.5*(f1L+f1R-d1)
            f(2,i)=0.5*(f2L+f2R-d2)
            f(3,i)=0.5*(f3L+f3R-d3)
         end do
!
! Update field
!
         do j=1, 3
         do i=1, n-1
            q(j,i)=q(j,i)+dt/dx*(f(j,i)-f(j,i+1))
         end do
         end do
!
! output
!
         if(time.ge.tout) then
            tout=tout+toutput
            ip=ip+1
            write(6,*) 'Time=',time
            write(ip,'(4e16.8)') (q(1,i),fnu(i),fnp(i),fnu(i)/fna(i),i=1,n-1)
         end if
      end do

      stop
	end program Prog7