		program rigorous2
! Quasi-1d nozzle flow, including a normal shock wave
! rigorous solution
! Copyright 2001 Toru Shimada
!
		implicit real*8 (a-h,o-z)
		real*8 M_in, M_ex, M1
		parameter (n=101)
		dimension x(n), r(n),u(n),p(n), A(n), amach(n)
		fnarea(xx)=1.+2.2*(xx-1.5)**2
!
		pL=147000. ! (Pa)
		pR=9.8e4
		TL=3500. ! (K)
		gam=1.4
		amol=28.96
		gasc=8314.511212/amol
		gm=gam-1.
		gp=gam+1.
		rL=pL/(gasc*TL)

		crit=1.e-10

		xmin=0.
		xmax=3.
		at=fnarea(1.5)

! Solve mach number at inlet

		aratio=fnarea(0.)/at
		em=0.1
		diff=1.e30
		do while(diff.gt.crit) 
			f$=(2./gp*(1.+0.5*gm*em**2))**(gp/gm)/em**2-aratio**2
			dfdm=4.**(gam/gm)/gp*(em**2-1.)/em**3*(0.5*(2.+gm*em**2)/gp)**(2./gm)
			diff=abs(f$/dfdm)
			em=em-f$/dfdm
		end do
		M_in = em
		write(6,*) 'M_in=',M_in
!
! stagnation point value at pre-shock region
		p01=pL*(1.+0.5*gm*M_in**2)**(gam/gm)
		r01=rL*(1.+0.5*gm*M_in**2)**(1./gm)
! Total temperature
		T0 = p01/(r01*gasc)
		write(6,*) 'T0=',T0

! Solve mach number at exit
		Aex = fnarea(3.0)
		M_ex = 0.1
		diff=1.e30
		do while(diff.gt.crit) 
			f$=(Aex*M_ex*pR)**2*(1.+0.5*gm*M_ex**2)**(2.*gam/gm) &
			  -(At*p01)**2*((2.+M_ex**2*gm)/gp)**(gp/gm)
			dfdm=(2.*M_ex*(-(At*p01)**2*((2.+M_ex**2*gm)/gp)**(gp/gm)*gp &
				+(Aex*pR)**2*(1.+0.5*gm*M_ex**2)**(2.*gam/gm)*(2.+M_ex**2*(-1.+3.*gam)))) &
				/(2.+M_ex**2*gm)
			diff=abs(f$/dfdm)
			M_ex=M_ex-f$/dfdm
		end do
		write(6,*) 'M_ex=',M_ex
!
! calculate post-shock values
		p02 = pR*(1.+0.5*gm*M_ex**2)**(gam/gm)
		r02 = p02 / (T0*gasc)
		At2 = p01*At / p02
!
! calculate pre-shock Mach number
		M1 = 2.0
		diff=1.e30
		do while(diff.gt.crit) 
			alpha=((M1**2*gp)/(2.+M1**2*gm))**(gam/gm)
			beta=(gp/(-gm+2.*gam*M1**2))**(1./gm)
			f$=p02/p01 - alpha*beta
			dadm=(4.*gam*((M1**2*gp)/(2.+M1**2*gm))**(2.+1./gm))/(M1**3*(gam**2-1.))
			dbdm=(-4.*M1*gam*(gp/(-gm+2.*gam*M1**2))**(1./gm))/(gm*(1.+(-1.+2.*M1**2)*gam))
			dfdm=-beta*dadm-alpha*dbdm
			diff=abs(f$/dfdm)
			M1=M1-f$/dfdm
		end do
		write(6,*) 'M1=',M1
!
! calculate the shock location
		As = (2.**(gp/(2.*gm))*At*((1. + (M1**2*gm)/2.)/gp)**(gp/(2.*gm)))/M1
		xs = 1.5 + sqrt((As-1.)/2.2)
		write(6,*) 'xs=',xs
!
! calc. all values
        n1 = int (float(n-1)*xs / 3.)+1 
		dx1=(xs-xmin)/float(n1-1)
		do i=1, n1-1
			x(i)=xmin+dx1*float(i-1)
			A(i)=fnarea(x(i))
		end do
		n2 = n1+1

		x(n1)=xs
		x(n2)=xs
		A(n1)=fnarea(xs)
		A(n2)=A(n1)

		dx2 = (xmax-xs) / float(n-n2)
		do i=n2+1, n
			x(i)=xs+dx2*float(i-n2)
			A(i)=fnarea(x(i))
		end do

		do i=1, n
			if(i.le.n1) then
				aratio=a(i)/at
			else
				aratio=a(i)/at2
			end if
!
			if(x(i).lt.1.5) then
				em=0.1
			else if(i.le.n1) then
				em=1.1
			else
				em=0.1
			end if

			diff=1.e30
			do while(diff.gt.crit) 
				f$=(2./gp*(1.+0.5*gm*em**2))**(gp/gm)/em**2-aratio**2
				dfdm=4.**(gam/gm)/gp*(em**2-1.)/em**3*(0.5*(2.+gm*em**2)/gp)**(2./gm)
				diff=abs(f$/dfdm)
				em=em-f$/dfdm
			end do
			amach(i)=em

			if(i.le.n1) then
				ptot=p01
				rtot=r01
			else
				ptot=p02
				rtot=r02
			end if

			p(i)=ptot/(1.+0.5*gm*amach(i)**2)**(gam/gm)
			r(i)=rtot/(1.+0.5*gm*amach(i)**2)**(1./gm)
			u(i)=sqrt(gam*p(i)/r(i))*amach(i)
		end do
!
! output 
		a1=sqrt(gam*p(1)/r(1))
		open (10,file='rigorous2.dat',status='unknown')
		write(10,'(5e16.8)') (x(i),r(i)/r(1),u(i)/a1,p(i)/p(1),amach(i),i=1,n)

		stop
		end