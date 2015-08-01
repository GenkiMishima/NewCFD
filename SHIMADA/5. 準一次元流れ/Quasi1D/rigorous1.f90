		program rigorous1
		implicit real*8 (a-h,o-z)
		parameter (n=101)
		dimension x(n), r(n),u(n),p(n), A(n), amach(n)
		fnarea(xx)=1.+2.2*(xx-1.5)**2
!
		pL=10.e6 ! (Pa)
		TL=3500. ! (K)
		gam=1.4
		amol=28.96
		gasc=8314.511212/amol
		gm=gam-1.
		gp=gam+1.
		rL=pL/(gasc*TL)

		xmin=0.
		xmax=3.
		dx=(xmax-xmin)/float(n-1)
		do i=1, n
			x(i)=xmin+dx*float(i-1)
			A(i)=fnarea(x(i))
		end do
		at=fnarea(1.5)

! Solve mach number

		crit=1.e-5

		do i=1, n
			aratio=a(i)/at
			if(x(i).lt.1.5) then
				em=0.1
			else
				em=1.1
			end if

			diff=1.e30
			do while(diff.gt.crit) 
				f$=(2./gp*(1.+0.5*gm*em**2))**(gp/gm)/em**2-aratio**2
				dfdm=4.**(gam/gm)/gp*(em**2-1.)/em**3*(0.5*(2.+gm*em**2)/gp)**(2./gm)
				diff=abs(f$/dfdm)
				em=em-f$/dfdm
			end do

			amach(i)=em
		end do
!
! stagnation point value
		p0=pL*(1.+0.5*gm*amach(1)**2)**(gam/gm)
		r0=rL*(1.+0.5*gm*amach(1)**2)**(1./gm)
!
! calc. all values

		do i=1, n
			p(i)=p0/(1.+0.5*gm*amach(i)**2)**(gam/gm)
			r(i)=r0/(1.+0.5*gm*amach(i)**2)**(1./gm)
			u(i)=sqrt(gam*p(i)/r(i))*amach(i)
		end do
!
! output 
		a1=sqrt(gam*p(1)/r(1))
		open (10,file='rigorous1.dat',status='unknown')
		write(10,'(5e16.8)') (x(i),r(i)/r(1),u(i)/a1,p(i)/p(1),amach(i),i=1,n)

		stop
		end
