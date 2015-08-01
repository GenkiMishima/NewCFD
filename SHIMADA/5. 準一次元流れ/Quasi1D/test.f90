		implicit real*8 (a-h,o-z)
		real*8 M_in
		pc=10.e6 ! (Pa)
		tc=3530.
		gam=1.165
		amol=29.66
		gasc=8314.511212/amol
		gm=gam-1.
		gp=gam+1.

		crit=1.e-10


! Solve mach number at inlet

		aratio=1.73
		em=1.0001
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
		pexit=pc/(1.+.5*gm*M_in**2)**(gam/gm)
		temp=tc/(1.+.5*gm*M_in**2)
		dens=pexit/(gasc*temp)
		sos=sqrt(gam*pexit/dens)
		velo=M_in*sos
		write(6,*) 'dens=',dens
		write(6,*) 'velo=',velo
		write(6,*) 'pexit=',pexit
		write(6,*) 'temp=',temp
!
		areaexit=atan2(0.,-1.)*(9.2e-3)**2
		write(6,*) 'areaexit=',areaexit
		flowmass=areaexit*dens*velo
		write(6,*) 'flowmass=',flowmass,'(kg/s)'

		em=1.
		M_in = em
		write(6,*) 'M_in=',M_in


		pexit=pc/(1.+.5*gm*M_in**2)**(gam/gm)
		temp=tc/(1.+.5*gm*M_in**2)
		dens=pexit/(gasc*temp)
		sos=sqrt(gam*pexit/dens)
		velo=M_in*sos
		write(6,*) 'dens=',dens
		write(6,*) 'velo=',velo
		write(6,*) 'pexit=',pexit
		write(6,*) 'temp=',temp
!
		areaexit=atan2(0.,-1.)*(9.2e-3)**2
		write(6,*) 'areaexit=',areaexit
		flowmass=areaexit*dens*velo
		write(6,*) 'flowmass=',flowmass,'(kg/s)'
!
		stop
		end