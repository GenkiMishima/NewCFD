	program area_distribution
	implicit real*8 (a-h,o-z)
	parameter (n=151)
	dimension x(n), y(n), a(n)
!
	pi=atan2(0.,-1.)
	xmin=0.
	xmax=3.
	dx=(xmax-xmin)/float(n-1)
	do i=1, n
		x(i)=xmin+dx*float(i-1)
		a(i)=1.+2.2*(x(i)-1.5)**2
		y(i)=sqrt(a(i)/pi)
	end do
!
	open (10,file='area.txt',status='unknown')
		write(10,'(3e16.8)') (x(i),y(i),a(i),i=1,n)
	close(10)
!
	stop
	end