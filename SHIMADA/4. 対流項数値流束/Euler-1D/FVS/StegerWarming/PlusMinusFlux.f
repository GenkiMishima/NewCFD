      program plusminusflux
      implicit real*8 (a-h,o-z)
      real*8 mach, m
      parameter (n=101)
      dimension mach(n), fp(3,n), fm(3,n)
*
      g=1.4
*
      ammax=2.
      ammin=-2.
      dm=(ammax-(ammin))/float(n-1)
      do i=1,n
         mach(i)=ammin+dm*float(i-1)
         m=mach(i)
         fp(1,i)=2.*g*m+abs(1.-m)+2.*(g-1.)*abs(m)+abs(1.+m)
         fp(2,i)=2.+2.*g*m**2+(m-1.)*abs(1.-m)+2.*m*(g-1.)*abs(m)
     &          +(1.+m)*abs(1.+m)
         fp(3,i)=(2.+(m-2)*m*(g-1.))*abs(1.-m)
     &          +2.*m*((2.+m**2*(g-1.))*g+m*(g-1.)**2*abs(m))
     &          +(2.+m*(2.+m)*(g-1.))*abs(m+1.)
*
         fm(1,i)=2.*g*m-abs(-1.+m)-2.*(g-1.)*abs(m)-abs(1.+m)
         fm(2,i)=2.+2.*g*m**2-(m-1.)*abs(-1.+m)-2.*m*(g-1.)*abs(m)
     &          -(1.+m)*abs(1.+m)
         fm(3,i)=(-2.-(m-2)*m*(g-1.))*abs(-1.+m)
     &          +2.*m*((2.+m**2*(g-1.))*g-m*(g-1.)**2*abs(m))
     &          -(2.+m*(2.+m)*(g-1.))*abs(m+1.)
      end do
*
      write(50,'(7e16.8)') (mach(i),fp(1,i),fp(2,i),fp(3,i)
     &                             ,fm(1,i),fm(2,i),fm(3,i),i=1,n)
      stop
      end
        
