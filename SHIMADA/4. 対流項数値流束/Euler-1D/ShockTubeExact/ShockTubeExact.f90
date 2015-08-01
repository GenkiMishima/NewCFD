!===============================
!  Program 8
!  Shock Tube Analytic Solution
!
!  T.Shimada 00/12/7
!===============================
!
      program Prog8
	  implicit real*8 (a-h,o-z)
      parameter (nex=101, nmax=200)
      real*8 m4, m1, IM3, IM4, IPlus
      dimension x(nmax),pr(nmax),ve(nmax),ro(nmax)
!
      data r4,u4,p4 /5.,0.,5./
      data r1,u1,p1 /1.,0.,1./
      data time /0.5/
      data xmin, xmax, xdia /-1.0, 1.0, 0./

      open(unit=11,file='pres25.dat')
      open(unit=12,file='dens25.dat')
      open(unit=13,file='velo25.dat')
      open(unit=14,file='mach25.dat')
	  
      gam=1.4
      gm=gam-1.
      gp=gam+1.

      a4=sqrt(gam*p4/r4)
      a1=sqrt(gam*p1/r1)

      m4=u4/a4
      p41=p4/p1
!
! evaluate shock speed
!
      m1=2.
      do loo=1, 100
         term1 = (2.*gam*m1**2-gm)/gp
         term2 = 1.+0.5*gm*(m4-a1/a4*(u1/a1+2./gp*(m1**2-1.)/m1))
         f     = p41-term1*term2**(-2.*gam/gm)
		 dt1dm = 4.*gam/gp*m1
		 dt2dm = -a1*(1.+m1**2)*gm/(a4*m1**2*gp)
		 df    = term2**(-(3.*gam+5.)/gm)*(-gm*term2*dt1dm+2.*gam*term1*dt2dm)/gm
         m1    = m1-f/df
      end do
      u_shock=u1+m1*a1
!
      p21=(1.-gam+2.*m1**2*gam)/(1.+gam)
      p2=p1*p21
      p3=p2
!
      r2=r1*m1**2*gp/(2.+m1**2*gm)
      u2=2.*a1/gp*(m1-1./m1)+u1
      u3=u2
!
      a3=a4+.5*gm*(u4-u3)
      r3=gam*p3/a3**2
!
! expansion fan
!
      u_ex_r=u3-a3
      u_ex_f=u4-a4
!
      x_ex_f=xdia+time*u_ex_f
      x(1)=xmin
      x(2)=x_ex_f
      ro(1)=r4
      pr(1)=p4
      ve(1)=u4
      ro(2)=r4
      pr(2)=p4
      ve(2)=u4
!
      IM3=2.*a3/gm-u3
      IM4=2.*a4/gm-u4
      IPlus=2.*a4/gm+u4
!
!    IM3 and IM4 are the (minus) Rieman Invariant, I-, in region 3 and region 4
!    IPlus is the (plus) Rieman Invariant, I+, constant across the expansion fan
!
      dq=(IM3-IM4)/float(nex-1)
      do i=1, nex
         q=IM4+dq*float(i-1)
         ic=i+2
         u=0.5*(IPlus-q)
         a=0.25*gm*(IPlus+q)
         u_ex=u-a
         x(ic)=xdia+u_ex*time

		 pr(ic)=p4*(a/a4)**(2.*gam/gm)
         ro(ic)=gam*pr(ic)/a**2
         ve(ic)=u
      end do
!
! contact surface
!
      ic=ic+1
      x(ic)=xdia+time*u3
      pr(ic)=p3
      ro(ic)=r3
      ve(ic)=u3
      ic=ic+1
      x(ic)=xdia+time*u3
      pr(ic)=p2
      ro(ic)=r2
      ve(ic)=u2
!
! shock wave
!
      ic=ic+1
      x(ic)=xdia+time*u_shock
      pr(ic)=p2
      ro(ic)=r2
      ve(ic)=u2
      ic=ic+1
      x(ic)=xdia+time*u_shock
      pr(ic)=p1
      ro(ic)=r1
      ve(ic)=u1
!
      ic=ic+1
      x(ic)=xmax
      pr(ic)=p1
      ro(ic)=r1
      ve(ic)=u1
!
      write(11,'(2e16.8)') (x(i),pr(i),i=1,ic)
      write(12,'(2e16.8)') (x(i),ro(i),i=1,ic)
      write(13,'(2e16.8)') (x(i),ve(i),i=1,ic)
      write(14,'(2e16.8)') (x(i),ve(i)/sqrt(gam*pr(i)/ro(i)),i=1,ic)

      close(11)
      close(12)
      close(13)
      close(14)
      stop
      end
