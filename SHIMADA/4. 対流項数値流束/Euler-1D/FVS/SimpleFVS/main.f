      Program Simple_FVS
      implicit real*8 (a-h,o-z)
      real*8 kappa
      parameter(n=101)
*
      dimension x(n),q(3,0:n),f(3,n)
      dimension r_L(n),u_L(n),p_L(n)
      dimension r_R(n),u_R(n),p_R(n)
*
      fnp(i)=(gam-1.)*(q(3,i)-0.5*q(2,i)**2/q(1,i))
      fnu(i)=q(2,i)/q(1,i)
      fnt(i)=fnp(i)/q(1,i)
*
      xmin=-50.
      xmax=50.
*
      dx=(xmax-xmin)/float(n-1)
      do i=1, n
         x(i)=xmin+dx*float(i-1)
      end do
*
      gam=1.4
      rL=1.
      uL=5.
      pL=1.
      rR=1.
      uR=-5.
      pR=1.
*
      xhalf=0.5*(xmin+xmax)
      do i=1, n-1
         xh=0.5*(x(i)+x(i+1))
         if(xh.lt.xhalf) then
            q(1,i)=rL
            q(2,i)=rL*uL
            q(3,i)=pL/(gam-1.)+0.5*rL*uL**2
         else
            q(1,i)=rR
            q(2,i)=rR*uR
            q(3,i)=pR/(gam-1.)+0.5*rR*uR**2
         end if
      end do
*
      q(1,0)=q(1,1)
      q(2,0)=q(2,1)
      q(3,0)=q(3,1)
      q(1,n)=q(1,n-1)
      q(2,n)=q(2,n-1)
      q(3,n)=q(3,n-1)
*
      tend=20.
      cfl=0.1
      kappa=-1.

*        
      time=0.
      ic=0
      ip=10
      write(ip,'(e16.8)') (0.5*(x(i)+x(i+1)),i=1,n-1)
*
      do while(time.lt.tend)
         umax=0.
         do i=1,n-1
            u=q(2,i)/q(1,i)
            a=sqrt((gam-1.)*gam*(-q(2,i)**2+2.*q(3,i)*q(1,i))
     &       /(2.*q(1,i)**2))            
            umax=max(umax,abs(u)+a)
         end do
         dt=cfl*dx/umax
*
         time=time+dt
         ic=ic+1
*
         do i=1, n-1
            r0=q(1,i-1)
            r1=q(1,i  )
            r2=q(1,i+1)
            u0=q(2,i-1)/q(1,i-1)
            u1=q(2,i  )/q(1,i  )
            u2=q(2,i+1)/q(1,i+1)
            p0=-(-1.+gam)*(q(2,i-1)**2-2.*q(3,i-1)*q(1,i-1))
     &         /(2.*q(1,i-1))
            p1=-(-1.+gam)*(q(2,i  )**2-2.*q(3,i  )*q(1,i  ))
     &         /(2.*q(1,i  ))
            p2=-(-1.+gam)*(q(2,i+1)**2-2.*q(3,i+1)*q(1,i+1))
     &         /(2.*q(1,i+1))
*
            dr_b=r1-r0
            dr_f=r2-r1
            del_r_r=0.5*((1.+kappa)*dr_b+(1.-kappa)*dr_f)
            del_r_l=0.5*((1.-kappa)*dr_b+(1.+kappa)*dr_f)
*
            du_b=u1-u0
            du_f=u2-u1
            del_u_r=0.5*((1.+kappa)*du_b+(1.-kappa)*du_f)
            del_u_l=0.5*((1.-kappa)*du_b+(1.+kappa)*du_f)
*
            dp_b=p1-p0
            dp_f=p2-p1
            del_p_r=0.5*((1.+kappa)*dp_b+(1.-kappa)*dp_f)
            del_p_l=0.5*((1.-kappa)*dp_b+(1.+kappa)*dp_f)
*
            del_r_lim=0.
            if(dr_b.ne.0.) then
               theta=dr_f/dr_b
               if(theta.gt.0.) then
                  dl=del_r_l*min(4.*theta/(1.-kappa+(1.+kappa)*theta)
     &               ,1.)
                  dr=del_r_r*min(4./(1.+kappa+(1.-kappa)*theta),1.)
                  del_r_lim=0.5*(sign(1.,dl)+sign(1.,dr))
     &                     *min(abs(dl),abs(dr))
               end if
            end if
*
            del_u_lim=0.
            if(du_b.ne.0.) then
               theta=du_f/du_b
               if(theta.gt.0.) then
                  dl=del_u_l*min(4.*theta/(1.-kappa+(1.+kappa)*theta)
     &               ,1.)
                  dr=del_u_r*min(4./(1.+kappa+(1.-kappa)*theta),1.)
                  del_u_lim=0.5*(sign(1.,dl)+sign(1.,dr))
     &                     *min(abs(dl),abs(dr))
               end if
            end if
*
            del_p_lim=0.
            if(dp_b.ne.0.) then
               theta=dp_f/dp_b
               if(theta.gt.0.) then
                  dl=del_p_l*min(4.*theta/(1.-kappa+(1.+kappa)*theta)
     &               ,1.)
                  dr=del_p_r*min(4./(1.+kappa+(1.-kappa)*theta),1.)
                  del_p_lim=0.5*(sign(1.,dl)+sign(1.,dr))
     &                     *min(abs(dl),abs(dr))
               end if
            end if
*
            r_R(i  )=r1-0.5*del_r_lim
            r_L(i+1)=r1+0.5*del_r_lim
            u_R(i  )=u1-0.5*del_u_lim
            u_L(i+1)=u1+0.5*del_u_lim
            p_R(i  )=p1-0.5*del_p_lim
            p_L(i+1)=p1+0.5*del_p_lim
*
         end do
*
         r_R(n)=rR
         u_R(n)=uR
         p_R(n)=pR
*
         r_L(1)=rL
         u_L(1)=uL
         p_L(1)=pL
*
         do i=1, n
            aL=sqrt(p_L(i)*gam/r_L(i))
            aR=sqrt(p_R(i)*gam/r_R(i))
            emL=u_L(i)/aL
            emR=u_R(i)/aR

            f1p=aL*r_L(i)/(4.*gam)
     &         *(2.*emL*gam+abs(1.-emL)
     &         +2.*(gam-1.)*abs(emL)+abs(1.+emL))
            f1m=aR*r_R(i)/(4.*gam)*(2.*emR*gam-abs(1.-emR)
     &         -2.*(gam-1.)*abs(emR)-abs(1.+emR))

            f2p=aL**2*r_L(i)/(4.*gam)
     &         *(2.+2.*emL**2*gam+(-1.+emL)*abs(1.-emL)
     &         +2.*emL*(gam-1.)*abs(emL)+(1.+emL)*abs(1.+emL))
            f2m=aR**2*r_R(i)/(4.*gam)
     &         *(2.+2.*emR**2*gam-(-1.+emR)*abs(1.-emR)
     &         -2.*emR*(gam-1.)*abs(emR)-(1.+emR)*abs(1.+emR))

            f3p=aL**3*r_L(i)/(8.*(gam-1.)*gam)
     &         *( (+2.+(-2.+emL)*emL*(gam-1.))*abs(1.-emL)
     &         +2.*emL*((2.+emL**2*(-1.+gam))*gam
     &         +emL*(-1.+gam)**2*abs(emL))
     &         +(2.+emL*(2.+emL)*(-1.+gam))*abs(1.+emL))
            f3m=aR**3*r_R(i)/(8.*(gam-1.)*gam)
     &         *( (-2.-(-2.+emR)*emR*(gam-1.))*abs(1.-emR)
     &         +2.*emR*((2.+emR**2*(-1.+gam))*gam
     &         -emR*(-1.+gam)**2*abs(emR))
     &         -(2.+emR*(2.+emR)*(-1.+gam))*abs(1.+emR))
            f(1,i)=f1p+f1m
            f(2,i)=f2p+f2m
            f(3,i)=f3p+f3m

c            f(1,i)=0.5*( r_L(i)*(u_L(i)+aL+abs(u_L(i)))
c     &                  +r_R(i)*(u_R(i)-aR-abs(u_R(i))) )
c
c            f(2,i)=0.5*( p_L(i)+u_L(i)**2*r_L(i)
c     &                  +p_R(i)+u_R(i)**2*r_R(i)
c     &                  +r_L(i)*u_L(i)*(aL+abs(u_L(i)))
c     &                  -r_R(i)*u_R(i)*(aR+abs(u_R(i))) )
c            f(3,i)=0.25*(2.*p_L(i)*u_L(i)*gam/(gam-1.)
c     &                  +2.*p_R(i)*u_R(i)*gam/(gam-1.)
c     &            +u_L(i)**3*r_L(i)+u_R(i)**3*r_R(i)
c     &            +1./(gam-1.)*((2.*p_L(i)+u_L(i)**2*(gam-1.)*r_L(i))
c     &            *(aL+abs(u_L(i))))
c     &            -1./(gam-1.)*((2.*p_R(i)+u_R(i)**2*(gam-1.)*r_R(i))
c     &            *(aR+abs(u_R(i)))))
         end do
*
         do i=1,n-1
            q(1,i)=q(1,i)+dt/dx*(f(1,i)-f(1,i+1))
            q(2,i)=q(2,i)+dt/dx*(f(2,i)-f(2,i+1))
            q(3,i)=q(3,i)+dt/dx*(f(3,i)-f(3,i+1))
         end do
      end do
*
      write(6,*) 'Time=',time
      write(ip+1,'(4e16.8)') (q(1,i),fnu(i),fnp(i),fnt(i),i=1,n-1)
      stop
      end 