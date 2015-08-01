!		Header File : roe2d.h
!
!		necessary local arrays are
!			delta_q(4), xinv(4,4), alpha(4), aramda(4), XM(4,4), fL(4), fR(4), dissip(4)
!
				HL=gam/gm*pL/rL+0.5*(uL**2+vL**2)
				HR=gam/gm*pR/rR+0.5*(uR**2+vR**2)
!
				rrl=sqrt(rL)
				rrr=sqrt(rR)
				fctL=rrl/(rrl+rrr)
				fctR=rrr/(rrl+rrr)

				uH=fctL*uL + fctR*uR
				vH=fctL*vL + fctR*vR
				HH=fctL*HL + fctR*HR

				aH=sqrt(gm*(HH-0.5*(uH**2+vH**2)))

				delta_q(1)=rR - rL
				delta_q(2)=rR*uR - rL*uL
				delta_q(3)=rR*vR - rL*vL
				delta_q(4)=pR/gm+0.5*rR*(uR**2+vR**2) - pL/gm-0.5*rL*(uL**2+vL**2)

				factor=0.5*gm/aH**2
				xinv(1,1)=-2.*aH**2*vH/gm
				xinv(1,2)=0.
				xinv(1,3)=2.*aH**2/gm
				xinv(1,4)=0.
				xinv(2,1)=-2.*HH+4.*aH**2/gm
				xinv(2,2)=2.*uH
				xinv(2,3)=2.*vH
				xinv(2,4)=-2.
				xinv(3,1)=HH+aH*(uH-aH)/gm
				xinv(3,2)=-uH-aH/gm
				xinv(3,3)=-vH
				xinv(3,4)=1.
				xinv(4,1)=HH-aH*(uH+aH)/gm
				xinv(4,2)=-uH+aH/gm
				xinv(4,3)=-vH
				xinv(4,4)=1.

				xinv(1:4,1:4)=xinv(1:4,1:4) * factor

				alpha=matmul(xinv,delta_q)

				aramda(1)=abs(uH)
				aramda(2)=abs(uH)
				aramda(3)=abs(uH-aH)
				aramda(4)=abs(uH+aH)
!
! Entropy Fix
!
				if(aramda(3).lt.delta*aH) aramda(3)=(aramda(3)**2+delta**2*aH**2)/(2.*delta*aH)
				if(aramda(4).lt.delta*aH) aramda(4)=(aramda(4)**2+delta**2*aH**2)/(2.*delta*aH)
!
				alpha(1)=aramda(1)*alpha(1)
				alpha(2)=aramda(2)*alpha(2)
				alpha(3)=aramda(3)*alpha(3)
				alpha(4)=aramda(4)*alpha(4)
!
				XM(1,1)=0.
				XM(1,2)=1.
				XM(1,3)=1.
				XM(1,4)=1.
				XM(2,1)=0.
				XM(2,2)=uH
				XM(2,3)=uH-aH
				XM(2,4)=uH+aH
				XM(3,1)=1.
				XM(3,2)=vH
				XM(3,3)=vH
				XM(3,4)=vH
				XM(4,1)=vH
				XM(4,2)=0.5*(uH**2+vH**2)
				XM(4,3)=HH-aH*uH
				XM(4,4)=HH+aH*uH
!
				dissip=matmul(XM,alpha)
!
				fL(1)=rL*uL
				fL(2)=rL*uL**2+pL
				fL(3)=rL*uL*vL
				fL(4)=rL*uL*(pL*gam/gm/rL+0.5*(uL**2+vL**2))
!
				fR(1)=rR*uR
				fR(2)=rR*uR**2+pR
				fR(3)=rR*uR*vR
				fR(4)=rR*uR*(pR*gam/gm/rR+0.5*(uR**2+vR**2))
!
				fL(1:4)=0.5*fL(1:4)-0.25*dissip(1:4)
				fR(1:4)=0.5*fR(1:4)-0.25*dissip(1:4)