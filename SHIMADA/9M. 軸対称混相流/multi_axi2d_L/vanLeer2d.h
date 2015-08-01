!		Header File : vanLeer2d.h
!
!		necessary local arrays are
!			fL(4), fR(4)
!
				HL=gam/gm*pL/rL+0.5*(uL**2+vL**2)
				HR=gam/gm*pR/rR+0.5*(uR**2+vR**2)
				aL=sqrt(gam*pL/rL)
				aR=sqrt(gam*pR/rR)
!
				emL=uL/aL
				emR=uR/aR
				emyL=vL/aL
				emyR=vR/aR
!
				if(emL.lt.-1.) then
					fL(1:4)=0.
				else if(emL.lt.1.) then
					fL(1)=0.25*(1.+emL)**2*rL*aL
					fL(2)=fL(1)*(2.+gm*emL)/gam*aL
					fL(3)=fL(1)*emyL*aL
					fL(4)=fL(1)*(emL+0.5*(emyL**2-1.)+1./gm)*aL**2
				else
					fL(1)=rL*uL
					fL(2)=rL*uL**2+pL
					fL(3)=rL*uL*vL
					fL(4)=rL*uL*HL
				end if
!
				if(1.le.emR) then
					fR(1:4)=0.
				else if(-1.le.emR) then
					fR(1)=-0.25*(-1.+emR)**2*rR*aR
					fR(2)=fR(1)*(-2.+gm*emR)/gam*aR
					fR(3)=fR(1)*emyR*aR
					fR(4)=fR(1)*(-emR+0.5*(emyR**2-1.)+1./gm)*aR**2
				else
					fR(1)=rR*uR
					fR(2)=rR*uR**2+pR
					fR(3)=rR*uR*vR
					fR(4)=rR*uR*HR
				end if
	
!
