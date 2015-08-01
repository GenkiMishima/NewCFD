!*
!* include file for AUSM_DV_ef3
!*
!* version 1.0 by t.shimada, August 23, 1994
!* version 1.1 by t.shimada, August 31, 1994
!*
!*     RL$,UL$,VL$,EIL$, RR$,UR$,VR$,EIR$ are given
!*     CS, SN, GAM, GM are given
!*
!*     Returns fn1$, fn2$, fn3$, fn4$.
!*
!*
			eiL=pL/rL/gm
            aL=sqrt(gam*pL/rL)
            enL=EIL+.5*(uL**2+vL**2)
            HL=enL+pL/rL
!
			eiR=pR/rR/gm
            aR=sqrt(gam*pR/rR)
            enR=EIR+.5*(uR**2+vR**2)
            HR=enR+pR/rR
!
            cm$=max(aL,aR)
!*
            aml$=uL/cm$
            amr$=uR/cm$
			gconst=sqrt(0.5*gam)
!*====
            arg=gconst*aml$
            b$$=abs(arg)
            ft4$=0.5*(1.+sign(1.,4.-b$$))
            c$$=exp(-b$$*b$$)
            t$$=1./(1.+.3275911*b$$)
            erf$$=sign(1.,arg)     &
	           *( 1.+ ft4$*(     &
              -(.254829592*t$$-.284496736*t$$*t$$ &
             +1.421413741*t$$**3 &
            -1.453152027*t$$**4+1.061405429*t$$**5)*c$$ ) )
!*====
            f_erf_l=1.+erf$$
            exp_l=c$$
!*====
            arg=gconst*amr$
            b$$=abs(arg)
            ft4$=0.5*(1.+sign(1.,4.-b$$))
            c$$=exp(-b$$*b$$)
            t$$=1./(1.+.3275911*b$$)
            erf$$=sign(1.,arg)     &           
			*( 1.+ ft4$*(     &           
			-(.254829592*t$$-.284496736*t$$*t$$     &           
			+1.421413741*t$$**3      &           
			-1.453152027*t$$**4+1.061405429*t$$**5)*c$$ ) )
!*====
            f_erf_r=1.-erf$$
            exp_r=c$$
!*
            prl_p$=0.5*pL*f_erf_l
            prr_m$=0.5*pR*f_erf_r
            p_12$ = prl_p$ + prr_m$
            pu_12$= uL*prl_p$ + uR*prr_m$
!*
            factor=1.
!c           factor=gconst
!            factor=0.78
			pi=atan2(0.,-1.)
			pi_const=sqrt(2./(gam*pi))
            alpha_l$=factor*gam*pL/rL/cm$*pi_const
            alpha_r$=factor*gam*pR/rR/cm$*pi_const
!*
!cL          alpha2_l$= 0.5*pL*cm$*pi_const / factor
!cL          alpha2_r$= 0.5*pR*cm$*pi_const / factor
!*
            unl_p$=0.5*(uL*f_erf_l+alpha_l$*exp_l)
            unr_m$=0.5*(uR*f_erf_r-alpha_r$*exp_r)
            ru_12$=unl_p$*rL+unr_m$*rR
!*
            s$=0.5*min(1.,10.*abs(pR-pL)/min(pL,pR))
!*
            fn1$=ru_12$
            f2_D=0.5*(ru_12$*(uL+uR)-abs(ru_12$)*(uR-uL))
            f2_V=unl_p$*rL*uL+unr_m$*rR*uR
            fn2$=(0.5+s$)*f2_V+(0.5-s$)*f2_D + p_12$
!c           fn2$=f2_V + p_12$
            fn3$=0.5*(ru_12$*(vL+vR)-abs(ru_12$)*(vR-vL))
            fn4$=0.5*(ru_12$*(HL+HR)-abs(ru_12$)*(HR-HL))
!cL          fn4$=0.5*(ru_12$*(enl$+enr$)-abs(ru_12$)*(enr$-enl$ ))
!cL   &          + pu_12$
!cL   &          +0.5*alpha2_l$*exp_l
!cL   &          -0.5*alpha2_r$*exp_r
!cL  This term (enl,enr) are not working good for double_ellipse
!cL  calculation.   August 31, 1994.
!cL So, this version is latest as ausm_DV_ef3.h
