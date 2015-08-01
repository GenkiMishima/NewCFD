		subroutine init_param
!
		include 'common.h'
!
		gasc=8314.511212/amol
		gm=gam-1.
!
		cp=gam/gm*gasc
		cv=1./gm*gasc
		Prandtl=4.*gam/(9.*gam-5.)
		alambda=amyu*cp/Prandtl
		beta=cc/cv

		gam_eff=(gam+aloadratio*beta)/(1.+aloadratio*beta)
		amol_eff=(1.+aloadratio)*amol
		gasc_eff=8314.511212/amol_eff
!
		totalh=gam_eff*gasc_eff*t_chamber/(gam_eff-1.)
		cdis=gam_eff/sqrt(totalh*(gam_eff-1.))* &
		     (2./(gam_eff+1.))**(0.5*(gam_eff+1.)/(gamm_eff-1.))
		un_in=athroat/aburn*gasc_eff*t_chamber*cdis/(1.+aloadratio)
!
		time = 0.
		ic = 0
		ipout = ipiter
!
		rho_amb    = pres_amb     / (gasc_eff * t_amb)
		rho_chamber= pres_chamber / (gasc_eff * t_chamber)
!
		rho_p_amb = rho_amb * aloadratio
		rho_p_chamber = rho_chamber * aloadratio
!
		return
		end