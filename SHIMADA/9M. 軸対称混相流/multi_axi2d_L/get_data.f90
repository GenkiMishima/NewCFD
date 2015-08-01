		Subroutine get_data
!
!			Copyright 2001, Toru SHIMADA
!
		include 'common.h'
!
		namelist /input/ CFL, maxiter, ipiter, phi, kappa, delta, itype, ilocal
		namelist /gas/ amol, gam, amyu, rho_amb, t_amb, pres_amb, u_amb, v_amb, &
					rho_chamber, t_chamber, pres_chamber, u_chamber, v_chamber
		namelist /particle/ sigmap, dp, cc, aloadratio, t_p_amb, u_p_amb, v_p_amb, &
						t_p_chamber, u_p_chamber, v_p_chamber
!
! Set Default values
!
      gam=1.211
	  amol=20.33					! g / mole
      amyu=7.e-5					! Pa-s
!
!     particle phase
!
      sigmap=3204.					! kg/m^3
	  dp = 1.e-6					! m
	  cc = 1380.					! J/kg-K
	  aloadratio=0.4


		read(10,input)
		read(10,gas)
		read(10,particle)
!
		write(6,input)
		write(6,gas)
		write(6,particle)

		write(60,input)
		write(60,gas)
		write(60,particle)
!
		if(itype.eq.0) then
			open(70,file='res.data',status='unknown')
		else
			open(70,file='res.data',status='unknown',access='APPEND')
		end if
!
		return
		end