		parameter (ni=101,nj=54)
		implicit real*8 (a-h,o-z)
		real*8 kappa, norm_i, norm_j, nu_a_g, nu_b_g, nu_a_p, nu_b_p
		common /f_array/ q_g(4,0:ni,0:nj), rhs_g(4,ni,nj), &
						 q_p(4,0:ni,0:nj), rhs_p(4,ni,nj), &
						 x(ni,nj), r(ni,nj), rave_i(ni,nj), rave_j(ni,nj), &
						 ds_i(ni,nj), ds_j(ni,nj), area(ni,nj), vol(ni,nj), &
						 norm_i(2,ni,nj), norm_j(2,ni,nj), dt(ni,nj), &
						 alpha_g(ni,nj), alpha_p(ni,nj), &
						 nu_a_g(ni,nj), nu_b_g(ni,nj), nu_a_p(ni,nj), nu_b_p(ni,nj), &
						 tau_v(ni,nj), tau_t(ni,nj)
						 
		common /i_array/ iflag(0:ni,0:nj)
		common /i_var/ ic,maxiter,ipiter,ipout, itype, ilocal
		common /f_var/ cfl,phi,kappa, delta, amol, gam, gm, rho_amb, t_amb, pres_amb, u_amb, v_amb, &
					   rho_chamber, t_chamber, pres_chamber, u_chamber, v_chamber,  &
					   gasc, time, residual, cp, cv, Prandtl, amyu, alambda, sigmap, dp, cc, beta, &
					   aloadratio, rho_p_amb, t_p_amb, u_p_amb, v_p_amb, &
					   rho_p_chamber, t_p_chamber, u_p_chamber, v_p_chamber, &
					   gam_eff, amol_eff, gasc_eff, aburn, athroat,  un_in

		common /cputime/ start_time, end_time, ellapse_time, step_cpu

!