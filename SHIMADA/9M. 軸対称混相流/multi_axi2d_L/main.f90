		program main
!
!		Multi-phase 		
!		Two-Dimensional Axisymmetric Euler Solver
!		for arbitrary geometry, structured grid
!			
!			LU-SGS time integration with local-time step used
!			Finite Volume Method
!			MUSCLE approach for higher-order spatial approximation
!				Limiting at most
!
!			Copyright 2003, Toru Shimada
!
		include 'common.h'
		real rtime, ta(2)
!
		call init_io
		call get_data
!
		call geomet
!
		call init_param
		call init_field
!
! Main Loop
!
		do while (ic.le.maxiter)
			rtime=etime(ta)
			start_time=ta(1)
!
			call time_step
!
			call I_boundary
			call I_Flux
!
			call J_boundary
			call J_Flux
!
			call Source
!
			call Point_Implicit
!
			call LU_SGS
!
			call update
!
!			if(mod(ic,ipout).eq.0) 	call out_field(1)
!
			call show_time_residual		
		end do
!
		call out_field(1)
		call out_field(2)

		close(50)
!
		stop
		end